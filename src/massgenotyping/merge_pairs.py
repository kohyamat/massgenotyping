from __future__ import annotations

import gzip
import os
import shutil
import signal
import tempfile
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from Bio import SeqIO
from natsort import natsorted
from tqdm import tqdm

from .argument_parser import get_args, print_args
from .base import MarkerData, count_records, read_fastx
from .trim_sequence import trim_low_qual, trim_primer
from .utils import run_NGmerge


class MergePairs(MarkerData):
    """
    Merge all paired-end reads in the input directory
    """

    def __init__(
        self,
        indir=None,
        outdir=None,
        trim_quality=True,
        quality_threshold=32,
        window_size=8,
        step_size=1,
        n_cpu=1,
        quiet=True,
        **kwargs,
    ):
        if indir:
            self.indir = Path(indir)
            self.subdirs = natsorted(
                set(
                    [
                        f.parent
                        for f in self.indir.glob("**/*_R1.fastq.gz")
                        if f.parent.name.find("unassigned") < 0
                    ]
                ),
                key=lambda x: x.name,
            )
        if outdir:
            self.outdir = Path(outdir)
        elif indir:
            self.outdir = Path(indir)
        self.trim_quality = trim_quality
        self.quality_threshold = quality_threshold
        self.window_size = window_size
        self.step_size = step_size
        self.n_cpu = n_cpu
        self.quiet = quiet
        if "path_marker_data" in kwargs and kwargs["path_marker_data"]:
            super().__init__(**kwargs)

    def iter_read_pairs(self, subdir):
        read1_list = natsorted(subdir.glob("*_R1.fastq.gz"), key=lambda x: x.name)
        read2_list = natsorted(subdir.glob("*_R2.fastq.gz"), key=lambda x: x.name)

        # Check number of files
        if len(read2_list) == 0:
            errmsg = "Read 2 files did not found in input directory"
            raise RuntimeError(errmsg)
        elif len(read1_list) != len(read2_list):
            errmsg = "Number of read1 files do not same as that of read2 files"
            raise RuntimeError(errmsg)

        for read1, read2 in zip(read1_list, read2_list):
            yield read1, read2

    def quality_trimming(self, read1, read2, **kwargs):
        # make temporarlly files to write trimmed records
        ftmp1, tmppath1 = tempfile.mkstemp()
        ftmp2, tmppath2 = tempfile.mkstemp()
        with os.fdopen(ftmp1, "w") as tmp:
            for rec in read_fastx(read1):
                rec_qt = trim_low_qual(rec, **kwargs)
                SeqIO.write(rec_qt, tmp, "fastq")
        with os.fdopen(ftmp2, "w") as tmp:
            for rec in read_fastx(read2):
                rec_qt = trim_low_qual(rec, **kwargs)
                SeqIO.write(rec_qt, tmp, "fastq")
        return tmppath1, tmppath2

    def merge(self, read_pair):
        read1, read2 = read_pair
        sd = read1.parent
        idx = read1.name.replace("_R1.fastq.gz", "")
        n_reads = count_records(read1)
        if idx != read2.name.replace("_R2.fastq.gz", ""):
            errmsg = "file indexes did not match between read1 and read2"
            raise RuntimeError(errmsg)

        if self.trim_quality:
            qt1, qt2 = self.quality_trimming(
                read1,
                read2,
                threshold=self.quality_threshold,
                window_size=self.window_size,
                step_size=self.step_size,
            )
            run_NGmerge(read1=qt1, read2=qt2, outdir=self.tmpdir, out_prefix=idx)
            # remove temporarily files
            os.remove(qt1)
            os.remove(qt2)
        else:
            run_NGmerge(read1=read1, read2=read2, outdir=self.tmpdir, out_prefix=idx)

        # post mergeing
        of0 = self.tmpdir.joinpath(idx + "_merged.fastq.gz")
        of1 = natsorted(self.tmpdir.glob(idx + "_notMerged*"), key=lambda x: x.name)

        # keep one of the unmerged reads if its quality is not too low
        n_merged = count_records(of0)
        if n_reads > n_merged:
            i = 0
            # if self.frag_len:
            #     frag_len = self.frag_len
            frag_len = []
            if n_merged > 0:
                for rec in read_fastx(of0):
                    frag_len.append(len(rec))
                    i += 1
                    if i == 50:
                        break
            rec_keep = []
            rec1_not_keep = []
            rec2_not_keep = []
            for recs in zip(read_fastx(of1[0]), read_fastx(of1[1])):
                recs[0].description = recs[0].description + "_notMerged1"
                recs[1].description = recs[1].description + "_notMerged2"
                recs[1].seq = recs[1].seq.reverse_complement()
                qual_rev = recs[1].letter_annotations["phred_quality"][::-1]
                recs[1].letter_annotations["phred_quality"] = qual_rev
                pr_matches = np.array([x.description.count(sd.name) for x in recs])
                seq_lens = np.array([len(x) for x in recs])
                qual_means = np.array(
                    [np.mean(x.letter_annotations["phred_quality"]) for x in recs]
                )
                if frag_len:
                    len_th = np.median(frag_len) * 0.8
                else:
                    len_th = 50

                if all(seq_lens > len_th):
                    if all(pr_matches == 2):
                        i = np.argmax(qual_means)
                    elif pr_matches[0] != pr_matches[1]:
                        i = np.argmax(pr_matches)
                    else:
                        i = np.argmax(qual_means)
                elif any(seq_lens > len_th):
                    i = np.argmax(seq_lens)
                else:
                    i = -1

                if i > -1 and qual_means[i] > 25:
                    rec_keep.append(recs[i])
                else:
                    rec1_not_keep.append(recs[0])
                    rec2_not_keep.append(recs[1])

            n_keep = len(rec_keep)
            if n_keep > 0:
                with gzip.open(of0, "at") as outfile:
                    SeqIO.write(rec_keep, outfile, "fastq")

                if len(rec1_not_keep) > 0:
                    with gzip.open(of1[0], "wt") as outfile:
                        SeqIO.write(rec1_not_keep, outfile, "fastq")
                    with gzip.open(of1[1], "wt") as outfile:
                        SeqIO.write(rec2_not_keep, outfile, "fastq")
                else:
                    for f in of1:
                        f.unlink()
        else:
            n_keep = 0
            for f in of1:
                f.unlink()

        return idx, n_reads, n_merged, n_keep

    def merge_all_pairs_in_subdir(self, sd):
        self.tmpdir = Path(tempfile.mkdtemp())
        n_pairs = len(list(sd.glob("*_R1.fastq.gz")))

        if self.dict_frag_len:
            self.frag_len = self.dict_frag_len[sd.name]
        else:
            self.frag_len = None

        if self.n_cpu == 1:
            if self.quiet:
                res = [self.merge(x) for x in self.iter_read_pairs(sd)]
            else:
                res = [
                    self.merge(x) for x in tqdm(self.iter_read_pairs(sd), total=n_pairs)
                ]
        else:
            original_sigint_handler = signal.getsignal(signal.SIGINT)
            pool = Pool(self.n_cpu)
            signal.signal(signal.SIGINT, original_sigint_handler)
            try:
                if self.quiet:
                    res = list(pool.imap(self.merge, self.iter_read_pairs(sd)))
                else:
                    res = list(
                        tqdm(
                            pool.imap(self.merge, self.iter_read_pairs(sd)),
                            total=n_pairs,
                        )
                    )
            except KeyboardInterrupt:
                print("Caught KeyboardInterrupt, terminating workers")
                pool.terminate()
                pool.join()
                shutil.rmtree(str(self.tmpdir))
                raise KeyboardInterrupt
            else:
                pool.close()
                pool.join()

        outpath0 = self.outdir.joinpath(sd.name)
        outpath1 = self.outdir.joinpath("result_summary", "merge_pairs")
        outpath2 = outpath1.joinpath(sd.name + "_merge_summary.csv")
        outpath0.mkdir(parents=True, exist_ok=True)
        outpath1.mkdir(parents=True, exist_ok=True)

        for f in self.tmpdir.iterdir():
            shutil.copy(str(f), str(outpath0))

        with outpath2.open("w") as outfile:
            outfile.write("ID,N_read_pairs,N_merged,N_kept_unmerged,")
            outfile.write("Percent_merged,Percent_remained\n")
            for idx, n0, n1, n2 in natsorted(res, key=lambda x: x[0]):
                outfile.write(
                    "{},{},{},{},{:.1f},{:.1f}\n".format(
                        idx, n0, n1, n2, n1 / n0 * 100, (n1 + n2) / n0 * 100
                    )
                )

        shutil.rmtree(str(self.tmpdir))

        if not self.quiet:
            sums = np.array(np.array(res)[:, 1:3], dtype=int).sum(axis=0)
            p_merged = np.divide(*sums[::-1])
            print("Percentage merged: {:.2%}".format(p_merged))

        return res

    def run(self):
        results = []
        if self.quiet:
            for sd in tqdm(self.subdirs, total=len(self.subdirs)):
                results += self.merge_all_pairs_in_subdir(sd)
        else:
            for sd in self.subdirs:
                print("." * shutil.get_terminal_size()[0])
                print("Processing {} ({})".format(sd.name, time.ctime()))
                results += self.merge_all_pairs_in_subdir(sd)
            print("")

        sums = np.array(np.array(results)[:, 1:3], dtype=int).sum(axis=0)
        p_merged = np.divide(*sums[::-1])
        print("Overall percentage merged: {:.2%}".format(p_merged))


class TrimPrimers(MarkerData):
    """Trim primer sequences from post-merging sequences"""

    def __init__(self, indir, max_mismatch=None, n_cpu=1, **kwargs):
        self.indir = Path(indir)
        self.max_mismatch = max_mismatch
        self.n_cpu = n_cpu
        super().__init__(**kwargs)

    def worker(self, seq_file):
        sd = seq_file.parent
        self.tmpdir.joinpath(sd.name).mkdir(exist_ok=True)
        outpath = self.tmpdir.joinpath(sd.name, seq_file.name)
        try:
            with gzip.open(outpath, "wt") as outfile:
                for record in read_fastx(seq_file):
                    records_pt = trim_primer(
                        record, self.dict_primer[sd.name], self.max_mismatch
                    )
                    SeqIO.write(records_pt, outfile, "fastq")
        except RuntimeError:
            pass

    def run(self):
        self.tmpdir = Path(tempfile.mkdtemp())
        g = "*_merged.fastq.gz"
        n_files = len(list(self.indir.glob("**/{}".format(g))))
        if self.n_cpu == 1:
            for seq_file in tqdm(self.indir.glob("**/{}".format(g)), total=n_files):
                self.worker(seq_file)
        else:
            original_sigint_handler = signal.getsignal(signal.SIGINT)
            pool = Pool(self.n_cpu)
            signal.signal(signal.SIGINT, original_sigint_handler)
            try:
                list(
                    tqdm(
                        pool.imap_unordered(
                            self.worker, self.indir.glob("**/{}".format(g))
                        ),
                        total=n_files,
                    )
                )
            except KeyboardInterrupt:
                print("Caught KeyboardInterrupt, terminating workers")
                pool.terminate()
                pool.join()
                shutil.rmtree(str(self.tmpdir))
                raise KeyboardInterrupt
            else:
                pool.close()
                pool.join()

        for f in self.tmpdir.glob("**/{}".format(g)):
            outpath = self.indir.joinpath(f.parent.name, f.name)
            shutil.copy(str(f), str(outpath))


def main(args):
    print("Start merging paired-end reads ... ({})".format(time.ctime()))
    t1 = time.time()

    if isinstance(args, dict):
        kwargs = args
    else:
        kwargs = vars(args)

    if "subcommand" in kwargs:
        kwargs.pop("subcommand")

    if not kwargs["quiet"]:
        print_args(kwargs)

    mgp = MergePairs(**kwargs)
    mgp.run()

    if kwargs["trim_primer"]:
        print("\nTrimming primer sequences ...")
        if kwargs["outdir"]:
            kwargs["indir"] = kwargs.pop("outdir")
        kwargs["max_mismatch"] = kwargs.pop("max_mismatch_primer")
        trp = TrimPrimers(**kwargs)
        trp.run()

    t2 = time.time()
    elapsed_time = time.strftime("%H:%M:%S", time.gmtime(t2 - t1))
    print("\nFinished (elapsed time: {})".format(elapsed_time))


if __name__ == "__main__":
    import sys

    argv = ["merge-pairs"]
    argv.extend(sys.argv[1:])
    main(get_args(argv))
