from __future__ import annotations

import gzip
import os
import signal
import tempfile
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from Bio import AlignIO, SeqIO
from tqdm import tqdm

from .argument_parser import get_args, print_args
from .base import count_records, count_uniq_seq, read_fastx
from .utils import run_mafft


class Denoise:
    def __init__(
        self,
        infile,
        a=0.1,
        b=0.9,
        outdir=None,
        out_suffix=None,
        localpair=False,
        globalpair=False,
        n_cpu=1,
        quiet=False,
        **kwargs,
    ):
        self.infile = infile
        self.a = a
        self.b = b
        self.outdir = outdir
        if out_suffix:
            self.out_suffix = out_suffix
        else:
            self.out_suffix = ""
        self.localpair = localpair
        self.globalpair = globalpair
        self.n_cpu = n_cpu
        self.quiet = quiet

    @staticmethod
    def denoise(
        seq_file,
        a=0.1,
        b=0.9,
        outdir=None,
        out_suffix="_denoised",
        exclude_outlier_in_length=False,
        localpair=False,
        globalpair=False,
    ):
        idx = Path(seq_file).name.split(".")[0].replace("_merged", "")
        if outdir:
            outpath = Path(outdir).joinpath(idx + out_suffix + ".fasta.gz")
        else:
            outpath = Path(seq_file).parent.joinpath(idx + out_suffix + ".fasta.gz")

        # count number of reads for unique sequences
        seq_count = count_uniq_seq(seq_file)
        if exclude_outlier_in_length:
            seq_count = Denoise.exclude_outliers_in_sequence_length(seq_count)
        n_reads = np.sum(list(seq_count.values()))

        # align unique sequences
        seqs = list(seq_count.keys())

        if len(seqs) == 1:
            with gzip.open(outpath, "wt") as outfile:
                for rec in read_fastx(seq_file):
                    SeqIO.write(rec, outfile, "fasta")
            return

        if localpair or globalpair:
            if localpair:
                opts = "--localpair --maxiterate 1000"
            else:
                opts = "--globalpair --maxiterate 1000"
            if len(seqs) > 50:
                align_tmp = run_mafft(seqs[:50], opts=opts)
                ftmp, tmppath = tempfile.mkstemp()
                try:
                    with os.fdopen(ftmp, "w") as outfile:
                        AlignIO.write(align_tmp, outfile, "fasta")
                    align = run_mafft(seqs[50:], add_to_existence_alignment=tmppath)
                finally:
                    os.remove(tmppath)
            else:
                align = run_mafft(seqs, opts=opts)
        else:
            align = run_mafft(seqs)

        align = np.array(align)
        align0 = np.repeat(align, tuple(seq_count.values()), 0)

        # calculate base frequencies for each position
        base_uniq = np.unique(align)
        base_frq = np.array([(align0 == x).sum(axis=0) for x in base_uniq]) / n_reads

        # the most frequent base for each position
        base_frq_max = base_frq.max(axis=0)
        most_frq = base_uniq[base_frq.argmax(axis=0)]

        # replace rare bases with the most frequent base for each presumably
        # homogeneous nucleotide postion
        J, K = np.where((base_frq > 0) & (base_frq < a) & (base_frq_max > b))
        for j, k in zip(J, K):
            align[:, k][align[:, k] == base_uniq[j]] = most_frq[k]

        # delete positions where the most freqeunt character is '-'
        align = np.delete(
            align, np.where((most_frq == "-") & (base_frq_max > b)), axis=1
        )

        # write denoised sequences in a file (FASTA format)
        mapping = {
            list(seq_count.keys())[i]: "".join(align[i, :]).upper().replace("-", "")
            for i in range(len(align))
        }

        with gzip.open(outpath, "wt") as outfile:
            for rec in read_fastx(seq_file):
                seq = str(rec.seq)
                desc = rec.description
                if seq in mapping.keys():
                    outfile.write(">{}\n{}\n".format(desc, mapping[seq]))

    @staticmethod
    def exclude_outliers_in_sequence_length(seq_count):
        """exclude sequences whose length is out of quantiles ± IQR * 1.5"""
        lengths, n_reads = zip(*[[len(s), c] for s, c in seq_count.items()])
        x = np.repeat(lengths, n_reads)
        if np.sum(n_reads) > 4 and len(lengths) > 4:
            quantiles = [np.percentile(x, q) for q in [25, 75]]
            IQR = quantiles[1] - quantiles[0]
            quartiles_ext = quantiles + np.array([-IQR * 1.5, IQR * 1.5])
            res = {}
            for (k, v), length in zip(seq_count.items(), lengths):
                if length >= quartiles_ext[0] and length <= quartiles_ext[1]:
                    res[k] = v
            return res
        else:
            return seq_count

    def iter_args(self):
        for f in self.infile:
            yield f, self.a, self.b, self.outdir, self.out_suffix

    def worker(self, args):
        try:
            if count_records(args[0]) > 1:
                if self.localpair:
                    self.denoise(*args, localpair=True)
                elif self.globalpair:
                    self.denoise(*args, globalpair=True)
                else:
                    self.denoise(*args)
            else:
                idx = Path(args[0]).name.split(".")[0]
                if args[3]:
                    outpath = Path(args[3]).joinpath(idx + args[4] + ".fasta.gz")
                else:
                    outpath = Path(args[0]).parent.joinpath(idx + args[4] + ".fasta.gz")
                with gzip.open(outpath, "wt") as outfile:
                    for rec in read_fastx(args[0]):
                        SeqIO.write(rec, outfile, "fasta")
        except RuntimeError:
            pass

    def run(self):
        if self.n_cpu == 1:
            if self.quiet:
                for args in self.iter_args():
                    self.worker(args)
            else:
                for args in tqdm(self.iter_args(), total=len(self.infile)):
                    self.worker(args)
        else:
            original_sigint_handler = signal.getsignal(signal.SIGINT)
            pool = Pool(self.n_cpu)
            signal.signal(signal.SIGINT, original_sigint_handler)
            try:
                list(
                    tqdm(
                        pool.imap_unordered(self.worker, self.iter_args()),
                        total=len(self.infile),
                    )
                )
            except KeyboardInterrupt:
                print("Caught KeyboardInterrupt, terminating workers")
                pool.terminate()
                pool.join()
                raise KeyboardInterrupt
            else:
                pool.close()
                pool.join()


def main(args):
    print("Start denoising ... ({})".format(time.ctime()))
    t1 = time.time()

    if isinstance(args, dict):
        kwargs = args
    else:
        kwargs = vars(args)

    if "subcommand" in kwargs:
        kwargs.pop("subcommand")

    if not kwargs["quiet"]:
        print_args(kwargs)

    inpath = kwargs.pop("inpath")
    glob_pattern = kwargs.pop("glob_pattern")

    if len(inpath) == 1 and Path(inpath[0]).is_dir():
        infile = list(Path(inpath[0]).glob("**/" + glob_pattern))
        if not infile:
            raise RuntimeError("No files match with the pattern")
    else:
        infile = inpath

    dns = Denoise(infile=infile, **kwargs)
    dns.run()

    t2 = time.time()
    elapsed_time = time.strftime("%H:%M:%S", time.gmtime(t2 - t1))
    print("\nFinished (elapsed time: {})".format(elapsed_time))


if __name__ == "__main__":
    import sys

    argv = ["denoise"]
    argv.extend(sys.argv[1:])
    main(get_args(argv))
