from __future__ import annotations

import gzip
import os
import re
import shutil
import signal
import subprocess
import tempfile
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from Bio import SeqIO
from fuzzysearch import find_near_matches
from natsort import natsorted
from tqdm import tqdm

from .argument_parser import get_args, print_args
from .base import (
    MarkerData,
    PrimerMatch,
    check_no_wrapped,
    count_records,
    guess_fmt,
    read_fastx,
    revc,
)
from .utils import common_prefix


class Demultiplex(MarkerData):
    def __init__(
        self,
        infile=None,
        outdir="output",
        out_prefix="",
        max_mismatch=[0.12, 0.14],
        n_cpu=1,
        trim_primer=False,
        quiet=False,
        *args,
        **kwargs,
    ):
        if infile:
            self.set_infile(infile)
        self.outdir = Path(outdir)
        self.out_prefix0 = out_prefix
        if self.out_prefix0:
            self.out_prefix0 += "_"
        self.set_out_prefix(out_prefix)
        self.max_mismatch = max_mismatch
        self.n_cpu = n_cpu
        self.trim_primer = trim_primer
        self.quiet = quiet
        self.outdir.mkdir(exist_ok=True)
        super().__init__(*args, **kwargs)
        if isinstance(list(self.dict_primer.values())[0], str):
            self.fwd_only = True
        else:
            self.fwd_only = False

    def set_infile(self, infile):
        if isinstance(infile, list) and len(infile) <= 2:
            self.infile = infile
            self.n_reads = count_records(infile[0])
            if len(infile) == 2:
                self.paired = True
                if self.n_reads != count_records(infile[1]):
                    errmsg = "The number of sequences differs between "
                    errmsg += "{} and {}".format(infile[0], infile[1])
                    raise RuntimeError(errmsg)
            else:
                self.paired = False
        else:
            errmsg = "infile must be a list of 1 or 2 file path(s)"
            raise ValueError(errmsg)

    def set_out_prefix(self, out_prefix):
        if out_prefix:
            if self.out_prefix0:
                self.out_prefix = self.out_prefix0 + out_prefix + "_"
            else:
                self.out_prefix = out_prefix + "_"
        else:
            self.out_prefix = self.out_prefix0

    def assign(self, i, record, pid=""):
        """
        Assign loci to sequences
        """
        if self.paired:
            rec1, rec2 = record
            seq1, seq2 = str(rec1.seq), str(rec2.seq)
        else:
            rec1 = record[0]
            seq1 = str(rec1.seq)

        if i in self.exact_matches1:
            match1 = self.exact_matches1[i]
            for j, m in enumerate(match1):
                if isinstance(m, str):
                    if m:
                        match1[j] = find_approximate_matches(
                            seq1, {m: self.dict_primer[m]}, 0, False, j
                        )
                    else:
                        match1[j] = find_approximate_matches(
                            seq1, self.dict_primer, self.max_mismatch, False, j
                        )
                elif m == -1:
                    match1[j] = find_approximate_matches(
                        seq1, self.dict_primer, self.max_mismatch, False, j
                    )
                else:
                    # if m is PrimerMatch
                    pass
        else:
            match1 = find_approximate_matches(
                seq1, self.dict_primer, self.max_mismatch, False
            )
            if self.fwd_only:
                match1 = [match1]

        match1 = [m if m else PrimerMatch("", -1, -1, 100) for m in match1]

        if self.paired:
            if i in self.exact_matches2:
                match2 = self.exact_matches2[i]
                for j, m in enumerate(match2):
                    if isinstance(m, str):
                        if m:
                            match2[j] = find_approximate_matches(
                                seq2, {m: self.dict_primer[m]}, 0, True, j
                            )
                        else:
                            match2[j] = find_approximate_matches(
                                seq2, self.dict_primer, self.max_mismatch, True, j
                            )
                    elif m == -1:
                        match2[j] = find_approximate_matches(
                            seq2, self.dict_primer, self.max_mismatch, True, j
                        )
                    else:
                        # if m is PrimerMatch
                        pass
            else:
                match2 = find_approximate_matches(
                    seq2, self.dict_primer, self.max_mismatch, True
                )
                if self.fwd_only:
                    match2 = [match2]

            match2 = [m if m else PrimerMatch("", -1, -1, 100) for m in match2]

        # Assign sequences
        if self.fwd_only:
            if self.paired:
                fwd_best = get_best_match([match1[0], match2[0]])
                primer_len = fwd_best.end - fwd_best.start
                short = len(rec1) < (primer_len + 10) and len(rec2) < (primer_len + 10)
            else:
                fwd_best = match1[0]
                primer_len = fwd_best.end - fwd_best.start
                short = len(rec1) < (primer_len + 10)

            if short:
                match = "unassigned_SS"
            elif fwd_best.match:
                match = fwd_best.match
            else:
                match = "unassigned_NF"
            assign = (fwd_best.match, match)
        else:
            if self.paired:
                fwd_best = get_best_match([match1[0], match2[1]])
                rev_best = get_best_match([match1[1], match2[0]])
                primer_len = fwd_best.end - fwd_best.start
                primer_len += rev_best.end - rev_best.start
                short = len(rec1) < (primer_len + 10) and len(rec2) < (primer_len + 10)
            else:
                fwd_best, rev_best = match1
                primer_len = fwd_best.end - fwd_best.start
                primer_len += rev_best.end - rev_best.start
                short = len(rec1) < (primer_len + 10)

            if short:
                match = "unassigned_SS"
            elif fwd_best.match == rev_best.match != "":
                match = fwd_best.match
            # elif match1[0].match == match2[0].match != "":
            #     match = match1[0].match
            # elif match1[0].match == match1[1].match != "":
            #     match = match1[0].match
            # elif match2[0].match == match2[1].match != "":
            #     match = match2[0].match
            # elif match1[1].match == match2[1].match != "":
            #     match = match1[1].match
            elif fwd_best.match == "" or rev_best.match == "":
                match = "unassigned_NF"
            else:
                match = "unassigned_WC"
            assign = (fwd_best.match, rev_best.match, match)

        # Trim primer sequences
        if self.trim_primer:
            if self.fwd_only:
                if match1[0].match != "":
                    rec1 = rec1[match1[0].end :]
                if self.paired and match2[0].match != "":
                    rec2 = rec2[: match2[0].start]
            else:
                if match1[0].match != "" and match1[1].match != "":
                    rec1 = rec1[match1[0].end : match1[1].start]
                elif match1[0].match != "":
                    rec1 = rec1[match1[0].end :]
                elif match1[1].match != "":
                    rec1 = rec1[: match1[1].start]
                if self.paired:
                    if match2[0].match != "" and match2[1].match != "":
                        rec2 = rec2[match2[0].end : match2[1].start]
                    elif match2[0].match != "":
                        rec2 = rec2[match2[0].end :]
                    elif match2[1].match != "":
                        rec2 = rec2[: match2[1].start]

        self.write_record_tmp(assign, record, pid)

        if match.find("unassigned_") > -1:
            return 0
        else:
            return 1

    def write_record_tmp(self, assign, record, pid=""):
        """
        write sequence records to temporally files
        """
        if self.fwd_only:
            fwd_name, match = assign
        else:
            fwd_name, rev_name, match = assign

        if self.paired:
            rec1, rec2 = record
        else:
            rec1 = record[0]

        tmppath0 = self.tmpdir.joinpath("tmp{}_assign.csv".format(pid))
        with tmppath0.open(mode="a") as outfile:
            if len(assign) == 2:
                outfile.write(",".join((rec1.id, fwd_name, match)) + "\n")
            else:
                outfile.write(",".join((rec1.id, fwd_name, rev_name, match)) + "\n")

        tmppath1 = self.tmpdir.joinpath(match, "tmp{}_R1.fastq".format(pid))
        with tmppath1.open(mode="a") as outfile:
            SeqIO.write(rec1, outfile, "fastq")

        if len(record) == 2:
            tmppath2 = self.tmpdir.joinpath(match, "tmp{}_R2.fastq".format(pid))
            with tmppath2.open(mode="a") as outfile:
                SeqIO.write(rec2, outfile, "fastq")

    def iter_args(self):
        for i, record in enumerate(zip(*[read_fastx(x) for x in self.infile])):
            yield i, record

    def worker(self, args):
        pid = os.getpid()
        return self.assign(*args, pid=pid)

    def run(self):
        # make temporally directories
        self.tmpdir = Path(tempfile.mkdtemp())
        for x in self.dict_primer:
            self.tmpdir.joinpath(x).mkdir(exist_ok=True)
        for x in ["unassigned_SS", "unassigned_WC", "unassigned_NF"]:
            self.tmpdir.joinpath(x).mkdir(exist_ok=True)

        # find exact matches
        self.exact_matches1 = find_exact_matches(self.infile[0], self.dict_primer)
        if len(self.infile) == 2:
            self.exact_matches2 = find_exact_matches(self.infile[1], self.dict_primer)

        # find approximate matches
        if self.n_cpu == 1:
            if self.quiet:
                res = [self.worker(args) for args in self.iter_args()]
            else:
                res = [
                    self.worker(args)
                    for args in tqdm(self.iter_args(), total=self.n_reads)
                ]
        else:
            if self.n_reads > 200:
                chunksize = int(self.n_reads / 200)
            else:
                chunksize = 1
            original_sigint_handler = signal.getsignal(signal.SIGINT)
            pool = Pool(self.n_cpu)
            signal.signal(signal.SIGINT, original_sigint_handler)
            try:
                if self.quiet:
                    res = list(
                        pool.map(
                            self.worker,
                            self.iter_args(),
                            chunksize=chunksize,
                        )
                    )
                else:
                    res = list(
                        tqdm(
                            pool.imap_unordered(
                                self.worker,
                                self.iter_args(),
                                chunksize=chunksize,
                            ),
                            total=self.n_reads,
                        )
                    )
            except KeyboardInterrupt:
                print("Caught KeyboardInterrupt, terminating workers")
                pool.terminate()
                pool.join()
                shutil.rmtree(self.tmpdir)
                raise KeyboardInterrupt
            else:
                pool.close()
                pool.join()

        if not self.quiet:
            mode = "paired" if len(self.infile) == 2 else "single"
            n_assined = np.array(res).sum()
            excl = list(self.tmpdir.glob("tmp*_assign.csv"))
            excl.extend(list(self.tmpdir.glob("unassigned_*")))
            n_loci = len(
                [
                    x
                    for x in self.tmpdir.iterdir()
                    if x not in excl and len(list(x.iterdir())) > 0
                ]
            )
            msg = "{} out of {} {}-end reads ({:.2%}) have been assigned to {} loci"
            msg = msg.format(
                n_assined, self.n_reads, mode, n_assined / self.n_reads, n_loci
            )
            print(msg)

        # Combine temporally files
        resdir = self.outdir.joinpath("result_summary", "demultiplex")
        resdir.mkdir(parents=True, exist_ok=True)
        outpath0 = resdir.joinpath("{}assign.csv.gz".format(self.out_prefix))
        ftmp0 = self.tmpdir.glob("tmp*_assign.csv")
        with gzip.open(outpath0, "wb") as outfile:
            for f in ftmp0:
                with f.open(mode="rb") as tmp:
                    shutil.copyfileobj(tmp, outfile, 1024 * 1024 * 10)

        for sd in self.tmpdir.iterdir():
            fname1 = "{}{}_R1.fastq.gz".format(self.out_prefix, sd.name)
            outpath1 = self.outdir.joinpath(sd.name, fname1)
            ftmp1 = natsorted(sd.glob("tmp*_R1.fastq"), key=lambda x: x.name)
            if len(ftmp1) > 0:
                self.outdir.joinpath(sd.name).mkdir(exist_ok=True)
                with gzip.open(outpath1, "wb") as outfile:
                    for f in ftmp1:
                        with f.open(mode="rb") as tmp:
                            shutil.copyfileobj(tmp, outfile, 1024 * 1024 * 10)

            if len(self.infile) == 2:
                fname2 = "{}{}_R2.fastq.gz".format(self.out_prefix, sd.name)
                outpath2 = self.outdir.joinpath(sd.name, fname2)
                ftmp2 = natsorted(sd.glob("tmp*_R2.fastq"), key=lambda x: x.name)
                if len(ftmp2) > 0:
                    with gzip.open(outpath2, "wb") as outfile:
                        for f in ftmp2:
                            with f.open(mode="rb") as tmp:
                                shutil.copyfileobj(tmp, outfile, 1024 * 1024 * 10)

        shutil.rmtree(self.tmpdir)
        return res


def _find_exact_matches(filepath, subseqs, subseq_lookup, fmt):
    filepath = Path(filepath)
    opts = "-n -o"
    if shutil.which("rg"):
        prog = "rg"
        if filepath.suffix == ".gz":
            opts = opts + " -z --column"
        else:
            opts = opts + " --column"
    else:
        if filepath.suffix == ".gz":
            prog = "zgrep"
        else:
            prog = "grep"
        opts = opts + " -E"

    if fmt == "fasta":
        k = 2
    elif fmt == "fastq":
        k = 4
    else:
        raise ValueError("'fmt' must be 'fasta' or 'fastq'")

    if any([x.find(":-)rev") > -1 for x in subseq_lookup.values()]):
        fwd_only = False
    else:
        fwd_only = True

    cmd0 = r"{} {} {} {}".format(prog, opts, subseqs, filepath).split()
    res = subprocess.Popen(cmd0, stdout=subprocess.PIPE)
    matches = {}
    for line in res.stdout.read().decode("utf-8").strip().split():
        splits = line.split(":")
        if len(splits) > 2:
            line_no, column, match = line.split(":")
        elif len(splits) > 1:
            line_no, match = splits
        else:
            match = splits[0]
        rec_no = int((int(line_no) - 2) / k)
        marker, match_type = subseq_lookup[match].split(":-)")
        if fwd_only:
            j = 0
        else:
            j = 1 if match_type.find("rc") > -1 else 0
        if rec_no not in matches:
            if fwd_only:
                matches[rec_no] = [""]
            else:
                matches[rec_no] = [""] * 2
        if matches[rec_no][j]:
            # if multiple matches
            matches[rec_no][j] = -1
        else:
            if len(splits) > 2:
                start = int(column) - 1
                end = start + len(match)
                matches[rec_no][j] = PrimerMatch(marker, start, end, 0)
            else:
                matches[rec_no][j] = marker

    return matches


def find_exact_matches(filepath, dict_subseqs, fmt="auto", read2=False):
    """
    Find exact matches using grep or ripgrep

    Parameters
    -----
    filepath : str or pathlib.PosixPath
        path to the input sequene file
    dict_subseqs : dict
        dict of sub-sequences to search
    fmt : ["fasta", "fastq", "auto"]
        file format (default: "auto")
    """

    if fmt == "auto":
        fmt = guess_fmt(filepath)
    elif fmt not in ["fasta", "fastq"]:
        raise ValueError("'fmt' must be 'fasta', 'fastq', or 'auto'")

    try:
        check_no_wrapped(filepath, fmt)
    except AssertionError:
        return {}

    subseq_lookup = {}
    for name, s in dict_subseqs.items():
        if isinstance(s, str):
            subseq_lookup[s] = name + ":-)fwd"
            subseq_lookup[revc(s)] = name + ":-)fwdrc"
        elif isinstance(s, list) and len(s) == 2:
            subseq_lookup[s[0]] = name + ":-)fwd"
            subseq_lookup[revc(s[0])] = name + ":-)fwdrc"
            subseq_lookup[s[1]] = name + ":-)rev"
            subseq_lookup[revc(s[1])] = name + ":-)revrc"
        else:
            raise ValueError("Unable to read 'dict_subseqs'")

    try:
        fwd_seqs, rev_seqs = zip(*(dict_subseqs.values()))
    except ValueError:
        fwd_seqs = list(dict_subseqs.values())
        rev_seqs = []

    if read2:
        if rev_seqs:
            subseqs = "|".join(list(rev_seqs) + [revc(s) for s in fwd_seqs])
        else:
            subseqs = "|".join([revc(s) for s in fwd_seqs])
    else:
        if rev_seqs:
            subseqs = "|".join(list(fwd_seqs) + [revc(s) for s in rev_seqs])
        else:
            subseqs = "|".join(list(fwd_seqs))

    return _find_exact_matches(filepath, subseqs, subseq_lookup, fmt)


def _find_approximate_matches(idx, subseq, seq, max_mismatch):
    if max_mismatch < 1:
        max_l_dist = int(round(len(subseq) * max_mismatch))
    else:
        max_l_dist = int(round(max_mismatch))

    matches = [
        PrimerMatch(idx, m.start, m.end, m.dist)
        for m in find_near_matches(subseq, seq, max_l_dist=max_l_dist)
    ]

    return get_best_match(matches)


def find_approximate_matches(seq, dict_subseqs, max_mismatch, read2=False, which=None):
    if isinstance(max_mismatch, list):
        if len(max_mismatch) == 1:
            mm1, mm2 = max_mismatch[0], max_mismatch[0]
        else:
            mm1, mm2 = max_mismatch
    elif isinstance(max_mismatch, (int, float)):
        mm1, mm2 = max_mismatch, max_mismatch
    else:
        raise ValueError

    if isinstance(list(dict_subseqs.values())[0], str):
        fwd_only = True
        matches1 = []
    else:
        fwd_only = False
        matches1, matches2 = [], []

    for name, s in dict_subseqs.items():
        if fwd_only:
            if read2:
                match1 = _find_approximate_matches(name, revc(s), seq, mm1)
            else:
                match1 = _find_approximate_matches(name, s, seq, mm1)
            if match1:
                matches1.append(match1)
        else:
            if read2:
                if which != 1:
                    match1 = _find_approximate_matches(name, s[1], seq, mm1)
                if which != 0:
                    match2 = _find_approximate_matches(name, revc(s[0]), seq, mm2)
            else:
                if which != 1:
                    match1 = _find_approximate_matches(name, s[0], seq, mm1)
                if which != 0:
                    match2 = _find_approximate_matches(name, revc(s[1]), seq, mm2)
            if which != 1 and match1:
                matches1.append(match1)
            if which != 0 and match2:
                matches2.append(match2)

    if fwd_only:
        return get_best_match(matches1)
    else:
        if which == 0:
            return get_best_match(matches1)
        elif which == 1:
            return get_best_match(matches2)
        else:
            return [get_best_match(matches1), get_best_match(matches2)]


def get_best_match(x):
    x = [x for x in x if isinstance(x, PrimerMatch)]
    if x:
        return min(x, key=lambda y: (y.dist, -(y.end - y.start)))
    else:
        return


def main(args):
    print("Start demultiplexing ... ({})".format(time.ctime()))
    t1 = time.time()

    if isinstance(args, dict):
        kwargs = args
    else:
        kwargs = vars(args)

    if "subcommand" in kwargs:
        kwargs.pop("subcommand")

    quiet = kwargs["quiet"]
    if not quiet:
        print_args(kwargs)

    inpath = kwargs.pop("inpath")
    glob_pattern = kwargs.pop("glob_pattern")
    pat = re.compile("{}.*".format(glob_pattern.strip("*").replace("*", ".*")))
    outdir = Path(kwargs["outdir"])

    dmp = Demultiplex(**kwargs)
    if not quiet:
        print("{} primer pairs are loaded successfully".format(len(dmp.dict_primer)))

    if len(inpath) == 1 and Path(inpath[0]).is_dir():
        indir = Path(inpath[0])
        files = indir.glob("{}".format(glob_pattern))
        idx_list = natsorted(set([pat.sub("", f.name) for f in files]))

        if len(idx_list) == 0:
            errmsg = "No files match with the pattern "
            errmsg += "'{}' found in {}".format(glob_pattern, indir)
            raise RuntimeError(errmsg)

        results = np.array([], dtype=int)
        if quiet:
            for idx in tqdm(idx_list, total=len(idx_list)):
                infile = natsorted(indir.glob("{}*".format(idx)), key=lambda x: x.name)
                dmp.set_infile(infile)
                dmp.set_out_prefix(idx)
                res = dmp.run()
                results = np.append(results, res)
        else:
            for idx in idx_list:
                print("." * shutil.get_terminal_size()[0])
                print("Processing {} ({}): ".format(idx, time.ctime()))
                infile = natsorted(indir.glob("{}*".format(idx)), key=lambda x: x.name)
                dmp.set_infile(infile)
                dmp.set_out_prefix(idx)
                res = dmp.run()
                results = np.append(results, res)

        n_reads = len(results)
        n_assined = results.sum()
        excldir = [x.name for x in outdir.glob("unassigned_*")] + ["result_summary"]
        n_loci = len(
            [x for x in outdir.glob("[!.*]*") if x.is_dir() and x.name not in excldir]
        )
        msg = "\nOverall, {} out of {} reads ({:.2%}) have been assigned to {} loci"
        msg = msg.format(n_assined, n_reads, n_assined / n_reads, n_loci)
        print(msg)

    else:
        dmp.quiet = False
        dmp.set_infile(inpath)
        if not kwargs["out_prefix"]:
            dmp.set_out_prefix(Path(common_prefix(inpath)).stem)
        results = np.array(dmp.run(), dtype=int)

    t2 = time.time()
    elapsed_time = time.strftime("%H:%M:%S", time.gmtime(t2 - t1))
    print("\nFinished (elapsed time: {})".format(elapsed_time))


if __name__ == "__main__":
    import sys

    argv = ["demultiplex"]
    argv.extend(sys.argv[1:])
    main(get_args(argv))
