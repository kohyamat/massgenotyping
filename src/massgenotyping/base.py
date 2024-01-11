from __future__ import annotations

import gzip
import shutil
import subprocess
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from Bio import SeqIO
from dataclasses_json import config, dataclass_json


@dataclass_json
@dataclass(frozen=True)
class PrimerMatch:
    match: str = field(metadata=config(field_name="match"))
    start: int = field(metadata=config(field_name="start"))
    end: int = field(metadata=config(field_name="end"))
    dist: int = field(metadata=config(field_name="dist"))


@dataclass_json
@dataclass(frozen=True)
class RepData:
    start: int = field(metadata=config(field_name="start"))
    end: int = field(metadata=config(field_name="end"))
    n_reps: int = field(metadata=config(field_name="n_reps"))
    rep_seq: str = field(metadata=config(field_name="rep_seq"))
    motif: str = field(metadata=config(field_name="motif"))
    motif_class: str = field(metadata=config(field_name="motif_class"))


@dataclass_json
@dataclass
class SeqData:
    id: str = field(metadata=config(field_name="id"))
    seq: str = field(metadata=config(field_name="seq"))
    counts: int = field(metadata=config(field_name="counts"))
    length: int = field(default=0)
    rep_data: list[RepData] = field(
        default_factory=list, metadata=config(field_name="rep_data")
    )
    samples: list[str] = field(
        default_factory=list, metadata=config(field_name="samples")
    )
    non_ssr_seq: str = field(default="", metadata=config(field_name="non_ssr_seq"))
    stutter_s: list[str] = field(
        default_factory=list, metadata=config(field_name="stutter_s")
    )
    stutter_l: list[str] = field(
        default_factory=list, metadata=config(field_name="stutter_l")
    )

    def __post_init__(self):
        self.seq = self.seq.upper()
        self.length = len(self.seq.replace("-", ""))


@dataclass_json
@dataclass
class GenData:
    id: str = field(metadata=config(field_name="id"))
    n_reads_all: int = field(metadata=config(field_name="n_reads_all"))
    genotype: list[str] = field(
        default_factory=list, metadata=config(field_name="genotype")
    )
    n_reads_each: dict[str, int] = field(
        default_factory=dict, metadata=config(field_name="n_reads_each")
    )

    def __post_init__(self):
        try:
            err = "Allele sets differed between genotype and n_reads_each"
            alleles = [x for x in self.genotype if x != "NA"]
            assert set(alleles) == set(self.n_reads_each.keys()), err
        except AssertionError:
            print(set(alleles), set(self.n_reads_each))


class MarkerData(object):
    def __init__(
        self,
        path_marker_data=None,
        path_fwd_primers=None,
        path_rev_primers=None,
        verbose=False,
        **kwargs,
    ):
        self.dict_primer = self.read_marker_file(
            path_marker_data, path_fwd_primers, path_fwd_primers, verbose
        )
        self.path_marker_data = path_marker_data
        if self.path_marker_data:
            self.set_frag_len()
            self.set_max_alleles()

    def set_frag_len(self):
        try:
            self.dict_frag_len = gen_dict_from_table(
                self.path_marker_data, "Name", "Frag_len"
            )
        except ValueError:
            self.dict_frag_len = None

    def set_max_alleles(self):
        try:
            self.dict_max_alleles = gen_dict_from_table(
                self.path_marker_data, "Name", "Max_alleles"
            )
        except ValueError:
            self.dict_max_alleles = None

    def read_marker_file(
        self,
        path_marker_data=None,
        path_fwd_primers=None,
        path_rev_primers=None,
        verbose=False,
    ):
        if path_marker_data:
            dict_primer = gen_dict_from_table(
                path_marker_data, key="Name", value=["Fwd", "Rev"]
            )
        elif path_fwd_primers and path_rev_primers:
            len_f = count_records(path_fwd_primers)
            len_r = count_records(path_rev_primers)
            if len_f != len_r:
                msg = "The number of sequences differs between "
                msg += "{} and {}".format(path_fwd_primers, path_rev_primers)
                raise RuntimeError(msg)

            dict_primer = {}
            msg = "\nThe names of foward and rev primers do not match. "
            msg += "The foward primer names are used as marker names:"
            for f, r in zip(read_fastx(path_fwd_primers), read_fastx(path_rev_primers)):
                if f.id != r.id and verbose:
                    if msg:
                        print(msg)
                        msg = ""
                    print("Fwd: {0}, Rev: {1} -----> {0}".format(f.id, r.id))
                dict_primer[f.id] = [str(f.seq), str(r.seq)]
        elif path_fwd_primers:
            msg = "'path_fwd_primers' not provided"
            raise ValueError(msg)
        elif path_rev_primers:
            msg = "'path_rev_primers' not provided"
            raise ValueError(msg)
        else:
            msg = "No information for primer sequences provided"
            raise ValueError(msg)

        return dict_primer


def revc(s: str) -> str:
    """Return the reverse compelement of given nucleotide sequence."""
    o = "ACGTUWSMKRYBDHVNZacgtuwsmkrybdhvnz-"
    c = "TGCAAWSKMYRVHDBNZtgcaawskmyrvhdbnz-"
    if len(set(s) & set(o)) > len(set(s)):
        errmsg = "invalid character was found in the sequeces"
        raise RuntimeError(errmsg)
    return s.translate(str.maketrans(o, c))[::-1]


def check_file(filepath: str | Path):
    """Check if a path exists and if it is a file."""
    if isinstance(filepath, str):
        filepath = Path(filepath)

    if not filepath.exists():
        errmsg = "File not found: {}".format(filepath)
        raise FileNotFoundError(errmsg)

    if not filepath.is_file():
        errmsg = "'{}' is not a file".format(filepath)
        raise RuntimeError(errmsg)


def check_no_wrapped(filepath: str | Path, fmt: str = "fastq"):
    """
    Check the input sequence file and raise an error if it is wrapped.

    Parameters
    ----------
    filepath :
        path to the sequence file
    fmt :
        file format ('fasta' or 'fastq')

    """
    if isinstance(filepath, str):
        filepath = Path(filepath)
    check_file(filepath)

    if filepath.suffix == ".gz":
        if shutil.which("rg"):
            prog0 = "rg -z"
            prog1 = "rg"
        else:
            prog0 = "zgrep"
            prog1 = "grep"
    else:
        if shutil.which("rg"):
            prog0 = prog1 = "rg"
        else:
            prog0 = prog1 = "grep"

    if fmt == "fasta":
        cmd0 = r"{} -n -m 20 ^> {}".format(prog0, filepath).split()
        cmd1 = r"cut -d: -f1".split()
        res = subprocess.Popen(cmd0, stdout=subprocess.PIPE)
        res = subprocess.Popen(cmd1, stdin=res.stdout, stdout=subprocess.PIPE)
        j, k = 1, 2

    elif fmt == "fastq":
        cmd0 = r"{} -A2 ^@ {}".format(prog0, filepath).split()
        cmd1 = r"{} -n -m 20 ^\+".format(prog1).split()
        cmd2 = r"cut -d: -f1".split()
        res = subprocess.Popen(cmd0, stdout=subprocess.PIPE)
        res = subprocess.Popen(cmd1, stdin=res.stdout, stdout=subprocess.PIPE)
        res = subprocess.Popen(cmd2, stdin=res.stdout, stdout=subprocess.PIPE)
        j, k = 3, 4

    line_nos = res.stdout.read().decode("utf-8").strip().split()
    if line_nos:
        line_nos = np.array(line_nos).astype(int)
        n = len(line_nos)
        assert np.all((line_nos - j) / k == np.arange(n))


def count_records(filepath: str | Path, fmt: str = "fastq", opts: str = "") -> int:
    """
    Count the number of sequence records in a fasta/fastq file.

    Parameters
    ----------
    filepath:
        path to the input fasta/fastq file
    fmt:
        file format (default: "fasta")
    opts:
        options for grep

    """
    if isinstance(filepath, str):
        filepath = Path(filepath)
    check_file(filepath)

    if filepath.suffix == ".gz":
        if shutil.which("rg"):
            prog0 = "rg -z"
            prog1 = "rg"
        else:
            prog0 = "zgrep"
            prog1 = "grep"
    else:
        if shutil.which("rg"):
            prog0 = "rg"
            prog1 = "rg"
        else:
            prog0 = "grep"
            prog1 = "grep"

    if fmt == "fasta":
        cmd0 = "{} -c ^> {} {}".format(prog0, filepath, opts).split()
        res0 = subprocess.Popen(cmd0, stdout=subprocess.PIPE)
        return int(res0.stdout.read().decode("utf-8").strip())

    elif fmt == "fastq":
        cmd0 = r"{} -A2 ^@ {} {}".format(prog0, filepath, opts).split()
        cmd1 = r"{} -c ^\+".format(prog1).split()
        res0 = subprocess.Popen(cmd0, stdout=subprocess.PIPE)
        res1 = subprocess.Popen(cmd1, stdin=res0.stdout, stdout=subprocess.PIPE)
        line_no = res1.stdout.read().decode("utf-8").strip()
        if line_no:
            return int(line_no)
        else:
            return 0

    else:
        errmsg = "fmt must be either 'fastq' or 'fasta'"
        raise ValueError(errmsg)


def read_fastx(filepath: str | Path, fmt: str = "auto") -> Iterator[SeqIO.SeqRecord]:
    """
    Read a fasta/fastq file and return a generator of Bio.Seq.SeqRecord objects.

    Parameters
    ----------
    filepath:
        path to the input fasta/fastq file
    fmt:
        'fasta', 'fastq' or 'auto' (default: "auto")

    See Also
    --------
    Bio.SeqIO.parse

    """
    if isinstance(filepath, str):
        filepath = Path(filepath)

    if fmt == "auto":
        fmt = guess_fmt(filepath)
    elif fmt not in ["fasta", "fastq"]:
        raise ValueError("'fmt' must be 'fasta', 'fastq', or 'auto'")

    if not count_records(filepath, fmt):
        errmsg = "No sequence records found in {}".format(filepath)
        raise RuntimeError(errmsg)

    if filepath.suffix == ".gz":
        with gzip.open(filepath, "rt") as handle:
            for record in SeqIO.parse(handle, fmt):
                yield record
    else:
        with filepath.open("r") as handle:
            for record in SeqIO.parse(handle, fmt):
                yield record


def count_uniq_seq(filepath: str | Path, read_count_in_id: bool = False, **kwargs):
    """
    Count the number of reads for each unique sequence.

    Parameters
    ----------
    filepath :
        path to the input fasta/fastq file
    fmt : []
        "fasta", "fastq" or "auto" (default: "auto")
    read_count_in_id:
        read counts in sequence id
    kwargs:
        keyward arguments

    """
    seq_count = {}
    for rec in read_fastx(filepath, **kwargs):
        seq = str(rec.seq).replace("-", "")
        if read_count_in_id:
            count = int(rec.id.split(":")[-1])
            idx = "_".join(rec.id.split("_")[:-1])
            if seq in seq_count.keys():
                seq_count[seq][0] += count
                seq_count[seq][1] += [idx]
            else:
                seq_count[seq] = [count, [idx]]
        else:
            if seq in seq_count.keys():
                seq_count[seq] += 1
            else:
                seq_count[seq] = 1

    if read_count_in_id:
        seq_count = dict(
            sorted(
                seq_count.items(),
                key=lambda x: (len(x[1][1]), x[1][0]),
                reverse=True,
            )
        )
    else:
        seq_count = dict(sorted(seq_count.items(), key=lambda x: x[1], reverse=True))

    return seq_count


def gen_dict_from_table(filepath, key, value, header=True, delimiter=","):
    """
    Generate a dict object from a tablular file.

    Parameters
    ----------
    filepath : str or Path
        path to the input tabular file
    key : int or str
        column index or name for dict key
    value : int or str or list
        column index or name for dict value
    header : bool
        Whether the input file contains a header row (default: True)
    delimiter : str
        delimiter character of the input file (default: ",")

    """
    filepath = Path(filepath)
    check_file(filepath)

    if not delimiter:
        if filepath.suffix == ".csv":
            delimiter = ","
        if filepath.suffix == ".tsv":
            delimiter = r"\t+"
        if filepath.suffix == ".txt":
            delimiter = None

    with filepath.open() as f:
        line = f.readline()
        if header:
            hd = line.strip().split(delimiter)
            line = f.readline()
        else:
            hd = []
        items_list = []
        while line:
            items_list.append([parse_item(x) for x in line.strip().split(delimiter)])
            line = f.readline().strip()

    if not hd and not isinstance(key, int):
        raise TypeError("key must be int when no header row in input file")
    else:
        if isinstance(key, int):
            idx0 = key
        else:
            try:
                idx0 = hd.index(key)
            except ValueError:
                raise ValueError("Column '{}' not found in {}".format(key, filepath))

    if isinstance(value, str) or isinstance(value, int):
        value = [value]
    elif not isinstance(value, list):
        raise TypeError("value must be a str, int, or list")

    if np.array(value).dtype == "int64":
        idx1 = value
    else:
        try:
            idx1 = [hd.index(v) for v in value]
        except ValueError:
            not_found = ", ".join([v for v in value if v not in hd])
            raise ValueError("Column '{}' not found in {}".format(not_found, filepath))

    if len(idx1) > 1:
        return {items[idx0]: [items[x] for x in idx1] for items in items_list}
    else:
        return {items[idx0]: items[idx1[0]] for items in items_list}


def guess_fmt(filepath):
    """
    Guess the file format (FASTA/FASTQ) of the input sequence file
    """
    filepath = Path(filepath)
    if filepath.name.find("fastq") > -1:
        fmt = "fastq"
    elif filepath.name.find("fasta") > -1:
        fmt = "fasta"
    elif count_records(filepath, "fastq", "-m 10"):
        fmt = "fastq"
    elif count_records(filepath, "fasta", "-m 10"):
        fmt = "fasta"
    else:
        raise RuntimeError("Unable to determine file format")

    return fmt


def parse_item(s):
    return int(s) if s.isdigit() else float(s) if isfloat(s) else s


def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False
