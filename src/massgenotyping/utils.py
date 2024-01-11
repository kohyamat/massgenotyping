from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from io import StringIO
from pathlib import Path

from Bio import AlignIO, Phylo, SeqIO


def common_prefix(x):
    prefix = ""
    i = 0
    while all([xx[i] == x[0][i] for xx in x]):
        prefix += x[0][i]
        i += 1
        if (i + 1) > min([len(xx) for xx in x]):
            break
    return prefix


def run_NGmerge(
    read1,
    read2,
    outdir=".",
    out_prefix="",
    min_overlap=20,
    max_mismatch=0.1,
    min_overlap_dovetailed=50,
):
    """
    Run NGmerge by using the subprocess module.

    Parameters
    -----
    read1 : str
        Input FASTQ file with reads from forward direction
    read2 : str
        Input FASTQ file with reads from reverse direction
    outdir : str
        Output directory (default: ".")
    out_prefix : str
        Prefix of output file
    min_overlap : int
        Minimum overlap of the paired-end reads (default 20)
    max_mismatch : float
        Mismatches to allow in the overlapped region (a fraction of the overlap length;
        default 0.10)
    min_overlap_dovetailed:
        Minimum overlap of dovetailed alignments (default 50)

    Notes
    -----
    [1] Gasper, "NGmerge: merging paired-end reads via novel empirically-derived models
    of sequencing errorsPartTree," BMC Bioinformatics, vol. 19, 536, 2018.
    """
    outdir = Path(outdir)

    if not shutil.which("NGmerge"):
        msg = "The exutable file 'NGmerge' not found in $PATH"
        raise RuntimeError(msg)

    if not out_prefix:
        out_prefix = Path(common_prefix([read1, read2])).name
        out_prefix = out_prefix.rstrip("_R")

    cmd = "NGmerge -1 {} -2 {} ".format(str(read1), str(read2))
    cmd += " -o {}".format(outdir.joinpath(out_prefix + "_merged.fastq.gz"))
    cmd += " -m {} -p {}".format(min_overlap, max_mismatch)
    cmd += " -d -e {}".format(min_overlap_dovetailed)
    cmd += " -f {} -z".format(outdir.joinpath(out_prefix + "_notMerged"))

    proc = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # stdout = proc.stdout.read().decode("utf-8")
    stderr = proc.stderr.read().decode("utf-8")

    if stderr:
        raise RuntimeError(stderr)


def run_flash(
    read1,
    read2,
    outdir=".",
    out_prefix="out",
    max_mismatch_density=0.2,
    frag_len=None,
    read_len=None,
    max_overlap=None,
    return_max_overlap=False,
    verbose=False,
):
    """
    Run FLUSH by using the subprocess module.

    Parameters
    -----
    read1 : str
        Input FASTQ file with reads from forward direction
    read2 : str
        Input FASTQ file with reads from reverse direction
    outdir : str
        Output directory (default: ".")
    out_prefix : str
        Prefix of output file

    Notes
    -----
    [1] Magoč & Salzberg, "FLASH: fast length adjustment of short reads to improve
    genome assemblies," Bioinformatics, vol. 27, pp. 2957–2963, 2011.
    """

    if not shutil.which("flash"):
        msg = "The exutable file 'flash' not found in $PATH"
        raise RuntimeError(msg)

    cmd0 = "flash {} {} -d {} -o {} -O -z -t 1 -x {}".format(
        str(read1), str(read2), outdir, out_prefix, max_mismatch_density
    )

    def increment_max_overlap(cmd, max_overlap_start=50):
        max_overlap = max_overlap_start
        while True:
            cmd = cmd + " -M %i" % round(max_overlap)
            proc = subprocess.Popen(
                cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stderr = proc.stderr.read().decode("utf-8")
            stdout = [x.decode("utf-8") for x in list(proc.stdout.readlines())]

            if stderr.find("--max-overlap") < 0:
                break
            max_overlap += 10

        return stdout, stderr, max_overlap

    if max_overlap:
        stdout, stderr, max_overlap = increment_max_overlap(cmd0, max_overlap)
    elif frag_len and read_len:
        cmd1 = cmd0 + " -r {:.0f} -f {:.0f} -s {:.0f}".format(
            read_len, frag_len, frag_len * 0.1
        )
        proc = subprocess.Popen(
            cmd1.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        stderr = proc.stderr.read().decode("utf-8")
        stdout = proc.stdout.read().decode("utf-8")

        if stderr.find("--max-overlap") > 0:
            stdout, stderr, max_overlap = increment_max_overlap(cmd0)
        else:
            max_overlap = None
    else:
        stdout, stderr, max_overlap = increment_max_overlap(cmd0)

    if stderr:
        raise RuntimeError(stderr)
    elif verbose:
        [print(line) for line in stdout]

    if return_max_overlap:
        return max_overlap


def run_mafft(
    obj, treeout=False, reorder=False, add_to_existence_alignment="", opts=""
):
    """
    Run MAFFT by using the subprocess module.
    Output is a Bio.Align.MultipleSeqAlignment object.

    Parameters
    -----
    obj : Bio.Align.MultipleSeqAlignment, list of Bio.Seq.SeqRecord or str
        input sequence alingment
    treeout : bool
        return guide tree (default: False)
    reorder : bool
        reorder output alingments (default: False)
    opts : str, optional
        other mafft options

    See Also
    --------
    Bio.Align.MultipleSeqAlignment
    Bio.Seq.SeqRecord

    Notes
    -----
    [1] Katoh & Toh, "PartTree: an algorithm to build an approximate tree from
    a large number of unaligned sequences (describes the PartTree algorithm),"
    Bioinformatics, vol. 23, pp. 372-374, 2007.

    [2] Katoh, Kuma, Toh & Miyata, "MAFFT version 5: improvement in accuracy
    of multiple sequence alignment (describes [ancestral versions of] the
    G-INS-i, L-INS-i and E-INS-i strategies)," Nucleic Acids Res.,
    vol. 33, pp. 511-518, 2005.

    [3] Katoh, Misawa, Kuma & Miyata, "MAFFT: a novel method for rapid
    multiple sequence alignment based on fast Fourier transform (describes
    the FFT-NS-1, FFT-NS-2 and FFT-NS-i strategies)," Nucleic Acids Res.,
    vol. 30, pp. 3059-3066, 2002.
    """
    if not shutil.which("mafft"):
        msg = "The exutable file 'mafft' not found in $PATH"
        raise RuntimeError(msg)

    ftmp, tmppath = tempfile.mkstemp()
    try:
        if isinstance(obj, AlignIO.MultipleSeqAlignment):
            with os.fdopen(ftmp, "w") as tmp:
                AlignIO.write(obj, tmp, "fasta")
        elif isinstance(obj, list) and isinstance(obj[0], SeqIO.SeqRecord):
            with os.fdopen(ftmp, "w") as tmp:
                SeqIO.write(obj, tmp, "fasta")
        elif isinstance(obj, list) and isinstance(obj[0], str):
            with os.fdopen(ftmp, "w") as tmp:
                for i, seq in enumerate(obj):
                    tmp.write(">seq{}\n{}\n".format(i, seq))
        else:
            os.close(ftmp)
            os.remove(tmppath)
            return
        if treeout:
            opts += " --treeout"
        if reorder:
            opts += " --reorder"

        if add_to_existence_alignment:
            if os.path.exists(add_to_existence_alignment):
                cmd = "mafft {} --add {} {}".format(
                    opts, tmppath, add_to_existence_alignment
                )
            else:
                "{} not exists".format(add_to_existence_alignment)
        else:
            cmd = "mafft {} {}".format(opts, tmppath)
        proc = subprocess.Popen(
            cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        stdout = proc.stdout.read().decode("utf-8")
        stderr = proc.stderr.read().decode("utf-8")

        if treeout:
            treefile = tmppath + ".tree"
            tree = Phylo.read(treefile, "newick")

    finally:
        os.remove(tmppath)
        if treeout:
            os.remove(treefile)

    if stdout:
        align = AlignIO.read(StringIO(stdout), "fasta")
        if treeout:
            return align, tree
        else:
            return align
    else:
        print(stderr)
        return


def make_blastdb(seq_file, title, out, dtype="nucl"):
    if not shutil.which("makeblastdb"):
        msg = "The exutable file 'makeblastdb' not found in $PATH"
        raise RuntimeError(msg)
    cmd = "makeblastdb -in {} -title {} ".format(seq_file, title)
    cmd += "-out {} -dbtype {}".format(out, dtype)
    proc = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stderr = proc.stderr.read().decode("utf-8")
    if stderr:
        raise RuntimeError(stderr)


def run_blastn(
    query,
    database,
    out,
    outfmt=10,
    max_target_seqs=1,
    gapopen=1,
    gapextend=1,
    penalty=-4,
    evalue="1.0E-10",
):
    if not shutil.which("blastn"):
        msg = "The exutable file 'blast' not found in $PATH"
        raise RuntimeError(msg)
    cmd = "blastn -query {} -db {} -out {} ".format(query, database, out)
    cmd += "-outfmt {} -max_target_seqs {} ".format(outfmt, max_target_seqs)
    cmd += "-evalue {}".format(evalue)
    # cmd += '-gapopen {} -gapextend {} '.format(gapopen, gapextend)
    # cmd += '-penalty {} -evalue {}'.format(penalty, evalue)
    proc = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stderr = proc.stderr.read().decode("utf-8")
    if stderr:
        r = "Warning: [blastn] Examining 5 or more matches is recommended"
        stderr = stderr.replace(r, "").strip()
        # if stderr:
        # raise RuntimeError(stderr)
