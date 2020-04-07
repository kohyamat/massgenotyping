import numpy as np
from Bio import SeqIO
from fuzzysearch import find_near_matches
from fuzzysearch.common import get_best_match_in_group

from .base import revc


def sliding_window(x, window_size=10, step_size=1):
    """
    Generate substrings by the sliding window algorithm

    Parameters
    -----
    x : str
        input string
    window_size : int
        window size (default: 10)
    step_size : int
        step size (default: 1)
    """
    n_chunks = ((len(x) - window_size) / step_size) + 1
    for i in range(0, int(n_chunks * step_size), step_size):
        yield x[i : i + window_size]


def trim_low_qual(
    seq_record, threshold=20, window_size=10, step_size=1, return_n_trim=False
):
    """
    Trim low quality bases in the 3'- side of the sequence based on the phred
    quality scores.

    Parameters
    -----
    seq_record : Bio.Seq.SeqRecord
        input sequence record with quality score
    threshold : int
        cut off threshold of the quality score (default: 20)
    window_size : int
        window size (default: 10)
    step_size : int
        step size (default: 1)
    """
    qual = np.array(seq_record.letter_annotations["phred_quality"])
    for i, q in enumerate(sliding_window(qual[::-1], window_size, step_size)):
        if q.mean() > threshold:
            break

    n_trim = i * step_size
    if n_trim:
        tr = seq_record[:-n_trim]
    else:
        tr = seq_record[:]

    if return_n_trim:
        return tr, n_trim
    else:
        return tr


def trim_primer(seq_record, primer_seqs, max_mismatch=0.14):
    """
    Trim primer sequences.

    Parameters
    -----
    seq_record : Bio.Seq.SeqRecord
        input sequence record
    primer_seqs : list
        list of the foward and reverse primer sequnces
    max_mismatch : float
        Maximum proportion of mismatches allowed for serching primer sequeces
        (range: 0-1; default: 0.14)
    """
    seq = seq_record.seq
    fwd, rev = primer_seqs
    rev_rc = revc(rev)
    len_fwd, len_rev = len(fwd), len(rev)

    if max_mismatch > 1:
        max_l_dist1 = max_l_dist2 = max_mismatch
    elif max_mismatch > 0:
        max_l_dist1 = round(len_fwd * max_mismatch)
        max_l_dist2 = round(len_rev * max_mismatch)
    else:
        print("max_mismatch must be a positive value")
        return

    m0 = find_near_matches(fwd, str(seq), max_l_dist=max_l_dist1)
    m1 = find_near_matches(rev_rc, str(seq), max_l_dist=max_l_dist2)

    if len(m0) > 0:
        match_fwd = get_best_match_in_group(m0)
    if len(m1) > 0:
        match_rev_rc = get_best_match_in_group(m1)

    if len(m0) > 0 and len(m1) > 0:
        tr = seq_record[match_fwd.end : match_rev_rc.start]
    elif len(m0) > 0:
        tr = seq_record[match_fwd.end :]
    elif len(m1) > 0:
        tr = seq_record[: match_rev_rc.start]
    else:
        tr = seq_record[:]

    return tr


if __name__ == "__main__":
    # test
    import time
    from utils import read_fastx

    seq_file = "../example/sequence_data/test1_S1_L001_R1_001.fastq.gz"

    # quality trimming
    opts = {"threshold": 20, "window_size": 10, "step_size": 1}
    t0 = time.time()
    records_tr, n_tr = zip(
        *[trim_low_qual(x, **opts, return_n_trim=True) for x in read_fastx(seq_file)]
    )
    print(
        "{} bases in all were trimmed across {} sequences".format(sum(n_tr), len(n_tr))
    )
    print("time elapsed: {:.5f}".format(time.time() - t0))

    # primer trimming
    primer_seqs = ["ATGCATGC", "AGAGGGGC"]
    rec = SeqIO.SeqRecord(seq="ATGCATGCAAAAAAAAAAGCCCCTCT")
    print(trim_primer(rec, primer_seqs, 0.12).seq)
