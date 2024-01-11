from __future__ import annotations

from collections.abc import Iterator

import numpy as np
from Bio import SeqIO
from fuzzysearch import find_near_matches
from fuzzysearch.common import get_best_match_in_group

from .base import revc


def sliding_window(
    x: list[int], window_size: int = 10, step_size: int = 1
) -> Iterator[list[int]]:
    """
    Generate substrings by the sliding window algorithm.

    Parameters
    ----------
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
    seq_record: SeqIO.SeqRecord,
    threshold: int = 20,
    window_size: int = 10,
    step_size: int = 1,
    return_n_trim: bool = False,
) -> SeqIO.SeqRecord | tuple[SeqIO.SeqRecord, int]:
    """
    Trim low quality bases in the 3'- side of the sequence based on the phred quality
    scores.

    Parameters
    ----------
    seq_record : Bio.Seq.SeqRecord
        input sequence record with quality scores
    threshold : int
        cut off threshold of the quality score (default: 20)
    window_size : int
        window size (default: 10)
    step_size : int
        step size (default: 1)
    return_n_trim : bool
        return number of bases trimmed (default: False)

    """
    qual = seq_record.letter_annotations["phred_quality"]
    for i, q in enumerate(sliding_window(qual[::-1], window_size, step_size)):
        if np.mean(q) > threshold:
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


def trim_primer(
    seq_record: SeqIO.SeqRecord,
    primer_seqs: list[str],
    max_mismatch: float | int = 0.14,
) -> SeqIO.SeqRecord:
    """
    Trim primer sequences.

    Parameters
    ----------
    seq_record : Bio.Seq.SeqRecord
        input sequence record
    primer_seqs : list
        list of the foward and reverse primer sequnces
    max_mismatch : float
        Maximum number (or proportion) of mismatches allowed for searching primer
        sequeces (default: 0.14)

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
        raise ValueError("max_mismatch must be a positive value")

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
