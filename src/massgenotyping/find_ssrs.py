from __future__ import annotations

import re
from collections.abc import Iterator
from itertools import product

from fuzzysearch.common import group_matches

from .base import RepData, revc

rep_seq_pat = re.compile("\\(([ACGT]+)\\)_{([0-9]+)}([ACTG]*)")


def rotate(string: str) -> Iterator[str]:
    """Rotate the characters in the input string one by one."""
    for i in range(len(string)):
        yield string[i:] + string[:i]


def motif_set(motif_class: str) -> Iterator[str]:
    """Return all sequences of the same motif class."""
    for i in rotate(motif_class):
        yield i
    for i in rotate(revc(motif_class)):
        yield i


def gen_motif_classes(min_motif_len: int = 2, max_motif_len: int = 6) -> Iterator[str]:
    """Generate repeat-motif classes for a given range of sequence lengths."""
    motifs = []
    for mlen in range(min_motif_len, max_motif_len + 1):
        for mcls_ in product("ACGT", repeat=mlen):
            mcls = "".join(mcls_)
            if len(set(mcls)) == 1:
                continue
            if mcls not in motifs:
                motifs.extend(list(motif_set(mcls)))
                yield mcls


def get_motif_class(motif: str) -> str:
    """Return the class of the given motif."""
    for mcls in gen_motif_classes(len(motif), len(motif) + 1):
        for m in motif_set(mcls):
            if m == motif:
                return mcls
    else:
        raise ValueError(
            "Unable to find the class of the given motif. "
            "Maybe it contains a character other than ['A', 'C', 'G', 'T']?"
        )


def remove_gaps(seq: str) -> tuple[str, list[int]]:
    """Remove gaps in the given sequence."""
    seqout, idx = "", []
    for i, char in enumerate(seq):
        if char != "-":
            seqout += char
            idx.append(i)
    return seqout, idx


def _find_ssrs(
    seq: str, motif: str, min_repeats: int = 3, max_interrupt: int = 0
) -> Iterator[tuple[int, int, int, str]]:
    seqng, idx = remove_gaps(seq.upper())
    mlen = len(motif)

    def find_in_index_range(
        string: str,
        substring: str,
        n: int,
        start_index: int,
        end_index: int | None = None,
    ) -> tuple[int, int, int]:
        try:
            start = string.find(substring * n, start_index, end_index)
        except AttributeError:
            raise ValueError("'seq' must be a str object")
        match = start
        while match > -1:
            n += 1
            match = string.find(substring * n, start_index, start + len(substring) * n)
        n_reps = n - 1
        return start, start + mlen * n_reps - 1, n_reps

    start, end, n_reps = find_in_index_range(seqng, motif, min_repeats, 0)
    while start > -1:
        rep_seq = "({})_{{{}}}".format(motif, n_reps)

        if max_interrupt:
            # 3'-end side
            start_new, end_new, n_reps_new = find_in_index_range(
                seqng, motif, 1, end + 1, end + 1 + max_interrupt + mlen
            )
            while start_new > -1:
                interrupt_seq = seqng[(end + 1) : start_new]
                rep_seq += "{}({})_{{{}}}".format(interrupt_seq, motif, n_reps_new)
                n_reps += n_reps_new + len(interrupt_seq) // mlen
                end = end_new
                start_new, end_new, n_reps_new = find_in_index_range(
                    seqng, motif, 1, end + 1, end + 1 + max_interrupt + mlen
                )

            # 5'-end side
            start_rev = len(seqng) - start
            start_new, end_new, n_reps_new = find_in_index_range(
                seqng[::-1],
                motif[::-1],
                1,
                start_rev,
                start_rev + max_interrupt + mlen,
            )
            while start_new > -1:
                interrupt_seq = seqng[::-1][start_rev:start_new][::-1]
                rep_seq_new = "({})_{{{}}}{}".format(motif, n_reps_new, interrupt_seq)
                rep_seq = rep_seq_new + rep_seq
                n_reps += n_reps_new + len(interrupt_seq) // mlen
                start = len(seqng) - (end_new + 1)
                start_rev = len(seqng) - start
                start_new, end_new, n_reps_new = find_in_index_range(
                    seqng[::-1],
                    motif[::-1],
                    1,
                    start_rev,
                    start_rev + max_interrupt + mlen,
                )

        yield idx[start], idx[end] + 1, n_reps, rep_seq
        start, end, n_reps = find_in_index_range(seqng, motif, min_repeats, end + 1)


def find_ssrs(
    seq: str,
    min_repeats: int = 3,
    motif: str | list[str] = [],
    motif_class: str | list[str] = [],
    min_motif_len: int = 2,
    max_motif_len: int = 6,
    max_interrupt: int = 0,
    start: int = 0,
    end: int | None = None,
    **kwargs,
) -> list[RepData] | None:
    """
    Find short sequence repeats in the given sequence string.

    Parameters
    ----------
    seq: str
        input sequence string
    min_repeats: int
        minimum number of repeats to search (default: 3)
    motif: str or list
        repeat motif to search (default: None)
    motif_class: str or list
        class of repeat motif to search (default: None)
    min_motif_len: int
        minimum length of repeat motif (default: 2)
    max_motif_len: int
        maximum length of repeat motif (default: 6)
    max_interrupt: int
        maximum length of interruption to allow (default: None)
    start: int
        starting postion where ssr needs to be search
    end: int
        ending postion where ssr needs to be search

    """
    subseq = seq[start:end]
    matches = []
    if motif:
        if isinstance(motif, list):
            motifs = motif
        elif isinstance(motif, str):
            motifs = [motif]
        else:
            raise TypeError("motif must be a str or list")

        motif_classes = [get_motif_class(m) for m in motifs]
        for m, mcls in zip(motifs, motif_classes):
            for s, e, n, rs in _find_ssrs(subseq, m, min_repeats, max_interrupt):
                matches.append(RepData(s + start, e + start, n, rs, m, mcls))
    else:
        if not motif_class:
            motif_classes = list(gen_motif_classes(min_motif_len, max_motif_len))
        elif isinstance(motif_class, list):
            motif_classes = motif_class
        elif isinstance(motif_class, str):
            motif_classes = [motif_class]
        else:
            raise TypeError("motif_class must be a str or list")

        for mcls in motif_classes:
            for m in motif_set(mcls):
                for s, e, n, rs in _find_ssrs(subseq, m, min_repeats, max_interrupt):
                    matches.append(RepData(s + start, e + start, n, rs, m, mcls))

    match_groups = group_matches(matches)
    best_matches = [get_longest_RepData(g) for g in match_groups if g]

    return list(sorted(best_matches, key=lambda m: m.start))


def get_longest_RepData(group: list[RepData]) -> RepData:
    """Return the best match RepData in the group."""
    return max(group, key=lambda x: (x.n_reps, -x.start))


def unfold_rep_seq(rep_seq: str) -> str:
    """Unfold the rep_seq string."""
    unfolded = ""
    matches = rep_seq_pat.findall(rep_seq)
    for m in matches:
        unfolded += m[0] * int(m[1]) + m[2]
    return unfolded
