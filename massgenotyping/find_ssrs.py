import re
from itertools import product

from fuzzysearch.common import group_matches

from .base import RepData, revc

rep_seq_pat = re.compile("\\(([ACGT]+)\\)_{([0-9]+)}([ACTG]*)")


def rotate(string):
    for i in range(len(string)):
        yield string[i:] + string[:i]


def motif_set(motif_class):
    for i in rotate(motif_class):
        yield i
    for i in rotate(revc(motif_class)):
        yield i


def gen_motif_classes(min_motif_len=2, max_motif_len=6):
    motifs = []
    for mlen in range(min_motif_len, max_motif_len + 1):
        for mcls in product(list("ACGT"), repeat=mlen):
            mcls = "".join(mcls)
            if len(set(mcls)) > 1:
                if mcls not in motifs:
                    motifs += list(motif_set(mcls))
                    yield mcls


def get_motif_class(motif):
    for mcls in gen_motif_classes(len(motif), len(motif) + 1):
        for m in motif_set(mcls):
            if m == motif:
                return mcls
    else:
        raise ValueError(
            "Unable to find the class of input motif. "
            "Maybe it contains a character other than A, C, G, T?"
        )


def get_best_match_in_group_RepData(group):
    try:
        return max(group, key=lambda x: (x.n_reps, -x.start))
    except ValueError:
        return


def remove_gaps(seq):
    seqout = ""
    idx = []
    for i, char in enumerate(seq):
        if char != "-":
            seqout += char
            idx.append(i)
    return seqout, idx


def _find_ssrs(seq, motif, min_repeats=3, max_interrupt=0):
    seqng, idx = remove_gaps(seq.upper())
    mlen = len(motif)

    def find_in_index_range(string, substring, n, start_index, end_index=None):
        try:
            match = string.find(substring * n, start_index, end_index)
        except AttributeError:
            raise ValueError("seq must be a str like object")
        start = match
        while match >= 0:
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
    seq,
    min_repeats=3,
    motif=None,
    motif_class=None,
    min_motif_len=2,
    max_motif_len=4,
    max_interrupt=None,
    return_only_best=False,
    start=0,
    end=None,
    **kwargs
):
    """
    Find short sequence repeats in the given sequence string

    Parameters
    -----
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
        maximum length of repeat motif (default: 5)
    max_interrupt: int
        maximum length of interruption to allow (default: None)
    return_only_best: bool
        remove overlapping matches and return the longest match (default: False)
    start: int
        starting postion where ssr needs to be search
    end: int
        ending postion where ssr needs to be search
    """
    subseq = seq[start:end]
    matches = []
    if not motif:
        if not motif_class:
            motif_classes = gen_motif_classes(min_motif_len, max_motif_len)
        elif isinstance(motif_class, str):
            motif_classes = [motif_class]
        elif not isinstance(motif_class, list):
            raise TypeError("motif_class must be a str or list")

        for mcls in motif_classes:
            for m in motif_set(mcls):
                for s, e, n, rs in _find_ssrs(subseq, m, min_repeats, max_interrupt):
                    matches.append(RepData(s + start, e + start, n, rs, m, mcls))
    else:
        if isinstance(motif_class, list):
            motifs = motif
        if isinstance(motif, str):
            motifs = [motif]
        else:
            raise TypeError("motif must be a str or list")

        motif_classes = [get_motif_class(m) for m in motifs]
        for m, mcls in zip(motifs, motif_classes):
            for s, e, n, rs in _find_ssrs(subseq, m, min_repeats, max_interrupt):
                matches.append(RepData(s + start, e + start, n, rs, m, mcls))

    match_groups = group_matches(matches)
    best_matches = [get_best_match_in_group_RepData(g) for g in match_groups]

    if return_only_best:
        return get_best_match_in_group_RepData(best_matches)
    else:
        return sorted(best_matches, key=lambda m: m.start)


def unfold_rep_seq(rep_seq):
    unfolded = ""
    matches = rep_seq_pat.findall(rep_seq)
    for m in matches:
        unfolded += m[0] * int(m[1]) + m[2]
    return unfolded


if __name__ == "__main__":
    # test
    import timeit

    seq = "ATT--AATATTA---TTATT-ATTTTTATTGCGGGGGGGGATTATTATTATT"
    motif = "ATT"
    rep = find_ssrs(seq, motif=motif, max_interrupt=3)
    print([r.rep_seq for r in rep])
    print([unfold_rep_seq(r.rep_seq) for r in rep])

    t1 = timeit.timeit(
        "find_ssrs(seq, motif=motif, return_only_best=True)",
        number=1000,
        globals=globals(),
    )
    print(find_ssrs(seq, motif=motif, return_only_best=True))
    print(t1)

    t2 = timeit.timeit(
        "find_ssrs(seq, return_only_best=True)", number=1000, globals=globals()
    )
    print(find_ssrs(seq, return_only_best=True))
    print(t2)
