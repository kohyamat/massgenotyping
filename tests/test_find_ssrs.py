from massgenotyping.base import RepData
from massgenotyping.find_ssrs import (
    find_ssrs,
    gen_motif_classes,
    get_longest_RepData,
    get_motif_class,
    motif_set,
    remove_gaps,
    rotate,
    unfold_rep_seq,
)


def test_rotate():
    ans = ["test", "estt", "stte", "ttes"]
    for x, y in zip(rotate("test"), ans):
        assert x == y


def test_motif_set():
    ans = ["AAG", "AGA", "GAA", "CTT", "TTC", "TCT"]
    for x, y in zip(motif_set("AAG"), ans):
        assert x == y


def test_gen_motif_classes():
    ans = ["AAC", "AAG", "AAT", "ACC", "ACG", "ACT", "AGC", "AGG", "ATC", "CCG"]
    for x, y in zip(gen_motif_classes(3, 3), ans):
        assert x == y


def test_get_motif_class():
    for test_string in ["AAG", "AGA", "GAA", "CTT", "TTC", "TCT"]:
        assert get_motif_class(test_string) == "AAG"


def test_remove_gaps():
    result = remove_gaps("ATT--AAT")
    assert result == ("ATTAAT", [0, 1, 2, 5, 6, 7])


class TestFindSSRs:
    seq = "GC-GC---GCTTAT-TATTATTACCTTACCCTTATT"

    def test_find_ssrs_default(self):
        result = find_ssrs(self.seq)
        expect = [
            RepData(0, 10, 3, rep_seq="(GC)_{3}", motif="GC", motif_class="CG"),
            RepData(10, 23, 4, rep_seq="(TTA)_{4}", motif="TTA", motif_class="AAT"),
        ]
        assert result == expect

    def test_find_ssrs_with_motif1(self):
        result = find_ssrs(self.seq, motif="GC")
        expect = [RepData(0, 10, 3, rep_seq="(GC)_{3}", motif="GC", motif_class="CG")]
        assert result == expect

    def test_find_ssrs_with_motif2(self):
        result = find_ssrs(self.seq, motif="TTA")
        expect = [
            RepData(10, 23, 4, rep_seq="(TTA)_{4}", motif="TTA", motif_class="AAT")
        ]
        assert result == expect

    def test_find_ssrs_with_motif3(self):
        result = find_ssrs(self.seq, motif="TAT")
        expect = [
            RepData(11, 21, 3, rep_seq="(TAT)_{3}", motif="TAT", motif_class="AAT")
        ]
        assert result == expect

    def test_find_ssrs_with_motif4(self):
        result = find_ssrs(self.seq, motif=["GC", "TTA"])
        expect = [
            RepData(0, 10, 3, rep_seq="(GC)_{3}", motif="GC", motif_class="CG"),
            RepData(10, 23, 4, rep_seq="(TTA)_{4}", motif="TTA", motif_class="AAT"),
        ]
        assert result == expect

    def test_find_ssrs_with_motif5(self):
        result = find_ssrs(self.seq, motif="AT")
        assert result == []

    def test_find_ssrs_with_motif_class1(self):
        result = find_ssrs(self.seq, motif_class="AAT")
        expect = [
            RepData(10, 23, 4, rep_seq="(TTA)_{4}", motif="TTA", motif_class="AAT")
        ]
        assert result == expect

    def test_find_ssrs_with_motif_class2(self):
        result = find_ssrs(self.seq, motif_class=["CG", "AAT"])
        expect = [
            RepData(0, 10, 3, rep_seq="(GC)_{3}", motif="GC", motif_class="CG"),
            RepData(10, 23, 4, rep_seq="(TTA)_{4}", motif="TTA", motif_class="AAT"),
        ]
        assert result == expect

    def test_find_ssrs_with_min_repeats(self):
        result = find_ssrs(self.seq, min_repeats=4)
        expect = [
            RepData(10, 23, 4, rep_seq="(TTA)_{4}", motif="TTA", motif_class="AAT")
        ]
        assert result == expect

    def test_find_ssrs_with_max_interrupt1(self):
        result = find_ssrs(self.seq, motif="TTA", max_interrupt=2)
        expect = [
            RepData(
                10,
                28,
                5,
                rep_seq="(TTA)_{4}CC(TTA)_{1}",
                motif="TTA",
                motif_class="AAT",
            )
        ]
        assert result == expect

    def test_find_ssrs_with_max_interrupt2(self):
        result = find_ssrs(self.seq, motif="TTA", max_interrupt=3)
        expect = [
            RepData(
                10,
                34,
                7,
                rep_seq="(TTA)_{4}CC(TTA)_{1}CCC(TTA)_{1}",
                motif="TTA",
                motif_class="AAT",
            )
        ]
        assert result == expect

    def test_find_ssrs_with_max_motif_len(self):
        result = find_ssrs(self.seq, max_motif_len=2)
        expect = [RepData(0, 10, 3, rep_seq="(GC)_{3}", motif="GC", motif_class="CG")]
        assert result == expect

    def test_find_ssrs_with_min_motif_len(self):
        result = find_ssrs(self.seq, min_motif_len=3)
        expect = [
            RepData(10, 23, 4, rep_seq="(TTA)_{4}", motif="TTA", motif_class="AAT")
        ]
        assert result == expect

    def test_find_ssrs_with_start(self):
        result = find_ssrs(self.seq, start=10)
        expect = [
            RepData(10, 23, 4, rep_seq="(TTA)_{4}", motif="TTA", motif_class="AAT")
        ]
        assert result == expect

    def test_find_ssrs_with_end(self):
        result = find_ssrs(self.seq, end=10)
        expect = [RepData(0, 10, 3, rep_seq="(GC)_{3}", motif="GC", motif_class="CG")]
        assert result == expect


class TestUnfoldRepSeq:
    def test_unfold_rep_seq(self):
        result = unfold_rep_seq("(TTA)_{10}")
        expect = "TTATTATTATTATTATTATTATTATTATTA"
        assert result == expect

    def test_unfold_rep_seq_with_interruption(self):
        result = unfold_rep_seq("(AT)_{3}AC(AT)_{5}")
        expect = "ATATATACATATATATAT"
        assert result == expect

    def test_unfold_rep_seq_with_two_interruption(self):
        result = unfold_rep_seq("(AT)_{3}AC(AT)_{5}CA(AT)_{2}")
        expect = "ATATATACATATATATATCAATAT"
        assert result == expect


def test_get_longest_RepData():
    group = [RepData(0, 12, 4, "", "", ""), RepData(0, 15, 5, "", "", "")]
    result = get_longest_RepData(group)
    assert result == group[1]
