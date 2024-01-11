from __future__ import annotations

import os
import re
import shutil
import tempfile
import time
from itertools import chain
from pathlib import Path

import Levenshtein
import matplotlib.pyplot as plt
import numpy as np
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, _DistanceMatrix
from fuzzysearch.common import group_matches
from matplotlib import colors, gridspec, patches
from natsort import natsorted

from .argument_parser import ask_user, get_args, print_args
from .base import SeqData, count_uniq_seq, read_fastx
from .find_ssrs import find_ssrs, get_longest_RepData, remove_gaps, unfold_rep_seq
from .utils import make_blastdb, run_mafft
from .variant_filter import pairwise_dist_Levenstein


class AlleleCheck(object):
    def __init__(
        self,
        infile=None,
        outdir=None,
        force_no_visual_check=False,
        quiet=False,
        **kwargs,
    ):
        self.infile = None
        self.outdir = None
        if infile:
            self.set_infile(infile)
        if outdir:
            self.outdir = Path(outdir)
        self.force_no_visual_check = force_no_visual_check
        self.quiet = quiet
        self.kwargs = kwargs
        self.seqdat = []
        self.align = []
        self.motifs = []
        self.ssr_regions = []

    def set_infile(self, infile):
        self.infile = Path(infile)
        self.indir = self.infile.resolve().parents[1]
        self.subdir = self.infile.resolve().parents[0]
        if not self.outdir:
            self.outdir = self.indir

    def add_rep_data(self):
        # clear rep_data
        for s in self.seqdat:
            s.rep_data = []

        # add rep_data to seqdat
        reps = [
            characterise_ssrs(self.align, r, m)
            for r, m in zip(self.ssr_regions, self.motifs)
        ]
        reps_ = [list(x) for x in zip(*reps)]
        for i, r in enumerate(reps_):
            self.seqdat[i].rep_data = r

    def init_allele_id(self, prefix):
        for s in self.seqdat:
            if not s.id:
                s.id = prefix + "*"

    def add_allele_id(self):
        name_list = []
        allele_group = {}
        for s in self.seqdat:
            if s.id.split("*")[1]:
                name_list.append(s.id)
                group = s.id.split("*")[1].split(":")[0]
                if s.non_ssr_seq not in allele_group.keys():
                    allele_group[s.non_ssr_seq] = group

        non_ssr_uniq = sorted(
            set(map(lambda x: x.non_ssr_seq, self.seqdat)),
            key=lambda x: -len([s for s in self.seqdat if s.non_ssr_seq == x]),
        )

        if allele_group:
            for x in non_ssr_uniq:
                if x not in allele_group.keys():
                    i = 1
                    while "{:02}".format(i) in allele_group.values():
                        i += 1
                    allele_group[x] = "{:02}".format(i)
        else:
            for i, x in enumerate(non_ssr_uniq):
                allele_group[x] = "{:02}".format(i + 1)

        if name_list:
            for s in self.seqdat:
                if not s.id.split("*")[1]:
                    s.id += allele_group[s.non_ssr_seq] + ":"
                    if s.rep_data:
                        s.id += "".join(
                            ["{:02}:".format(r.n_reps) for r in s.rep_data if r]
                        )
                    i = 1
                    while s.id + "{:02}".format(i) in name_list:
                        i += 1
                    s.id += "{:02}".format(i)
        else:
            for s in self.seqdat:
                s.id += allele_group[s.non_ssr_seq] + ":"
                if s.rep_data:
                    s.id += "".join(
                        ["{:02}:".format(r.n_reps) for r in s.rep_data if r]
                    )

            name_uniq = set(map(lambda x: x.id, self.seqdat))
            group = [[y for y in self.seqdat if y.id == x] for x in name_uniq]
            for g in group:
                g = sorted(g, key=l_dist_from_simple_repeats)
                for i, al in enumerate(g):
                    al.id += "{:02}".format(i + 1)

    def visual_check(self):
        while True:
            vca = VisualAlleleCheck(suptitle=self.subdir.name, **self.__dict__)
            if self.seq_discarded:
                while self.seq_discarded:
                    vca.selected.append(self.seq_discarded.pop())
                vca.change_rect_color()
            vca.show()
            if len(vca.selected) > 0:
                ds = ", ".join([self.seqdat[i].id for i in sorted(vca.selected)])
                msg = "Would you like to discard {} sequence{} ({})?".format(
                    len(vca.selected), "s" if len(vca.selected) > 1 else "", ds
                )
                if ask_user(msg, default="y", quit=True):
                    self.seqdat = [
                        self.seqdat[i]
                        for i in range(len(self.seqdat))
                        if i not in vca.selected
                    ]
                    self.align = AlignIO.MultipleSeqAlignment(
                        [
                            self.align[i]
                            for i in range(len(self.align))
                            if i not in vca.selected
                        ]
                    )
                    # remove gaps
                    self.align = remove_gap_pos(self.align)

                    self.ssr_regions, self.motifs = find_variable_ssrs(
                        self.align, **self.kwargs
                    )

                    # add rep_data to seqdat
                    self.add_rep_data()

                    if len(self.align) > 1:
                        self.tree = construct_tree(
                            self.align, self.ssr_regions, self.motifs
                        )
                    else:
                        self.tree = None
                    print("Reconstructing phylogeny ...")
            else:
                msg = "Keep all sequences and write results?"
                if ask_user(msg, default="y", quit=True):
                    break

    def write_results(self):
        # write the sequence alignment of all alleles to FASTA file
        with self.outpath0.open("w") as outfile:
            for s, a in zip(self.seqdat, self.align):
                outfile.write(">{}\n{}\n".format(s.id, a.seq))

        # write seqdat to json file
        json_str = ",\n".join(
            [s.to_json(indent=4, ensure_ascii=False) for s in self.seqdat]
        )
        json_str = "[\n{}\n]".format(json_str)
        with self.outpath1.open("w") as outfile:
            outfile.write(json_str)

    def make_database(self):
        """
        make a blast database
        """
        blastdb = str(self.subdir.joinpath(self.subdir.name + "_allele_blastdb"))
        ftmp, tmppath = tempfile.mkstemp()
        try:
            with os.fdopen(ftmp, "w") as tmpfile:
                for s in self.seqdat:
                    tmpfile.write(">{}\n{}\n".format(s.id, s.seq))
            make_blastdb(tmppath, self.subdir.name, blastdb)
        finally:
            os.remove(tmppath)

    def run(self):
        if not self.infile:
            raise RuntimeError("infile is not set")

        self.outpath0 = self.subdir.joinpath(self.subdir.name + "_allele_aligned.fasta")
        self.outpath1 = self.subdir.joinpath(self.subdir.name + "_allele_data.json")

        seq_count = count_uniq_seq(self.infile, read_count_in_id=True)
        self.seqdat = [SeqData("", s, c, samples=n) for s, (c, n) in seq_count.items()]

        self.seq_discarded = []
        if self.outpath0.exists():
            msg = "Output file {} already exists. Would you like to merge?"
            msg = msg.format(self.outpath0.name)
            try:
                if ask_user(msg, default="n", overwrite=True):
                    seq_kept, seq_kept_id = zip(
                        *[
                            [str(rec.seq).replace("-", "").upper(), rec.id]
                            for rec in read_fastx(self.outpath0)
                        ]
                    )
                else:
                    print("Skip {}".format(self.infile.name))
                    return
            except Exception:
                pass
            else:
                msg = "Would you like to keep previous allele names?"
                keep_name = ask_user(msg, default="y", quit=False)
                for i, s in enumerate(self.seqdat):
                    if s.seq in seq_kept:
                        if keep_name:
                            s.id = seq_kept_id[seq_kept.index(s.seq)]
                    else:
                        self.seq_discarded.append(i)

                self.align = run_mafft(
                    [s.seq for s in self.seqdat],
                    add_to_existence_alignment=str(self.outpath0),
                )
                self.align = self.align[-len(self.seqdat) :]

        if not self.align:
            self.align = run_mafft(
                [s.seq for s in self.seqdat],
                reorder=False,
                opts="--globalpair --maxiterate 1000",
            )

        # find variable SSRs
        self.ssr_regions, self.motifs = find_variable_ssrs(self.align, **self.kwargs)

        # add rep_data to seqdat
        self.add_rep_data()

        # add non-SSR sequences
        for s, a in zip(self.seqdat, self.align):
            s.non_ssr_seq = get_non_ssr_seq(str(a.seq), self.ssr_regions)

        # add allele names
        self.init_allele_id(prefix=self.subdir.name)
        self.add_allele_id()

        # construct an allele tree
        if len(self.align) > 1:
            self.tree = construct_tree(self.align, self.ssr_regions, self.motifs)
        else:
            self.tree = None

        # add lists of potential stutter sequences
        all_seqs = [s.seq for s in self.seqdat]
        for s in self.seqdat:
            le, mo = potential_stutter_sequences(s)
            s.stutter_s = [
                self.seqdat[all_seqs.index(x)].id
                for x in set(le).intersection(set(all_seqs))
            ]
            s.stutter_l = [
                self.seqdat[all_seqs.index(x)].id
                for x in set(mo).intersection(set(all_seqs))
            ]

        if not self.force_no_visual_check:
            self.visual_check()
        self.write_results()
        self.make_database()


class VisualAlleleCheck(object):
    """
    Draw an interactive figure from the input sequence alignment
    """

    def __init__(self, seqdat, align, tree, ssr_regions, motifs, suptitle="", **kwargs):
        if tree:
            self.order = [
                int(c.name.split("_")[-1].replace("seq", ""))
                for c in tree.get_terminals()
            ]
        else:
            self.order = list(range(len(seqdat)))
        self.text_list = []
        for o in self.order:
            s = seqdat[o]
            n_sample = len(s.samples)
            if n_sample > 20:
                samples = ", ".join(s.samples[:20])
                samples += ", ..."
            else:
                samples = ", ".join(s.samples)
            if n_sample == 1:
                samples += " ({} sample)".format(n_sample)
            else:
                samples += " ({} samples)".format(n_sample)

            text = r"$\bf{" + s.id + "}$"
            text += "\n{} bp; {} reads in total; ".format(s.length, s.counts)
            text += "{}\n".format(
                "/".join(
                    [re.sub("(_{[0-9]*})", "$\\1$", r.rep_seq) for r in s.rep_data if r]
                )
            )
            text += samples
            self.text_list.append(text)
        self.draw_figure(seqdat, align, tree, ssr_regions, motifs, suptitle)
        self.selected = []
        self.idx_tmp = -1
        self.show_info(self.idx_tmp)
        self.subplot = None

    def draw_figure(self, seqdat, align, tree, ssr_regions=[], motifs=[], suptitle=""):
        self.fig = plt.figure(figsize=(16, 7.5))

        gs0 = gridspec.GridSpec(1, 1)
        gs0.update(top=0.94, bottom=0.88)
        gs1 = gridspec.GridSpec(1, 6)
        gs1.update(top=0.85, bottom=0.06, left=0.04, right=0.84, wspace=1)
        gs2 = gridspec.GridSpec(1, 1)
        gs2.update(top=0.85, bottom=0.06, left=0.86, right=0.96)

        self.ax0 = self.fig.add_subplot(gs0[:, :])
        self.ax1 = self.fig.add_subplot(gs1[:, :1])
        self.ax2 = self.fig.add_subplot(gs1[:, 1:])
        self.ax3 = self.fig.add_subplot(gs2[:, :])

        self.ax0.set_axis_off()

        # phylogeny
        if tree:
            Phylo.draw(tree, do_show=False, axes=self.ax1, label_func=lambda x: None)
            self.ax1.set_xlim((0, max(list(tree.depths().values()))))

        # sequence alignment
        align_arr = np.array(align)[self.order, :]
        if np.any(align_arr == "-"):
            charlist = ["a", "c", "g", "t", "-"]
            choicelist = [0, 1, 2, 3, 4]
            cpal = ["#9b59b6", "#3498db", "#95a5a6", "#34495e", "#ffffff"]
        else:
            charlist = ["a", "c", "g", "t"]
            choicelist = [0, 1, 2, 3]
            cpal = ["#9b59b6", "#3498db", "#95a5a6", "#34495e"]

        condlist = [align_arr == x for x in charlist]
        cm = colors.ListedColormap(cpal)

        m, n = align_arr.shape
        x = np.arange(0, n + 1)
        y = np.arange(0, m + 1) + 0.5
        X, Y = np.meshgrid(x, y)
        Z = np.select(condlist, choicelist).astype(int)

        self.ax2.pcolor(X, Y, Z, cmap=cm)
        self.ax2.hlines(
            y=np.arange(0, m) + 0.5, xmin=0, xmax=n + 1, color=(1, 1, 1, 0.2), lw=2
        )

        # variable sites
        n_uniq = np.array([len(set(list(align[:, j]))) for j in range(n)])
        c = ["#e74c3c" if x > 1 else "#34495e" for x in n_uniq]
        s = [6 if x > 1 else 0.25 for x in n_uniq]
        self.ax2.scatter(np.arange(0, n, 1), [m + 1] * n, s=s, c=c)

        self.ax2.legend(
            [patches.Patch(color=cpal[i]) for i in range(4)],
            ["A", "C", "G", "T"],
            loc="upper right",
            bbox_to_anchor=(1, -0.04),
            ncol=4,
            borderaxespad=0.0,
            frameon=False,
            handlelength=1.2,
            handletextpad=0.6,
            columnspacing=0.6,
        )

        self.ax2.set_xlim(-1, n + 1)
        self.ax2.set_xticks(np.arange(0, n, 20))
        self.ax2.set_yticks(np.arange(1, m + 3))
        self.ytl2 = [re.sub("^.*\\*", "*", seqdat[o].id) for o in self.order] + [
            "variable site",
            "SSR region",
        ]
        self.ax2.set_yticklabels(self.ytl2)

        # ssr region
        if ssr_regions and motifs:
            pad = n * 0.01
            height = 0.8
            for sr, mot in zip(ssr_regions, motifs):
                paths = polygon_path(sr, m + 2, pad, height)
                polygon = plt.Polygon(paths, fc="#2ecc71", alpha=0.8)
                self.ax2.add_patch(polygon)
                self.ax2.text(np.sum(sr) / 2, m + 2, s=mot, va="center", ha="center")

        opts = dict(
            width=n,
            height=1,
            color=(1, 1, 1, 0),
            edgecolor=(1, 1, 1, 0),
            linewidth=2,
            zorder=5,
            picker=1,
        )

        self.rects = [
            x.get_children()[0]
            for x in [self.ax2.barh(i, **opts) for i in range(1, m + 1)]
        ]

        # barchert of the number of samples
        n_samples = [len(seqdat[o].samples) for o in self.order]
        self.ax3.barh(range(1, m + 1), width=n_samples, align="center")

        self.ax1.get_xaxis().tick_bottom()
        self.ax2.get_xaxis().tick_bottom()
        self.ax3.get_xaxis().tick_bottom()
        self.ax1.axes.get_yaxis().set_visible(False)
        self.ax2.yaxis.set_ticks_position("none")
        self.ax3.axes.get_yaxis().set_visible(False)
        for s in ["top", "right", "left"]:
            self.ax1.spines[s].set_visible(False)
            self.ax2.spines[s].set_visible(False)
            self.ax3.spines[s].set_visible(False)

        self.ax2.set_ylim((0.2, m + 2.8))
        self.ax1.set_ylim(self.ax2.get_ylim())
        self.ax3.set_ylim(self.ax2.get_ylim())

        self.ax1.set_xlabel("Branch length")
        self.ax2.set_xlabel("Nucleotide position (bp)")
        self.ax3.set_xlabel("Number of samples")

        self.fig.suptitle("{}".format(suptitle))

    def on_pick(self, event):
        artist = event.artist
        mouseevent = event.mouseevent

        if self.subplot == 2:
            idx = self.rects.index(artist)

            if mouseevent.button == 1:
                self.reset_info()
                if self.idx_tmp != idx:
                    self.idx_tmp = idx
                else:
                    self.idx_tmp = -1
                self.show_info(self.idx_tmp)

            elif mouseevent.button == 3:
                if self.order[idx] in self.selected:
                    self.selected.remove(self.order[idx])
                else:
                    self.selected.append(self.order[idx])
                self.change_rect_color()
            event.canvas.draw()

    def on_key(self, event):
        if event.key in ["up", "ctrl+k"]:
            if self.idx_tmp == len(self.rects) - 1:
                self.idx_tmp = -1
            else:
                self.idx_tmp += 1

        elif event.key in ["down", "ctrl+j"]:
            if self.idx_tmp == -1:
                self.idx_tmp = len(self.rects) - 1
            else:
                self.idx_tmp -= 1

        elif event.key in ["right", "ctrl+l"]:
            if self.idx_tmp > -1:
                if self.order[self.idx_tmp] not in self.selected:
                    self.selected.append(self.order[self.idx_tmp])

                self.change_rect_color()

        elif event.key in ["left", "ctrl+h"]:
            if self.idx_tmp > -1:
                if self.order[self.idx_tmp] in self.selected:
                    self.selected.remove(self.order[self.idx_tmp])

                self.change_rect_color()

        if self.idx_tmp > -2:
            self.reset_info()
            self.show_info(self.idx_tmp)

        event.canvas.draw()

    def change_rect_color(self):
        selected_idx = [self.order.index(i) for i in self.selected]
        for i, rect in enumerate(self.rects):
            if i in selected_idx:
                rect.set_edgecolor("mediumvioletred")
                rect.set_zorder(10)
            else:
                rect.set_edgecolor((1, 1, 1, 0))
                rect.set_zorder(5)

    def show_info(self, idx):
        if idx == -1:
            s = "Left click on a sequence  (or press \u2191/\u2193 key) to "
            s += "show information \n right click (or press \u2190/\u2192"
            s += r" key) to select the sequence to $\bf{discard}$"
            self.txt0 = self.ax0.text(
                0.5, 0.6, s=s, color="grey", va="center", ha="center", wrap=True
            )
            return

        rect = self.rects[idx]
        rect.set_facecolor((0, 0, 0, 0.2))
        self.txt0 = self.ax0.text(
            0.5,
            0.6,
            s=self.text_list[idx],
            color="k",
            va="center",
            ha="center",
            wrap=True,
        )
        ytl2 = [
            r"$\bf{" + ytl + "}$" if i == idx else ytl
            for i, ytl in enumerate(self.ytl2)
        ]
        self.ax2.set_yticklabels(ytl2)

    def reset_info(self):
        if self.txt0:
            self.txt0.remove()
            self.txt0 = ""
        for rect in self.rects:
            rect.set_facecolor((1, 1, 1, 0))
        self.ax2.set_yticklabels(self.ytl2)

    def enter_axes(self, event):
        x0, y0 = event.inaxes.get_position().min
        if (x0, y0) == tuple(self.ax0.get_position().min):
            self.subplot = 0
        elif (x0, y0) == tuple(self.ax1.get_position().min):
            self.subplot = 1
        elif (x0, y0) == tuple(self.ax2.get_position().min):
            self.subplot = 2
        elif (x0, y0) == tuple(self.ax3.get_position().min):
            self.subplot = 3

    def leave_axes(self, event):
        self.subplot = None

    def connect(self):
        self.cidpick = self.fig.canvas.mpl_connect("pick_event", self.on_pick)
        self.cidkey = self.fig.canvas.mpl_connect("key_press_event", self.on_key)
        self.cidenteraxes = self.fig.canvas.mpl_connect(
            "axes_enter_event", self.enter_axes
        )
        self.cidleaveaxes = self.fig.canvas.mpl_connect(
            "axes_leave_event", self.leave_axes
        )

    def disconnect(self):
        """disconnect all the stored connection ids"""
        self.fig.canvas.mpl_disconnect(self.cidpick)
        self.fig.canvas.mpl_disconnect(self.cidkey)
        self.fig.canvas.mpl_disconnect(self.cidenteraxes)
        self.fig.canvas.mpl_disconnect(self.cidleaveaxes)

    def show(self):
        self.connect()
        plt.show()
        self.disconnect()


def polygon_path(ssr_range, y, pad=2, height=0.8):
    x0, x1 = ssr_range[0] + pad, ssr_range[1] - pad
    y0, y1 = y - height / 2, y + height / 2
    paths = [
        (x0, y0),
        (x0 - pad, (y0 + y1) / 2),
        (x0, y1),
        (x1, y1),
        (x1 + pad, (y0 + y1) / 2),
        (x1, y0),
    ]
    return paths


def l_dist_from_simple_repeats(allele):
    return sum(
        [
            Levenshtein.distance(*[unfold_rep_seq(r.rep_seq), r.motif * r.n_reps])
            for r in allele.rep_data
            if r
        ]
    )


def construct_tree(align, ssr_regions, motifs, weights=[1, 0.1]):
    """
    Construct an upgma tree based on a pairwise Levenshtein distance matrix.
    For each pairwise comparison, the Levenshtein distances are calculated
    for sequences of non-SSR and SSR regions separately, and the weighted
    sum of them are used as the distance to construct an upgma tree. By default,
    weights for non-SSR and SSR regions are 1 and 0.1, respectively. In SSR
    regions, one repeat difference is considered to be one edit distance.

    Parameters
    ----------
    align: Bio.AlignIO.MultipleSeqAlignment
        input sequence alignment
    ssr_regions: list of tuple
        start and end positions of SSR regions in the alignment
    motifs: list
        repeat motifs
    weights: list
        weights for non-SSR and SSR regions to culculate pairwise distances
        (default: [1, 0.1])
    """
    non_ssr_seqs = []
    ssr_seqs = []
    for a in align:
        seq = str(a.seq.upper())
        ssr_idx = np.array(list(chain(*[list(range(*x)) for x in ssr_regions])))
        non_ssr_idx = list(set(range(len(seq))) - set(ssr_idx))
        seq_arr = np.array(list(seq))
        non_ssr_seq = "".join(seq_arr[non_ssr_idx])
        non_ssr_seqs.append(non_ssr_seq)

        ssr_seq = ""
        for rr, mot in zip(ssr_regions, motifs):
            ssr_seq += seq[rr[0] : rr[1]].replace("-", "").replace(mot, "x")
        ssr_seqs.append(ssr_seq)

    mat1 = pairwise_dist_Levenstein(non_ssr_seqs)
    mat2 = pairwise_dist_Levenstein(ssr_seqs)

    mat = [
        list(np.array(i) * weights[0] + np.array(j) * weights[1])
        for i, j in zip(mat1, mat2)
    ]
    names = ["seq{}".format(i) for i in range(len(align))]
    dmat = _DistanceMatrix(names, mat)

    constructor = DistanceTreeConstructor()
    return constructor.upgma(dmat)


def remove_gap_pos(align):
    align_new = align[:, :0]
    for n in range(align.get_alignment_length()):
        char_uniq = set(align[:, n])
        if len(char_uniq) == 1 and "-" in char_uniq:
            pass
        else:
            align_new += align[:, n : (n + 1)]
    return align_new


def potential_stutter_sequences(allele):
    seq = allele.seq

    def more_or_less_repeat(rep_seq):
        pat_digit = re.compile("(\\d+)(?!\\d)")
        one_less_rep = []
        one_more_rep = []
        n_matches = len(pat_digit.findall(rep_seq))
        for i in range(1, n_matches + 1):
            j = -1 if i < 2 else i - 1
            ss = pat_digit.sub(lambda x: str(int(x.group(0)) - 1), rep_seq, i)
            ss = pat_digit.sub(lambda x: str(int(x.group(0)) + 1), ss, j)
            one_less_rep.append(ss)
            ss = pat_digit.sub(lambda x: str(int(x.group(0)) + 1), rep_seq, i)
            ss = pat_digit.sub(lambda x: str(int(x.group(0)) - 1), ss, j)
            one_more_rep.append(ss)
        return one_less_rep, one_more_rep

    seq_one_less_rep = []
    seq_one_more_rep = []
    for r in allele.rep_data:
        one_less_rep, one_more_rep = more_or_less_repeat(r.rep_seq)
        seq_one_less_rep.extend(
            [seq[: r.start] + unfold_rep_seq(x) + seq[r.end :] for x in one_less_rep]
        )
        seq_one_more_rep.extend(
            [seq[: r.start] + unfold_rep_seq(x) + seq[r.end :] for x in one_more_rep]
        )
    return seq_one_less_rep, seq_one_more_rep


def find_variable_ssrs(align, min_variants=3, **kwargs):
    """
    Find variable SSR regions from a multiple sequence alignment

    Parameters
    ----------
    align: Bio.AlignIO.MultipleSeqAlignment
        input alignment
    min_variants: int
        Minimum number of variants for valiable SSR regions
    **kwargs
        The keyward arguments are used for find_ssr()

    See Also
    ----------
    find_ssrs()
    """
    matches = [find_ssrs(str(a.seq), max_interrupt=0, **kwargs) for a in align]
    match_groups = group_matches(chain(*matches))
    ssr_regions = []
    motifs = []
    for group in match_groups:
        if len(group) >= min_variants:
            starts, ends = list(zip(*[[rep.start, rep.end] for rep in group]))
            ssr_regions.append((np.min(starts), np.max(ends)))
            motifs.append(get_longest_RepData(group).motif)

    return ssr_regions, motifs


def characterise_ssrs(align, ssr_region, motif, max_interrupt=None):
    start, end = ssr_region
    if not max_interrupt:
        max_interrupt = len(motif)

    reps = []
    for a in align:
        seq = str(a.seq)
        seq_ng, idx_ng = remove_gaps(seq)
        ssr_region_ng = [i for i, j in enumerate(idx_ng) if j >= start and j <= end]
        if ssr_region_ng:
            rep = find_ssrs(
                seq_ng,
                motif=motif,
                min_repeats=2,
                start=ssr_region_ng[0],
                end=ssr_region_ng[-1],
                max_interrupt=len(motif),
            )
            best = get_longest_RepData(rep)
            reps.append(best)
        else:
            reps.append(None)
    return reps


def get_non_ssr_seq(seq, ssr_regions=[]):
    seqout = ""
    init = 0
    for start, end in ssr_regions:
        seqout += seq[init:start]
        init = end
    seqout += seq[init:]
    return seqout.upper().replace("-", "")


def main(args):
    print("Start allele checking ... ({})".format(time.ctime()))
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

    indir = Path(kwargs.pop("indir"))
    glob_pattern = kwargs.pop("glob_pattern")
    files = natsorted(indir.glob("*/" + glob_pattern), key=lambda x: x.name)

    subdir = kwargs.pop("subdir")
    if subdir:
        files = [f for f in files if f.parent.name in subdir]

    for f in files:
        if not quiet:
            print("." * shutil.get_terminal_size()[0])
            print("Processing {} ({}): ".format(f.name, time.ctime()))
        ach = AlleleCheck(infile=f, **kwargs)
        ach.run()

    t2 = time.time()
    elapsed_time = time.strftime("%H:%M:%S", time.gmtime(t2 - t1))
    print("\nFinished (elapsed time: {})".format(elapsed_time))


if __name__ == "__main__":
    import sys

    argv = ["allele-check"]
    argv.extend(sys.argv[1:])
    main(get_args(argv))
