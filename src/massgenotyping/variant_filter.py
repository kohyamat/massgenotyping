from __future__ import annotations

import gzip
import os
import re
import shutil
import signal
import sys
import tempfile
import time
from itertools import chain
from multiprocessing import Manager, Pool
from pathlib import Path

import Levenshtein
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib import colors, get_backend, gridspec, patches, ticker
from natsort import natsorted
from tqdm import tqdm

from .argument_parser import ask_user, get_args, print_args
from .base import MarkerData, SeqData, count_records, count_uniq_seq, read_fastx
from .find_ssrs import find_ssrs, get_longest_RepData
from .utils import run_mafft


class VariantFilter(MarkerData):
    def __init__(
        self,
        indir,
        outdir=None,
        subdir=[],
        glob_pattern="*_denoised.fasta.gz",
        threshold=0.25,
        max_alleles=2,
        min_l_dist=3,
        min_reads=5,
        force_visual_check=False,
        force_no_visual_check=False,
        n_cpu=1,
        quiet=False,
        **kwargs,
    ):
        self.indir = Path(indir)
        self.subdirs = natsorted(
            set(
                [
                    f.parent
                    for f in self.indir.glob("**/" + glob_pattern)
                    if f.parent.name.find("unassigned") < 0
                ]
            ),
            key=lambda x: x.name,
        )
        if subdir:
            self.subdirs = [sd for sd in self.subdirs if sd.name in subdir]
        self.outdir = Path(outdir) if outdir else Path(indir)
        self.glob_pattern = glob_pattern
        self.threshold = threshold
        self.min_l_dist = min_l_dist
        self.min_reads = min_reads
        self.force_visual_check = force_visual_check
        self.force_no_visual_check = force_no_visual_check
        self.n_cpu = n_cpu
        self.quiet = quiet

        if "path_marker_data" in kwargs and kwargs["path_marker_data"]:
            super().__init__(**kwargs)
        else:
            self.dict_max_alleles = {}

        if not self.dict_max_alleles:
            self.dict_max_alleles = {sd.name: max_alleles for sd in self.subdirs}

    def filter(self, seq_file):
        marker = seq_file.parent.name
        idx = re.sub(
            "{}.*".format(self.glob_pattern.strip("*").replace("*", ".*")),
            "",
            seq_file.name,
        )
        idx = idx.replace("_" + marker, "")

        outpath1 = self.outdir.joinpath(
            "result_summary", "variant_filter", marker, idx + ".csv"
        )

        if not self.overwrite and outpath1.exists():
            return 0

        if count_records(seq_file, "fasta") < self.min_reads:
            return 0

        seqdat = [
            SeqData("Seq{}".format(i), s, c)
            for i, (s, c) in enumerate(count_uniq_seq(seq_file).items())
        ]

        n_reads = np.array([x.counts for x in seqdat])

        meet = n_reads > n_reads.max() * self.threshold
        seq_meet = [s.seq for s, m in zip(seqdat, meet) if m]
        n_alleles = self.dict_max_alleles[marker]

        if self.force_no_visual_check:
            keep = meet.astype(int)
        else:
            if meet.sum() <= n_alleles:
                keep = meet.astype(int)
            elif meet.sum() > n_alleles:
                if np.any([x not in self.uniq_seqs for x in seq_meet]):
                    keep = None
                else:
                    keep = meet.astype(int)
                    keep[n_alleles:] = 0
            else:
                l_dists = np.array(list(chain(*pairwise_dist_Levenstein(seq_meet))))
                if np.any(l_dists < self.min_l_dist) > 0:
                    keep = None
                else:
                    keep = meet.astype(int)

        if keep is not None:
            keep = meet.astype(int)
            self.uniq_seqs.update(seq_meet)
            self.write_table(seqdat, keep, outpath1)
            self.write_sequences_tmp(idx, seqdat, keep)
            return 0
        else:
            if self.n_cpu > 1:
                self.mp_queue.append([seq_file, seqdat, meet])
            else:
                self.queue.append((seq_file, seqdat, meet))
            return 1

    def visual_check(self, seq_file=None, seqdat=None, meet=None):
        keep = None
        idx = re.sub(
            "{}.*".format(self.glob_pattern.strip("*").replace("*", ".*")),
            "",
            seq_file.name,
        )
        marker = seq_file.parent.name
        outpath1 = self.outdir.joinpath(
            "result_summary", "variant_filter", marker, idx + ".csv"
        )

        if not self.overwrite and outpath1.exists():
            return

        n_alleles = self.dict_max_alleles[marker]
        if self.force_visual_check:
            if count_records(seq_file, "fasta") < self.min_reads:
                return

            seqdat = [
                SeqData("Seq{}".format(i), s, c)
                for i, (s, c) in enumerate(count_uniq_seq(seq_file).items())
            ]

        else:
            seq_meet = [s.seq for s, m in zip(seqdat, meet) if m]
            new = set(seq_meet).difference(self.uniq_seqs)
            if not new:
                if meet.sum() > n_alleles:
                    keep = meet.astype(int)
                    keep[n_alleles:] = 0
                else:
                    keep = meet.astype(int)

        in_list = [True if s.seq in self.uniq_seqs else False for s in seqdat]

        for s in seqdat:
            try:
                s.rep_data = get_longest_RepData(find_ssrs(s.seq))
            except ValueError:
                pass

        if keep is None:
            while True:
                vsc = VisualCheck(seqdat, suptitle=idx, in_list=in_list)
                vsc.show()
                if vsc.selected != []:
                    break
                elif ask_user("No sequence was selected. Proceed anyway?"):
                    break
            keep = [1 if i in vsc.selected else 0 for i in range(len(seqdat))]
            seq_keep = [s.seq for s, k in zip(seqdat, keep) if k == 1]
            self.uniq_seqs.update(seq_keep)
            vsc.disconnect()

        if np.sum(keep) > 0:
            self.write_table(seqdat, keep, outpath1)
            self.write_sequences_tmp(idx, seqdat, keep)

    def visual_check_force(self, seq_file):
        self.visual_check(seq_file=seq_file)

    def write_table(self, seqdat, keep, outpath):
        with outpath.open("w") as outfile:
            outfile.write("ID,Keep,Counts,Length,Seq\n")
            for i, (s, k) in enumerate(zip(seqdat, keep)):
                results = [
                    "seq{}".format(i),
                    str(k),
                    str(s.counts),
                    str(s.length),
                    s.seq,
                ]
                outfile.write(",".join(results) + "\n")

    def write_sequences_tmp(self, idx, seqdat, keep):
        pid = os.getpid()
        tmppath = self.tmpdir.joinpath("tmp{}.fasta".format(pid))
        with tmppath.open("a") as outfile:
            for i, (s, k) in enumerate(zip(seqdat, keep)):
                if k == 1:
                    outfile.write(">{}_seq{}:{}\n{}\n".format(idx, i, s.counts, s.seq))

    def process_all_files_in_subdir(self, subdir):
        files = natsorted(subdir.glob(self.glob_pattern), key=lambda x: x.name)
        outpath0 = self.outdir.joinpath(subdir.name, subdir.name + "_filtered.fasta.gz")
        outpath0.parent.mkdir(parents=True, exist_ok=True)

        resdir = self.outdir.joinpath("result_summary", "variant_filter", subdir.name)
        resdir.mkdir(parents=True, exist_ok=True)

        self.tmpdir = Path(tempfile.mkdtemp())

        if not self.quiet:
            print("." * shutil.get_terminal_size()[0])
            print("Processing {} ({}): ".format(subdir.name, time.ctime()))

        self.overwrite = False
        self.uniq_seqs = set()
        if outpath0.exists():
            msg = "Output file '{}' already exists. Overwrite?"
            if ask_user(msg.format(str(outpath0)), "n"):
                outpath0.unlink()
                self.overwrite = True
            else:
                self.uniq_seqs.update([str(x.seq) for x in read_fastx(outpath0)])

        if self.force_visual_check:
            for f in files:
                self.visual_check_force(f)
        else:
            self.queue = []
            if self.n_cpu == 1:
                if self.quiet:
                    res = [self.filter(f) for f in files]
                else:
                    res = [self.filter(f) for f in tqdm(files, total=len(files))]
            else:
                with Manager() as manager:
                    self.mp_queue = manager.list()
                    original_sigint_handler = signal.getsignal(signal.SIGINT)
                    pool = Pool(self.n_cpu)
                    signal.signal(signal.SIGINT, original_sigint_handler)
                    try:
                        if self.quiet:
                            res = list(pool.imap(self.filter, files))
                        else:
                            res = list(
                                tqdm(
                                    pool.imap(self.filter, files),
                                    total=len(files),
                                )
                            )
                    except KeyboardInterrupt:
                        print("Caught KeyboardInterrupt, terminating workers")
                        pool.terminate()
                        pool.join()
                        shutil.rmtree(str(self.tmpdir))
                        raise KeyboardInterrupt
                    else:
                        pool.close()
                        pool.join()
                        while self.mp_queue:
                            self.queue.append(self.mp_queue.pop())

            if np.sum(res) > 0:
                self.queue = self.queue[::-1]
                if not self.quiet:
                    print("\nStart visual-checking procedure ...")
                while self.queue:
                    self.visual_check(*self.queue.pop())

        for ftmp in self.tmpdir.glob("tmp*.fasta"):
            with gzip.open(outpath0, "wb") as outfile:
                with ftmp.open(mode="rb") as tmp:
                    shutil.copyfileobj(tmp, outfile, 1024 * 1024 * 10)
        shutil.rmtree(str(self.tmpdir))

        if not self.quiet and subdir != self.subdirs[-1] and np.sum(res) > 0:
            if not ask_user("Done. Proceed to the next?", default="y", quit=False):
                sys.exit()

    def run(self):
        for sd in self.subdirs:
            self.process_all_files_in_subdir(sd)


class VisualCheck(object):
    def __init__(self, seqdat, align=None, suptitle="", in_list=None):
        self.draw_figure(seqdat, align, suptitle)

        self.labels = []
        for s in seqdat:
            lab = r"{}; {} bp; {} reads".format(s.id, s.length, s.counts)
            if s.rep_data:
                if isinstance(s.rep_data, list):
                    rs = [re.sub("(_{[0-9]*})", "$\\1$", r.rep_seq) for r in s.rep_data]
                    lab += r"; {}".format(", ".join(rs))
                else:
                    lab += r"; {}".format(
                        re.sub("(_{[0-9]*})", "$\\1$", s.rep_data.rep_seq)
                    )
            self.labels.append(lab)

        if in_list:
            x = []
            y = []
            for i, rect in enumerate(self.rects1):
                if in_list[i]:
                    x0, y0 = rect.get_xy()
                    h, w = rect.get_height(), rect.get_width()
                    x.append(x0 + w * 0.5)
                    y.append(y0 + h * 0.5)
            scatter = self.ax1.plot(x, y, "o", zorder=100, color="#34495e")
            self.ax1.legend(
                scatter,
                ("already in the list",),
                loc="upper right",
                bbox_to_anchor=(1, 0.9),
                fontsize=10,
                frameon=False,
                handletextpad=0.1,
            )

        self.selected = []
        self.subplot = None
        self.txt0 = ""
        self.ann1 = ""
        self.idx_tmp = -1
        self.text_selected()

    def draw_figure(self, seqdat, align=None, suptitle=""):
        self.fig = plt.figure(figsize=(10, 8))
        gs0 = gridspec.GridSpec(1, 1)
        gs0.update(top=0.94, bottom=0.88)

        gs1 = gridspec.GridSpec(2, 1)
        gs1.update(top=0.85, bottom=0.06, hspace=0.3)

        self.ax0 = self.fig.add_subplot(gs0[:, :])
        self.ax1 = self.fig.add_subplot(gs1[0, :])
        self.ax2 = self.fig.add_subplot(gs1[1, :])

        self.ax0.set_axis_off()

        # Barplot of the distribution of sequence lengths
        cols = sns.cubehelix_palette(
            len(seqdat), start=2.4, rot=0.2, dark=0.4, reverse=True
        )

        lengths = np.array([x.length for x in seqdat])
        bottoms = {b: 0 for b in set(lengths)}

        for i, s in enumerate(seqdat):
            self.ax1.bar(
                s.length,
                s.counts,
                bottom=bottoms[s.length],
                color=cols[i],
                edgecolor="white",
                linewidth=1,
                picker=1,
            )
            bottoms[s.length] += s.counts

        self.ax1.set_title("Distribution of sequence lengths of unique variants")
        self.ax1.set_xlabel("Sequence length (bp)")
        self.ax1.set_ylabel("Number of reads")
        self.ax1.set_xlim(lengths.min() - 4, lengths.max() + 4)
        self.ax1.set_ylim(0, self.ax1.get_ylim()[1] * 1.2)
        self.ax1.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        self.ax1.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

        # sequence alignment
        if not align:
            align = run_mafft([x.seq for x in seqdat[:10]])
        align_arr = np.array(align)
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
        y = np.arange(0, m + 1) - 0.5
        X, Y = np.meshgrid(x, y)
        Z = np.select(condlist, choicelist).astype(int)

        self.ax2.pcolor(X, Y, Z, cmap=cm)
        self.ax2.hlines(
            y=np.arange(0, m - 1) + 0.5, xmin=0, xmax=n + 1, color=(1, 1, 1, 0.2), lw=3
        )

        # variable sites
        n_uniq = np.array([len(set(list(align[:, j]))) for j in range(n)])
        c = ["#e74c3c" if x > 1 else "#34495e" for x in n_uniq]
        s = [6 if x > 1 else 0.25 for x in n_uniq]
        self.ax2.scatter(np.arange(0, n, 1), [-1] * n, s=s, c=c)

        self.ax2.set_title(
            "Sequence alignment (sorted in decreasing order of frequencies)"
        )
        self.ax2.set_xlabel("Nucleotide position (bp)")

        self.ax2.set_xlim(-0.5, n + 0.5)
        self.ax2.set_xticks(np.arange(0, n, 20))
        self.ax2.set_ylim(-1.5, m)
        self.ax2.set_yticks(np.arange(-1, m))
        self.ytl2 = ["variable site"] + [re.sub(r"^.*\*", "*", s.id) for s in seqdat]
        self.ax2.set_yticklabels(self.ytl2)

        self.ax2.spines["right"].set_visible(False)
        self.ax2.spines["top"].set_visible(False)
        self.ax2.spines["left"].set_visible(False)
        self.ax2.invert_yaxis()

        self.ax2.legend(
            [patches.Patch(color=cpal[i]) for i in range(4)],
            ["A", "C", "G", "T"],
            loc="upper right",
            bbox_to_anchor=(1, -0.08),
            ncol=4,
            borderaxespad=0.0,
            frameon=False,
            handlelength=1.2,
            handletextpad=0.6,
            columnspacing=0.6,
        )

        self.rects1 = self.get_rects(self.ax1)

        kwargs = dict(
            width=self.ax2.get_xlim()[1] - 0.5,
            height=1,
            color=(1, 1, 1, 0),
            edgecolor=(1, 1, 1, 0),
            linewidth=3,
            zorder=1,
            picker=1,
        )
        self.rects2 = [
            x.get_children()[0] for x in [self.ax2.barh(i, **kwargs) for i in range(m)]
        ]

        xl0, xl1 = self.ax1.get_xlim()
        self.reset_button = self.ax1.text(
            xl0 + (xl1 - xl0) * 0.96,
            self.ax1.get_ylim()[1] * 0.96,
            s="Reset",
            ha="right",
            va="top",
            bbox=dict(boxstyle="round", ec=(1.0, 0.5, 0.5), fc=(1.0, 0.8, 0.8)),
            picker=1,
        )

        self.cols_org = []
        self.zorder_org = []
        for container in self.ax1.containers:
            for child in container.get_children():
                if isinstance(child, plt.Rectangle):
                    self.cols_org.append(child.get_fc())
                    self.zorder_org.append(child.get_zorder())

        self.fig.suptitle("{}".format(suptitle))

    def on_pick(self, event):
        artist = event.artist
        mouseevent = event.mouseevent

        if self.subplot is None:
            return

        if self.subplot == 1 and isinstance(artist, plt.Text):
            # click on the reset button
            self.reset_color_all()
            self.reset_info()
            self.selected = []
            self.text_selected()

        else:
            if self.subplot == 1:
                idx = self.rects1.index(artist)
            else:
                idx = self.rects2.index(artist)

            if mouseevent.button == 1:
                self.reset_info()
                if self.idx_tmp != idx:
                    self.idx_tmp = idx
                else:
                    self.idx_tmp = -1
                self.show_info(self.idx_tmp)

            elif mouseevent.button == 3:
                if idx not in self.selected:
                    self.selected.append(idx)
                    self.selected.sort()
                else:
                    self.selected.remove(idx)
                    self.reset_color_deselected(idx)

                self.change_color_selected()
                self.text_selected()

        event.canvas.draw()

    def on_key(self, event):
        if event.key in ["down", "ctrl+j"]:
            if self.idx_tmp == len(self.rects2) - 1:
                self.idx_tmp = -1
            else:
                self.idx_tmp += 1

        elif event.key in ["up", "ctrl+k"]:
            if self.idx_tmp == -1:
                self.idx_tmp = len(self.rects2) - 1
            else:
                self.idx_tmp -= 1

        elif event.key in ["right", "ctrl+l"]:
            if self.idx_tmp > -1:
                if self.idx_tmp not in self.selected:
                    self.selected.append(self.idx_tmp)
                    self.selected.sort()

            self.change_color_selected()
            self.text_selected()

        elif event.key in ["left", "ctrl+h"]:
            if self.idx_tmp > -1:
                if self.idx_tmp in self.selected:
                    self.selected.remove(self.idx_tmp)
                    self.reset_color_deselected(self.idx_tmp)

            self.change_color_selected()
            self.text_selected()

        elif event.key in ["ctrl+c"]:
            self.disconnect()
            sys.exit()

        if self.idx_tmp > -2:
            self.reset_info()
            self.show_info(self.idx_tmp)

        event.canvas.draw()

    def get_rects(self, ax):
        rects = []
        for container in ax.containers:
            for child in container.get_children():
                if isinstance(child, plt.Rectangle):
                    rects.append(child)
        return rects

    def text_selected(self):
        if self.txt0:
            self.txt0.remove()
        if len(self.selected) > 0:
            s = r"$\bf{Selected}$" + "\n"
            if self.labels:
                s += ", ".join(
                    [self.labels[x].split(";")[0].rstrip() for x in self.selected]
                )
            else:
                s += ", ".join(["Seq{}".format(x) for x in self.selected])
            c = "k"
        else:
            s = "Left click on a sequence  (or press \u2191/\u2193 key) to "
            s += "show information \n right click (or press \u2190/\u2192"
            s += r" key) to select the sequence to $\bf{keep}$"
            c = "grey"

        self.txt0 = self.ax0.text(
            0.5, 0.6, s=s, color=c, va="center", ha="center", wrap=True
        )

    def show_info(self, idx):
        if idx == -1:
            return
        else:
            # ax1
            rect1 = self.rects1[idx]
            rect1.set_edgecolor("black")
            rect1.set_linewidth(2)
            rect1.set_zorder(max(self.zorder_org) + 1)
            label = self.labels[idx] if self.labels else "Seq{}".format(idx)
            x0, y0 = rect1.get_xy()
            h, w = rect1.get_height(), rect1.get_width()
            x, y = x0 + w * 0.5, y0 + h * 0.5
            xl0, xl1 = self.ax1.get_xlim()
            xytext = (xl0 + (xl1 - xl0) * 0.04, self.ax1.get_ylim()[1] * 0.96)
            ap = dict(
                arrowstyle="-",
                color="k",
                patchB=rect1,
                shrinkB=0,
                linewidth=1.2,
                connectionstyle="arc3,rad=0.2",
            )
            bb = dict(boxstyle="round", fc="cyan", ec=(1, 1, 1, 0), linewidth=1.2)
            self.ann1 = self.ax1.annotate(
                label,
                xy=(x, y),
                xytext=xytext,
                color="navy",
                va="top",
                arrowprops=ap,
                bbox=bb,
            )

            # ax2
            if idx < 10:
                rect2 = self.rects2[idx]
                rect2.set_facecolor((0, 0, 0, 0.2))
                ytl2 = [
                    r"$\bf{" + ytl + "}$" if i == (idx + 1) else ytl
                    for i, ytl in enumerate(self.ytl2)
                ]
                self.ax2.set_yticklabels(ytl2)

    def reset_info(self):
        if self.ann1:
            self.ann1.remove()
            self.ann1 = ""
        for rect1, z in zip(self.rects1, self.zorder_org):
            rect1.set_zorder(z)
            rect1.set_edgecolor("white")
            rect1.set_linewidth(1)
        for rect2 in self.rects2:
            rect2.set_facecolor((1, 1, 1, 0))
        self.ax2.set_yticklabels(self.ytl2)

    def change_color_selected(self):
        if len(self.selected) > 0:
            for idx in set(self.selected):
                self.rects1[idx].set_facecolor("mediumvioletred")
                if idx < 10:
                    self.rects2[idx].set_edgecolor("mediumvioletred")
                    self.rects2[idx].set_zorder(100)

    def reset_color_deselected(self, idx):
        self.rects1[idx].set_facecolor(self.cols_org[idx])
        if idx < 10:
            self.rects2[idx].set_edgecolor((1, 1, 1, 0))
            self.rects2[idx].set_zorder(1)

    def reset_color_all(self):
        if len(self.selected) > 0:
            for i in self.selected:
                self.rects1[i].set_facecolor(self.cols_org[i])
                self.rects2[i].set_edgecolor((1, 1, 1, 0))
                self.rects2[i].set_zorder(1)

    def enter_axes(self, event):
        x0, y0 = event.inaxes.get_position().min
        if (x0, y0) == tuple(self.ax1.get_position().min):
            self.subplot = 1
        elif (x0, y0) == tuple(self.ax2.get_position().min):
            self.subplot = 2

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
        try:
            self.connect()
            plt.show()
        finally:
            self.disconnect()
            plt.close()


def pairwise_dist_Levenstein(seq_list):
    n = len(seq_list)
    dist_mat = [
        [Levenshtein.distance(seq_list[i], seq_list[j]) for j in range(i + 1)]
        for i in range(n)
    ]
    return dist_mat


def main(args):
    print("Start variant filtering ... ({})".format(time.ctime()))
    t1 = time.time()

    if isinstance(args, dict):
        kwargs = args
    else:
        kwargs = vars(args)

    if "subcommand" in kwargs:
        kwargs.pop("subcommand")

    if get_backend() == "agg":
        kwargs["force_visual_check"] = False
        kwargs["force_no_visual_check"] = True

    quiet = kwargs["quiet"]
    if not quiet:
        print_args(kwargs)

    vrf = VariantFilter(**kwargs)
    vrf.run()

    t2 = time.time()
    elapsed_time = time.strftime("%H:%M:%S", time.gmtime(t2 - t1))
    print("\nFinished (elapsed time: {})".format(elapsed_time))


if __name__ == "__main__":
    argv = ["filter"]
    argv.extend(sys.argv[1:])
    main(get_args(argv))
