from __future__ import annotations

import copy
import json
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

import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from matplotlib import get_backend
from natsort import natsorted
from tqdm import tqdm

from .argument_parser import ask_user, get_args, print_args
from .base import GenData, SeqData, count_records, gen_dict_from_table, read_fastx
from .utils import run_blastn
from .variant_filter import VariantFilter, VisualCheck, pairwise_dist_Levenstein


class AlleleCall(VariantFilter):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        if self.glob_pattern.find("fastq") > -1:
            self.fmt = "fastq"
        elif self.glob_pattern.find("fasta") > -1:
            self.fmt = "fasta"
        else:
            self.fmt = "auto"

    def filter(self, seq_file):
        n_reads = count_records(seq_file, self.fmt)
        if not n_reads:
            return 1

        marker = seq_file.parent.name
        pat_rm = "{}.*".format(self.glob_pattern.strip("*").replace("*", ".*"))
        idx = re.sub(pat_rm, "", seq_file.name)
        idx = idx.replace("_" + marker, "")

        # Run blastn against allele database
        ftmp, tmppath = tempfile.mkstemp()
        try:
            with os.fdopen(ftmp, "w") as tmpfile:
                SeqIO.write(list(read_fastx(seq_file)), tmpfile, "fasta")
            blast_out = self.resdir.joinpath(idx + "_blast.csv")
            run_blastn(tmppath, self.blastdb, str(blast_out))
        finally:
            os.remove(tmppath)

        # Read results of the blast search
        hits = read_blast_results(blast_out)
        candidates = [x for x in copy.deepcopy(self.alleledat) if x.id in hits]
        for x in candidates:
            x.counts = hits[x.id]

        # Filtering
        filtered = []
        visual_check = False
        groups = group_alleles(candidates)
        for g in groups:
            counts = [x.counts for x in g]
            keep = [x for x in g if x.counts > max(counts) * self.threshold]
            if len(keep) > self.max_alleles:
                visual_check = True
            elif len(keep) == self.max_alleles:
                n_reps = np.array([[y.n_reps for y in x.rep_data] for x in keep])
                neighbor = np.where(
                    np.sum(np.abs(np.diff(n_reps, axis=0)), axis=1) <= 1
                )[0]
                if len(neighbor) > 0:
                    visual_check = True
                # if len(neighbor) == 1:
                #     if cs[neighbor[0] + 1] - cs[neighbor[0]] > 0:
                #         visual_check = True
            filtered.extend(keep)

        counts = [x.counts for x in filtered]
        filtered = [x for x in filtered if x.counts > max(counts) * 0.05]
        filtered = sorted(filtered, key=lambda x: x.counts)[::-1]

        selected_idx = [candidates.index(x) for x in filtered]
        if len(selected_idx) > self.max_alleles:
            selected_idx = selected_idx[: self.max_alleles]

        n_filtered = len(filtered)
        if n_filtered > self.max_alleles:
            visual_check = True
        elif n_filtered < self.max_alleles:
            if sum([x.counts for x in filtered]) < self.min_reads:
                # rests to missing
                selected_idx.extend([-1] * (self.max_alleles - n_filtered))
            else:
                if self.max_alleles == 2:
                    # homozygote
                    selected_idx = selected_idx * 2
                else:
                    visual_check = True
        else:
            # check if candidate alleles are close to each other
            seqs = [x.seq for x in filtered]
            l_dists = np.array(list(chain(*pairwise_dist_Levenstein(seqs))))
            if np.any(l_dists <= 3):
                visual_check = True

        n_reads_each = {
            candidates[i].id: candidates[i].counts for i in set(selected_idx) if i > -1
        }
        genotype = [candidates[i].id if i > -1 else "NA" for i in selected_idx]

        if self.force_no_visual_check:
            if self.n_cpu > 1:
                self.mp_genotypes.append(GenData(idx, n_reads, genotype, n_reads_each))
            else:
                self.genotypes.append(GenData(idx, n_reads, genotype, n_reads_each))
            return 0
        elif visual_check or self.force_visual_check:
            if self.n_cpu > 1:
                self.mp_queue.append((idx, n_reads, candidates, selected_idx))
            else:
                self.queue.append((idx, n_reads, candidates, selected_idx))
            return 1
        else:
            if self.n_cpu > 1:
                self.mp_genotypes.append(GenData(idx, n_reads, genotype, n_reads_each))
            else:
                self.genotypes.append(GenData(idx, n_reads, genotype, n_reads_each))
            return 0

    def visual_check(self, idx, n_reads, candidates, selected_idx):
        while True:
            vsl = VisualSelection(
                candidates, suptitle=idx, max_alleles=self.max_alleles
            )
            while selected_idx:
                vsl.selected.append(selected_idx.pop())
            vsl.text_selected()
            vsl.change_color_selected()
            vsl.show()
            vsl.disconnect()
            if vsl.selected:
                if len(vsl.selected) == self.max_alleles:
                    break
                else:
                    msg = (
                        "The number of alleles selected ({}) is less than expected "
                        "({}).\nDo you want the rest of the alleles to be missing?"
                    )
                    msg = msg.format(len(vsl.selected), self.max_alleles)
                    if ask_user(msg, default="n"):
                        vsl.selected.extend(
                            [-1] * (self.max_alleles - len(vsl.selected))
                        )
                        break
                    else:
                        continue
            elif ask_user("No allele was selected. Proceed anyway?"):
                break

        genotype = [candidates[i].id if i > -1 else "NA" for i in vsl.selected]
        n_reads_each = {
            candidates[i].id: candidates[i].counts for i in set(vsl.selected) if i > -1
        }

        genotype = [candidates[i].id for i in vsl.selected]
        n_reads_each = {
            candidates[i].id: candidates[i].counts for i in set(vsl.selected)
        }
        if len(genotype) < self.max_alleles:
            genotype.extend(["NA"] * (self.max_alleles - len(genotype)))

        self.genotypes.append(GenData(idx, n_reads, genotype, n_reads_each))

    def write_genotype_table(self):
        # write genotype list to json file
        self.genotypes = natsorted(self.genotypes, key=lambda x: x.id)
        json_str = ",\n".join(
            [x.to_json(indent=4, ensure_ascii=False) for x in self.genotypes]
        )
        json_str = "[\n{}\n]".format(json_str)
        with self.outpath.open("w") as outfile:
            outfile.write(json_str)

    def process_all_files_in_subdir(self, subdir):
        if not self.quiet:
            print("." * shutil.get_terminal_size()[0])
            print("Processing {} ({}): ".format(subdir.name, time.ctime()))

        self.max_alleles = self.dict_max_alleles[subdir.name]
        self.resdir = self.indir.joinpath("result_summary", "allele_call", subdir.name)
        self.resdir.mkdir(parents=True, exist_ok=True)

        self.outpath = subdir.joinpath(subdir.name + "_genotype.json")
        self.blastdb = str(subdir.joinpath(subdir.name + "_allele_blastdb"))
        path_allele_data = subdir.joinpath(subdir.name + "_allele_data.json")

        if path_allele_data.exists():
            self.alleledat = read_allele_data(path_allele_data)
        else:
            print("Allele database file not found")
            return

        if self.outpath.exists():
            msg = "Output file '{}' already exists. Overwrite?"
            msg = msg.format(str(self.outpath))
            if ask_user(msg, default="n"):
                self.outpath.unlink()
                self.overwrite = True
            else:
                return

        self.queue = []
        self.genotypes = []
        files = natsorted(subdir.glob(self.glob_pattern), key=lambda x: x.name)
        if self.n_cpu == 1:
            if self.quiet:
                res = [self.filter(f) for f in files]
            else:
                res = [self.filter(f) for f in tqdm(files, total=len(files))]
        else:
            with Manager() as manager:
                self.mp_queue = manager.list()
                self.mp_genotypes = manager.list()
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
                    while self.mp_genotypes:
                        self.genotypes.append(self.mp_genotypes.pop())

        if np.sum(res) > 0:
            self.queue = natsorted(self.queue, key=lambda x: x[0])[::-1]
            if not self.quiet:
                print("\nStart visual-checking ...")
            while self.queue:
                self.visual_check(*self.queue.pop())

        self.write_genotype_table()

        if not self.force_no_visual_check:
            if subdir != self.subdirs[-1] and np.sum(res) > 0:
                if not ask_user("Done. Proceed to the next?", default="y", quit=False):
                    sys.exit()

    def run(self):
        for sd in self.subdirs:
            self.process_all_files_in_subdir(sd)


class VisualSelection(VisualCheck):
    def __init__(self, *args, max_alleles=2, **kwargs):
        super().__init__(*args, **kwargs)
        self.numbers = []
        self.max_alleles = max_alleles
        self.ax1.set_title("Distribution of Sequence Lengths of Alleles")

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
            if self.numbers:
                for n in self.numbers:
                    n.remove()
                self.numbers = []

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
                if idx in self.selected:
                    if len(self.selected) < self.max_alleles:
                        self.selected.append(idx)
                        self.selected.sort()
                    elif len(self.selected) == self.max_alleles:
                        self.selected = [x for x in self.selected if x != idx]
                        self.reset_color_deselected(idx)
                else:
                    if len(self.selected) < self.max_alleles:
                        self.selected.append(idx)
                        self.selected.sort()
                    else:
                        pass

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
                if self.idx_tmp in self.selected:
                    if len(self.selected) < self.max_alleles:
                        self.selected.append(self.idx_tmp)
                        self.selected.sort()
                else:
                    if len(self.selected) < self.max_alleles:
                        self.selected.append(self.idx_tmp)
                        self.selected.sort()

            self.change_color_selected()
            self.text_selected()

        elif event.key in ["left", "ctrl+h"]:
            if self.idx_tmp > -1:
                if self.idx_tmp in self.selected:
                    self.selected.remove(self.idx_tmp)
                    if self.selected.count(self.idx_tmp) == 0:
                        self.reset_color_deselected(self.idx_tmp)

            self.change_color_selected()
            self.text_selected()

        elif event.key in ["ctrl+c"]:
            sys.exit()

        if self.idx_tmp > -2:
            self.reset_info()
            self.show_info(self.idx_tmp)

        event.canvas.draw()

    def change_color_selected(self):
        if self.numbers:
            for n in self.numbers:
                n.remove()
            self.numbers = []

        if len(self.selected) > 0:
            for idx in set(self.selected):
                rect1 = self.rects1[idx]
                rect1.set_facecolor("mediumvioletred")
                x0, y0 = rect1.get_xy()
                h, w = rect1.get_height(), rect1.get_width()
                x, y = x0 + w * 0.5, y0 + h * 0.5
                n = self.ax1.text(
                    x,
                    y,
                    s=str(self.selected.count(idx)),
                    zorder=100,
                    color="k",
                    fontweight="bold",
                    va="center",
                    ha="center",
                )
                self.numbers.append(n)

                if idx < 10:
                    self.rects2[idx].set_edgecolor("mediumvioletred")
                    self.rects2[idx].set_zorder(100)

    def reset_color_deselected(self, idx):
        self.rects1[idx].set_facecolor(self.cols_org[idx])
        if idx < 10:
            self.rects2[idx].set_edgecolor((1, 1, 1, 0))
            self.rects2[idx].set_zorder(1)


def read_allele_data(json_file):
    with open(json_file, "r") as infile:
        json_data = json.load(infile)
        return [SeqData.from_dict(i) for i in json_data]


def read_blast_results(filepath):
    """
    Read blast results and return counts of blast search hits
    """
    hits = {}
    with open(filepath) as infile:
        line = infile.readline()
        while line:
            hit = line.strip().split(",")[1]
            if hit in hits:
                hits[hit] += 1
            else:
                hits[hit] = 1
            line = infile.readline()
    return dict(sorted(hits.items(), key=lambda x: x[1])[::-1])


def group_alleles(alleledat):
    groups = []
    for al in alleledat:
        if not groups:
            groups.append([al])
            continue

        stutters = al.stutter_s + al.stutter_l

        overlapping_groups = [
            g for g in groups if np.any([x in [y.id for y in g] for x in stutters])
        ]
        if not overlapping_groups:
            groups.append([al])
        elif len(overlapping_groups) == 1:
            overlapping_groups[0].append(al)
        else:
            new_group = [al]
            for group in overlapping_groups:
                for al in group:
                    if al not in new_group:
                        new_group.append(al)
            groups = [g for g in groups if g not in overlapping_groups]
            groups.append(new_group)
    return groups


def read_genotype_data(json_file):
    with open(json_file, "r") as infile:
        json_data = json.load(infile)
        return [GenData.from_dict(i) for i in json_data]


def make_genotype_table(indir, path_marker_data=None, **kwargs):
    indir = Path(indir)

    path_genotype_files = natsorted(indir.glob("**/*_genotype.json"))
    all_markers = [f.parent.name for f in path_genotype_files]

    if path_marker_data:
        marker_group = gen_dict_from_table(path_marker_data, "Name", "Group")
        groups = natsorted(set(marker_group.values()))

        if len(groups) == 1 and not marker_group[groups[0]]:
            marker_group = {k: "all" for k in marker_group.keys()}
            groups = natsorted(set(marker_group.values()))

    else:
        marker_group = {k: "all" for k in all_markers}

    sample_list = {}
    genotypes = {}
    for f in path_genotype_files:
        marker = f.name.replace("_genotype.json", "")
        group = marker_group[marker]
        if not group:
            group = "all"
        genotype = read_genotype_data(f)
        if group in genotypes:
            genotypes[group].update({marker: genotype})
        else:
            genotypes[group] = {marker: genotype}
        for g in genotype:
            idx = g.id.replace("_" + marker, "")
            if group in sample_list:
                sample_list[group].update([idx])
            else:
                sample_list[group] = set([idx])

    for group in genotypes:
        outpath = indir.joinpath(group + "_genotype_table.csv")

        markers = []
        sample_ids = natsorted(sample_list[group])
        gen_table = {i: [] for i in sample_ids}

        for marker, genotype in genotypes[group].items():
            markers.append(marker)
            path_allele_data = list(
                indir.glob("**/{}_allele_data.json".format(marker))
            )[0]
            allele_data = read_allele_data(path_allele_data)
            alleles = [x.id for x in allele_data]
            for g in genotype:
                genid = [alleles.index(x) + 1 if x != "NA" else 0 for x in g.genotype]
                genstr = "".join(["{:03}".format(x) for x in genid])
                gen_table[g.id].append(genstr)

        with outpath.open("w") as outfile:
            outfile.write(",{}\n".format(",".join(markers)))
            for k, v in gen_table.items():
                outfile.write("{},{}\n".format(k, ",".join(v)))


def main(args):
    print("Start allele calling ... ({})".format(time.ctime()))
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

    if get_backend() == "agg":
        kwargs["force_no_visual_check"] = True

    quiet = kwargs["quiet"]
    if not quiet:
        print_args(kwargs)

    acl = AlleleCall(**kwargs)
    acl.run()

    make_genotype_table(**kwargs)

    t2 = time.time()
    elapsed_time = time.strftime("%H:%M:%S", time.gmtime(t2 - t1))
    print("\nFinished (elapsed time: {})".format(elapsed_time))


if __name__ == "__main__":
    argv = ["allele-call"]
    argv.extend(sys.argv[1:])
    main(get_args(argv))
