from __future__ import annotations

import sys
from pathlib import Path

from matplotlib import get_backend

from .argument_parser import get_args
from .base import SeqData, count_uniq_seq
from .find_ssrs import find_ssrs, get_longest_RepData
from .variant_filter import VisualCheck


class ShowAlignment(VisualCheck):
    def __init__(self, seqdat, suptitle):
        super().__init__(seqdat, suptitle=suptitle)
        self.reset_button.remove()

    def on_pick(self, event):
        artist = event.artist
        mouseevent = event.mouseevent

        if self.subplot is None:
            return

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

        elif event.key in ["ctrl+c"]:
            self.disconnect()
            sys.exit()

        if self.idx_tmp > -2:
            self.reset_info()
            self.show_info(self.idx_tmp)

        event.canvas.draw()

    def text_selected(self):
        if self.txt0:
            self.txt0.remove()
        s = "Left click on a sequence (or press \u2191/\u2193 key) to show information"
        self.txt0 = self.ax0.text(
            0.5, 0.6, s=s, color="grey", va="center", ha="center", wrap=True
        )


def main(args):
    if get_backend == "agg":
        sys.exit()

    for seq_file in args.infile:
        seqdat = [
            SeqData("Seq{}".format(i), s, c)
            for i, (s, c) in enumerate(count_uniq_seq(seq_file).items())
        ]
        for s in seqdat:
            s.rep_data = get_longest_RepData(find_ssrs(s.seq))

        sal = ShowAlignment(seqdat, Path(seq_file).name)
        sal.show()


if __name__ == "__main__":
    argv = ["show-alignment"]
    argv.extend(sys.argv[1:])
    main(get_args(argv))
