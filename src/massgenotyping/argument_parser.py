from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path
from typing import Literal

here = Path(__file__).parent.resolve()

version: dict[str, str | dict] = {}
with (here / "__about__.py").open() as fp:
    exec(fp.read(), version)


def get_args(argv=None):
    # top-level parser
    parser = argparse.ArgumentParser(
        prog="mgt",
        description="massgenotyping v{}".format(version["__version__"]),
        epilog=(
            "show subcommand help:\n"
            "  %(prog)s SUBCOMMAND -h\n\n"
            "Any feedback would be greatly appreciated!\n"
            "Homepage: https://github.com/kohyamat/massgenotyping"
        ),
        formatter_class=MyFormatter,
    )

    subparsers = parser.add_subparsers(
        dest="subcommand",
        metavar="SUBCOMMAND",
    )

    # parser for the "demultiplex" command
    parser_a = subparsers.add_parser(
        "demultiplex",
        help="Demultiplex raw amplicon sequences",
        description="Demultiplex raw amplicon sequences based on primer sequences",
        usage=(
            "%(prog)s [options] read1.fastq read2.fastq -m marker_data.csv (paired-end)"
            "\n       %(prog)s [options] read1.fastq -m marker_data.csv (single-end)"
            "\n       %(prog)s [options] read1.fastq read2.fastq -f primer_fwd.fasta "
            "-r primer_rev.fasta"
            "\n       %(prog)s [options] INDIR -g '*_R[12]*' -m marker_data.csv "
            "(sequential processing)"
        ),
        formatter_class=MyFormatter2,
    )
    parser_a.add_argument(
        "inpath",
        nargs="+",
        metavar="INPATH",
        help=(
            "Path to one or two sequence files in the FASTQ format (support gzip "
            "compression), or the directory containing the sequence files. If the "
            "input consists of two files, the input data is considered as paired-end; "
            "if it consists of a single file, it is considered as single-end. If the "
            "input is a directory, all sequence files in the directory are "
            "automatically processed. In this case, the --glob-pattern option is used "
            "to find sequence files in the directory."
        ),
        action=MaxTwoArgs,
    )
    parser_a.add_argument(
        "-g",
        "--glob-pattern",
        default="*_R[12]_*",
        metavar="STR",
        help=(
            "A wildcard pattern used to find sequence files in the input directory. "
            "The file name string before the pattern will be used as the prefix of "
            "the output files. (default: '*_R[12]_*')"
        ),
    )
    parser_a.add_argument(
        "-m",
        "--path-marker-data",
        default=None,
        metavar="FILE",
        help=(
            "Path to the CSV file of the marker information table which contains the "
            "fields of marker names and sequences of forward and reverse primers. "
            "If this option is used, -f and -r arguments are ignored."
        ),
    )
    parser_a.add_argument(
        "-f",
        "--path-fwd-primers",
        default=None,
        metavar="FILE",
        help="Path to the sequence file of forward primers in the FASTA format",
    )
    parser_a.add_argument(
        "-r",
        "--path-rev-primers",
        default=None,
        metavar="FILE",
        help="Path to the sequence file of reverse primers in the FASTA format",
    )
    parser_a.add_argument(
        "-o",
        "--outdir",
        default="project",
        metavar="DIR",
        help="Path to the output directory (default: %(default)s)",
    )
    parser_a.add_argument(
        "-O",
        "--out-prefix",
        default="",
        metavar="STR",
        help="Prefix for the output file name (default: None)",
    )
    parser_a.add_argument(
        "-M",
        "--max-mismatch",
        type=check_positive,
        default=[0.12, 0.14],
        metavar="INT or FLOAT",
        nargs="+",
        help=(
            "Maximum number (or propotion) of mismatches allowed for searching primer "
            "sequences. Primer sequences at the 5'- and 3'-ends will be searched with "
            "different mismatch toleraces when two space-separated values are given. "
            "(default: 0.12 0.14)"
        ),
        action=MaxTwoArgs,
    )
    parser_a.add_argument(
        "--n-cpu",
        type=int,
        default=1,
        metavar="INT",
        help="Number of processers (default: %(default)i)",
    )
    parser_a.add_argument(
        "--trim-primer",
        action="store_true",
        help="Trim primer sequences",
    )
    parser_a.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Reduce stdout messages",
    )

    # parser for the "merge-pairs" command
    parser_b = subparsers.add_parser(
        "merge-pairs",
        help="Merge paired-end reads",
        description=(
            "Merge paired-end reads using NGmerge "
            "(Gaspar 2018, BMC Bioinformatics, 19(1):536)"
        ),
        usage="%(prog)s [options] INDIR",
    )
    parser_b.add_argument(
        "indir",
        metavar="INDIR",
        help=(
            "Path to the input directory that contains sequence files of paierd-end "
            "reads with suffixes '_R1.fastq.gz' and '_R2.fastq.gz'. Normally, the "
            "output directory created by the 'demultiplex' program"
        ),
    )
    parser_b.add_argument(
        "-o",
        "--outdir",
        default=None,
        metavar="DIR",
        help=(
            "Path to the output directory. Output files will be created within the "
            "input directory when the --outdir option was not used. (default: "
            "%(default)s)"
        ),
    )
    parser_b.add_argument(
        "-t",
        "--quality-threshold",
        type=float,
        default=32,
        metavar="INT",
        help=(
            "A cutoff threshold for quality trimming. The 3' end of the sequence reads "
            "will be trimmed before merging read pairs when the value > 0. (default: "
            "%(default)s)"
        ),
    )
    parser_b.add_argument(
        "-w",
        "--window-size",
        type=int,
        default=8,
        metavar="INT",
        help="Window size for quality trimming (default: %(default)i)",
    )
    parser_b.add_argument(
        "-s",
        "--step-size",
        type=int,
        default=1,
        metavar="INT",
        help="Step size for quality trimming (default: %(default)i)",
    )
    parser_b.add_argument(
        "-p",
        "--trim-primer",
        action="store_true",
        help=(
            "Trim primer sequences after merging read pairs. Use this option with "
            "-m (or -f & -r) option. (default: False)"
        ),
    )
    parser_b.add_argument(
        "-m",
        "--path-marker-data",
        default=None,
        metavar="FILE",
        help=(
            "Path to the CSV file of the marker information table which contains the "
            "fields of marker names and sequences of forward/reverse primers. "
            "-f and -r arguments will be ignored if this options is used."
        ),
    )
    parser_b.add_argument(
        "-f",
        "--path-fwd-primers",
        default=None,
        metavar="FILE",
        help="Path to the sequence file of forward primers in the FASTA format",
    )
    parser_b.add_argument(
        "-r",
        "--path-rev-primers",
        default=None,
        metavar="FILE",
        help="Path to the sequence file of reverse primers in the FASTA format",
    )
    parser_b.add_argument(
        "-M",
        "--max-mismatch-primer",
        type=check_positive,
        default=0.12,
        metavar="NUM",
        help=(
            "Maximum number (or proportion) of mismatches allowed for searching "
            "primer sequences (default: %(default)s)"
        ),
    )
    parser_b.add_argument(
        "--n-cpu",
        type=int,
        default=1,
        metavar="INT",
        help="Number of processers (default: %(default)i)",
    )
    parser_b.add_argument(
        "-q", "--quiet", action="store_true", help="Reduce stdout messages"
    )

    # parser for the "denoise" command
    parser_c = subparsers.add_parser(
        "denoise",
        help="Reduce putative sequence errors",
        description=(
            "Replace rare bases, which probably be caused due to errors in base "
            "calling in sequencing or PCR, with the most frequent base at each "
            "presumably homogeneous nucleotide position across the multiple "
            "sequence alingment. The program MAFFT is required."
        ),
        usage=(
            "%(prog)s [options] FILE [FILE, ...]"
            "\n       %(prog)s [options] INDIR -g '*_merged.fastq.gz*'"
        ),
    )
    parser_c.add_argument(
        "inpath",
        metavar="INPATH",
        nargs="+",
        help=(
            "Path to the sequences file(s) in FASTA/FASTQ format (support gzip "
            "compression), or a directory  which contains sequence files."
            "If the input is a directory, the --glob-pattern option is used to find "
            "sequence files"
        ),
    )
    parser_c.add_argument(
        "-g",
        "--glob-pattern",
        default="*_merged.fastq.gz",
        metavar="STR",
        help=(
            "A wildcard pattern to find out sequence files in the input directory. "
            "The character string of the file name before the pattern will be used "
            "as the prefix of the output files. (default: '*_merged.fastq.gz')"
        ),
    )
    parser_c.add_argument(
        "-o",
        "--outdir",
        default=None,
        metavar="DIR",
        help=(
            "Path to the output directory. By default, output files will be created "
            "within the input directory."
        ),
    )
    parser_c.add_argument(
        "-O",
        "--out-suffix",
        metavar="STR",
        default="_denoised",
        help="Suffix for the output file name (default: %(default)s)",
    )
    parser_c.add_argument(
        "-a",
        type=check_0to1,
        default=0.1,
        metavar="FLOAT",
        help=(
            "A threshold for detecting erroneous bases (default: %(default)s). Bases "
            "that are less frequent than the threshold are considered erroneous and "
            "are replaced by the most frequent base. Replacements are only performed "
            "at nucleotide positions that meet the condition of the -b option."
        ),
    )
    parser_c.add_argument(
        "-b",
        type=check_0to1,
        default=0.9,
        metavar="FLOAT",
        help=(
            "A threshold for detectng presumably homogeneous nucleotide positions "
            "(default: %(default)s). Nucleotide positions are consider as homogeneous "
            "if the frequency of the most frequenct base is higher than the threshold."
        ),
    )
    parser_c.add_argument(
        "--localpair",
        action="store_true",
        help=(
            "Use the --localpair option of MAFFT to align sequences. "
            "More accurate but slower"
        ),
    )
    parser_c.add_argument(
        "--globalpair",
        action="store_true",
        help=(
            "Use the --globalpair option of MAFFT to align sequences. "
            "More accurate but slower"
        ),
    )
    parser_c.add_argument(
        "--n-cpu",
        type=int,
        metavar="INT",
        default=1,
        help="Number of processers (default: %(default)i)",
    )
    parser_c.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Reduce stdout messages",
    )

    # parser for the "filter" command
    parser_d = subparsers.add_parser(
        "filter",
        help="Filtering out erroreous sequence variants",
        description=(
            "Filtering erroreous sequence variants such as 'stutter' fragments "
            "and PCR chimeras, and extract putative aleleles"
        ),
        usage="%(prog)s [options] INDIR",
    )
    parser_d.add_argument(
        "indir",
        metavar="INDIR",
        help="Path to the input directory",
    )
    parser_d.add_argument(
        "-g",
        "--glob-pattern",
        default="*_denoised.fasta.gz",
        metavar="STR",
        help=(
            "A wildcard pattern to find out sequence files in the input directory. "
            "The character string of the file name before the pattern will be used "
            "as the prefix of the output files. (default: '*_denoised.fasta.gz')"
        ),
    )
    parser_d.add_argument(
        "-s",
        "--subdir",
        default=None,
        metavar="DIR",
        nargs="+",
        help="Names of subdirectories (i.e. marker names) to be processed",
    )
    parser_d.add_argument(
        "-o",
        "--outdir",
        default=None,
        metavar="DIR",
        help=(
            "Path to the output directory. Output files will be created within the "
            "input directory when the --outdir option was not used. (default: "
            "%(default)s)"
        ),
    )
    parser_d.add_argument(
        "-t",
        "--threshold",
        type=check_0to1,
        default=0.25,
        metavar="FLOAT",
        help=(
            "A cutoff threshold for rare sequence variants. The sequence variants "
            "will be discarded when thier relative frequencies to the most frequent "
            "variant are lower than the threshold. (default: %(default)s)"
        ),
    )
    parser_d.add_argument(
        "-n",
        "--max_alleles",
        type=int,
        default=2,
        metavar="INT",
        help=(
            "Maximum number of alleles amplified in each marker. If it differs between "
            "markers, use --path-marker-data option for set the value for each marker. "
            "If number of unique sequence variants was greater the expected number of "
            "alleles, the visual-checking procedure will be launch. "
            "(default: %(default)i)"
        ),
    )
    parser_d.add_argument(
        "-m",
        "--path-marker-data",
        default=None,
        metavar="FILE",
        help=(
            "Path to the CSV file of the marker information table which contains the "
            "fields of marker names and maximum number of alleles for each marker"
        ),
    )
    parser_d.add_argument(
        "-l",
        "--min_l_dist",
        type=int,
        default=3,
        metavar="INT",
        help=(
            "Minimum Levenshtein distance between seqeunce variants to consider them "
            "as different alleles automatically. If two or more seqeunce variants are "
            "very close to each other and at least one of them is newly detected, the "
            "visual-checking procedure will be launch. (default: %(default)s)"
        ),
    )
    parser_d.add_argument(
        "-M",
        "--min_reads",
        type=int,
        default=6,
        help="Minimum number of reads required for proccessing (default: %(default)s)",
    )
    exclusive_group_d = parser_d.add_mutually_exclusive_group(required=False)
    exclusive_group_d.add_argument(
        "--force-visual-check",
        action="store_true",
        help="Perform visual checking for all samples",
    )
    exclusive_group_d.add_argument(
        "--force-no-visual-check",
        action="store_true",
        help=(
            "Do not perform visual checking. This will keep all sequence variants that "
            "meet the critera"
        ),
    )
    parser_d.add_argument(
        "--n-cpu",
        type=int,
        default=1,
        metavar="INT",
        help="Number of processers (default: %(default)i)",
    )
    parser_d.add_argument(
        "-q", "--quiet", action="store_true", help="Reduce stdout messages"
    )

    # parser for the "allele-check" command
    parser_e = subparsers.add_parser(
        "allele-check",
        help="Check allele candidates and create an allele database",
        description=(
            (
                "Check the sequence alingment of allele candidates and create an "
                "allele database"
            )
        ),
        usage="%(prog)s [options] INDIR",
    )
    parser_e.add_argument(
        "indir",
        metavar="INDIR",
        help="Path to the input directory",
    )
    parser_e.add_argument(
        "-s",
        "--subdir",
        metavar="DIR",
        nargs="+",
        help="Names of subdirectories (i.e. marker names) to be processed",
    )
    parser_e.add_argument(
        "-g",
        "--glob-pattern",
        default="*_filtered.fasta.gz",
        metavar="STR",
        help=(
            "A wildcard pattern to find out sequence files in the input directory. "
            "The character string of the file name before the pattern will be used "
            "as the prefix of the output files. (default: "
            "'*_filtered.fasta.gz')"
        ),
    )
    parser_e.add_argument(
        "-o",
        "--outdir",
        default=None,
        metavar="DIR",
        help=(
            "Path to the output directory. Output files will be created within the "
            "input directory when the --outdir option was not used. (default: "
            "%(default)s)"
        ),
    )
    parser_e.add_argument(
        "--min-repeats",
        default=3,
        metavar="INT",
        help="Minimum numbers of repeats needs to be checked (default: %(default)s)",
    )
    parser_e.add_argument(
        "--min-motif-len",
        default=2,
        metavar="INT",
        help=(
            "Minimum length of the repeat motif needs to be checked "
            "(default: %(default)s)"
        ),
    )
    parser_e.add_argument(
        "--max-motif-len",
        default=4,
        metavar="INT",
        help=(
            "Maximum length of the repeat motif needs to be checked "
            "(default: %(default)s)"
        ),
    )
    parser_e.add_argument(
        "--force-no-visual-check",
        action="store_true",
        help=("Do not perform visual checking"),
    )
    parser_e.add_argument(
        "-q", "--quiet", action="store_true", help="Reduce stdout messages"
    )

    # parser for the "allele-call" command
    parser_f = subparsers.add_parser(
        "allele-call",
        help="Assign alleles to raw amplicon sequences",
        description=(
            "Perform BLASTn search against the database created for each marker and "
            "assign alleles to the raw sequencing data."
        ),
        usage="%(prog)s [options]",
    )
    parser_f.add_argument(
        "indir",
        metavar="INDIR",
        help="Path to the input directory",
    )
    parser_f.add_argument(
        "-g",
        "--glob-pattern",
        default="*_merged.fastq.gz",
        metavar="PATTERN",
        help=(
            "A wildcard pattern to find sequence files in the input directory"
            " (default: '*_merged.fastq.gz')"
        ),
    )
    parser_f.add_argument(
        "-s",
        "--subdir",
        metavar="DIR",
        nargs="+",
        help="Names of subdirectories (i.e. marker names) to be processed",
    )
    parser_f.add_argument(
        "-o",
        "--outdir",
        default=None,
        metavar="DIR",
        help=(
            "Path to the output directory. Output files will be created within the "
            "input directory when the --outdir option was not used."
            " (default: %(default)s)"
        ),
    )
    parser_f.add_argument(
        "-t",
        "--threshold",
        type=check_0to1,
        default=0.4,
        metavar="FLOAT",
        help=(
            "A cutoff threshold for filtering out stutter sequences. "
            "(default: %(default)s)"
        ),
    )
    parser_f.add_argument(
        "-n",
        "--max_alleles",
        type=int,
        default=2,
        metavar="INT",
        help=(
            "Maximum number of alleles amplified in each marker. If it differs between "
            "markers, use --path-marker-data option for set the value for each marker. "
            "If number of unique sequence variants was greater the expected number of "
            "alleles, the visual-checking procedure will be launch. "
            "(default: %(default)i)"
        ),
    )
    parser_f.add_argument(
        "-m",
        "--path-marker-data",
        default=None,
        metavar="FILE",
        help=(
            "Path to the CSV file of the marker information table which contains the "
            "fields of marker names and maximum number of alleles for each marker"
        ),
    )
    parser_f.add_argument(
        "-M",
        "--min_reads",
        type=int,
        default=6,
        help="Minimum number of reads required for proccessing (default: %(default)s)",
    )
    exclusive_group_f = parser_f.add_mutually_exclusive_group(required=False)
    exclusive_group_f.add_argument(
        "--force-visual-check",
        action="store_true",
        help="Perform visual checking for all samples",
    )
    exclusive_group_f.add_argument(
        "--force-no-visual-check",
        action="store_true",
        help=("Do not perform visual checking."),
    )
    parser_f.add_argument(
        "--n-cpu",
        type=int,
        default=1,
        metavar="INT",
        help="Number of processers (default: %(default)i)",
    )
    parser_f.add_argument(
        "-q", "--quiet", action="store_true", help="Reduce stdout messages"
    )

    # parser for the "show-alignment" command
    parser_h = subparsers.add_parser(
        "show-alignment",
        help="Show the sequence alingment",
        description="Show the sequence alingment of the input sequence file",
        usage="%(prog)s FILE [FILE, ...]",
    )
    parser_h.add_argument(
        "infile",
        metavar="INFILE",
        nargs="+",
        help="Sequences file in FASTA/FASTQ format (support gzip compression)",
    )

    if argv:
        args = parser.parse_args(argv)
    else:
        args = parser.parse_args()

    if args.subcommand == "demultiplex":
        # Check marker information file
        check_marker_data_in_args(args)

        # Check output directory
        if Path(args.outdir).exists():
            msg = "Output directory '{}' already exists. Proceed anyway?"
            try:
                if ask_user(msg.format(args.outdir), "y", False, True):
                    pass
                else:
                    sys.exit()
            except Exception:
                shutil.rmtree(args.outdir)

    elif args.subcommand == "merge-pairs":
        # Check marker information file
        if args.trim_primer:
            check_marker_data_in_args(args)

    return args


class MyFormatter(argparse.RawTextHelpFormatter):
    """
    Corrected _max_action_length for the indenting of subactions
    (https://stackoverflow.com/questions/32888815/max-help-position-\
    is-not-works-in-python-argparse-library)
    """

    def add_argument(self, action):
        if action.help is not argparse.SUPPRESS:
            # find all invocations
            get_invocation = self._format_action_invocation
            invocations = [get_invocation(action)]
            current_indent = self._current_indent
            for subaction in self._iter_indented_subactions(action):
                # compensate for the indent that will be added
                indent_chg = self._current_indent - current_indent
                added_indent = "x" * indent_chg
                invocations.append(added_indent + get_invocation(subaction))
            # print('inv', invocations)

            # update the maximum item length
            invocation_length = max([len(s) for s in invocations])
            action_length = invocation_length + self._current_indent
            self._action_max_length = max(self._action_max_length, action_length)

            # add the item to the list
            self._add_item(self._format_action, [action])


class MyFormatter2(argparse.HelpFormatter):
    """
    modify _format_action_invocation()
    show only two args_strings when nargs == "+"
    """

    def _format_action_invocation(self, action):
        if not action.option_strings:
            default = self._get_default_metavar_for_positional(action)
            (metavar,) = self._metavar_formatter(action, default)(1)
            return metavar

        else:
            parts = []

            if action.nargs == 0:
                parts.extend(action.option_strings)

            elif action.nargs == "+":
                default = self._get_default_metavar_for_optional(action)
                args_string = self._format_args(action, default)
                args_string = args_string.replace(" ...", "")
                for option_string in action.option_strings:
                    parts.append("{} {}".format(option_string, args_string))

            else:
                default = self._get_default_metavar_for_optional(action)
                args_string = self._format_args(action, default)
                for option_string in action.option_strings:
                    parts.append("{} {}".format(option_string, args_string))

            return ", ".join(parts)


class MaxTwoArgs(argparse.Action):
    def __call__(self, parser, args, values, option_string=None):
        if len(values) > 2:
            errmsg = "'{}' requires 1 or 2 argument(s)".format(self.dest)
            raise ValueError(errmsg)
        setattr(args, self.dest, values)


def check_positive(value):
    fvalue = float(value)
    if fvalue < 0:
        raise argparse.ArgumentTypeError("the value must be positive")
    return fvalue


def check_0to1(value):
    fvalue = float(value)
    if fvalue < 0 and fvalue > 1:
        raise argparse.ArgumentTypeError("the value must be between 0 to 1")
    return fvalue


def check_marker_data_in_args(args):
    if args.path_marker_data or (args.path_fwd_primers and args.path_rev_primers):
        pass
    else:
        errmsg = (
            "The information of primer sequences must be provided. "
            "Use -m or -f and -r option(s) to set primer sequences."
        )
        sys.exit(errmsg)


def print_args(args) -> None:
    if not isinstance(args, dict):
        args = vars(args)
    print("\nOPTIONS")
    for k, v in args.items():
        if isinstance(v, list):
            n_args = len(v)
            if n_args > 2 and Path(v[0]).is_file():
                v = ", ".join([str(x) for x in v[:2]])
                v += ", ... ({} files)".format(n_args)
            else:
                v = ", ".join([str(x) if isinstance(x, float) else x for x in v])
        print("{:.<25}: {}".format(k, v))
    print("")


def gen_incomplete_string(string: str) -> list[str]:
    return [string[:i] for i in range(1, len(string) + 1)]


def ask_user(
    msg: str,
    default: Literal["y", "n", "q", "o"],
    quit: bool = True,
    overwrite: bool = False,
) -> bool:
    y = gen_incomplete_string("yes")
    n = gen_incomplete_string("no")
    q = gen_incomplete_string("quit")
    o = gen_incomplete_string("overwrite")

    choice = [y, n]
    if quit:
        choice.append(q)
    if overwrite:
        choice.append(o)

    choice_init = [i[0] for i in choice]
    choice_full = [i[-1] for i in choice]
    choice_s = [i.upper() if i == default else i for i in choice_init]
    choice_s_full = ["[{}]{}".format(i, j[1:]) for i, j in zip(choice_s, choice_full)]

    for i in choice:
        if i[0] == default:
            i.append("")
            break
    else:
        errmsg = "'default' must be one of the following: {}".format(str(choice_init))
        raise ValueError(errmsg)

    prompt = "{} [{}]: ".format(msg, "/".join(choice_s))

    while True:
        try:
            input_str = input(prompt).lower()
            match = [c for c in choice if input_str in c]
            if match:
                if match[0] == y:
                    return True
                elif match[0] == n:
                    return False
                elif match[0] == q:
                    sys.exit()
                elif match[0] == o:
                    raise Exception("overwrite")
            else:
                prompt = "input {}: ".format("/".join(choice_s_full))
        except KeyboardInterrupt:
            sys.exit()


if __name__ == "__main__":
    args = get_args()
    print_args(args)
