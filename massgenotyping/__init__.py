from .argument_parser import get_args

__version__ = "0.1.0"


def main():
    args = get_args()
    if args.subcommand == "demultiplex":
        from . import demultiplex

        demultiplex.main(args)

    elif args.subcommand == "merge-pairs":
        from . import merge_pairs

        merge_pairs.main(args)

    elif args.subcommand == "denoise":
        from . import denoise

        denoise.main(args)

    elif args.subcommand == "filter":
        from . import variant_filter

        variant_filter.main(args)

    elif args.subcommand == "allele-check":
        from . import allele_check

        allele_check.main(args)

    elif args.subcommand == "allele-call":
        from . import allele_call

        allele_call.main(args)

    elif args.subcommand == "show-alignment":
        from . import show_alignment

        show_alignment.main(args)
