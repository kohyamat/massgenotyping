from massgenotyping.argument_parser import get_args

__version__ = "0.1.1"


def main():
    args = get_args()
    if args.subcommand == "demultiplex":
        from massgenotyping import demultiplex

        demultiplex.main(args)

    elif args.subcommand == "merge-pairs":
        from massgenotyping import merge_pairs

        merge_pairs.main(args)

    elif args.subcommand == "denoise":
        from massgenotyping import denoise

        denoise.main(args)

    else:
        from matplotlib import get_backend

        if get_backend() == "agg":
            print(
                (
                    "\033[33mWARNING: No display found. "
                    "Skip visual checking procedures\033[0m"
                )
            )

        if args.subcommand == "filter":
            from massgenotyping import variant_filter

            variant_filter.main(args)

        elif args.subcommand == "allele-check":
            from massgenotyping import allele_check

            allele_check.main(args)

        elif args.subcommand == "allele-call":
            from massgenotyping import allele_call

            allele_call.main(args)

        elif args.subcommand == "show-alignment":
            from massgenotyping import show_alignment

            show_alignment.main(args)
