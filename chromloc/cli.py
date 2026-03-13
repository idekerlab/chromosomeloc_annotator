import argparse
import sys
from typing import Optional, Sequence
import ndex2.client
from .annotate import Config, run_update


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Annotate HCX hierarchy nodes with chromosome counts and pie style."
        )
    )
    parser.add_argument(
        dest="input",
        default="-",
        help=(
            "Input CX2 file."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        default="-",
        help=(
            "Output JSON file for updated network. Use '-' (default) to write to stdout."
        ),
    )
    parser.add_argument(
        "--chrom-map",
        dest="chrom_map",
        default=None,
        help="Path to hgnc gene→chromosome JSON/TSV map (default: packaged map if present).",
    )

    parser.add_argument(
        "--interaction-uuid-attr",
        dest="interaction_uuid_attr",
        default="HCX::interactionNetworkUUID",
        help="Node attribute that stores the interaction network NDEx UUID.",
    )
    parser.add_argument(
        "--interaction-uuid-network-attr",
        dest="interaction_uuid_network_attr",
        default="HCX::interactionNetworkUUID",
        help="Network attribute fallback for interaction network UUID.",
    )
    parser.add_argument(
        "--gene-list-delim",
        dest="gene_list_delim",
        default=",",
        help="Delimiter used when gene attributes contain multiple genes in a string.",
    )
    parser.add_argument(
        "--species",
        dest="species",
        default="human",
        help="Species code (default human). Controls chromosome list naming.",
    )
    parser.add_argument(
        "--ndex-server",
        dest="ndex_server",
        default="https://public.ndexbio.org",
        help="NDEx server URL for fetching interaction networks.",
    )
    parser.add_argument(
        "--ndex-username",
        dest="ndex_username",
        default=None,
        help="NDEx username (optional, for private networks).",
    )
    parser.add_argument(
        "--ndex-password",
        dest="ndex_password",
        default=None,
        help="NDEx password (optional, for private networks).",
    )

    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """
    Console script entry point.

    This is wired up in setup.py as:

    """
    if argv is None:
        argv = sys.argv[1:]

    parser = _build_arg_parser()
    args = parser.parse_args(argv)

    config = Config(
        client=ndex2.client.Ndex2(args.ndex_server, args.ndex_username, args.ndex_password),
        chrom_map_path=args.chrom_map,
        interaction_uuid_attr=args.interaction_uuid_attr,
        interaction_uuid_network_attr=args.interaction_uuid_network_attr,
        gene_list_delim=args.gene_list_delim,
        species=args.species,
        ndex_server=args.ndex_server,
        ndex_username=args.ndex_username,
        ndex_password=args.ndex_password,
    )

    # Determine output stream
    if args.output == "-":
        out_stream = sys.stdout
    else:
        out_stream = open(args.output, "w", encoding="utf-8")

    try:
        run_update(
            input_source=args.input,
            output_stream=out_stream,
            config=config,
        )
    finally:
        if args.output != "-":
            out_stream.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
