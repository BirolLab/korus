#!/usr/bin/env python3
"""Merge single-sample KORUS matrices."""

import argparse
import math
from pathlib import Path
import sys


def parse_arguments():
    """Parse matrix paths and output options."""
    parser = argparse.ArgumentParser(
        description=(
            "Merge single-sample matrices with matching window sizes. "
            "Preserve input sample order and first-seen coordinate order; "
            "fill missing values with NA."
        ),
        epilog=(
            "Examples:\n"
            "  python3 aggregateMatrices.py sample1.tsv sample2.tsv -o merged.tsv\n"
            "  python3 aggregateMatrices.py --matrix-list matrices.txt -o merged.tsv\n\n"
            "Each matrix must have a header: chr window_start:<size> <sample>.\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("matrices", nargs="*", type=Path, help="Input matrix paths")
    parser.add_argument(
        "--matrix-list", type=Path, metavar="FILE",
        help="File with input matrices (one per line)",
    )
    parser.add_argument(
        "-o", "--output", default="-", metavar="FILE",
        help="Output TSV path (default: standard output; '-' also means stdout)",
    )
        
    return parser, parser.parse_args()


def collect_paths(args):
    """Collect input paths, resolving list entries relative to the list file."""
    paths = list(args.matrices)
    if args.matrix_list is not None:
        with args.matrix_list.open('r', encoding="utf-8") as handle:
            for line in handle:
                entry = line.strip()
                if entry:
                    path = Path(entry)
                    if not path.is_absolute():
                        path = args.matrix_list.parent / path
                    paths.append(path)
    if not paths:
        raise ValueError("Provide at least one matrix path or a nonempty --matrix-list.")
    return paths


def read_header(handle, path):
    """Validate a matrix header and return its window size and sample name."""
    header = handle.readline().strip().split("\t")
    if (
        len(header) != 3
        or header[0] != "chr"
        or not header[1].startswith("window_start:")
    ):
        raise ValueError(
            f"{path}:1: Expected header 'chr window_start:<size> <sample>'."
        )
    size = header[1].partition(":")[2]
    return int(size), header[2]


def read_rows(handle, path, sample_index, rows):
    """Validate rows and add values to the merged coordinate dictionary."""
    seen = set()
    for line_number, line in enumerate(handle, start=2):
        fields = line.strip().split("\t")
        if len(fields) != 3:
            raise ValueError(f"{path}:{line_number}: Expected exactly three columns.")
        chromosome, window, value = fields
        coordinate = (chromosome, str(int(window)))
        if coordinate in seen:
            raise ValueError(
                f"{path}:{line_number}: Duplicate coordinate {chromosome}, {window}."
            )
        seen.add(coordinate)
        rows.setdefault(coordinate, {})[sample_index] = value
    if not seen:
        raise ValueError(f"{path}: Matrix contains no data rows.")


def merge_matrices(paths):
    """Read matrices and require distinct samples and matching sizes."""
    samples = []
    seen_samples = set()
    rows = {}
    window_size = None
    for path in paths:
        with path.open('r', encoding="utf-8") as handle:
            size, sample = read_header(handle, path)
            if window_size is not None and size != window_size:
                raise ValueError(
                    f"{path}:1: Window size {size} does not match {window_size} "
                    f"from {paths[0]}."
                )
            if sample in seen_samples:
                raise ValueError(f"{path}:1: Duplicate sample name '{sample}'.")
            window_size = size
            read_rows(handle, path, len(samples), rows)
            samples.append(sample)
            seen_samples.add(sample)
    return window_size, samples, rows


def write_matrix(handle, window_size, samples, rows):
    """Write a merged TSV, including NA for absent coordinates."""
    handle.write("\t".join(["chr", f"window_start:{window_size}", *samples]) + "\n")
    for coordinate, values in rows.items():
        fields = [*coordinate, *(values.get(index, "NA") for index in range(len(samples)))]
        handle.write("\t".join(fields) + "\n")


def main():
    """Validate all inputs before opening the output and report useful errors."""
    parser, args = parse_arguments()
    
    if not args.matrices and not args.matrix_list:
        print("Please either supply matrix files names as positional arguments or using `--matrix-list`")
        parser.print_help()
        parser.exit()
    
    paths = collect_paths(args)
    window_size, samples, rows = merge_matrices(paths)
    if args.output == "-":
        write_matrix(sys.stdout, window_size, samples, rows)
    else:
        with Path(args.output).open("w", encoding="utf-8") as handle:
            write_matrix(handle, window_size, samples, rows)



if __name__ == "__main__":
    main()
