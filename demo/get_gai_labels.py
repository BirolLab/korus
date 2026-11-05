#!/usr/bin/env python3
"""Report the GAI label with the highest cluster fraction for each sample."""

import argparse



def get_gai_labels(input_file):
    """Yield (sample ID, GAI label) pairs from an open text file.

    Values after the first semicolon (such as p-values) are ignored.
    If fractions tie, the first cluster in the table header wins.
    """
    cluster_labels = {}
    reading_mapping = False
    cluster_names = None

    for line in input_file:
        line = line.strip()
        if not line:
            continue

        if line == "--- Cluster Label Mapping ---":
            reading_mapping = True
            continue

        if reading_mapping:
            if set(line) == {"-"}:
                reading_mapping = False
                continue
            cluster, label = line.split(":")

            cluster_labels[cluster.strip()] = label.strip()
            continue

        fields = line.split("\t")
        if fields[0] == "ID":
            cluster_names = fields[1:]
            continue

        fractions = [float(value.split(";")[0]) for value in fields[1:]]

        highest_index = max(range(len(fractions)), key=fractions.__getitem__)
        yield fields[0], cluster_labels[cluster_names[highest_index]]



def main():
    """Parse a results file and write sample IDs and GAI labels to stdout."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_file", help="Tab-delimited cluster results file")
    args = parser.parse_args()

    with open(args.input_file, 'r', encoding="utf-8") as input_file:
        for id, gai in get_gai_labels(input_file):
            print(f"{id}\t{gai}")



if __name__ == "__main__":
    main()
