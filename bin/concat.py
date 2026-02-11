#!/usr/bin/env python3

"""Concatenate all sample tab files into a single TSV using concat.R."""

import argparse
import os
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description="Concatenate sample tab files")
    parser.add_argument("--reference", required=True, help="Reference BED file")
    parser.add_argument("--input", action="append", required=True,
                        help="Input .tab file(s). Can be specified multiple times.")
    parser.add_argument("--output", required=True, help="Output TSV file")
    args = parser.parse_args()

    # Pegasus stages concat.R into the job's working directory
    concat_r = os.path.join(os.getcwd(), "concat.R")

    # concat.R takes: ref_bed output_file tab_file1 tab_file2 ...
    cmd = ["Rscript", concat_r, args.reference, args.output] + args.input

    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)

    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
