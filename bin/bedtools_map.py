#!/usr/bin/env python3

"""Map read counts to gene features using bedtools map."""

import argparse
import os
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description="Map counts to gene features with bedtools")
    parser.add_argument("--input", required=True, help="Input bedGraph coverage file")
    parser.add_argument("--reference", required=True, help="Reference BED file with gene features")
    parser.add_argument("--output", required=True, help="Output mapped counts file")
    parser.add_argument("--operation", required=True, choices=["sum", "count"],
                        help="Bedtools map operation: sum (total counts) or count (unique positions)")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    cmd = (
        f"bedtools map -a {args.reference} -b {args.input} "
        f"-c 4 -o {args.operation} > {args.output}"
    )

    print(f"Running: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)

    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
