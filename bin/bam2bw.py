#!/usr/bin/env python3

"""Generate BigWig coverage tracks from BAM files."""

import argparse
import os
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description="Generate BigWig from BAM")
    parser.add_argument("--input", required=True, help="Input BAM file")
    parser.add_argument("--output", required=True, help="Output BigWig file")
    parser.add_argument("--threads", type=int, default=2, help="Number of threads")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    # Step 1: Index the BAM file
    index_cmd = f"samtools index -@ {args.threads} {args.input}"
    print(f"Running: {index_cmd}")
    result = subprocess.run(index_cmd, shell=True, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)

    # Step 2: Generate BigWig
    bw_cmd = f"bamCoverage -b {args.input} -o {args.output}"
    print(f"Running: {bw_cmd}")
    result = subprocess.run(bw_cmd, shell=True, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)

    print(f"Output: {args.output}")
    print(f"BAM index: {args.input}.bai")


if __name__ == "__main__":
    main()
