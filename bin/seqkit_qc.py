#!/usr/bin/env python3

"""Generate quality control statistics using seqkit."""

import argparse
import os
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description="QC statistics with seqkit")
    parser.add_argument("--mode", required=True, choices=["fastq", "bam"],
                        help="QC mode: fastq (seqkit stats) or bam (seqkit bam)")
    parser.add_argument("--output", required=True, help="Output QC report file")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    parser.add_argument("--input", action="append", required=True,
                        help="Input file(s). Can be specified multiple times.")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    input_files = " ".join(args.input)

    if args.mode == "fastq":
        cmd = f"seqkit stats -j {args.threads} {input_files} --quiet -a -T > {args.output} 2>&1"
    else:
        cmd = f"seqkit bam -j {args.threads} -s {input_files} --quiet > {args.output} 2>&1"

    print(f"Running: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)

    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
