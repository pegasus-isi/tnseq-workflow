#!/usr/bin/env python3

"""Align reads to reference genome using BWA-MEM and sort with samtools."""

import argparse
import os
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description="Align reads with BWA-MEM")
    parser.add_argument("--input", required=True, help="Input FASTQ file (.fq.gz)")
    parser.add_argument("--output", required=True, help="Output sorted BAM file")
    parser.add_argument("--reference", required=True, help="Reference FASTA file")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    cmd = (
        f"bwa mem -t {args.threads} {args.reference} {args.input} | "
        f"samtools sort -@2 - > {args.output}"
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
