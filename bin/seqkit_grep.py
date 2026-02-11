#!/usr/bin/env python3

"""Filter FASTQ reads by transposon sequence using seqkit grep."""

import argparse
import os
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description="Filter reads by transposon sequence")
    parser.add_argument("--input", required=True, help="Input FASTQ file (.fq.gz)")
    parser.add_argument("--output", required=True, help="Output filtered FASTQ file")
    parser.add_argument("--pattern", required=True, help="Transposon static region sequence")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    # zcat | seqkit grep
    cmd = (
        f"zcat < {args.input} | "
        f"seqkit grep --threads {args.threads} --by-seq --ignore-case "
        f"--pattern {args.pattern} -o {args.output}"
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
