#!/usr/bin/env python3

"""Generate genome coverage (5' positions) using bedtools genomecov."""

import argparse
import os
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description="Generate genome coverage with bedtools")
    parser.add_argument("--input", required=True, help="Input BAM file")
    parser.add_argument("--output", required=True, help="Output bedGraph coverage file")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    cmd = f"bedtools genomecov -5 -bg -ibam {args.input} > {args.output}"

    print(f"Running: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)

    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
