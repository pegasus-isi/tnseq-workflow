#!/usr/bin/env python3

"""Generate individual sample count tab files using tab.R."""

import argparse
import os
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description="Generate tab files from mapped data")
    parser.add_argument("--input", required=True, help="Input mapped counts file")
    parser.add_argument("--output", required=True, help="Output tab file")
    args = parser.parse_args()

    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # Pegasus stages tab.R into the job's working directory
    tab_r = os.path.join(os.getcwd(), "tab.R")

    cmd = ["Rscript", tab_r, "--data", args.input, "--out", args.output]

    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)

    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
