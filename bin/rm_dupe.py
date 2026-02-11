#!/usr/bin/env python3

"""Remove PCR duplicates using Je markdupes with UMI awareness."""

import argparse
import os
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description="Remove PCR duplicates using Je")
    parser.add_argument("--input", required=True, help="Input BAM file")
    parser.add_argument("--output", required=True, help="Output deduplicated BAM file")
    parser.add_argument("--metrics", required=True, help="Output deduplication metrics file")
    parser.add_argument("--je-jar", required=True, help="Path to je_1.2_bundle.jar")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    os.makedirs(os.path.dirname(args.metrics), exist_ok=True)

    cmd = [
        "java", "-Xmx5500m", "-jar", args.je_jar,
        "markdupes",
        f"I={args.input}",
        f"O={args.output}",
        f"M={args.metrics}",
        "MM=0",
        "REMOVE_DUPLICATES=TRUE",
    ]

    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)

    print(f"Output: {args.output}")
    print(f"Metrics: {args.metrics}")


if __name__ == "__main__":
    main()
