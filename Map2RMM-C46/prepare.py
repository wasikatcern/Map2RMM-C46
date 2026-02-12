#!/usr/bin/env python3
import ROOT
import numpy as np
import csv
import sys

def main():
    if len(sys.argv) < 3:
        print("Usage:  python prepare.py  input.root  output.csv")
        sys.exit(1)

    infile  = sys.argv[1]
    outfile = sys.argv[2]

    print(f"[INFO] Opening ROOT file: {infile}")

    f = ROOT.TFile(infile)
    if not f or f.IsZombie():
        print(f"[ERROR] Cannot open file {infile}")
        sys.exit(1)

    t = f.Get("inputNN")
    if not t:
        print("[ERROR] TTree 'inputNN' not found in file.")
        sys.exit(1)

    # Verify branches
    if not hasattr(t, "c46"):
        print("[ERROR] Branch 'c46' not found. This ROOT file is not RMM-C46 format.")
        sys.exit(1)

    print(f"[INFO] Events in tree: {t.GetEntries()}")

    # Build CSV header
    header = ["id"] + [f"C46_{i+1}" for i in range(46)]

    print(f"[INFO] Writing CSV: {outfile}")
    with open(outfile, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)

        for i, event in enumerate(t):
            c46_vec = list(event.c46)
            if len(c46_vec) != 46:
                print(f"[WARNING] event {i}: c46 vector has size {len(c46_vec)} (expected 46)")

            row = [int(event.id)] + c46_vec
            writer.writerow(row)

            if (i+1) % 5000 == 0:
                print(f"[INFO] Processed {i+1} events...")

    print(f"[INFO] Done. CSV saved to {outfile}")


if __name__ == "__main__":
    main()

