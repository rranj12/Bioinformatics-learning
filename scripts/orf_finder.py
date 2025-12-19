from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import csv
import os
import matplotlib.pyplot as plt

STOP_CODONS = {"TAA", "TAG", "TGA"}
START_CODON = "ATG"


def find_orfs(seq_str: str, min_length: int, strand: str):
    s = seq_str.upper()
    orfs = []

    for frame in range(3):
        start_pos = None  # start of current ORF (first ATG after last stop)
        i = frame

        while i + 3 <= len(s):
            codon = s[i:i+3]

            if codon in STOP_CODONS:
                # if we had a start, close ORF at this stop
                if start_pos is not None:
                    end = i + 3
                    length = end - start_pos
                    if length >= min_length:
                        orfs.append({
                            "start": start_pos,
                            "end": end,
                            "length": length,
                            "frame": frame,
                            "strand": strand,
                            "dna": s[start_pos:end],
                        })
                # reset after stop codon
                start_pos = None

            elif codon == START_CODON:
                # only take the FIRST start in an open segment
                if start_pos is None:
                    start_pos = i

            i += 3

    return orfs


def translate_dna(dna: str, table: int = 11):
    return str(Seq(dna).translate(table=table, to_stop=True))

def map_rev_coords(orfs, genome_len: int):
    mapped = []
    for o in orfs:
        m = dict(o)
        m["start"], m["end"] = genome_len - o["end"], genome_len - o["start"]
        mapped.append(m)
    return mapped

def main():
    parser = argparse.ArgumentParser(description="ORF Finder")
    parser.add_argument("--fasta", required=True, help="Input FASTA file")
    parser.add_argument("--min-length", type=int, default=100, help="Minimum ORF length (bp)")
    parser.add_argument("--out", default="results/orfs", help="Output prefix (default: results/orfs)")
    parser.add_argument("--table", type=int, default=11, help="Translation table (default 11 for bacteria/E. coli)")
    args = parser.parse_args()

    record = next(SeqIO.parse(args.fasta, "fasta"))
    seq = record.seq
    genome_len = len(seq)

    # forward
    plus = find_orfs(str(seq), args.min_length, strand="+")
    for o in plus:
        o["protein"] = translate_dna(o["dna"], table=args.table)

    # reverse complement
    rc = seq.reverse_complement()
    minus = find_orfs(str(rc), args.min_length, strand="-")
    minus = map_rev_coords(minus, genome_len)
    for o in minus:
        # Need DNA sequence in original orientation for translation?
        # easiest: translate from the reverse-complement ORF DNA itself
        # We'll reconstruct by taking forward genome slice and revcomp it.
        dna_forward_slice = seq[o["start"]:o["end"]]          # slice on forward genome
        dna_on_minus = dna_forward_slice.reverse_complement() # actual coding DNA on minus strand
        o["dna"] = str(dna_on_minus).upper()
        o["protein"] = translate_dna(o["dna"], table=args.table)

    orfs = plus + minus
    orfs.sort(key=lambda x: (x["start"], x["strand"], x["frame"]))

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    csv_path = args.out + ".csv"
    png_path = args.out + "_length_hist.png"

    # write CSV
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["start", "end", "length", "frame", "strand", "protein"])
        for o in orfs:
            w.writerow([o["start"], o["end"], o["length"], o["frame"], o["strand"], o["protein"]])

    # histogram
    lengths = [o["length"] for o in orfs]
    plt.figure()
    plt.hist(lengths, bins=50)
    plt.xlabel("ORF length (bp)")
    plt.ylabel("Count")
    plt.title(f"ORF length distribution (min={args.min_length})")
    plt.tight_layout()
    plt.savefig(png_path, dpi=200)
    plt.close()

    print(f"Loaded {record.id}")
    print(f"Genome length: {genome_len} bp")
    print(f"Found {len(orfs)} ORFs (min_length={args.min_length})")
    print(f"Wrote {csv_path}")
    print(f"Wrote {png_path}")

if __name__ == "__main__":
    main()
