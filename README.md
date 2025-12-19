# Day 1: Genomic Data Formats & File Handling

## Summary
- Set up `compbio` conda environment with core bioinformatics tools
- Downloaded E. coli reference genome (NC_000913.3) and parsed FASTA
- Computed length, GC%, nucleotide counts, quality plot
  
## Learnigs
- FASTA stores sequences; FASTQ stores sequences + per-base quality scores
- Phred scores map to error probabilities 
- Quality tends to degrade toward the end of reads (platform-dependent)

## Lingering questions
1. How do tools reliably detect Phred33 vs Phred64 in the wild?

