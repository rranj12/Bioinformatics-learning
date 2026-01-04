# Genomic Data Formats & File Handling

## Summary
- Downloaded E. coli reference genome (NC_000913.3) and parsed FASTA
- Computed length, GC%, nucleotide counts, quality plot
  
## Learnigs
- FASTA stores sequences; FASTQ stores sequences + per-base quality scores
- Phred scores map to error probabilities 
- Quality tends to degrade toward the end of reads (platform-dependent)


# Gene Finding & ORF Analysis

## Summary
- Wrote a custom ORF finder in Python using Biopython
- Parsed a full bacterial genome (E. coli, NC_000913.3)
- Identified open reading frames across all six reading frames
- Filtered ORFs by minimum length and exported results to CSV
- Generated basic visualizations of ORF length distributions
  
## ORF-Finding Script 
- Scans both forward and reverse-complement strands
- Considers all three reading frames per strand and detects start/stop codons
- Tracks genomic coordinates correctly (including reverse-strand mapping)
- Filters ORFs based on minimum length
- outputs structured results and a histogram of ORF lengths


# Sequence Alignment & BLAST
I used BLAST to perform a comparative analysis between two related bacterial genomes.

## Summary: 
- Learned core sequence alignment concepts (local vs global, heuristic vs exact)
- Set up and ran BLAST locally
- Downloaded and compared closely related bacterial genomes
- Parsed BLAST output into a structured format for analysis
- Explored conserved regions between E. coli and Shigella

## BLAST Alignment Steps:
- Downloaded E. coli and Shigella reference genomes from NCBI
- Built a local BLAST nucleotide database
- Ran blastn to align sequences
- Exported results in tabular format
- Analyzed percent identity, alignment length, and significance in a notebook


# Alignment Theory & RNA-seq with STAR

## Summary
- Studied alignment theory in more depth
- Built a STAR reference genome index
- Aligned real RNA-seq reads from a public dataset
- Generated gene-level expression counts and QC metrics

## Pipeline steps
- Downloaded RNA-seq reads from SRA
- Downloaded reference genome (chr 22) + GTF annotations
- Built a STAR genome index
- Aligned reads using STAR
- Inspected alignment logs and splice junction outputs
- Generated a gene expression matrix

