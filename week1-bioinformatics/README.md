# Bioinformatics Day 1: Genomic Data Formats & File Handling 
## Summary: 
- Set up 'compbio' conda environment with core bioinformatics tools 
- Downloaded E.Cole reference genome (NC_00913.3) and parsed FASTA
- Computed length, GC%, nucleotide counts
- Built a FASTQ quality analysis plot/analysis (per-position mean)

## Takeaways
- FASTA stores sequences; FASTQ stores sequences + per-base Phred quality scores 
- Phred scores map to error probabilities 
- Quality tends to degrade toward end of reads 

## Lingering Questions: 
- how to reliably detect Phred 33 vs Phred 64 in the wild? 

