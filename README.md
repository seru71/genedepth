# genedepth

Usage:

 1) Download NCBIrefSeq transcript table to `refseq.txt`
	
 2) Remove non-standard chromosomes and sort by chr and position, e.g.
 
    `grep -vP 'chr[0-9a-zA-Z]+_' refseq.txt | sort -k3,3V -k5,5n -k6,6n > refseq.sorted.txt`
    
 3) Extract exons and sort by coordinate
 
	`cat refseq.sorted.txt | ./genedepth.py | bedtools sort -i - > refseq.exons.bed`
	
 4) Calculate depth
 
	`samtools depth -a BAM | ./gene_depth.py refseq.exons.bed OUTPUT_PREFIX`
	
 5) Mean, median, and fraction of bp covered by 10, 20 and 30x will be provided in two files:
	`OUTPUT_PREFIX.exon_stats`
	`OUTPUT_PREFIX.transcript_stats`
	

 TODOs
	Stats for transcripts present on two chromosomes (e.g. X and Y) are not separated and may be wrong
