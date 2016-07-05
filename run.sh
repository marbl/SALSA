#!/bin/bash

# $1 : bam file
# $2 : assembled contigs

echo $1
echo $2
fai='.fai'
echo $2$fai
#first get contig lengths from assembled contigs
echo 'finding contig lengths....'
cleaned_assembly='cleaned.fa'
samtools faidx $2

cut -f 1,2 $2$fai > contig_lengths 

#convert alignment bam file to bedfile
echo 'converting from bam to bed format....'
bamToBed -i $1 > alignment_unique.bed
#first break contigs using .cpp file
echo 'finding misassemblies using mapping of reads to contigs ......'
./break_contigs

#write broken contigs to a new fasta file using pythonfile
echo 'writing corrected assembly......'
python break_contigs.py -b breakpoints -a $2 -o $cleaned_assembly

#now create links
echo 'generating links between contigs.....'
python create_links.py > new_links
#sort file containgng scores
sort -k 3 -g -r new_links > new_links_sorted

#use these score to links_to_graph.py and get scaffolds
echo 'generating and writing out scaffolds....'
python links_to_graph.py -c $cleaned_assembly