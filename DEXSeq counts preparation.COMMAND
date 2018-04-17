# The addresses should be returned in the "Environment" panel
# Set the working directory in Terminal to the folder that contains your .gtf and .bam files
# The following command entered into the terminal turns a genome annotation file from a .gtf to a .gff file (mind the spaces):
python /~returned_address~/library/DEXSeq/python_scripts/dexseq_prepare_annotation.py name_of.gtf desired_name_of.gff

# Use samtools or pysam to sort .bam files if this hasn't been done already.
# Next, call the following python module to get the .bam data aligned to exons in the .gff file:

python /~returned_address~/library/DEXSeq/python_scripts/dexseq_count.py your.gff -p yes -r pos -s no -f bam file.bam desired_name.txt
