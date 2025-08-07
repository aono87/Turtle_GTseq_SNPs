#!/bin/bash
# START WITH THIS FILE IN /home/$user/gtseq/$sp/$run/raw_fastq/$project/


# run in directory containing fastq.gz 
#**# indicates that the information needs to be edited for the code to work properly

#**# add a run number example : ms30
run=ms48

#**# add a species identifyer 
sp=dcor

#**# add a species project, this will be used for the directories use on Ahi
project=stx_hatch
plate=
#**# add user information, this should be your ahi user name, again used for directories used on ahi
user=afrey

#add a reference file name
ref=/home/$user/gtseq/ref/DcGTseqPanel_205_112122.fasta
ref2=/home/$user/gtseq/ref/dc_msat_primerF.txt
ref3=/home/$user/gtseq/ref/dc205_targetregions.bed

LIST=/home/$user/gtseq/ref/fprim.list

#add run directory
loc=/home/$user/gtseq/$sp/$run/raw_fastq/$project/
loc2=/home/$user/gtseq/$sp/$run/megasat/msat_fastq/$project/
loc3=/home/$user/gtseq/$sp/$run/trimmed_fastq/$project/
loc4=/home/$user/gtseq/$sp/$run/bamfiles/$project/
loc5=/home/$user/gtseq/$sp/$run/metadata/trim_stat_files/$project/
loc6=/home/$user/gtseq/$sp/$run/metadata/run/$project/
loc7=/home/$user/gtseq/$sp/$run/bamfiles/$project/depth/
loc8=/home/$user/gtseq/$sp/$run/raw_fastq/$project/unmapped/

# index the reference fasta file
bwa index $ref
samtools faidx $ref
# this is really fast for a small file, less than a second


for file in *.fastq.gz
# files are labeled as labID_sample#_L001_R1_001.fastq.gz

do

echo $file

#extracts the first section (labid) of the filename separated by "_"
sample=`echo $file | cut -f1 -d "."`
echo $sample
#counts the reads in the orginal fastq.gz file and outputs to a text file
#echo "read count"
#zcat "$sample"_*.fastq.gz | echo $((`wc -l`/4))
#zcat "$sample"_*.fastq.gz  | echo $sample $((`wc -l`/4)) >>  "$loc6"readcounts_"$run"_"$project".txt
#renames the file to eliminate everything but the z#
#mv ./"$sample"_*.fastq.gz "$sample".fastq.gz


#exec < $LIST
	
#	while read LINE
#	do
#		zgrep -o "$line" ./$sample.fastq.gz | wc -l >> $loc7"$sample"_fprim_ct.txt
	
#	done

#trims the adaptors from the fastq.gz
bbduk.sh in=$sample.fastq.gz ref=/home/$user/gtseq/ref/adaptorGT.fasta k=12 stats=$loc5"$sample"_stats.txt out=$loc3"$sample"_trimmed.fastq.gz outm=$loc8"$sample"_trash.fastq.gz
#removes the original fastq.gz file

echo $sample
#counts the reads in the trimmed fastq.gz file and outputs to a text file
#echo "trimmed read count"
#zcat $loc3"$sample"_trimmed.fastq.gz | echo $((`wc -l`/4))
#zcat $loc3"$sample"_trimmed.fastq.gz | echo $sample $((`wc -l`/4)) >> "$loc6"trimmedreadcounts_"$run"_"$project".txt
# generates the sam file
#bwa mem -M -v 3 -t 15 -R "@RG\tID:$sample\tLB:amplicon\tPL:ILLUMINA\tSM:"$sp$sample"_"$project"" $ref $loc3"$sample"_trimmed.fastq.gz > $loc4"$sample"_mapped2ref.sam
# generates the bam file and sorts it
#/usr/local/samtools-1.2/bin/samtools view -bh -@ 6 $loc4"$sample"_mapped2ref.sam | /usr/local/samtools-1.2/bin/samtools sort -@ 6 - $loc4"$sample"_mapped2ref_sorted
#index bam
#/usr/local/samtools-1.2/bin/samtools index $loc4"$sample"_mapped2ref_sorted.bam
#counts mapped reads and outputs to a text file
#echo $sample
#echo "mapped read count"
#samtools view -F 0x904 -c $loc4"$sample"_mapped2ref_sorted.bam
#printf "`echo $sample` `samtools view -F 0x904 -c $loc4"$sample"_mapped2ref_sorted.bam`\n" >>"$loc6"mapped2ref_"$run"_"$project".txt
#gets locus read depth
samtools depth -b /home/afrey/gtseq/ref/dc205_targetregions.bed $loc4"$sample"_mapped2ref_sorted.bam > $loc7"$sample".coverage

#extract msat reads
#/usr/local/samtools-1.2/bin/samtools view -f 4 $loc4"$sample"_mapped2ref_sorted.bam > $loc4"$sample"_unmapped.sam
#cat $loc4"$sample"_unmapped.sam | fgrep -f $ref2 > $loc4"$sample"_msat.sam
#cat $loc4"$sample"_msat.sam | grep -v @ | awk '{print "@"$1"\n"$10"\n+\n"$11}' > $loc2"$sample"_msat.fastq
#cat $loc2"$sample"_msat.fastq | echo $((`wc -l`/4))
#cat $loc2"$sample"_msat.fastq | echo $sample $((`wc -l`/4)) >> "$loc6"Msat_readcounts_"$run"_"$project".txt

#creating unmapped.fastq
#cat $loc4"$sample"_unmapped.sam | grep -v @ | awk '{print "@"$1"\n"$10"\n+\n"$11}' > $loc8"$sample"_unmapped.fastq

done


# tarball up the sams for easy transfer (for use in microhaplotype R package).

#tar -cvf - "$loc4"*_mapped2ref.sam | gzip -9 > $loc6"$sp"_"$run"_"$project"_sam_output.tar.gz


