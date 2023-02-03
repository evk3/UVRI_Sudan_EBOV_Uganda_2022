#!/bin/bash -l
## Grid Engine Example Job Script  

# -- Begin SGE embedded arguments --
#$ -V
#Pass all environment variables to job
#$ -cwd
#Use current working directory

#$ -N Pakistan
# Name of script

#$ -j y
#Combine standard error and output files.

#$-q short.q
#Use the short.q queue, and not any other queue.

#$-pe smp 2
#Ask for a parallel environment for multi-threading.

#bad/funky nodes
#$-l h=!'node228.hpc.biotech.cdc.gov'&!'node229.hpc.biotech.cdc.gov'&!'node230.hpc.biotech.cdc.gov'

# -- End SGE embedded arguments --

module load bwa/0.7.17
module load samtools/1.9
module load picard/2.21.1
module load BEDTools/2.27.1
module load cutadapt/2.3
module load prinseq/0.20.3
module load htslib/1.9
module load bowtie2/2.3.5.1
module load gatk/4.1.7.0
module load htslib/1.10
module load gcc/9.2.0

# create temp directory for work on /scratch

mkdir -p /scicomp/scratch/evk3/EBOV/mapping/
scratch='/scicomp/scratch/evk3/EBOV/mapping/'

echo "Hostname:" $HOSTNAME
echo "SGE Value " $SGE_TASK_ID

OUTPUT_PATH=/scicomp/home/evk3/Diagnostics/2022_Ebola_Sudan_Uganda/2022003623_3624_4105/
echo "Output path: " $OUTPUT_PATH


REFERENCE=./Nakisamata_reference.fasta
cp $REFERENCE $scratch/

echo "Reference: " $REFERENCE

SEEDFILE=./file_names.txt
file_num=$(awk "NR==$SGE_TASK_ID" $SEEDFILE)

echo "First File" $file_num

sample_num=$(sed -r 's/(.*)_.*_.*_.*_.*$/\1/g' <<< "$file_num")
echo "Sample number is: " $sample_num

L1_READ1=$file_num
L1_READ2=$(awk "NR==($SGE_TASK_ID + 1)" $SEEDFILE)

echo $L1_READ1
echo $L1_READ2


echo "Gunzipping now!"
gunzip -c $L1_READ1 > "$scratch"/"$sample_num"_R1_cutadapt.fastq
gunzip -c $L1_READ2 > "$scratch"/"$sample_num"_R2_cutadapt.fastq

#Remove low quality reads:
echo "starting printseq-lite"
prinseq-lite -fastq "$scratch"/"$sample_num"_R1_cutadapt.fastq -fastq2 "$scratch"/"$sample_num"_R2_cutadapt.fastq -min_qual_mean 25 -trim_qual_right 20 -min_len 50 -out_good "$scratch"/"$sample_num"_trimmed

echo "Indexing reference sequence using bowtie2"
bwa index "$scratch"/"$REFERENCE" -p "$scratch"/"$REFERENCE"


echo "Mapping reads to reference genome"
bwa mem -t $NSLOTS "$scratch"/"$REFERENCE" "$scratch"/"$sample_num"_trimmed_1.fastq "$scratch"/"$sample_num"_trimmed_2.fastq > "$scratch"/"$sample_num"_reads.sam


echo "Starting samtools - convert SAM to BAM"
samtools view -S -b -o "$scratch"/"$sample_num"_reads.bam "$scratch"/"$sample_num"_reads.sam


# Previously generated samtools index of reference genome.  Generates *.fai file and only need to do 1X.
# This step may not be necessary with the samtool pipeline below:
echo "Indexing Ebo genome with samtools"
samtools faidx "$scratch"/"$REFERENCE"

echo "Starting samtools sort BAM file"
samtools sort -@ $NSLOTS "$scratch"/"$sample_num"_reads.bam -o "$scratch"/"$sample_num"_reads.sorted.bam

echo "Starting samtools index BAM file"
samtools index "$scratch"/"$sample_num"_reads.sorted.bam


JAVA_OPTS='-Xmx50g'
TMP_DIR=/tmp


#Copy only mapped bases, and save:
samtools view -b -F 4 "$scratch"/"$sample_num"_reads.sorted.bam > "$scratch"/"$sample_num"_reads-mapped.sorted.bam


#Make intermediate1 fasta:
echo "Making intermediate1 fasta!"
bcftools mpileup -A -d 6000000 -B -Q 0 -Ov -f "$scratch"/"$REFERENCE" "$scratch"/"$sample_num"_reads-mapped.sorted.bam | bcftools call --ploidy 1 -mv -Ov -o "$scratch"/"$sample_num"_reads-mapped.sorted.bam.vcf

gatk IndexFeatureFile --input "$scratch"/"$sample_num"_reads-mapped.sorted.bam.vcf

gatk CreateSequenceDictionary -R "$scratch"/"$REFERENCE"

gatk FastaAlternateReferenceMaker -R "$scratch"/"$REFERENCE" -O "$scratch"/"$sample_num"_intermediate1_no_contigs.fasta -V "$scratch"/"$sample_num"_reads-mapped.sorted.bam.vcf

#**********************************************************************************************************

echo "Performing mapping to Self #1"

echo "Indexing reference sequence using bowtie2"
bowtie2-build-s "$scratch"/"$sample_num"_intermediate1_no_contigs.fasta "$scratch"/"$sample_num"_intermediate1_myIndex 

echo "Mapping reads to reference genome"
bowtie2-align-s -I 0 -X 800 -p 32 --sensitive -q -x "$scratch"/"$sample_num"_intermediate1_myIndex -1 "$scratch"/"$sample_num"_trimmed_1.fastq -2 "$scratch"/"$sample_num"_trimmed_2.fastq -S "$scratch"/"$sample_num"_intermediate1.sam

# Previously generated samtools index of reference genome.  Generates *.fai file and only need to do 1X.
# This step may not be necessary with the samtool pipeline below:
echo "Indexing inetermediate1 genome with samtools"
samtools faidx "$scratch"/"$sample_num"_intermediate1_no_contigs.fasta

echo "Starting samtools - convert SAM to BAM"
samtools view -S -b -o "$scratch"/"$sample_num"_intermediate1.bam "$scratch"/"$sample_num"_intermediate1.sam

echo "Starting samtools sort BAM file"
samtools sort -@ $NSLOTS "$scratch"/"$sample_num"_intermediate1.bam -o "$scratch"/"$sample_num"_intermediate1.sorted.bam

echo "Starting samtools index BAM file"
samtools index "$scratch"/"$sample_num"_intermediate1.sorted.bam


JAVA_OPTS='-Xmx50g'
TMP_DIR=/tmp


#Copy only mapped bases, and save:
samtools view -b -F 4 "$scratch"/"$sample_num"_intermediate1.sorted.bam > "$scratch"/"$sample_num"_intermediate1-mapped.sorted.bam


#Make intermediate1 fasta:
echo "Making intermediate1 fasta!"
bcftools mpileup -A -d 6000000 -B -Q 0 -Ov -f "$scratch"/"$sample_num"_intermediate1_no_contigs.fasta "$scratch"/"$sample_num"_intermediate1-mapped.sorted.bam | bcftools call --ploidy 1 -mv -Ov -o "$scratch"/"$sample_num"_intermediate1-mapped.sorted.bam.vcf

gatk IndexFeatureFile --input "$scratch"/"$sample_num"_intermediate1-mapped.sorted.bam.vcf

gatk CreateSequenceDictionary -R "$scratch"/"$sample_num"_intermediate1_no_contigs.fasta

gatk FastaAlternateReferenceMaker -R "$scratch"/"$sample_num"_intermediate1_no_contigs.fasta -O "$scratch"/"$sample_num"_intermediate2_no_contigs.fasta -V "$scratch"/"$sample_num"_intermediate1-mapped.sorted.bam.vcf

#**********************************************************************************************************

echo "Performing mapping to Self #2"

echo "Indexing reference sequence using bowtie2"
bowtie2-build-s "$scratch"/"$sample_num"_intermediate2_no_contigs.fasta "$scratch"/"$sample_num"_intermediate2_myIndex 

echo "Mapping reads to reference genome"
bowtie2-align-s -I 0 -X 800 -p 32 --sensitive -q -x "$scratch"/"$sample_num"_intermediate2_myIndex -1 "$scratch"/"$sample_num"_trimmed_1.fastq -2 "$scratch"/"$sample_num"_trimmed_2.fastq -S "$scratch"/"$sample_num"_intermediate2.sam

# Previously generated samtools index of reference genome.  Generates *.fai file and only need to do 1X.
# This step may not be necessary with the samtool pipeline below:
echo "Indexing inetermediate1 genome with samtools"
samtools faidx "$scratch"/"$sample_num"_intermediate2_no_contigs.fasta

echo "Starting samtools - convert SAM to BAM"
samtools view -S -b -o "$scratch"/"$sample_num"_intermediate2.bam "$scratch"/"$sample_num"_intermediate2.sam

echo "Starting samtools sort BAM file"
samtools sort -@ $NSLOTS "$scratch"/"$sample_num"_intermediate2.bam -o "$scratch"/"$sample_num"_intermediate2.sorted.bam

echo "Starting samtools index BAM file"
samtools index "$scratch"/"$sample_num"_intermediate2.sorted.bam


JAVA_OPTS='-Xmx50g'
TMP_DIR=/tmp


#Copy only mapped bases, and save:
samtools view -b -F 4 "$scratch"/"$sample_num"_intermediate2.sorted.bam > "$scratch"/"$sample_num"_intermediate2-mapped.sorted.bam

samtools index "$scratch"/"$sample_num"_intermediate2-mapped.sorted.bam

#Make intermediate1 fasta:
echo "Making final fasta!"

echo "Making consensus fasta!"

samtools mpileup -r 1 -A -aa -d 6000000 -B -Q 0 -f "$scratch"/"$sample_num"_intermediate2_no_contigs.fasta "$scratch"/"$sample_num"_intermediate2-mapped.sorted.bam | /scicomp/home-pure/evk3/setup/ivar-master_1.3/src/ivar consensus -p "$scratch"/"$sample_num".consensus -m 2 -n N


#**********************************************************************************************************


#copy results from node /scratch/evk3/ebo back to home dir

cp "$scratch"/"$sample_num"_intermediate2_no_contigs.fasta "$OUTPUT_PATH"
cp "$scratch"/"$sample_num"_intermediate2-mapped.sorted.bam "$OUTPUT_PATH"
cp "$scratch"/"$sample_num".consensus.fa "$OUTPUT_PATH"
cp "$scratch"/"$sample_num".consensus.qual.txt "$OUTPUT_PATH"


module unload bwa/0.7.17
module unload samtools/1.9
module unload picard/2.21.1
module unload BEDTools/2.27.1
module unload cutadapt/2.3
module unload prinseq/0.20.3
module unload htslib/1.9
module unload bowtie2/2.3.5.1
module unload gatk/4.1.7.0
module unload htslib/1.10
module unload gcc/9.2.0

echo "Script finish"

