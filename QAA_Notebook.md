# Quality Assessment Assignment:

The objectives of this assignment are:
- to use existing tools for quality assessment and adaptor trimming
- compare the quality assessments to those from your own software
- demonstrate your ability to summarize other important information about this RNA-Seq data set.

### Data:
My files:
1_2A_control_S1_L008    17_3E_fox_S13_L008

1_2A_control_S1_L008_R1_001.fastq.gz    
1_2A_control_S1_L008_R2_001.fastq.gz

17_3E_fox_S13_L008_R1_001.fastq.gz
17_3E_fox_S13_L008_R2_001.fastq.gz

## PART ONE: READ QUAITY SCORE DISTRIBUTIONS

1. Using FastQC via the command line on Talapas, produce plots of quality score distributions for R1 and R2 reads. 
Also, produce plots of the per-base N content, and comment on whether or not they are consistent with the quality score plots.

### FASTQC: QAA_P1.srun

#!/bin/bash
#SBATCH --account=bgmp 
#SBATCH --partition=bgmp 
#SBATCH --job-name=QAA_P1
#SBATCH --output=QAA_P1_%j.out
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=8 
#SBATCH --time=1:00:00 
#SBATCH --cpus-per-task=1

module load fastqc/0.11.5
/usr/bin/time -v zcat /projects/bgmp/shared/2017_sequencing/demultiplexed/1_2A_control_S1_L008_R1_001.fastq.gz | fastqc stdin -o S1_R1 
/usr/bin/time -v zcat /projects/bgmp/shared/2017_sequencing/demultiplexed/1_2A_control_S1_L008_R2_001.fastq.gz | fastqc stdin -o S1_R2

/usr/bin/time -v zcat /projects/bgmp/shared/2017_sequencing/demultiplexed/17_3E_fox_S13_L008_R1_001.fastq.gz | fastqc stdin -o S13_R1
/usr/bin/time -v zcat /projects/bgmp/shared/2017_sequencing/demultiplexed/17_3E_fox_S13_L008_R2_001.fastq.gz | fastqc stdin -o S13_R2

### JOB OUTPUT FILE: QAA_P1_16172097.out
* These were super fast! Each one took 10-16 sec

2. Run your quality score plotting script from your Demultiplexing assignment from Bi622. 
Describe how the FastQC quality score distribution plots compare to your own. 
If different, propose an explanation. Also, does the runtime differ? If so, why?

### P1_Dist.srun

#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --job-name=P1_Dist     
#SBATCH --output=P1_Dist_%j.out
#SBATCH --time=8:00:00
#SBATCH --nodes=1               
#SBATCH --ntasks-per-node=8   
#SBATCH --cpus-per-task=1

/usr/bin/time -v ./P1_Dist.py -f /projects/bgmp/shared/2017_sequencing/demultiplexed/1_2A_control_S1_L008_R1_001.fastq.gz -p 101 -o S1_R1_hist.png -r Read1

/usr/bin/time -v ./P1_Dist.py -f /projects/bgmp/shared/2017_sequencing/demultiplexed/1_2A_control_S1_L008_R2_001.fastq.gz -p 101 -o S1_R2_hist.png -r Read1

/usr/bin/time -v ./P1_Dist.py -f /projects/bgmp/shared/2017_sequencing/demultiplexed/17_3E_fox_S13_L008_R1_001.fastq.gz -p 101 -o S13_R1_hist.png -r Read1

/usr/bin/time -v ./P1_Dist.py -f /projects/bgmp/shared/2017_sequencing/demultiplexed/17_3E_fox_S13_L008_R2_001.fastq.gz -p 101 -o S13_R2_hist.png -r Read2

### JOB OUTPUT FILE: P1_Dist_16172572.out
*These runs were much slower, each one took ~6 minutes.

## PART TWO: ADAPTER TRIMMING COMPARISON

5. Using cutadapt, properly trim adapter sequences from your assigned files. Be sure to read how to use cutadapt. Use default settings. What proportion of reads (both forward and reverse) were trimmed?

Create environment: 
conda create --name QAA_env
conda install python=3.9.5
conda install -c bioconda cutadapt
conda install -c bioconda trimmomatic

cutadapt --version (should be 3.4)
trimmomatic -version (should be 0.39)

Read 1 Adapter: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
Read 2 Adapter: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

### RUN CUTADAPT:

#!/bin/bash
#SBATCH --account=bgmp 
#SBATCH --partition=bgmp 
#SBATCH --job-name=cutadapt
#SBATCH --output=cutadapt_%j.out
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=8 
#SBATCH --time=8:00:00 
#SBATCH --cpus-per-task=1

/usr/bin/time -v cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -j 0 -o S1_R1.fastq -p S1_R2.fastq /projects/bgmp/shared/2017_sequencing/demultiplexed/1_2A_control_S1_L008_R1_001.fastq.gz /projects/bgmp/shared/2017_sequencing/demultiplexed/1_2A_control_S1_L008_R2_001.fastq.gz
/usr/bin/time -v cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -j 0 -o S13_R1.fastq -p S13_R2.fastq /projects/bgmp/shared/2017_sequencing/demultiplexed/17_3E_fox_S13_L008_R1_001.fastq.gz /projects/bgmp/shared/2017_sequencing/demultiplexed/17_3E_fox_S13_L008_R2_001.fastq.gz

### JOB OUTPUT FILE: cutadapt_16175059.out
System time (seconds): 10.07

6. Use Trimmomatic to quality trim your reads. Specify the following, in this order:

LEADING: quality of 3
TRAILING: quality of 3
SLIDING WINDOW: window size of 5 and required quality of 15
MINLENGTH: 35 bases

### RUN TRIMMOMATIC: Trimmomatic.srun

#!/bin/bash
#SBATCH --account=bgmp 
#SBATCH --partition=bgmp 
#SBATCH --job-name=trimmomatic
#SBATCH --output=trimmomatic_%j.out
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=8 
#SBATCH --time=1:00:00 
#SBATCH --cpus-per-task=1

/usr/bin/time -v trimmomatic PE -threads 4 S1_R1.fastq S1_R2.fastq S1_R1.trimmed.fastq S1_R1un.trimmed.fastq S1_R2.trimmed.fastq S1_R2un.trimmed.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:35
/usr/bin/time -v trimmomatic PE -threads 4 S13_R1.fastq S13_R2.fastq S13_R1.trimmed.fastq S13_R1un.trimmed.fastq S13_R2.trimmed.fastq S13_R2un.trimmed.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:35

### JOB OUTPUT FILE: trimmomatic_16181904.out

Took trimmed files and turned them into TSVs. 
cat S1_R1.trimmed.fastq | sed -n "2~4p" | awk '{print length($0)}' > S1_R1.trimmed.nsv

## PART 3: ALIGNMENT AND STRAND SPECIFICITY

9. Find publicly available mouse genome fasta files (Ensemble release 104) and generate an alignment database from them. Align the reads to your mouse genomic database using a splice-aware aligner. Use the settings specified in PS8 from Bi621.

went to Ensembl > Downloads > FTP Downloads > Mouse > DNA FASTA & GTF FASTA

Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
Mus_musculus.GRCm39.104.gtf.gz

Install STAR: conda install star -c bioconda
INSTALL NUMPY: conda install numpy -c bioconda
INSTALL PYSAM: conda install pysam -c bioconda
pip install HTSeq

### GENERATE DB: STAR_generate_db.srun

#!/bin/bash
#SBATCH --account=bgmp 
#SBATCH --partition=bgmp
#SBATCH --job-name=STAR_generate_db
#SBATCH --output=STAR_%j.out
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --time=2:00:00 
#SBATCH --cpus-per-task=8 

/usr/bin/time -v STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir /projects/bgmp/nkearns/bioinformatics/Bi623/QAA/Mus_musculus.GRCm39.dna.ens104.STAR_2.7.9a \
--genomeFastaFiles /projects/bgmp/nkearns/bioinformatics/Bi623/QAA/Mus_musculus.GRCm39.dna.primary_assembly.fa \
--sjdbGTFfile /projects/bgmp/nkearns/bioinformatics/Bi623/QAA/Mus_musculus.GRCm39.104.gtf
--outFileNamePrefix Align_OUTPUT

### JOB OUTPUT FILE: 
STAR_16199053.out
System time (seconds): 54.65

### ALIGN TO DATABASE: STAR_Align.srun

#!/bin/bash
#SBATCH --account=bgmp 
#SBATCH --partition=bgmp
#SBATCH --job-name=STAR_align
#SBATCH --output=STAR_align_%j.out
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --time=2:00:00 
#SBATCH --cpus-per-task=8 

/usr/bin/time -v STAR --runThreadN 8 \
--runMode alignReads \
--outFilterMultimapNmax 3 \
--outSAMunmapped Within KeepPairs \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--readFilesCommand zcat \
--readFilesIn /projects/bgmp/nkearns/bioinformatics/Bi623/QAA/S1_R1.trimmed.fastq.gz /projects/bgmp/nkearns/bioinformatics/Bi623/QAA/S1_R2.trimmed.fastq.gz \
--genomeDir /projects/bgmp/nkearns/bioinformatics/Bi623/QAA/Mus_musculus.GRCm39.dna.ens104.STAR_2.7.9a \
--outFileNamePrefix S1_Align_to_ref

### JOB OUTPUT FILE: STAR_align_16206461.out

10. USE SCRIPT FROM PS8 TO COUNT MAPPED AND UNMAPPED READS

./ParseSam.py -f S1_Align_to_refAligned.out.sam

1_2A_control_S1:
Mapped: 15627450
Unmapped: 305236

./ParseSam.py -f S13_Align_to_refAligned.out.sam

17_3E_fox_S13:
Mapped: 21532850
Unmapped: 948692

11. Count reads that map to features using htseq-count. You should run htseq-count twice: once with --stranded=yes and again with --stranded=no. Use default parameters otherwise.

### HTSEQ.srun

#!/bin/bash
#SBATCH --account=bgmp 
#SBATCH --partition=bgmp
#SBATCH --job-name=htseq
#SBATCH --output=HTSEQ_%j.out
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --time=8:00:00 
#SBATCH --cpus-per-task=8 

/usr/bin/time -v htseq-count -f sam --stranded=yes --samout=S1_htseq_stranded_out  /projects/bgmp/nkearns/bioinformatics/Bi623/QAA/S1_Alignment/S1_Align_sorted /projects/bgmp/nkearns/bioinformatics/Bi623/QAA/Mus_musculus.GRCm39.104.gtf > S1.Stranded_htseq.out
/usr/bin/time -v htseq-count -f sam --stranded=no --samout=S1_htseq_out  /projects/bgmp/nkearns/bioinformatics/Bi623/QAA/S1_Alignment/S1_Align_sorted /projects/bgmp/nkearns/bioinformatics/Bi623/QAA/Mus_musculus.GRCm39.104.gtf > S1.NotSranded_htseq.out
/usr/bin/time -v htseq-count -f sam --stranded=yes --samout=S13_htseq_stranded_out  /projects/bgmp/nkearns/bioinformatics/Bi623/QAA/S13_Alignment/S13_Align_sorted /projects/bgmp/nkearns/bioinformatics/Bi623/QAA/Mus_musculus.GRCm39.104.gtf > S13.Stranded_htseq.out
/usr/bin/time -v htseq-count -f sam --stranded=no --samout=S13_htseq_out  /projects/bgmp/nkearns/bioinformatics/Bi623/QAA/S13_Alignment/S13_Align_sorted /projects/bgmp/nkearns/bioinformatics/Bi623/QAA/Mus_musculus.GRCm39.104.gtf > S13.NotStranded_htseq.out

### JOB OUTPUT FILE: HTSEQ_16215111.out

S1_stranded
no_feature    7163565 = 94.5 %

S1_NOT_stranded
no_feature    408922 = 5.5% %

S13_stranded
no_feature    9781740 = 90.8 %

S13_NOT_stranded
no_feature    994711 = 1.8 %








