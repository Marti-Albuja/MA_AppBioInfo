# Set the trace
set -uex

# SRR number
SRR=SRR28572035

# Number of reads to sample
N=10000

# The output read names
R1=reads/${SRR}_1.fastq
R2=reads/${SRR}_2.fastq

# Trimmed read names
T1=reads/${SRR}_1.trimmed.fastq
T2=reads/${SRR}_2.trimmed.fastq

# The reads directory
RDIR=reads

# The reports directory
PDIR=reports

#
#-----NOTHING NEEDS TO BE CHANGED BELOW THIS LINE-----
#

# Make the necessary directories
mkdir -p ${RDIR} ${PDIR}

# Download the FASTQ file
fastq-dump -X ${N} --split-files -O ${RDIR} ${SRR} 

# Run fastqc
fastqc -q -o ${PDIR} ${R1} ${R2}

# Run fastp and trim for quality and read length
fastp --cut_right -f 30 -F 30 -T 80 -i ${R1} -I ${R2} -o ${T1} -O ${T2}

# Run fastqc
fastqc -q -o ${PDIR} ${T1} ${T2}

