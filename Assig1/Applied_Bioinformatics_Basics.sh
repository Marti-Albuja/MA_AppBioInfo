#Applied Bioinformatics

#To activate through conda or micromamba
conda activate bioinfo
micromamba activate bioinfo

#Create a folder
mkdir -p work # -p means that you can create intermediate folders inside an original folder

# Listing file name paths
find .

# Create an environment
micromamba create -y -n heatshock
micromamba create --name bedops
micromamba create -y --name foo python=3.8
micromamba create -n bioinfo -y python=3.10

# Install the software
micromamba install fastqc bwa bcftools trimmomatic

#List the current environments
micromamba env list

#To list the channel priority
micromamba config list

#To change the order of the channels
micromamba config prepend channels conda-forge

#To remove and install things inside an environment
micromamba remove bcftools
micromamba install bcftools

#Update packages inside an environment
micromamba update --all
micromamba update toolname

#Procedure to work with micromamba
micromamba create --name mynewthing -y
micromamba activate mynewthing
micromamba install toolname

# Enable the bioconda channels.
micromamba config prepend channels conda-forge
micromamba config append channels bioconda

# Set the channel priority.
micromamba config set channel_priority strict

#To know where you are
pwd

#List files inside a directory
ls

#Change directory
cd name_directory
cd #takes you to home directory

#To go back
cd ..

#Manual of a tool
man nametool

#Create files
touch file_name

#Remove directories
rmdir directory_name

#Viewing files
more
less
cat 

#Editing with nano

nano text_name

#To append text to a file
echo "text" >> file_name

#To count lines in a file
wc -l file_name

#To count words in a file
wc -w file_name

#To count characters in a file
wc -c file_name

#To eliminate or don't include something when using grep
grep -v "text" file_name

#To sort things in reverse order
sort -r file_name

#To sort things numerically
sort -rn file_name
sort-uniq-count-rank

#Stats of a file
seqkit stats file_name

#Index a fasta file
samtools faidx file_name

#CDS are the coding sequences where you can visualize the start and end codons. Genes also have untranslated sequence.

