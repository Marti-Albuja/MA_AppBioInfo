#Applied Bioinformatics

#To activate through conda or micromamba
conda activate bioinfo
micromamba activate bioinfo

#Create a folder
mkdir -p work

# Listing file name paths
find .

# Create an environment
micromamba create -y -n heatshock
micromamba create --name bedops
micromamba create -y --name foo python=3.8

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
