# coffeeSubGenome

# Running this pipeline

## 1. Clone this repository to your working space

```
git clone https://github.com/tellmathieu/coffeeSubGenome.git
```

## 2. Set up an environment to run this pipeline
- a) install mamba if you don't already have it - you'll have to go through the installation process or do a command like `module load mamba`
```
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```
- b) Initialize mamba and set up an environment - this can take 5-10 minute
```
source ~/.bashrc #or you can log out and then log back in
mamba create -c conda-forge -c bioconda -y -n subgenome augustus braker3 snakemake ucsc-liftover conda-forge::ruby ucsc-ldhggene ucsc-pslcdnafilter ucsc-fasplit parallel ucsc-blat busco
mamba activate subgenome 
```

## 3. Edit 1snake.smk to your variables and filepaths

```
#######################
#     Variables       #
#######################

# Fill in this info - with the provided folder from Thales 
# - I renamed the subgenome folders to match my naming convention (subspecies variable)
subspecies=['subEugenioides', 'subCanephora']
subspeciesConversion='/home/tmathieu/testingCofSM/subGenomeConversion.txt'
cdsTranscriptome='/home/tmathieu/coffeaArabica/0.80.colapsed.all.transcriptomes.fasta.transdecoder.cds'
prevGenomeLink='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/713/225/GCF_003713225.1_Cara_1.0/GCF_003713225.1_Cara_1.0_genomic.fna.gz'
prevGenomeShorthand='Cara_1.0'
hifiGenomeLink='https://bioinformatics.psb.ugent.be/gdb/coffea_arabica/Cara_r6_tgs15k_ChrNames_racon_w3kb_chrs_NoOrganelles.flip.fix.fa.gz'
hifiGenomeShorthand='Cara_hifi'
providedFolder='/home/tmathieu/coffeaArabica/RNAseqBasedPrediction' #I got this file from Thales
chrConversion='/home/tmathieu/testingCofSM/chrConversion.txt'
environmentFolder='/home/tmathieu/miniforge3/envs/liftnew'
numGenesTest=3000 #chosen by Thales - see his script
```

## 4. Run pipeline - this will likely take more than a week unless you allocate more resources to it

```
snakemake --snakefile  1snake.smk -j 20  -p 
```


