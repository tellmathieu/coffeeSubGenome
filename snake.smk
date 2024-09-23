import os, sys, glob

#######################
# Snakemake adaptation of Thales Henrique Cherubino Ribeiro's
# Thesis Chapter 1 - Starting at the augustus steps
# By Tell Mathieu
# Starts with Thales bonafide.gb file and other prep files he made. 
# He sent these to me and they are in the provided folder
# Just testing it on the updated genome from this paper: 
# 	https://www.nature.com/articles/s41588-024-01695-w
#######################

#######################
#     Variables       #
#######################

# Fill in this info - with the provided folder from Thales 
# - I renamed the subgenome folders to match my naming convention (subspecies variable)
subspecies=['subEugenioides', 'subCanephora']
subspeciesFastaChr=['EE', 'CC']
providedFolder='/home/tmathieu/coffeaArabica/RNAseqBasedPrediction'
cdsTranscriptome='0.80.colapsed.all.transcriptomes.fasta.transdecoder.cds'
prevGenomeLink='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/713/225/GCF_003713225.1_Cara_1.0/GCF_003713225.1_Cara_1.0_genomic.fna.gz'
prevGenomeShorthand = 'Cara_1.0'
hifiGenomeLink='https://bioinformatics.psb.ugent.be/gdb/coffea_arabica/Cara_r6_tgs15k_ChrNames_racon_w3kb_chrs_NoOrganelles.flip.fix.fa.gz'
hifiGenomeShorthand = 'Cara_hifi'

###############################
# Don't Edit after this point #
###############################

projectFolder=os.getcwd()

prevGenomeName = prevGenomeLink.split('/')[-1]
hifiGenomeName = hifiGenomeLink.split('/')[-1]

prevGenomeUnzipped = prevGenomeName.split('.')[0:-1]
hifiGenomeUnzipped = hifiGenomeName.split('.')[0:-1]

prevGenomeCleaned = prevGenomeShorthand + '_cleaned.fa'
hifiGenomeCleaned = hifiGenomeShorthand + '_cleaned.fa'



#probably unnecessary - check later
binPath='/data/data2/homes/tmathieu/miniforge3/envs/prediction/bin'
cfgPath='/data/data2/homes/tmathieu/miniforge3/envs/prediction/config'
fasta='/home/tmathieu/coffeaArabica/genome/EEsubgenome.chr.fa'
subGenomeFolder='/home/tmathieu/coffeaArabica/RNAseqBasedPrediction/ceu.subgenome'
cfgPath='/data/data2/homes/tmathieu/miniforge3/envs/prediction'

#cfgPath=$(echo $binPath | rev | cut -d"/" -f2- | rev)

rule all:
	input:
		prevGenomeCleaned = os.path.join(projectFolder, 'genome', prevGenomeCleaned),
		hifiGenomeCleaned = os.path.join(projectFolder, 'genome', hifiGenomeCleaned)

rule downloadGenomes:
	params:
		prevGenomeLink = prevGenomeLink,
		hifiGenomeLink = hifiGenomeLink,
		genomeDir = os.path.join(projectFolder, 'genome'),
		prevGenomeName = prevGenomeName,
		hifiGenomeName = hifiGenomeName
	output:
		prevGenomeFile = os.path.join(projectFolder, 'genome', prevGenomeName),
		hifiGenomeFile = os.path.join(projectFolder, 'genome', hifiGenomeName)
	shell: '''
		# download genomes used - original (used by Thales) and new (from the nature article)
		wget {params.prevGenomeLink}
		wget {params.hifiGenomeLink}

		# make directory and move genomes into it
		mkdir -p {params.genomeDir}
		mv {params.prevGenomeName} {output.prevGenomeFile}
		mv {params.hifiGenomeName} {output.hifiGenomeFile}
	'''

rule unzipGenomes:
	input:
		prevGenomeFile = os.path.join(projectFolder, 'genome', prevGenomeName),
		hifiGenomeFile = os.path.join(projectFolder, 'genome', hifiGenomeName)
	output:
		prevGenomeFileUnzipped = os.path.join(projectFolder, 'genome', prevGenomeUnzipped),
		hifiGenomeFileUnzipped = os.path.join(projectFolder, 'genome', hifiGenomeUnzipped)
	shell: '''
		# unzip prevGenome
		gunzip -k {input.prevGenomeFile}

		# unzip hifiGenome
		gunzip -k {input.hifiGenomeFile}
	'''

rule cleanGenomes:
	input:
		prevGenomeFileUnzipped = os.path.join(projectFolder, 'genome', prevGenomeUnzipped),
		hifiGenomeFileUnzipped = os.path.join(projectFolder, 'genome', hifiGenomeUnzipped)
	output:
		prevGenomeCleaned = os.path.join(projectFolder, 'genome', prevGenomeCleaned),
		hifiGenomeCleaned = os.path.join(projectFolder, 'genome', hifiGenomeCleaned)
	shell: '''
		# clean prevGenome
		awk '{{print $1}}' {input.prevGenomeFileUnzipped} > temp_prev_genome.fa
		sed -e 's/\(^>.*$\)/#\1#/' temp_prev_genome.fa | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' > {output.prevGenomeCleaned}

		# clean hifiGenome
		awk '{{print $1}}' {input.hifiGenomeFileUnzipped} > temp_hifi_genome.fa
		sed -e 's/\(^>.*$\)/#\1#/' temp_hifi_genome.fa | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' > {output.hifiGenomeCleaned}

		# remove temporary files
		rm temp*
	'''

rule prep:
	input:
		subspecies = subspecies
	output:


	shell:'''
	# Making directory for later steps
	mkdir -p $subspecies

	# Augustus script to set up the software for using an organism that isn't prebuilt
	new_species.pl --species=$subspecies
	'''
