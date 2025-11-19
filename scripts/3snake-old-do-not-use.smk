import os, sys, glob

#######################
# Snakemake adaptation of Thales Henrique Cherubino Ribeiro's
# Thesis Chapter 1 - Starting at the augustus steps
# By Tell Mathieu
# Starts with Thales bonafide.gb file and other prep files he made. 
# He sent these to me and they are in the provided folder
# Just testing it on the updated genome from this paper: 
# 	https://www.nature.com/articles/s41588-024-01695-w
# Biological question: Is there subgenome dominance in Coffea arabica? Is one of the subgenomes losing universal orthologs faster than others?
# Collaborators: Thales Henrique Cherubino Ribeiro, Blake Meyers
#######################

#######################
#     Variables       #
#######################

# Fill in this info - with the provided folder from Thales 
# - I renamed the subgenome folders to match my naming convention (subspecies variable)
subspecies=['subEugenioides', 'subCanephora']
subspeciesConversion='/quobyte/bcmeyersgrp/tell/testingCofSM/subGenomeConversion.txt'
cdsTranscriptome='/quobyte/bcmeyersgrp/tell/coffee/coffeaArabica/0.80.colapsed.all.transcriptomes.fasta.transdecoder.cds'
prevGenomeLink='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/713/225/GCF_003713225.1_Cara_1.0/GCF_003713225.1_Cara_1.0_genomic.fna.gz'
prevGenomeShorthand='Cara_1' # technically 1.0, but removing the dot made this simpler to run
hifiGenomeLink='https://bioinformatics.psb.ugent.be/gdb/coffea_arabica/Cara_r6_tgs15k_ChrNames_racon_w3kb_chrs_NoOrganelles.flip.fix.fasta.gz'
hifiGenomeShorthand='Cara_hifi'
providedFolder='/quobyte/bcmeyersgrp/tell/coffee/coffeaArabica/RNAseqBasedPrediction' #I got this file from Thales
chrConversion='/quobyte/bcmeyersgrp/tell/testingCofSM/chrConversion.txt'
environmentFolder='/home/tmathieu/.conda/envs/liftnew'
numGenesTest=3000 #chosen by Thales - see his script

###############################
# Don't Edit after this point #
###############################

projectFolder=os.getcwd()

genomeList = [prevGenomeShorthand, hifiGenomeShorthand]
prevGenomeGZ = prevGenomeShorthand + '.fa.gz'
prevGenomeCleaned = prevGenomeShorthand + '_cleaned.fa'
hifiGenomeGZ = hifiGenomeShorthand + '.fa.gz'
hifiGenomeCleaned = hifiGenomeShorthand + '_cleaned.fa'

FLO = ['CITATION.cff','README.md','Rakefile','gff_compare.rb','gff_helpers.rb','gff_longest_transcripts.rb','gff_recover.rb','gff_remove_feats.rb','opts_example.yaml','scripts/all_mrna.rb','scripts/all_mrna_cds.rb','scripts/clean_NCBI_gff.sh','scripts/how_many_split_across_scaffolds.rb','scripts/install.sh','scripts/longest_mrna.rb','scripts/longest_mrna_cds.rb','scripts/qa.rb','scripts/remove_mitochondrial_annotations.sh','scripts/trim_NCBI_seqids_to_accession_numbers.sh']

rule all:
	input: 
		#os.path.join(projectFolder, 'busco.pdf'),
		os.path.join(projectFolder, 'buscoUTR.pdf')

rule downloadGenomes:
	params:
		prevGenomeLink = prevGenomeLink,
		hifiGenomeLink = hifiGenomeLink,
		genomeDir = os.path.join(projectFolder, 'genome')
	output:
		prevGenomeFile = os.path.join(projectFolder, 'genome', prevGenomeGZ),
		hifiGenomeFile = os.path.join(projectFolder, 'genome', hifiGenomeGZ)
	shell: '''
		# make directory if not already present
		mkdir -p {params.genomeDir}

		# download genomes used - original (used by Thales) and new (from the nature article)
		wget -O - {params.prevGenomeLink} > {output.prevGenomeFile}
		wget -O - {params.hifiGenomeLink} > {output.hifiGenomeFile}
	'''

rule unzipGenomes:
	input:
		genomeFile = os.path.join(projectFolder, 'genome', '{genome}.fa.gz')
	output:
		genomeFileUnzipped = os.path.join(projectFolder, 'genome', '{genome}.fa')
	shell: '''
		# unzip prevGenome
		gunzip -k {input.genomeFile}
	'''

rule cleanGenomes:
	input:
		genomeFileUnzipped = os.path.join(projectFolder, 'genome', '{genome}.fa')
	output:
		genomeCleaned = os.path.join(projectFolder, 'genome', '{genome}_cleaned.fa')
	shell: '''
		# clean genome
		sed -e 's/\\(^>.*$\\)/#\\1#/' {input.genomeFileUnzipped} | tr -d "\\r" | tr -d "\\n" | sed -e 's/$/#/' | tr "#" "\\n" | sed -e '/^$/d' > temp_genome_{wildcards.genome}.fa
		awk '{{print $1}}' temp_genome_{wildcards.genome}.fa > {output.genomeCleaned}

		# remove temporary files
		rm temp_genome_{wildcards.genome}.fa
	'''

rule convertChrGenome:
	input:
		chrConversion = chrConversion,
		genomeCleaned = os.path.join(projectFolder, 'genome', '{genome}_cleaned.fa')
	output:
		genomeConverted = os.path.join(projectFolder, 'genome', '{genome}_converted.fa')
	shell: '''
		# get # of lines in the conversion document
		lines=$(cat {input.chrConversion} | wc -l)
		
		# set a counting var
		line=1

		# make copy of cleaned genome - so nothing happens to original
		cp {input.genomeCleaned} temp_{wildcards.genome}_chr.fa

		# iterate through the lines of the conversion doc
		while [ $line -le $lines ]
			do
				# get old chr name
				{wildcards.genome}oldChrName=$(awk '{{ print $1 }}' {input.chrConversion} | awk -v var1="$line" 'NR==var1')
				
				# get new chr name
				{wildcards.genome}newChrName=$(awk '{{ print $2 }}' {input.chrConversion} | awk -v var1="$line" 'NR==var1')
				
				# replace old with new
				sed -i -e "s/${wildcards.genome}oldChrName/${wildcards.genome}newChrName/g" temp_{wildcards.genome}_chr.fa
				
				# add one to line var so we go to next line
				let line++
			done

		# rename temp file to final file
		mv temp_{wildcards.genome}_chr.fa {output.genomeConverted}
	'''

rule rmUnknownsGenome:
	input:
		chrConversion = chrConversion,
		genomeConverted = os.path.join(projectFolder, 'genome', '{genome}_converted.fa')
	output:
		genomeFinal = os.path.join(projectFolder, 'genome', '{genome}_final.fa')
	shell: '''
		# get # of lines in the conversion document
		lines=$(cat {input.chrConversion} | wc -l)

		# set a counting var
		line=1

		# iterate through the lines of the conversion doc
		while [ $line -le $lines ]
			do
				# get new chr name
				{wildcards.genome}newChrName=$(awk '{{ print $2 }}' {input.chrConversion} | awk -v var1="$line" 'NR==var1')
				
				# only add fasta record to temp doc if its on our list
				sed -n "/>${wildcards.genome}newChrName/{{p;n;p}}" {input.genomeConverted} >> temp_{wildcards.genome}_rep.fa
				
				# add one to line var so we go to next line
				let line++
			done

		# rename temp file to final file
		mv temp_{wildcards.genome}_rep.fa {output.genomeFinal}
	'''

rule separateSubGenomes:
	input:
		genomeFinal = expand(os.path.join(projectFolder, 'genome', '{genome}_final.fa'), zip, genome=genomeList),
		subspeciesConversion = subspeciesConversion
	params:
		subgenomeDir = os.path.join(projectFolder, '{subgenome}'),
		genomeDir = os.path.join(projectFolder, 'genome'),
	output:
		subGenomeFinal = os.path.join(projectFolder, '{subgenome}', '{subgenome}_{genome}_final.fa')
	shell: '''
		# make directory if not already present
		mkdir -p {params.subgenomeDir}

		# grab chr name - I add the subgenome name to the variables, because I problems in the past when running in parallel
		{wildcards.subgenome}chrName=$(grep {wildcards.subgenome} {input.subspeciesConversion} | awk '{{ print $1 }}')

		# put only the grabbed chr names in a new fasta file - separating the genomes
		grep -A 1 "sg${wildcards.subgenome}chrName" {params.genomeDir}/{wildcards.genome}_final.fa > {output.subGenomeFinal}

	'''

# Next we have to combine the hints files that Thales generated 
# and then lift (get new coodinates on the new genome) them 
# from the old genome to the new genome
# The files he generated: prot.hints, rnaseq.gff, introns.f.gff
# There wasn't a provided CDNA hints file, so I recreated it using his cds transcriptome file
# the next few steps are for processing the CDNA hints that weren't included in the provided folder from Thales
#align transcriptome with blat
rule blatSubGenomeCDNAGFF:
	input:
		cdsTranscriptome = cdsTranscriptome,
		subGenomeSource = os.path.join(projectFolder, '{subgenome}', '{subgenome}_' + prevGenomeShorthand + '_final.fa')
	params:
		subgenomeProcessingDir = os.path.join(projectFolder, '{subgenome}', 'processing')
	output:
		os.path.join(projectFolder, '{subgenome}', 'processing', 'cdna_{subgenome}.psl')
	shell: '''
		blat -noHead -minIdentity=92 \
			{input.subGenomeSource} \
			{input.cdsTranscriptome} \
			{output}
	'''

#filter
rule filterSubGenomeCDNAGFF:
	input:
		os.path.join(projectFolder, '{subgenome}', 'processing', 'cdna_{subgenome}.psl')
	output:
		os.path.join(projectFolder, '{subgenome}', 'processing', 'cdna_{subgenome}.f.psl')
	shell: '''
		pslCDnaFilter -minId=0.9 -localNearBest=0.005 \
			-ignoreNs -bestOverlap \
			{input} \
			{output}
	'''

#sort
rule sortSubGenomeCDNAGFF:
	input:
		os.path.join(projectFolder, '{subgenome}', 'processing', 'cdna_{subgenome}.f.psl')
	output:
		os.path.join(projectFolder, '{subgenome}', 'processing', 'cdna_{subgenome}.fs.psl')
	shell: '''
		cat {input} \
			| sort -n -k 16,16 \
			| sort -s -k 14,14 \
			> {output}
	'''

#convert to hints file
rule convertSubGenomeCDNAGFF:
	input:
		os.path.join(projectFolder, '{subgenome}', 'processing', 'cdna_{subgenome}.fs.psl')
	output:
		os.path.join(projectFolder, '{subgenome}', 'processing', 'cdna_{subgenome}.hints')
	shell: '''
		blat2hints.pl \
			--in={input} \
			--out={output} \
			--minintronlen=35 \
			--trunkSS
	'''

rule getSubGenomeHintsGFF:
	input:
		protHints = os.path.join(providedFolder, '{subgenome}', 'prot.hints'),
		cdnaHints = os.path.join(projectFolder, '{subgenome}', 'processing', 'cdna_{subgenome}.hints'),
		rnaseqHints = os.path.join(providedFolder, '{subgenome}', 'rnaseq.gff'),
		intronsHints = os.path.join(providedFolder, '{subgenome}', 'introns.f.gff')
	output:
		subGenomeHints = os.path.join(projectFolder, '{subgenome}', '{subgenome}_hints.gff')
	shell: '''
		# concatenating all the subgenome hints
		cat {input.protHints} \
		{input.cdnaHints} \
		{input.rnaseqHints} \
		{input.intronsHints} \
		> {output.subGenomeHints}
	'''

rule convertChrHints:
	input:
		chrConversion = chrConversion,
		subGenomeHints = os.path.join(projectFolder, '{subgenome}', '{subgenome}_hints.gff')
	output:
		convertedHints = os.path.join(projectFolder, '{subgenome}', '{subgenome}_hints_converted.gff')
	shell: '''
		# this process is similar to rules: convertChrGenome and rmUnknownsGenome
		lines=$(cat {input.chrConversion} | wc -l)
		line=1
		cp {input.subGenomeHints} temp_{wildcards.subgenome}_chr_hints.gff

		while [ $line -le $lines ]
			do
				oldChrName=$(awk '{{ print $1 }}' {input.chrConversion} | awk -v var1="$line" 'NR==var1')
				newChrName=$(awk '{{ print $2 }}' {input.chrConversion} | awk -v var1="$line" 'NR==var1')
				sed -i -e "s/$oldChrName/$newChrName/g" temp_{wildcards.subgenome}_chr_hints.gff
				let line++
			done

		mv temp_{wildcards.subgenome}_chr_hints.gff {output.convertedHints}
	'''

#this is to download a pipeline for lifting the hints to the new genome for when we run the file augustus
rule downloadFlo:
	output:
		os.path.join(projectFolder, 'flo.tar.gz')
	shell: '''
		# download flo
		wget -c https://github.com/yeban/flo/archive/master.tar.gz -O flo.tar.gz
	'''

rule gunzipFlo:
	input:
		os.path.join(projectFolder, 'flo.tar.gz')
	output:
		expand(os.path.join(projectFolder, 'flo-master', '{flo}'),flo=FLO)
	shell: '''
		# uncompress flo
		tar xvf flo.tar.gz
		wait
	'''

rule moveFlo:
	input:
		expand(os.path.join(projectFolder, 'flo-master', '{flo}'),flo=FLO)
	output:
		expand(os.path.join(projectFolder, 'flo', '{flo}'),flo=FLO)
	shell: '''
		# change name of uncompressed flo file - flo requires this
		if [ -d flo/ ]; then
			echo "flo exists - removing"
			rm -rdf flo
			# make dir
			mkdir -p flo
			mv flo-master/* flo/
		else
			# make dir
			mkdir -p flo
			mv flo-master/* flo/
		fi

		# remove old folder
		rm -rdf flo-master
	'''

rule installFlo:
	input:
		expand(os.path.join(projectFolder, 'flo', '{flo}'),flo=FLO),
		installFloScript = os.path.join(projectFolder, 'flo', 'scripts', 'install.sh')
	output:
		floInstallDoneTxt = os.path.join(projectFolder, 'floInstallDone.txt')
	shell: '''
		# run install script
		bash {input.installFloScript}

		# create file that says this process is done - it was easier this way than finding the new files
		echo Flo Install Done! > {output.floInstallDoneTxt}
	'''

rule chooseOptionsForFlo:
	input:
		expand(os.path.join(projectFolder, 'flo', '{flo}'),flo=FLO),
		floInstallDoneTxt = os.path.join(projectFolder, 'floInstallDone.txt'),
		subGenomeSource = os.path.join(projectFolder, '{subgenome}', '{subgenome}_' + prevGenomeShorthand +'_final.fa'),
		subGenomeTarget = os.path.join(projectFolder, '{subgenome}', '{subgenome}_' + hifiGenomeShorthand + '_final.fa'),
		convertedHints = os.path.join(projectFolder, '{subgenome}', '{subgenome}_hints_converted.gff')
	params:
		projectFolder = projectFolder
	output:
		optYaml = os.path.join(projectFolder, '{subgenome}', 'flo_opts.yaml')
	shell: '''
		# creating yaml file using echo
		echo ":add_to_path:" > {output.optYaml}
		echo "  - '{params.projectFolder}/ext/kent/bin'" >> {output.optYaml}
		echo "  - '{params.projectFolder}/ext/parallel-20150722/src'" >> {output.optYaml}
		echo "  - '{params.projectFolder}/ext/genometools-1.5.6/bin'" >> {output.optYaml}
		echo ":source_fa: '{input.subGenomeSource}'" >> {output.optYaml}
		echo ":target_fa: '{input.subGenomeTarget}'" >> {output.optYaml}
		echo ":processes: '20'" >> {output.optYaml}
		echo ":blat_opts: '-fastMap -tileSize=12 -minIdentity=98'"  >> {output.optYaml}
		echo ":lift:" >> {output.optYaml}
		echo "  - '{input.convertedHints}'" >> {output.optYaml}
	'''

# running the flo pipeline on both subgenomes
# had to use a script file for this step, since you need to change directories
rule runLiftOverSubGenomes:
	threads: 20
	input:
		optYaml = expand(os.path.join(projectFolder, '{subgenome}', 'flo_opts.yaml'), subgenome = subspecies),
		floInstallDoneTxt = os.path.join(projectFolder, 'floInstallDone.txt')
	params:
		runFloScript = os.path.join(projectFolder, 'scripts', 'runFlo.sh'),
		projectFolder = projectFolder,
		subgenomes = subspecies
	output:
		liftOverCompleteTxt = os.path.join(projectFolder, 'LiftOverComplete.txt'),
	shell: '''
		# run flo
		bash {params.runFloScript} {params.subgenomes} {params.projectFolder} ../flo/Rakefile

		# show finished
		echo Flo Liftovers complete! > {output.liftOverCompleteTxt}
	'''	

rule copyLiftedHints:
	input:
		liftOverCompleteTxt = os.path.join(projectFolder, 'LiftOverComplete.txt')
	params:
		liftedHints = os.path.join(projectFolder, '{subgenome}', 'run', '{subgenome}_nogenes_hints', 'lifted.gff3')
	output:
		renamedLiftedHints = os.path.join(projectFolder, '{subgenome}', '{subgenome}_lifted_hints.gff3')
	shell: '''
		# renaming lifted hints to include subgenome name
		cp {params.liftedHints} {output.renamedLiftedHints}
	'''

rule newSpecies:
	threads: 100
	input:
		renamedLiftedHints = os.path.join(projectFolder, '{subgenome}', '{subgenome}_lifted_hints.gff3')
	output:
		parametersCFG = os.path.join(environmentFolder, 'config', 'species', '{subgenome}', '{subgenome}_parameters.cfg'),
		weightMatrix = os.path.join(environmentFolder, 'config', 'species', '{subgenome}', '{subgenome}_weightmatrix.txt'),
		metapars = os.path.join(environmentFolder, 'config', 'species', '{subgenome}', '{subgenome}_metapars.cfg'),
		newSpeciesCreatedTxt = os.path.join(projectFolder, '{subgenome}', '{subgenome}_New_Species_Setup.txt')
	shell:'''
		# Augustus script to set up the software for using an organism that isn't prebuilt - start here
		new_species.pl --species={wildcards.subgenome}

		echo "New species setup for {wildcards.subgenome}" > {output.newSpeciesCreatedTxt}
	'''

# this is where we use the provided bonafide.gb file from Thales -- I did not recreate this file
# He got this by following:
# "Predicting Genes in Single Genomes with AUGUSTUS" 
# by Katharina J. Hoff and Mario Stanke

# ********** Steps He Did ************ - Do check this on your own to confirm
# Generating Training Gene Structures From Short Read RNA-seq Data -> Main Protocol (from fastqs and fasta to bonafide.gb)
# Alternate Protocol 2 with steps 3 and 4 of Alternate Protocol 1 - generating training gene structures from ESTs
# Support Protocol 2

# ************** I start at this point in the process using Thales' files ******************
# Alternate Protocol 4 - refining the parameters to train augustus
# first times we run etraining are to adjust parameters which we do in the next few steps
rule initialETraining:
	threads: 20
	input:
		newSpeciesCreatedTxt = os.path.join(projectFolder, '{subgenome}', '{subgenome}_New_Species_Setup.txt'),
		bonafideGB = os.path.join(providedFolder, '{subgenome}', 'bonafide.gb')
	params:
		trainingDir = os.path.join(projectFolder, '{subgenome}', 'training')
	output:
		bonafideOut = os.path.join(projectFolder, '{subgenome}', 'training', '{subgenome}_bonafide.out'),
		cpBonafideGB = os.path.join(projectFolder, '{subgenome}', 'bonafide.gb')
	shell: '''
		mkdir -p {params.trainingDir}
		cp {input.bonafideGB} {output.cpBonafideGB}

		# run initial etraining - to get some parameters we need
		etraining \
			--species={wildcards.subgenome} \
			{output.cpBonafideGB} \
			&> {output.bonafideOut}
	'''

# Refining the bonafide.gb to remove bad training genes
# Thales already did this, but I'm doing this as a check
# Alternate Protocol 4
rule getErrorsForRefineParameters:
	threads: 20
	input:
		bonafideOut = os.path.join(projectFolder, '{subgenome}', 'training', '{subgenome}_bonafide.out'),
		bonafideGB = os.path.join(projectFolder, '{subgenome}', 'bonafide.gb'),
	params:
		parametersCFG = os.path.join(environmentFolder, 'config', 'species', '{subgenome}', '{subgenome}_parameters.cfg'),
		refineParametersScript = os.path.join(projectFolder, 'scripts', 'refineParameters.sh')
	output:
		initialTrainingCheck = os.path.join(projectFolder, '{subgenome}', 'training', '{subgenome}_error_check.txt'),
		badList = os.path.join(projectFolder, '{subgenome}', 'training', 'bad.lst'),
		bonafideFGB = os.path.join(projectFolder, '{subgenome}', 'training', '{subgenome}_bonafide.f.gb')

	shell: '''
		# runs through a few refining steps
		bash {params.refineParametersScript} \
			{wildcards.subgenome} \
			{input.bonafideOut} \
			{input.bonafideGB} \
			{output.initialTrainingCheck} \
			{params.parametersCFG} \
			{output.badList} \
			{output.bonafideFGB} \
	'''

rule createTrainingSet:
	threads: 20
	input:
		initialTrainingCheck = os.path.join(projectFolder, '{subgenome}', 'training', '{subgenome}_error_check.txt')
	params:
		numGenesTest = numGenesTest,
		badList = os.path.join(projectFolder, '{subgenome}', 'training', 'bad.lst'),
		bonafideFGB = os.path.join(projectFolder, '{subgenome}', 'training', '{subgenome}_bonafide.f.gb')
	output:
		testGB = os.path.join(projectFolder, '{subgenome}', 'training', 'test.gb'),
		trainGB = os.path.join(projectFolder, '{subgenome}', 'training', 'train.gb')
	shell: '''
		# randomly splits bonafide file so we have a training set and a testing set
		randomSplit.pl {params.bonafideFGB} {params.numGenesTest}
		mv {params.bonafideFGB}.test {output.testGB}
		mv {params.bonafideFGB}.train {output.trainGB}
	'''

rule eTrain:
	threads: 100
	input:
		testGB = os.path.join(projectFolder, '{subgenome}', 'training', 'test.gb'),
		trainGB = os.path.join(projectFolder, '{subgenome}', 'training', 'train.gb')
	output:
		trainOut = os.path.join(projectFolder, '{subgenome}', 'training', 'train.out')
	shell: '''
		# etraining with the training set
		etraining \
			--species={wildcards.subgenome} \
			{input.trainGB} \
			&> {output.trainOut}
	'''

rule changeConfig:
	threads: 100
	input:
		trainOut = os.path.join(projectFolder, '{subgenome}', 'training', 'train.out'),
	params:
		parametersCFG = os.path.join(environmentFolder, 'config', 'species', '{subgenome}', '{subgenome}_parameters.cfg')
	output:
		paramChangeTxt = os.path.join(projectFolder, '{subgenome}', 'training', '{subgenome}_paramChange.txt')
	shell: '''
		#changing tag percentage
		tail -6 {input.trainOut} | head -1 | awk '{{ print $3 }}' | sed "s/(//g" | sed "s/)//g" > {output.paramChangeTxt}
		grep /Constant/amberprob {params.parametersCFG} >> {output.paramChangeTxt}
		{wildcards.subgenome}tag=$(awk 'FNR==1 {{ print $0 }}' {output.paramChangeTxt})
		{wildcards.subgenome}origCfgTag=$(awk 'FNR==2 {{ print $0 }}' {output.paramChangeTxt})
		echo ${wildcards.subgenome}origCfgTag | sed "s/[0-9].[0-9][0-9]/${wildcards.subgenome}tag/g" >> {output.paramChangeTxt}
		{wildcards.subgenome}cfgTag=$(awk 'FNR==3 {{ print $0 }}' {output.paramChangeTxt})
		sed -i "s|${wildcards.subgenome}origCfgTag|${wildcards.subgenome}cfgTag|g" {params.parametersCFG}

		#changing taa percentage
		tail -5 {input.trainOut} | head -1 | awk '{{ print $3 }}' | sed "s/(//g" | sed "s/)//g" >> {output.paramChangeTxt}
		grep /Constant/ochreprob {params.parametersCFG} >> {output.paramChangeTxt}
		{wildcards.subgenome}taa=$(awk 'FNR==4 {{ print $0 }}' {output.paramChangeTxt})
		{wildcards.subgenome}origCfgTaa=$(awk 'FNR==5 {{ print $0 }}' {output.paramChangeTxt})
		echo ${wildcards.subgenome}origCfgTaa | sed "s/[0-9].[0-9][0-9]/${wildcards.subgenome}taa/g" >> {output.paramChangeTxt}
		{wildcards.subgenome}cfgTaa=$(awk 'FNR==6 {{ print $0 }}' {output.paramChangeTxt})
		sed -i "s|${wildcards.subgenome}origCfgTaa|${wildcards.subgenome}cfgTaa|g" {params.parametersCFG}

		#changing tga percentage
		tail -4 {input.trainOut} | head -1 | awk '{{ print $3 }}' | sed "s/(//g" | sed "s/)//g" >> {output.paramChangeTxt}
		grep /Constant/opalprob {params.parametersCFG} >> {output.paramChangeTxt}
		{wildcards.subgenome}tga=$(awk 'FNR==7 {{ print $0 }}' {output.paramChangeTxt})
		{wildcards.subgenome}origCfgTga=$(awk 'FNR==8 {{ print $0 }}' {output.paramChangeTxt})
		echo ${wildcards.subgenome}origCfgTga | sed "s/[0-9].[0-9][0-9]/${wildcards.subgenome}tga/g" >> {output.paramChangeTxt}
		{wildcards.subgenome}cfgTga=$(awk 'FNR==9 {{ print $0 }}' {output.paramChangeTxt})
		sed -i "s|${wildcards.subgenome}origCfgTga|${wildcards.subgenome}cfgTga|g" {params.parametersCFG}
	'''

rule runTestAugustus:
	threads: 100
	input:
		testGB = os.path.join(projectFolder, '{subgenome}', 'training', 'test.gb'),
		paramChangeTxt = os.path.join(projectFolder, '{subgenome}', 'training', '{subgenome}_paramChange.txt')
	output:
		testOut = os.path.join(projectFolder, '{subgenome}', 'training', 'test.out')
	shell: '''
		# running augustus with new parameters
		augustus \
			--species={wildcards.subgenome} \
			{input.testGB} \
			&> {output.testOut}
	'''

rule optimizeUTRs:
	threads: 100
	input:
		testOut = os.path.join(projectFolder, '{subgenome}', 'training', 'test.out'),
		trainGB = os.path.join(projectFolder, '{subgenome}', 'training', 'train.gb')
	params:
		augustusConfig = os.path.join(environmentFolder, 'config'),
		UTRmetapars = os.path.join(environmentFolder, 'config', 'species', '{subgenome}', '{subgenome}_metapars.utr.cfg')
	output:
		optimizeReport = os.path.join(projectFolder, '{subgenome}', 'training', 'optimizeReport.UTR.out')
	shell: '''
		# refine augustus parameters for UTRs
		optimize_augustus.pl  \
			--cpus={threads} \
			--species={wildcards.subgenome} \
			--kfold={threads} \
			{input.trainGB} \
			--trainOnlyUtr=1 \
			--AUGUSTUS_CONFIG_PATH={params.augustusConfig} \
			--metapars={params.UTRmetapars} \
			| tee {output.optimizeReport}
	'''

rule testUTR:
	threads: 100
	input:
		testGB = os.path.join(projectFolder, '{subgenome}', 'training', 'test.gb'),
		optimizeReport = os.path.join(projectFolder, '{subgenome}', 'training', 'optimizeReport.UTR.out')
	output:
		testUTRFinalOut = os.path.join(projectFolder, '{subgenome}', 'training', 'test.UTR.Final.out')
	shell: '''
		# testing with UTRs on to see if accuracy or precision got better
		augustus \
			--species={wildcards.subgenome} \
			{input.testGB} \
			--UTR=on \
			--print_utr=on \
			> {output.testUTRFinalOut}
	'''

rule optimize:
	threads: 100
	input:
		testUTROut = os.path.join(projectFolder, '{subgenome}', 'training', 'test.UTR.Final.out'),
		trainGB = os.path.join(projectFolder, '{subgenome}', 'training', 'train.gb')
	params:
		augustusConfig = os.path.join(environmentFolder, 'config')
	output:
		optimizeReport = os.path.join(projectFolder, '{subgenome}', 'training', 'optimizeReport.out')
	shell: '''
		# refine augustus parameters for everything but UTRS
		optimize_augustus.pl \
			--cpus={threads} \
			--species={wildcards.subgenome} \
			--kfold={threads} \
			{input.trainGB} \
			--trainOnlyUtr=0 \
			--AUGUSTUS_CONFIG_PATH={params.augustusConfig} \
			| tee {output.optimizeReport}
	'''

rule testFinal:
	threads: 100
	input:
		testGB = os.path.join(projectFolder, '{subgenome}', 'training', 'test.gb'),
		optimizeReport = os.path.join(projectFolder, '{subgenome}', 'training', 'optimizeReport.out')
	output:
		testFinalOut = os.path.join(projectFolder, '{subgenome}', 'training', 'test.noUTR.Final.out')
	shell: '''
		# final testing run of augusutus
		augustus \
			--species={wildcards.subgenome} \
			{input.testGB} \
			--UTR=off \
			> {output.testFinalOut}
	'''

rule trainingTestResults:
	threads: 100
	input:
		testOut = os.path.join(projectFolder, '{subgenome}', 'training', 'test.out'),
		testFinalOut = os.path.join(projectFolder, '{subgenome}', 'training', 'test.noUTR.Final.out'),
		testUTRFinalOut = os.path.join(projectFolder, '{subgenome}', 'training', 'test.UTR.Final.out')
	output:
		trainingSummaryReport = os.path.join(projectFolder, '{subgenome}', 'training', 'trainingSummaryReport.txt')
	shell: '''
		# taking results from the tests and putting them into the same file
		echo "******************* {wildcards.subgenome} ************************" > {output.trainingSummaryReport}
		echo "*******************************************************" >> {output.trainingSummaryReport}
		echo "************** Original Training Output ***************" >> {output.trainingSummaryReport}
		echo "*******************************************************" >> {output.trainingSummaryReport}
		grep -A 40 "Evaluation of gene prediction" {input.testOut} >> {output.trainingSummaryReport}

		echo "*******************************************************" >> {output.trainingSummaryReport}
		echo "************** Final Training Output ******************" >> {output.trainingSummaryReport}
		echo "*******************************************************" >> {output.trainingSummaryReport}
		grep -A 40 "Evaluation of gene prediction" {input.testFinalOut} >> {output.trainingSummaryReport}

		echo "*******************************************************" >> {output.trainingSummaryReport}
		echo "************ Final UTR Training Output ****************" >> {output.trainingSummaryReport}
		echo "*******************************************************" >> {output.trainingSummaryReport}
		grep -A 40 "Evaluation of gene prediction" {input.testUTRFinalOut} >> {output.trainingSummaryReport}
	'''

rule runMainAugustusAnalysisUTR:
	threads: 100
	input:
		trainingSummaryReport = os.path.join(projectFolder, '{subgenome}', 'training', 'trainingSummaryReport.txt'),
		subGenomeHifiFinal = os.path.join(projectFolder, '{subgenome}', '{subgenome}_' + hifiGenomeShorthand + '_final.fa'),
		renamedLiftedHints = os.path.join(projectFolder, '{subgenome}', '{subgenome}_lifted_hints.gff3')
	params:
		extCFG= os.path.join(environmentFolder, 'config', 'extrinsic', 'extrinsic.M.RM.E.W.P.cfg')
	output:
		cpExtCFG = os.path.join(projectFolder, '{subgenome}', 'training', '{subgenome}.UTR.extrinsic.M.RM.E.W.P.cfg'),
		errFile = os.path.join(projectFolder, '{subgenome}', 'training', 'errors_{subgenome}.UTR.err'),
		predictedGFF = os.path.join(projectFolder, '{subgenome}', 'training', 'predicted.UTR.C.{subgenome}.gff')
	shell: '''
		# copying extrinsic augustus file - not sure why we do this
		cp {params.extCFG} {output.cpExtCFG}

		# running Augustus with hints on - official run
		augustus \
			--species={wildcards.subgenome} \
			--UTR=on \
			--extrinsicCfgFile={output.cpExtCFG} \
			--allow_hinted_splicesites=atac \
			{input.subGenomeHifiFinal} \
			--codingseq=on \
			--protein=on \
			--outfile={output.predictedGFF} \
			--progress=true \
			--genemodel=complete \
			--hintsfile={input.renamedLiftedHints} \
			--errfile={output.errFile}
	'''

rule runMainAugustusAnalysisNOUTR:
	threads: 100
	input:
		trainingSummaryReport = os.path.join(projectFolder, '{subgenome}', 'training', 'trainingSummaryReport.txt'),
		subGenomeHifiFinal = os.path.join(projectFolder, '{subgenome}', '{subgenome}_' + hifiGenomeShorthand + '_final.fa'),
		renamedLiftedHints = os.path.join(projectFolder, '{subgenome}', '{subgenome}_lifted_hints.gff3')
	params:
		extCFG= os.path.join(environmentFolder, 'config', 'extrinsic', 'extrinsic.M.RM.E.W.P.cfg')
	output:
		cpExtCFG = os.path.join(projectFolder, '{subgenome}', 'training', '{subgenome}.extrinsic.M.RM.E.W.P.cfg'),
		errFile = os.path.join(projectFolder, '{subgenome}', 'training', 'errors_{subgenome}.err'),
		predictedGFF = os.path.join(projectFolder, '{subgenome}', 'training', 'predicted.noUTR.C.{subgenome}.gff')
	shell: '''
		# copying extrinsic augustus file - not sure why we do this
		cp {params.extCFG} {output.cpExtCFG}

		# running Augustus with hints on - official run
		augustus \
			--species={wildcards.subgenome} \
			--UTR=off \
			--extrinsicCfgFile={output.cpExtCFG} \
			--allow_hinted_splicesites=atac \
			{input.subGenomeHifiFinal} \
			--codingseq=on \
			--protein=on \
			--outfile={output.predictedGFF} \
			--progress=true \
			--genemodel=complete \
			--hintsfile={input.renamedLiftedHints} \
			--errfile={output.errFile}
	'''

rule getAnnoFasta:
	input:
		subGenomeHifiFinal = os.path.join(projectFolder, '{subgenome}', '{subgenome}_' + hifiGenomeShorthand + '_final.fa'),
		predictedGFF = os.path.join(projectFolder, '{subgenome}', 'training', 'predicted.noUTR.C.{subgenome}.gff')
	output:
		predictedAA = os.path.join(projectFolder, '{subgenome}', 'training', 'predicted.noUTR.C.{subgenome}.aa')
	shell: '''
		# get protein annotations
		getAnnoFasta.pl \
			--seqfile {input.subGenomeHifiFinal} \
			{input.predictedGFF}
	'''

rule renamingPredictedAA:
	input:
		predictedAA = os.path.join(projectFolder, '{subgenome}', 'training', 'predicted.noUTR.C.{subgenome}.aa')
	output:
		renamedPredictedAA = os.path.join(projectFolder, '{subgenome}', 'training', '{subgenome}.noUTR.Proteins.aa')
	shell: '''
		# adding subgenome name to protein prediction
		sed "s/>/>{wildcards.subgenome}_/g" {input.predictedAA} > {output.renamedPredictedAA}
	'''

rule BUSCOproteinAnalyis:
	threads: 20
	input:
		renamedPredictedAA = os.path.join(projectFolder, '{subgenome}', 'training', '{subgenome}.noUTR.Proteins.aa')
	params:
		buscoDir = os.path.join(projectFolder, '{subgenome}busco')
	output:
		outputSummary = os.path.join(projectFolder, '{subgenome}busco', 'short_summary.specific.eudicots_odb10.{subgenome}busco.txt'),
		full_table = os.path.join(projectFolder, '{subgenome}busco', 'run_eudicots_odb10', 'full_table.tsv')
	shell: '''
		mkdir -p {params.buscoDir}

		busco \
			-f \
			-i {input.renamedPredictedAA} \
			-l eudicots_odb10 \
			-o {wildcards.subgenome}busco \
			-m prot \
			-c {threads} \
			--long
		'''

rule parseBUSCO:
	input:
		outputSummary = os.path.join(projectFolder, '{subgenome}busco', 'short_summary.specific.eudicots_odb10.{subgenome}busco.txt'),
		full_table = os.path.join(projectFolder, '{subgenome}busco', 'run_eudicots_odb10', 'full_table.tsv')
	params:
		buscoDir = os.path.join(projectFolder, '{subgenome}busco',)
	output:
		buscoCSV = os.path.join(projectFolder, '{subgenome}busco', 'busco.csv')
	shell: '''
		echo "Strain,Complete_single_copy,Complete_duplicated,Fragmented,Missing" > {output.buscoCSV}
		cat {input.outputSummary} | grep "(S)" | awk -v strain="{wildcards.subgenome}" '{{print strain","$1}}' > {params.buscoDir}/complete_single.txt
		cat {input.outputSummary} | grep "(D)" | awk '{{print $1}}' > {params.buscoDir}/complete_duplicated.txt
		cat {input.outputSummary} | grep "(F)" | awk '{{print $1}}' > {params.buscoDir}/fragmented.txt
		cat {input.outputSummary} | grep "(M)" | awk '{{print $1}}' > {params.buscoDir}/missing.txt
		paste -d "," {params.buscoDir}/complete_single.txt {params.buscoDir}/complete_duplicated.txt {params.buscoDir}/fragmented.txt {params.buscoDir}/missing.txt >> {output.buscoCSV}
		rm {params.buscoDir}/complete_single.txt {params.buscoDir}/complete_duplicated.txt {params.buscoDir}/fragmented.txt {params.buscoDir}/missing.txt
	'''

rule combineBUSCO:
	input:
		buscoCSV = expand(os.path.join(projectFolder, '{subgenome}busco', 'busco.csv'), subgenome = subspecies)
	output:
		combinedBuscoCSV = os.path.join(projectFolder, 'combinedBusco.csv')
	shell: '''
		echo "Strain,Complete_single_copy,Complete_duplicated,Fragmented,Missing" > {output.combinedBuscoCSV}
		cat {input.buscoCSV} | sort | uniq -u >> {output.combinedBuscoCSV}
	'''

rule buscoRscript:
	input:
		combinedBuscoCSV = os.path.join(projectFolder, 'combinedBusco.csv')
	params:
		buscoRscript = os.path.join(projectFolder, 'scripts', 'busco.R')
	output:
		buscoPDF = os.path.join(projectFolder, 'busco.pdf')
	shell: '''
		Rscript {params.buscoRscript} \
			{input.combinedBuscoCSV} \
			{output.buscoPDF}
	'''

rule getAnnoFastaUTR:
	input:
		subGenomeHifiFinal = os.path.join(projectFolder, '{subgenome}', '{subgenome}_' + hifiGenomeShorthand + '_final.fa'),
		predictedGFF = os.path.join(projectFolder, '{subgenome}', 'training', 'predicted.UTR.C.{subgenome}.gff')
	output:
		predictedAA = os.path.join(projectFolder, '{subgenome}', 'training', 'predicted.UTR.C.{subgenome}.aa')
	shell: '''
		# get protein annotations
		getAnnoFasta.pl \
			--seqfile {input.subGenomeHifiFinal} \
			{input.predictedGFF}
	'''

rule renamingPredictedAAUTR:
	input:
		predictedAA = os.path.join(projectFolder, '{subgenome}', 'training', 'predicted.UTR.C.{subgenome}.aa')
	output:
		renamedPredictedAA = os.path.join(projectFolder, '{subgenome}', 'training', '{subgenome}.UTR.Proteins.aa')
	shell: '''
		# adding subgenome name to protein prediction
		sed "s/>/>{wildcards.subgenome}_/g" {input.predictedAA} > {output.renamedPredictedAA}
	'''


rule BUSCOproteinAnalyisUTR:
	threads: 20
	input:
		renamedPredictedAA = os.path.join(projectFolder, '{subgenome}', 'training', '{subgenome}.UTR.Proteins.aa')
	params:
		buscoDir = os.path.join(projectFolder, '{subgenome}buscoUTR')
	output:
		outputSummary = os.path.join(projectFolder, '{subgenome}buscoUTR', 'short_summary.specific.eudicots_odb10.{subgenome}busco.txt'),
		full_table = os.path.join(projectFolder, '{subgenome}buscoUTR', 'run_eudicots_odb10', 'full_table.tsv')
	shell: '''
		mkdir -p {params.buscoDir}

		busco \
			-f \
			-i {input.renamedPredictedAA} \
			-l eudicots_odb10 \
			-o {wildcards.subgenome}buscoUTR \
			-m prot \
			-c {threads} \
			--long
		'''

rule parseBUSCOUTR:
	input:
		outputSummary = os.path.join(projectFolder, '{subgenome}buscoUTR', 'short_summary.specific.eudicots_odb10.{subgenome}busco.txt'),
		full_table = os.path.join(projectFolder, '{subgenome}buscoUTR', 'run_eudicots_odb10', 'full_table.tsv')
	params:
		buscoDir = os.path.join(projectFolder, '{subgenome}buscoUTR',)
	output:
		buscoCSV = os.path.join(projectFolder, '{subgenome}buscoUTR', 'busco.csv')
	shell: '''
		echo "Strain,Complete_single_copy,Complete_duplicated,Fragmented,Missing" > {output.buscoCSV}
		cat {input.outputSummary} | grep "(S)" | awk -v strain="{wildcards.subgenome}" '{{print strain","$1}}' > {params.buscoDir}/complete_single.txt
		cat {input.outputSummary} | grep "(D)" | awk '{{print $1}}' > {params.buscoDir}/complete_duplicated.txt
		cat {input.outputSummary} | grep "(F)" | awk '{{print $1}}' > {params.buscoDir}/fragmented.txt
		cat {input.outputSummary} | grep "(M)" | awk '{{print $1}}' > {params.buscoDir}/missing.txt
		paste -d "," {params.buscoDir}/complete_single.txt {params.buscoDir}/complete_duplicated.txt {params.buscoDir}/fragmented.txt {params.buscoDir}/missing.txt >> {output.buscoCSV}
		rm {params.buscoDir}/complete_single.txt {params.buscoDir}/complete_duplicated.txt {params.buscoDir}/fragmented.txt {params.buscoDir}/missing.txt
	'''

rule combineBUSCOUTR:
	input:
		buscoCSV = expand(os.path.join(projectFolder, '{subgenome}buscoUTR', 'busco.csv'), subgenome = subspecies)
	output:
		combinedBuscoCSV = os.path.join(projectFolder, 'combinedBuscoUTR.csv')
	shell: '''
		echo "Strain,Complete_single_copy,Complete_duplicated,Fragmented,Missing" > {output.combinedBuscoCSV}
		cat {input.buscoCSV} | sort | uniq -u >> {output.combinedBuscoCSV}
	'''

rule buscoRscriptUTR:
	input:
		combinedBuscoCSV = os.path.join(projectFolder, 'combinedBuscoUTR.csv')
	params:
		buscoRscript = os.path.join(projectFolder, 'scripts', 'busco.R')
	output:
		buscoPDF = os.path.join(projectFolder, 'buscoUTR.pdf')
	shell: '''
		Rscript {params.buscoRscript} \
			{input.combinedBuscoCSV} \
			{output.buscoPDF}
	'''