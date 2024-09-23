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
		expand(os.path.join(projectFolder, '{subgenome}', 'training', 'predicted.C.{subgenome}.aa'), subgenome = subspecies),

rule downloadGenomes:
	params:
		prevGenomeLink = prevGenomeLink,
		hifiGenomeLink = hifiGenomeLink,
		genomeDir = os.path.join(projectFolder, 'genome')
	output:
		prevGenomeFile = os.path.join(projectFolder, 'genome', prevGenomeGZ),
		hifiGenomeFile = os.path.join(projectFolder, 'genome', hifiGenomeGZ)
	shell: '''
		# make directory
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
		lines=$(cat {input.chrConversion} | wc -l)
		line=1
		cp {input.genomeCleaned} temp_{wildcards.genome}_chr.fa

		while [ $line -le $lines ]
			do
				oldChrName=$(awk '{{ print $1 }}' {input.chrConversion} | awk -v var1="$line" 'NR==var1')
				newChrName=$(awk '{{ print $2 }}' {input.chrConversion} | awk -v var1="$line" 'NR==var1')
				sed -i -e "s/$oldChrName/$newChrName/g" temp_{wildcards.genome}_chr.fa
				let line++
			done

		mv temp_{wildcards.genome}_chr.fa {output.genomeConverted}
	'''

rule rmUnknownsGenome:
	input:
		chrConversion = chrConversion,
		genomeConverted = os.path.join(projectFolder, 'genome', '{genome}_converted.fa')
	output:
		genomeFinal = os.path.join(projectFolder, 'genome', '{genome}_final.fa')
	shell: '''
		lines=$(cat {input.chrConversion} | wc -l)
		line=1

		while [ $line -le $lines ]
			do
				newChrName=$(awk '{{ print $2 }}' {input.chrConversion} | awk -v var1="$line" 'NR==var1')
				sed -n "/>$newChrName/{{p;n;p}}" {input.genomeConverted} >> temp_{wildcards.genome}_rep.fa
				let line++
			done

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
		mkdir -p {params.subgenomeDir}

		chrName=$(grep {wildcards.subgenome} {input.subspeciesConversion} | awk '{{ print $1 }}')

		grep -A 1 "sg$chrName" {params.genomeDir}/{wildcards.genome}_final.fa > {output.subGenomeFinal}

	'''

# Nexts we have to combine the hints files that Thales generated 
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
		wget -c https://github.com/yeban/flo/archive/master.tar.gz -O flo.tar.gz
	'''

rule gunzipFlo:
	input:
		os.path.join(projectFolder, 'flo.tar.gz')
	output:
		expand(os.path.join(projectFolder, 'flo-master', '{flo}'),flo=FLO)
	shell: '''
		tar xvf flo.tar.gz
		wait
	'''

rule moveFlo:
	input:
		expand(os.path.join(projectFolder, 'flo-master', '{flo}'),flo=FLO)
	output:
		expand(os.path.join(projectFolder, 'flo', '{flo}'),flo=FLO)
	shell: '''
		mkdir -p flo

		mv flo-master/* flo/

		rmdir flo-master
	'''

rule installFlo:
	input:
		expand(os.path.join(projectFolder, 'flo', '{flo}'),flo=FLO),
		installFloScript = os.path.join(projectFolder, 'flo', 'scripts', 'install.sh')
	output:
		floInstallDoneTxt = os.path.join(projectFolder, 'floInstallDone.txt')
	shell: '''
		bash {input.installFloScript}

		echo Flo Install Done! > {output.floInstallDoneTxt}
	'''

# this step was something the flo instructions said to do - removing genes from the hints files I got earlier
rule removeGenesFromHintsGFF:
	input:
		convertedHints = os.path.join(projectFolder, '{subgenome}', '{subgenome}_hints_converted.gff'),
		floInstallDoneTxt = os.path.join(projectFolder, 'floInstallDone.txt')
	params:
		floPath = os.path.join(projectFolder, 'flo')
	output:
		subGenomeNoGenesHints = os.path.join(projectFolder, '{subgenome}', '{subgenome}_nogenes_hints.gff')
	shell: '''
		{params.floPath}/gff_remove_feats.rb \
		gene \
		{input.convertedHints} \
		> {output.subGenomeNoGenesHints}
	'''

rule chooseOptionsForFlo:
	input:
		expand(os.path.join(projectFolder, 'flo', '{flo}'),flo=FLO),
		subGenomeSource = os.path.join(projectFolder, '{subgenome}', '{subgenome}_' + prevGenomeShorthand +'_final.fa'),
		subGenomeTarget = os.path.join(projectFolder, '{subgenome}', '{subgenome}_' + hifiGenomeShorthand + '_final.fa'),
		subGenomeNoGenesHints = os.path.join(projectFolder, '{subgenome}', '{subgenome}_nogenes_hints.gff')
	params:
		projectFolder = projectFolder
	output:
		optYaml = os.path.join(projectFolder, '{subgenome}', 'flo_opts.yaml')
	shell: '''
		echo ":add_to_path:" > {output.optYaml}
		echo "  - '{params.projectFolder}/ext/kent/bin'" >> {output.optYaml}
		echo "  - '{params.projectFolder}/ext/parallel-20150722/src'" >> {output.optYaml}
		echo "  - '{params.projectFolder}/ext/genometools-1.5.6/bin'" >> {output.optYaml}
		echo ":source_fa: '{input.subGenomeSource}'" >> {output.optYaml}
		echo ":target_fa: '{input.subGenomeTarget}'" >> {output.optYaml}
		echo ":processes: '20'" >> {output.optYaml}
		echo ":blat_opts: '-fastMap -tileSize=12 -minIdentity=98'"  >> {output.optYaml}
		echo ":lift:" >> {output.optYaml}
		echo "  - '{input.subGenomeNoGenesHints}'" >> {output.optYaml}
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
		bash {params.runFloScript} {params.subgenomes} {params.projectFolder} ../flo/Rakefile

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
		cp {params.liftedHints} {output.renamedLiftedHints}
	'''

rule newSpecies:
	threads: 20
	input:
		renamedLiftedHints = os.path.join(projectFolder, '{subgenome}', '{subgenome}_lifted_hints.gff3')
	output:
		parametersCFG = os.path.join(environmentFolder, 'config', 'species', '{subgenome}', '{subgenome}_parameters.cfg'),
		weightMatrix = os.path.join(environmentFolder, 'config', 'species', '{subgenome}', '{subgenome}_weightmatrix.txt'),
		metapars = os.path.join(environmentFolder, 'config', 'species', '{subgenome}', '{subgenome}_metapars.cfg'),
		newSpeciesCreatedTxt = os.path.join(projectFolder, '{subgenome}', '{subgenome}_New_Species_Setup.txt')
	shell:'''
		# Augustus script to set up the software for using an organism that isn't prebuilt
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
		randomSplit.pl {params.bonafideFGB} {params.numGenesTest}
		mv {params.bonafideFGB}.test {output.testGB}
		mv {params.bonafideFGB}.train {output.trainGB}
	'''

rule eTrain:
	threads: 20
	input:
		testGB = os.path.join(projectFolder, '{subgenome}', 'training', 'test.gb'),
		trainGB = os.path.join(projectFolder, '{subgenome}', 'training', 'train.gb')
	output:
		trainOut = os.path.join(projectFolder, '{subgenome}', 'training', 'train.out')
	shell: '''
		etraining \
			--species={wildcards.subgenome} \
			{input.trainGB} \
			&> {output.trainOut}
	'''

rule changeConfig:
	threads: 20
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
	threads: 20
	input:
		testGB = os.path.join(projectFolder, '{subgenome}', 'training', 'test.gb'),
		paramChangeTxt = os.path.join(projectFolder, '{subgenome}', 'training', '{subgenome}_paramChange.txt')
	output:
		testOut = os.path.join(projectFolder, '{subgenome}', 'training', 'test.out')
	shell: '''
		augustus \
			--species={wildcards.subgenome} \
			{input.testGB} \
			&> {output.testOut}
	'''

rule runUTRTestAugustus:
	threads: 20
	input:
		testGB = os.path.join(projectFolder, '{subgenome}', 'training', 'test.gb'),
		testOut = os.path.join(projectFolder, '{subgenome}', 'training', 'test.out')
	output:
		testOut = os.path.join(projectFolder, '{subgenome}', 'training', 'test.utr.out')
	shell: '''
		augustus \
			--species={wildcards.subgenome} \
			{input.testGB} \
			--UTR=on \
			--print_utr=on \
			&> {output.testOut}
	'''

rule optimizeUTRs:
	threads: 20
	input:
		testOut = os.path.join(projectFolder, '{subgenome}', 'training', 'test.utr.out'),
		trainGB = os.path.join(projectFolder, '{subgenome}', 'training', 'train.gb')
	params:
		augustusConfig = os.path.join(environmentFolder, 'config'),
		UTRmetapars = os.path.join(environmentFolder, 'config', 'species', '{subgenome}', '{subgenome}_metapars.utr.cfg')
	output:
		optimizeReport = os.path.join(projectFolder, '{subgenome}', 'training', 'optimizeReport.UTR.out')
	shell: '''
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

rule trainAfterUTROptimization:
	threads: 20
	input: 
		optimizeReport = os.path.join(projectFolder, '{subgenome}', 'training', 'optimizeReport.UTR.out'),
		trainGB = os.path.join(projectFolder, '{subgenome}', 'training', 'train.gb')
	output:
		finalTrainOut = os.path.join(projectFolder, '{subgenome}', 'training', 'train.UTR.out')
	shell: '''
		etraining \
			--species={wildcards.subgenome} \
			{input.trainGB} \
			&> {output.finalTrainOut}
	'''

rule testUTR:
	threads: 20
	input:
		testGB = os.path.join(projectFolder, '{subgenome}', 'training', 'test.gb'),
		finalTrainOut = os.path.join(projectFolder, '{subgenome}', 'training', 'train.UTR.out')
	output:
		testUTRFinalOut = os.path.join(projectFolder, '{subgenome}', 'training', 'test.UTR.Final.out')
	shell: '''
		augustus \
			--species={wildcards.subgenome} \
			{input.testGB} \
			--UTR=on \
			--print_utr=on \
			> {output.testUTRFinalOut}
	'''

rule optimize:
	threads: 20
	input:
		testUTROut = os.path.join(projectFolder, '{subgenome}', 'training', 'test.UTR.Final.out'),
		trainGB = os.path.join(projectFolder, '{subgenome}', 'training', 'train.gb')
	params:
		augustusConfig = os.path.join(environmentFolder, 'config')
	output:
		optimizeReport = os.path.join(projectFolder, '{subgenome}', 'training', 'optimizeReport.out')
	shell: '''
		optimize_augustus.pl \
			--cpus={threads} \
			--species={wildcards.subgenome} \
			--kfold={threads} \
			{input.trainGB} \
			--trainOnlyUtr=0 \
			--AUGUSTUS_CONFIG_PATH={params.augustusConfig} \
			| tee {output.optimizeReport}
	'''

rule trainAfterOptimization:
	threads: 20
	input: 
		optimizeReport = os.path.join(projectFolder, '{subgenome}', 'training', 'optimizeReport.out'),
		trainGB = os.path.join(projectFolder, '{subgenome}', 'training', 'train.gb')
	output:
		finalTrainOut = os.path.join(projectFolder, '{subgenome}', 'training', 'train.afterOpt.out')
	shell: '''
		etraining \
			--species={wildcards.subgenome} \
			{input.trainGB} \
			&> {output.finalTrainOut}
	'''

rule testFinal:
	threads: 20
	input:
		testGB = os.path.join(projectFolder, '{subgenome}', 'training', 'test.gb'),
		testUTRFinalOut = os.path.join(projectFolder, '{subgenome}', 'training', 'train.afterOpt.out')
	output:
		testFinalOut = os.path.join(projectFolder, '{subgenome}', 'training', 'test.noUTR.Final.out')
	shell: '''
		augustus \
			--species={wildcards.subgenome} \
			{input.testGB} \
			--UTR=off \
			> {output.testFinalOut}
	'''

rule trainingTestResults:
	threads: 20
	input:
		testOut = os.path.join(projectFolder, '{subgenome}', 'training', 'test.out'),
		testUTROut = os.path.join(projectFolder, '{subgenome}', 'training', 'test.utr.out'),
		testFinalOut = os.path.join(projectFolder, '{subgenome}', 'training', 'test.noUTR.Final.out'),
		testUTRFinalOut = os.path.join(projectFolder, '{subgenome}', 'training', 'test.UTR.Final.out')
	output:
		trainingSummaryReport = os.path.join(projectFolder, '{subgenome}', 'training', 'trainingSummaryReport.txt')
	shell: '''
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
		echo "********** Optimized UTR Training Output **************" >> {output.trainingSummaryReport}
		echo "*******************************************************" >> {output.trainingSummaryReport}
		grep -A 40 "Evaluation of gene prediction" {input.testUTROut} >> {output.trainingSummaryReport}

		echo "*******************************************************" >> {output.trainingSummaryReport}
		echo "************ Final UTR Training Output ****************" >> {output.trainingSummaryReport}
		echo "*******************************************************" >> {output.trainingSummaryReport}
		grep -A 40 "Evaluation of gene prediction" {input.testUTRFinalOut} >> {output.trainingSummaryReport}
	'''

rule runMainAugustusAnalysis:
	threads: 20
	input:
		trainingSummaryReport = os.path.join(projectFolder, '{subgenome}', 'training', 'trainingSummaryReport.txt'),
		subGenomeHifiFinal = os.path.join(projectFolder, '{subgenome}', '{subgenome}_' + hifiGenomeShorthand + '_final.fa'),
		renamedLiftedHints = os.path.join(projectFolder, '{subgenome}', '{subgenome}_lifted_hints.gff3')
	params:
		extCFG= os.path.join(environmentFolder, 'config', 'extrinsic', 'extrinsic.M.RM.E.W.P.cfg')
	output:
		cpExtCFG = os.path.join(projectFolder, '{subgenome}', 'training', '{subgenome}.extrinsic.M.RM.E.W.P.cfg'),
		errFile = os.path.join(projectFolder, '{subgenome}', 'training', 'errors_{subgenome}.err'),
		predictedGFF = os.path.join(projectFolder, '{subgenome}', 'training', 'predicted.C.{subgenome}.gff')
	shell: '''
		cp {params.extCFG} {output.cpExtCFG}

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

rule getAnnoFasta:
	input:
		subGenomeHifiFinal = os.path.join(projectFolder, '{subgenome}', '{subgenome}_' + hifiGenomeShorthand + '_final.fa'),
		predictedGFF = os.path.join(projectFolder, '{subgenome}', 'training', 'predicted.C.{subgenome}.gff')
	output:
		predictedAA = os.path.join(projectFolder, '{subgenome}', 'training', 'predicted.C.{subgenome}.aa')
	shell: '''
		getAnnoFasta.pl \
			--seqfile {input.subGenomeHifiFinal} \
			{input.predictedGFF}
	'''

rule renamingPredictedAA:
	input:
		predictedAA = os.path.join(projectFolder, '{subgenome}', 'training', 'predicted.C.{subgenome}.aa')
	output:
		renamedPredictedAA = os.path.join(projectFolder, '{subgenome}', 'training', '{subgenome}.Proteins.aa')
	shell: '''
		sed "s/>/>SubC.{wildcards.subgenome}_/g" {input.predictedAA} > {output.renamedPredictedAA}
	'''
