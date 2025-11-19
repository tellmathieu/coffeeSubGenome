import os, sys, glob
from multiprocessing import cpu_count

#######################
# BRAKER3 Version of Thales thesis project of exploring sub-genome dominance in coffee arabica
# By Tell Mathieu
# 
# Thales paper: https://www.biorxiv.org/content/10.1101/2023.03.08.531780v1.full
# 
# Updated genome from this paper: 
# 	https://www.nature.com/articles/s41588-024-01695-w
# Biological question: Is there subgenome dominance in Coffea arabica? Is one of the subgenomes losing universal orthologs faster than others?
# Collaborators: Thales Henrique Cherubino Ribeiro, Blake Meyers
#######################

#######################
#        Input        #
#######################

srafile = "input/sra.txt"
javaBin = "/quobyte/bcmeyersgrp/tell/.conda/envs/java-1-8/bin"
gushrBin = "/quobyte/bcmeyersgrp/tell/braker/gushr/GUSHR"
genemarkBin = "/quobyte/bcmeyersgrp/tell/braker/GeneMark-ETP/bin"

# also put your subgenome fastas in a 
# folder called "genome" in your project folder

#######################
#     Variables       #
#######################

subgenome = ['subCanephora', 'subEugenioides']
#put subgenome fastas in this folder following correct naming convention
genomeDir = os.path.join(projectFolder, 'genome')
genome = os.path.join(genomeDir, '{subgenome}_Cara_hifi_final.fa')

#######################
#   Pre-SnakeMake     #
#######################

# set working directory, in case we need it
projectFolder=os.getcwd()

# get sra ids from srafile
with open(srafile, 'r') as file:
    a = list(file)

sra = [item.rstrip() for item in a]

sralist = '' + sra[0]
for item in sra[1:]:
	sralist = sralist + ',' + item

### File names

fastqDir = os.path.join(projectFolder, 'fastq')
splitFastqDir = os.path.join(projectFolder, 'splitFastq')
trimmedFastqDir = os.path.join(projectFolder, 'trimmedFastq')
processDir = os.path.join(projectFolder, 'processing')
bamDir = os.path.join(processDir, 'bam')
braker = os.path.join(projectFolder, '{subgenome}','braker_utr.gtf')
aa = os.path.join(projectFolder, '{subgenome}','braker_utr.aa')
indexBuilt = os.path.join(genomeDir, '{subgenome}_index_built.txt')
sam = os.path.join(processDir, '{sra}-{subgenome}.sam')
unsortedbam = os.path.join(processDir, '{sra}-{subgenome}.unsorted.bam')
bam = os.path.join(bamDir, '{sra}-{subgenome}.bam')
buscoDir = os.path.join(projectFolder, '{subgenome}buscoUTR')
outputSummary = os.path.join(projectFolder, '{subgenome}buscoUTR', 'short_summary.specific.eudicots_odb10.{subgenome}buscoUTR.txt'),
full_table = os.path.join(projectFolder, '{subgenome}buscoUTR', 'run_eudicots_odb10', 'full_table.tsv')
buscoCSV = os.path.join(projectFolder, '{subgenome}buscoUTR', 'busco.csv')
combinedBuscoCSV = os.path.join(projectFolder, 'combinedBuscoUTR.csv')
buscoRscript = os.path.join(projectFolder, 'scripts', 'busco.R')
buscoPDF = os.path.join(projectFolder, 'buscoUTR.pdf')

bamlist = '' + sra[0]
for item in sra[1:]:
	bamlist = bamlist + ',' + item


rule all:
	input: 
		buscoPDF,

rule downloadFastq:
	params:
		fastqDir = fastqDir
	output:
		fastq = os.path.join(fastqDir, '{sra}.fastq')
	shell:'''
		mkdir -p {params.fastqDir}

		fastq-dump --outdir {params.fastqDir} --accession {wildcards.sra}
	'''

rule splitFastq:
	input:
		fastq = os.path.join(fastqDir, '{sra}.fastq')
	params:
		splitFastqDir = splitFastqDir
	output:
		r1 = os.path.join(trimmedFastqDir, '{sra}_R1.fastq'),
		r2 = os.path.join(trimmedFastqDir, '{sra}_R2.fastq')
	shell: '''
		mkdir -p {params.splitFastqDir}

		seqtk seq -1 {input.fastq} > {output.r1}
		seqtk seq -2 {input.fastq} > {output.r2}
	'''

rule trimFastq:
	input:
		r1 = os.path.join(trimmedFastqDir, '{sra}_R1.fastq'),
		r2 = os.path.join(trimmedFastqDir, '{sra}_R2.fastq')
	params:
		trimmedFastqDir = trimmedFastqDir
	output:
		trimmed1 = os.path.join(trimmedFastqDir, '{sra}_R1_trimmed.fastq'),
		trimmed2 = os.path.join(trimmedFastqDir, '{sra}_R2_trimmed.fastq')
	shell: '''
		mkdir -p {params.trimmedFastqDir}

		fastp -i {input.r1} -I {input.r2} -o {output.trimmed1} -O {output.trimmed2} --detect_adapter_for_pe
	'''

rule build:
	threads: cpu_count()
	conda: 'envs/hisat2.yaml'
	input:
		genome = genome
	params:
		genomeDir = genomeDir
	output:
		indexBuilt = indexBuilt
	shell: '''
		hisat2-build -p {threads} {input.genome} {params.genomeDir}/{wildcards.subgenome}

		echo Index Built > {output.indexBuilt}
	'''

rule mapBam:
	threads: max(2, cpu_count() // 3)
	conda: 'envs/hisat2.yaml'
	input:
		indexBuilt = indexBuilt,
		read1 = os.path.join(trimmedFastqDir, '{sra}_R1_trimmed.fastq'),
		read2 = os.path.join(trimmedFastqDir, '{sra}_R2_trimmed.fastq')
	params:
		genomeDir = genomeDir,
		processDir = processDir
	output:
		sam = temp(sam)
	shell: '''
		mkdir -p {params.processDir}

		hisat2 -p {threads} -x {params.genomeDir}/{wildcards.subgenome} -1 {input.read1} -2 {input.read2} -S {output.sam}
	'''


rule samToBam:
	threads: max(2, cpu_count() // 3)
	conda: 'envs/samtools.yaml'
	input:
		sam = sam
	output:
		unsortedbam = temp(unsortedbam)
	shell: '''
			samtools view \
				-@ {threads} \
				-Su \
				{input.sam} \
				> {output.unsortedbam}
		'''

rule sortBam:
	threads: max(2, cpu_count() // 3)
	conda: 'envs/samtools.yaml'
	input:
		unsortedbam = unsortedbam
	params:
		bamDir = bamDir
	output:
		bam = bam
	shell: '''
			mkdir -p {params.bamDir}

			samtools sort \
    			-@ {threads} \
    			{input.unsortedbam} \
    			-o {output.bam}
		'''

rule braker:
	threads: cpu_count()
	input:
		genome = genome,
		bam = lambda wildcards: expand(bam, subgenome=wildcards.subgenome, sra=sra)
	params:	
		javaBin = javaBin,
		gushrBin = gushrBin,
		genemarkBin = genemarkBin
	output:
		braker = braker,
		aa = aa
	shell: '''
		mkdir -p {wildcards.subgenome}

		bamFiles=( {input.bam} )
		bamInput=$(IFS=, ; echo "${{bamFiles[@]}}")

		echo $bamInput

		braker.pl --species={wildcards.subgenome} \
			--genome={input.genome} \
			--bam="$bamInput" \
			--UTR=on \
			--threads {threads} \
			--workingdir={wildcards.subgenome} \
			--verbosity=3 \
			--JAVA_PATH={params.javaBin} \
			--GUSHR_PATH={params.gushrBin} \
			--GENEMARK_PATH={params.genemarkBin}
	'''


rule BUSCOproteinAnalyisUTR:
	threads: 20
	input:
		aa = aa
	params:
		buscoDir = buscoDir
	output:
		outputSummary = outputSummary,
		full_table = full_table
	shell: '''
		mkdir -p {params.buscoDir}

		busco \
			-f \
			-i {input.aa} \
			-l eudicots_odb10 \
			-o {wildcards.subgenome}buscoUTR \
			-m prot \
			-c {threads} \
			--long
		'''

rule parseBUSCOUTR:
	input:
		outputSummary = outputSummary,
		full_table = full_table
	params:
		buscoDir = buscoDir
	output:
		buscoCSV = buscoCSV
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
		buscoCSV = expand(buscoCSV, subgenome = subgenome)
	output:
		combinedBuscoCSV = combinedBuscoCSV
	shell: '''
		echo "Strain,Complete_single_copy,Complete_duplicated,Fragmented,Missing" > {output.combinedBuscoCSV}
		cat {input.buscoCSV} | sort | uniq -u >> {output.combinedBuscoCSV}
	'''

rule buscoRscriptUTR:
	input:
		combinedBuscoCSV = combinedBuscoCSV
	params:
		buscoRscript = buscoRscript
	output:
		buscoPDF = buscoPDF
	shell: '''
		Rscript {params.buscoRscript} \
			{input.combinedBuscoCSV} \
			{output.buscoPDF}
	'''
		