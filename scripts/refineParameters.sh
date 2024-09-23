#!/bin/sh

subgenome=$1
bonafideOut=$2
bonafideGB=$3
initialTrainingCheck=$4
parametersCFG=$5
badList=$6
bonafideFGB=$7

errors=$(grep -c "Variable stopCodonExcludedFromCDS set right" $bonafideOut)
totalLoci=$(grep -c LOCUS $bonafideGB)

errorRate=$((($errors/$totalLoci)*100))

echo "Initial Check for Bad Training Genes" > $initialTrainingCheck
echo "# of errors stopCodonExcludedFromCDS:" >> $initialTrainingCheck
echo $errors >> $initialTrainingCheck
echo "total # of loci in bonafide.gb:" >> $initialTrainingCheck
echo $totalLoci >> $initialTrainingCheck

# In Alternate Protocol 4, authors stated that stopCodonExcludedFromCDS 
# did not need to be set to TRUE unless the number of errors 
# wasn't a minority.
# I set this as 20%, but this would be up for debate, a correct percentage wasn't indicated by the authors.
if [ $errorRate -gt 20 ]; then
	echo "Need to refine parameters" >> $initialTrainingCheck
	sed -i 's/stopCodonExcludedFromCDS false/stopCodonExcludedFromCDS true/g' $parametersCFG
	echo "stopCodonExcludedFromCDS is changed to TRUE" >> $initialTrainingCheck
else
	echo "Good to go for stopCodonExcludedFromCDS!" >> $initialTrainingCheck
fi

# 2nd part of refining parameters - already done by Thales, but included here to make sure
etraining --species=$subgenome $bonafideGB 2>&1 \
	| grep "in sequence" \
	| perl -pe 's/.*n sequence (\S+):.*/$1/' \
	| sort -u > $badList
filterGenes.pl $badList $bonafideGB > $bonafideFGB

remainingLoci=$(grep -c LOCUS $bonafideGB $bonafideFGB)
echo $remainingLoci >> $initialTrainingCheck
