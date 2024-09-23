#!/bin.sh

subgenome1=$1
subgenome2=$2
projectFolder=$3
Rakefile=$4

echo $subgenome1
cd $projectFolder/$subgenome1
rake -f $Rakefile

wait

echo $subgenome2
cd $projectFolder/$subgenome2
rake -f $Rakefile

wait
