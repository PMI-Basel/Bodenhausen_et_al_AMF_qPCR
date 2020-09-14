#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --mem=20g
#SBATCH --output=ASV.out
#SBATCH --error=ASV.error
#SBATCH --job-name=ASV_clustering
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=jan.waelchli@students.unibe.ch
#SBATCH --mail-type=ALL

#load modules
module add UHTS/Quality_control/fastqc/0.11.7
module add UHTS/Quality_control/cutadapt/2.5

 #--------------------------------------------------------------------
 #Jan Wälchli | 06.11.2019
 #--------------------------------------------------------------------

#running time notification
echo 'Start script'

#get the runs
runs=$(awk '{print $3}' design.tab | uniq)

## --------------------------------------------------------------------
## A | Quality Control - FastQC
## --------------------------------------------------------------------

#create a folder
rm -r  ../../4_output/qc 2> /dev/null #suppress error message
mkdir ../../4_output/qc

#quality control
fastqc -t 20 -k 0 -q ../../2_data/* -o ../../4_output/qc 2> /dev/null

#remove no longer needed files
rm ../../4_output/qc/*.zip

#running time notification
echo 'A - Quality Control done'


## --------------------------------------------------------------------
## B | Primer Files
## --------------------------------------------------------------------

#create a folder
rm -r  ./primer 2> /dev/null
mkdir ./primer

#create the primer files
for run in ${runs}; do

	#padding sequence to remove
	padding=$(grep ${run} design.tab | awk '{print $9}' | uniq)
	if [ ${padding} = 'F' ]; then padding=''; fi

	#forward primers
	grep ${run} design.tab | awk '{print $4, $6}' > ./primer/${run}_F_primers.txt
	cat ./primer/${run}_F_primers.txt | \
	tr ' ' '\n' | tr '_' '-' | \
	sed 's/^/\^/g' | \
	sed 's'/'^\^F'/'>F'/'g' | \
	sed 's'/'^'${padding}/''/'g' > ./primer/${run}_F_primers.fasta
	
	#reverse primers
	grep ${run} design.tab | awk '{print $5, $7}' > ./primer/${run}_R_primers.txt
	cat ./primer/${run}_R_primers.txt | \
	tr ' ' '\n' | tr '_' '-' | \
	sed 's'/'^R'/'>R'/'g' | \
	sed 's'/'^'${padding}/''/'g' > ./primer/${run}_R_primers.fasta
	
	
done

#remove no longer needed files
rm ./primer/*.txt

#running time notification
echo 'B - Primer Files done'


## ---------------------------------------------------------------------
## C | Primer Spliting
## ---------------------------------------------------------------------

for run in ${runs}; do
	cutadapt \
		-e 0.01 --no-indels \
		-g file:./primer/${run}_F_primers.fasta \
		-G file:./primer/${run}_R_primers.fasta \
		-o ./primer/${run}_{name2}-{name1}-r1_.fastq \
		-p ./primer/${run}_{name2}-{name1}-r2_.fastq \
		../../2_data/${run}_r1.fastq.gz ../../2_data/${run}_r2.fastq.gz \
		--discard-untrimmed
done

#running time notification
echo 'C - Primer Spliting done'


## ---------------------------------------------------------------------
## D | Clean up
## ---------------------------------------------------------------------

#mv primer lists
mkdir ./primer/primer_lists
mv ./primer/*.fasta ./primer/primer_lists

#move files with primer combinations not in the desing file

#forward primer files
awk '{print $3,$5,$4}' design.tab | \
tr ' _' '-' | \
sed 's/-R/_R/g' | \
sed 's/$/-r1_.fastq/g' | \
uniq > files_to_keep.txt

#reverse primer files
awk '{print $3,$5,$4}' design.tab | \
tr ' _' '-' | \
sed 's/-R/_R/g' | \
sed 's/$/-r2_.fastq/g' | \
uniq >> files_to_keep.txt

#files to move
touch files_to_move.txt
cd primer
mkdir ./unused_files
for i in *.fastq; do
	if ! grep -qxFe "$i" ../files_to_keep.txt; then
		echo "$i" >> ../files_to_move.txt
	fi
done

#move the files
while read line; do
	mv ${line} unused_files
done < ../files_to_move.txt
cd ..

#sort files by run
for run in ${runs}; do
	rm -r ./primer/${run} 2> /dev/null
	mkdir ./primer/${run}
	mv ./primer/${run}*.fastq ./primer/${run}
done

#sort runs by taxa
rm -r ./primer/bacteria ./primer/fungi 2> /dev/null
mkdir ./primer/bacteria ./primer/fungi

for run in ${runs}; do
	taxa=$(grep -m 1 ${run} design.tab | awk '{print $2}')
	if [ ${taxa} = 'b' ]; then mv ./primer/${run} ./primer/bacteria/${run}; fi
	if [ ${taxa} = 'f' ]; then mv ./primer/${run} ./primer/fungi/${run}; fi
done

#running time notification
echo 'D - Clean up done'
echo 'End script'
