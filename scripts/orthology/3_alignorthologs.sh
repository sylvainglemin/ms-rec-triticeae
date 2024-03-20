#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=2G

# Script that aligns the sequence of a focal species with its ortholog in Hordeum
# then compute the dS between the two sequences


############ Input ##########
BestBlast=$1
Spe=$2
OrthoSpe=$3
OrthoCDS=$4
ContigSeq=$5
OutDir=$6
OutDir=${OutDir}/${Spe}
OutDirDS=$7

nbLine=$(wc -l $BestBlast | cut -d' ' -f1)
#############################

## create tmp directory
if [ ! -d "tmp_blastbest_aln" ]; then
        mkdir tmp_blastbest_aln
fi

## Create output directory
if [ ! -d "$OutDir" ]; then
        mkdir $OutDir
fi


source /softs/local/env/envblast-2.7.1.sh

## For each contigs of the species

for i in $(seq 2 $nbLine)
do

   line=$(sed -n ${i}p $BestBlast)
   #get contig name
   Contig=$(echo $line | perl -lane 'print "$F[0]"')
   echo $Contig
   #get ortho gene
   Ortho=$(echo $line | perl -lane 'print "$F[1]"')
   echo $Ortho

   #get sequence of contig
   grep -w -A1 $Contig $ContigSeq > tmp_blastbest_aln/${Spe}_${Contig}_seq.fasta

   #get sequence of ortho gene
   if [ ! -f "tmp_blastbest_aln/${Ortho}_${OrthoSpe}_seq.fasta" ]; then
      grep -A1 $Ortho $OrthoCDS > tmp_blastbest_aln/${Ortho}_${OrthoSpe}_seq.fasta
   fi

   #BlastN
   source /local/env/envblast-2.9.0.sh
   blastn -query tmp_blastbest_aln/${Spe}_${Contig}_seq.fasta \
   -subject tmp_blastbest_aln/${Ortho}_${OrthoSpe}_seq.fasta \
   -out tmp_blastbest_aln/${Spe}_${Contig}_blast_${OrthoSpe}.txt \
   -outfmt '6 qseqid sseqid qlen slen qstart qend sstart send evalue length pident' \
   -max_hsps 4

   #get best hit CDS for aln
   BestHitLine=$(sed -n 1p tmp_blastbest_aln/${Spe}_${Contig}_blast_${OrthoSpe}.txt)
   BestHit=$(echo $BestHitLine | perl -lane 'print "$F[1]"')
   grep -A1 "$BestHit " tmp_blastbest_aln/${Ortho}_${OrthoSpe}_seq.fasta >> tmp_blastbest_aln/${Spe}_${Contig}_seq.fasta

   #Align with MACSE
   source /local/env/envjava-1.8.0.sh
   java -jar macse_v2.03.jar -prog alignSequences -seq tmp_blastbest_aln/${Spe}_${Contig}_seq.fasta \
   -out_NT ${OutDir}/${Spe}_${Contig}_${OrthoSpe}_${Ortho}_alnNT.fasta \
   -out_AA ${OutDir}/${Spe}_${Contig}_${OrthoSpe}_${Ortho}_alnAA.fasta
   
   ## Modif aln file
   ALN=${OutDir}/${Spe}_${Contig}_${OrthoSpe}_${Ortho}_alnNT.fasta
   sed -i "/${BestHit}/c >${BestHit}" $ALN
   sed -i 's/!/-/g' $ALN

   ##Calcul ds
   #Prep codeml config file
   echo "seqfile = ${ALN}" > tmp_blastbest_aln/${Spe}_${Contig}_codeml.ctl
   echo "outfile = tmp_blastbest_aln/${Spe}_${Contig}_codeml_out.txt" >> tmp_blastbest_aln/${Spe}_${Contig}_codeml.ctl
   cat codeml.txt >> tmp_blastbest_aln/${Spe}_${Contig}_codeml.ctl
   
   source /local/env/envpaml-4.9.sh
   codeml tmp_blastbest_aln/${Spe}_${Contig}_codeml.ctl

   DSline=$(grep "dS =" tmp_blastbest_aln/${Spe}_${Contig}_codeml_out.txt)
   DS=$(echo $DSline | cut -d'=' -f7)

   echo $Contig $Ortho $DS $DSline >> ${OutDirDS}/${Spe}_${OrthoSpe}_ds.txt

   ##Remove tmp files
   #rm tmp_blastbest_aln/${Spe}_${Contig}_seq.fasta
   #rm tmp_blastbest_aln/${Spe}_${Contig}_blast_${OrthoSpe}.txt
   rm tmp_blastbest_aln/${Spe}_${Contig}_codeml.ctl
   #rm tmp_blastbest_aln/${Spe}_${Contig}_codeml_out.txt
   rm ${OutDir}/${Spe}_${Contig}_${OrthoSpe}_${Ortho}_alnNT.fasta
   rm ${OutDir}/${Spe}_${Contig}_${OrthoSpe}_${Ortho}_alnNAA.fasta
done