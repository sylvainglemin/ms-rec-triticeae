#!bin/sh

# Create the list of focal species fasta files ####
# To be done in the folder with multifasta files (one per species) with all reference sequences
ls > focal_species_fastafile_list.txt



# Create the databases ####

# Hordeum database
makeblastdb -in Hordeum_vulgare.Hv_IBSC_PGSB_v2.cdna.all.fa -parse_seqids -dbtype nucl
# Focal species database
for species in $(cat focal_species_fastafile_list.txt); do makeblastdb -in ${species} -parse_seqids -dbtype nucl;done


# Blast sequences ####

# Note that for focal2hordeum we use the qlen option but for the hordeum2focal we use the slen
# The reason is that the focal corresponds to gene fragments
for species in $(cat focal_species_fastafile_list.txt);
do echo $species;
output=${species}_focal2hordeum;
blastn -db Hordeum_vulgare.Hv_IBSC_PGSB_v2.cdna.all.fa -query ${species} -out ${output}.temp -evalue 1e-10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen gaps";
# Sort blast output for the R script to sort and clean blast outputs
sort -k1,1 -k2,2 -k7,7n ${output}.temp > ${output}.blast.txt;
rm ${output}.temp;
output=${species}_hordeum2focal;	
blastn -db ${species} -query Hordeum_vulgare.Hv_IBSC_PGSB_v2.cdna.all.fa -out ${output}.temp -evalue 1e-10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen gaps";
sort -k1,1 -k2,2 -k7,7n ${output}.temp > ${output}.blast.txt;
rm ${output}.temp;
done;
