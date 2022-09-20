Aim: To create a subset of a blast database like NCBI's nt, so that the pipeline can be tested quickly

Steps:
1. Gather fasta sequences for test assemblies
2. Run blastn for each test assembly against nt
3. Gather hits within nt and create reduced blastdbs using taxids of hits


## 1. Gather fasta sequences for test assemblies
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/910/594/005/GCA_910594005.1_mCerEla1.1/GCA_910594005.1_mCerEla1.1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/922/990/625/GCA_922990625.1_mMelMel3.1_maternal_haplotype_uncurated/GCA_922990625.1_mMelMel3.1_maternal_haplotype_uncurated_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/905/404/285/GCA_905404285.1_ilEraDefo1.1/GCA_905404285.1_ilEraDefo1.1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/927/399/515/GCA_927399515.1_gfLaeSulp1.1/GCA_927399515.1_gfLaeSulp1.1_genomic.fna.gz
```
## 2. Run blastn for each test assembly against nt
```
for a in *gz
do 
  blastn -task megablast \
    -query <(gunzip -c $a) \
    -db /path/to/blastdb/nt \
    -outfmt "6 qseqid staxids bitscore std" \
    -max_target_seqs 10 \
    -max_hsps 1 \
    -evalue 1.0e-10 \
    -num_threads 8 \
    -lcase_masking \
    -dust "20 64 1" \
  > $a.blastn.nt.out.raw
done
```
## 3. Gather hits within nt and create reduced blastdbs using taxids of hits
```
for a in mCerEla1.1 mMelMel3.1 ilEraDefo1.1 gfLaeSulp1.1
do
  mkdir nt_$a
  cut -f5 *_${a}_*blastn.nt.out.raw | sort | uniq > $a.hits.txt
  blastdbcmd -entry_batch $a.hits.txt -db /path/to/blastdb/nt > $a.hits.fasta
  blastdbcmd -entry_batch $a.hits.txt -db /path/to/blastdb/nt -outfmt '%a %T' > $a.hits.taxidmap
  makeblastdb -dbtype nucl -in $a.hits.fasta -taxid_map $a.hits.taxidmap -out nt_$a/nt_$a -parse_seqids
done
```
