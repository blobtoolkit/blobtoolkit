# Uniprot

Based on Uniprot reference proteomes downloaded June 2021. Commands below were run from a directory containing files for the full database.

Install dependencies

```
conda install -c conda-forge -c bioconda diamond seqtk
```

Fetch sequence names from public BTK results

```
mkdir -p reduced && cd reduced

echo "gfLaeSulp1_1 gfLaeSulp1.1
CAJQFK01.1 ilEraDefo1.1
CAJVBL01 mCerEla1.1
CAKLPM01 mMelMel3.1" > assembly.names

while read BTKID TOLID; do
  for TAXRULE in buscogenes buscoregions; do
    curl -s "https://blobtoolkit.genomehubs.org/api/v1/field/${BTKID}/${TAXRULE}_positions" | grep -oP "tr\|[\w|]+" | sort -u > ${TOLID}.${TAXRULE}.ids
  done
done < assembly.names
```

Generate a set of files containing all matching sequences

```
cat *.ids | sort -u > all.ids

seqtk subseq ../reference_proteomes.fasta.gz all.ids | gzip -c > all.fasta.gz

printf "accession\taccession.version\ttaxid\tgi\n" > all.taxid_map

grep -Ff <(cut -d"|" -f 2 all.ids) ../reference_proteomes.taxid_map >> all.taxid_map

diamond makedb -p 4 --in all.fasta.gz --taxonmap all.taxid_map --taxonnodes ../../taxdump_2021_06/nodes.dmp -d all.dmnd
```

Generate files for each assembly

```
for IDS in *s.ids; do
  DB=${IDS%*.ids}
  seqtk subseq all.fasta.gz $IDS | gzip -c > $DB.fasta.gz
  printf "accession\taccession.version\ttaxid\tgi\n" > $DB.taxid_map
  grep -Ff <(cut -d"|" -f 2 $IDS) all.taxid_map >> $DB.taxid_map
  diamond makedb -p 4 --in $DB.fasta.gz --taxonmap $DB.taxid_map --taxonnodes ../../taxdump_2021_06/nodes.dmp -d $DB.dmnd
done
```
