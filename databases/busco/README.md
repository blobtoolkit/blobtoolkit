Aim: To create busco subsets with only a few busco genes per lineage for pipeline testing

Create a list of the lineages for which you want to make subsets:
```
echo 'carnivora_odb10
laurasiatheria_odb10
eutheria_odb10
mammalia_odb10
tetrapoda_odb10
vertebrata_odb10
metazoa_odb10
eukaryota_odb10' > lineages.txt
```

Now, for each lineage, copy a subset of files over, and a subset of genes for each lineage:

```
D=/path/to/downloaded/buscov5/lineages

cat lineages.txt | while read L
do
  N=`grep number_of_BUSCOs $D/$L/dataset.cfg | cut -f2 -d"="`
  C=`printf %.0f $(echo "$N *0.05" | bc -l)`

  # 0.05 = 5% of busco genes. Change if you want a higher percentage of genes

  rm -rf  $L
  mkdir -p $L && cd $L
  ls $D/$L/hmms | shuf | head -n $C | cut -f1 -d. > $L.subset
  mkdir -p hmms prfl
  cat $L.subset | while read a; do rsync -a $D/$L/hmms/$a.hmm  hmms/; done
  cat $L.subset | while read a; do rsync -a $D/$L/prfl/$a.prfl prfl/; done
  for a in ancestral ancestral_variants; do cat $D/$L/$a | paste - - | fgrep -f $L.subset - | sed 's/\t/\n/g' > $a; done
  for a in refseq_db.faa; do cat $D/$L/$a | paste - - | fgrep -f $L.subset - | sed 's/\t/\n/g' > $a; done
  for f in lengths_cutoff scores_cutoff; do fgrep -f $L.subset $D/$L/$f > $f; done
  cp $D/$L/dataset.cfg ./
  echo $L
  echo $N
  echo $C
  cd ..
done
```
