# README.md - Sequences_aa_CAZy

## Overview

Sequences for Cel5Z obtained from CAZy on 6 November 2015 (CAZy release 2015-11-03)


***GH5 family ****

Concatenated sequences:

```
$ ls -1 *.fasta
Dickeya_chrysanthemi_PY35.fasta
Dickeya_dadantii_3937.fasta
Dickeya_dadantii_Ech586.fasta
Dickeya_dadantii_Ech703.fasta
Dickeya_zeae_EC1.fasta
Dickeya_zeae_Ech1591.fasta
concatenated.fasta
$ cat *.fasta > bogus_cat.fasta
```


The above sequences were downloaded in fasta format one-by-one from caZy.

#### Concatenated genome sequences from Dickeya:


$ cd /Users/eirinixemantilotou/Documents/PhD/PhD_year1\ /Bioinformatics/Identifying_Dickeya_enzymes/IBioIC_Dickeya_data/features_aa 


```
SL_RW192_renamed_genecall_translations.fasta
CSL_RW240_renamed_genecall_translations.fasta
DW_0440_renamed_genecall_translations.fasta
GBBC2039_renamed_genecall_translations.fasta
GBBC2040_renamed_genecall_translations.fasta
IPO_2222_renamed_genecall_translations.fasta
IPO_980_renamed_genecall_translations.fasta
MK10_renamed_genecall_translations.fasta
MK16_renamed_genecall_translations.fasta
MK19_renamed_genecall_translations.fasta
MK7_renamed_genecall_translations.fasta
NCPPB_2511_renamed_genecall_translations.fasta
NCPPB_2538_renamed_genecall_translations.fasta
NCPPB_2976_renamed_genecall_translations.fasta
NCPPB_3274_renamed_genecall_translations.fasta
NCPPB_3531_renamed_genecall_translations.fasta
NCPPB_3532_renamed_genecall_translations.fasta
NCPPB_3533_renamed_genecall_translations.fasta
NCPPB_3534_renamed_genecall_translations.fasta
NCPPB_3537_renamed_genecall_translations.fasta
NCPPB_402_renamed_genecall_translations.fasta
NCPPB_453_renamed_genecall_translations.fasta
NCPPB_516_renamed_genecall_translations.fasta
NCPPB_569_renamed_genecall_translations.fasta
NCPPB_898_renamed_genecall_translations.fasta

```

**$ cat *.fasta > concatenated.fasta**

By using the above command I merged all the Dickeya genomes in one file.

The next step was to create a database of the 24 Dickeya genomes. In order to achieve that I used the following scrip:

```
$ makeblastdb -in concatenated.fasta -dbtype prot -title " concatenated database" -out concatenated_genome

```

In order to compare the 6 genes against the 24 gnomes I wrote a script for blastp (compare amino acid sequence)
The script is the above:

```

$ blastp -query concatenated.fasta -db conncatenated_genome -out concatenated_blastp.out

```

After creating a datable for the 24 genomes 3 files were created :

```
conncatenated_genome.phr














