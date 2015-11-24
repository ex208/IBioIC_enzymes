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
conncatenated_genome.phrconncatenated_genome.pinconncatenated_genome.psq
```
The final file name output is:
```concatenated_blastp.out
```
However Blast includes other formats which they might 
I attempted to download all the sequences from PL1 family (in which PelE belongs to) by using CAZY. However, I realised that it is getting a bit harder as Dickeya genus has quite a few representations in this family. Further studying to understand how Cazy works. Studying to understand what kind of data I will need to collect. 
As long as it is difficult to collect raw data for Dickeya representation in each individual family from CAZy the goal is to find a way where i can easily collet and download all the protein sequences from either each individual family or each individual super-family or even the the whole Dickeya CAZyome. 
Developing a scrip where can automatically can lead you to download each family  would be an option
Using CAZYnes Analysis Toolkit (CAT)  might be another option [link](http://mothra.ornl.gov/cgi-bin/cat/cat.cgi)
dbCAN [link](http://csbl.bmb.uga.edu/dbCAN/)
[AHV.dk](http://www.ahv.dk/index.php/cazy-extract-cazyome)
Using CATlook I downloaded each indivualy family and made separate files for each of them containing the protein sequence in a fasta format.
One important observation is that the website is not updated according to CAZy. 
For example CAZy for GH5 family has 6 Dickeya representations while CAT gave me 4 as result!!!
Dickeya zeae EC1[link](https://open.library.ubc.ca/cIRcle/collections/ubclibraryandarchives/2689/items/1.0077841#downloadfiles)