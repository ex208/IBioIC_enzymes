#Identification of RBBH 
######8-12 February 2016


* I created copies of the ***.fasta*** files from:


````
/Users/eirinixemantilotou/Documents/PhD/PhD_year1/Bioinformatics/Identifying_Dickeya_enzymes/IBioIC_Dickeya_data/features_aa

````

A new folder which is named **Identification_RBBH** was created and the .fasta copy files were stored within the folder: 

````
/Users/eirinixemantilotou/Documents/PhD/PhD_year1/Bioinformatics/Identifying_Dickeya_enzymes/Identification_RBBH

````

The folder includes:

````
1. CSL_RW192_renamed_genecall_translations.fasta
2. CSL_RW240_renamed_genecall_translations.fasta
3. DW_0440_renamed_genecall_translations.fasta
4. GBBC2039_renamed_genecall_translations.fasta
5. GBBC2040_renamed_genecall_translations.fasta
6. IPO_2222_renamed_genecall_translations.fasta
7. IPO_980_renamed_genecall_translations.fasta
8. MK10_renamed_genecall_translations.fasta
9. MK16_renamed_genecall_translations.fasta
10.MK19_renamed_genecall_translations.fasta
11.MK7_renamed_genecall_translations.fasta
12.NCPPB_2511_renamed_genecall_translations.fasta
13.NCPPB_2538_renamed_genecall_translations.fasta
14.NCPPB_2976_renamed_genecall_translations.fasta
15.NCPPB_3274_renamed_genecall_translations.fasta
16.NCPPB_3531_renamed_genecall_translations.fasta
17.NCPPB_3532_renamed_genecall_translations.fasta
18.NCPPB_3533_renamed_genecall_translations.fasta
19.NCPPB_3534_renamed_genecall_translations.fasta
20.NCPPB_3537_renamed_genecall_translations.fasta
21.NCPPB_402_renamed_genecall_translations.fasta
22.NCPPB_453_renamed_genecall_translations.fasta
23.NCPPB_516_renamed_genecall_translations.fasta
24.NCPPB_569_renamed_genecall_translations.fasta
25.NCPPB_898_renamed_genecall_translations.fasta
````
 * The next step was to browse through the command line within the directory I keep the translations:

 /Users/eirinixemantilotou/Documents/PhD/PhD_year1/Bioinformatics/Identifying_Dickeya_enzymes/Identification_RBBH
 
 and to open a jupiter notebook which I named Identification_RBBH.ipynb
 
* I viewed the content of the filed by using **!ls** to make sure that I am within the right directory and I have access to the files I will need. 

* The step was to create databases for each file. In order to achieve that from the python notebook I imported the subprocess module and I used the command 

"makeblastdb -in %s -dbtype prot -title %s -out %s" % (i,i,i)

After exerting the command, I checked within the directory I keep the data and I noticed that most of the databases had failed to be created. According to the command line there was no residues given. 

* In order to fix the problematic files I used another jupyter notebook which is named: filtering_FASTA_sequences_for_nonzero_length and is stored at 
 /Users/eirinixemantilotou/Documents/PhD/PhD_year1/Bioinformatics/Identifying_Dickeya_enzymes/Identification_RBBH. 
 
 Please read the notebook for further information.
 
 * The new/fixed files were stored at a new folder: /Users/eirinixemantilotou/Documents/PhD/PhD_year1/Bioinformatics/Identifying_Dickeya_enzymes/Identification_RBBH_2

* A new jupiter notebook was created under the /Users/eirinixemantilotou/Documents/PhD/PhD_year1/Bioinformatics/Identifying_Dickeya_enzymes/Identification_RBBH_2 directory and is named Identification_RBBH_2.ipynb

* I listed the content of the directory again to make sure I am working with the new files

* I run "makeblastdb -in %s -dbtype prot -title %s -out %s" % (i,i,i)
and 3 files for each genome was created within the same folder within the /Users/eirinixemantilotou/Documents/PhD/PhD_year1/Bioinformatics/Identifying_Dickeya_enzymes/Identification_RBBH_2  directory 

* I then run "blastp -query %s -db %s -out %s_%s.tab -outfmt 7" %(i,j,i.split('.')[0],j.split('.')[0]) in order to execute blast for each genome against every existing database of each other genome. 

* The procedure was quite time-consuming!

* After the command has run for all different combination I should end up with 24 .tab files for each genome which I can use to identify the RBBH. 







