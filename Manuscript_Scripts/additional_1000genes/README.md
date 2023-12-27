# Add 1000 genes due to the reviwers' comments:
The input files are the fasta and sub files, Justin gave to to me:
It is in the Dropbox/BAD_Mutations-LRT/Data/1000genes.tar.gz.

But in reality, only 979 genes because 21 of them don't have SNPs.

Here is Justin's email:
===
I've got files ready for BAD_Mutations
There are 979 fasta and peptide sequences and positions files
There are 40760 SNPs in 1001genomes_1000genes.update.csv

21 of the random 1000 genes selected have no snp calls. Probably dup or repeat genes where confidence is bad. Two are mitochondrial. So total is actually 979 genes.

Reminder: AT1G05260.1 had to be converted to AT1G05260_1 for hyphy which does not like "." in tree files. File names remain unchanged.
We should use original genome versions based on last years analysis, rather than download latest genomes and annotations. I have these if needed.

I will being working on PPH and SIFT output. You can start on BAD_Mutations

Justin

===

## Processes:

### Step1: rename the fasta and pos files:

```
rename .1.pos _1.pos *.1.pos
 
rename .1.fasta _1.fasta *.1.fasta

```
### Step2: create the file list for alignments:

The files are here:

[fasta.list00](https://github.com/lilei1/BAD_mutation_Ath/blob/master/additional_1000genes/input/fasta.list00)

[fasta.list01](https://github.com/lilei1/BAD_mutation_Ath/blob/master/additional_1000genes/input/fasta.list00)

Then run [run_alignment_list00.job](https://github.com/lilei1/BAD_mutation_Ath/blob/master/additional_1000genes/jobs/run_alignment_list00.job)
and run [run_alignment_list01.job](https://github.com/lilei1/BAD_mutation_Ath/blob/master/additional_1000genes/jobs/run_alignment_list01.job)

### step3: create file list for prediction:
sorted_MSA.list00

sorted_MSA.list01

sorted_fasta_tree.list00

sorted_fasta_tree.list01

sorted_sub.list01

sorted_sub_tree.list00

sorted_tree.list00

sorted_tree.list01

Then run [run_prediction_mesabi00.job](https://github.com/lilei1/BAD_mutation_Ath/blob/master/additional_1000genes/jobs/run_prediction_mesabi00.job)

run [run_prediction_mesabi01.job](https://github.com/lilei1/BAD_mutation_Ath/blob/master/additional_1000genes/jobs/run_prediction_mesabi01.job)

### step 4: create the ANNOVAR format file to run compile:

```
./run_creat_annova.sh /panfs/roc/groups/9/morrellp/llei/Deleterious_mutation_project/LRT_BAD_mutation/A_thaliana_BAD_Mutation/Additional1000genes/file_list/sorted_sub_tree.list00 >list00_annova
```

### step 5: run compile to get the effect table and logistic regression values:

```
module load bad_mutations/1.0
 
 ./BAD_Mutations.py -v DEBUG compile -P /panfs/roc/groups/9/morrellp/llei/Deleterious_mutation_project/LRT_BAD_mutation/A_thaliana_BAD_Mutation/Additional1000genes/out/list00 -S /panfs/roc/groups/9/morrellp/llei/Deleterious_mutation_project/LRT_BAD_mutation/A_thaliana_BAD_Mutation/Additional1000genes/annova_like_file/sorted_list00_annova
 
 ./BAD_Mutations.py -v DEBUG compile -P /panfs/roc/groups/9/morrellp/llei/Deleterious_mutation_project/LRT_BAD_mutation/A_thaliana_BAD_Mutation/Additional1000genes/out/list01 -S /panfs/roc/groups/9/morrellp/llei/Deleterious_mutation_project/LRT_BAD_mutation/A_thaliana_BAD_Mutation/Additional1000genes/annova_like_file/sorted_list01_annova

```
### step 6: combined all the files and sort:

```
cat /panfs/roc/groups/9/morrellp/llei/Deleterious_mutation_project/LRT_BAD_mutation/A_thaliana_BAD_Mutation/Additional1000genes/out/list01/Combined_Report.txt /panfs/roc/groups/9/morrellp/llei/Deleterious_mutation_project/LRT_BAD_mutation/A_thaliana_BAD_Mutation/Additional1000genes/out/list00/Combined_Report.txt >972_Combined_Report.txt

```
The final file is [here](https://github.com/lilei1/BAD_mutation_Ath/blob/master/additional_1000genes/output/972_Combined_Report.txt).


