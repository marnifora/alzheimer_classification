### Scripts documentation
<br></br>
#### Preprocessing vcf files downloaded from the database
<br></br>
#### prepreparing.sh

Preprocessing of downloaded vcf files. Unpacking, filtering SNPs, running vcf_to_matrix.py for every given chromosome 
(see below).

Multiprocessing for Bash used in this file is available thanks to script job_pool.sh created by Vince Tse with changes 
by Geoff Clements (https://github.com/vincetse/shellutils/blob/master/job_pool.sh).

Selection of SNPs variants from vcf files was carried out using GATK tool (https://software.broadinstitute.org/gatk/).

##### Input:
- {istart}{chr}{iend}.vcf.gz.tar files containing WGS information for each chromosome
- {reference}.fasta file with reference genome for given vcf files - default named "human_g1k_v37.fasta"

##### Output:
- {base}_chr{chr_number}_SNPs.vcf - vcf file with information about SNPs from chr {chr_number} for every patient, 
{base} is a name of database where the data comes from
- genome_stats.txt - file with information about number of patients and SNPs on each chromosome
If -matrix or -stats options are on there is more output from different scripts (see cmd parameters section).

##### Cmd parameters:
- -chr [NR] - NR is number of chromosome to analyse
- -all - all chromosomes (for human: from 1 to 23) given to analyse
- -x - last chromosome is named 'X' despite '23'
- -from [NR] - the beginning of the scope of chromosomes to analyse
- -to [NR] - the end of the scope of chromosomes to analyse
- -tar - unpacking files from .tar
- -gz - unpacking files from .gz
- -snp - vcf files filtering for SNPs after unpacking
- -stats - running vcf_stats.py for every vcf file (additional statistical analysis)
- -matrix - running vcf_to_matrix.py for every vcf file (rewriting vcf files into numpy matrices and txt files)
- -base [DATABASE] - set name of database where data comes from
- -indir [DIR] – input directory
- -outdir [DIR] – output directory
- -gatkdir [DIR] - directory to folder with GATK tool
- -reference [DIR+NAME] - name of the file with reference genome in fasta format
- -name [NAME] - pattern of name of vcf file with genetic data, if the name is different for every chromosome it should 
be marked as in the example: 'name{chr}rest' where {chr} is the number of the chromosome

<br></br>
#### vcf_stats.py

Making basic statistics of given vcf file. Writing number of SNPs, number of patients and their IDs into file.

##### Input:
- {input}.vcf file to analyse

##### Output:
- {input}_stats.txt - file which contains statistics for vcf file, {input} is the same string as in the input file

##### Cmd parameters:
- -input [NAME] - set path+name of the input file


<br></br>
#### vcf_to_matrix.py

Processing given vcf file. Writing order of patients IDs into txt file, writing order of SNPs and their values into txt 
file, writing genetic information into numpy matrix.

##### Input:
- {database}_chr{chr_number}_SNPsMNPs.vcf - vcf file from database named {database} for given chromosome contains only 
SNPs and MNPs

##### Output:
- pid_chr{chr_number}.txt - order of patients ID numbers in data for given chromosome (each line is one ID number), 
example:
```
    1280
    
    1901
    
    2300
```    
- snps_chr{chr_number}.txt - information about SNPs on given chromosome: each line is one SNP, first column contains 
position of given SNP on the chromosome, second column is the reference nucleobase in this position, third column 
contains alternative nucleobases occurred in the given position, example:
```
    1234    A    C
    
    2344    C    G,A,T
    
    3071    C    T,A
```    
- matrix_chr{chr_number}.npy - information about SNPs (each column is one SNP) from given chromosome for all patient 
(each row is one patient), example:

```
  [ [ [0,1]  ,  [-1,-1]  ,  [1,2] ], 
  
    [ [1,1]  ,   [2,3]   ,  [0,1] ], 
    
    [ [0,0]  ,  [-1,-1]  ,  [2,2] ]  ]
``` 

To sum up all given examples for given chromosome:

Element (0,0) in the 3. table which is '0,1' means that SNP on the position 1234 ((0,0) element of the snps_chr table)
 for patient whose ID=1280 (0 line of pid_chr file) is A in the first allele (beacuse first number is 0) and C in the 
 second allele (beacuse second number is 1).
 
Element (0,1) in the 3. table, which is '-1,-1', means that for SNP on the position 2344 ((1,0) element of the snps_
chr table) for the same patient (because it is the same row) there is no data. 

Element (1,1) in the 3. table, which is '2,3' means that for SNP on the position 2344 ((1,0) element of the snps_chr
table)) for patient 1901 (1 line of the pid_chr file) is A in the first allele (because first number is 2) and T in
the second allele (beacuse second number is 3).

##### Cmd parameters:
- -chr [NR] - number of chromosome which vcf file is given to analyse
- -input [NAME] - set path+name of the input file
- -output [DIR] – output directory


<br></br>
#### Making X and Y matrices needed to classification

<br></br>
#### make_pid-diagnoses.py

Unification all pid_chr{chr_number}.txt files in one pid_chr.txt file (after checking whether the order is the same in 
every file). Making file diagnoses.txt with diagnose for every patient in the same order as in the pid_chr.txt file.

For every database there should be different function of mapping diagnoses! Because different data sets have different 
coding for the diagnoses.

##### Input:
- pid_chr{chr_number}.txt - file with order of patients IDs in data for {chr_number} chromosome, output file from 
vcf_to_matrix.py
- {diagnose_file} - file or files with diagnoses assigned to every patient

##### Output:
- pid_chr.txt - list of patients IDs for every chromosome, order of patients is universal, the same for every file from 
rosmap dataset
- diagnoses.txt - list of diagnoses (AD/NL/DIF/NN) for every patient, order of diagnoses is like in the pid_chr.txt

##### Cmd parameters:
- -indir [DIR] - input directory
- -outdir [DIR] - output directory
- -dataset [NAME] - name of data set for which diagnoses should be established (e.g. adni or rosmap)
- -diagdir [DIR] - directory to a folder with file(s) based on which diagnoses can be established


<br></br>
#### makeY.py

Making csv file with Y matrix needed for classification process (it contains information about classes to which the 
objects are assigned). Writing numbers of patients with diagnoses different than NL/AD into txt file. Updating 
genome_stats.txt file - add to every line number of patients with diagnosis NL or AD (who can be used in further 
analysis).

##### Input:
- pid_chr.txt - list of patients' IDs, output of make_pit-diagnoses.py
- diagnoses.txt - diagnoses for all patients in the order corresponding to the order from pid_chr file. Possible 
diagnoses:
	AD - Alzheimer's disease
	DIF - different diagnosis (for example Mild Cognitive Impairment)
	NL - normal
	NN - no data
- genome_stats.txt - file with information about number of patients and SNPs on each chromosome

##### Output:
- Y_chr.csv - list of diagnoses coded as 0 (NL) or 1 (AD) in the order corresponding to the order in pid_chr file. 
Patients with DIF diagnosis was excluded from this file.
- dif_chr.txt - list of IDs of patients having other diagnosis than NL/AD and number of line in which they are written 
in pid_chr file
- genome_stats.txt (UPDATED)

##### Cmd parameters:
- -indir [DIR] - input directory
- -outdir [DIR] - output directory


<br></br>
#### makeX_pooling.sh

Running function makeX.py for every given chromosome.

Multiprocessing for Bash used in this file is available thanks to script job_pool.sh created by Vince Tse with changes 
by Geoff Clements (https://github.com/vincetse/shellutils/blob/master/job_pool.sh).

##### Input: 
see input for makeX.py script (below)

##### Output:
see output for makeX.py script (below)

##### Cmd parameters: 
- -chr [NR] - NR is number of chromosome to analyze
- -all - all chromosomes (for human: from 1 to 23) given to analyze
- -from [NR] - the beginning of the scope of chromosomes to analyze
- -to [NR] - the end of the scope of chromosomes to analyze
- -indir [DIR] - input directory
- -outdir [DIR] - output directory


<br></br>
#### makeX.py

Selection of patients with appropriate diagnosis (rejection of diagnoses different than AD/NL – see makeY.py above). 
Making X matrix (based on selected patients) for each chromosome. Writing them to csv file. Matrix X is needed for 
classification process, it contains information about every SNP (columns of matrix) for every patient (rows of matrix).

##### Input:
- dif_chr.txt - list of patients with DIF diagnosis, output from makeY.py.
- matrix_chr{chr_number}.npy - numpy matrix, output from vcf_to_matrix.py.

##### Output:
- X_chr{chr_number}.csv - csv table, where columns are SNPs, rows are patients. It contains information from 
matrix_chr{chr_number}.npy but only for first allele (one number in each position of the table).
- X_chr{chr_number}_nodif.csv - similar to X_chr{chr_number}.csv file but without patients with DIF diagnosis.

##### Cmd parameters:
- -chr [NR] - NR is number of chromosome to analyze
- -indir [DIR] - input directory
- -outdir [DIR] - output directory


<br></br><br></br>
#### Selection of the most important SNPs (Boruta algorithm), building classifier, carrying out the classification

<br></br>
#### boruta_classification.py

Running Boruta algorithm on the given data set(s), building classifier based on the chosen by Boruta SNPs from given 
train set, conducting classification based on the given test set.

##### Input:
- genome_stats.txt - file with information about number of patients and SNPs on each chromosome
- X_chr{chr_number}_nodif.csv - csv table, where columns are SNPs, rows are patients, output of makeX.py
- Y_chr.csv - list of diagnoses for each patient, output of makeY.py
###### Optional input:
- {subset}_snps_chr{chr_number}.txt - list of SNPs from chromosome {chr_number}, belonging to {subset} (see section 
“Subsets of SNPs” below)

##### Output:
- bestsnps_chr{chr_number}\_{perc}\_{run_number}.txt - list of chosen SNPs from chromosome {chr} by Boruta with parameter 
perc={perc} in the run {run_number}
- testpat_{run_number}.txt – drawn test set of patients from the run {run_number}
- class_scores_{run_number}.txt – result of classification conducted in the run {run_number}
- X_train_genome_{perc}_{run_number}.npy
- y_train_genome_{perc}_{run_number}.npy
- X_test_genome_{perc}_{run_number}.npy
- y_test_genome_{perc}_{run_number}.npy
- X_train_chr{chr_number}\_{perc}\_{run_number}.npy
- y_train_{run_number}.npy
- X_test_chr{chr_number}\_{perc}\_{run_number}.npy
- y_test_{run_number}.npy

Make in first run (another runs write in addition):
- boruta_runs.txt -  list of conducted runs of Boruta analysis with their parameters
- class_runs.txt – list of conducted runs of classification with their parameters




##### Cmd parameters:
- -dataset [NAME] [DIR] - set of data which should be used to Boruta analysis
- -testset [NAME] [DIR] - set of data which should be used as test data for classification
- -test [SIZE] - size of subset of elements which should be used as test set for classification, default = 0
- -outdir [DIR] - directory for output files, if not given directory of first set in dataset is taken as outdir
- -perc [VALUE] - value(s) of perc parameter of Boruta analysis, if more than one given it should be written as a list 
(e.g. '-perc [80,90,100]'), default = 90
- -classperc [VALUE] - value of perc parameter of Boruta analysis on which result the classification should be based, if 
not given first value from perc list is set as class_perc
- -r [SIZE] - size of the window of Boruta (how many of SNPs Boruta analyze in one round), default = 5000
- -chr [RANGE] - set of chromosomes to analyze from given data set, e.g. '-chr 1-12,14,18-20', default = ‘1-23’
- -class - if only classification should be run
- -boruta - if only Boruta analysis should be run
- -run [VALUE] - number of run of both Boruta and classification analyses, default is the next number from run files
- -borutarun [VALUE] - number of run of Boruta analysis
- -classrun [VALUE] - number of run of classification analysis
- -fixed - if number of run or list of chromosomes shouldn't be changed even is it was run before
- -subset [NAME] - name of SNPs subset which should be under consideration during analysis, e.g. '-subset shared'
- -subsetrun [VALUE] - number of run of subset establishing
- -cont - if continuation of started Boruta analysis should be run (for different chromosomes)
- -makeY - if y matrices should be established before classification

##### Parameters required for different type of analyses:
1. Boruta + classification
    - -dataset [NAME] [DIR]
    - -testset [NAME] [DIR] or -test [SIZE]
2. Only Boruta
    - -boruta
    - -dataset [NAME] [DIR] (only one data set should be given)
3. Only classification
    - -class
    - -borutarun [VALUE]
    - -testset [NAME] [DIR] or -test [SIZE]
4. Continuation of previous started Boruta analysis:
    - -cont
    - -boruta
    - -chr [RANGE]
    - -borutarun [VALUE]



<br></br><br></br>
#### Additional analyzes

<br></br>
#### SNPs subsets

##### shared_snps.py

Subset 'shared' from the given group of data sets contains SNPs which occur in every of them.

1. shared_snps_chr{chr_number}.txt - list of SNPs from chromosome {chr_number} which get through the pruning AND occur 
also in data from other cohort. Output of compare_snps_pruned.py

<br></br><br></br>