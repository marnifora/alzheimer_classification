# Alzheimer classification


<br></br>
### The aim of this project is to classify patients into groups healthy/ill based on their genetic information.
<br></br>

#### Versions

GNU bash 4.4.19(1)

GATK 4.0.10.1 (https://software.broadinstitute.org/gatk/)

Python 3.6.6
- boruta 0.1.5
- numpy 1.15.4
- pandas 0.23.4
- scikit-learn 0.20.1

<br></br>
#### How to use

Processing vcf files:

- prepreparing.sh
- make_pid-diagnoses.py
- makeY.py
- makeX_pooling.sh

Selection of attributes and classification:

- boruta_classification.py
<br></br>
#### Steps of an examplary analysis
Using prepared testing data (400 patients, ~38k SNPs) from "./testing/files" directory.
- prepare csv matrices from vcf files
```
./prepreparing.sh -all -tar -gz -stats -matrix -base test -vcf test_chr_SNPs.vcf -dir ${PWD}/testing/
```
- prepare file with diagnoses
```
python make_pid-diagnoses.py -dir ${PWD}/testing/ -diagdir ${PWD}/testing/diagnoses/ -dataset test
```
- build Y vectors (containing diagnosis for each patient)
```
python makeY.py -dir ${PWD}/testing/
```
- build X matrices
```
./makeX_pooling.sh -all -dir ${PWD}/testing/
```
- run boruta in the correct way (train/test split before selection of important SNPs)
```
python boruta_classification.py -boruta -dataset test ${PWD}/testing/ -borutarun 1 -test 0.1
```
- run boruta in the wrong way (no train/test slit before selection of important SNPs)
```
python boruta_classification.py -boruta -dataset test ${PWD}/testing/ -borutarun 2
```
- build random classifier based on SNPs selected from the first run of boruta, proceed classification
```
python boruta_classification.py -class -dataset test ${PWD}/testing/ -borutarun 1 -classrun 1
```
- build random classifier based on SNPs selected from the second run of boruta, proceed classification
```
python boruta_classification.py -class -dataset test ${PWD}/testing/ -borutarun 2 -classrun 2 -test 0.1
```
Results of the correct classification procedure are saved in "./testing/dataset/boruta/class_results_1.txt".
Results of the wrong classification procedure are saved in "./testing/dataset/boruta/class_results_2.txt".

<br></br>
#### Description

Project initially started as the part of my Bachelor's degree thesis, named "Classification of patients with Alzheimer's 
disease based on DNA polymorphisms". Later it has developed into bigger project of patients classification based on WGS 
and GWAS data from ADNI and Rosmap consortia.

##### Used data
Three sets of data have been used:
1. Whole genome sequencing (WGS) data of 486 patients (235 cases, 251 controls), obtained from ADNI consortium 
(https://adni.loni.usc.edu/).
2. WGS data of 1033 patients (530 cases, 503 controls), obtained from Rosmap project 
(https://www.synapse.org/#!Synapse:syn10901595).
3. Data from Genome-wide association study (GWAS) of 432 patients, obtained from ADNI consortium.


##### Basic steps of analysis
The basic analysis of each data set can be described by following steps:

1. Rewriting genetic data into matrices and information about patients and diagnoses into text files.
2. Division of patients into training and testing set.
3. Selection of the most important SNPs by **Boruta algorithm** based on training set of patients.
4. Training the **Random Forest** classifier based on selected SNPs.
5. Testing the classifier based on testing set of patients.

##### Boruta algorithm

**Boruta algorithm** has been developed by Miron B. Kursa and Witold R. Rudnicki ("Feature Selection with the Boruta Package" - 
Journal of Statistical Software, Vol 36 (2010), https://www.jstatsoft.org/article/view/v036i11). In this project Boruta
implementation for Python made by Daniel Homola 
(http://danielhomola.com/2015/05/08/borutapy-an-all-relevant-feature-selection-method/) was used.

Boruta is a feature selection method. It is designed as a wrapper around a Random Forest classification algorithm. 
It iteratively removes the features which are proved by a statistical test to be less relevant than random probes. 

##### Additional analysis

- analysis of subset of SNPs (e.g. SNPs shared between two data sets)
- removing outlier patients


<br></br>
