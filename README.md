Project
=======


This repostitory contains everything for the BS32010 project including raw data, images created and scripts used.

Information on the scripts used and where to find them:

R scripts

affycode.R - This was used to carry out the differential expression analysis on the microarray data
RNAseq.R - This was used to carry out differential expression analysis on the RNA-Seq data
Script to get plots.R - Created the scatterplots used in the report

Phyloscript.R - This script was used to carry out a model test, create trees, and calculate regression. Names were                        modified for each gene.
model table.R - was used to write the results from the model test to a file.
projectstuff.R - This script plots the residuals together


Results

All of the residuals.png, trees.png, bootstrap.png, and modeltest.csv files contain the results of the phylogenetic analysis carried out on each gene.
The ROS1.png, ROS1Normalisd.png etc. show the results from script to get plots.R
densityRAW.pdf and densityRMA.pdf contain the histograms of before and after normalisation respectively.
sumRAW and sumRMA.pdf contain the boxplots of before and after normalisation respectively. 
annotResults.txt is the unsorted results from the RNASeq.R analysis
celResults.txt is the unsorted results file from the affycode.R analysis
top6diffexp.txt is the top 6 results from annotResults.txt when sorted by adjusted p value
top6affy.txt is the top 6 results from celResults.txt when sorted by adjusted p value
