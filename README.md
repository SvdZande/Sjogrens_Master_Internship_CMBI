# A Multimodal Analysis of B Lymphocytes in Sjögren's Syndrome
This repository contains the full analysis I performed on multimodal datasets from patients with Sjögrens syndrome for my Masters internship (2023).
The internship was conducted at the department of Molecular BioSciences (formerly known as the CMBI) at the Radboudumc in Nijmegen, the Netherlands.
For any questions regarding the repository or its contents, you can mail me at sidney.vanderzande@ru.nl.

This internship was performed under the supervision of prof. dr. Martijn Huijnen and dr. Prashant Singh, in collaboration with the Department of Rheumatology at the Radboudumc.
A paper published previously by our group, which has some overlap with this project, can be found under the DOI 10.1136/annrheumdis-2021-221604.

## Docker Image
For all analyses, I used the Bioconductor v.3.17 docker image available at https://hub.docker.com/r/bioconductor/bioconductor_docker. Given that Singularity is installed on your cluster node, it can be installed using the command 
```
singularity pull docker://bioconductor/bioconductor_docker:RELEASE_3_17
```
If you have any questions on Singularity or Docker, or the CMBI in general, the wiki might be of help: https://gitlab.cmbi.umcn.nl/cmbi/general/-/wikis/home

##Data
Most of the data I used as input for my scripts is present in this GitHub repository. The cellranger filtered_contig_annotations are present for each donor. I also added the autoreactive CDR3 sequences found for each donor in the **Autoreactive_sequences** folder.

## Scripts
The version number of each script is indicated in the file name. If, for any reason, you would like to request a previous version of any script mentioned here, and you do not have access to the Radboudumc Octarine cluster, please email me at the email address mentioned at the beginning.
Most of these packages use functions that I moved to the self-designed __sjogren__ package. This package can be downloaded off of the Octarine cluster or off of this GitHub repository. I will not be further developing this package after completion of the project, but I am happy to answer any questions you may have on the functions that are in it.

The main script I used in my analysis is the _full_analysis.Rmd_ script. This script contains the main single-cell analysis performed on all samples of SjS patients.

For B cell clonotyping, I used the _clonotype_calling_final.Rmd_ script.

For running the Immcantation analysis framework, I used the _immcantation.Rmd_ script.

For comparisons to healthy tonsil data from public datasets, I used the _tonsil_b_analysis.Rmd_ script.

For cell-cell interaction analyses using cellchat, I used the _cellchat.Rmd_ script.

A script for the comparison of SjS data to spleen data is also included in this folder, however, this script was never used in any outputs, nor is it tested properly.

Hope this was informative!

Sidney
