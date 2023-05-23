# Sjogrens_Master_Internship_CMBI
This repository contains the full analysis I performed on multimodal datasets from patients with Sjögrens syndrome for my Masters internship (2023).
The internship was conducted at the department of Molecular BioSciences (formrly known as the CMBI) at the Radboudumc in Nijmegen, the Netherlands.
For any questions regarding the repository or its contents, you can mail me at sidney.vanderzande@ru.nl.

## Docker Image
For all analyses I used the Bioconductor v.3.16 docker image available at https://hub.docker.com/r/bioconductor/bioconductor_docker. Given that Singularity is installed on your cluster node, it can be installed using the command 
```
singularity pull docker://bioconductor/bioconductor_docker:RELEASE_3_16
```


## Scripts
The version number of each script is indicated in the file name. If, for any reason, you would like to request a previous version of any script mentioned here, and you do not have access to the Radboudumc Octarine cluster, please email me at the email address mentioned at the beginning.
Most of these packages use functions that I moved to the self-designed sjogren package. This package can be downloaded off of the Octarine cluster. Later on, I will make this package available to GitHub.

The main script I used in my analysis is the full_analysis.R script. This script contains the main single-cell analysis performed on all samples of SjS patients.

For B cell clonotyping, I used the clonotype_calling_final.R script.

For running the Immcantation analysis framework, I used the immcantation.R script.

For comparisons to healthy tonsil data from public datasets, I used the tonsil_b_analysis.R script.
