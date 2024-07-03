# Data Collection Scripts
## downloadTFs.sh
- This script downloads all the human transcription factors (identified by Lambert _et al._, (2018)) FASTA sequences.
- This script uses Bash Shell scripting 
- For the script to work ensure you have downloaded the genenames of the transcription factors from this URL: https://humantfs.ccbr.utoronto.ca/download.php and save it in the same directory as the script.
- The script can now be run.

## downloadEP300.sh
- This script downloads the EP300 protein's FASTA sequence.
- This script uses Bash Shell scripting
- There are no requirements for this script to work, it can be run straight away. 

# Data preprocessing Scripts
## DBD_identification.py
- This script identifies the DNA binding domains (DBDs) in each FASTA file for the TFs via local search of downloaded alignment files and/or UniProt search based on specified keywords. It outputs the locations of the DBDs and a .txt file taht contains names of any TFs that had no DBD identified. 
- This script uses python
- Run the **downloadZipfile.py** script to download the neccesary file.
- You can then successfully run this script to identify the DBDs.

## AD_identification.py
- This script identifies the activation domains (ADs) in eaxh FASTA file for the TFs using data that has been downloaded from the Alerasool _et al._, (2022) study.
- This script uses Python
- Run the **downloadADdata.py** script to download the necessary data in the same same directory as the AD_identification.py script.
- You can then successfully run this script to identify the ADs. 

# AlphaFold Structrual Predictions
## Customfile.py
- This script creates the .txt file needed to run AlphaFold2 predictions on all of the TFs on the command line interfact (CLI).
- This script uses python
- Run this script given you have already ran the **DBD_identification.py** script.
- After this script run the **filterCustomfile.py** script to filter any predictions that may be insignificant.

## AD_Customfile.py
- This script creates the .txt file needed to run AlphaFold2 predictions on the command line interfact (CLI).
- This script uses python
- Run this script given you have already ran the **AD_identification.py** script.
- This script predicts the AD regions of the TFs with annotated ADs against only the Taz2 domain of the p300 protein
  
## submit_alphapulldown_alignments.sh
- This scripts automates the submission of each TF for AlphaPulldown to produce a MSA to the HPC.
- This script uses Bash Shell scripting
- You must have all of the TF FASTA sequences in one directory in order for this script to run successfully.

## submit_alphapulldown_AD_alignments.sh
- This scripts automates the submission of each TF (with an annotated AD) for AlphaPulldown to produce a MSA to the HPC.
- This script uses Bash Shell scripting
- You must have all of the TFs with annotated ADs FASTA sequences in one directory in order for this script to run successfully.

## submit_alphafold_predictions.sh
- This scripts automates the submission of the TF AlphaFold predictions to the HPC.
- This script utilises Bash Shell scripting
- You need to run the **submit_alphapulldown_alignments.sh** script prior to running this script.

## submit_alphafold_AD_predictions.sh
- This scripts automates the submission of the TFs (with annotated ADs) AlphaFold predictions againt p300 core domain to the HPC.
- This script utilises Bash Shell scripting
- You need to run the **submit_alphapulldown_AD_alignments.sh** script prior to running this script. 
