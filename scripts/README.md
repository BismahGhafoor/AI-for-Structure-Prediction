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

## dbd_identification.py
- This script identifies the DNA bidning domains (DBDs) in each FASTA file for the TFs via local search of downloaded alignment files and/or UniProt search based on specified keywords. It outputs the locations of the DBDs and a .txt file taht contains names of any TFs that had no DBD identified. 
- This script uses python (PYTHON VERSION)
- Run the **downloadZipfile.py** script to download the neccesary file.
- You can then successfully run this script to identify the DBD.
