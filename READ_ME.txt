Script #1:

- AIM: The overall aim is to modify the summary file which functions as a "road map" for the analysis. 
  * Removing lines where H chain and L chain has the same ID. 
  * Removing duplicate lines that belong to the same pdb. This is because some files contain multiple biological units. Continue with the line where complex has lowest b factor. 
  * Extracting fasta sequences of antibody variable domains.
  * Counting the number of amino acids in the different antibody domains of each chain. 

- Takes "Summary_all.tsv", "Thera_all.tsv" and pdb files located in "imgt_all_clean" folder. 
- The pdb folder path and the input files are hard coded after defining the functions (approx. line 171). 
