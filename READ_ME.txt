SCRIPT #1:

- AIM: The overall aim is to modify the summary file which functions as a "road map" for the analysis. 
  * Removing lines where H chain and L chain has the same ID. 
  * Removing duplicate lines that belong to the same pdb. This is because some files contain multiple biological units. Continue with the line where complex has lowest b factor. 
  * Extracting fasta sequences of antibody variable domains.
  * Counting the number of amino acids in the different antibody domains of each chain. 

- Takes "Summary_all.tsv", "Thera_all.tsv" and pdb files located in "imgt_all_clean" folder. 
- The pdb folder path and the input files are hard coded after defining the functions (approx. line 171). 

SCRIPT #2:

- AIM: The overall aim is to modify the summary file
 * Exclude lines that represent homologous antibodies based on the clustering analysis. 
 * The script should only exclude antibodies where both chains are placed in the same clusters. E.g. if Ab1 and Ab2 both have H chain placed in cluster 1 and L chain placed in cluster 2 then only one of the structures should stay. If Ab1 has H chain in cluster 1 and L chain in cluster 2 and Ab2 has H chain in cluster 3 and L chain in cluster 2 then both structures should stay. 
 
- Takes 'Summary_all_sorted.csv' hard coded on line 3
- Takes "Cluster.txt" hard coded on line 36. The cluster file is made with CD-HIT and the fasta file generated from script #1. 
