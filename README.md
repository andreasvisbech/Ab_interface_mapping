

# DSSP
DSSP is required for script 3. The latest version is dssp 4, but biopython, which we use to run dssp with, doesnt support version 4 yet (08/07/2022). So for now, you must use dssp 3 (dssp 2 will also work). You can edit script 3 to point to a specific dssp, but by default it will just try to run "dssp".  
Dssp can either be compiled from source or installed via your linux package manager, if available. (dssp 4 is available in ubuntu 21.10+22.04, dssp3 is available in 18.04+20.04)
- [dssp 3 github](https://github.com/cmbi/dssp)
- [dssp 4 github](https://github.com/PDB-REDO/dssp)

# Script 1

- AIM: The overall aim is to modify the summary file which functions as a "road map" for the analysis.
	- Removing lines where H chain and L chain has the same ID.
	- Removing duplicate lines that belong to the same pdb. This is because some files contain multiple biological units. Continue with the line where complex has lowest b factor.
	- Extracting fasta sequences of antibody variable domains.
	- Counting the number of amino acids in the different antibody domains of each chain.

- Takes "Summary_all.tsv" and pdb files located in "imgt_all_clean" folder.

# Script 2

- AIM: The overall aim is to modify the summary file
	- Exclude lines that represent homologous antibodies based on the clustering analysis.
	- The script should only exclude antibodies where both chains are placed in the same clusters. E.g. if Ab1 and Ab2 both have H chain placed in cluster 1 and L chain placed in cluster 2 then only one of the structures should stay. If Ab1 has H chain in cluster 1 and L chain in cluster 2 and Ab2 has H chain in cluster 3 and L chain in cluster 2 then both structures should stay.

- Takes 'Summary_all_sorted.csv' hard coded on line 3
- Takes "Cluster.txt" hard coded on line 36. The cluster file is made with CD-HIT and the fasta file generated from script #1.

# Script 3

- AIM:
	- Takes summary file from script #2 and pdb files. 
    - Then uses summary file as road map and looks through all pdb files for atom-atom contacts.

- path for pdb files is hard coded on line 314
- file for dssp analysis specified on line 395
- summary file is hard coded on line 680

# Script 4
- AIM:
	- Takes the file from script #3 containing contact points and calculates summary statistics. 

  
