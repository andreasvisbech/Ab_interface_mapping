import pandas as pd
from Bio.PDB import *
import statistics

# Loading functions
def get_Ab_type(VH_id, VL_id):
    
    func_dict = {}
    
    if VH_id.find('nan') < 0 and VL_id.find('nan') < 0:
        func_var = 'Fv'
        #func_dict['Ab_ids'] = [VH_id , VL_id]
        
    elif VH_id.find('nan') < 0 and VL_id.find('nan') >= 0:
        func_var = 'VH sdAb'
        #func_dict['Ab_ids'] = [VH_id]
        
    elif VH_id.find('nan') >= 0 and VL_id.find('nan') < 0:
        func_var = 'VL sdAb'
        #func_dict['Ab_ids'] = [VL_id]
        
    return func_var

def get_VHVL_seq(structure, V_id):
    
    seq = ''
    
    # The function extracts the amino acid sequence of the VH or VL domain.
    for Ab_resi_func in structure[0][V_id]:
        
        #Check that the residue is an amino acid
        if any(str(Ab_resi_func.get_resname()) == x for x in aa_list):
            
            # Check if the aa belongs to V domain if the residue ID is smaller or equal to 128
            if int(Ab_resi_func.get_id()[1]) <= 128:
                
                Ab_aa_3letter = str(Ab_resi_func.get_resname())
                Ab_aa = aa_dict[Ab_aa_3letter]
                seq = seq + Ab_aa
                
    return seq
            
# Setting path for the pdb files
path = './imgt_all_clean/'

# Initialize variables for run
aa_list = aa_list = ['ARG','HIS','LYS','ASP','GLU','SER','THR','ASN','GLN','CYS','GLY','PRO','ALA','ILE','LEU','MET','PHE','TRP','TYR','VAL']
aa_dict = {'ARG':'R' , 'HIS':'H', 'LYS':'K', 'ASP':'D', 'GLU':'E', 'SER':'S', 'THR':'T', 'ASN':'N', 'GLN':'Q', 'CYS':'C', 'GLY':'G', 'PRO':'P', 'ALA':'A', 'VAL':'V', 'ILE':'I', 'LEU':'L', 'MET':'M', 'PHE':'F', 'TYR':'Y', 'TRP':'W'}

# Specify the summary file
data = pd.read_csv('Summary_all.tsv', sep='\t',header=0)
thera_data = pd.read_csv('Thera_all.tsv', sep='\t',header=0)
thera_pdb_list = thera_data['pdb'].tolist()

# Load parser
parser = PDBParser()

# Load empty string for storing fasta in
fasta = ''

### Writing all pdb IDs to a list
pdb_list = data['pdb'].tolist()

### Small code block to remove all antibodies where antigen is N/A or H chain and L chain are the same i.e. scFvs where H chain and L chain are not properly annotated
drop_list = []
for a in range(len(pdb_list)):

    if str(data.iloc[a][4]) == 'nan' or str(data.iloc[a][1]) == str(data.iloc[a][2]):
        drop_list.append(a)

    else:
        None          
data = data.drop(labels=drop_list, axis=0)

# Create for loop with the purpose of removing duplicate biological units in the crystal structure.
# The loop should go over the individual complexes and only continue with the complex that has lowest average b factor
pdb_list_unique = data['pdb'].drop_duplicates().tolist()

pd_list = []

for b in range(len(pdb_list_unique)):
    
    print(b)
    
    # Loading in the pdb file
    #pdb2 = str(path + str(data.iloc[a][0]) + '.pdb')
    #structure = parser.get_structure("structure", pdb2)
    
    ### Loading the pdb structure into python
    pdb1 = str(path + pdb_list_unique[b] + '.pdb')
    structure = parser.get_structure("structure", pdb1)
    
    # Getting data from dataframe relating to the specific pdb
    my_pd = data[data['pdb'] == pdb_list_unique[b]]
    
    # Check how many rows in the remaining data rows contain the pdb id.
    # Since each row in summary files represents a single biological unit multiple rows indicate the crystal
    # contains more than one biological unit. If only one row is present there is no redundant data in the crystal.
    if my_pd.shape[0] == 1:
        pd_list.append(my_pd.iloc[0].tolist())
        
    elif my_pd.shape[0] > 1:
        
        bfactor_avr = []
        
        for c in range(my_pd.shape[0]):
            
            Ag_chain_list = []
            Ab_chain_list = []
            
            VH_id = str(my_pd.iloc[c][1])
            VL_id = str(my_pd.iloc[c][2])
            
            if len(str(my_pd.iloc[c][4])) == 1:                  ###Append the antigen chains to the chain list
                Ag_chain_list.append(str(my_pd.iloc[c][4][0]))
            elif len(str(my_pd.iloc[c][4])) == 5:
                Ag_chain_list.append(str(my_pd.iloc[c][4][0]))
                Ag_chain_list.append(str(my_pd.iloc[c][4][4]))
            elif len(str(my_pd.iloc[c][4])) == 9:
                Ag_chain_list.append(str(my_pd.iloc[c][4][0]))
                Ag_chain_list.append(str(my_pd.iloc[c][4][4]))
                Ag_chain_list.append(str(my_pd.iloc[c][4][8]))
            
            if VH_id.find('nan') < 0 and VL_id.find('nan') < 0:         ### Append the Ab chains to the chain list
                Ab_type = 'Fv'
                Ab_chain_list.append(str(my_pd.iloc[c][1]))
                Ab_chain_list.append(str(my_pd.iloc[c][2]))
                
            elif VH_id.find('nan') < 0 and VL_id.find('nan') >= 0:
                Ab_type = 'VH sdAb'
                Ab_chain_list.append(str(my_pd.iloc[c][1]))
                
            elif VH_id.find('nan') >= 0 and VL_id.find('nan') < 0:
                Ab_type = 'VL sdAb'
                Ab_chain_list.append(str(my_pd.iloc[c][2]))
            
            bfactor_list = []
            
            # Go over all atoms in antibody and get bfactors for each atom
            for chain_ID in Ab_chain_list:
                for residue_A in structure[0][chain_ID]:
                    if any(str(residue_A.get_resname()) == x for x in aa_list): # Only include if the residue is an amino acid
                        
                        for atom_A in residue_A:
                            bfactor_list.append(atom_A.get_bfactor())

            # Go over all atoms in antigen and get bfactors for each atom
            for chain_ID in Ag_chain_list:
                for residue_B in structure[0][chain_ID]:
                    if any(str(residue_B.get_resname()) == x for x in aa_list):         ### Only include if the residue is an amino acid
                        for atom_B in residue_B:
                            bfactor_list.append(atom_B.get_bfactor())
    
            bfactor_avr.append(statistics.mean(bfactor_list))

        bfactor_min_idx = bfactor_avr.index(min(bfactor_avr))
        pd_list.append(my_pd.iloc[bfactor_min_idx].tolist())

data = pd.DataFrame(data=pd_list, columns=data.columns.values)          

# Creating for loop for extracting fasta sequence of the antibody V domains. Only amino acids with ID below 128 is included. 
for a in range(data.shape[0]):
    
    print(str(a) + ' --- ' + str(data.iloc[a][0]))
    
    VH_id = str(data.iloc[a][1])
    VL_id = str(data.iloc[a][2])
    
    # Loading in the pdb file
    pdb2 = str(path + str(data.iloc[a][0]) + '.pdb')
    structure = parser.get_structure("structure", pdb2)
    
    # Getting the antibody type and the associated chain ids for antibody chains
    VH_id = str(data.iloc[a][1])
    VL_id = str(data.iloc[a][2])
    Ab_type = get_Ab_type(VH_id, VL_id)
    
    if Ab_type == 'Fv':
        VH_seq = get_VHVL_seq(structure, VH_id)
        VL_seq = get_VHVL_seq(structure, VL_id)
        
        fasta = fasta + '>' + str(data.iloc[a][0]) + ', VH' + '\n' + VH_seq + '\n'
        fasta = fasta + '>' + str(data.iloc[a][0]) + ', VL' + '\n' + VL_seq + '\n'
        
    elif Ab_type == 'VH sdAb':
        VH_seq = get_VHVL_seq(structure, VH_id)
        
        fasta = fasta + '>' + str(data.iloc[a][0]) + ', VH' + '\n' + VH_seq + '\n'
        
    elif Ab_type == 'VL sdAb':
        VL_seq = get_VHVL_seq(structure, VL_id)
        
        fasta = fasta + '>' + str(data.iloc[a][0]) + ', VH' + '\n' + VL_seq + '\n'
        
file2_out = open("Fasta.txt","a")
file2_out.write(fasta)


### Create for loop for counting total CDR amino acids in each structure
CDR_tot_list = []
VH_FR1_list = []
VH_CDR1_list = []
VH_FR2_list = []
VH_CDR2_list = []
VH_FR3_list = []
VH_CDR3_list = []
VH_FR4_list = []
VL_FR1_list = []
VL_CDR1_list = []
VL_FR2_list = []
VL_CDR2_list = []
VL_FR3_list = []
VL_CDR3_list = []
VL_FR4_list = []

for e in range(data.shape[0]):

    CDR_tot = 0
    
    VH_FR1 = 0
    VH_CDR1 = 0
    VH_FR2 = 0
    VH_CDR2 = 0
    VH_FR3 = 0
    VH_CDR3 = 0
    VH_FR4 = 0
    
    VL_FR1 = 0
    VL_CDR1 = 0
    VL_FR2 = 0
    VL_CDR2 = 0
    VL_FR3 = 0
    VL_CDR3 = 0
    VL_FR4 = 0
    
    pdb3 = str('./imgt_all_clean/' + str(data['pdb'][e]) + '.pdb')
    structure = parser.get_structure("structure", pdb3)
    
    #print('Counting CDR of: ' + str(data3['pdb'][e]) + ' (' + str(e) + '/' + str(data3.shape[0]) + ')')
    
    Ab_chain_list = []
    
    VH_id = str(data.iloc[e][1])
    VL_id = str(data.iloc[e][2])
    
    if VH_id.find('nan') < 0 and VL_id.find('nan') < 0:         ### Append the Ab chains to the chain list
        Ab_type = 'Fv'
        VH_id = str(data.iloc[e][1])
        VL_id = str(data.iloc[e][2])
        Ab_chain_list.append(VH_id)
        Ab_chain_list.append(VL_id)
        
        for residue_C in structure[0][VH_id]:
            if any(str(residue_C.get_resname()) == x for x in aa_list):         ### Only include if the residue is an amino acid   
                    
                if residue_C.get_id()[1] < 27:
                    VH_FR1 = VH_FR1+1
                elif 27 <= residue_C.get_id()[1] <= 38:
                    VH_CDR1 = VH_CDR1+1
                elif 39 <= residue_C.get_id()[1] <= 55:
                    VH_FR2 = VH_FR2+1
                elif 56 <= residue_C.get_id()[1] <= 65:
                    VH_CDR2 = VH_CDR2+1
                elif 66 <= residue_C.get_id()[1] <= 104:
                    VH_FR3 = VH_FR3+1
                elif 105 <= residue_C.get_id()[1] <= 117:
                    VH_CDR3 = VH_CDR3+1
                elif 118 <= residue_C.get_id()[1] <= 128:
                    VH_FR4 = VH_FR4+1
                    
        for residue_C in structure[0][VL_id]:
            if any(str(residue_C.get_resname()) == x for x in aa_list):         ### Only include if the residue is an amino acid   
                    
                if residue_C.get_id()[1] < 27:
                    VL_FR1 = VL_FR1+1
                elif 27 <= residue_C.get_id()[1] <= 38:
                    VL_CDR1 = VL_CDR1+1
                elif 39 <= residue_C.get_id()[1] <= 55:
                    VL_FR2 = VL_FR2+1
                elif 56 <= residue_C.get_id()[1] <= 65:
                    VL_CDR2 = VL_CDR2+1
                elif 66 <= residue_C.get_id()[1] <= 104:
                    VL_FR3 = VL_FR3+1
                elif 105 <= residue_C.get_id()[1] <= 117:
                    VL_CDR3 = VL_CDR3+1
                elif 118 <= residue_C.get_id()[1] <= 128:
                    VL_FR4 = VL_FR4+1
        
    elif VH_id.find('nan') < 0 and VL_id.find('nan') >= 0:
        Ab_type = 'VH sdAb'
        VH_id = str(data.iloc[e][1])
        Ab_chain_list.append(VH_id)
        
        for residue_C in structure[0][VH_id]:
            if any(str(residue_C.get_resname()) == x for x in aa_list):         ### Only include if the residue is an amino acid   
                    
                if residue_C.get_id()[1] < 27:
                    VH_FR1 = VH_FR1+1
                elif 27 <= residue_C.get_id()[1] <= 38:
                    VH_CDR1 = VH_CDR1+1
                elif 39 <= residue_C.get_id()[1] <= 55:
                    VH_FR2 = VH_FR2+1
                elif 56 <= residue_C.get_id()[1] <= 65:
                    VH_CDR2 = VH_CDR2+1
                elif 66 <= residue_C.get_id()[1] <= 104:
                    VH_FR3 = VH_FR3+1
                elif 105 <= residue_C.get_id()[1] <= 117:
                    VH_CDR3 = VH_CDR3+1
                elif 118 <= residue_C.get_id()[1] <= 128:
                    VH_FR4 = VH_FR4+1
                
    elif VH_id.find('nan') >= 0 and VL_id.find('nan') < 0:
        Ab_type = 'VL sdAb'
        VL_id = str(data.iloc[e][2])
        Ab_chain_list.append(VL_id)
        
        for residue_C in structure[0][VL_id]:
            if any(str(residue_C.get_resname()) == x for x in aa_list):         ### Only include if the residue is an amino acid   
                    
                if residue_C.get_id()[1] < 27:
                    VL_FR1 = VL_FR1+1
                elif 27 <= residue_C.get_id()[1] <= 38:
                    VL_CDR1 = VL_CDR1+1
                elif 39 <= residue_C.get_id()[1] <= 55:
                    VL_FR2 = VL_FR2+1
                elif 56 <= residue_C.get_id()[1] <= 65:
                    VL_CDR2 = VL_CDR2+1
                elif 66 <= residue_C.get_id()[1] <= 104:
                    VL_FR3 = VL_FR3+1
                elif 105 <= residue_C.get_id()[1] <= 117:
                    VL_CDR3 = VL_CDR3+1
                elif 118 <= residue_C.get_id()[1] <= 128:
                    VL_FR4 = VL_FR4+1
        
    for chain_ID in Ab_chain_list:
        for residue_C in structure[0][chain_ID]:
            if any(str(residue_C.get_resname()) == x for x in aa_list):         ### Only include if the residue is an amino acid                        
                
                #VHVL_id = str(residue_C.get_id()[1]) + str(residue_C.get_id()[2])
                #VHVL_id_list.append(VHVL_id.replace(" ", ""))
                #VHVL_id_list_number.append(float(residue_C.get_id()[1]))
                #VHVL_id_list_letter.append(str(residue_C.get_id()[2]))
                
                if 27 <= residue_C.get_id()[1] <= 38 or 56 <= residue_C.get_id()[1] <= 65 or 105 <= residue_C.get_id()[1] <= 117:
                    CDR_tot = CDR_tot + 1
                    
    CDR_tot_list.append(CDR_tot)
    VH_FR1_list.append(VH_FR1)
    VH_CDR1_list.append(VH_CDR1)
    VH_FR2_list.append(VH_FR2)
    VH_CDR2_list.append(VH_CDR2)
    VH_FR3_list.append(VH_FR3)
    VH_CDR3_list.append(VH_CDR3)
    VH_FR4_list.append(VH_FR4)
    VL_FR1_list.append(VL_FR1)
    VL_CDR1_list.append(VL_CDR1)
    VL_FR2_list.append(VL_FR2)
    VL_CDR2_list.append(VL_CDR2)
    VL_FR3_list.append(VL_FR3)
    VL_CDR3_list.append(VL_CDR3)
    VL_FR4_list.append(VL_FR4)

CDR_data = pd.DataFrame({'Total CDR residues':CDR_tot_list , 
                         '#aa in VH FR1 (ref)':VH_FR1_list,
                         '#aa in VH CDR1 (ref)':VH_CDR1_list,
                         '#aa in VH FR2 (ref)':VH_FR2_list,
                         '#aa in VH CDR2 (ref)':VH_CDR2_list,
                         '#aa in VH FR3 (ref)':VH_FR3_list,
                         '#aa in VH CDR3 (ref)':VH_CDR3_list,
                         '#aa in VH FR4 (ref)':VH_FR4_list,
                         '#aa in VL FR1 (ref)':VL_FR1_list,
                         '#aa in VL CDR1 (ref)':VL_CDR1_list,
                         '#aa in VL FR2 (ref)':VL_FR2_list,
                         '#aa in VL CDR2 (ref)':VL_CDR2_list,
                         '#aa in VL FR3 (ref)':VL_FR3_list,
                         '#aa in VL CDR3 (ref)':VL_CDR3_list,
                         '#aa in VL FR4 (ref)':VL_FR4_list}
                          )
data4 = pd.concat([data,CDR_data] , axis=1)
data4.to_csv('Summary_all_sorted.csv', sep=';',  index=False)
