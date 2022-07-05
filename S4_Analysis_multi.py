import pandas as pd
from multiprocessing import Pool
from multiprocessing import Manager

# Getting a list of the PDB names
#pdb_list = data['PDB'].tolist()
#pdb_list = list(dict.fromkeys(pdb_list))

#aa_list = ['ARG','HIS','LYS','ASP','GLU','SER','THR','ASN','GLN','CYS','GLY','PRO','ALA','ILE','LEU','MET','PHE','TRP','TYR','VAL']
#Ab_domains = ['VH FR1' , 'VH CDR1', 'VH FR2', 'VH CDR2', 'VH FR3', 'VH CDR3', 'VH FR4', 'VL FR1' , 'VL CDR1', 'VL FR2', 'VL CDR2', 'VL FR3', 'VL CDR3', 'VL FR4']

# Set number of multiprocesses
mul = 8

# Define functions
def get_Ab_type(df):
    
    # Get a list of the antibody types
    Ab_type_list = df['Ab type'].tolist()
    Ab_type_list = list(dict.fromkeys(Ab_type_list))

    if len(Ab_type_list) > 1:
        print(Ab_type_list)
        print('The pdb entry contains more than one Ab type!?')
        
    return Ab_type_list[0]

def get_target_type(df):
    # Get a list of the antibody types
    target_type_list = df['Antigen type'].tolist()
    target_type_list = list(dict.fromkeys(target_type_list))
    
    if len(target_type_list) > 1:
        print(target_type_list)
        print('The pdb entry contains more than one target type!?')
        
    return target_type_list[0]

def unique_contact_resi(df, Ab_resi_specific_filter_list):
    
    # Function counts total CDR contacts by dropping rows that are duplicate on antibody chain, residue ID, residue name and CDR. 
    func_df = df.drop_duplicates(subset=Ab_resi_specific_filter_list, keep='last')
    
    return func_df.shape[0]

def unique_contact_atom(df, Ab_atom_specific_filter_list):
    
    # Function counts total CDR contacts by dropping rows that are duplicate on antibody chain, residue ID, residue name and CDR. 
    func_df = df.drop_duplicates(subset=Ab_atom_specific_filter_list, keep='last')
    
    return func_df.shape[0]

def epitope_resi_contact(df, Ag_resi_specific_filter_list):

    func_df = df.drop_duplicates(subset=Ag_resi_specific_filter_list, keep='last')
    
    return func_df.shape[0]

def epitope_atom_contact(df, Ag_atom_specific_filter_list):
    
    func_df = df.drop_duplicates(subset=Ag_atom_specific_filter_list, keep='last')
    
    return func_df.shape[0]

def Ab_aa_count(df, residue, Ab_resi_specific_filter_list):
    
    # The function calculates the count of a specific amino acid in the paratope data

    # Get unique residues on the antibody
    func_df = df.drop_duplicates(subset=Ab_resi_specific_filter_list, keep='last')
    
    # Get all the rows where the amino acid is the same as specified by user. 
    func_df2 = func_df.loc[func_df['Ab resi aa'] == residue]
    
    return func_df2.shape[0]
    
def Ag_aa_count(df, residue, Ag_resi_specific_filter_list):
    
    # The function calculates the count of a specific amino acid in the epitope data

    # Get unique residues on the antigen
    func_df = df.drop_duplicates(subset=Ag_resi_specific_filter_list, keep='last')
    
    # Get all the rows where the amino acid is the same as specified by user. 
    func_df2 = func_df.loc[func_df['Ag resi aa'] == residue]
    
    return func_df2.shape[0]

def total_resi_in_each_domain(df, Ab_domain, Ab_resi_specific_filter_list):
    
    # Drop duplicates in the data to only get unique residues
    df_NR = df.drop_duplicates(subset=Ab_resi_specific_filter_list, keep='last')
    
    # The chain (VH or VL) is specified as the first two letters in the domain (e.g. "VH FR1"
    chain = Ab_domain[0:2]
    
    # The domain, e.g. CDR1 is specified similar to above only it is the last part of the string.
    domain = Ab_domain[3:]
    
    data_slice1 = df_NR[df_NR['Ab chain type (VH/VL)'] == chain]
    
    data_slice2 = data_slice1[data_slice1['Ab domain'] == domain]
    
    domain_count = data_slice2.shape[0]
    
    return domain_count
    
def total_atom_in_each_domain(df, Ab_domain, Ab_atom_specific_filter_list):
    
    # Drop duplicates in the data to only get unique residues
    df_NR = df.drop_duplicates(subset=Ab_atom_specific_filter_list, keep='last')
    
    # The chain (VH or VL) is specified as the first two letters in the domain (e.g. "VH FR1"
    chain = Ab_domain[0:2]
    
    # The domain, e.g. CDR1 is specified similar to above only it is the last part of the string.
    domain = Ab_domain[3:]
    
    data_slice1 = df_NR[df_NR['Ab chain type (VH/VL)'] == chain]
    
    data_slice2 = data_slice1[data_slice1['Ab domain'] == domain]
    
    domain_count = data_slice2.shape[0]
    
    return domain_count

def total_contacts_in_domain(df, Ab_domain):
    
    # The chain (VH or VL) is specified as the first two letters in the domain (e.g. "VH FR1"
    chain = Ab_domain[0:2]
    
    # The domain, e.g. CDR1 is specified similar to above only it is the last part of the string.
    domain = Ab_domain[3:]
    
    data_slice1 = data[data['Ab chain type (VH/VL)'] == chain]
    
    data_slice2 = data_slice1[data_slice1['Ab domain'] == domain]
    
    domain_count = data_slice2.shape[0]
    
    return domain_count

def aa_freq_in_domains(df, Ab_domain, aa, Ab_resi_specific_filter_list):
    
    # The chain (VH or VL) is specified as the first two letters in the domain (e.g. "VH FR1"
    chain = Ab_domain[0:2]
    
    # The domain, e.g. CDR1 is specified similar to above only it is the last part of the string.
    domain = Ab_domain[3:]
    
    df_NR = df.drop_duplicates(subset=Ab_resi_specific_filter_list)
    
    data_slice1 = df_NR[df_NR['Ab chain type (VH/VL)'] == chain]
    
    data_slice2 = data_slice1[data_slice1['Ab domain'] == domain]
    
    data_slice3 = data_slice2[data_slice2['Ab resi aa'] == aa]
    
    aa_count = data_slice3.shape[0]
    
    #aa_count = Ab_aa_count(data_slice, aa, Ab_resi_specific_filter_list)
    
    return aa_count
    
def get_chain_organisms(df):
    
    H_chain_organism = df['H chain organism'].tolist()
    L_chain_organism = df['L chain organism'].tolist()
    
    H_chain_organism = list(dict.fromkeys(H_chain_organism))
    L_chain_organism = list(dict.fromkeys(L_chain_organism))
    
    if len(H_chain_organism) > 1:
        print('More than one H chain organism found!')
    
    if len(L_chain_organism) > 1:
        print('More than one L chain organism found')
    
    return H_chain_organism[0], L_chain_organism[0]

def get_domain_lengths_ref(df, Ab_domain):
    
    count = df['#aa in ' + Ab_domain + ' (ref)'].tolist()
    count = list(dict.fromkeys(count))
    
    return count[0]


   
# Define master dict
#master_dict = {}
#master_dict['pdb_name'] = []
#master_dict['Ab type'] = []
#master_dict['PDB file group'] = []
#master_dict['unique #aa in Ab'] = []
#master_dict['unique #atom in Ab'] = []
#master_dict['unique #aa in Ag'] = []
#master_dict['unique #atom in Ag'] = []

#for a in aa_list:
#    master_dict['#'+str(a)+' in Ab'] = []

#for a in aa_list:
#    master_dict['#'+str(a)+' in Ag'] = []  

#for a in aa_list:
#    for b in Ab_domains:
#        master_dict['#'+str(a)+' in '+str(b)] = []


def main_func(my_tuple):
    
    # Defining local variables
    Ab_resi_specific_filter_list = ['PDB', 'Ab chain type (VH/VL)', 'Ab resi ID', 'Ab resi aa']
    Ab_atom_specific_filter_list = ['PDB', 'Ab chain type (VH/VL)', 'Ab resi ID', 'Ab resi aa', 'Ab atom']
    Ag_resi_specific_filter_list = ['PDB', 'Ag chain ID', 'Ag resi ID', 'Ag resi aa']
    Ag_atom_specific_filter_list = ['PDB', 'Ag chain ID', 'Ag resi ID', 'Ag resi aa', 'Ag atom']
    Ab_domains = ['VH FR1' , 'VH CDR1', 'VH FR2', 'VH CDR2', 'VH FR3', 'VH CDR3', 'VH FR4', 'VL FR1' , 'VL CDR1', 'VL FR2', 'VL CDR2', 'VL FR3', 'VL CDR3', 'VL FR4']
    aa_list = ['ARG','HIS','LYS','ASP','GLU','SER','THR','ASN','GLN','CYS','GLY','PRO','ALA','ILE','LEU','MET','PHE','TRP','TYR','VAL']
    
    # Defining lists for output
    master_dict = {}
    master_dict['Antigen type'] = []
    master_dict['PDB'] = []
    master_dict['Ab type'] = []
    master_dict['Hchain species'] = []
    master_dict['Lchain species'] = []
    master_dict['No. of total contacts'] = []
    master_dict['No. paratope residues'] = []
    master_dict['No. paratope atoms'] = []
    master_dict['No. epitope residues'] = []
    master_dict['No. epitope atoms'] = []
    
    master_dict['#Resi in VH'] = []
    master_dict['#Resi in VL'] = []
    
    for a in Ab_domains:
        master_dict['#Resi in ' + str(a)] = []
        
    for a in Ab_domains:
        master_dict['#Resi in ' + str(a) + ' (ref)'] = []
    
    for a in Ab_domains:
        master_dict['#Atom in ' + str(a)] = []
    
    for a in Ab_domains:
        master_dict['Total contacts in ' + str(a)] = []
    
    for a in aa_list:
        master_dict['#'+str(a)+' in Ab'] = []
        master_dict['#'+str(a)+' in VH'] = []
        master_dict['#'+str(a)+' in VL'] = []
        master_dict['#'+str(a)+' in Ag'] = []
        
    for a in aa_list:
        for b in Ab_domains:
            master_dict['#'+str(a)+' in '+str(b)] = []
    
    pdb_id = my_tuple[1]
    master_dict['PDB'].append(pdb_id)
    #pdb_list.append(pdb_id)
    print(pdb_id)
    
    #print(line + '/'+ str(len(pdb_list))+ ' (' + str(pdb_list[i]) + ')')
    
    data_new = data[data['PDB'] == pdb_id]
    data_new_ref = data_ref[data_ref['PDB'] == pdb_id] 
    
    # Get the count of residues in the different domains in the antibody
    data_new_VH = data_new.loc[data_new['Ab chain type (VH/VL)'] == 'VH']
    data_new_VL = data_new.loc[data_new['Ab chain type (VH/VL)'] == 'VL']
    
    #pdb_name = pdb_id
    #master_dict['pdb_name'].append(pdb_name)
    
    # Getting the antibody type for the specific pdb
    Ab_type = get_Ab_type(data_new)
    master_dict['Ab type'].append(Ab_type)
    #Ab_type_list.append(Ab_type)
    #master_dict['Ab type'].append(Ab_type)
    
    # Getting the target type for the specific pdb. This can be protein or peptide
    target_type = get_target_type(data_new)
    master_dict['Antigen type'].append(target_type)
    #master_dict['PDB file group'].append(target_type)
    
    # Getting species for each of the chains
    H_chain_organism, L_chain_organism = get_chain_organisms(data_new)
    master_dict['Hchain species'].append(H_chain_organism)
    master_dict['Lchain species'].append(L_chain_organism)
    
    # Get the total number of unique contacts in the given structure
    total_contacts = data_new.shape[0]
    master_dict['No. of total contacts'].append(total_contacts)
    
    # Get the total number of unique antibody amino acid residues in the interface
    Ab_contact_resi = unique_contact_resi(data_new, Ab_resi_specific_filter_list)
    master_dict['No. paratope residues'].append(Ab_contact_resi)
    
    # Get the total number of unique antibody amino acid residues in the interface in each of the chains
    Ab_contact_resi_VH = unique_contact_resi(data_new_VH, Ab_resi_specific_filter_list)
    master_dict['#Resi in VH'].append(Ab_contact_resi_VH)
    Ab_contact_resi_VL = unique_contact_resi(data_new_VL, Ab_resi_specific_filter_list)
    master_dict['#Resi in VL'].append(Ab_contact_resi_VL)
    
    # Get the total number of unique antibody atoms in the interface
    Ab_contact_atom = unique_contact_atom(data_new, Ab_atom_specific_filter_list)
    master_dict['No. paratope atoms'].append(Ab_contact_atom)
    #Ab_contact_atom_list.append(Ab_contact_atom)
    #master_dict['unique #atom in Ab'].append(Ab_contact_atom)
    
    # Get total number of unique residues on the epitope
    epi_contact_resi = epitope_resi_contact(data_new, Ag_resi_specific_filter_list)
    master_dict['No. epitope residues'].append(epi_contact_resi)
    #Ag_contact_resi_list.append(epi_contact_resi)
    #master_dict['unique #aa in Ag'].append(epi_contact_resi)
    
    # Get total number of unique atoms on the epitope
    epi_contact_atom = epitope_atom_contact(data_new, Ag_atom_specific_filter_list)
    master_dict['No. epitope atoms'].append(epi_contact_atom)
    #Ag_contact_atom_list.append(epi_contact_atom)
    #master_dict['unique #atom in Ag'].append(epi_contact_atom)
    
    # Get total count of residues in antibody paratope
    for a in aa_list:
        aa_count = Ab_aa_count(data_new, a, Ab_resi_specific_filter_list)
        master_dict['#'+str(a)+' in Ab'].append(aa_count)
        
        aa_count_VH = Ab_aa_count(data_new_VH, a, Ab_resi_specific_filter_list)
        master_dict['#'+str(a)+' in VH'].append(aa_count_VH)
        
        aa_count_VL = Ab_aa_count(data_new_VL, a, Ab_resi_specific_filter_list)
        master_dict['#'+str(a)+' in VL'].append(aa_count_VL)
        
    # Get total count of residues in the antigen epitope
    for a in aa_list:
        aa_count = Ag_aa_count(data_new, a, Ag_resi_specific_filter_list)
        master_dict['#' + str(a) + ' in Ag'].append(aa_count)
    
    # Get total count of unique residues in each of the paratope domains
    for a in Ab_domains:
        domain_count = total_resi_in_each_domain(data_new, a, Ab_resi_specific_filter_list)
        master_dict['#Resi in ' + str(a)].append(domain_count)
        
    # Get total count of residues in the reference data i.e. the domain lengths irrespective of if they are contact or not. 
    for a in Ab_domains:
        domain_length = get_domain_lengths_ref(data_new, a)
        master_dict['#Resi in ' + str(a) + ' (ref)'].append(domain_length)
    
    # Get total count of unique atoms in each of the paratope domains   
    for a in Ab_domains:
        domain_count = total_atom_in_each_domain(data_new, a, Ab_atom_specific_filter_list)
        master_dict['#Atom in ' + str(a)].append(domain_count)
    
    # Get total count of contacts in each of the paratope domains  
    for a in Ab_domains:
        domain_count = total_contacts_in_domain(data_new, a)
        master_dict['Total contacts in ' + str(a)].append(domain_count)
    
    # Get count of the different amino acids in the different antibody domains of both VH and VL
    for a in Ab_domains:
        #chain = a[:1]
        #domain = a[len(a)-3:]
        for b in aa_list:
            my_var = aa_freq_in_domains(data_new, a, b, Ab_resi_specific_filter_list)
            master_dict['#'+str(b)+' in '+str(a)].append(my_var)
            
    
    df_out = pd.DataFrame.from_dict(master_dict)
    
    return df_out

# Loading the csv file containing data
# Oh no... I am defining dataframes as global variables... Bad idea...? But they are not changed anywhere in the main function, only used for making sub-dataframes
data_ref = pd.read_csv('Output_ref.csv',delimiter=';')
data = pd.read_csv('Output_contact.csv',delimiter=';')
pdb_data = data['PDB'].drop_duplicates()

if __name__ == '__main__':

    pool = Pool(mul)
    results = pool.map(main_func, pdb_data.iteritems())
    
    results_analysis = pd.concat(results)
    results_analysis.to_csv('Output_Analysis.csv', sep=';')


