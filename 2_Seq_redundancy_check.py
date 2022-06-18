import pandas as pd

Summary_data = pd.read_csv('Summary_all_sorted.csv', sep=';',header=0)

file1_list = []

pdb_list = []
Resolution_list = []
Affinity_list = []
Ab_type_list = []
Cluster_ID = -1

for a in range(Summary_data.shape[0]):
    pdb_list.append(Summary_data.iloc[a][0])
    Resolution_list.append(Summary_data.iloc[a][16])
    Affinity_list.append(Summary_data.iloc[a][25])
    
    VH_id = str(Summary_data.iloc[a][1])
    VL_id = str(Summary_data.iloc[a][2])
    
    if VH_id != 'nan' and VL_id != 'nan':
        Ab_type = 'Fv'
    elif VH_id != 'nan' and VL_id == 'nan':
        Ab_type = 'VH sdAb' 
    elif VH_id == 'nan' and VL_id != 'nan':
        Ab_type = 'VL sdAb'
    else:
        Ab_type = 'N/A'
        print('Ab type not recognized!')

    Ab_type_list.append(Ab_type)

# Make a list with empty lists for storing cluster IDs in. Each list corresponds to the specific pdb identifier
Cluster_ID_list = [ [] for _ in range(len(pdb_list)) ]

with open("Cluster.txt", "r") as a_file:

    # Go through all lines of cluster file outputted from CD-HIT
    for line in a_file:
        if line.count('>Cluster'):
            Cluster_ID = Cluster_ID + 1             # If new cluster is identified add to counter

        # In the specific line go over all the pdb entries and see if any of them match.
        # If the pdb entry is found on the line this specific cluster id is added to the list with index mathing
        # the pdb.
        for i in range(len(pdb_list)):
            if line.count(pdb_list[i]) > 0:
                Cluster_ID_list[i].append(Cluster_ID)

Used_cluster_combi = []
Final_list = []

# The Cluster_ID_list contains a list of lists. Each of these lists are representing a unique pdb. 
# Each element represent the alignment cluster that the antibody chains from the pdb falls into
Cluster_ID_list_unique = []
[Cluster_ID_list_unique.append(x) for x in Cluster_ID_list if x not in Cluster_ID_list_unique]

# Go through all unique cluster id combinations. For each combination look through the full list of combinations and store the index of this combination. 
# The point is to see if this specific combination of clusters appear more than once. The reason for doing this is to allow inclusion of common light chain antibodies. 
# So even if a VH1 show similarity to another VH2 in a given cluster. If VL1 does not share similarity with VL2 in a corresponding cluster then the antibody is considered unique. 
for j in range(len(Cluster_ID_list_unique)):
    
    Cluster_combi = Cluster_ID_list_unique[j]
    index_list = []
    
    for a in range(0,len(Cluster_ID_list)):
            if Cluster_ID_list[a] == Cluster_combi:
                index_list.append(a)
    
    # Only one entry in the index_list show the antibody does not share similarity with any other Abs in the data. 
    if len(index_list) == 1:
        Final_list.append(pdb_list[index_list[0]])
    
    # If more than one index is found redundant antibodies have been found. 
    # If this is the case antibodies with affinity values are prioritized and otherwise we just take the structure with best resolution.
    elif len(index_list) > 1:
        data_with_aff = []
        for b in index_list:
            if Affinity_list[b] != 'None':
                data_with_aff.append(b)
                
        if len(data_with_aff) == 0:
            candidate_index = 0
            Min_resolution = 10
            for c in index_list:
                if Resolution_list[c] < Min_resolution:
                    Min_resolution = Resolution_list[c]
                    candidate_index = c
        
            Final_list.append(pdb_list[candidate_index])
        
        elif len(data_with_aff) == 1:
            Final_list.append(pdb_list[data_with_aff[0]])
            
        elif len(data_with_aff) > 1:
            candidate_index = 0
            Min_resolution = 10
            for d in data_with_aff:
                if Resolution_list[d] < Min_resolution:
                    Min_resolution = Resolution_list[d]
                    candidate_index = d
        
            Final_list.append(pdb_list[candidate_index])

for i in range(Summary_data.shape[0]):
    if Final_list.count(Summary_data['pdb'][i]) == 1:
        file1_list.append(Summary_data.iloc[i].tolist())
        
data2 = pd.DataFrame(data=file1_list, columns=Summary_data.columns.values)
data2.to_csv('Summary_all_sorted_nonredundant.csv', sep=';' ,index=False)
