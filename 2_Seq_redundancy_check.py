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
    
    if VH_id.find('nan') < 0 and VL_id.find('nan') < 0:
        Ab_type = 'Fv'
    elif VH_id.find('nan') < 0 and VL_id.find('nan') >= 0:
        Ab_type = 'VH sdAb' 
    elif VH_id.find('nan') >= 0 and VL_id.find('nan') < 0:
        Ab_type = 'VL sdAb'  
    Ab_type_list.append(Ab_type)

Cluster_ID_list = [ [] for _ in range(len(pdb_list)) ]

with open("Cluster.txt", "r") as a_file:
    for line in a_file:
        if line.count('>Cluster'):
            Cluster_ID = Cluster_ID + 1
            
        for i in range(len(pdb_list)):
            if line.count(pdb_list[i]) > 0:
                Cluster_ID_list[i].append(Cluster_ID)

Used_cluster_combi = []
Final_list = []

Cluster_ID_list_unique = []
[Cluster_ID_list_unique.append(x) for x in Cluster_ID_list if x not in Cluster_ID_list_unique]

for j in range(len(Cluster_ID_list_unique)):
    
    Cluster_combi = Cluster_ID_list_unique[j]
    index_list = []
    
    for a in range(0,len(Cluster_ID_list)):
            if Cluster_ID_list[a] == Cluster_combi:
                index_list.append(a)
    
    if len(index_list) == 1:
        Final_list.append(pdb_list[index_list[0]])
    
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