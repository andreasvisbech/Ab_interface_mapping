from pandas import read_feather, DataFrame
from pathlib import Path


def main(summary_data_path: Path, cluster_path: Path, output_path: Path):
	summary_data = read_feather(summary_data_path)

	file1_list = []

	pdb_list = []
	resolution_list = []
	affinity_list = []
	ab_type_list = []
	cluster_id = -1

	for a in range(summary_data.shape[0]):
		pdb_list.append(summary_data.iloc[a][0])
		resolution_list.append(summary_data.iloc[a][16])
		affinity_list.append(summary_data.iloc[a][25])

		vh_id = str(summary_data.iloc[a][1])
		vl_id = str(summary_data.iloc[a][2])

		if vh_id != 'nan' and vl_id != 'nan':
			Ab_type = 'Fv'
		elif vh_id != 'nan' and vl_id == 'nan':
			Ab_type = 'VH sdAb'
		elif vh_id == 'nan' and vl_id != 'nan':
			Ab_type = 'VL sdAb'
		else:
			Ab_type = 'N/A'
			print('Ab type not recognized!')

		ab_type_list.append(Ab_type)

	# Make a list with empty lists for storing cluster IDs in. Each list corresponds to the specific pdb identifier
	cluster_id_list = [[] for _ in range(len(pdb_list))]

	with open(cluster_path) as a_file:
		# Go through all lines of cluster file outputted from CD-HIT
		for line in a_file:
			if line.count('>Cluster'):
				cluster_id += 1  # If new cluster is identified add to counter

			# In the specific line go over all the pdb entries and see if any of them match.
			# If the pdb entry is found on the line this specific cluster id is added to the list with index mathing
			# the pdb.
			for i in range(len(pdb_list)):
				if line.count(pdb_list[i]) > 0:
					cluster_id_list[i].append(cluster_id)

	final_list = []

	# The cluster_id_list contains a list of lists. Each of these lists are representing a unique pdb.
	# Each element represent the alignment cluster that the antibody chains from the pdb falls into
	cluster_id_list_unique = []
	[cluster_id_list_unique.append(x) for x in cluster_id_list if x not in cluster_id_list_unique]

	# Go through all unique cluster id combinations. For each combination look through the full list of combinations and store the index of this combination.
	# The point is to see if this specific combination of clusters appear more than once. The reason for doing this is to allow inclusion of common light chain antibodies.
	# So even if a VH1 show similarity to another VH2 in a given cluster. If VL1 does not share similarity with VL2 in a corresponding cluster then the antibody is considered unique.
	for j in range(len(cluster_id_list_unique)):
		cluster_combi = cluster_id_list_unique[j]
		index_list = []

		for a in range(0, len(cluster_id_list)):
			if cluster_id_list[a] == cluster_combi:
				index_list.append(a)

		# Only one entry in the index_list show the antibody does not share similarity with any other Abs in the data.
		if len(index_list) == 1:
			final_list.append(pdb_list[index_list[0]])

		# If more than one index is found redundant antibodies have been found.
		# If this is the case antibodies with affinity values are prioritized and otherwise we just take the structure with best resolution.
		elif len(index_list) > 1:
			data_with_aff = []
			for b in index_list:
				if affinity_list[b] != 'None':
					data_with_aff.append(b)

			if len(data_with_aff) == 0:
				candidate_index = 0
				min_resolution = 10
				for c in index_list:
					if resolution_list[c] < min_resolution:
						min_resolution = resolution_list[c]
						candidate_index = c

				final_list.append(pdb_list[candidate_index])

			elif len(data_with_aff) == 1:
				final_list.append(pdb_list[data_with_aff[0]])

			elif len(data_with_aff) > 1:
				candidate_index = 0
				min_resolution = 10
				for d in data_with_aff:
					if resolution_list[d] < min_resolution:
						min_resolution = resolution_list[d]
						candidate_index = d

				final_list.append(pdb_list[candidate_index])

	for i in range(summary_data.shape[0]):
		if final_list.count(summary_data['pdb'][i]) == 1:
			file1_list.append(summary_data.iloc[i].tolist())

	data2 = DataFrame(data=file1_list, columns=summary_data.columns.values)
	data2.to_feather(output_path / 'Summary_all_sorted_nonredundant.fea.zst', compression='zstd')


if __name__ == '__main__':
	main(Path("test/work_dir/Summary_all_sorted.fea.zst"), Path("test/work_dir/Cluster.txt"), Path("."))
