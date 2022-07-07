from pathlib import Path
import argparse
from pandas import read_parquet, DataFrame
from shared_funcs import get_ab_type


def main(summary_data_path: Path, cluster_path: Path, output_path: Path, csv_output: bool):
	summary_data = read_parquet(summary_data_path)
	file1_list = []
	pdb_list = []
	resolution_list = []
	affinity_list = []
	ab_type_list = []
	cluster_id = -1

	for i, row in summary_data.iterrows():
		pdb_list.append(row.pdb)
		resolution_list.append(row.resolution)
		affinity_list.append(row.affinity)
		ab_type = get_ab_type(row[['Hchain', 'Lchain']].isna())
		ab_type_list.append(ab_type)

	# Make a list with empty lists for storing cluster IDs in. Each list corresponds to the specific pdb identifier
	cluster_id_list = [[] for _ in range(len(pdb_list))]

	with open(cluster_path) as a_file:
		# Go through all lines of cluster file outputted from CD-HIT
		for line in a_file:
			if line.count('>Cluster'):
				cluster_id += 1  # If new cluster is identified add to counter

			# In the specific row go over all the pdb entries and see if any of them match.
			# If the pdb entry is found on the row this specific cluster id is added to the
			# list with index mathing the pdb.
			for i in range(len(pdb_list)):
				if line.count(pdb_list[i]) > 0:
					cluster_id_list[i].append(cluster_id)

	final_list = []

	# The cluster_id_list contains a list of lists. Each of these lists are representing a unique pdb.
	# Each element represent the alignment cluster that the antibody chains from the pdb falls into
	cluster_id_list_unique = []

	for x in cluster_id_list:
		if x not in cluster_id_list_unique:
			cluster_id_list_unique.append(x)

	# Go through all unique cluster id combinations. For each combination look through the full list of combinations and store the index of this combination.
	# The point is to see if this specific combination of clusters appear more than once. The reason for doing this is to allow inclusion of common light chain antibodies.
	# So even if a VH1 show similarity to another VH2 in a given cluster. If VL1 does not share similarity with VL2 in a corresponding cluster then the antibody is considered unique.
	for cluster_combi in cluster_id_list_unique:
		index_list = [i for i, a in enumerate(cluster_id_list) if a == cluster_combi]

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

	for i, row in summary_data.iterrows():
		if final_list.count(row.pdb) == 1:
			file1_list.append(row)

	data2 = DataFrame(file1_list)
	data2.to_parquet(output_path / 'Summary_all_sorted_nonredundant.parquet')
	if csv_output:
		data2.to_csv(output_path / 'Summary_all_sorted_nonredundant.parquet')


if __name__ == '__main__':
	aparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	aparser.add_argument("-s", default=Path('test/work_dir/Summary_all_sorted.parquet'), type=Path, help="Path to Summary_all_sorted.parquet")
	aparser.add_argument("-c", default=Path('test/work_dir/Cluster.txt'), type=Path, help="Path to Cluster.txt")
	aparser.add_argument("-o", default=Path('.'), type=Path, help="Path to write output to")
	aparser.add_argument("--csv", action='store_true', help="Also output csv format")
	args = aparser.parse_args()
	main(args.s, args.c, args.o, args.csv)
