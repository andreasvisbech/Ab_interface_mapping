from multiprocessing import Pool
import pandas as pd
import argparse
from pathlib import Path
from shared_funcs import aa_list, domain_list


def get_ab_type(df):
	"""
	Get a list of the antibody types
	:param df: A dataframe containing "Ab type" column
	:return: The singular unique value of "Ab type"
	"""
	ab_type_list = df['Ab type'].unique()

	if len(ab_type_list) > 1:
		print(ab_type_list)
		print('The pdb entry contains more than one Ab type!?')

	return ab_type_list[0]


def get_target_type(df):
	"""
	Get a list of the antibody types
	:param df:
	:return:
	"""
	target_type_list = df['Antigen type'].unique()

	if len(target_type_list) > 1:
		print(target_type_list)
		print('The pdb entry contains more than one target type!?')

	return target_type_list[0]


def unique_contact_resi(df, ab_resi_specific_filter_list):
	"""
	Counts total CDR contacts by dropping rows that are duplicate on antibody chain, residue ID, residue name and CDR.
	:param df:
	:param ab_resi_specific_filter_list:
	:return:
	"""
	return len(df.drop_duplicates(subset=ab_resi_specific_filter_list))


def unique_contact_atom(df, ab_atom_specific_filter_list):
	"""
	Function counts total CDR contacts by dropping rows that are duplicate on antibody chain, residue ID, residue name and CDR.
	:param df:
	:param ab_atom_specific_filter_list:
	:return:
	"""
	return len(df.drop_duplicates(subset=ab_atom_specific_filter_list, keep='last'))


def epitope_resi_contact(df, ag_resi_specific_filter_list):
	"""

	:param df:
	:param ag_resi_specific_filter_list:
	:return:
	"""
	return len(df.drop_duplicates(subset=ag_resi_specific_filter_list, keep='last'))


def epitope_atom_contact(df, ag_atom_specific_filter_list):
	"""

	:param df:
	:param ag_atom_specific_filter_list:
	:return:
	"""
	return len(df.drop_duplicates(subset=ag_atom_specific_filter_list, keep='last'))


def Ab_aa_count(df, residue, ab_resi_specific_filter_list):
	"""
	The function calculates the count of a specific amino acid in the paratope data
	:param df:
	:param residue:
	:param ab_resi_specific_filter_list:
	:return:
	"""
	# Get unique residues on the antibody
	func_df = df.drop_duplicates(subset=ab_resi_specific_filter_list, keep='last')
	# Get all the rows where the amino acid is the same as specified by user.
	func_df2 = func_df.loc[func_df['Ab resi aa'] == residue]

	return len(func_df2)


def Ag_aa_count(df, residue, ag_resi_specific_filter_list):
	"""
	The function calculates the count of a specific amino acid in the epitope data
	:param df:
	:param residue:
	:param ag_resi_specific_filter_list:
	:return:
	"""
	# Get unique residues on the antigen
	func_df = df.drop_duplicates(subset=ag_resi_specific_filter_list, keep='last')

	# Get all the rows where the amino acid is the same as specified by user.
	func_df2 = func_df.loc[func_df['Ag resi aa'] == residue]

	return len(func_df2)


def total_in_each_domain(df, ab_domain, ab_specific_filter_list):
	"""

	:param df:
	:param ab_domain:
	:param ab_specific_filter_list:
	:return:
	"""
	# Drop duplicates in the data to only get unique residues
	df_nr = df.drop_duplicates(subset=ab_specific_filter_list, keep='last')
	# The chain (VH or VL) is specified as the first two letters in the domain (e.g. "VH FR1"
	chain = ab_domain[0:2]
	# The domain, e.g. CDR1 is specified similar to above only it is the last part of the string.
	domain = ab_domain[3:]
	data_slice1 = df_nr[df_nr['Ab chain type (VH/VL)'] == chain]
	data_slice2 = data_slice1[data_slice1['Ab domain'] == domain]
	return len(data_slice2)


def total_contacts_in_domain(data, ab_domain):
	"""

	:param data:
	:param ab_domain:
	:return:
	"""
	# The chain (VH or VL) is specified as the first two letters in the domain (e.g. "VH FR1"
	chain = ab_domain[0:2]
	# The domain, e.g. CDR1 is specified similar to above only it is the last part of the string.
	domain = ab_domain[3:]
	data_slice1 = data[data['Ab chain type (VH/VL)'] == chain]
	data_slice2 = data_slice1[data_slice1['Ab domain'] == domain]

	return len(data_slice2)


def aa_freq_in_domains(df, ab_domain, aa, ab_resi_specific_filter_list):
	"""

	:param df:
	:param ab_domain:
	:param aa:
	:param ab_resi_specific_filter_list:
	:return:
	"""
	# The chain (VH or VL) is specified as the first two letters in the domain (e.g. "VH FR1"
	chain = ab_domain[0:2]
	# The domain, e.g. CDR1 is specified similar to above only it is the last part of the string.
	domain = ab_domain[3:]
	df_nr = df.drop_duplicates(subset=ab_resi_specific_filter_list)
	data_slice1 = df_nr[df_nr['Ab chain type (VH/VL)'] == chain]
	data_slice2 = data_slice1[data_slice1['Ab domain'] == domain]
	data_slice3 = data_slice2[data_slice2['Ab resi aa'] == aa]
	return len(data_slice3)


def get_chain_organisms(df):
	"""

	:param df:
	:return:
	"""
	h_chain_organism = df['H chain organism'].unique()
	l_chain_organism = df['L chain organism'].unique()
	if len(h_chain_organism) > 1:
		print('More than one H chain organism found!')

	if len(l_chain_organism) > 1:
		print('More than one L chain organism found')

	return h_chain_organism[0], l_chain_organism[0]


def get_domain_lengths_ref(df, ab_domain):
	"""

	:param df:
	:param ab_domain:
	:return:
	"""
	count = df['#aa in ' + ab_domain + ' (ref)'].unique()
	return count[0]


def pool_runner(grouped_by_pdb_df):
	"""

	:param grouped_by_pdb_df:
	:return:
	"""
	ab_resi_specific_filter_list = ['PDB', 'Ab chain type (VH/VL)', 'Ab resi ID', 'Ab resi aa']
	ab_atom_specific_filter_list = ['PDB', 'Ab chain type (VH/VL)', 'Ab resi ID', 'Ab resi aa', 'Ab atom']
	ag_resi_specific_filter_list = ['PDB', 'Ag chain ID', 'Ag resi ID', 'Ag resi aa']
	ag_atom_specific_filter_list = ['PDB', 'Ag chain ID', 'Ag resi ID', 'Ag resi aa', 'Ag atom']

	# Defining lists for output
	master_dict = {
		'Antigen type': [], 'PDB': [], 'Ab type': [], 'Hchain species': [], 'Lchain species': [],
		'No. of total contacts': [], 'No. paratope residues': [], 'No. paratope atoms': [],
		'No. epitope residues': [], 'No. epitope atoms': [], '#Resi in VH': [], '#Resi in VL': []
	}

	for a in domain_list:
		master_dict['#Resi in ' + str(a)] = []
		master_dict['#Resi in ' + str(a) + ' (ref)'] = []
		master_dict['#Atom in ' + str(a)] = []
		master_dict['Total contacts in ' + str(a)] = []

	for a in aa_list:
		master_dict['#' + str(a) + ' in Ab'] = []
		master_dict['#' + str(a) + ' in VH'] = []
		master_dict['#' + str(a) + ' in VL'] = []
		master_dict['#' + str(a) + ' in Ag'] = []
		for b in domain_list:
			master_dict['#' + str(a) + ' in ' + str(b)] = []

	pdb_id = grouped_by_pdb_df[0]
	master_dict['PDB'].append(pdb_id)
	print(pdb_id)

	data_new = grouped_by_pdb_df[1][grouped_by_pdb_df[1]['kilde'] == 'data']

	# Get the count of residues in the different domains in the antibody
	data_new_vh = data_new.loc[data_new['Ab chain type (VH/VL)'] == 'VH']
	data_new_vl = data_new.loc[data_new['Ab chain type (VH/VL)'] == 'VL']

	# Getting the antibody type for the specific pdb
	ab_type = get_ab_type(data_new)
	master_dict['Ab type'].append(ab_type)

	# Getting the target type for the specific pdb. This can be protein or peptide
	target_type = get_target_type(data_new)
	master_dict['Antigen type'].append(target_type)

	# Getting species for each of the chains
	h_chain_organism, l_chain_organism = get_chain_organisms(data_new)
	master_dict['Hchain species'].append(h_chain_organism)
	master_dict['Lchain species'].append(l_chain_organism)

	# Get the total number of unique contacts in the given structure
	total_contacts = len(data_new)
	master_dict['No. of total contacts'].append(total_contacts)

	# Get the total number of unique antibody amino acid residues in the interface
	ab_contact_resi = unique_contact_resi(data_new, ab_resi_specific_filter_list)
	master_dict['No. paratope residues'].append(ab_contact_resi)

	# Get the total number of unique antibody amino acid residues in the interface in each of the chains
	ab_contact_resi_vh = unique_contact_resi(data_new_vh, ab_resi_specific_filter_list)
	master_dict['#Resi in VH'].append(ab_contact_resi_vh)
	ab_contact_resi_vl = unique_contact_resi(data_new_vl, ab_resi_specific_filter_list)
	master_dict['#Resi in VL'].append(ab_contact_resi_vl)

	# Get the total number of unique antibody atoms in the interface
	ab_contact_atom = unique_contact_atom(data_new, ab_atom_specific_filter_list)
	master_dict['No. paratope atoms'].append(ab_contact_atom)

	# Get total number of unique residues on the epitope
	epi_contact_resi = epitope_resi_contact(data_new, ag_resi_specific_filter_list)
	master_dict['No. epitope residues'].append(epi_contact_resi)

	# Get total number of unique atoms on the epitope
	epi_contact_atom = epitope_atom_contact(data_new, ag_atom_specific_filter_list)
	master_dict['No. epitope atoms'].append(epi_contact_atom)

	# Get total count of residues in antibody paratope
	for a in aa_list:
		aa_count = Ab_aa_count(data_new, a, ab_resi_specific_filter_list)
		master_dict['#' + str(a) + ' in Ab'].append(aa_count)

		aa_count_vh = Ab_aa_count(data_new_vh, a, ab_resi_specific_filter_list)
		master_dict['#' + str(a) + ' in VH'].append(aa_count_vh)

		aa_count_vl = Ab_aa_count(data_new_vl, a, ab_resi_specific_filter_list)
		master_dict['#' + str(a) + ' in VL'].append(aa_count_vl)

	# Get total count of residues in the antigen epitope
	for a in aa_list:
		aa_count = Ag_aa_count(data_new, a, ag_resi_specific_filter_list)
		master_dict['#' + str(a) + ' in Ag'].append(aa_count)

	# Get total count of unique residues in each of the paratope domains
	for a in domain_list:
		domain_count = total_in_each_domain(data_new, a, ab_resi_specific_filter_list)
		master_dict['#Resi in ' + str(a)].append(domain_count)

	# Get total count of residues in the reference data i.e. the domain lengths irrespective of if they are contact or not.
	for a in domain_list:
		domain_length = get_domain_lengths_ref(data_new, a)
		master_dict['#Resi in ' + str(a) + ' (ref)'].append(domain_length)

	# Get total count of unique atoms in each of the paratope domains
	for a in domain_list:
		domain_count = total_in_each_domain(data_new, a, ab_atom_specific_filter_list)
		master_dict['#Atom in ' + str(a)].append(domain_count)

	# Get total count of contacts in each of the paratope domains
	for a in domain_list:
		domain_count = total_contacts_in_domain(data_new, a)
		master_dict['Total contacts in ' + str(a)].append(domain_count)

	# Get count of the different amino acids in the different antibody domains of both VH and VL
	for a in domain_list:
		# chain = a[:1]
		# domain = a[len(a)-3:]
		for b in aa_list:
			my_var = aa_freq_in_domains(data_new, a, b, ab_resi_specific_filter_list)
			master_dict['#' + str(b) + ' in ' + str(a)].append(my_var)

	df_out = pd.DataFrame.from_dict(master_dict)

	return df_out


def main(output_ref: Path, output_contact: Path, output_path: Path, threads: int, csv_output: bool):
	"""

	:param output_ref:
	:param output_contact:
	:param output_path:
	:param threads:
	:param csv_output:
	:return:
	"""
	data = pd.read_parquet(output_contact)
	ref = pd.read_parquet(output_ref)
	data['kilde'] = 'data'
	ref['kilde'] = 'ref'

	results = Pool(threads).map(pool_runner, pd.concat([data, ref]).groupby('PDB'))

	results_analysis = pd.concat(results)
	results_analysis.to_parquet(output_path / 'Output_Analysis.parquet')

	if csv_output:
		results_analysis.to_csv(output_path / 'Output_Analysis.csv')


if __name__ == '__main__':
	aparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	aparser.add_argument("-r", default=Path("test/work_dir/Output_ref.parquet"), type=Path,	 help="Path to script3s Output_ref.parquet")
	aparser.add_argument("-c", default=Path("test/work_dir/Output_contact.parquet"), type=Path, help="Path to script3s Output_contact.parquet")
	aparser.add_argument("-o", default=Path('.'), type=Path, help="Path to write output to")
	aparser.add_argument("-t", default=4, type=int, help="Number of threads to use")
	aparser.add_argument("--csv", action='store_true', help="Also output csv format")
	args = aparser.parse_args()

	main(args.r, args.c, args.o, args.t, args.csv)
