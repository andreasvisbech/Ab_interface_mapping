import pandas as pd
import numpy as np
from datetime import datetime
from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.DSSP import DSSP
from multiprocessing import Pool
from pathlib import Path
from itertools import repeat
import argparse
from shared_funcs import aa_list, get_ab_type


def get_Ag_ids(row):
	"""
	The function looks at the length of the input and gets the chain ids.

	:param row:
	:return:
	"""
	#  If crystal contains 1 antigen chain the length will be 1.
	#  If two antigen chains the length will be 5.
	#  If 3 antigen chain length will be 4
	func_list = []
	if len(row['antigen_chain']) == 1:
		func_list.append(row['antigen_chain'][0])
	elif len(row['antigen_chain']) == 5:
		func_list.append(row['antigen_chain'][0])
		func_list.append(row['antigen_chain'][4])
	elif len(row['antigen_chain']) == 9:
		func_list.append(row['antigen_chain'][0])
		func_list.append(row['antigen_chain'][4])
		func_list.append(row['antigen_chain'][8])
	elif len(row['antigen_chain']) == 13:
		func_list.append(row['antigen_chain'][0])
		func_list.append(row['antigen_chain'][4])
		func_list.append(row['antigen_chain'][8])
		func_list.append(row['antigen_chain'][12])

	return func_list


def get_chain_seq(structure, chain_id):
	"""
	The function gets the full sequence of the antibody chain in question.

	:param structure:
	:param chain_id:
	:return:
	"""

	ppb = PPBuilder()
	func_var = ''
	try:
		for pp in ppb.build_peptides(structure[0][chain_id]):
			func_var += str(pp.get_sequence())
	except KeyError:
		func_var = 'N/A'
	return func_var


def get_Ab_domain(resi):
	"""
	The function determines which antibody domain the residue belongs to (FR1, CDR1, FR2, CDR2, FR3, CDR3 or FR4).
	The numbering comes from https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html

	:param resi:
	:return:
	"""

	if 27 > resi.get_id()[1]:
		cdr_value = 'FR1'
	elif 27 <= resi.get_id()[1] <= 38:
		cdr_value = 'CDR1'
	elif 38 < resi.get_id()[1] < 56:
		cdr_value = 'FR2'
	elif 56 <= resi.get_id()[1] <= 65:
		cdr_value = 'CDR2'
	elif 65 < resi.get_id()[1] < 105:
		cdr_value = 'FR3'
	elif 105 <= resi.get_id()[1] <= 117:
		cdr_value = 'CDR3'
	elif 117 < resi.get_id()[1] <= 128:
		cdr_value = 'FR4'
	else:
		cdr_value = 'N/A'

	return cdr_value


def domain_count(structure, ab_type, vh_id, vl_id):
	out_dict = {}

	vh_fr1_ref_count = vh_cdr1_ref_count = vh_fr2_ref_count = vh_cdr2_ref_count = 0
	vh_fr3_ref_count = vh_cdr3_ref_count = vh_fr4_ref_count = 0

	vl_fr1_ref_count = vl_cdr1_ref_count = vl_fr2_ref_count = vl_cdr2_ref_count = 0
	vl_fr3_ref_count = vl_cdr3_ref_count = vl_fr4_ref_count = 0

	if ab_type == 'Fv':
		for ab_resi in structure[0][vh_id]:
			cdr_value = get_Ab_domain(ab_resi)
			if cdr_value == 'FR1':
				vh_fr1_ref_count += 1
			elif cdr_value == 'CDR1':
				vh_cdr1_ref_count += 1
			elif cdr_value == 'FR2':
				vh_fr2_ref_count += 1
			elif cdr_value == 'CDR2':
				vh_cdr2_ref_count += 1
			elif cdr_value == 'FR3':
				vh_fr3_ref_count += 1
			elif cdr_value == 'CDR3':
				vh_cdr3_ref_count += 1
			elif cdr_value == 'FR4':
				vh_fr4_ref_count += 1
		for ab_resi in structure[0][vl_id]:
			cdr_value = get_Ab_domain(ab_resi)
			if cdr_value == 'FR1':
				vl_fr1_ref_count += 1
			elif cdr_value == 'CDR1':
				vl_cdr1_ref_count += 1
			elif cdr_value == 'FR2':
				vl_fr2_ref_count += 1
			elif cdr_value == 'CDR2':
				vl_cdr2_ref_count += 1
			elif cdr_value == 'FR3':
				vl_fr3_ref_count += 1
			elif cdr_value == 'CDR3':
				vl_cdr3_ref_count += 1
			elif cdr_value == 'FR4':
				vl_fr4_ref_count += 1
	elif ab_type == 'VH sdAb':
		for ab_resi in structure[0][vh_id]:
			cdr_value = get_Ab_domain(ab_resi)
			if cdr_value == 'FR1':
				vh_fr1_ref_count += 1
			elif cdr_value == 'CDR1':
				vh_cdr1_ref_count += 1
			elif cdr_value == 'FR2':
				vh_fr2_ref_count += 1
			elif cdr_value == 'CDR2':
				vh_cdr2_ref_count += 1
			elif cdr_value == 'FR3':
				vh_fr3_ref_count += 1
			elif cdr_value == 'CDR3':
				vh_cdr3_ref_count += 1
			elif cdr_value == 'FR4':
				vh_fr4_ref_count += 1
	elif ab_type == 'VL sdAb':
		for ab_resi in structure[0][vl_id]:
			cdr_value = get_Ab_domain(ab_resi)
			if cdr_value == 'FR1':
				vl_fr1_ref_count += 1
			elif cdr_value == 'CDR1':
				vl_cdr1_ref_count += 1
			elif cdr_value == 'FR2':
				vl_fr2_ref_count += 1
			elif cdr_value == 'CDR2':
				vl_cdr2_ref_count += 1
			elif cdr_value == 'FR3':
				vl_fr3_ref_count += 1
			elif cdr_value == 'CDR3':
				vl_cdr3_ref_count += 1
			elif cdr_value == 'FR4':
				vl_fr4_ref_count += 1

	out_dict['vh_fr1_ref_count'] = vh_fr1_ref_count
	out_dict['vh_cdr1_ref_count'] = vh_cdr1_ref_count
	out_dict['vh_fr2_ref_count'] = vh_fr2_ref_count
	out_dict['vh_cdr2_ref_count'] = vh_cdr2_ref_count
	out_dict['vh_fr3_ref_count'] = vh_fr3_ref_count
	out_dict['vh_cdr3_ref_count'] = vh_cdr3_ref_count
	out_dict['vh_fr4_ref_count'] = vh_fr4_ref_count

	out_dict['vl_fr1_ref_count'] = vl_fr1_ref_count
	out_dict['vl_cdr1_ref_count'] = vl_cdr1_ref_count
	out_dict['vl_fr2_ref_count'] = vl_fr2_ref_count
	out_dict['vl_cdr2_ref_count'] = vl_cdr2_ref_count
	out_dict['vl_fr3_ref_count'] = vl_fr3_ref_count
	out_dict['vl_cdr3_ref_count'] = vl_cdr3_ref_count
	out_dict['vl_fr4_ref_count'] = vl_fr4_ref_count

	return out_dict


def get_resi_id(resi):
	"""
	This function determines the residue id based on the pdb structure supplied.
	It includes both the number (index 1) and a potential letter (index 2) which can be present in the pdb.
	Lastly the function strips any whitespace.

	:param resi:
	:return:
	"""
	return f"{resi.get_id()[1]}{resi.get_id()[2]}".strip()


def get_SS(dssp_dict, chain_id, resi):
	return dssp_dict[chain_id, (' ', resi.get_id()[1], resi.get_id()[2])][2]


def get_atom_coord(atom):
	return [atom.get_coord()[0], atom.get_coord()[1], atom.get_coord()[2]]


def get_resi_ASA(dssp_dict, chain_id, resi):
	"""
	Function calculates available surface area of the residue using the dssp

	:param dssp_dict:
	:param chain_id:
	:param resi:
	:return:
	"""
	return dssp_dict[chain_id, (' ', resi.get_id()[1], resi.get_id()[2])][3]


def Ab_raw_extract(structure, index):
	out_list = []

	try:
		for func_resi in structure[0][index]:
			# Check that the residue is an amino acid
			if func_resi.get_resname() in aa_list:
				func_var = func_resi.get_resname()
				func_var2 = str(func_resi.get_id()[1]) + str(func_resi.get_id()[2])
				func_var2 = str(func_var2).strip()
				func_var3 = func_var + func_var2
				out_list.append(func_var3)
	except KeyError:
		out_list.append('N/A')

	return out_list


def pool_runner(row: pd.Series, pdb_path: Path):
	row = row[1]  # row = (id, actual row)
	contact_cutoff = 5

	pdb_type_list = []
	pdb_list = []
	Ab_type_list = []
	Ab_chain_seq_list = []
	Ab_chain_type_list = []
	Ab_chain_id_list = []
	Ab_domain_list = []
	Ab_resi_id_list = []
	Ab_resi_aa_list = []
	Ab_atom_list = []
	Ab_atom_coord_list = []
	Ab_SS_list = []
	Ab_resi_ASA_list = []
	Ab_raw_extract_list = []
	Ag_chain_id_list = []
	Ag_resi_id_list = []
	Ag_resi_aa_list = []
	Ag_atom_list = []
	Ag_atom_coord_list = []
	Ag_SS_list = []
	Ag_resi_ASA_list = []
	Contact_distance_list = []
	resolution_list = []
	H_chain_organism_list = []
	L_chain_organism_list = []

	VH_FR1_ref_list = []
	VH_CDR1_ref_list = []
	VH_FR2_ref_list = []
	VH_CDR2_ref_list = []
	VH_FR3_ref_list = []
	VH_CDR3_ref_list = []
	VH_FR4_ref_list = []
	VL_FR1_ref_list = []
	VL_CDR1_ref_list = []
	VL_FR2_ref_list = []
	VL_CDR2_ref_list = []
	VL_FR3_ref_list = []
	VL_CDR3_ref_list = []
	VL_FR4_ref_list = []

	pdb_type_list_raw = []
	pdb_list_raw = []
	Ab_type_list_raw = []
	Ab_chain_type_list_raw = []
	Ab_chain_id_list_raw = []
	Ab_domain_list_raw = []
	Ab_resi_id_list_raw = []
	Ab_resi_aa_list_raw = []
	Ab_atom_list_raw = []
	Ab_atom_coord_list_raw = []

	now = datetime.now()
	start = now.strftime("%H:%M:%S")

	print('Running: ', row.pdb, ' --- start: ', str(start))

	# Loading in the pdb file
	pdb1 = pdb_path / f"{row.pdb}.pdb"
	parser = PDBParser()
	structure = parser.get_structure("structure", pdb1)
	dssp_dict = DSSP(structure[0], pdb1, dssp="./dssp")
	pdb_type = row['Antigen group']  # Get antigen type
	resolution = row['resolution']  # Get structure resolution
	vh_organism = row['heavy_species']  # Get structure H chain organism
	vl_organism = row['light_species']  # Get structure L chain organism
	pdb_id = row.pdb
	ag_id_list = get_Ag_ids(row)  # Get a list of antigen ids in the structure

	# Getting the antibody type and the associated chain ids for antibody chains
	vh_id = row.Hchain
	vl_id = row.Lchain
	ab_type = get_ab_type(row[['Hchain', 'Lchain']].isna())
	if ab_type == 'Fv':
		ab_id_list = [vh_id, vl_id]
	elif ab_type == 'VH sdAb':
		ab_id_list = [vh_id]
	elif ab_type == 'VL sdAb':
		ab_id_list = [vl_id]
	else:
		ab_id_list = []

	# Create for loop that iterates through all antibody chains
	for Ab_idx in ab_id_list:
		ab_chain_seq = get_chain_seq(structure, Ab_idx)  # Get sequence of Ab chain and put into chain id list
		ab_raw = Ab_raw_extract(structure, Ab_idx)  # Get raw extract of the antibody chain
		# Get the antibody chain type i.e. if it is heavy chain or light chain
		if Ab_idx == vh_id:
			ab_chain_type = 'VH'
		elif Ab_idx == vl_id:
			ab_chain_type = 'VL'
		else:
			ab_chain_type = 'N/A'

		# Create for loop to go over all residues from the Ab_idx chain
		for Ab_resi in structure[0][Ab_idx]:
			# Check that the residue is an amino acid
			if Ab_resi.get_resname() in aa_list:
				ab_domain = get_Ab_domain(Ab_resi)  # Get the Ab domain for the residue
				ab_resi_id = get_resi_id(Ab_resi)  # Get the id of the antibody residue
				ab_resi_aa = Ab_resi.get_resname()  # Get the aa of the specific Ab residue
				try:  # Get residue secondary structure
					ab_resi_ss = get_SS(dssp_dict, Ab_idx, Ab_resi)
				except KeyError:
					ab_resi_ss = 'Missing atoms in residue'
					continue

				ab_resi_asa = get_resi_ASA(dssp_dict, Ab_idx, Ab_resi)  # Get antibody available surface area for the residue
				for Ab_atom in Ab_resi:  # Create for loop that goes through all atoms in the Ab residue
					ab_atom_id = Ab_atom.get_id()  # Get Ab atom
					ab_atom_coord = get_atom_coord(Ab_atom)  # Get Ab atom coords

					# Save the raw reference data irrespective of whether it is contact or not
					pdb_type_list_raw.append(pdb_type)
					pdb_list_raw.append(pdb_id)
					Ab_type_list_raw.append(ab_type)
					Ab_chain_type_list_raw.append(ab_chain_type)
					Ab_chain_id_list_raw.append(Ab_idx)
					Ab_domain_list_raw.append(ab_domain)
					Ab_resi_id_list_raw.append(ab_resi_id)
					Ab_resi_aa_list_raw.append(ab_resi_aa)
					Ab_atom_list_raw.append(ab_atom_id)
					Ab_atom_coord_list_raw.append(ab_atom_coord)

					# Create for loops that go through all atoms in antigen
					for Ag_idx in ag_id_list:
						for Ag_resi in structure[0][Ag_idx]:
							if Ag_resi.get_resname() in aa_list:  # Only consider amino acids
								ag_resi_aa = Ag_resi.get_resname()  # Get the aa of the specific Ab residue
								ag_resi_id = get_resi_id(Ag_resi)  # Get the id of the antibody residue
								try:  # Get residue secondary structure
									ag_resi_ss = get_SS(dssp_dict, Ag_idx, Ag_resi)
								except KeyError:
									ag_resi_ss = 'Missing atoms in residue'
									continue

								# Get antibody available surface area for the residue
								ag_resi_asa = get_resi_ASA(dssp_dict, Ag_idx, Ag_resi)
								for Ag_atom in Ag_resi:
									ag_atom_id = Ag_atom.get_id()  # Get Ag atom
									ag_atom_coord = get_atom_coord(Ag_atom)  # Get Ab atom coords
									# Determine distance between antibody atom and antigen atom and see if the distance
									# is below user specified cutoff i.e. if the atom-atom pair is considered a contact
									atom_distance = np.linalg.norm(Ab_atom.coord - Ag_atom.coord)

									if atom_distance <= contact_cutoff and ab_atom_id[0] != 'H' and ag_atom_id[0] != 'H':
										pdb_type_list.append(pdb_type)
										pdb_list.append(pdb_id)
										Ab_type_list.append(ab_type)
										Ab_chain_seq_list.append(ab_chain_seq)
										Ab_chain_type_list.append(ab_chain_type)
										Ab_chain_id_list.append(Ab_idx)
										Ab_domain_list.append(ab_domain)
										Ab_resi_id_list.append(ab_resi_id)
										Ab_resi_aa_list.append(ab_resi_aa)
										Ab_atom_list.append(ab_atom_id)
										Ab_atom_coord_list.append(ab_atom_coord)
										Ab_SS_list.append(ab_resi_ss)
										Ab_resi_ASA_list.append(ab_resi_asa)
										Ab_raw_extract_list.append(ab_raw)

										Ag_chain_id_list.append(Ag_idx)
										Ag_resi_id_list.append(ag_resi_id)
										Ag_resi_aa_list.append(ag_resi_aa)
										Ag_atom_list.append(ag_atom_id)
										Ag_atom_coord_list.append(ag_atom_coord)
										Ag_SS_list.append(ag_resi_ss)
										Ag_resi_ASA_list.append(ag_resi_asa)

										Contact_distance_list.append(atom_distance)
										resolution_list.append(resolution)
										H_chain_organism_list.append(vh_organism)
										L_chain_organism_list.append(vl_organism)

										VH_FR1_ref_list.append(row['#aa in VH FR1 (ref)'])
										VH_CDR1_ref_list.append(row['#aa in VH CDR1 (ref)'])
										VH_FR2_ref_list.append(row['#aa in VH FR2 (ref)'])
										VH_CDR2_ref_list.append(row['#aa in VH CDR2 (ref)'])
										VH_FR3_ref_list.append(row['#aa in VH FR3 (ref)'])
										VH_CDR3_ref_list.append(row['#aa in VH CDR3 (ref)'])
										VH_FR4_ref_list.append(row['#aa in VH FR4 (ref)'])
										VL_FR1_ref_list.append(row['#aa in VL FR1 (ref)'])
										VL_CDR1_ref_list.append(row['#aa in VL CDR1 (ref)'])
										VL_FR2_ref_list.append(row['#aa in VL FR2 (ref)'])
										VL_CDR2_ref_list.append(row['#aa in VL CDR2 (ref)'])
										VL_FR3_ref_list.append(row['#aa in VL FR3 (ref)'])
										VL_CDR3_ref_list.append(row['#aa in VL CDR3 (ref)'])
										VL_FR4_ref_list.append(row['#aa in VL FR4 (ref)'])

	# Writing a dataframe object for outputting. The dataframe is made from the lists that are filled above.
	df_contact = pd.DataFrame(list(zip(
		pdb_type_list, pdb_list, Ab_type_list, Ab_chain_seq_list, Ab_chain_type_list, Ab_chain_id_list, Ab_domain_list,
		Ab_resi_id_list, Ab_resi_aa_list, Ab_atom_list, Ab_atom_coord_list, Ab_SS_list, Ab_resi_ASA_list,
		Ab_raw_extract_list, Ag_chain_id_list, Ag_resi_id_list, Ag_resi_aa_list, Ag_atom_list, Ag_atom_coord_list,
		Ag_SS_list, Ag_resi_ASA_list, Contact_distance_list, resolution_list, H_chain_organism_list,
		L_chain_organism_list, VH_FR1_ref_list, VH_CDR1_ref_list, VH_FR2_ref_list, VH_CDR2_ref_list, VH_FR3_ref_list,
		VH_CDR3_ref_list, VH_FR4_ref_list, VL_FR1_ref_list, VL_CDR1_ref_list, VL_FR2_ref_list, VL_CDR2_ref_list,
		VL_FR3_ref_list, VL_CDR3_ref_list, VL_FR4_ref_list
	)), columns=['Antigen type', 'PDB', 'Ab type', 'Ab chain seq', 'Ab chain type (VH/VL)', 'Ab chain ID', 'Ab domain',
				 'Ab resi ID', 'Ab resi aa', 'Ab atom', 'Ab atom coord', 'Ab secondary structure', 'Ab ASA',
				 'Ab raw extract', 'Ag chain ID', 'Ag resi ID', 'Ag resi aa', 'Ag atom', 'Ag atom coord',
				 'Ag secondary structure', 'Ag ASA', 'Contact distance', 'Structure resolution', 'H chain organism',
				 'L chain organism', '#aa in VH FR1 (ref)', '#aa in VH CDR1 (ref)', '#aa in VH FR2 (ref)',
				 '#aa in VH CDR2 (ref)', '#aa in VH FR3 (ref)', '#aa in VH CDR3 (ref)', '#aa in VH FR4 (ref)',
				 '#aa in VL FR1 (ref)', '#aa in VL CDR1 (ref)', '#aa in VL FR2 (ref)', '#aa in VL CDR2 (ref)',
				 '#aa in VL FR3 (ref)', '#aa in VL CDR3 (ref)', '#aa in VL FR4 (ref)'])

	df_ref = pd.DataFrame(list(zip(
		pdb_type_list_raw, pdb_list_raw, Ab_type_list_raw, Ab_chain_type_list_raw, Ab_chain_id_list_raw,
		Ab_domain_list_raw, Ab_resi_id_list_raw, Ab_resi_aa_list_raw, Ab_atom_list_raw, Ab_atom_coord_list_raw )),
		columns=[
			'Antigen type', 'PDB', 'Ab type', 'Ab chain type (VH/VL)', 'Ab chain ID', 'Ab domain', 'Ab resi ID',
			'Ab resi aa', 'Ab atom', 'Ab atom coord'
		])

	return df_contact, df_ref


def main(summary_data_path: Path, pdb_path: Path, output_path: Path, threads: int, csv_output: bool):
	# Loading in the data from the summary file
	data = pd.read_parquet(summary_data_path)
	print(f"Splitting work-load on {threads} threads")
	results = Pool(threads).starmap(pool_runner, zip(data.iterrows(), repeat(pdb_path)))
	print(f"Outputting results to {output_path}")
	results_contact = pd.concat([x[0] for x in results]).reset_index(drop=True)
	results_contact.to_parquet(output_path / 'Output_contact.parquet')

	results_ref = pd.concat([x[1] for x in results]).reset_index(drop=True)
	results_ref.to_parquet(output_path / 'Output_ref.parquet')

	if csv_output:
		results_contact.to_csv(output_path / 'Output_contact.csv')
		results_ref.to_csv(output_path / 'Output_ref.csv')


if __name__ == '__main__':
	aparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	aparser.add_argument("-s", default=Path('test/work_dir/Summary_all_sorted_nonredundant.parquet'), type=Path, help="Path to Summary_all_sorted_nonredundant.parquet")
	aparser.add_argument("-p", default=Path("test/work_dir/imgt_all_clean_align"), type=Path, help="Path folder with pdb files")
	aparser.add_argument("-o", default=Path('.'), type=Path, help="Path to write output to")
	aparser.add_argument("-t", default=4, type=int, help="Number of threads to use")
	aparser.add_argument("--csv", action='store_true', help="Also output csv format")
	aparser.add_argument("--add_headers", action="store_true", help="Add headers to the pdb files. Needed for dssp version 4.")
	args = aparser.parse_args()
	if args.add_headers:
		from shared_funcs import prepend_headers

		prepend_headers(args.p)
	main(args.s, args.p, args.o, args.t, args.csv)
