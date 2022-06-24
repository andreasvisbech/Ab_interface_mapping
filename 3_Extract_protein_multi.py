# Import modules 
import pandas as pd
import numpy as np
from datetime import datetime
from Bio.PDB import *
from Bio import SeqIO
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
from Bio.PDB.DSSP import DSSP

from multiprocessing import Pool

# Set number of multiprocesses
mul = 16

parser = PDBParser()
ppb = PPBuilder()


# Loading functions
def get_Ag_ids(line):
	# The function looks at the length of the input and gets the chain ids.

	data = str(line['antigen_chain'])

	func_list = []

	# If crystal contains 1 antigen chain the length will be 1. If two antigen chains the length will be 5. If 3 antigen chain length will be 4

	if len(data) == 1:
		func_list.append(data[0])
	elif len(data) == 5:
		func_list.append(data[0])
		func_list.append(data[4])
	elif len(data) == 9:
		func_list.append(data[0])
		func_list.append(data[4])
		func_list.append(data[8])
	elif len(data) == 13:
		func_list.append(data[0])
		func_list.append(data[4])
		func_list.append(data[8])
		func_list.append(data[12])

	return func_list


def get_Ab_type(VH_id, VL_id):
	func_dict = {}

	if VH_id.find('nan') < 0 and VL_id.find('nan') < 0:
		func_var = 'Fv'
		# func_dict['Ab_ids'] = [VH_id , VL_id]

	elif VH_id.find('nan') < 0 and VL_id.find('nan') >= 0:
		func_var = 'VH sdAb'
		# func_dict['Ab_ids'] = [VH_id]

	elif VH_id.find('nan') >= 0 and VL_id.find('nan') < 0:
		func_var = 'VL sdAb'
		# func_dict['Ab_ids'] = [VL_id]

	else:
		func_var = 'N/A'

	return func_var


def get_chain_seq(structure, chain_id):
	# The function gets the full sequence of the antibody chain in question.

	func_var = ''
	try:
		for pp in ppb.build_peptides(structure[0][chain_id]):
			func_var = func_var + str(pp.get_sequence())
	except KeyError:
		func_var = 'N/A'
	return func_var


def get_Ab_domain(resi):
	# The function determines which antibody domain the residue belongs to (FR1, CDR1, FR2, CDR2, FR3, CDR3 or FR4).
	# The numbering comes from https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html

	if 27 > resi.get_id()[1]:
		CDR_value = 'FR1'
	elif 27 <= resi.get_id()[1] and 38 >= resi.get_id()[1]:
		CDR_value = 'CDR1'
	elif 38 < resi.get_id()[1] and 56 > resi.get_id()[1]:
		CDR_value = 'FR2'
	elif 56 <= resi.get_id()[1] and 65 >= resi.get_id()[1]:
		CDR_value = 'CDR2'
	elif 65 < resi.get_id()[1] and 105 > resi.get_id()[1]:
		CDR_value = 'FR3'
	elif 105 <= resi.get_id()[1] and 117 >= resi.get_id()[1]:
		CDR_value = 'CDR3'
	elif 117 < resi.get_id()[1] <= 128:
		CDR_value = 'FR4'
	else:
		CDR_value = 'N/A'

	return CDR_value


def domain_count(structure, Ab_type, VH_id, VL_id, domain_list):
	out_dict = {}

	VH_FR1_ref_count = 0
	VH_CDR1_ref_count = 0
	VH_FR2_ref_count = 0
	VH_CDR2_ref_count = 0
	VH_FR3_ref_count = 0
	VH_CDR3_ref_count = 0
	VH_FR4_ref_count = 0

	VL_FR1_ref_count = 0
	VL_CDR1_ref_count = 0
	VL_FR2_ref_count = 0
	VL_CDR2_ref_count = 0
	VL_FR3_ref_count = 0
	VL_CDR3_ref_count = 0
	VL_FR4_ref_count = 0

	if Ab_type == 'Fv':

		for Ab_resi in structure[0][VH_id]:

			CDR_value = get_Ab_domain(Ab_resi)
			if CDR_value == 'FR1':
				VH_FR1_ref_count = VH_FR1_ref_count + 1
			elif CDR_value == 'CDR1':
				VH_CDR1_ref_count = VH_CDR1_ref_count + 1
			elif CDR_value == 'FR2':
				VH_FR2_ref_count = VH_FR2_ref_count + 1
			elif CDR_value == 'CDR2':
				VH_CDR2_ref_count = VH_CDR2_ref_count + 1
			elif CDR_value == 'FR3':
				VH_FR3_ref_count = VH_FR3_ref_count + 1
			elif CDR_value == 'CDR3':
				VH_CDR3_ref_count = VH_CDR3_ref_count + 1
			elif CDR_value == 'FR4':
				VH_FR4_ref_count = VH_FR4_ref_count + 1

		for Ab_resi in structure[0][VL_id]:

			CDR_value = get_Ab_domain(Ab_resi)
			if CDR_value == 'FR1':
				VL_FR1_ref_count = VL_FR1_ref_count + 1
			elif CDR_value == 'CDR1':
				VL_CDR1_ref_count = VL_CDR1_ref_count + 1
			elif CDR_value == 'FR2':
				VL_FR2_ref_count = VL_FR2_ref_count + 1
			elif CDR_value == 'CDR2':
				VL_CDR2_ref_count = VL_CDR2_ref_count + 1
			elif CDR_value == 'FR3':
				VL_FR3_ref_count = VL_FR3_ref_count + 1
			elif CDR_value == 'CDR3':
				VL_CDR3_ref_count = VL_CDR3_ref_count + 1
			elif CDR_value == 'FR4':
				VL_FR4_ref_count = VL_FR4_ref_count + 1

	elif Ab_type == 'VH sdAb':

		for Ab_resi in structure[0][VH_id]:

			CDR_value = get_Ab_domain(Ab_resi)
			if CDR_value == 'FR1':
				VH_FR1_ref_count = VH_FR1_ref_count + 1
			elif CDR_value == 'CDR1':
				VH_CDR1_ref_count = VH_CDR1_ref_count + 1
			elif CDR_value == 'FR2':
				VH_FR2_ref_count = VH_FR2_ref_count + 1
			elif CDR_value == 'CDR2':
				VH_CDR2_ref_count = VH_CDR2_ref_count + 1
			elif CDR_value == 'FR3':
				VH_FR3_ref_count = VH_FR3_ref_count + 1
			elif CDR_value == 'CDR3':
				VH_CDR3_ref_count = VH_CDR3_ref_count + 1
			elif CDR_value == 'FR4':
				VH_FR4_ref_count = VH_FR4_ref_count + 1

	elif Ab_type == 'VL sdAb':

		for Ab_resi in structure[0][VL_id]:

			CDR_value = get_Ab_domain(Ab_resi)
			if CDR_value == 'FR1':
				VL_FR1_ref_count = VL_FR1_ref_count + 1
			elif CDR_value == 'CDR1':
				VL_CDR1_ref_count = VL_CDR1_ref_count + 1
			elif CDR_value == 'FR2':
				VL_FR2_ref_count = VL_FR2_ref_count + 1
			elif CDR_value == 'CDR2':
				VL_CDR2_ref_count = VL_CDR2_ref_count + 1
			elif CDR_value == 'FR3':
				VL_FR3_ref_count = VL_FR3_ref_count + 1
			elif CDR_value == 'CDR3':
				VL_CDR3_ref_count = VL_CDR3_ref_count + 1
			elif CDR_value == 'FR4':
				VL_FR4_ref_count = VL_FR4_ref_count + 1

	out_dict['VH_FR1_ref_count'] = VH_FR1_ref_count
	out_dict['VH_CDR1_ref_count'] = VH_CDR1_ref_count
	out_dict['VH_FR2_ref_count'] = VH_FR2_ref_count
	out_dict['VH_CDR2_ref_count'] = VH_CDR2_ref_count
	out_dict['VH_FR3_ref_count'] = VH_FR3_ref_count
	out_dict['VH_CDR3_ref_count'] = VH_CDR3_ref_count
	out_dict['VH_FR4_ref_count'] = VH_FR4_ref_count

	out_dict['VL_FR1_ref_count'] = VL_FR1_ref_count
	out_dict['VL_CDR1_ref_count'] = VL_CDR1_ref_count
	out_dict['VL_FR2_ref_count'] = VL_FR2_ref_count
	out_dict['VL_CDR2_ref_count'] = VL_CDR2_ref_count
	out_dict['VL_FR3_ref_count'] = VL_FR3_ref_count
	out_dict['VL_CDR3_ref_count'] = VL_CDR3_ref_count
	out_dict['VL_FR4_ref_count'] = VL_FR4_ref_count

	return out_dict


def get_resi_id(resi):
	# The function determines the residue id based on the pdb structure supplied.
	# The function includes both the number (index 1) and a potential letter (index 2) which can be present in the pdb.
	# Lastly the function strips any whitespace.

	func_var = str(resi.get_id()[1]) + str(resi.get_id()[2])
	func_var = str(func_var).strip()

	return func_var


def get_SS(dssp_dict, chain_id, resi):
	func_var = dssp_dict[chain_id, (' ', resi.get_id()[1], resi.get_id()[2])][2]

	return func_var


def get_atom_coord(atom):
	func_list = []

	func_list.append(atom.get_coord()[0])
	func_list.append(atom.get_coord()[1])
	func_list.append(atom.get_coord()[2])

	return func_list


def get_resi_ASA(dssp_dict, chain_id, resi):
	# Function calculates available surface area of the residue using the dssp

	func_var = dssp_dict[chain_id, (' ', resi.get_id()[1], resi.get_id()[2])][3]

	return func_var


def Ab_raw_extract(structure, index, aa_list):
	out_list = []

	try:
		for func_resi in structure[0][index]:

			# Check that the residue is an amino acid
			if any(str(func_resi.get_resname()) == x for x in aa_list):
				func_var = func_resi.get_resname()

				func_var2 = str(func_resi.get_id()[1]) + str(func_resi.get_id()[2])
				func_var2 = str(func_var2).strip()

				func_var3 = func_var + func_var2

				out_list.append(func_var3)

	except KeyError:
		out_list.append('N/A')

	return out_list


# def writing_raw(pdb_id, Ab_type, Ab_chain_type, Ab_domain, Ab_resi_id, Ab_resi_aa, Ab_atom_id, Ab_atom_coord):
def writing_raw(a, pdb_type, pdb_id, Ab_type, Ab_chain_type, Ab_idx, Ab_domain, Ab_resi_aa, Ab_resi_id, Ab_atom_id,
				Ab_atom_coord):
	if any(str(Ab_resi_aa) == x for x in aa_list) and Ab_atom_id[0] != 'H':
		master_dict_raw['pdb_type_list_raw'].append(pdb_type)
		master_dict_raw['pdb_raw'].append(pdb_id)
		master_dict_raw['Ab_type_list_raw'].append(Ab_type)
		master_dict_raw['Ab_chain_type_list_raw'].append(Ab_chain_type)
		master_dict_raw['Ab_chain_id_list_raw'].append(Ab_idx)
		master_dict_raw['Ab_domain_list_raw'].append(Ab_domain)
		master_dict_raw['Ab_resi_aa_list_raw'].append(Ab_resi_aa)
		master_dict_raw['Ab_resi_id_list_raw'].append(Ab_resi_id)
		master_dict_raw['Ab_atom_list_raw'].append(Ab_atom_id)
		master_dict_raw['Ab_atom_coord_list_raw'].append(Ab_atom_coord)

		master_dict_raw['#aa in VH FR1 (reference data)'] = VH_FR1_ref_list[a]
		master_dict_raw['#aa in VH CDR1 (reference data)'] = VH_CDR1_ref_list[a]
		master_dict_raw['#aa in VH FR2 (reference data)'] = VH_FR2_ref_list[a]
		master_dict_raw['#aa in VH CDR2 (reference data)'] = VH_CDR2_ref_list[a]
		master_dict_raw['#aa in VH FR3 (reference data)'] = VH_FR3_ref_list[a]
		master_dict_raw['#aa in VH CDR3 (reference data)'] = VH_CDR3_ref_list[a]
		master_dict_raw['#aa in VH FR4 (reference data)'] = VH_FR4_ref_list[a]
		master_dict_raw['#aa in VL FR1 (reference data)'] = VL_FR1_ref_list[a]
		master_dict_raw['#aa in VL CDR1 (reference data)'] = VL_CDR1_ref_list[a]
		master_dict_raw['#aa in VL FR2 (reference data)'] = VL_FR2_ref_list[a]
		master_dict_raw['#aa in VL CDR2 (reference data)'] = VL_CDR2_ref_list[a]
		master_dict_raw['#aa in VL FR3 (reference data)'] = VL_FR3_ref_list[a]
		master_dict_raw['#aa in VL CDR3 (reference data)'] = VL_CDR3_ref_list[a]
		master_dict_raw['#aa in VL FR4 (reference data)'] = VL_FR4_ref_list[a]


def set_parameters():
	# Distance cutoff
	contact_cutoff = 5

	# Setting path for the pdb files
	path = '/work1/avima/Ab_interface/imgt_all_clean_align/'

	# Initialize variables for run
	aa_list = aa_list = ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS', 'GLY', 'PRO', 'ALA',
						 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'TYR', 'VAL']

	domain_list = ['VH FR1', 'VH CDR1', 'VH FR2', 'VH CDR2', 'VH FR3', 'VH CDR3', 'VH FR4', 'VL FR1', 'VL CDR1',
				   'VL FR2', 'VL CDR2', 'VL FR3', 'VL CDR3', 'VL FR4']

	return contact_cutoff, path, aa_list, domain_list


def main_func(line):
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

	# Setting path for the pdb files
	contact_cutoff, path, aa_list, domain_list = set_parameters()

	# Getting object from the line
	line = line[1]

	now = datetime.now()
	start = now.strftime("%H:%M:%S")

	print('Running: ' + str(line['pdb']) + ' (' + str(line['my_index']) + ') --- start: ' + str(start))

	# Loading in the pdb file
	pdb1 = str(path + line['pdb'] + '.pdb')
	structure = parser.get_structure("structure", pdb1)

	# Getting dssp dict
	dssp_dict = DSSP(structure[0], pdb1, dssp='/work1/avima/Ab_interface/Data_extract/dssp-2.0.4-linux-amd64')

	# Get antigen type
	pdb_type = line['Antigen group']

	# Get structure resolution
	resolution = line['resolution']

	# Get structure H chain organism
	VH_organism = line['heavy_species']

	# Get structure L chain organism
	VL_organism = line['light_species']

	# Save pdb id
	pdb_id = line['pdb']

	# Get a list of antigen ids in the structure
	Ag_id_list = get_Ag_ids(line)

	# Getting the antibody type and the associated chain ids for antibody chains
	VH_id = str(line['Hchain'])
	VL_id = str(line['Lchain'])
	Ab_type = get_Ab_type(VH_id, VL_id)
	if Ab_type == 'Fv':
		Ab_id_list = [VH_id, VL_id]
	elif Ab_type == 'VH sdAb':
		Ab_id_list = [VH_id]
	elif Ab_type == 'VL sdAb':
		Ab_id_list = [VL_id]
	else:
		Ab_id_list = []

	# Create for loop that iterates through all antibody chains
	for Ab_idx in Ab_id_list:

		# Get sequence of Ab chain and put into chain id list
		Ab_chain_seq = get_chain_seq(structure, Ab_idx)

		# Get raw extract of the antibody chain
		Ab_raw = Ab_raw_extract(structure, Ab_idx, aa_list)

		# Get the antibody chain type i.e. if it is heavy chain or light chain
		if Ab_idx == VH_id:
			Ab_chain_type = 'VH'
		elif Ab_idx == VL_id:
			Ab_chain_type = 'VL'
		else:
			Ab_chain_type = 'N/A'

		# Create for loop to go over all residues from the Ab_idx chain
		for Ab_resi in structure[0][Ab_idx]:

			# Check that the residue is an amino acid
			if any(str(Ab_resi.get_resname()) == x for x in aa_list):

				# Get the Ab domain for the residue
				Ab_domain = get_Ab_domain(Ab_resi)

				# Get the id of the antibody residue
				Ab_resi_id = get_resi_id(Ab_resi)

				# Get the aa of the specific Ab residue
				Ab_resi_aa = Ab_resi.get_resname()

				# Get residue secondary structure
				try:
					Ab_resi_SS = get_SS(dssp_dict, Ab_idx, Ab_resi)
				except KeyError:
					Ab_resi_SS = 'Missing atoms in residue'
					continue

				# Get antibody available surface area for the residue
				Ab_resi_ASA = get_resi_ASA(dssp_dict, Ab_idx, Ab_resi)

				# Create for loop that goes through all atoms in the Ab residue
				for Ab_atom in Ab_resi:

					# Get Ab atom
					Ab_atom_id = Ab_atom.get_id()

					# Get Ab atom coords
					Ab_atom_coord = get_atom_coord(Ab_atom)

					# Save the raw reference data irrespective of whether it is contact or not
					pdb_type_list_raw.append(pdb_type)
					pdb_list_raw.append(pdb_id)
					Ab_type_list_raw.append(Ab_type)
					Ab_chain_type_list_raw.append(Ab_chain_type)
					Ab_chain_id_list_raw.append(Ab_idx)
					Ab_domain_list_raw.append(Ab_domain)
					Ab_resi_id_list_raw.append(Ab_resi_id)
					Ab_resi_aa_list_raw.append(Ab_resi_aa)
					Ab_atom_list_raw.append(Ab_atom_id)
					Ab_atom_coord_list_raw.append(Ab_atom_coord)

					# Create for loops that go through all atoms in antigen
					for Ag_idx in Ag_id_list:

						for Ag_resi in structure[0][Ag_idx]:

							if any(str(Ag_resi.get_resname()) == x for x in aa_list):  ### Only consider amino acids

								# Get the aa of the specific Ab residue
								Ag_resi_aa = Ag_resi.get_resname()

								# Get the id of the antibody residue
								Ag_resi_id = get_resi_id(Ag_resi)

								# Get residue secondary structure
								try:
									Ag_resi_SS = get_SS(dssp_dict, Ag_idx, Ag_resi)
								except KeyError:
									Ag_resi_SS = 'Missing atoms in residue'
									continue

								# Get antibody available surface area for the residue
								Ag_resi_ASA = get_resi_ASA(dssp_dict, Ag_idx, Ag_resi)

								for Ag_atom in Ag_resi:

									# Get Ag atom
									Ag_atom_id = Ag_atom.get_id()

									# Get Ab atom coords
									Ag_atom_coord = get_atom_coord(Ag_atom)

									# Determine distance between antibody atom and antigen atom and see if the distance is below user specified cutoff i.e. if the atom-atom pair is considered a contact
									atom_distance = np.linalg.norm(Ab_atom.coord - Ag_atom.coord)

									if atom_distance <= contact_cutoff and Ab_atom_id[0] != 'H' and Ag_atom_id[
										0] != 'H':
										pdb_type_list.append(pdb_type)
										pdb_list.append(pdb_id)
										Ab_type_list.append(Ab_type)
										Ab_chain_seq_list.append(Ab_chain_seq)
										Ab_chain_type_list.append(Ab_chain_type)
										Ab_chain_id_list.append(Ab_idx)
										Ab_domain_list.append(Ab_domain)
										Ab_resi_id_list.append(Ab_resi_id)
										Ab_resi_aa_list.append(Ab_resi_aa)
										Ab_atom_list.append(Ab_atom_id)
										Ab_atom_coord_list.append(Ab_atom_coord)
										Ab_SS_list.append(Ab_resi_SS)
										Ab_resi_ASA_list.append(Ab_resi_ASA)
										Ab_raw_extract_list.append(Ab_raw)

										Ag_chain_id_list.append(Ag_idx)
										Ag_resi_id_list.append(Ag_resi_id)
										Ag_resi_aa_list.append(Ag_resi_aa)
										Ag_atom_list.append(Ag_atom_id)
										Ag_atom_coord_list.append(Ag_atom_coord)
										Ag_SS_list.append(Ag_resi_SS)
										Ag_resi_ASA_list.append(Ag_resi_ASA)

										Contact_distance_list.append(atom_distance)
										resolution_list.append(resolution)
										H_chain_organism_list.append(VH_organism)
										L_chain_organism_list.append(VL_organism)

										VH_FR1_ref_list.append(line['#aa in VH FR1 (ref)'])
										VH_CDR1_ref_list.append(line['#aa in VH CDR1 (ref)'])
										VH_FR2_ref_list.append(line['#aa in VH FR2 (ref)'])
										VH_CDR2_ref_list.append(line['#aa in VH CDR2 (ref)'])
										VH_FR3_ref_list.append(line['#aa in VH FR3 (ref)'])
										VH_CDR3_ref_list.append(line['#aa in VH CDR3 (ref)'])
										VH_FR4_ref_list.append(line['#aa in VH FR4 (ref)'])
										VL_FR1_ref_list.append(line['#aa in VL FR1 (ref)'])
										VL_CDR1_ref_list.append(line['#aa in VL CDR1 (ref)'])
										VL_FR2_ref_list.append(line['#aa in VL FR2 (ref)'])
										VL_CDR2_ref_list.append(line['#aa in VL CDR2 (ref)'])
										VL_FR3_ref_list.append(line['#aa in VL FR3 (ref)'])
										VL_CDR3_ref_list.append(line['#aa in VL CDR3 (ref)'])
										VL_FR4_ref_list.append(line['#aa in VL FR4 (ref)'])

	# Writing a dataframe object for outputting. The dataframe is made from the lists that are filled above.
	df_contact = pd.DataFrame(list(zip(pdb_type_list,
									   pdb_list,
									   Ab_type_list,
									   Ab_chain_seq_list,
									   Ab_chain_type_list,
									   Ab_chain_id_list,
									   Ab_domain_list,
									   Ab_resi_id_list,
									   Ab_resi_aa_list,
									   Ab_atom_list,
									   Ab_atom_coord_list,
									   Ab_SS_list,
									   Ab_resi_ASA_list,
									   Ab_raw_extract_list,
									   Ag_chain_id_list,
									   Ag_resi_id_list,
									   Ag_resi_aa_list,
									   Ag_atom_list,
									   Ag_atom_coord_list,
									   Ag_SS_list,
									   Ag_resi_ASA_list,
									   Contact_distance_list,
									   resolution_list,
									   H_chain_organism_list,
									   L_chain_organism_list,
									   VH_FR1_ref_list,
									   VH_CDR1_ref_list,
									   VH_FR2_ref_list,
									   VH_CDR2_ref_list,
									   VH_FR3_ref_list,
									   VH_CDR3_ref_list,
									   VH_FR4_ref_list,
									   VL_FR1_ref_list,
									   VL_CDR1_ref_list,
									   VL_FR2_ref_list,
									   VL_CDR2_ref_list,
									   VL_FR3_ref_list,
									   VL_CDR3_ref_list,
									   VL_FR4_ref_list
									   )),
							  columns=['Antigen type',
									   'PDB',
									   'Ab type',
									   'Ab chain seq',
									   'Ab chain type (VH/VL)',
									   'Ab chain ID',
									   'Ab domain',
									   'Ab resi ID',
									   'Ab resi aa',
									   'Ab atom',
									   'Ab atom coord',
									   'Ab secondary structure',
									   'Ab ASA',
									   'Ab raw extract',
									   'Ag chain ID',
									   'Ag resi ID',
									   'Ag resi aa',
									   'Ag atom',
									   'Ag atom coord',
									   'Ag secondary structure',
									   'Ag ASA',
									   'Contact distance',
									   'Structure resolution',
									   'H chain organism',
									   'L chain organism',
									   '#aa in VH FR1 (ref)',
									   '#aa in VH CDR1 (ref)',
									   '#aa in VH FR2 (ref)',
									   '#aa in VH CDR2 (ref)',
									   '#aa in VH FR3 (ref)',
									   '#aa in VH CDR3 (ref)',
									   '#aa in VH FR4 (ref)',
									   '#aa in VL FR1 (ref)',
									   '#aa in VL CDR1 (ref)',
									   '#aa in VL FR2 (ref)',
									   '#aa in VL CDR2 (ref)',
									   '#aa in VL FR3 (ref)',
									   '#aa in VL CDR3 (ref)',
									   '#aa in VL FR4 (ref)'
									   ])

	df_ref = pd.DataFrame(list(zip(pdb_type_list_raw,
								   pdb_list_raw,
								   Ab_type_list_raw,
								   Ab_chain_type_list_raw,
								   Ab_chain_id_list_raw,
								   Ab_domain_list_raw,
								   Ab_resi_id_list_raw,
								   Ab_resi_aa_list_raw,
								   Ab_atom_list_raw,
								   Ab_atom_coord_list_raw
								   )),
						  columns=['Antigen type',
								   'PDB',
								   'Ab type',
								   'Ab chain type (VH/VL)',
								   'Ab chain ID',
								   'Ab domain',
								   'Ab resi ID',
								   'Ab resi aa',
								   'Ab atom',
								   'Ab atom coord'
								   ])

	#return df_contact
	return df_ref


# Loading in the data from the summary file
data = pd.read_csv('/work1/avima/Ab_interface/Data_extract/Summary_all_sorted_nonredundant_index.csv', sep=';', header=0)

if __name__ == '__main__':
	pool = Pool(mul)
	results = pool.map(main_func, data.iterrows())

	#results_contact = pd.concat(results)
	#results_contact.to_csv('Output_contact.csv', sep=';')

	results_ref = pd.concat(results)
	results_ref.to_csv('Output_ref.csv', sep=';')




