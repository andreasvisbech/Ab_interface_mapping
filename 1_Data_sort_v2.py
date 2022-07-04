import statistics

import pandas as pd
from Bio.PDB import PDBParser

from shared_funcs import get_ab_type, aa_list, aa_dict


def get_vhvl_seq(structure, v_id):
	seq = ''

	# The function extracts the amino acid sequence of the VH or VL domain.
	for Ab_resi_func in structure[0][v_id]:
		# Check that the residue is an amino acid
		if str(Ab_resi_func.get_resname()) in aa_list:

			# Check if the aa belongs to V domain if the residue ID is smaller or equal to 128
			if int(Ab_resi_func.get_id()[1]) <= 128:
				ab_aa_3letter = str(Ab_resi_func.get_resname())
				ab_aa = aa_dict[ab_aa_3letter]
				seq += ab_aa

	return seq


def main(summary_data, pdb_path):
	parser = PDBParser()  # Load parser
	fasta = ''  # Load empty string for storing fasta in
	pdb_list = summary_data['pdb'].tolist()  # Writing all pdb IDs to a list

	# Small code block to remove all antibodies where antigen is N/A or H chain and L chain are the same i.e. scFvs where H chain and L chain are not properly annotated
	drop_list = []
	for a in range(len(pdb_list)):
		if str(summary_data.iloc[a][4]) == 'nan' or str(summary_data.iloc[a][1]) == str(summary_data.iloc[a][2]):
			drop_list.append(a)

	summary_data = summary_data.drop(labels=drop_list, axis=0)

	# Create for loop with the purpose of removing duplicate biological units in the crystal structure.
	# The loop should go over the individual complexes and only continue with the complex that has lowest average b factor
	pdb_list_unique = summary_data['pdb'].drop_duplicates().tolist()

	pd_list = []

	for b in range(len(pdb_list_unique)):
		print(b)
		# Loading the pdb structure into python
		pdb1 = str(pdb_path + pdb_list_unique[b] + '.pdb')
		structure = parser.get_structure("structure", pdb1)

		# Getting summary_data from dataframe relating to the specific pdb
		my_pd = summary_data[summary_data['pdb'] == pdb_list_unique[b]]

		# Check how many rows in the remaining summary_data rows contain the pdb id.
		# Since each row in summary files represents a single biological unit multiple rows indicate the crystal
		# contains more than one biological unit. If only one row is present there is no redundant summary_data in the crystal.
		if my_pd.shape[0] == 1:
			pd_list.append(my_pd.iloc[0].tolist())

		elif my_pd.shape[0] > 1:
			bfactor_avr = []

			for c in range(my_pd.shape[0]):
				Ag_chain_list = []
				Ab_chain_list = []

				VH_id = str(my_pd.iloc[c][1])
				VL_id = str(my_pd.iloc[c][2])

				if len(str(my_pd.iloc[c][4])) == 1:  # Append the antigen chains to the chain list
					Ag_chain_list.append(str(my_pd.iloc[c][4][0]))
				elif len(str(my_pd.iloc[c][4])) == 5:
					Ag_chain_list.append(str(my_pd.iloc[c][4][0]))
					Ag_chain_list.append(str(my_pd.iloc[c][4][4]))
				elif len(str(my_pd.iloc[c][4])) == 9:
					Ag_chain_list.append(str(my_pd.iloc[c][4][0]))
					Ag_chain_list.append(str(my_pd.iloc[c][4][4]))
					Ag_chain_list.append(str(my_pd.iloc[c][4][8]))

				if VH_id.find('nan') < 0 and VL_id.find('nan') < 0:  # Append the Ab chains to the chain list
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
						if any(str(residue_A.get_resname()) == x for x in aa_list):  # Only include if the residue is an amino acid
							for atom_A in residue_A:
								bfactor_list.append(atom_A.get_bfactor())

				# Go over all atoms in antigen and get bfactors for each atom
				for chain_ID in Ag_chain_list:
					for residue_B in structure[0][chain_ID]:
						if any(str(residue_B.get_resname()) == x for x in aa_list):  # Only include if the residue is an amino acid
							for atom_B in residue_B:
								bfactor_list.append(atom_B.get_bfactor())

				bfactor_avr.append(statistics.mean(bfactor_list))

			bfactor_min_idx = bfactor_avr.index(min(bfactor_avr))
			pd_list.append(my_pd.iloc[bfactor_min_idx].tolist())

	summary_data = pd.DataFrame(data=pd_list, columns=summary_data.columns.values)

	# Creating for loop for extracting fasta sequence of the antibody V domains. Only amino acids with ID below 128 is included.
	for a in range(summary_data.shape[0]):
		print(str(a) + ' --- ' + str(summary_data.iloc[a][0]))
		#    VH_id = str(summary_data.iloc[a][1])
		#    VL_id = str(summary_data.iloc[a][2])

		# Loading in the pdb file
		pdb2 = str(pdb_path + str(summary_data.iloc[a][0]) + '.pdb')
		structure = parser.get_structure("structure", pdb2)

		# Getting the antibody type and the associated chain ids for antibody chains
		VH_id = str(summary_data.iloc[a][1])
		VL_id = str(summary_data.iloc[a][2])
		Ab_type = get_ab_type(VH_id, VL_id)

		if Ab_type == 'Fv':
			VH_seq = get_vhvl_seq(structure, VH_id)
			VL_seq = get_vhvl_seq(structure, VL_id)

			fasta = fasta + '>' + str(summary_data.iloc[a][0]) + ', VH' + '\n' + VH_seq + '\n'
			fasta = fasta + '>' + str(summary_data.iloc[a][0]) + ', VL' + '\n' + VL_seq + '\n'

		elif Ab_type == 'VH sdAb':
			VH_seq = get_vhvl_seq(structure, VH_id)

			fasta = fasta + '>' + str(summary_data.iloc[a][0]) + ', VH' + '\n' + VH_seq + '\n'

		elif Ab_type == 'VL sdAb':
			VL_seq = get_vhvl_seq(structure, VL_id)

			fasta = fasta + '>' + str(summary_data.iloc[a][0]) + ', VH' + '\n' + VL_seq + '\n'

	with open("fasta.fa", "a") as file2_out:
		file2_out.write(fasta)

	# Create for loop for counting amino acids in the different domains
	resi_tot_list = []
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

	for e in range(summary_data.shape[0]):
		resi_tot = 0
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

		pdb3 = str(pdb_path + str(summary_data['pdb'][e]) + '.pdb')
		structure = parser.get_structure("structure", pdb3)

		Ab_chain_list = []

		VH_id = str(summary_data.iloc[e][1])
		VL_id = str(summary_data.iloc[e][2])

		# If a VH id is given in the summary file count residues in the different domains
		if VH_id != 'nan':
			for residue_C in structure[0][VH_id]:
				if any(str(residue_C.get_resname()) == x for x in
					   aa_list):  # Only include if the residue is an amino acid

					if residue_C.get_id()[1] < 27:
						resi_tot += 1
						VH_FR1 += 1
					elif 27 <= residue_C.get_id()[1] <= 38:
						resi_tot += 1
						VH_CDR1 += 1
						CDR_tot += 1
					elif 39 <= residue_C.get_id()[1] <= 55:
						resi_tot += 1
						VH_FR2 += 1
					elif 56 <= residue_C.get_id()[1] <= 65:
						resi_tot += 1
						VH_CDR2 += 1
						CDR_tot += 1
					elif 66 <= residue_C.get_id()[1] <= 104:
						resi_tot += 1
						VH_FR3 += 1
					elif 105 <= residue_C.get_id()[1] <= 117:
						resi_tot += 1
						VH_CDR3 += 1
						CDR_tot += 1
					elif 118 <= residue_C.get_id()[1] <= 128:
						resi_tot += 1
						VH_FR4 += 1

		# If a VL id is given in the summary file count residues in the different domains
		if VL_id != 'nan':
			for residue_D in structure[0][VL_id]:
				if any(str(residue_D.get_resname()) == x for x in
					   aa_list):  # Only include if the residue is an amino acid
					if residue_D.get_id()[1] < 27:
						resi_tot += 1
						VL_FR1 += 1
					elif 27 <= residue_D.get_id()[1] <= 38:
						resi_tot += 1
						VL_CDR1 += 1
						CDR_tot += 1
					elif 39 <= residue_D.get_id()[1] <= 55:
						resi_tot += 1
						VL_FR2 += 1
					elif 56 <= residue_D.get_id()[1] <= 65:
						resi_tot += 1
						VL_CDR2 += 1
						CDR_tot += 1
					elif 66 <= residue_D.get_id()[1] <= 104:
						resi_tot += 1
						VL_FR3 += 1
					elif 105 <= residue_D.get_id()[1] <= 117:
						resi_tot += 1
						VL_CDR3 += 1
						CDR_tot += 1
					elif 118 <= residue_D.get_id()[1] <= 128:
						resi_tot += 1
						VL_FR4 += 1

		resi_tot_list.append(resi_tot)
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

	CDR_data = pd.DataFrame({
		'Total residues in Ab V domain': resi_tot_list,
		'Total CDR residues': CDR_tot_list,
		'#aa in VH FR1 (ref)': VH_FR1_list,
		'#aa in VH CDR1 (ref)': VH_CDR1_list,
		'#aa in VH FR2 (ref)': VH_FR2_list,
		'#aa in VH CDR2 (ref)': VH_CDR2_list,
		'#aa in VH FR3 (ref)': VH_FR3_list,
		'#aa in VH CDR3 (ref)': VH_CDR3_list,
		'#aa in VH FR4 (ref)': VH_FR4_list,
		'#aa in VL FR1 (ref)': VL_FR1_list,
		'#aa in VL CDR1 (ref)': VL_CDR1_list,
		'#aa in VL FR2 (ref)': VL_FR2_list,
		'#aa in VL CDR2 (ref)': VL_CDR2_list,
		'#aa in VL FR3 (ref)': VL_FR3_list,
		'#aa in VL CDR3 (ref)': VL_CDR3_list,
		'#aa in VL FR4 (ref)': VL_FR4_list}
	)
	data4 = pd.concat([summary_data, CDR_data], axis=1)
	data4.to_csv('Summary_all_sorted.csv', sep=';', index=False)


if __name__ == "__main__":
	main(pd.read_csv('test/test_data/Script1/Summary_all.tsv', sep='\t', header=0), "test/test_data/Script1/imgt_all_clean/")
