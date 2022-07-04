from statistics import mean
from pathlib import Path
from pandas import read_csv, DataFrame, concat
from Bio.PDB import PDBParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from shared_funcs import get_ab_type, aa_list, aa_dict


def get_vhvl_seq(structure, v_id: str):
	seq = ''

	# The function extracts the amino acid sequence of the VH or VL domain.
	for Ab_resi_func in structure[0][v_id]:
		# Check that the residue is an amino acid
		if Ab_resi_func.get_resname() in aa_list:

			# Check if the aa belongs to V domain if the residue ID is smaller or equal to 128
			if int(Ab_resi_func.get_id()[1]) <= 128:
				ab_aa_3letter = Ab_resi_func.get_resname()
				ab_aa = aa_dict[ab_aa_3letter]
				seq += ab_aa

	return seq


def create_fasta(summary_data: DataFrame, parser: PDBParser, pdb_path: Path, output_path: Path):
	"""
	Extracting fasta sequence of the antibody V domains. Only amino acids with ID below 128 is included.

	:param summary_data:
	:param parser:
	:param pdb_path:
	:param output_path:
	:return: Nothing
	"""
	sequences = []
	for i, row in summary_data.iterrows():
		print(i, ' --- ', row.pdb)

		# Loading in the pdb file
		pdb_file = pdb_path / f"{row.pdb}.pdb"
		structure = parser.get_structure("structure", pdb_file)

		# Getting the antibody type and the associated chain ids for antibody chains
		vh_id = row.Hchain
		vl_id = row.Lchain
		ab_type = get_ab_type(row[['Hchain', 'Lchain']].isna())

		if ab_type == 'Fv':
			vh_seq = get_vhvl_seq(structure, vh_id)
			vl_seq = get_vhvl_seq(structure, vl_id)
			sequences.append(SeqRecord(Seq(vh_seq), id=f"{row.pdb}_VH", description=""))
			sequences.append(SeqRecord(Seq(vl_seq), id=f"{row.pdb}_VL", description=""))
		elif ab_type == 'VH sdAb':
			vh_seq = get_vhvl_seq(structure, vh_id)
			sequences.append(SeqRecord(Seq(vh_seq), id=f"{row.pdb}_VH", description=""))
		elif ab_type == 'VL sdAb':
			vl_seq = get_vhvl_seq(structure, vl_id)
			sequences.append(SeqRecord(Seq(vl_seq), id=f"{row.pdb}_VL", description=""))

	with open(output_path / "fasta.fa", "w") as output_handle:
		SeqIO.write(sequences, output_handle, "fasta")


def main(summary_data_path: Path, pdb_path: Path, output_path: Path):
	summary_data = read_csv(summary_data_path, sep='\t')
	# Remove all antibodies where antigen is N/A or H chain and L chain are the same i.e. scFvs
	# where H chain and L chain are not properly annotated
	drop_list = []
	for i, row in summary_data.iterrows():
		if row.isna().antigen_chain or (row.Hchain == row.Lchain):
			drop_list.append(i)
	summary_data = summary_data[~summary_data.index.isin(drop_list)]

	# Remove duplicate biological units in the crystal structure.
	# Loop over the individual complexes and only continue with the complex that has the lowest average b factor
	pdb_list_unique = summary_data.pdb.unique()

	pd_list = []
	parser = PDBParser()  # Load parser

	for unique_pdb_name in pdb_list_unique:
		print(unique_pdb_name)
		# Loading the pdb structure into python
		pdb1 = pdb_path / f"{unique_pdb_name}.pdb"
		structure = parser.get_structure("structure", pdb1)

		# Getting summary_data from dataframe relating to the specific pdb
		my_pd = summary_data[summary_data['pdb'] == unique_pdb_name]

		# Check how many rows in the remaining summary_data rows contain the pdb id.
		# Since each row in summary files represents a single biological unit multiple rows indicate the crystal
		# contains more than one biological unit.
		# If only one row is present there is no redundant summary_data in the crystal.
		if len(my_pd) == 1:
			pd_list.append(my_pd.iloc[0])
		else:
			bfactor_avr = []

			for i, row in my_pd:
				ag_chain_list = []
				ab_chain_list = []

				if len(row.antigen_chain) == 1:  # Append the antigen chains to the chain list
					ag_chain_list.append(row.antigen_chain[0])
				elif len(row.antigen_chain) == 5:
					ag_chain_list.append(row.antigen_chain[0])
					ag_chain_list.append(row.antigen_chain[4])
				elif len(row.antigen_chain) == 9:
					ag_chain_list.append(row.antigen_chain[0])
					ag_chain_list.append(row.antigen_chain[4])
					ag_chain_list.append(row.antigen_chain[8])

				ab_type = get_ab_type(row[['Hchain', 'Lchain']].isna())
				if ab_type == 'Fv':
					ab_chain_list.append(row.Hchain)
					ab_chain_list.append(row.Lchain)
				elif ab_type == 'VH sdAb':
					ab_chain_list.append(row.Hchain)
				elif ab_type == 'VL sdAb':
					ab_chain_list.append(row.Lchain)

				bfactor_list = []

				# Go over all atoms in antibody and get bfactors for each atom
				for chain_ID in ab_chain_list:
					for residue_A in structure[0][chain_ID]:
						if residue_A.get_resname() in aa_list:  # Only include if the residue is an amino acid
							for atom_A in residue_A:
								bfactor_list.append(atom_A.get_bfactor())

				# Go over all atoms in antigen and get bfactors for each atom
				for chain_ID in ag_chain_list:
					for residue_B in structure[0][chain_ID]:
						if residue_B.get_resname() in aa_list:  # Only include if the residue is an amino acid
							for atom_B in residue_B:
								bfactor_list.append(atom_B.get_bfactor())

				bfactor_avr.append(mean(bfactor_list))

			bfactor_min_idx = bfactor_avr.index(min(bfactor_avr))
			pd_list.append(my_pd.iloc[bfactor_min_idx])

	summary_data = DataFrame(pd_list)

	# Create for loop for counting amino acids in the different domains
	resi_tot_list = []
	cdr_tot_list = []
	vh_fr1_list = []
	vh_cdr1_list = []
	vh_fr2_list = []
	vh_cdr2_list = []
	vh_fr3_list = []
	vh_cdr3_list = []
	vh_fr4_list = []
	vl_fr1_list = []
	vl_cdr1_list = []
	vl_fr2_list = []
	vl_cdr2_list = []
	vl_fr3_list = []
	vl_cdr3_list = []
	vl_fr4_list = []

	for i, row in summary_data.iterrows():
		resi_tot = 0
		cdr_tot = 0

		vh_fr1 = 0
		vh_cdr1 = 0
		vh_fr2 = 0
		vh_cdr2 = 0
		vh_fr3 = 0
		vh_cdr3 = 0
		vh_fr4 = 0

		vl_fr1 = 0
		vl_cdr1 = 0
		vl_fr2 = 0
		vl_cdr2 = 0
		vl_fr3 = 0
		vl_cdr3 = 0
		vl_fr4 = 0

		pdb3 = pdb_path / f"{row.pdb}.pdb"
		structure = parser.get_structure("structure", pdb3)

		vh_id = row.Hchain
		vl_id = row.Lchain

		# If a VH id is given in the summary file count residues in the different domains
		if not row.isna().Hchain:
			for residue_C in structure[0][vh_id]:
				if residue_C.get_resname() in aa_list:  # Only include if the residue is an amino acid
					if residue_C.get_id()[1] < 27:
						resi_tot += 1
						vh_fr1 += 1
					elif 27 <= residue_C.get_id()[1] <= 38:
						resi_tot += 1
						vh_cdr1 += 1
						cdr_tot += 1
					elif 39 <= residue_C.get_id()[1] <= 55:
						resi_tot += 1
						vh_fr2 += 1
					elif 56 <= residue_C.get_id()[1] <= 65:
						resi_tot += 1
						vh_cdr2 += 1
						cdr_tot += 1
					elif 66 <= residue_C.get_id()[1] <= 104:
						resi_tot += 1
						vh_fr3 += 1
					elif 105 <= residue_C.get_id()[1] <= 117:
						resi_tot += 1
						vh_cdr3 += 1
						cdr_tot += 1
					elif 118 <= residue_C.get_id()[1] <= 128:
						resi_tot += 1
						vh_fr4 += 1

		# If a VL id is given in the summary file count residues in the different domains
		if not row.isna().Lchain:
			for residue_D in structure[0][vl_id]:
				if residue_D.get_resname() in aa_list:  # Only include if the residue is an amino acid
					if residue_D.get_id()[1] < 27:
						resi_tot += 1
						vl_fr1 += 1
					elif 27 <= residue_D.get_id()[1] <= 38:
						resi_tot += 1
						vl_cdr1 += 1
						cdr_tot += 1
					elif 39 <= residue_D.get_id()[1] <= 55:
						resi_tot += 1
						vl_fr2 += 1
					elif 56 <= residue_D.get_id()[1] <= 65:
						resi_tot += 1
						vl_cdr2 += 1
						cdr_tot += 1
					elif 66 <= residue_D.get_id()[1] <= 104:
						resi_tot += 1
						vl_fr3 += 1
					elif 105 <= residue_D.get_id()[1] <= 117:
						resi_tot += 1
						vl_cdr3 += 1
						cdr_tot += 1
					elif 118 <= residue_D.get_id()[1] <= 128:
						resi_tot += 1
						vl_fr4 += 1

		resi_tot_list.append(resi_tot)
		cdr_tot_list.append(cdr_tot)
		vh_fr1_list.append(vh_fr1)
		vh_cdr1_list.append(vh_cdr1)
		vh_fr2_list.append(vh_fr2)
		vh_cdr2_list.append(vh_cdr2)
		vh_fr3_list.append(vh_fr3)
		vh_cdr3_list.append(vh_cdr3)
		vh_fr4_list.append(vh_fr4)
		vl_fr1_list.append(vl_fr1)
		vl_cdr1_list.append(vl_cdr1)
		vl_fr2_list.append(vl_fr2)
		vl_cdr2_list.append(vl_cdr2)
		vl_fr3_list.append(vl_fr3)
		vl_cdr3_list.append(vl_cdr3)
		vl_fr4_list.append(vl_fr4)

	cdr_data = DataFrame({
		'Total residues in Ab V domain': resi_tot_list, 'Total CDR residues': cdr_tot_list,
		'#aa in VH FR1 (ref)': vh_fr1_list, '#aa in VH CDR1 (ref)': vh_cdr1_list, '#aa in VH FR2 (ref)': vh_fr2_list,
		'#aa in VH CDR2 (ref)': vh_cdr2_list, '#aa in VH FR3 (ref)': vh_fr3_list, '#aa in VH CDR3 (ref)': vh_cdr3_list,
		'#aa in VH FR4 (ref)': vh_fr4_list, '#aa in VL FR1 (ref)': vl_fr1_list, '#aa in VL CDR1 (ref)': vl_cdr1_list,
		'#aa in VL FR2 (ref)': vl_fr2_list, '#aa in VL CDR2 (ref)': vl_cdr2_list, '#aa in VL FR3 (ref)': vl_fr3_list,
		'#aa in VL CDR3 (ref)': vl_cdr3_list, '#aa in VL FR4 (ref)': vl_fr4_list}
	)
	data4 = concat([summary_data, cdr_data], axis=1)
	data4.to_feather(output_path / 'Summary_all_sorted.fea.zst', compression='zstd')
	create_fasta(summary_data, parser, pdb_path, output_path)


if __name__ == "__main__":
	main(
		summary_data_path=Path('test/test_data/Script1/Summary_all.tsv'),
		pdb_path=Path("test/test_data/Script1/imgt_all_clean/"),
		output_path=Path("./")
	)
