from pathlib import Path
import pandas as pd

aa_list = ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS', 'GLY', 'PRO', 'ALA', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'TYR', 'VAL']
aa_dict = {'ARG': 'R', 'HIS': 'H', 'LYS': 'K', 'ASP': 'D', 'GLU': 'E', 'SER': 'S', 'THR': 'T', 'ASN': 'N', 'GLN': 'Q', 'CYS': 'C', 'GLY': 'G', 'PRO': 'P', 'ALA': 'A', 'VAL': 'V', 'ILE': 'I', 'LEU': 'L', 'MET': 'M', 'PHE': 'F', 'TYR': 'Y', 'TRP': 'W'}
domain_list = ['VH FR1', 'VH CDR1', 'VH FR2', 'VH CDR2', 'VH FR3', 'VH CDR3', 'VH FR4', 'VL FR1', 'VL CDR1', 'VL FR2', 'VL CDR2', 'VL FR3', 'VL CDR3', 'VL FR4']


def get_ab_type(is_id_na: pd.Series):
	"""
	Determine if both, one or neither of the variable chains have an id
	:param is_id_na: the Hchain and Lchain run through isna
	:return:
	"""
	vh_id = is_id_na.Hchain
	vl_id = is_id_na.Lchain
	if not vh_id and not vl_id:
		return 'Fv'
	elif not vh_id and vl_id:
		return 'VH sdAb'
	elif vh_id and not vl_id:
		return 'VL sdAb'
	else:
		print('Ab type not recognized!')
		return 'N/A'


def prepend_headers(pdb_path: Path):
	if not isinstance(pdb_path, Path):
		pdb_path = Path(pdb_path)
	headers = pd.read_parquet("headers.parquet")

	for pdb_file_path in pdb_path.glob("*.pdb"):
		with open(pdb_file_path, 'r+') as pdb_file:
			c = pdb_file.read()
			if "CRYST" in c:  # check if header already exists
				print(f"Not writing header to {pdb_file_path}. Header already exists.")
			else:
				pdb_file.seek(0)
				pdb_file.write(headers.loc[Path(pdb_file_path).stem, 'header'] + c)
