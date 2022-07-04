aa_list = ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS', 'GLY', 'PRO', 'ALA', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'TYR', 'VAL']
aa_dict = {'ARG': 'R', 'HIS': 'H', 'LYS': 'K', 'ASP': 'D', 'GLU': 'E', 'SER': 'S', 'THR': 'T', 'ASN': 'N', 'GLN': 'Q', 'CYS': 'C', 'GLY': 'G', 'PRO': 'P', 'ALA': 'A', 'VAL': 'V', 'ILE': 'I', 'LEU': 'L', 'MET': 'M', 'PHE': 'F', 'TYR': 'Y', 'TRP': 'W'}


def get_ab_type(vh_id, vl_id):
	if vh_id.find('nan') < 0 and vl_id.find('nan') < 0:
		func_var = 'Fv'

	elif vh_id.find('nan') < 0 and vl_id.find('nan') >= 0:
		func_var = 'VH sdAb'

	elif vh_id.find('nan') >= 0 and vl_id.find('nan') < 0:
		func_var = 'VL sdAb'
	else:
		func_var = 'N/A'

	return func_var
