class TestScript1:
	def test_output(self, tmp_path):
		from S1_sort_data import main
		from pandas import read_csv
		from pathlib import Path
		from filecmp import cmp
		main(
			summary_data=read_csv('test/test_data/Script1/Summary_all.tsv', sep='\t'),
			pdb_path=Path("test/test_data/Script1/imgt_all_clean/"),
			output_path=Path(tmp_path)
		)

		assert cmp(tmp_path / "fasta.fa", "test/test_data/Script1/fasta.fa")
		assert cmp(tmp_path / "Summary_all_sorted.csv", "test/test_data/Script1/Summary_all_sorted.csv")
