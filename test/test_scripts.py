from filecmp import cmp
from pathlib import Path


class TestScripts:
	def test_script1(self, tmp_path):
		from S1_sort_data import main
		main(
			summary_data_path=Path('test/work_dir/Summary_all.tsv'),
			pdb_path=Path("test/work_dir/imgt_all_clean/"),
			output_path=Path(tmp_path)
		)

		assert cmp(tmp_path / "fasta.fa", "test/work_dir/fasta.fa")
		assert cmp(tmp_path / "Summary_all_sorted.fea.zst", "test/work_dir/Summary_all_sorted.fea.zst")

	def test_script2(self, tmp_path):
		from S2_check_seq_redundancy import main
		main(Path("test/work_dir/Summary_all_sorted.fea.zst"), Path("test/work_dir/Cluster.txt"), Path(tmp_path))
		assert cmp(tmp_path / "Summary_all_sorted_nonredundant.fea.zst", "test/work_dir/Summary_all_sorted_nonredundant.fea.zst")
