from filecmp import cmp
from pathlib import Path
from pandas import read_feather, testing


class TestScripts:
	def test_script1(self, tmp_path):
		from S1_sort_data import main
		main(
			summary_data_path=Path('test/work_dir/Summary_all.tsv'),
			pdb_path=Path("test/work_dir/imgt_all_clean/"),
			output_path=Path(tmp_path)
		)

		assert cmp(tmp_path / "fasta.fa", "test/work_dir/fasta.fa")
		test = read_feather(tmp_path / "Summary_all_sorted.fea.zst")
		verified = read_feather("test/work_dir/Summary_all_sorted.fea.zst")
		testing.assert_frame_equal(test, verified)

	def test_script2(self, tmp_path):
		from S2_check_seq_redundancy import main
		main(Path("test/work_dir/Summary_all_sorted.fea.zst"), Path("test/work_dir/Cluster.txt"), Path(tmp_path))
		test = read_feather(tmp_path / "Summary_all_sorted_nonredundant.fea.zst")
		verified = read_feather("test/work_dir/Summary_all_sorted_nonredundant.fea.zst")
		testing.assert_frame_equal(test, verified)
