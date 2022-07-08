from filecmp import cmp
from pathlib import Path
from pandas import read_parquet, testing


class TestScripts:
	def test_script1(self, tmp_path):
		from S1_sort_data import main
		main(
			summary_data_path=Path('test/work_dir/Summary_all.tsv'),
			pdb_path=Path("test/work_dir/imgt_all_clean/"),
			output_path=Path(tmp_path),
			csv_output=False
		)

		assert cmp(tmp_path / "fasta.fa", "test/work_dir/fasta.fa")
		test = read_parquet(tmp_path / "Summary_all_sorted.parquet")
		verified = read_parquet("test/work_dir/Summary_all_sorted.parquet")
		testing.assert_frame_equal(test, verified)

	def test_script2(self, tmp_path):
		from S2_check_seq_redundancy import main
		main(Path("test/work_dir/Summary_all_sorted.parquet"), Path("test/work_dir/Cluster.txt"), Path(tmp_path), False)
		test = read_parquet(tmp_path / "Summary_all_sorted_nonredundant.parquet")
		verified = read_parquet("test/work_dir/Summary_all_sorted_nonredundant.parquet")
		testing.assert_frame_equal(test, verified)

	def test_script3(self, tmp_path):
		from S3_extract_protein_multi import main
		main(
			summary_data_path=Path('test/work_dir/Summary_all_sorted_nonredundant.parquet'),
			pdb_path=Path("test/work_dir/imgt_all_clean_align"),
			output_path=Path(tmp_path),
			threads=4,
			csv_output=False
		)
		test_contact = read_parquet(tmp_path / "Output_contact.parquet")
		verified_contact = read_parquet("test/work_dir/Output_contact.parquet")
		testing.assert_frame_equal(test_contact, verified_contact)
		test_ref = read_parquet(tmp_path / "Output_ref.parquet")
		verified_ref = read_parquet("test/work_dir/Output_ref.parquet")
		testing.assert_frame_equal(test_ref, verified_ref)
