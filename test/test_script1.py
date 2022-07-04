import unittest
from importlib import import_module

class Script1(unittest.TestCase):
	def compare_output(self):
		import_module("1_Data_sort_v2.main")
		self.assertEqual(True, False)  # add assertion here


if __name__ == '__main__':
	unittest.main()
