'''
Perform unit tests for essential based methods

Author: Kewalin Samart
Created date: 06/09/2023
Last Modified: 06/11/2023
'''
import os, sys, io
cur_path = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, cur_path+"/..")

import unittest
import unittest.mock
from src.utilities import *

class Kmer_tests(unittest.TestCase):
    def test_getkmers1(self):
        '''
        test kmers length of 3
        '''
        input_seq = "IAMCUTE"
        k = 3
        kmer_list = list(get_kmers(input_seq, k))
        expected = ["IAM","AMC","MCU","CUT","UTE"]
        self.assertEqual(kmer_list.sort(),expected.sort())

    def test_getkmers2(self):
        '''
        test kmers length of 4
        '''
        input_seq = "IAMCCUTE"
        k = 4
        kmer_list = list(get_kmers(input_seq, k))
        expected = ["IAM","AMCC","MCCU","CCUT","CUTE"]
        self.assertEqual(kmer_list.sort(),expected.sort())

    def test_getkmer3(self):
        '''
        test kmers length of 9; k > length of the sequence
        '''
        input_seq = "IAMCCUTE"
        k = 9
        with unittest.mock.patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            get_kmers(input_seq, k)
            self.assertEqual(
                mock_stdout.getvalue(),
                "kmer must be shorter than the sequence length\n"  
            )

class Numkey_tests(unittest.TestCase):
    def test_get_numerical_rep1(self):
        '''
        test a sequence with A T C G
        '''
        input_seq = "ATCGATCGGGT"
        numkey = get_numerical_rep(input_seq)
        expected = 888235
        self.assertEqual(numkey,expected)
        

    def test_get_numerical_rep2(self):
        '''
        test a sequence with A T C G and other characters
        '''
        input_seq = "NNNNNAT"
        numkey = get_numerical_rep(input_seq)
        expected = 3
        self.assertEqual(numkey,expected)

    def test_get_numerical_rep3(self):
        '''
        test an empty sequence
        '''
        input_seq = ""
        numkey = get_numerical_rep(input_seq)
        expected = None
        self.assertEqual(numkey,expected)
        

class Jaccard_tests(unittest.TestCase):

    def test_jaccard1(self):
        '''
        test jaccard with overlaps
        '''
        A = [1,2,3,4]
        B = [3,4,5,6,7,8]
        jaccard_res = jaccard(A,B)
        expected = 0.25
        self.assertEqual(jaccard_res,expected)

    def test_jaccard2(self):
        '''
        test jaccard with non overlaps
        '''
        A = [1,2,3,4]
        B = [100,89]
        jaccard_res = jaccard(A,B)
        expected = 0
        self.assertEqual(jaccard_res,expected)

if __name__ == '__main__':
    unittest.main()
