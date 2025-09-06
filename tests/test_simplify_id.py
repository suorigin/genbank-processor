#!/usr/bin/env python3
"""
Tests for simplify_id functionality
"""

import os
import tempfile
import unittest
from genbank_processor.main import simplify_sequence_id, batch_simplify_ids

class TestSimplifyID(unittest.TestCase):
    """Test simplify_id functionality"""
    
    def setUp(self):
        """Set up test environment"""
        self.test_dir = tempfile.mkdtemp()
        
    def tearDown(self):
        """Clean up test environment"""
        import shutil
        shutil.rmtree(self.test_dir)
    
    def test_simplify_sequence_id(self):
        """Test simplifying sequence IDs in a FASTA file"""
        # Create a test FASTA file
        test_fasta = os.path.join(self.test_dir, "test.fasta")
        output_fasta = os.path.join(self.test_dir, "test_simplified.fasta")
        
        with open(test_fasta, 'w') as f:
            f.write(""">rps12 ribosomal protein S12 | PP988495.1 | join{[71952:72066](-), [100414:100646](-), [99852:99878](-)}
ATGCCAACTATTAAACAACTTATTAGAAATACAAGACAGCCAATCAGAAATGTCACGAAA
>rps12 ribosomal protein S12 | PP988495.1 | join{[71952:72066](-), [100414:100646](-), [99852:99878](-)}
ATGCCAACTATTAAACAACTTATTAGAAATACAAGACAGCCAATCAGAAATGTCACGAAA
""")
        
        processed_sequences = simplify_sequence_id(test_fasta, output_fasta)
        
        self.assertEqual(processed_sequences, 2)
        self.assertTrue(os.path.exists(output_fasta))
        
        # Check if the output file contains simplified IDs
        with open(output_fasta, 'r') as f:
            content = f.read()
            self.assertIn("rps12.1_PP988495.1", content)
            self.assertIn("rps12.2_PP988495.1", content)

if __name__ == "__main__":
    unittest.main()
