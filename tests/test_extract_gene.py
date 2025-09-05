#!/usr/bin/env python3
"""
Tests for extract_gene functionality
"""

import os
import tempfile
import unittest
from genbank_processor.main import extract_sequences_from_genbank, batch_extract_sequences

class TestExtractGene(unittest.TestCase):
    """Test extract_gene functionality"""
    
    def setUp(self):
        """Set up test environment"""
        self.test_dir = tempfile.mkdtemp()
        
    def tearDown(self):
        """Clean up test environment"""
        import shutil
        shutil.rmtree(self.test_dir)
    
    def test_extract_sequences(self):
        """Test extracting sequences from a GenBank file"""
        # This test requires a sample GenBank file
        # You should add a sample.gb file to tests/test_data/
        test_gb_file = "tests/test_data/sample.gb"
        
        if os.path.exists(test_gb_file):
            output_dir = os.path.join(self.test_dir, "output")
            os.makedirs(output_dir, exist_ok=True)
            
            cds_count, protein_count, accession = extract_sequences_from_genbank(
                test_gb_file, output_dir, "both"
            )
            
            self.assertIsInstance(cds_count, int)
            self.assertIsInstance(protein_count, int)
            self.assertIsInstance(accession, str)
            
            # Check if output files were created
            cds_file = os.path.join(output_dir, f"{accession}_cds.fasta")
            protein_file = os.path.join(output_dir, f"{accession}_protein.fasta")
            
            if cds_count > 0:
                self.assertTrue(os.path.exists(cds_file))
            
            if protein_count > 0:
                self.assertTrue(os.path.exists(protein_file))

if __name__ == "__main__":
    unittest.main()
