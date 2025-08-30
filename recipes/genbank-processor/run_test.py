#!/usr/bin/env python3

import subprocess
import sys

def test_cli():
    """Test that the command line interface works"""
    # Test English help
    result = subprocess.run([
        sys.executable, '-m', 'genbank_processor.main', '--help'
    ], capture_output=True, text=True)
    
    assert result.returncode == 0
    assert 'GenBank Processor' in result.stdout
    print("CLI test passed!")

if __name__ == "__main__":
    test_cli()
