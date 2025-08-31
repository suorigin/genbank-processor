#!/usr/bin/env python3
"""
Basic tests for GenBank Processor
"""

def test_import():
    """Test that the package can be imported"""
    try:
        import genbank_processor
        assert True
    except ImportError as e:
        assert False, f"Failed to import genbank_processor: {e}"

def test_help_command():
    """Test that the help command works"""
    import subprocess
    import sys
    
    try:
        # Test command line help
        result = subprocess.run([
            sys.executable, '-m', 'genbank_processor.main', '--help'
        ], capture_output=True, text=True, timeout=30)
        
        # Check if the command ran successfully
        assert result.returncode == 0, f"Command failed with return code {result.returncode}"
        
        # Check if the expected text is in the output
        assert 'GenBank Processor' in result.stdout or 'genbank-processor' in result.stdout
        
    except subprocess.TimeoutExpired:
        assert False, "Command timed out"
    except Exception as e:
        assert False, f"Unexpected error: {e}"

def test_dummy():
    """A simple dummy test to ensure at least one test runs"""
    assert True
