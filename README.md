# GenBank Processor

A comprehensive command-line tool for processing, filtering, and downloading organelle genomic data from NCBI GenBank. This tool provides a streamlined workflow for researchers working with complete  organelle genome data.

## Dependencies

- Python 3.8+
- requests >= 2.25.0

## Features

- **Filtering**: Extract complete organelle genome records from raw data files
- **Conversion**: Transform filtered data into TSV format for easy analysis
- **Download**: Batch download GenBank files using accession numbers
- **Pipeline**: Complete workflow from raw data to downloaded GenBank files
- **Multilingual Support**: Help documentation available in English and Chinese

## Installation
### Prerequisites

- `Python` 3.8 or higher
- `requests` library

### From source
```bash
git clone https://github.com/yourusername/genbank-processor.git
cd genbank-processor
pip install .
# or
python setup.py install
```

## Usage

GenBank Processor provides several subcommands for different processing tasks:

### Filter complete genomes

```bash
genbank-processor filter -i input.txt -o filtered.txt
```

### Convert to TSV format

```bash
genbank-processor convert -i filtered.txt -o genomes.tsv
```

### Download GenBank files

```bash
genbank-processor download -i genomes.tsv -o downloads/ --delay 1.5 --retries 5
```

### Run complete pipeline

```bash
genbank-processor pipeline -i raw_data.txt -t genomes.tsv -o downloads/ --delay 1.0 --retries 3
```

## Command Reference

### Global Options

- `-h [en|zh], --help [en|zh]`: Show help message in specified language (default: English)
- `-v, --verbose`: Enable verbose output

### Filter Command

Extracts complete genome records from input data.

```
genbank-processor filter -i INPUT -o OUTPUT
```

- `-i, --input`: Input file path (required)
- `-o, --output`: Output file path (required)

### Convert Command

Converts filtered data to TSV format with three columns: Species Information, Sequence Information, and Accession Number.

```
genbank-processor convert -i INPUT -o OUTPUT
```

- `-i, --input`: Input file path (required)
- `-o, --output`: Output TSV file path (required)

### Download Command

Downloads GenBank files using accession numbers from a TSV file.

```
genbank-processor download -i INPUT -o OUTPUT [--delay DELAY] [--retries RETRIES]
```

- `-i, --input`: Input TSV file path (required)
- `-o, --output`: Output directory path (required)
- `--delay`: Delay between requests in seconds (default: 1.0)
- `--retries`: Number of retry attempts for failed downloads (default: 3)

### Pipeline Command

Runs the complete processing pipeline: filter → convert → download.

```
genbank-processor pipeline -i INPUT -o OUTPUT [-t TSV] [--delay DELAY] [--retries RETRIES]
```

- `-i, --input`: Raw input file path (required)
- `-o, --output`: Output directory path (required)
- `-t, --tsv`: Intermediate TSV file path (default: genomes.tsv)
- `--delay`: Delay between requests in seconds (default: 1.0)
- `--retries`: Number of retry attempts for failed downloads (default: 3)

## Input Format

The tool expects input files with the following format:

```
Species description line (must contain "complete genome")
Sequence information line
Accession number line

[Additional records separated by blank lines]
```

Example:
```
1.Homo sapiens complete genome
Sequence information for Homo sapiens
NC_000001.1 GI:XXXXXXXX

2.Mus musculus complete genome
Sequence information for Mus musculus
NC_000067.6 GI:XXXXXXXX
```

## Output Format

### Filter Output

The filter command produces output with the same format as the input, but only containing records that include "complete genome" in the species description.

### Convert Output

The convert command produces a TSV file with the following columns:
- Species Information
- Sequence Information
- Accession Number

### Download Output

The download command saves individual GenBank files (.gb format) named by their accession numbers in the specified output directory.

## Multilingual Support

GenBank Processor provides help documentation in both English and Chinese. To view help in a specific language:

```bash
# English help
genbank-processor -h en --help

# Chinese help
genbank-processor -h zh --help
```

The tool automatically detects the user's preferred language based on:
1. Command line arguments
2. GENBANK_LANG environment variable
3. System locale settings

## Examples

There is a Summary_result.txt file in the tests/test_data folder for demonstration purposes.
You can switch to the tests/test_data folder to conduct test learning.
```bash
cd ./tests/test_data
```

### Basic Usage

```bash
# Filter complete genomes from raw data
genbank-processor filter -i Summary_result.txt -o complete_genomes.txt

# Convert to TSV format
genbank-processor convert -i complete_genomes.txt -o genomes.tsv

# Download GenBank files
genbank-processor download -i genomes.tsv -o gb_files/ --delay 2.0
```

### Complete Pipeline

```bash
# Run the complete pipeline with verbose output
genbank-processor pipeline -i Summary_result.txt -o gb_files/ -t genomes.tsv --delay 2.0 
```

### Using Different Languages

```bash
# View help in Chinese
genbank-processor -h zh --help

# View filter command help in English
genbank-processor filter -h en --help
```

## Error Handling

The tool includes robust error handling with:
- Automatic retries for failed downloads (configurable with --retries)
- Graceful handling of network timeouts
- Clear error messages for common issues
- Proper cleanup on interruption (Ctrl+C)


## Contributing

- Contributions are welcome! Please feel free to submit a Pull Request.
- For more information, you can follow the user:(https://www.zhihu.com/people/su-dian-dian-87-36)
- Author’s email: [1436636379@qq.com]

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Support

If you encounter any issues or have questions, please open an issue on the GitHub repository.

## Acknowledgments

- NCBI for providing the GenBank database and E-utilities API
- The bioinformatics community for valuable feedback and testing

---

*GenBank Processor - Streamlining genomic data processing for researchers*
