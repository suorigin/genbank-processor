# GenBank Processor

A comprehensive command-line tool for processing, filtering, downloading, and renaming organelle genomic data from NCBI GenBank. This tool provides a streamlined workflow for researchers working with complete  organelle genome data.

## Dependencies

- Python 3.8+
- requests >= 2.25.0
- biopython >= 1.79
- pandas >= 1.3.0

## Features

- **Filtering**: Extract complete organelle genome records from raw data files
- **Conversion**: Transform filtered data into TSV format for easy analysis
- **Download**: Batch download GenBank files using accession numbers
- **Pipeline**: Complete workflow from raw data to downloaded GenBank files
- **Rename gbfiles**: Rename GenBank files based on a mapping table
- **Extractcp**: Extract chloroplast genome sequences and rename them
- **Extract Gene**: Extract CDS and protein sequences from GenBank files
- **Simplify ID**: Simplify FASTA sequence IDs to gene_short_accession format
- **Multilingual Support**: Help documentation available in English and Chinese

## Installation
### Prerequisites

- `Python` 3.8 or higher
- `requests`, `biopython`, and `pandas` libraries

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

### Rename GenBank files

```bash
genbank-processor renamegbfiles -i gb_files -r replace.txt -o renamed_files
```

### Extract chloroplast sequences

```bash
genbank-processor extractcp -i gb_files -o chloroplast_sequences -r replace.txt
```

### Extract CDS and protein sequences

```bash
genbank-processor extract_gene -i gb_files -o gene_sequences -t both
```

### Simplify sequence IDs

```bash
genbank-processor simplify_id -i gene_sequences -o simplified_sequences
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

### Renamegbfiles Command

Renames GenBank files based on a mapping table.

```
genbank-processor renamegbfiles -i INPUT -r REPLACE -o OUTPUT
```

- `-i, --input`: Input directory containing GB files (required)
- `-r, --replace`: Mapping table file path (TSV format) (required)
- `-o, --output`: Output directory path (required)

### Extractcp Command

Extracts chloroplast genome sequences from GenBank files and optionally renames them.

```
genbank-processor extractcp -i INPUT -o OUTPUT [-r REPLACE]
```

- `-i, --input`: Input directory containing GB files (required)
- `-o, --output`: Output directory path (required)
- `-r, --replace`: Mapping table file path (TSV format) (optional)

### Extract Gene Command

Extracts CDS and/or protein sequences from GenBank files.

```
genbank-processor extract_gene -i INPUT -o OUTPUT [-t TYPE] [-r RENAME]
```

- `-i, --input`: Input directory containing GB files (required)
- `-o, --output`: Output directory path (required)
- `-t, --type`: Extraction type: cds, protein, or both (default: both)
- `-r, --rename`: Mapping table file path for renaming output files (optional)

### Simplify ID Command

Simplifies FASTA sequence IDs to gene_short_accession format.

```
genbank-processor simplify_id -i INPUT_DIR -o OUTPUT_DIR
# or
genbank-processor simplify_id -I INPUT_FILE -O OUTPUT_FILE
```

- `-i, --input-dir`: Input directory containing FASTA files
- `-o, --output-dir`: Output directory path
- `-I, --input-file`: Input FASTA file path
- `-O, --output-file`: Output FASTA file path

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

The renamegbfiles and extractcp commands require a mapping table in TSV format with two columns:

Example:
```
Accession	Name
PP988495.1	DASZ
PV915199.1	WYLC
PQ510294.1	GL
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

### Renamegbfiles Output

The renamegbfiles command saves renamed GenBank files in the specified output directory, using the names from the mapping table.

### Extractcp Output

The extractcp command saves FASTA files containing chloroplast genome sequences, optionally renamed according to the mapping table.

### Extract Gene Output

The extract_gene command saves FASTA files containing CDS and/or protein sequences, with file names following the pattern `{accession}_{type}.fasta` (e.g., `PP988495.1_cds.fasta`).

### Simplify ID Output

The simplify_id command saves FASTA files with simplified sequence IDs in the format `gene_short.accession` (e.g., `rps12.1_PP988495.1`).

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

# Rename GenBank files
genbank-processor renamegbfiles -i gb_files/ -r replace.txt -o renamed_files/

# Extract chloroplast sequences
genbank-processor extractcp -i gb_files/ -o chloroplast_sequences/ -r replace.txt

# Extract CDS and protein sequences
genbank-processor extract_gene -i gb_files/ -o gene_sequences/ -t both

# Simplify sequence IDs
genbank-processor simplify_id -i gene_sequences/ -o simplified_sequences/
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
- Author's email: [1436636379@qq.com]

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Support

If you encounter any issues or have questions, please open an issue on the GitHub repository.

## Acknowledgments

- NCBI for providing the GenBank database and E-utilities API
- The bioinformatics community for valuable feedback and testing

---

*GenBank Processor - Streamlining genomic data processing for researchers*
```
