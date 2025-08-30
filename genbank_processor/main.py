#!/usr/bin/env python3
"""
GenBank Processor - A comprehensive tool for processing genomic data from NCBI
"""

import argparse
import os
import sys
import locale
import time
import requests
from urllib.parse import urlencode
from urllib.error import HTTPError

class HelpSystem:
    """Multilingual help system"""
    
    @staticmethod
    def get_help_texts(lang=None):
        """Get help texts in specified language"""
        if lang is None:
            lang = HelpSystem.detect_language()
        
        texts = {
            'en': {
                'description': "GenBank Processor - A tool for genomic data processing",
                'filter_help': "Filter complete genomes data from raw data",
                'convert_help': "Convert filtered data to TSV format",
                'download_help': "Download GenBank files from NCBI",
                'pipeline_help': "Run the complete processing pipeline",
                'input_help': "Input file path",
                'output_help': "Output file path",
                'delay_help': "Delay between requests in seconds",
                'retries_help': "Number of retry attempts",
                'tsv_help': "Intermediate TSV file path",
                'help_help': "Show help message in specified language (en/zh, default: en)",
                'verbose_help': "Enable verbose output",
                'examples': "\nExamples:\n  genbank-processor filter -i input.txt -o output.txt\n  genbank-processor pipeline -i data.txt -o results/\n",
            },
            'zh': {
                'description': "GenBank Processor - 基因组数据处理工具",
                'filter_help': "从原始数据中过滤，保留完整基因组信息",
                'convert_help': "将过滤后的数据转换为TSV格式",
                'download_help': "从NCBI下载GenBank文件",
                'pipeline_help': "运行完整处理流程",
                'input_help': "输入文件路径",
                'output_help': "输出文件路径",
                'delay_help': "请求之间的延迟时间（秒）",
                'retries_help': "重试次数",
                'tsv_help': "中间TSV文件路径",
                'help_help': "以指定语言显示帮助信息 (英文/中文，默认：英文)",
                'verbose_help': "显示详细输出",
                'examples': "\n示例:\n  genbank-processor filter -i input.txt -o output.txt\n  genbank-processor pipeline -i data.txt -o results/\n",
            }
        }
        
        return texts.get(lang, texts['en'])
    
    @staticmethod
    def detect_language():
        """Detect user's preferred language for help"""
        # 1. Check command line arguments
        for i, arg in enumerate(sys.argv):
            if arg in ['-h', '--help'] and i + 1 < len(sys.argv):
                lang = sys.argv[i + 1].lower()
                if lang in ['zh', 'cn', 'chinese']:
                    return 'zh'
                elif lang in ['en', 'english']:
                    return 'en'
        
        # 2. Check environment variable
        lang_env = os.environ.get('GENBANK_LANG', '').lower()
        if lang_env in ['zh', 'cn', 'chinese']:
            return 'zh'
        elif lang_env in ['en', 'english']:
            return 'en'
        
        # 3. Check system language
        try:
            # Try to get current locale setting
            current_locale = locale.getlocale()
            if current_locale and current_locale[0] and 'zh' in current_locale[0].lower():
                return 'zh'
            
            # If not set, check LANG environment variable
            lang_env = os.environ.get('LANG', '').lower()
            if lang_env and 'zh' in lang_env:
                return 'zh'
                
        except Exception:
            # Fall back to English if any exception occurs
            pass
        
        # 4. Default to English
        return 'en'

def filter_complete_genomes(input_file, output_file, verbose=False):
    """Filter records containing complete genomes"""
    if verbose:
        print(f"Filtering file: {input_file}")
    
    try:
        with open(input_file, 'r', encoding='utf-8') as f:
            content = f.read()
        
        units = content.strip().split('\n\n')
        filtered_units = []
        
        for unit in units:
            lines = unit.split('\n')
            if len(lines) >= 3 and "complete genome" in lines[0].lower():
                filtered_units.append(unit)
        
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write('\n\n'.join(filtered_units))
        
        if verbose:
            print(f"Filtering completed! Found {len(filtered_units)} complete genomes")
            print(f"Results saved to: {output_file}")
            
        return len(filtered_units)
        
    except Exception as e:
        print(f"Error during filtering: {str(e)}")
        sys.exit(1)

def convert_to_table(input_file, output_file, verbose=False):
    """Convert filtered data to TSV format"""
    if verbose:
        print(f"Converting file: {input_file}")
    
    try:
        with open(input_file, 'r', encoding='utf-8') as f:
            content = f.read().strip()
        
        units = content.split('\n\n')
        table_data = []
        
        for unit in units:
            lines = unit.split('\n')
            if len(lines) >= 3 and "complete genome" in lines[0].lower():
                species_info = lines[0].strip()
                sequence_info = lines[1].strip()
                accession_info = lines[2].strip()
                table_data.append([species_info, sequence_info, accession_info])
        
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("Species Information\tSequence Information\tAccession Number\n")
            for row in table_data:
                f.write("\t".join(row) + "\n")
        
        if verbose:
            print(f"Conversion completed! Generated {len(table_data)} rows of data")
            print(f"Results saved to: {output_file}")
            
        return len(table_data)
        
    except Exception as e:
        print(f"Error during conversion: {str(e)}")
        sys.exit(1)

def extract_accession_numbers(tsv_file):
    """Extract accession numbers from TSV file"""
    accessions = []
    try:
        with open(tsv_file, 'r', encoding='utf-8') as f:
            next(f)  # Skip header
            for line in f:
                columns = line.strip().split('\t')
                if len(columns) >= 3:
                    accession_field = columns[2]
                    accession = accession_field.split()[0]
                    accessions.append(accession)
        return accessions
    except Exception as e:
        print(f"Error extracting accession numbers: {str(e)}")
        sys.exit(1)

def download_genbank(accession, output_dir, retries=3, delay=2, verbose=False):
    """Download GenBank data for a single accession number"""
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    
    params = {
        'db': 'nucleotide',
        'id': accession,
        'rettype': 'gb',
        'retmode': 'text'
    }
    
    url = f"{base_url}?{urlencode(params)}"
    
    for attempt in range(retries):
        try:
            response = requests.get(url)
            response.raise_for_status()
            
            output_path = os.path.join(output_dir, f"{accession}.gb")
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(response.text)
            
            if verbose:
                print(f"Successfully downloaded: {accession}")
            return True
            
        except (HTTPError, requests.RequestException) as e:
            if verbose:
                print(f"Failed to download {accession} (attempt {attempt+1}/{retries}): {str(e)}")
            time.sleep(delay)
    
    if verbose:
        print(f"Unable to download {accession}, skipping...")
    return False

def download_genbank_files(tsv_file, output_dir, delay=1, retries=3, verbose=False):
    """Batch download GenBank files"""
    if verbose:
        print(f"Starting download of GenBank files from {tsv_file}")
    
    try:
        os.makedirs(output_dir, exist_ok=True)
        accessions = extract_accession_numbers(tsv_file)
        
        if not accessions:
            print("No valid accession numbers found")
            return 0
        
        if verbose:
            print(f"Found {len(accessions)} accession numbers, starting download...")
        
        success_count = 0
        for i, accession in enumerate(accessions):
            if verbose:
                print(f"Processing ({i+1}/{len(accessions)}): {accession}")
            
            if download_genbank(accession, output_dir, retries, delay, verbose):
                success_count += 1
            
            time.sleep(delay)  # Respect NCBI rate limits
        
        if verbose:
            print(f"Download completed! Successfully downloaded {success_count}/{len(accessions)} files")
            print(f"Files saved to: {output_dir}")
        
        return success_count
        
    except Exception as e:
        print(f"Error during download: {str(e)}")
        sys.exit(1)

def create_parser(lang='en'):
    """Create command line parser"""
    texts = HelpSystem.get_help_texts(lang)
    
    # Create main parser
    parser = argparse.ArgumentParser(
        description=texts['description'],
        epilog=texts['examples'],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False  # Disable default help
    )
    
    # Add help parameter
    parser.add_argument(
        '-h', '--help', 
        nargs='?', 
        const='en',  # Default value
        choices=['en', 'zh'],
        help=texts['help_help']
    )
    
    # Add verbose output parameter
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help=texts['verbose_help']
    )
    
    # Add subcommands
    subparsers = parser.add_subparsers(dest='command', help='available commands')
    
    # filter command
    filter_parser = subparsers.add_parser('filter', help=texts['filter_help'])
    filter_parser.add_argument('-i', '--input', required=True, help=texts['input_help'])
    filter_parser.add_argument('-o', '--output', required=True, help=texts['output_help'])
    
    # convert command
    convert_parser = subparsers.add_parser('convert', help=texts['convert_help'])
    convert_parser.add_argument('-i', '--input', required=True, help=texts['input_help'])
    convert_parser.add_argument('-o', '--output', required=True, help=texts['output_help'])
    
    # download command
    download_parser = subparsers.add_parser('download', help=texts['download_help'])
    download_parser.add_argument('-i', '--input', required=True, help=texts['input_help'])
    download_parser.add_argument('-o', '--output', required=True, help=texts['output_help'])
    download_parser.add_argument('--delay', type=float, default=1.0, help=texts['delay_help'])
    download_parser.add_argument('--retries', type=int, default=3, help=texts['retries_help'])
    
    # pipeline command
    pipeline_parser = subparsers.add_parser('pipeline', help=texts['pipeline_help'])
    pipeline_parser.add_argument('-i', '--input', required=True, help=texts['input_help'])
    pipeline_parser.add_argument('-t', '--tsv', default='genomes.tsv', help=texts['tsv_help'])
    pipeline_parser.add_argument('-o', '--output', required=True, help=texts['output_help'])
    pipeline_parser.add_argument('--delay', type=float, default=1.0, help=texts['delay_help'])
    pipeline_parser.add_argument('--retries', type=int, default=3, help=texts['retries_help'])
    
    return parser

def main():
    """Main function"""
    # Detect language
    help_lang = HelpSystem.detect_language()
    
    # Create parser
    parser = create_parser(help_lang)
    
    # Parse arguments
    args, remaining = parser.parse_known_args()
    
    # Handle help request
    if hasattr(args, 'help') and args.help:
        # If language is specified, recreate parser with that language
        if args.help in ['en', 'zh']:
            help_lang = args.help
            parser = create_parser(help_lang)
        
        # Print help information
        parser.print_help()
        return
    
    if not hasattr(args, 'command') or not args.command:
        parser.print_help()
        return
    
    try:
        if args.command == 'filter':
            print(f"Executing filter command")
            print(f"Input file: {args.input}")
            print(f"Output path: {args.output}")
            
            # Call the actual filtering function
            count = filter_complete_genomes(args.input, args.output, args.verbose)
            if count > 0:
                print(f"Successfully filtered {count} complete genome records")
            else:
                print("No complete genome records found")
            
        elif args.command == 'convert':
            print(f"Executing convert command")
            print(f"Input file: {args.input}")
            print(f"Output path: {args.output}")
            
            # Call the actual conversion function
            count = convert_to_table(args.input, args.output, args.verbose)
            if count > 0:
                print(f"Successfully converted {count} records")
            else:
                print("No records found for conversion")
            
        elif args.command == 'download':
            print(f"Executing download command")
            print(f"Input file: {args.input}")
            print(f"Output directory: {args.output}")
            print(f"Request delay: {args.delay} seconds")
            print(f"Retry attempts: {args.retries}")
            
            # Call the actual download function
            success_count = download_genbank_files(args.input, args.output, args.delay, args.retries, args.verbose)
            if success_count > 0:
                print(f"Successfully downloaded {success_count} files")
            else:
                print("Download failed")
            
        elif args.command == 'pipeline':
            print(f"Executing complete pipeline")
            print(f"Input file: {args.input}")
            print(f"Intermediate TSV file: {args.tsv}")
            print(f"Output directory: {args.output}")
            print(f"Request delay: {args.delay} seconds")
            print(f"Retry attempts: {args.retries}")
            
            # Create temporary file
            temp_file = "filtered_temp.txt"
            
            # Step 1: Filter
            filter_count = filter_complete_genomes(args.input, temp_file, args.verbose)
            if filter_count == 0:
                print("No complete genomes found, terminating pipeline")
                os.remove(temp_file)
                return
            
            # Step 2: Convert
            convert_count = convert_to_table(temp_file, args.tsv, args.verbose)
            
            # Step 3: Download
            download_count = download_genbank_files(args.tsv, args.output, args.delay, args.retries, args.verbose)
            
            # Clean up temporary file
            os.remove(temp_file)
            
            if args.verbose:
                print(f"Pipeline completed! Filtered: {filter_count}, Converted: {convert_count}, Downloaded: {download_count}")
    
    except KeyboardInterrupt:
        print("\nOperation interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
