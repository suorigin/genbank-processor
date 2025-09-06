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
import shutil
import re
import pandas as pd
from urllib.parse import urlencode
from urllib.error import HTTPError
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

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
                'renamegbfiles_help': "Rename GenBank files based on mapping table",
                'extractcp_help': "Extract complete genome sequences from GB files and rename IDs",
                'pipeline_help': "Run the complete processing pipeline",
                'extract_gene_help': "Extract CDS and/or protein sequences from GenBank files",
                'simplify_id_help': "Simplify FASTA sequence IDs to gene_short_accession format",
                'input_help': "Input file path",
                'output_help': "Output file path",
                'delay_help': "Delay between requests in seconds",
                'retries_help': "Number of retry attempts",
                'tsv_help': "Intermediate TSV file path",
                'mapping_help': "Mapping table file path (TSV format)",
                'type_help': "Extraction type: cds, protein or both (default: both)",
                'help_help': "Show help message in specified language (en/zh, default: en)",
                'verbose_help': "Enable verbose output",
                'examples': "\nExamples:\n  genbank-processor filter -i input.txt -o output.txt\n  genbank-processor pipeline -i data.txt -o results/\n  genbank-processor renamegbfiles -i gb_files -r mapping.tsv -o renamed_files\n  genbank-processor extractcp -i gb_files -o chloroplast_sequences -r mapping.tsv\n  genbank-processor extract_gene -i gb_files -o gene_sequences -t both\n  genbank-processor simplify_id -i gene_sequences -o simplified_sequences\n",
            },
            'zh': {
                'description': "GenBank Processor - 基因组数据处理工具",
                'filter_help': "从原始数据中过滤，保留完整基因组信息",
                'convert_help': "将过滤后的数据转换为TSV格式",
                'download_help': "从NCBI下载GenBank文件",
                'renamegbfiles_help': "根据映射表重命名GenBank文件",
                'extractcp_help': "从GB文件中提取完整基因组序列并重命名ID",
                'pipeline_help': "运行完整处理流程",
                'extract_gene_help': "从GenBank文件中提取CDS和/或蛋白质序列",
                'simplify_id_help': "简化FASTA序列ID为基因简称_登录号格式",
                'input_help': "输入文件路径",
                'output_help': "输出文件路径",
                'delay_help': "请求之间的延迟时间（秒）",
                'retries_help': "重试次数",
                'tsv_help': "中间TSV文件路径",
                'mapping_help': "映射表文件路径 (TSV格式)",
                'type_help': "提取类型: cds, protein 或 both (默认: both)",
                'help_help': "以指定语言显示帮助信息 (英文/中文，默认：英文)",
                'verbose_help': "显示详细输出",
                'examples': "\n示例:\n  genbank-processor filter -i input.txt -o output.txt\n  genbank-processor pipeline -i data.txt -o results/\n  genbank-processor renamegbfiles -i gb_files -r mapping.tsv -o renamed_files\n  genbank-processor extractcp -i gb_files -o chloroplast_sequences -r mapping.tsv\n  genbank-processor extract_gene -i gb_files -o gene_sequences -t both\n  genbank-processor simplify_id -i gene_sequences -o simplified_sequences\n",
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

def load_mapping(mapping_file, verbose=False):
    """Load mapping table file"""
    mapping = {}
    try:
        with open(mapping_file, 'r', encoding='utf-8') as f:
            # Skip header
            next(f)
            for line_num, line in enumerate(f, 2):  # Start counting from line 2
                line = line.strip()
                if not line:  # Skip empty lines
                    continue
                    
                parts = line.split('\t')
                if len(parts) < 2:
                    if verbose:
                        print(f"Warning: Line {line_num} format incorrect, skipped: {line}")
                    continue
                    
                old_name, new_name = parts[0], parts[1]
                if old_name in mapping:
                    if verbose:
                        print(f"Warning: Duplicate old name '{old_name}' at line {line_num}, will use the last occurrence")
                
                mapping[old_name] = new_name
                
        if verbose:
            print(f"Successfully loaded mapping table: {mapping_file}")
            print(f"Found {len(mapping)} mapping relationships")
        return mapping
    except Exception as e:
        print(f"Error: Unable to load mapping table - {str(e)}")
        return None

def renamegbfiles_gb_files(input_dir, output_dir, mapping, verbose=False):
    """Rename GB files"""
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Statistics results
    results = {
        'processed': 0,
        'success': [],
        'skipped': [],
        'failed': []
    }
    
    # Traverse all GB files in the input directory
    for filename in os.listdir(input_dir):
        if not filename.endswith(".gb"):
            continue
            
        results['processed'] += 1
        input_path = os.path.join(input_dir, filename)
        
        # Extract filename prefix (without extension)
        file_prefix = os.path.splitext(filename)[0]
        
        # Check if the file prefix is in the mapping table
        if file_prefix in mapping:
            new_name = mapping[file_prefix]
            output_filename = f"{new_name}.gb"
            output_path = os.path.join(output_dir, output_filename)
            
            try:
                # Copy file to output directory and rename
                shutil.copy2(input_path, output_path)
                results['success'].append((filename, output_filename))
                if verbose:
                    print(f"✓ Successfully renamed: {filename} -> {output_filename}")
            except Exception as e:
                error_msg = f"Failed to copy file: {str(e)}"
                results['failed'].append((filename, error_msg))
                if verbose:
                    print(f"✗ Failed: {filename} | Error: {error_msg}")
        else:
            results['skipped'].append(filename)
            if verbose:
                print(f"Warning: Prefix '{file_prefix}' of {filename} not found in mapping table, skipped")
    
    return results

def print_renamegbfiles_summary(results, output_dir):
    """Print rename processing summary"""
    print("\n" + "=" * 60)
    print("Rename Processing Summary")
    print("=" * 60)
    print(f"Total processed: {results['processed']} files")
    print(f"Successfully renamed: {len(results['success'])} files")
    print(f"Skipped: {len(results['skipped'])} files")
    print(f"Failed: {len(results['failed'])} files")
    print(f"Results saved to: {output_dir}")
    
    # Display successfully renamed files
    if results['success']:
        print("\nSuccessfully renamed files:")
        for old_name, new_name in results['success']:
            print(f"  {old_name} -> {new_name}")
    
    # Display skipped files
    if results['skipped']:
        print("\nSkipped files (not found in mapping table):")
        for filename in results['skipped']:
            print(f"  {filename}")
    
    # Display failed files
    if results['failed']:
        print("\nFailed files:")
        for filename, error in results['failed']:
            print(f"  {filename}: {error}")

def extract_chloroplast_sequence(gb_file):
    """
    Extract chloroplast genome sequence from a single GenBank file
    """
    records = list(SeqIO.parse(gb_file, "genbank"))
    
    # Check if chloroplast sequence is found
    chloroplast_found = False
    for record in records:
        # Look for chloroplast features
        if "chloroplast" in record.description.lower() or \
           any("chloroplast" in feat.qualifiers.get("organelle", [""])[0].lower() 
              for feat in record.features if feat.type == "source"):
            
            chloroplast_found = True
            # Get accession number and species name
            accession = record.id
            organism = record.annotations.get("organism", "Unknown_organism")
            
            # Create sequence description
            description = f"{accession}"
            
            # Return sequence record
            return SeqRecord(record.seq, id=accession, description=description)
    
    # If no chloroplast sequence is found
    if not chloroplast_found:
        print(f"Warning: No chloroplast sequence found in {os.path.basename(gb_file)}")
        return None

def load_mapping_fasta(mapping_file, verbose=False):
    """Load mapping table file for FASTA renaming"""
    mapping = {}
    try:
        with open(mapping_file, 'r', encoding='utf-8') as f:
            # Skip possible header line (try to detect)
            first_line = f.readline().strip()
            
            # Check if first line contains header features (like "Accession" or "Name")
            if "accession" in first_line.lower() and "name" in first_line.lower():
                # This is a header, start reading from next line
                pass
            else:
                # This is not a header, reset file pointer
                f.seek(0)
            
            # Read all lines
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:  # Skip empty lines
                    continue
                    
                parts = line.split('\t')
                if len(parts) < 2:
                    if verbose:
                        print(f"Warning: Line {line_num} format incorrect, skipped: {line}")
                    continue
                    
                old_name, new_name = parts[0], parts[1]
                if old_name in mapping:
                    if verbose:
                        print(f"Warning: Duplicate old name '{old_name}' at line {line_num}, will use the last occurrence")
                
                mapping[old_name] = new_name
                
        if verbose:
            print(f"Successfully loaded mapping table: {mapping_file}")
            print(f"Found {len(mapping)} mapping relationships")
        return mapping
    except Exception as e:
        print(f"Error: Unable to load mapping table - {str(e)}")
        return None

def rename_fasta_files(input_dir, output_dir, mapping, verbose=False):
    """
    Rename FASTA files and replace sequence IDs
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    processed = 0
    success = 0
    skipped = 0
    failed = 0
    
    # Traverse all FASTA files in the input directory
    for filename in os.listdir(input_dir):
        if not filename.endswith(".fasta"):
            continue
            
        processed += 1
        input_path = os.path.join(input_dir, filename)
        
        # Extract filename prefix (without extension)
        file_prefix = os.path.splitext(filename)[0]
        
        # Check if the file prefix is in the mapping table
        if file_prefix in mapping:
            species_name = mapping[file_prefix]
            if verbose:
                print(f"Processing: {filename} -> {species_name}")
            
            try:
                # Read FASTA file
                records = list(SeqIO.parse(input_path, "fasta"))
                
                # Create new sequence records (replace ID)
                new_records = []
                for record in records:
                    # Create new record, keeping only species name as ID
                    new_record = SeqRecord(
                        seq=record.seq,
                        id=species_name,
                        description="",  # Clear description information
                    )
                    new_records.append(new_record)
                
                # Create output filename
                output_filename = f"{species_name}.fasta"
                output_path = os.path.join(output_dir, output_filename)
                
                # Save new file
                SeqIO.write(new_records, output_path, "fasta")
                success += 1
                if verbose:
                    print(f"✓ Saved: {output_filename}")
                
            except Exception as e:
                error_msg = f"Failed to process: {filename} | Error: {str(e)}"
                if verbose:
                    print(f"✗ {error_msg}")
                failed += 1
        else:
            skipped += 1
            if verbose:
                print(f"Warning: Prefix '{file_prefix}' of {filename} not found in mapping table, skipped")
    
    return {
        'processed': processed,
        'success': success,
        'skipped': skipped,
        'failed': failed
    }

def print_extractcp_summary(results, output_dir):
    """Print extractcp processing summary"""
    print("\n" + "=" * 60)
    print("Extract and Rename Processing Summary")
    print("=" * 60)
    print(f"Total processed: {results['processed']} files")
    print(f"Successfully extracted: {results['success']} files")
    print(f"Skipped: {results['skipped']} files")
    print(f"Failed: {results['failed']} files")
    print(f"Results saved to: {output_dir}")

def batch_extract_chloroplast_genomes(input_dir, output_dir, mapping_file=None, verbose=False):
    """
    Batch process all GenBank files in a folder, optionally rename them
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # If renaming is needed, create temporary directory
    temp_dir = None
    if mapping_file:
        temp_dir = os.path.join(output_dir, "temp_extracted")
        os.makedirs(temp_dir, exist_ok=True)
        extract_dir = temp_dir
    else:
        extract_dir = output_dir
    
    # Statistics
    processed = 0
    success = 0
    skipped = 0
    failed = 0
    
    # Traverse all files in the input directory
    for filename in os.listdir(input_dir):
        if filename.lower().endswith((".gb", ".gbk", ".genbank")):
            input_path = os.path.join(input_dir, filename)
            processed += 1
            
            try:
                # Extract chloroplast sequence
                chloroplast_record = extract_chloroplast_sequence(input_path)
                
                if chloroplast_record:
                    # Create output filename
                    output_filename = f"{chloroplast_record.id}.fasta"
                    output_path = os.path.join(extract_dir, output_filename)
                    
                    # Save as FASTA
                    SeqIO.write(chloroplast_record, output_path, "fasta")
                    if verbose:
                        print(f"✓ Successfully extracted: {filename} -> {output_filename}")
                    success += 1
                else:
                    skipped += 1
                    if verbose:
                        print(f"Warning: No chloroplast sequence found in {filename}")
            except Exception as e:
                if verbose:
                    print(f"✗ Failed to process: {filename} | Error: {str(e)}")
                failed += 1
    
    # If renaming is needed, process renaming
    if mapping_file and temp_dir:
        if verbose:
            print("\n" + "=" * 50)
            print("Starting to rename extracted sequence files...")
        
        # Load mapping table
        mapping = load_mapping_fasta(mapping_file, verbose)
        if mapping:
            rename_results = rename_fasta_files(temp_dir, output_dir, mapping, verbose)
            
            # Clean up temporary directory
            shutil.rmtree(temp_dir)
            
            return {
                'processed': processed,
                'success': rename_results['success'],
                'skipped': skipped + rename_results['skipped'],
                'failed': failed + rename_results['failed']
            }
        else:
            print("Error: Unable to load mapping table, renaming operation cancelled")
            # Move temporary files to output directory
            for filename in os.listdir(temp_dir):
                shutil.move(os.path.join(temp_dir, filename), os.path.join(output_dir, filename))
            shutil.rmtree(temp_dir)
            
            return {
                'processed': processed,
                'success': success,
                'skipped': skipped,
                'failed': failed
            }
    else:
        return {
            'processed': processed,
            'success': success,
            'skipped': skipped,
            'failed': failed
        }

def extract_sequences_from_genbank(gb_file, output_dir, extract_type="both"):
    """
    从单个GenBank文件中提取CDS和/或蛋白质序列
    
    参数:
        gb_file: GenBank文件路径
        output_dir: 输出目录
        extract_type: 提取类型，可以是 "cds", "protein" 或 "both"
    
    返回:
        tuple: (成功提取的CDS数量, 成功提取的蛋白质数量, 登录号)
    """
    # 获取输入文件前缀（不含扩展名）作为登录号
    accession = os.path.splitext(os.path.basename(gb_file))[0]
    
    cds_records = []  # 存储所有CDS序列记录
    protein_records = []  # 存储所有蛋白质序列记录
    
    try:
        # 解析GenBank文件
        for record in SeqIO.parse(gb_file, "genbank"):
            # 遍历所有特征
            for feature in record.features:
                # 只处理CDS特征
                if feature.type == "CDS":
                    # 获取CDS标识信息
                    locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                    gene = feature.qualifiers.get("gene", [""])[0]
                    protein_id = feature.qualifiers.get("protein_id", [""])[0]
                    product = feature.qualifiers.get("product", [""])[0]
                    
                    # 创建有意义的ID
                    if gene:
                        seq_id = gene
                    elif locus_tag:
                        seq_id = locus_tag
                    elif protein_id:
                        seq_id = protein_id
                    else:
                        # 如果没有标识，使用位置信息
                        location = str(feature.location)
                        seq_id = f"CDS_{location.replace(':', '_').replace('[', '').replace(']', '')}"
                    
                    # 创建描述
                    description = f"{product} | {record.id} | {feature.location}"
                    if gene and gene != seq_id:
                        description = f"{gene}: {description}"
                    
                    # 提取CDS序列
                    if extract_type in ["cds", "both"]:
                        try:
                            cds_seq = feature.extract(record.seq)
                            cds_record = SeqRecord(
                                Seq(cds_seq),
                                id=seq_id,
                                description=description
                            )
                            cds_records.append(cds_record)
                        except Exception as e:
                            print(f"警告: 无法提取CDS序列 {seq_id}: {str(e)}")
                    
                    # 提取蛋白质序列
                    if extract_type in ["protein", "both"] and "translation" in feature.qualifiers:
                        try:
                            protein_seq = feature.qualifiers["translation"][0]
                            protein_record = SeqRecord(
                                Seq(protein_seq),
                                id=seq_id,
                                description=f"{description} | protein"
                            )
                            protein_records.append(protein_record)
                        except Exception as e:
                            print(f"警告: 无法提取蛋白质序列 {seq_id}: {str(e)}")
        
        # 保存序列到文件
        if extract_type in ["cds", "both"] and cds_records:
            cds_output_file = os.path.join(output_dir, f"{accession}_cds.fasta")
            SeqIO.write(cds_records, cds_output_file, "fasta")
            print(f"成功提取 {len(cds_records)} 条CDS序列到: {cds_output_file}")
        
        if extract_type in ["protein", "both"] and protein_records:
            protein_output_file = os.path.join(output_dir, f"{accession}_protein.fasta")
            SeqIO.write(protein_records, protein_output_file, "fasta")
            print(f"成功提取 {len(protein_records)} 条蛋白质序列到: {protein_output_file}")
        
        return len(cds_records), len(protein_records), accession
    
    except Exception as e:
        print(f"处理文件 {gb_file} 时出错: {str(e)}")
        return 0, 0, None

def batch_extract_sequences(input_dir, output_dir, extract_type="both", mapping_file=None):
    """
    批量处理目录中的所有GenBank文件
    
    参数:
        input_dir: 输入目录
        output_dir: 输出目录
        extract_type: 提取类型，可以是 "cds", "protein" 或 "both"
        mapping_file: 映射文件路径（可选）
    
    返回:
        tuple: (处理的文件数, 成功提取CDS的文件数, 成功提取蛋白质的文件数)
    """
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)
    
    # 读取映射文件（如果提供）
    rename_map = {}
    if mapping_file:
        try:
            # 读取映射文件，不假设列名
            df = pd.read_csv(mapping_file, sep='\t')
            
            # 假设第一列是原始名称，第二列是新名称
            if len(df.columns) >= 2:
                rename_map = dict(zip(df.iloc[:, 0], df.iloc[:, 1]))
                print(f"已加载 {len(rename_map)} 个重命名映射")
            else:
                print("警告: 映射文件必须至少包含两列，将跳过重命名")
        except Exception as e:
            print(f"警告: 无法读取映射文件 {mapping_file}: {str(e)}，将跳过重命名")
    
    # 统计信息
    processed_files = 0
    successful_cds_extractions = 0
    successful_protein_extractions = 0
    total_cds = 0
    total_proteins = 0
    
    # 存储所有提取的文件，用于后续重命名
    extracted_files = []
    
    # 遍历输入目录中的所有GenBank文件
    for filename in os.listdir(input_dir):
        if filename.lower().endswith(('.gb', '.gbk', '.genbank')):
            gb_file = os.path.join(input_dir, filename)
            print(f"\n处理文件: {filename}")
            cds_count, protein_count, accession = extract_sequences_from_genbank(gb_file, output_dir, extract_type)
            
            if cds_count > 0:
                successful_cds_extractions += 1
                total_cds += cds_count
                extracted_files.append(f"{accession}_cds.fasta")
            
            if protein_count > 0:
                successful_protein_extractions += 1
                total_proteins += protein_count
                extracted_files.append(f"{accession}_protein.fasta")
            
            processed_files += 1
    
    # 如果有映射文件，执行重命名
    if rename_map and extracted_files:
        rename_count = rename_extracted_files(output_dir, rename_map, extracted_files)
        print(f"\n重命名完成: 成功重命名 {rename_count} 个文件")
    
    print(f"\n处理完成!")
    print(f"共处理 {processed_files} 个文件")
    print(f"成功提取CDS的文件: {successful_cds_extractions}，共 {total_cds} 条CDS序列")
    print(f"成功提取蛋白质的文件: {successful_protein_extractions}，共 {total_proteins} 条蛋白质序列")
    
    return processed_files, successful_cds_extractions, successful_protein_extractions

def rename_extracted_files(output_dir, rename_map, extracted_files):
    """
    根据映射文件重命名已提取的文件
    
    参数:
        output_dir: 输出目录
        rename_map: 映射字典 {accession: species}
        extracted_files: 已提取的文件列表
    
    返回:
        int: 成功重命名的文件数量
    """
    rename_count = 0
    
    for filename in extracted_files:
        # 从文件名中提取登录号
        if filename.endswith('_cds.fasta'):
            accession = filename.replace('_cds.fasta', '')
            file_type = 'cds'
        elif filename.endswith('_protein.fasta'):
            accession = filename.replace('_protein.fasta', '')
            file_type = 'protein'
        else:
            continue
        
        # 检查是否有对应的映射
        if accession in rename_map:
            species = rename_map[accession]
            old_path = os.path.join(output_dir, filename)
            new_filename = f"{species}_{file_type}.fasta"
            new_path = os.path.join(output_dir, new_filename)
            
            # 重命名文件
            try:
                os.rename(old_path, new_path)
                print(f"重命名: {filename} -> {new_filename}")
                rename_count += 1
            except Exception as e:
                print(f"重命名文件 {filename} 时出错: {str(e)}")
    
    return rename_count

def simplify_sequence_id(input_path, output_path):
    """
    简化序列ID，只保留基因简称和登录号，并处理重复ID
    
    参数:
        input_path: 输入文件路径
        output_path: 输出文件路径
    
    返回:
        int: 处理的序列数量
    """
    processed_sequences = 0
    
    try:
        # 读取输入文件
        records = list(SeqIO.parse(input_path, "fasta"))
        
        # 用于跟踪基因简称的出现次数
        gene_count = defaultdict(int)
        gene_accession_map = {}
        
        # 第一遍：统计每个基因简称的出现次数
        for record in records:
            parts = record.description.split('|')
            
            if len(parts) >= 2:
                # 提取基因简称 (第一部分，空格前的部分)
                gene_part = parts[0].strip()
                gene_short = gene_part.split()[0] if ' ' in gene_part else gene_part
                
                # 提取登录号 (第二部分)
                accession = parts[1].strip()
                
                # 创建基础ID
                base_id = f"{gene_short}_{accession}"
                
                # 统计基因简称的出现次数
                gene_count[gene_short] += 1
                
                # 存储基因简称和登录号的映射
                gene_accession_map[record.id] = (gene_short, accession, base_id)
        
        # 第二遍：为重复的基因简称添加序号
        gene_occurrence = defaultdict(int)
        
        for record in records:
            if record.id in gene_accession_map:
                gene_short, accession, base_id = gene_accession_map[record.id]
                
                # 如果这个基因简称出现多次，添加序号
                if gene_count[gene_short] > 1:
                    gene_occurrence[gene_short] += 1
                    new_id = f"{gene_short}.{gene_occurrence[gene_short]}_{accession}"
                else:
                    new_id = base_id
                
                record.id = new_id
                record.description = new_id  # 清空描述，只保留ID
            
            processed_sequences += 1
        
        # 保存处理后的序列
        SeqIO.write(records, output_path, "fasta")
        print(f"成功简化 {processed_sequences} 条序列ID，保存到: {output_path}")
        
        return processed_sequences
    
    except Exception as e:
        print(f"处理文件 {input_path} 时出错: {str(e)}")
        return 0

def batch_simplify_ids(input_dir, output_dir):
    """
    批量处理目录中的所有FASTA文件
    
    参数:
        input_dir: 输入目录
        output_dir: 输出目录
    
    返回:
        int: 处理的文件数量
    """
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)
    
    processed_files = 0
    total_sequences = 0
    
    # 遍历输入目录中的所有FASTA文件
    for filename in os.listdir(input_dir):
        if filename.lower().endswith(('.fasta', '.fas', '.fa')):
            input_file = os.path.join(input_dir, filename)
            output_file = os.path.join(output_dir, filename)
            
            print(f"处理文件: {filename}")
            sequence_count = simplify_sequence_id(input_file, output_file)
            
            if sequence_count > 0:
                processed_files += 1
                total_sequences += sequence_count
    
    print(f"\n处理完成!")
    print(f"共处理 {processed_files} 个文件，{total_sequences} 条序列")
    
    return processed_files

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
    
    # renamegbfiles command
    renamegbfiles_parser = subparsers.add_parser('renamegbfiles', help=texts['renamegbfiles_help'])
    renamegbfiles_parser.add_argument('-i', '--input', required=True, help=texts['input_help'])
    renamegbfiles_parser.add_argument('-r', '--replace', required=True, help=texts['mapping_help'])
    renamegbfiles_parser.add_argument('-o', '--output', required=True, help=texts['output_help'])
    
    # extractcp command
    extractcp_parser = subparsers.add_parser('extractcp', help=texts['extractcp_help'])
    extractcp_parser.add_argument('-i', '--input', required=True, help=texts['input_help'])
    extractcp_parser.add_argument('-o', '--output', required=True, help=texts['output_help'])
    extractcp_parser.add_argument('-r', '--replace', help=texts['mapping_help'])
    
    # pipeline command
    pipeline_parser = subparsers.add_parser('pipeline', help=texts['pipeline_help'])
    pipeline_parser.add_argument('-i', '--input', required=True, help=texts['input_help'])
    pipeline_parser.add_argument('-t', '--tsv', default='genomes.tsv', help=texts['tsv_help'])
    pipeline_parser.add_argument('-o', '--output', required=True, help=texts['output_help'])
    pipeline_parser.add_argument('--delay', type=float, default=1.0, help=texts['delay_help'])
    pipeline_parser.add_argument('--retries', type=int, default=3, help=texts['retries_help'])
    
    # extract_gene command
    extract_gene_parser = subparsers.add_parser('extract_gene', help=texts['extract_gene_help'])
    extract_gene_parser.add_argument('-i', '--input', required=True, help=texts['input_help'])
    extract_gene_parser.add_argument('-o', '--output', required=True, help=texts['output_help'])
    extract_gene_parser.add_argument('-t', '--type', choices=['cds', 'protein', 'both'], 
                                   default='both', help=texts['type_help'])
    extract_gene_parser.add_argument('-r', '--rename', help=texts['mapping_help'])
    
    # simplify_id command
    simplify_id_parser = subparsers.add_parser('simplify_id', help=texts['simplify_id_help'])
    
    # simplify_id parameters (directory or file)
    simplify_group = simplify_id_parser.add_mutually_exclusive_group(required=True)
    simplify_group.add_argument('-i', '--input-dir', help=texts['input_help'])
    simplify_group.add_argument('-I', '--input-file', help=texts['input_help'])
    
    simplify_output_group = simplify_id_parser.add_mutually_exclusive_group(required=True)
    simplify_output_group.add_argument('-o', '--output-dir', help=texts['output_help'])
    simplify_output_group.add_argument('-O', '--output-file', help=texts['output_help'])
    
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
                
        elif args.command == 'renamegbfiles':
            print(f"Executing renamegbfiles command")
            print(f"Input directory: {args.input}")
            print(f"Mapping file: {args.replace}")
            print(f"Output directory: {args.output}")
            
            # Check if input directory exists
            if not os.path.isdir(args.input):
                print(f"Error: Input directory '{args.input}' does not exist")
                return
            
            # Check if mapping file exists
            if not os.path.isfile(args.replace):
                print(f"Error: Mapping file '{args.replace}' does not exist")
                return
            
            # Load mapping table
            mapping = load_mapping(args.replace, args.verbose)
            if mapping is None:
                return
            
            # Execute renaming
            results = renamegbfiles_gb_files(args.input, args.output, mapping, args.verbose)
            
            # Print summary (always show summary, regardless of verbose mode)
            print_renamegbfiles_summary(results, args.output)
        
        elif args.command == 'extractcp':
            print(f"Executing extractcp command")
            print(f"Input directory: {args.input}")
            print(f"Output directory: {args.output}")
            if args.replace:
                print(f"Mapping file: {args.replace}")
            
            # Check if input directory exists
            if not os.path.isdir(args.input):
                print(f"Error: Input directory '{args.input}' does not exist")
                return
            
            # If mapping file is provided, check if it exists
            if args.replace and not os.path.isfile(args.replace):
                print(f"Error: Mapping file '{args.replace}' does not exist")
                return
            
            # Execute extraction
            results = batch_extract_chloroplast_genomes(args.input, args.output, args.replace, args.verbose)
            
            # Print summary
            print_extractcp_summary(results, args.output)
            
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
        
        elif args.command == 'extract_gene':
            print(f"Executing extract_gene command")
            print(f"Input directory: {args.input}")
            print(f"Output directory: {args.output}")
            print(f"Extraction type: {args.type}")
            if args.rename:
                print(f"Mapping file: {args.rename}")
            
            # Check if input directory exists
            if not os.path.isdir(args.input):
                print(f"Error: Input directory '{args.input}' does not exist")
                return
            
            # If mapping file is provided, check if it exists
            if args.rename and not os.path.isfile(args.rename):
                print(f"Error: Mapping file '{args.rename}' does not exist")
                return
            
            # Execute extraction
            processed_files, successful_cds, successful_protein = batch_extract_sequences(
                args.input, args.output, args.type, args.rename
            )
            
            print(f"Extraction completed! Processed {processed_files} files")
            
        elif args.command == 'simplify_id':
            print(f"Executing simplify_id command")
            
            # Handle directory processing
            if hasattr(args, 'input_dir') and hasattr(args, 'output_dir') and args.input_dir and args.output_dir:
                print(f"Input directory: {args.input_dir}")
                print(f"Output directory: {args.output_dir}")
                
                # Check if input directory exists
                if not os.path.isdir(args.input_dir):
                    print(f"Error: Input directory '{args.input_dir}' does not exist")
                    return
                
                # Execute batch simplification
                processed_files = batch_simplify_ids(args.input_dir, args.output_dir)
                print(f"Simplification completed! Processed {processed_files} files")
            
            # Handle single file processing
            elif hasattr(args, 'input_file') and hasattr(args, 'output_file') and args.input_file and args.output_file:
                print(f"Input file: {args.input_file}")
                print(f"Output file: {args.output_file}")
                
                # Check if input file exists
                if not os.path.isfile(args.input_file):
                    print(f"Error: Input file '{args.input_file}' does not exist")
                    return
                
                # Execute single file simplification
                processed_sequences = simplify_sequence_id(args.input_file, args.output_file)
                print(f"Simplification completed! Processed {processed_sequences} sequences")
            
            else:
                print("Error: Invalid parameters for simplify_id command")
                return
    
    except KeyboardInterrupt:
        print("\nOperation interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
