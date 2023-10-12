#!/usr/bin/env python


import os
import subprocess

import gzip

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO

import qiime2

import pandas as pd
import biom

import skbio
from skbio import read
from skbio import TreeNode


class UTILS:
    
    def make_directory(self, dirname: str) -> None:
        
        """ make_directory
        
        make a new directory

        Args:
            dirname (str): directory path and name to make

        Returns:
            None
        """
        
        if not os.path.exists(dirname):
            os.mkdir(path=dirname)
            
        return None
    
    
    def run_the_command(self, command: list) -> None:
        
        """ run_the_command
        
        running the command provided
        
        Args:
            command (list): a list of command and its associated options
            
        Returns:
            None
        """
        
        subprocess.call(command, shell=True)
        
        
    def load_dataframe(self, 
                       table: str, 
                       **kwargs) -> pd.DataFrame:
        
        """ load_dataframe
        
        load tab or csv table file to the pandas dataframe
        
        Args:
            table (str): a table file separated with tab or comma
            sep (str): delimiter to use [default : ',']
            header (int): list of int [default : None]
            index_col (int, str, sequence of int / str, or False): columns to use as the index column
            usecols (list-like or callable): returns a subset of the columns
            prefix (str): add prefix to the column numbers when no header
            dtype (type): Type name or dict of column
            engine (str): parser engine to use, available : 'c', 'python', 'pyarrow'
            skiprows (list-like, int or callable): line numbers to skip or number of lines to skip
            verbose (bool): indicate number of NA values placed in non-numeric columns
            skip_blank_lines (bool): if True, skip over blank lines rather than interpreting as NaN values
            chunksize (int) : return TextFileReader object for iteration.
            thousands (str) : thousands separator
            encoding (str) : encoding to use for UTF when reading/writing (ex. 'utf-8')
            low_memory (bool) : if True, use the lower memory
            other options can be found from the pandas.read_csv()
            
        Returns:
            dataframe (pd.DataFrame): a pandas dataframe converted from a table provided
        """
        
        dataframe = pd.read_csv(
            table, 
            **kwargs)
        
        return dataframe
    
    
    def convert_biom(
        self, 
        biom_table: str) -> pd.DataFrame:
        
        """ load_biom
        
        loading biom formatted file
        
        Args:
            biom_table (str): a biom-formatted table to convert

        Returns:
            cvt_dataframe (pd.DataFrame): a converted table from the biom file
        """
        
        table = biom.load_table(biom_table)
        cvt_table = table.to_tsv()
        
        converted = {}
        columns = None
        values = []
        for line in cvt_table.split('\n'):
            if line.startswith('# Constructed'):
                pass
            
            else:
                if line.startswith('#OTU'):
                    columns = line.split('\t')
                
                else:
                    values.append(line.split('\t'))
                    
        converted = pd.DataFrame(values, columns = columns)
        converted = converted.set_index('#OTU ID')
        converted = converted.astype(float)
        return converted
    
    
    def get_samplenames(self, metadata: str) -> list:
        
        """ get_samplenames
        
        getting samplenames from metadata
        
        Args:
            metadata (str): path and name of metadata file

        Returns:
            samplenames (list): a list of samplenames
        """
        
        meta = self.load_dataframe(
            table=metadata, 
            sep='\t',
            header=0, 
            skiprows=lambda x:x==1)
        samplenames = meta['#SampleID'].tolist()
        
        return samplenames
    
    
    def artifact_to_dataframe(self, 
                              artifact: qiime2.Artifact) -> pd.DataFrame:
        
        """ artifact_to_dataframe
        
        converting qiime2 artifact (table format) to the dataframe
        
        Args:
            artifact (qiime2.Artifact): a qza formatted table

        Returns:
            dataframe (pandas.DataFrame): a converted dataframe
        """
        
        dataframe = artifact.table.view(qiime2.Metadata).to_dataframe()
        return dataframe
    
    
    def path2abspath(self, path: str):
        
        """ path2abspath
        
        converting relative path to absolute path
        
        Args:
            path (str): path

        Returns:
            abspath (str): absolute path
        """
        
        abspath = os.path.abspath(path)
        if not path == abspath:
            return abspath
        return path
    
    
    # export qiime2 artifact
    def export_q2_artifact(self, 
                           artifact: qiime2.Artifact,
                           output: str) -> None:
        
        """ export_q2_artifact
        
        exporting QIIME2 Artifact
        
        Args:
            artifact (qiime2.Artifact): QIIME2 artifact to export
            output (str): An output directory path and name
        
        Returns:
            None
        """
        
        artifact.export_data(output_dir=output)
        
        
    # load tree file
    def load_tree(
        self, 
        tree: str) -> skbio.TreeNode:
        
        return read(tree, format='newick', into=TreeNode)
    
    
class SEQ:
    
    # check if the file compressed with gzip
    def is_gzip(self,
                filename: str) -> bool:
        
        """ is_gzip
        
        checking if the file is compressed with gzip format
        
        Args:
            filename (str): a path and name of file

        Returns:
            (bool): if True, file is compressed by gzip
        """
        
        try:
            with open(filename, 'rb') as file:
                return file.read(2) == b'\x1f\x8b'
        except:
            return False
    
    
    # parsing fastq file
    def parse_fastq(self, 
                    fastq: str) -> tuple:
        
        """ parse_fastq
        
        parsing fastq file (available with fastq and fastq.gz)

        Raises:
            fastq (str): a file path and name of fastq

        Returns:
            read_ids (list): a list of read ids
            reads (list): a list of parsed fastq
            quals (list): a list of fastq qualities
        """
        
        read_ids = []
        reads = []
        quals = []
        if self.is_gzip(filename=fastq):
            with gzip.open(fastq, 'rt') as in_file:
                for record in SeqIO.parse(in_file, 'fastq'):
                    read_id = str(record.id)
                    read_seq = str(record.seq)
                    qual = record.letter_annotations['phred_quality']
                    
                    read_ids.append(read_id)
                    reads.append(read_seq)
                    quals.append(qual)
                    
        else:
            try:
                with open(fastq, 'r') as in_file:
                    for record in SeqIO.parse(in_file, 'fastq'):
                        read_id = str(record.id)
                        read_seq = str(record.seq)
                        qual = record.letter_annotations['phred_quality']
                        
                        read_ids.append(read_id)
                        reads.append(read_seq)
                        quals.append(qual)
                        
            except:
                raise FileNotFoundError(f'{os.path.basename(fastq)} does not exist')
            
        return (read_ids, reads, quals)
    
    
    
    # parsing fasta file
    def parse_fasta(self,
                    fasta: str) -> dict:
        
        """ parse_fasta
            
        parsing fasta file (available with fasta and fasta.gz)

        Raises:
            fasta (str): a file path and name of fasta

        Returns:
            reads (dict): a dictionary of parsed fasta
        """
        
        reads = {}
        if self.is_gzip(filename=fasta):
            with gzip.open(fasta, 'rt') as in_file:
                for record in SeqIO.parse(in_file, 'fasta'):
                    reads[record.id] = [str(record.seq)]
        else:
            try:
                for record in SeqIO.parse(fasta, 'fasta'):
                    reads[record.id] = [str(record.seq)]
                    
            except:
                raise FileNotFoundError(f'{os.path.basename(fasta)} does not exist')
            
        return reads