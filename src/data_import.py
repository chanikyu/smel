#!/usr/bin/env python


import os
import sys

from src import UTILS

import qiime2


class IMPORT:
    
    def __init__(self, dirname: str) -> None:
        self.ut = UTILS()
        self.dirname = dirname
        self.paired = self.check_paired_end()
        
        
    def check_paired_end(self) -> bool:
        
        """ check_paired_end
        
        checking if the sequencing files are paired end or not
        
        Args:
            dirname (str): a directory name that fastq files exist

        Returns:
            (bool): if paired end, returns True
        """
        
        fwd_files = sorted(
            list(
                filter(
                    lambda x:(x.endswith('_1.fastq.gz') or x.endswith('_R1.fastq.gz') or x.endswith('_R1_001.fastq.gz')), 
                    os.listdir(path=self.dirname))))
        if len(fwd_files) != 0:
            return True
        return False
    
    
    def check_file_exists(
        self, 
        samplename: str, 
        filename: str) -> bool:
        
        """ check_file_exists
        
        checking if the files exist
        
        Args:
            samplename (str): samplename
            filename (str): filename

        Returns:
            (bool): if samplename in filename, returns True
        """
        
        if samplename in filename:
            return True
        else:
            return False
        
    
    def create_manifest(
        self, 
        metadata: str,
        output: str) -> None:
        
        
        """ create_manifest
        
        creating manifest file
        
        Args:
            dirname (str): a directory path and name including rawdata
            metadata (str): a metadata file path and name
            output (str): an output manifest file path and name
        
        """
        
        samplenames = self.ut.get_samplenames(metadata=metadata)
        
        with open(output, 'a') as manifest:
            
            manifest.write('sample-id,absolute-filepath,direction\n')
            
            # paired end
            if self.paired:
                
                fwd_files = sorted(list(filter(
                    lambda x:x.endswith('_1.fastq.gz') or x.endswith('_R1.fastq.gz'),
                    os.listdir(path=self.dirname))))
                rev_files = sorted(list(filter(
                    lambda x:x.endswith('_2.fastq.gz') or x.endswith('_R2.fastq.gz'),
                    os.listdir(path=self.dirname))))
                
                for samplename in samplenames:
                    
                    exist = []
                    
                    if (f'{samplename}_1.fastq.gz' in fwd_files) or (f'{samplename}_R1.fastq.gz'):
                        exist.append(True)
                    else:
                        exist.append(False)
                        
                    if (f'{samplename}_2.fastq.gz' in rev_files) or (f'{samplename}_R2.fastq.gz'):
                        exist.append(True)
                    else:
                        exist.append(False)

                    if all(exist):
                        
                        if f'{samplename}_1.fastq.gz' in fwd_files:
                            fwd = f'{samplename}_1.fastq.gz'
                        elif f'{samplename}_R1.fastq.gz' in fwd_files:
                            fwd = f'{samplename}_R1.fastq.gz'
                            
                        if f'{samplename}_2.fastq.gz' in rev_files:
                            rev = f'{samplename}_2.fastq.gz'
                        elif f'{samplename}_R2.fastq.gz' in rev_files:
                            rev = f'{samplename}_R2.fastq.gz'
                            
                        dirpath = self.ut.path2abspath(path=self.dirname)
                            
                        fwd_abs_path = os.path.join(dirpath, fwd)
                        rev_abs_path = os.path.join(dirpath, rev)
                        
                        manifest.write(f'{samplename},{fwd_abs_path},forward\n')
                        manifest.write(f'{samplename},{rev_abs_path},reverse\n')
                    
                    # file missing
                    else:
                        print(f'\033[91mFile Missing : {samplename}\033[0m')
                        is_continue = str(input('Continue to the next process? [y/n]'))
                        if is_continue.lower() == 'y':
                            pass
                        else:
                            sys.exit(1)
                        
            # single end
            else:
                files = sorted(list(filter(
                    lambda x:x.endswith('.fastq.gz'),
                    os.listdir(path=self.dirname))))
                
                for samplename in samplenames:
                    
                    if (f'{samplename}.fastq.gz' in files):
                        dirpath = self.ut.path2abspath(path=self.dirname)
                        abs_path = os.path.join(dirpath, f'{samplename}.fastq.gz')
                        manifest.write(f'{samplename},{abs_path},forward\n')
                    
                    else:
                        print(f'\033[91mFile Missing : {samplename}\033[0m')
                        is_continue = str(input('Continue to the next process? [y/n]'))
                        if is_continue.lower() == 'y':
                            pass
                        else:
                            sys.exit(1)
                        
        return None
                        
                        
    # load sequence artifact
    def load_seq_artifact(
        self, 
        manifest: str) -> qiime2.sdk.result.Artifact:
        
        """ load_artifact
        
        loading sequence artifact via manifest
        
        Args:
            manifest (str): a manifest file to load files

        Returns:
            imported (qiime2.sdk.result.Artifact): imported data with QIIME2 Artifact format
        """
        
        from qiime2 import Artifact
        
        if self.paired:
            import_type = 'SampleData[PairedEndSequencesWithQuality]'
            view_type = 'PairedEndFastqManifestPhred33'
        else:
            import_type = 'SampleData[SequencesWithQuality]'
            view_type = 'SingleEndFastqManifestPhred33'
        
        imported = Artifact.import_data(
            type=import_type, 
            view=manifest, 
            view_type=view_type)
        
        return imported
    
    
    # load metadata artifact
    def load_metadata_artifact(
        self, 
        metadata: str) -> qiime2.sdk.result.Artifact:
        
        """ load_metadata_artifact
        
        loading metadata artifact
        
        Args:
            metadata (str): a metadata file

        Returns:
            (qiime2.sdk.result.Artifact) : imported metadata with QIIME2 Artifact format
        """
        
        from qiime2 import Metadata
        
        return Metadata.load(metadata)