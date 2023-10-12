#!/usr/bin/env python


import qiime2
from qiime2.plugins.cutadapt.methods import trim_paired

from Bio.Seq import Seq


class ADAPTER:
    def __init__(self, 
                 imported_data: qiime2.Artifact) -> None:
        self.imported_data = imported_data
        
        
        
    def trim_paired_sequences(self,
                              threads: int,
                              adapter_front: str,
                              adapter_tail: str,
                              error_rate: float = 0.1,
                              indels: bool = True,
                              times: int = 1,
                              overlap: int = 3,
                              match_read_wildcards: bool = False,
                              match_adapter_wildcards: bool = True,
                              minimum_length: int = 1,
                              discard_untrimmed: bool = True,
                              max_expected_errors: float = None,
                              max_n: float = None,
                              quality_cutoff_5end: int = 0,
                              quality_cutoff_3end: int = 0,
                              quality_base: int = 33) -> qiime2.Artifact:
        
        """ trim_paired_sequences
        
        trimming adapter sequences from the original sequences
        
        Args:
            threads (int): number of CPU cores to use
            adapter_front (str): Sequences of an adapter ligated to the front. 
                                 If multiple adapter sequences, separate sequence by comma.
            adapter_tail (str): Sequences of an adapter ligated to the tail. 
                                If multiple adapter sequences, separate sequence by comma.
            error_rate (float): Maximum allowed error rate [default: 0.1]
            indels (bool): Allow insertions or deletions of bases when matching adapters [default: True]
            times (int): Remove multiple occurrences of an adapter if it is repeated, up to `times` times. [default: 1]
            overlap (int): Require at least `overlap` bases of overlap between read and adapter for an adapter to be found.
                           [default: 3]
            match_read_wildcards (bool): Interpret IUPAC wildcards (e.g., N) in reads. [default: False]
            match_adapter_wildcards (bool): Interpret IUPAC wildcards (e.g., N) in adapters. [default: True]
            minimum_length (int): Discard reads shorter than specified value. 
                                  Note, the cutadapt default of 0 has been overridden, 
                                  because that value produces empty sequence records. [default: 1]
            discard_untrimmed (bool): Discard reads in which no adapter was found. [default: True]
            max_expected_errors (float): Discard reads that exceed maximum expected erroneous nucleotides. [default: None]
            max_n (float): Discard reads with more than COUNT N bases. 
                           If COUNT_or_FRACTION is a number between 0 and 1, 
                           it is interpreted as a fraction of the read length. [default: None]
            quality_cutoff_5end (int): Trim nucleotides with Phred score quality lower than threshold from 5 prime end. [default: 0]
            quality_cutoff_3end (int): Trim nucleotides with Phred score quality lower than threshold from 3 prime end. [default: 0]
            quality_base (int): How the Phred score is encoded (33 or 64). [default: 33]

        Returns:
            trimmed (SampleData[PairedEndSequencesWithQuality]): The resulting trimmed sequences.
        """
        
        
        front_f = adapter_front.split(',')
        adapter_f = adapter_tail.split(',')
        
        
        trimmed = trim_paired(
            demultiplexed_sequences = self.imported_data,
            cores = threads,
            adapter_f = adapter_f,
            front_f = [f'^{f}' for f in front_f],
            adapter_r = [str(Seq(f).reverse_complement()) for f in front_f],
            front_r = [f'^{str(Seq(f).reverse_complement())}' for f in adapter_f],
            error_rate = error_rate,
            indels = indels,
            times = times,
            overlap = overlap,
            match_read_wildcards = match_read_wildcards,
            match_adapter_wildcards = match_adapter_wildcards,
            minimum_length = minimum_length,
            discard_untrimmed = discard_untrimmed,
            max_expected_errors = max_expected_errors,
            max_n = max_n,
            quality_cutoff_5end = quality_cutoff_5end,
            quality_cutoff_3end = quality_cutoff_3end,
            quality_base = quality_base)
        
        return trimmed