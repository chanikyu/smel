#!/usr/bin/env python


import os
from tqdm import tqdm

import numpy as np
import statsmodels.api as sm

import qiime2
from qiime2.plugins.dada2.actions import denoise_paired
from qiime2.plugins.dada2.actions import denoise_ccs


class DENOISE:
    
    def __init__(self, 
                 demuxed_data: qiime2.Artifact) -> None:
        self.demuxed_data = demuxed_data
        
        
    # calculate truncation length
    # calculate mean quality scores at each positions  
    def mean_quality_scores_at_each_position(self,
                                             qscores: list) -> list:
        
        """ __mean_quality_scores_at_each_position
        
        calculate mean quality scores at each position
        
        Args:
            qscores (list): quality scores of reads

        Returns:
            mean_results (list): mean quality score list
        """
        
        max_len = np.max([len(qscore) for qscore in qscores])
        
        sums = [0] * max_len
        counts = [0] * max_len
        
        for sublist in tqdm(qscores, 
                            desc='Calculated Mean Quality',
                            ascii=True,
                            leave=False, 
                            position=0):
            for i in range(len(sublist)):
                sums[i] += sublist[i]
                counts[i] += 1
                    
        # calculate mean quality scores
        mean_results = []
        for vals in zip(sums, counts):
            mean_results.append(vals[0] / vals[1])
            
        return mean_results
    
    
    # LOESS and Cook's distance
    def loess_to_cooks_distance(self, 
                                  means : list, 
                                  frac : float = 0.2) -> tuple:
        
        """ __loess_to_cooks_distance
        
        fits LOESS model and calculates Cook's distance for the detection of outliers.
        the outliers refer to the trimming positions.
        
        Args:
            means (list): list of mean quality scores at each position
            frac (float): fraction of the data used when estimating each y-value [values between 0 - 1]
            
        Returns:
            loess (numpy.ndarray): the result of LOESS fit
            cooks_distance (numpy.ndarray): the result of Cook's distance
        
        """
        
        loess = sm.nonparametric.lowess(means, range(len(means)), frac=float(frac))
        # print(type(loess))
    
        infl = sm.OLS(loess[:, 1], sm.add_constant(means)).fit().get_influence()
        cooks_distance, _ = infl.cooks_distance
        # print(type(cooks_distance))
        
        return (loess, cooks_distance)
    
    
    # determine start and end position
    def determine_start_end_position(self, 
                                     cooks_distance: np.ndarray, 
                                     means : list,
                                     threshold : float = 0.05):
        
        """ determine_start_end_position
        
        determines and suggests start and end positions to trim.
        
        Args:
            cooks_distance (numpy.ndarray): Cook's distance
            means (list): list of mean qualities
            threshold (float): the cutoff value to determine outliers from Cook's distance [default : 0.05]
            
        Returns:
            * trim_start (int): start trimming position
            * trim_end (int): end trimming position
        """
    
        outliers = np.where(cooks_distance > float(threshold))[0]
        
        if len(outliers) == 0:
            trim_start = 0
            trim_end = len(means)
        else:
            front_outliers = []
            tail_outliers = []

            for out in outliers:
                if out < 150:
                    front_outliers.append(out)
                else:
                    tail_outliers.append(out)

            if len(front_outliers) != 0:
                trim_start = np.max(front_outliers)
            else:
                trim_start = 0

            if len(tail_outliers) != 0:
                trim_end = np.min(tail_outliers)
            else:
                trim_end = len(means)

        return (trim_start, trim_end)
    
    
    # DADA2 - paired end
    def runDADA2_paired(
        self,
        trunc_len_f: int,
        trunc_len_r: int,
        trim_left_f: int = 0,
        trim_left_r: int = 0,
        max_ee_f: float = 2.0,
        max_ee_r: float = 2.0,
        trunc_q: int = 2,
        min_overlap: int = 12,
        pooling_method: str = 'independent',
        chimera_method: str = 'consensus',
        min_fold_parent_over_abundance: float = 1.0,
        allow_one_off: bool = False,
        n_threads: int = 1,
        n_reads_learn: int = 1000000,
        hashed_feature_ids: bool = True) -> tuple:
        
        """ runDADA2_paired
        
        running DADA2 with paired end sequences
        
        Args:
            trunc_len_f (int): Position at which forward read sequences should be truncated due to decrease in quality. 
                               This truncates the 3' end of the of the input sequences, 
                               which will be the bases that were sequenced in the last cycles. 
                               Reads that are shorter than this value will be discarded. 
                               After this parameter is applied there must still be at least a 12 nucleotide
                               overlap between the forward and reverse reads. If 0 is provided, 
                               no truncation or length filtering will be performed
            trunc_len_r (int): Position at which reverse read sequences should be truncated due to decrease in quality. 
                               This truncates the 3' end of the of the input sequences, 
                               which will be the bases that were sequenced in the last cycles. 
                               Reads that are shorter than this value will be discarded. 
                               After this parameter is applied there must still be at least a 12 nucleotide
                               overlap between the forward and reverse reads. If 0 is provided, 
                               no truncation or length filtering will be performed
            trim_left_f (int): Position at which forward read sequences should be trimmed due to low quality. 
                               This trims the 5' end of the input sequences, 
                               which will be the bases that were sequenced in the first cycles. [default = 0]
            trim_left_r (int): Position at which reverse read sequences should be trimmed due to low quality. 
                               This trims the 5' end of the input sequences, 
                               which will be the bases that were sequenced in the first cycles. [default = 0]
            max_ee_f (float): Forward reads with number of expected errors higher than this value will be discarded. [default = 2.0]
            max_ee_r (float): Reverse reads with number of expected errors higher than this value will be discarded. [default = 2.0]
            trunc_q (int): Reads are truncated at the first instance of a quality score less than or equal to this value. 
                           If the resulting read is then shorter than `trunc_len_f` or `trunc_len_r` 
                           (depending on the direction of the read) it is discarded. [default = 2]
            min_overlap (int): The minimum length of the overlap required for merging the forward and reverse reads. [default = 12]
            pooling_method (str): The method used to pool samples for denoising. [default = independent]
                                    - "independent": Samples are denoised indpendently. 
                                    - "pseudo": The pseudo-pooling method is used to approximate pooling of samples. 
                                                In short, samples are denoised independently once, 
                                                ASVs detected in at least 2 samples are recorded, 
                                                and samples are denoised independently a second time, but this time
                                                with prior knowledge of the recorded ASVs and thus higher sensitivity to those ASVs.
            chimera_method (str): The method used to remove chimeras. [default = consensus]
                                    - "none": No chimera removal is performed. 
                                    - "pooled": All reads are pooled prior to chimera detection.
                                    - "consensus": Chimeras are detected in samples individually, 
                                                   and sequences found chimeric in a sufficient fraction of samples are removed.
            min_fold_parent_over_abundance (float): The minimum abundance of potential parents of a sequence being tested as chimeric,
                                                    expressed as a fold-change versus the abundance of the sequence being tested. 
                                                    Values should be greater than or equal to 1
                                                    (i.e. parents should be more abundant than the sequence being tested).
                                                    This parameter has no effect if chimera_method is "none". [default = 1.0]
            allow_one_off (bool): Bimeras that are one-off from exact are also identified if the
                                  `allow_one_off` argument is TrueIf True, a sequence will be identified
                                  as bimera if it is one mismatch or indel away from an exact bimera. [default = False]
            n_threads (int): The number of threads to use for multithreaded processing. 
                             If 0 is provided, all available cores will be used. [default = 1]
            n_reads_learn (int): The number of reads to use when training the error model. 
                                 Smaller numbers will result in a shorter run time but a less reliable error model. [default = 1000000]
            hashed_feature_ids (bool): If true, the feature ids in the resulting table will be presented as
                                       hashes of the sequences defining each feature. The hash will always be
                                       the same for the same sequence so this allows feature tables to be
                                       merged across runs of this method. You should only merge tables if the
                                       exact same parameters are used for each run. [default = True]

        Returns:
            table (FeatureTable[Frequency]): the resulting feature table.
            representative_sequences (FeatureData[Sequence]): The resulting feature sequences. 
                                                              Each feature in the feature table will
                                                              be represented by exactly one sequence, 
                                                              and these sequences will be the
                                                              joined paired-end sequences.
            denoising_stats (SampleData[DADA2Stats]) : Statistics of DADA2 denoising
        """
        
        table, rep_seqs, stats = denoise_paired(
            demultiplexed_seqs = self.demuxed_data,
            trunc_len_f = trunc_len_f,
            trunc_len_r = trunc_len_r,
            trim_left_f = trim_left_f,
            trim_left_r = trim_left_r,
            max_ee_f = max_ee_f,
            max_ee_r = max_ee_r,
            trunc_q = trunc_q,
            min_overlap = min_overlap,
            pooling_method = pooling_method,
            chimera_method = chimera_method,
            min_fold_parent_over_abundance = min_fold_parent_over_abundance,
            allow_one_off = allow_one_off,
            n_threads = n_threads,
            n_reads_learn = n_reads_learn,
            hashed_feature_ids = hashed_feature_ids
        )
        
        return (table, rep_seqs, stats)
    
    
    # DADA2 - CCS
    def runDADA2_CCS(
        self,
        front: str,
        adapter: str,
        max_mismatch: int = 2,
        indels: bool = False,
        trunc_len: int = 0,
        trim_left: int = 0,
        max_ee: float = 2.0,
        trunc_q: int = 2,
        min_len: int = 20,
        max_len: int = 0,
        pooling_method: str = 'independent',
        chimera_method: str = 'consensus',
        min_fold_parent_over_abundance: float = 3.5,
        allow_one_off: bool = False,
        n_threads: int = 1,
        n_reads_learn: int = 1000000,
        hashed_feature_ids: bool = True) -> tuple:
        
        """ runDADA2_CCS
        
        running DADA2 with CCS sequences
        
        Args:
            front (str): Sequence of an adapter ligated to the 5' end. 
                         The adapter and any preceding bases are trimmed. 
                         Can contain IUPAC ambiguous nucleotide codes. 
                         Note, primer direction is 5' to 3'. 
                         Primers are removed before trim and filter step. 
                         Reads that do not contain the primer are discarded. 
                         Each read is re-oriented if the reverse complement of the read is a better match to the provided primer sequence. 
                         This is recommended for PacBio CCS reads, which come in a random mix of forward and reverse-complement orientations.
            adapter (str): Sequence of an adapter ligated to the 3' end. 
                           The adapter and any preceding bases are trimmed. 
                           Can contain IUPAC ambiguous nucleotide codes. 
                           Note, primer direction is 5' to 3'. 
                           Primers are removed before trim and filter step. 
                           Reads that do not contain the primer are discarded.
            max_mismatch (int): The number of mismatches to tolerate when matching reads to primer sequences [default = 2]
            indels (bool): Allow insertions or deletions of bases when matching adapters. 
                           Note that primer matching can be significantly slower, currently about 4x slower [default = False]
            trunc_len (int): Position at which sequences should be truncated due to decrease in quality. 
                             This truncates the 3' end of the of the input sequences, 
                             which will be the bases that were sequenced in the last cycles. 
                             Reads that are shorter than this value will be discarded. 
                             If 0 is provided, no truncation or length filtering will be performed. 
                             Note: Since Pacbio CCS sequences were normally with very high quality scores, 
                             there is no need to truncate the Pacbio CCS sequences. [default = 0]
            trim_left (int): Position at which sequences should be trimmed due to low quality. 
                             This trims the 5' end of the of the input sequences, 
                             which will be the bases that were sequenced in the first cycles. [default = 0]
            max_ee (float): Reads with number of expected errors higher than this value will be discarded. [default = 2.0]
            trunc_q (int): Reads are truncated at the first instance of a quality score less than or equal to this value. 
                           If the resulting read is then shorter than `trunc_len`, it is discarded. [default = 2]
            min_len (int): Remove reads with length less than minLen. 
                           minLen is enforced after trimming and truncation. 
                           For 16S Pacbio CCS, suggest 1000. [default = 20]
            max_len (int): Remove reads prior to trimming or truncation which are longer than this value. 
                           If 0 is provided no reads will be removed based on length. 
                           For 16S Pacbio CCS, suggest 1600. [default = 0]
            pooling_method (str): The method used to pool samples for denoising. [default = independent]
                                    - "independent": Samples are denoised indpendently. 
                                    - "pseudo": The pseudo-pooling method is used to approximate pooling of samples. 
                                    In short, samples are denoised independently once, 
                                    ASVs detected in at least 2 samples are recorded, 
                                    and samples are denoised independently a second time, 
                                    but this time with prior knowledge of the recorded ASVs and thus higher sensitivity to those ASVs.
            chimera_method (str): The method used to remove chimeras. [default = consensus]
                                    - "none": No chimera removal is performed. 
                                    - "pooled": All reads are pooled prior to chimera detection.
                                    - "consensus": Chimeras are detected in samples individually, 
                                                   and sequences found chimeric in a sufficient fraction of samples are removed.
            min_fold_parent_over_abundance (float): The minimum abundance of potential parents of a sequence being tested as chimeric, 
                                                    expressed as a fold-change versus the abundance of the sequence being tested. 
                                                    Values should be greater than or equal to 1 
                                                    (i.e. parents should be more abundant than the sequence being tested).
                                                    Suggest 3.5. This parameter has no effect if chimera_method is "none". [default = 3.5]
            allow_one_off (bool): Bimeras that are one-off from exact are also identified if the `allow_one_off` argument is True. 
                                  If True, a sequence will be identified as bimera if it is one mismatch or indel away from an exact bimera. [default = False]
            n_threads (int): The number of threads to use for multithreaded processing. 
                             If 0 is provided, all available cores will be used. [default = 1]
            n_reads_learn (int): The number of reads to use when training the error model. 
                                 Smaller numbers will result in a shorter run time but a less reliable error model. [default = 1000000]
            hashed_feature_ids (bool): If true, the feature ids in the resulting table will be presented as hashes of the sequences defining each feature. 
                                       The hash will always be the same for the same sequence so this allows feature tables to be
                                       merged across runs of this method. You should only merge tables 
                                       if the exact same parameters are used for each run. [default = True]

        Returns:
            table (FeatureTable[Frequency]): the resulting feature table.
            representative_sequences (FeatureData[Sequence]): The resulting feature sequences. 
                                                              Each feature in the feature table will
                                                              be represented by exactly one sequence, 
                                                              and these sequences will be the
                                                              joined paired-end sequences.
            denoising_stats (SampleData[DADA2Stats]) : Statistics of DADA2 denoising
        """
        
        table, rep_seqs, stats = denoise_ccs(
            demultiplexed_seqs = self.demuxed_data,
            front = front,
            adapter = adapter,
            max_mismatch = max_mismatch,
            indels = indels,
            trunc_len = trunc_len,
            trim_left = trim_left,
            max_ee = max_ee,
            trunc_q = trunc_q,
            min_len = min_len,
            max_len = max_len,
            pooling_method = pooling_method,
            chimera_method = chimera_method,
            min_fold_parent_over_abundance = min_fold_parent_over_abundance,
            allow_one_off = allow_one_off,
            n_threads = n_threads,
            n_reads_learn = n_reads_learn,
            hashed_feature_ids = hashed_feature_ids
        )
        
        return (table, rep_seqs, stats)