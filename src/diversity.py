#!/usr/bin/env python


import os

import pandas as pd
import numpy as np

import qiime2
from qiime2.plugins.diversity.pipelines import core_metrics_phylogenetic

import conorm

from skbio.stats.composition import clr
from sklearn.preprocessing import StandardScaler

from skbio.diversity import alpha_diversity
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa

from src import UTILS
from src import SEQ


class DIVERSITY:
    def __init__(self, biom: str) -> None:
        self.util = UTILS()
        self.biom = biom
    
    
    # calculate the rarefaction depth
    def calculate_rarefaction_depth(
        self, 
        depth: int = 1000, 
        fail_ratio: float = 0.3) -> tuple:
        
        """ calculate_rarefaction_depth
        
        calculating rarefaction depth from the feature table
        
        Args:
            depth (int): a criterion of sampling depth
            fail_ratio (float): a criterion how many samples are failed

        Returns:
            depth (int): a calculated sampling depth
            the_ratio (float): the ratio of number of samples failed
        """
        
        cvt_biom = self.util.convert_biom(biom_table=self.biom)
        cvt_biom = cvt_biom.astype(float)
        cvt_biom = cvt_biom.astype(int)
        
        total_reads = cvt_biom.sum(axis=0)
        num_samples = len(total_reads)
        min_num_reads = np.min(total_reads)
        
        failed_samples = len(total_reads[total_reads <= depth])
        the_ratio = (failed_samples / num_samples)
        
        if the_ratio >= fail_ratio:
            print(ValueError(f'    \033[91mToo many samples are not passed the threshold\033[0m'))
            print(ValueError(f'    \033[91mSkip this step\033[0m\n'))
            return (0, 0)
        
        else:
            if min_num_reads <= depth:
                depth = 500
                    
            else:
                depth = min_num_reads
            
        return (depth, the_ratio)
    
    
    # alpha & beta diversity with rarefaction
    def diversity_with_rarefaction(
        self, 
        feature_table: qiime2.Artifact,
        rooted_tree: qiime2.Artifact,
        metadata: qiime2.Artifact,
        depth: int,
        replacement: bool = False,
        threads: 'int or str' = 1) -> tuple:
        
        
        """ diversity_with_rarefaction
        
        analyze alpha and beta diversity with rarefaction
        
        Args:
            feature_table (FeatureTable[Frequency]): The feature table containing the samples over which diversity metrics should be computed.
            rooted_tree (Phylogeny[Rooted]): Phylogenetic tree containing tip identifiers that correspond to the feature identifiers in the table. 
                                             This tree can contain tip ids that are not present in the table, 
                                             but all feature ids in the table must be present in this tree.
            metadata (Metadata): The sample metadata to use in the emperor plots
            depth (int): The total frequency that each sample should be rarefied to prior to computing diversity metrics.
            replacement (bool): Rarefy with replacement by sampling from the multinomial distribution instead of rarefying without replacement. 
                                [default = False]
            threads (int / str): The number of concurrent jobs or CPU threads to use in performing this calculation. 
                                 Individual methods will create jobs/threads as implemented in q2-diversity-lib dependencies. 
                                 May not exceed the number of available physical cores. 
                                 If n_jobs_or_threads = 'auto', one thread/job will be created for each identified CPU core on the host.
                                 [default = 1]

        Returns:
            _type_: _description_
        """
        
        
        rarefied_table, faith_pd_vector, observed_features_vector, shannon_vector, \
            evenness_vector, unweighted_unifrac_distance_matrix, weighted_unifrac_distance_matrix, \
                jaccard_distance_matrix, bray_curtis_distance_matrix, unweighted_unifrac_pcoa_results, \
                    weighted_unifrac_pcoa_results, jaccard_pcoa_results, bray_curtis_pcoa_results, \
                        unweighted_unifrac_emperor, weighted_unifrac_emperor, jaccard_emperor, bray_curtis_emperor = core_metrics_phylogenetic(
                            table = feature_table,
                            phylogeny = rooted_tree,
                            sampling_depth = depth,
                            metadata = metadata,
                            with_replacement = replacement,
                            n_jobs_or_threads = threads
                        )
                        
        return (rarefied_table, faith_pd_vector, observed_features_vector, shannon_vector, 
                evenness_vector, unweighted_unifrac_distance_matrix, weighted_unifrac_distance_matrix,
                jaccard_distance_matrix, bray_curtis_distance_matrix, unweighted_unifrac_pcoa_results,
                weighted_unifrac_pcoa_results, jaccard_pcoa_results, bray_curtis_pcoa_results, 
                unweighted_unifrac_emperor, weighted_unifrac_emperor, jaccard_emperor, bray_curtis_emperor)
    
    
    # normalization
    def tss_norm(self) -> pd.DataFrame:
        
        """ tss_norm
        
        total sum of squares normalization

        Returns:
            normalized_biom (pd.DataFrame): normalized table via total sum of squares
        """
        
        cvt_biom = self.util.convert_biom(biom_table=self.biom)
        cvt_biom = cvt_biom.astype(float)
        cvt_biom = cvt_biom.astype(int)
        
        total_sum_of_squares = np.sum(cvt_biom ** 2, axis=0)
        normalized_biom = cvt_biom / np.sqrt(total_sum_of_squares)
        
        return normalized_biom
    
    
    def clr_norm(self) -> pd.DataFrame:
        
        """ clr_norm
        
        clr normalization

        Returns:
            normalized_biom (pd.DataFrame): normalized table via clr
        """
        
        cvt_biom = self.util.convert_biom(biom_table=self.biom)
        cvt_biom = cvt_biom.astype(float)
        cvt_biom = cvt_biom.astype(int)
        
        columns = cvt_biom.columns
        index = cvt_biom.index
        
        normalized_biom = clr(cvt_biom)
        normalized_biom = pd.DataFrame(
            normalized_biom, 
            index=index, 
            columns=columns)
        
        return normalized_biom
    
    
    def tmm_norm(
        self, 
        trim_lfc: float = 0.3, 
        trim_mag: float = 0.05, 
        index_ref: 'float or str' = None) -> tuple:
        
        
        """ tmm_norm
        
        Trimmed Mean of M-values (TMM) normalization
        
        Args:
            trim_lfc (float): Quantile cutoff for M_g (logfoldchanges). [default : 0.3]
            trim_mag (float): Quantile cutoff for A_g (log magnitude). [default : 0.05]
            index_ref (float or str): Reference index or column name to use as reference in the TMM algorithm. [default : None]

        Returns:
            normalized_biom (pd.DataFrame): normalized table via TMM
            norm_factors (pd.DataFrame): the norm factors
        """
        
        cvt_biom = self.util.convert_biom(biom_table=self.biom)
        cvt_biom = cvt_biom.astype(float)
        cvt_biom = cvt_biom.astype(int)
        
        normalized_biom, norm_factors = conorm.tmm(
            data=cvt_biom,
            trim_lfc=trim_lfc,
            trim_mag=trim_mag,
            index_ref=index_ref,
            return_norm_factors=True
        )
        
        return (normalized_biom, norm_factors)
    
    
    def __getting_read_length__(
        self, 
        fasta: str) -> pd.DataFrame:
        
        
        """ __getting_read_length__
        
        getting length of each read
        
        Args:
            fasta (str): a fasta file path and name

        Returns:
            read_length_df (pd.DataFrame): a dataframe with calculated read length
        """
        
        seq = SEQ()
        reads = seq.parse_fasta(fasta=fasta)
        
        read_ids = []
        read_length = []
        for read_id, read_seq in reads.items():
            read_ids.append(read_id)
            read_length.append(len(read_seq))
            
        read_length_df = pd.DataFrame(read_length).T
        read_length_df.columns = read_ids
        read_length_df = read_length_df.T
        read_length_df.columns = ['length']
        
        return read_length_df
    
    
    def getmm_norm(
        self, 
        fasta: str) -> pd.DataFrame:
        
        """ getmm_norm
        
        getmm normalization
        
        Args:
            fasta (str): a fasta file path and name

        Returns:
            normalized_biom (pd.DataFrame): a tmm dataframe corrected by read length
        """
        
        read_length = self.__getting_read_length__(fasta=fasta)
        
        cvt_biom = self.util.convert_biom(biom_table=self.biom)
        cvt_biom = cvt_biom.astype(float)
        cvt_biom = cvt_biom.astype(int)
        
        table_with_length = pd.concat([cvt_biom, read_length], 
                                      axis=1)
        
        normalized_biom = pd.DataFrame()
        for col in table_with_length.columns:
            if not col == 'length':
                normed = table_with_length[col] / table_with_length['length']
                normalized_biom = pd.concat([normalized_biom, normed], 
                                            axis=1)
        
        normalized_biom.columns = table_with_length.columns[:-1]
        
        return normalized_biom
    
    
    # standardization
    def standardization(self) -> pd.DataFrame:
        
        
        """ standardization
        
        standardization

        Returns:
            standardized_biom (pd.DataFrame): a standardized dataframe
        """
        
        cvt_biom = self.util.convert_biom(biom_table=self.biom)
        cvt_biom = cvt_biom.astype(float)
        cvt_biom = cvt_biom.astype(int)
        
        scaler = StandardScaler()
        standardized = scaler.fit_transform(cvt_biom)
        
        standardized_biom = pd.DataFrame(standardized, 
                                         columns=cvt_biom.columns, 
                                         index=cvt_biom.index)
        
        return standardized_biom
    
    
    # alpha diversity
    def calculate_alpha_diversity(
        self, 
        normed_table: pd.DataFrame, 
        rooted_tree: str) -> pd.DataFrame:
        
        
        """ calculate_alpha_diversity
        
        calculate the alpha diversity using various metrics
        
        Args:
            normed_table (pd.DataFrame): a normalized table
            rooted_tree (str): a tree file path and name

        Returns:
            alpha_results (pd.DataFrame): a result of alpha diversity analysis
        """
        
        
        samplenames = list(normed_table.columns)
        asv_ids = list(normed_table.index)
        
        indices = ['chao1', 'faith_pd', 'goods_coverage', 
                   'mcintosh_e', 'observed_otus', 'pielou_e',
                   'shannon', 'simpson']
        
        alpha_results = pd.DataFrame()
        t_normed_table = normed_table.T
        
        for idx in indices:
            if not idx == 'faith_pd':
                alpha_table = alpha_diversity(
                    metric=idx, 
                    counts=t_normed_table, 
                    ids=samplenames)
                alpha_df = pd.DataFrame(alpha_table, columns=[idx])
                alpha_results = pd.concat([alpha_results, alpha_df], 
                                          axis=1)
                
            else:
                tree = self.util.load_tree(tree=rooted_tree)
                alpha_table = alpha_diversity(
                    metric=idx, 
                    counts=t_normed_table, 
                    ids=samplenames,
                    otu_ids=asv_ids,
                    tree=tree)
                alpha_df = pd.DataFrame(alpha_table, columns=[idx])
                alpha_results = pd.concat([alpha_results, alpha_df], 
                                          axis=1)
                
        return alpha_results
                
    
    # beta diversity
    def calculate_beta_diversity(
        self, 
        normed_table: pd.DataFrame,
        rooted_tree: str) -> tuple:
        
        metrics = ['unweighted_unifrac', 'weighted_unifrac', 
                   'bray_curtis', 'jaccard']
        
        
        for metric in metrics:
            if 'unifrac' in metric:
                tree = self.util.load_tree(tree=rooted_tree)
                dist_matrix = beta_diversity(
                    metric=metric, 
                    counts=normed_table.T,
                    ids=list(normed_table.columns), 
                    tree=tree,
                    otu_ids=list(normed_table.index))
                dist_matrix_df = dist_matrix.to_data_frame()
                
                pcoa_result = pcoa(dist_matrix_df)
                pcoa_result_values = pcoa_result.samples
                pcoa_result_values.index = list(dist_matrix_df.columns)
                pcoa_proportion_explained = pcoa_result.proportion_explained
                pcoa_eigen_values = pcoa_result.eigvals
                
                return (pcoa_result_values, pcoa_proportion_explained, pcoa_eigen_values)
                
            else:
                dist_matrix = beta_diversity(
                    metric=metric, 
                    counts=normed_table.T,
                    ids=list(normed_table.columns))
                dist_matrix_df = dist_matrix.to_data_frame()
                
                pcoa_result = pcoa(dist_matrix_df)
                pcoa_result_values = pcoa_result.samples
                pcoa_result_values.index = list(dist_matrix_df.columns)
                pcoa_proportion_explained = pcoa_result.proportion_explained
                pcoa_eigen_values = pcoa_result.eigvals
                
                return (pcoa_result_values, pcoa_proportion_explained, pcoa_eigen_values)