#!/usr/bin/env python

import os

import qiime2
import qiime2.plugins.feature_classifier.actions as feature_classifier
from qiime2.plugins.taxa.methods import filter_table
from qiime2.plugins.taxa.methods import filter_seqs


class CLASSIFICATION:
    def __init__(self,
                 reads: qiime2.Artifact,
                 classifier: str):
        
        self.reads = reads
        self.classifier = classifier
        
        
    # load classifier
    def load_classifier(self) -> qiime2.Artifact:
        
        """ load_classifier
        
        load classifier

        Returns:
            classifier_file (qiime2.Artifact): a taxonomic classifier
        """
        
        classifier_filename = ''
        
        if self.classifier.lower() == 'greengenes2':
            classifier_filename = os.path.join(os.environ['CONDA_PREFIX'], 'classifiers', 'gg_2022_10_backbone_full_length.nb.qza')
        elif self.classifier.lower() == 'silva':
            classifier_filename = os.path.join(os.environ['CONDA_PREFIX'], 'classifiers', 'silva-138-99-nb-classifier.qza')
            
        # print(classifier_filename)
        
        classifier_file = qiime2.Artifact.load(classifier_filename)
        return classifier_file
    
    
    # taxonomy classification
    def taxonomy_classification(
        self,
        reads_per_batch: 'int or str' = 'auto',
        n_jobs: int = 1,
        pre_dispatch: str = '2*n_jobs',
        confidence: 'float or str' = 0.7,
        read_orientation: str = 'auto') -> qiime2.Artifact:
        
        """ taxonomy_classification
        
        taxonomic classification of reads
        
        Args:
            reads_per_batch (int | str): Number of reads to process in each batch. 
                                         If "auto", this parameter is autoscaled to min (number of query sequences / n_jobs, 20000).
                                         [default = auto]
            n_jobs (int): The maximum number of concurrently worker processes. 
                          If -1 all CPUs are used. 
                          If 1 is given, no parallel computing code is used at all, which is useful for debugging. 
                          For n_jobs below -1, (n_cpus + 1 + n_jobs) are used. 
                          Thus for n_jobs = -2, all CPUs but one are used.
                          [default = 1]
            pre_dispatch (str): "all" or expression, as in "3*n_jobs". 
                                The number of batches (of tasks) to be pre-dispatched. [default = 2*n_jobs]
            confidence (float | str): Confidence threshold for limiting taxonomic depth. 
                                      Set to "disable" to disable confidence calculation, 
                                      or 0 to calculate confidence but not apply it to limit the taxonomic depth of the assignments.
                                      [default = 0.7]
            read_orientation (str): Direction of reads with respect to reference sequences. [default = auto]
                                        - "same" : will cause reads to be classified unchanged; 
                                        - "reverse-complement" : will cause reads to be reversed and complemented prior to classification. 
                                        - "auto" : will autodetect orientation based on the confidence estimates for the first 100 reads.

        Returns:
            taxonomy (FeatureData[Taxonomy]): taxonomic classified data
        """
        
        classifier_file = self.load_classifier()
        
        taxonomy = feature_classifier.classify_sklearn(
            classifier = classifier_file,
            reads = self.reads,
            reads_per_batch = reads_per_batch,
            n_jobs = n_jobs,
            pre_dispatch = pre_dispatch,
            confidence = confidence,
            read_orientation = read_orientation)
        
        return taxonomy
    
    
    # filter features from table by taxonomy
    def filter_table_by_taxon(
        self, 
        table: qiime2.Artifact,
        taxonomy: qiime2.Artifact,
        include: str = None,
        exclude: str = None,
        query_delimiter: str = ',',
        mode: str = 'contains') -> qiime2.Artifact:
        
        """ filter_table_by_taxon
        
        filtering features based on taxonomy from table
        
        Args:
            table (FeatureTable[Frequency]): Feature table to be filtered
            taxonomy (FeatureData[Taxonomy]): Taxonomic annotations for features in the provided feature table. 
                                              All features in the feature table must have a corresponding taxonomic annotation. 
                                              Taxonomic annotations for features that are not present in the feature table will be ignored.
            include (str): One or more search terms that indicate which taxa should be included in the resulting table. 
                           If providing more than one term, terms should be delimited by the query-delimiter character. 
                           By default, all taxa will be included. [default = None]
            exclude (str): One or more search terms that indicate which taxa should be excluded from the resulting table. 
                           If providing more than one term, terms should be delimited by the query-delimiter character. 
                           By default, no taxa will be excluded. [default = None]
            query_delimiter (str): The string used to delimit multiple search terms provided to include or exclude. 
                                   This parameter should only need to be modified if the default delimiter (a comma) is used in the provided taxonomic annotations.
                                   [default = ,]
            mode (str): Mode for determining if a search term matches a taxonomic annotation. [default = contains]
                            - "contains" requires that the annotation has the term as a substring;
                            - "exact" requires that the annotation is a perfect match to a search term.

        Returns:
            filtered_table (FeatureTable[Frequency]): The taxonomy-filtered feature table
        """
        
        filtered_table = filter_table(
            table = table,
            taxonomy = taxonomy,
            include = include,
            exclude = exclude,
            query_delimiter = query_delimiter,
            mode = mode)
        
        return filtered_table
    
    
    # filter features from sequences by taxonomy
    def filter_seqs_by_taxon(
        self,
        sequences: qiime2.Artifact,
        taxonomy: qiime2.Artifact,
        include: str = None,
        exclude: str = None,
        query_delimiter: str = ',',
        mode: str = 'contains') -> qiime2.Artifact:
        
        """ filter_seqs_by_taxon
        
        filtering features based on taxonomy from sequences
        
        Args:
            sequences (FeatureData[Sequence]): Feature sequences to be filtered
            taxonomy (FeatureData[Taxonomy]): Taxonomic annotations for features in the provided feature table. 
                                              All features in the feature table must have a corresponding taxonomic annotation. 
                                              Taxonomic annotations for features that are not present in the feature table will be ignored.
            include (str): One or more search terms that indicate which taxa should be included in the resulting table. 
                           If providing more than one term, terms should be delimited by the query-delimiter character. 
                           By default, all taxa will be included. [default = None]
            exclude (str): One or more search terms that indicate which taxa should be excluded from the resulting table. 
                           If providing more than one term, terms should be delimited by the query-delimiter character. 
                           By default, no taxa will be excluded. [default = None]
            query_delimiter (str): The string used to delimit multiple search terms provided to include or exclude. 
                                   This parameter should only need to be modified if the default delimiter (a comma) is used in the provided taxonomic annotations.
                                   [default = ,]
            mode (str): Mode for determining if a search term matches a taxonomic annotation. [default = contains]
                            - "contains" requires that the annotation has the term as a substring;
                            - "exact" requires that the annotation is a perfect match to a search term.

        Returns:
            filtered_sequences (FeatureData[Sequence]): The taxonomy-filtered feature sequences
        """
        
        filtered_sequences = filter_seqs(
            sequences = sequences,
            taxonomy = taxonomy,
            include = include,
            exclude = exclude,
            query_delimiter = query_delimiter,
            mode = mode)
        
        return filtered_sequences