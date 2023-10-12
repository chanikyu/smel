#!/usr/bin/env python


import qiime2
from qiime2.plugins.phylogeny.pipelines import align_to_tree_mafft_fasttree


class PHYLOGENY:
    def __init__(self, 
                 sequence: qiime2.Artifact,
                 threads: int) -> None:
        self.sequence = sequence
        self.threads = threads
    
    
    # alignment to fasttree
    def aln2fasttree(
        self,
        mask_max_gap_frequency: float = 1.0,
        mask_min_conservation: float = 0.4,
        parttree: bool = False) -> tuple:
        
        """ aln2fasttree
        
        alignment of reads and generation of phylogenetic tree
        
        Args:
            mask_max_gap_frequency (float): The maximum relative frequency of gap characters in a column for the column to be retained. 
                                            This relative frequency must be a number between 0.0 and 1.0 (inclusive), 
                                            where 0.0 retains only those columns without gap characters, 
                                            and 1.0 retains all columns regardless of gap character frequency. 
                                            This value is used when masking the aligned sequences. [default = 1.0]
            mask_min_conservation (float): The minimum relative frequency of at least one non-gap character in a column for that column to be retained. 
                                           This relative frequency must be a number between 0.0 and 1.0 (inclusive). 
                                           For example, if a value of 0.4 is provided, a column will only be retained 
                                           if it contains at least one character that is present in at least 40% of the sequences.
                                           This value is used when masking the aligned sequences. [default = 0.4]
            parttree (bool): This flag is required if the number of sequences being aligned are larger than 1000000. [default = False]

        Returns:
            alignment (FeatureData[AlignedSequence]): The aligned sequences
            masked_alignment (FeatureData[AlignedSequence]): The masked alignment
            tree (Phylogeny[Unrooted]): The unrooted phylogenetic tree
            rooted_tree (Phylogeny[Rooted]): The rooted phylogenetic tree
        """
        
        alignment, masked_alignment, tree, rooted_tree = align_to_tree_mafft_fasttree(
            sequences = self.sequence,
            n_threads = self.threads,
            mask_max_gap_frequency = mask_max_gap_frequency,
            mask_min_conservation = mask_min_conservation,
            parttree = parttree)
        
        return (alignment, masked_alignment, tree, rooted_tree)