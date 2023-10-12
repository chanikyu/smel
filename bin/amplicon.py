#!/usr/bin/env python

import os
import sys
sys.path.append(os.path.dirname(
    os.path.abspath(os.path.dirname(__file__))))

import numpy as np

import wget
import shutil

import pyfiglet
import qiime2

from termcolor import colored
from termcolor import cprint

from src import UTILS
from src import SEQ
from src import IMPORT
from src import ADAPTER
from src import DENOISE
from src import CLASSIFICATION
from src import PHYLOGENY
from src import DIVERSITY


# download classifiers
def download(classifier: str) -> None:
    
    silva_url = 'https://data.qiime2.org/2023.7/common/silva-138-99-nb-classifier.qza'
    greengenes2_url = 'https://data.qiime2.org/classifiers/greengenes/gg_2022_10_backbone_full_length.nb.qza'
    
    
    output = os.path.join(os.environ['CONDA_PREFIX'], 'classifiers')
    util = UTILS()
    
    util.make_directory(dirname=output)
    
    if classifier.lower() == 'silva':
        wget.download(silva_url, 
                      out=output)
        
    elif classifier.lower() == 'greengenes2':
        wget.download(greengenes2_url,
                      out=output)


# paired end reads analysis
def run(dirname: str, 
        metadata: str,
        adapter_front: str,
        adapter_tail: str,
        classifier: str,
        sampling_depth: int,
        threads: int,
        ccs: bool) -> None:
    
    print()
    
    if ccs:
        starting_text = f'Amplicon Analysis with PacBio CCS'
    else:
        starting_text = f'Amplicon Analysis with Paired End'
    
    cprint(pyfiglet.figlet_format('=' * 15, font='slant'), 'green')
    cprint(pyfiglet.figlet_format(starting_text, font='slant'), 'green')
    cprint(pyfiglet.figlet_format('=' * 15, font='slant'), 'green')
    
    print()
    
    cprint('Creating Manifest File', 'green')
    imp = IMPORT(dirname=dirname)
    if not os.path.exists(os.path.join(os.getcwd(), 'manifest')):
        imp.create_manifest(metadata=metadata, 
                            output=os.path.join(os.getcwd(), 'manifest'))
        print('  -> Saved manifest file\n')
    else:
        print('  -> Already created manifest file\n')
    
    
    cprint('Importing Sequences to QIIME2 Artifact', 'green')
    util = UTILS()
    if not os.path.exists(os.path.join(os.getcwd(), 'q2_files')):
        util.make_directory(dirname=os.path.join(os.getcwd(), 'q2_files'))
    
    if not os.path.exists(os.path.join(os.getcwd(), 'q2_files', 'imported.qza')):
        imported_data = imp.load_seq_artifact(manifest=os.path.join(os.getcwd(), 'manifest'))
        imported_data.save(os.path.join(os.getcwd(), 'q2_files', 'imported.qza'))
        print('  -> Imported sequence data')
        print(f'''  -> Saved as "{os.path.join(os.getcwd(), 'q2_files', 'imported.qza')}"\n''')
    else:
        imported_data = qiime2.Artifact.load(os.path.join(os.getcwd(), 'q2_files', 'imported.qza'))
        print('  -> Already imported sequence data\n')
        
    
    cprint('ASV Denoising with DADA2 :', 'green', end=' ')
    if ccs:
        if not os.path.exists(os.path.join(os.getcwd(), 'q2_files', 'dada2_table.qza')) \
            and not os.path.exists(os.path.join(os.getcwd(), 'q2_files', 'dada2_seqs.qza')):
                
            cprint('PacBio CCS', 'red')
            
            seq = SEQ()
            denoi = DENOISE(demuxed_data=imported_data)
            
            table, rep_seqs, stats = denoi.runDADA2_CCS(
                front=adapter_front,
                adapter=adapter_tail,
                max_mismatch=2,
                indels=False,
                trunc_len=0,
                trim_left=0,
                max_ee=2.0,
                trunc_q=2,
                min_len=1000,
                max_len=1600,
                pooling_method='independent',
                chimera_method='consensus',
                min_fold_parent_over_abundance=3.5,
                allow_one_off=False,
                n_threads=threads,
                n_reads_learn=1000000,
                hashed_feature_ids=True)
            
            table.save(os.path.join(os.getcwd(), 'q2_files', 'dada2_table.qza'))
            util.export_q2_artifact(
                artifact=table, 
                output=os.path.join(os.getcwd(), 'dada2.out'))
            cvt_biom = util.convert_biom(
                biom_table=os.path.join(os.getcwd(), 'dada2.out', 'feature-table.biom'))
            cvt_biom.to_csv(
                os.path.join(os.getcwd(), 'dada2.out', 'feature-table.txt'), 
                sep='\t',
                header=True,
                index=True)
            print(f'''  -> Saved as "{os.path.join(os.getcwd(), 'q2_files', 'dada2_table.qza')}"''')
            
            rep_seqs.save(os.path.join(os.getcwd(), 'q2_files', 'dada2_seqs.qza'))
            util.export_q2_artifact(
                artifact=rep_seqs, 
                output=os.path.join(os.getcwd(), 'dada2.out'))
            print(f'''  -> Saved as "{os.path.join(os.getcwd(), 'q2_files', 'dada2_seqs.qza')}"''')
            
            stats.save(os.path.join(os.getcwd(), 'q2_files', 'dada2_stats.qza'))
            util.export_q2_artifact(
                artifact=stats,
                output=os.path.join(os.getcwd(), 'dada2.out'))
            print(f'''  -> Saved as "{os.path.join(os.getcwd(), 'q2_files', 'dada2_stats.qza')}"''')
            
            print('  -> Finished DADA2\n')
    
    else:
        if not os.path.exists(os.path.join(os.getcwd(), 'q2_files', 'dada2_table.qza')) \
            and not os.path.exists(os.path.join(os.getcwd(), 'q2_files', 'dada2_seqs.qza')):
                
            cprint('Paired-End Reads', 'red')
            
            print('  -> Calculating Trimming Positions')
            
            seq = SEQ()
            denoi = DENOISE(demuxed_data=imported_data)
            
            fwd_ids = []
            fwd_seqs = []
            fwd_quals = []
            
            manifest = util.load_dataframe(table=os.path.join(os.getcwd(), 'manifest'))
            for idx in range(0, len(manifest), 2):
                fwd = manifest.iloc[idx]['absolute-filepath']
                rev = manifest.iloc[idx + 1]['absolute-filepath']
                
                print(f'  -> Forward read: {os.path.basename(fwd)}, Reverse read: {os.path.basename(rev)}')
                
                fwd_read_id, fwd_read, fwd_qual = seq.parse_fastq(fastq=fwd)
                
                fwd_ids.append(fwd_read_id)
                fwd_seqs.append(fwd_read)
                for q in fwd_qual:
                    fwd_quals.append(q)
                
            print(f'  -> A Total {len(fwd_quals):,} Reads to Calculate')
            print(f'  -> A Total {np.sum([len(q) for q in fwd_quals]):,} Bases to Calculate')
            
            fwd_mean_results = denoi.mean_quality_scores_at_each_position(qscores = fwd_quals)
            
            fwd_loess, fwd_cooks_dist = denoi.loess_to_cooks_distance(
                means = fwd_mean_results, 
                frac = 0.2)
            # print(type(fwd_loess))
            # print(fwd_loess)
            
            fwd_trim_pos, fwd_trunc_pos = denoi.determine_start_end_position(
                cooks_distance = fwd_cooks_dist,
                means = fwd_mean_results,
                threshold = 0.05)
            
            print()
            print(f'  -> The trimming position of read is {fwd_trim_pos}')
            print(f'  -> The truncation position of read is {fwd_trunc_pos // 2}\n')
        
            table, rep_seqs, stats = denoi.runDADA2_paired(
                trunc_len_f=fwd_trunc_pos // 2,
                trunc_len_r=fwd_trunc_pos // 2,
                trim_left_f=fwd_trim_pos,
                trim_left_r=fwd_trim_pos,
                max_ee_f=2.0,
                max_ee_r=2.0,
                trunc_q=15,
                min_overlap=12,
                pooling_method='independent',
                chimera_method='consensus',
                min_fold_parent_over_abundance=1.0,
                allow_one_off=False,
                n_threads=threads,
                n_reads_learn=1000000,
                hashed_feature_ids=True)
            
            table.save(os.path.join(os.getcwd(), 'q2_files', 'dada2_table.qza'))
            util.export_q2_artifact(
                artifact=table, 
                output=os.path.join(os.getcwd(), 'dada2.out'))
            cvt_biom = util.convert_biom(
                biom_table=os.path.join(os.getcwd(), 'dada2.out', 'feature-table.biom'))
            cvt_biom.to_csv(
                os.path.join(os.getcwd(), 'dada2.out', 'feature-table.txt'), 
                sep='\t',
                header=True,
                index=True)
            print(f'''  -> Saved as "{os.path.join(os.getcwd(), 'q2_files', 'dada2_table.qza')}"''')
            
            rep_seqs.save(os.path.join(os.getcwd(), 'q2_files', 'dada2_seqs.qza'))
            util.export_q2_artifact(
                artifact=rep_seqs, 
                output=os.path.join(os.getcwd(), 'dada2.out'))
            print(f'''  -> Saved as "{os.path.join(os.getcwd(), 'q2_files', 'dada2_seqs.qza')}"''')
            
            stats.save(os.path.join(os.getcwd(), 'q2_files', 'dada2_stats.qza'))
            util.export_q2_artifact(
                artifact=stats,
                output=os.path.join(os.getcwd(), 'dada2.out'))
            print(f'''  -> Saved as "{os.path.join(os.getcwd(), 'q2_files', 'dada2_stats.qza')}"''')
            
            print('  -> Finished DADA2\n')
            
        else:
            table = qiime2.Artifact.load(os.path.join(os.getcwd(), 'q2_files', 'dada2_table.qza'))
            rep_seqs = qiime2.Artifact.load(os.path.join(os.getcwd(), 'q2_files', 'dada2_seqs.qza'))
            print('\n  -> Already denoised sequences\n')
    
    
    cprint('Taxonomic Classification', 'green')
    classification = CLASSIFICATION(
        reads=rep_seqs, 
        classifier=classifier)
    
    if not os.path.exists(os.path.join(os.getcwd(), 'q2_files', 'taxonomy.qza')):
        taxonomy, = classification.taxonomy_classification(
            reads_per_batch='auto',
            n_jobs=threads,
            pre_dispatch='2*n_jobs',
            confidence=0.7,
            read_orientation='auto')
        taxonomy.save(os.path.join(os.getcwd(), 'q2_files', 'taxonomy.qza'))
        util.export_q2_artifact(
            artifact=taxonomy, 
            output=os.path.join(os.getcwd(), 'taxonomy.out'))
        print('  -> Reads are classified')
        print(f'''  -> Saved as "{os.path.join(os.getcwd(), 'q2_files', 'taxonomy.qza')}"''')
        
    else:
        taxonomy = qiime2.Artifact.load(os.path.join(os.getcwd(), 'q2_files', 'taxonomy.qza'))
        print('  -> Already classified sequences\n')
    
    
    cprint('Filtration of Mitochondria and Chloroplast', 'green')
    if not os.path.exists(os.path.join(os.getcwd(), 'q2_files', 'filtered_dada2_table.qza')) and not os.path.exists(os.path.join(os.getcwd(), 'q2_files', 'filtered_dada2_seqs.qza')):
        filtered_table = classification.filter_table_by_taxon(
            table=table,
            taxonomy=taxonomy,
            include=None,
            exclude='Mitochondria,Chloroplast,Plastid',
            query_delimiter=',',
            mode='contains').filtered_table
        filtered_sequences = classification.filter_seqs_by_taxon(
            sequences=rep_seqs,
            taxonomy=taxonomy,
            include=None,
            exclude='Mitochondria,Chloroplast,Plastid',
            query_delimiter=',',
            mode='contains').filtered_sequences
        
        print('  -> ASV Frequency Table and Sequences are Filtered by Mitochondria and Chloroplast (Plastid)')
        filtered_table.save(os.path.join(os.getcwd(), 'q2_files', 'filtered_dada2_table.qza'))
        util.export_q2_artifact(
            artifact=filtered_table,
            output=os.path.join(os.getcwd(), 'tax_filtered.out'))
        cvt_biom2 = util.convert_biom(
            biom_table=os.path.join(os.getcwd(), 'tax_filtered.out', 'feature-table.biom'))
        cvt_biom2.to_csv(
            os.path.join(os.getcwd(), 'tax_filtered.out', 'feature-table.txt'),
            sep='\t', 
            header=True, 
            index=True)
        print(f'''  -> Saved as "{os.path.join(os.getcwd(), 'q2_files', 'filtered_dada2_table.qza')}"''')
        
        filtered_sequences.save(os.path.join(os.getcwd(), 'q2_files', 'filtered_dada2_seqs.qza'))
        util.export_q2_artifact(
            artifact=filtered_sequences,
            output=os.path.join(os.getcwd(), 'tax_filtered.out'))
        print(f'''  -> Saved as "{os.path.join(os.getcwd(), 'q2_files', 'filtered_dada2_seqs.qza')}"''')
        
        print()
        
    else:
        filtered_table = qiime2.Artifact.load(os.path.join(os.getcwd(), 'q2_files', 'filtered_dada2_table.qza'))
        filtered_sequences = qiime2.Artifact.load(os.path.join(os.getcwd(), 'q2_files', 'filtered_dada2_seqs.qza'))
        print('  -> Already filtered taxa from table and sequences\n')
        
        
    cprint('Generate a Phylogenetic Tree for Diversity Analyses', 'green')
    if not os.path.exists(os.path.join(os.getcwd(), 'q2_files', 'alignment.qza')) and not os.path.exists(os.path.join(os.getcwd(), 'q2_files', 'rooted_tree.qza')):
        phylo = PHYLOGENY(
            sequence=filtered_sequences, 
            threads=threads)
        
        alignment, masked_alignment, tree, rooted_tree = phylo.aln2fasttree(
            mask_max_gap_frequency=1.0,
            mask_min_conservation=0.4,
            parttree=True)
        
        print('  -> Reads are aligned and phylogenetic tree generated')
        
        alignment.save(os.path.join(os.getcwd(), 'q2_files', 'alignment.qza'))
        util.export_q2_artifact(
            artifact=alignment,
            output=os.path.join(os.getcwd(), 'phylogeny.out'))
        print(f'''  -> Saved as "{os.path.join(os.getcwd(), 'q2_files', 'alignment.qza')}"''')
        
        masked_alignment.save(os.path.join(os.getcwd(), 'q2_files', 'masked_alignment.qza'))
        util.export_q2_artifact(
            artifact=masked_alignment,
            output=os.path.join(os.getcwd(), 'phylogeny.out'))
        print(f'''  -> Saved as "{os.path.join(os.getcwd(), 'q2_files', 'masked_alignment.qza')}"''')
        
        tree.save(os.path.join(os.getcwd(), 'q2_files', 'tree.qza'))
        util.export_q2_artifact(
            artifact=tree,
            output=os.path.join(os.getcwd(), 'phylogeny.out'))
        print(f'''  -> Saved as "{os.path.join(os.getcwd(), 'q2_files', 'tree.qza')}"''')
        
        rooted_tree.save(os.path.join(os.getcwd(), 'q2_files', 'rooted_tree.qza'))
        util.export_q2_artifact(
            artifact=rooted_tree,
            output=os.path.join(os.getcwd(), 'phylogeny.out'))
        print(f'''  -> Saved as "{os.path.join(os.getcwd(), 'q2_files', 'rooted_tree.qza')}"\n''')
        
    else:
        alignment = qiime2.Artifact.load(os.path.join(os.getcwd(), 'q2_files', 'alignment.qza'))
        masked_alignment = qiime2.Artifact.load(os.path.join(os.getcwd(), 'q2_files', 'masked_alignment.qza'))
        tree = qiime2.Artifact.load(os.path.join(os.getcwd(), 'q2_files', 'tree.qza'))
        rooted_tree = qiime2.Artifact.load(os.path.join(os.getcwd(), 'q2_files', 'rooted_tree.qza'))
        
        print('  -> Already conducted alignment and generation of phylogenetic tree\n')
        
    
    cprint('Diversity Analysis', 'green')
    if not os.path.exists(os.path.join(os.getcwd(), 'diversity.out')):
        
        util.make_directory(dirname=os.path.join(os.getcwd(), 'diversity.out'))
        biom = os.path.join(os.getcwd(), 'tax_filtered.out', 'feature-table.biom')
        div = DIVERSITY(biom=biom)
        
        print('  1. Diversity Analysis with Rarefied Sequences')
        if not os.path.exists(os.path.join(os.getcwd(), 'diversity.out', 'rarefaction')):
            
            util.make_directory(dirname=os.path.join(os.getcwd(), 'diversity.out', 'rarefaction'))
            
            depth, failed_ratio = div.calculate_rarefaction_depth(
                depth=sampling_depth,
                fail_ratio=0.3)
            
            if not depth == 0 and not failed_ratio == 0:
                print(f'    -> Calculated sampling depth as {depth:,}')
                print(f'    -> {failed_ratio * 100:.2f}% of samples are not passed\n')
                
                _metadata_ = qiime2.Metadata.load(metadata)
                
                results = div.diversity_with_rarefaction(
                    feature_table=filtered_table,
                    rooted_tree=rooted_tree,
                    metadata=_metadata_,
                    depth=depth,
                    replacement=False,
                    threads=threads)
                    
                result_names = ['rarefied_table', 'faith_pd_vector', 'observed_features_vector', 'shannon_vector', 
                                'evenness_vector', 'unweighted_unifrac_distance_matrix', 'weighted_unifrac_distance_matrix', 
                                'jaccard_distance_matrix', 'bray_curtis_distance_matrix', 'unweighted_unifrac_pcoa_results', 
                                'weighted_unifrac_pcoa_results', 'jaccard_pcoa_results', 'bray_curtis_pcoa_results', 
                                'unweighted_unifrac_emperor', 'weighted_unifrac_emperor', 'jaccard_emperor', 'bray_curtis_emperor']
                
                for result_name, result_value in zip(result_names, results):
                    result_value.save(os.path.join(os.getcwd(), 'q2_files', f'{result_name}.qza'))
                    
                    if not '_emperor' in result_name:
                        util.export_q2_artifact(
                            artifact=result_value, 
                            output=os.path.join(os.getcwd(), 'diversity.out', 'rarefaction'))
                        
                        if 'vector' in result_name:
                            shutil.move(src=os.path.join(os.getcwd(), 'diversity.out', 'rarefaction', 'alpha-diversity.tsv'),
                                        dst=os.path.join(os.getcwd(), 'diversity.out', 'rarefaction', f'{result_name}.tsv'))
                            
                        elif 'matrix' in result_name:
                            shutil.move(src=os.path.join(os.getcwd(), 'diversity.out', 'rarefaction', 'distance-matrix.tsv'),
                                        dst=os.path.join(os.getcwd(), 'diversity.out', 'rarefaction', f'{result_name}.tsv'))
                            
                        elif 'pcoa' in result_name:
                            shutil.move(src=os.path.join(os.getcwd(), 'diversity.out', 'rarefaction', 'ordination.txt'),
                                        dst=os.path.join(os.getcwd(), 'diversity.out', 'rarefaction', f'{result_name}.txt'))
                    
                    print(f'''    -> Saved as "{os.path.join(os.getcwd(), 'q2_files', f'{result_name}.qza')}"''')
                    print()
                
        else:
            print('    -> Alpha and beta diversity analyses with rarefaction were already performed\n')
        
        
        print('  2. Diversity Analysis without Rarefied Sequences: Normalization and Standardization')
        if not os.path.exists(os.path.join(os.getcwd(), 'diversity.out', 'normalization')):
            
            util.make_directory(dirname=os.path.join(os.getcwd(), 'diversity.out', 'normalization'))
            
            # tss
            try:
                tss_normed = div.tss_norm()
                tss_normed.to_csv(os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'tss_normalized.csv'), 
                                index=True, 
                                header=True)
                print(f'''    -> Saved as "{os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'tss_normalized.csv')}"''')
                
                tss_alpha_diversity_result = div.calculate_alpha_diversity(
                    normed_table=tss_normed, 
                    rooted_tree=os.path.join(os.getcwd(), 'phylogeny.out', 'tree.nwk'))
                tss_alpha_diversity_result.to_csv(os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'tss_alpha_diversity.csv'),
                                                index=True,
                                                header=True)
                print(f'''    -> Saved as "{os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'tss_alpha_diversity.csv')}"''')
                
                
                tss_pcoa_vals, tss_pcoa_proportion, tss_pcoa_eigvals = div.calculate_beta_diversity(
                    normed_table=tss_normed,
                    rooted_tree=os.path.join(os.getcwd(), 'phylogeny.out', 'tree.nwk')
                )
                tss_pcoa_vals.to_csv(os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'tss_pc_values.csv'), 
                                    index=True, header=True)
                print(f'''    -> Saved as "{os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'tss_pc_values.csv')}"''')
                tss_pcoa_proportion.to_csv(os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'tss_pcoa_proportion_explained.csv'), 
                                    index=True, header=True)
                print(f'''    -> Saved as "{os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'tss_pcoa_proportion_explained.csv')}"''')
                tss_pcoa_eigvals.to_csv(os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'tss_pcoa_eigen_values.csv'), 
                                    index=True, header=True)
                print(f'''    -> Saved as "{os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'tss_pcoa_eigen_values.csv')}"\n''')
            except:
                cprint('    -> ValueError emerged. Skip TSS normalization step.', 'red')
            
            
            # clr
            try:
                clr_normed = div.clr_norm()
                clr_normed.to_csv(os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'clr_normalized.csv'), 
                                index=True, 
                                header=True)
                print(f'''    -> Saved as "{os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'clr_normalized.csv')}"''')
                
                clr_alpha_diversity_result = div.calculate_alpha_diversity(
                    normed_table=clr_normed, 
                    rooted_tree=os.path.join(os.getcwd(), 'phylogeny.out', 'tree.nwk'))
                clr_alpha_diversity_result.to_csv(os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'clr_alpha_diversity.csv'),
                                                index=True,
                                                header=True)
                print(f'''    -> Saved as "{os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'clr_alpha_diversity.csv')}"''')
                
                clr_pcoa_vals, clr_pcoa_proportion, clr_pcoa_eigvals = div.calculate_beta_diversity(
                    normed_table=clr_normed,
                    rooted_tree=os.path.join(os.getcwd(), 'phylogeny.out', 'tree.nwk')
                )
                clr_pcoa_vals.to_csv(os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'clr_pc_values.csv'), 
                                    index=True, header=True)
                print(f'''    -> Saved as "{os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'clr_pc_values.csv')}"''')
                clr_pcoa_proportion.to_csv(os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'clr_pcoa_proportion_explained.csv'), 
                                    index=True, header=True)
                print(f'''  -> Saved as "{os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'clr_pcoa_proportion_explained.csv')}"''')
                clr_pcoa_eigvals.to_csv(os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'clr_pcoa_eigen_values.csv'), 
                                    index=True, header=True)
                print(f'''    -> Saved as "{os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'clr_pcoa_eigen_values.csv')}"\n''')
            except:
                cprint('  -> ValueError emerged. Skip this CLR normalization step.', 'red')
            
            
            # tmm
            try:
                tmm_normed, norm_factors = div.tmm_norm(
                    trim_lfc=0.3, 
                    trim_mag=0.05)
                tmm_normed.to_csv(os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'tmm_normalized.csv'), 
                                index=True, 
                                header=True)
                print(f'''    -> Saved as "{os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'tmm_normalized.csv')}"''')
                
                tmm_alpha_diversity_result = div.calculate_alpha_diversity(
                    normed_table=tmm_normed, 
                    rooted_tree=os.path.join(os.getcwd(), 'phylogeny.out', 'tree.nwk'))
                tmm_alpha_diversity_result.to_csv(os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'tmm_alpha_diversity.csv'),
                                                index=True,
                                                header=True)
                print(f'''    -> Saved as "{os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'tmm_alpha_diversity.csv')}"''')
                
                tmm_pcoa_vals, tmm_pcoa_proportion, tmm_pcoa_eigvals = div.calculate_beta_diversity(
                    normed_table=tmm_normed,
                    rooted_tree=os.path.join(os.getcwd(), 'phylogeny.out', 'tree.nwk')
                )
                tmm_pcoa_vals.to_csv(os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'tmm_pc_values.csv'), 
                                    index=True, header=True)
                print(f'''    -> Saved as "{os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'tmm_pc_values.csv')}"''')
                tmm_pcoa_proportion.to_csv(os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'tmm_pcoa_proportion_explained.csv'), 
                                    index=True, header=True)
                print(f'''    -> Saved as "{os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'tmm_pcoa_proportion_explained.csv')}"''')
                tmm_pcoa_eigvals.to_csv(os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'tmm_pcoa_eigen_values.csv'), 
                                    index=True, header=True)
                print(f'''    -> Saved as "{os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'tmm_pcoa_eigen_values.csv')}"\n''')
            except:
                cprint('    -> ValueError emerged. Skip this TMM normalization step.', 'red')
            
            
            # getmm
            try:
                getmm_normed = div.getmm_norm(fasta=os.path.join(os.getcwd(), 'tax_filtered.out', 'dna-sequences.fasta'))
                getmm_normed.to_csv(os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'getmm_normalized.csv'), 
                                index=True, 
                                header=True)
                print(f'''  -> Saved as "{os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'getmm_normalized.csv')}"''')
                
                getmm_alpha_diversity_result = div.calculate_alpha_diversity(
                    normed_table=getmm_normed, 
                    rooted_tree=os.path.join(os.getcwd(), 'phylogeny.out', 'tree.nwk'))
                getmm_alpha_diversity_result.to_csv(os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'getmm_alpha_diversity.csv'),
                                                    index=True,
                                                    header=True)
                print(f'''    -> Saved as "{os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'getmm_alpha_diversity.csv')}"''')
                
                getmm_pcoa_vals, getmm_pcoa_proportion, getmm_pcoa_eigvals = div.calculate_beta_diversity(
                    normed_table=getmm_normed,
                    rooted_tree=os.path.join(os.getcwd(), 'phylogeny.out', 'tree.nwk')
                )
                getmm_pcoa_vals.to_csv(os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'getmm_pc_values.csv'), 
                                    index=True, header=True)
                print(f'''    -> Saved as "{os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'getmm_pc_values.csv')}"''')
                getmm_pcoa_proportion.to_csv(os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'getmm_pcoa_proportion_explained.csv'), 
                                    index=True, header=True)
                print(f'''    -> Saved as "{os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'getmm_pcoa_proportion_explained.csv')}"''')
                getmm_pcoa_eigvals.to_csv(os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'getmm_pcoa_eigen_values.csv'), 
                                    index=True, header=True)
                print(f'''    -> Saved as "{os.path.join(os.getcwd(), 'diversity.out', 'normalization', 'getmm_pcoa_eigen_values.csv')}"\n''')
            except:
                cprint('    -> ValueError emerged. Skip this GeTMM normalization step.', 'red')
            
        
        else:
            print('    -> Alpha and beta diversity analyses without rarefaction were already performed\n')
            
    else:
        print('    -> Alpha and beta diversity analyses were already performed\n')
        
        
    cprint('All process has been completed.', 'green')