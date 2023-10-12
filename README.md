# smel

A pipeline tool for automatic analysis of 16S amplicon sequences and metagenomes.
 * Amplicon analysis pipeline is based on the [QIIME2](https://qiime2.org/).
    * both paired end reads and pacbio ccs reads are available
 * Shotgun metagenome analysis pipeline is not ready.

 ### version
 v0.0.3

## simple usage
* for getting help
```shell
smel -h
or
smel --help
```

* to download classifiers
```shell
smel amplicon download -c classifier
```

* for amplicon sequences
```shell
# for paired end reads
smel amplicon run -d rawdata_directory -m metadata.txt

# for pacbio ccs reads
smel amplicon run -d rawdata_directory -m metadata.txt --adapter_front forward_adapter_sequence --adapter_rev reverse_adapter_sequence --ccs
```

* for shotgun metagenome (coming soon)
```shell
smel shotgun ???
```

## installation with conda
* first download **`yml`** file which is located in **`yaml`** directory from this repository.
    * [smel-py38-linux-conda.yml](https://drive.google.com/file/d/1Peretn-BN0Awo6TdBm6_NAfVq8n7c6P-/view?usp=sharing)
    * [smel-py38-osx-conda.yml](https://drive.google.com/file/d/1YcT5CJlmoGtO0C2QHiXeNsANfEPYuebc/view?usp=sharing)
```shell
# linux
conda env create -n smel --file smel-py38-linux-conda.yml

# macOS
conda env create -n smel --file smel-py38-linux-conda.yml
```

* if you have already installed QIIME2 with conda,
```shell
pip install smel
```

## Amplicon analysis
#### input file
* rawdata_directory
    * includes rawdata fastq files

* metadata.txt
    * QIIME2 metadata format

    | #SampleID | site | pH | $\cdots$ |
    | --------- | ---------- | ---------- | --------- |
    | #q2:types | categorical | numeric | $\cdots$ |
    | sample_1 | europe | 5.4 | $\cdots$ |
    | sample_2 | asia | 7.4 | $\cdots$ |

    * columns are tab-seperated

#### output files
* *.out
    * dada2.out
        * `dna-sequences.fasta` : denoised sequences
        * `feature-table.biom` : biom format of ASV table
        * `feature-table.txt` : txt format of ASV table
        * `stats.tsv` : DADA2 statistics
    * diversity.out
        * normalization : normalized without rarefaction (CLR/GeTMM/TMM/TSS)
            * `[clr/getmm/tmm/tss]_alpha_diversity.csv` : alpha diversity table
            * `[clr/getmm/tmm/tss]_normalized.csv` : normalized table
            * `[clr/getmm/tmm/tss]_pc_values.csv` : PC values
            * `[clr/getmm/tmm/tss]_pcoa_eigen_values.csv` : eigen values
            * `[clr/getmm/tmm/tss]_pcoa_proportion_explained.csv` : proportion explained ratio
        * rarefaction : normalized with rarefaction
            * `[bray_curtis/jaccard/unweighted_unifrac/weighted_unifrac]_distance_matrix.tsv` : distance matrix table
            * `[bray_curtis/jaccard/unweighted_unifrac/weighted_unifrac]_pcoa_results.txt` : PCoA results
            * `[evenness/faith_pd/observed_features/shannon]_vector.tsv` : alpha diversity results
            * `feature-table.biom` : rarefied ASV table
    * phylogeny.out
        * `aligned-dna-sequences.fasta` : multisequence aligned sequences
        * `tree.nwk` : newick format of rooted tree file
    * tax_filtered.out
        * `dna-sequences.fasta` : chloroplast and mitochondria filtered sequences
        * `feature-table.biom` : biom format of chloroplast and mitochondria filtered ASV table
        * `feature-table.txt` : txt format of chloroplast and mitochondria filtered ASV table
    * taxonomy.out
        * `taxonomy.tsv` : result of taxonomic assignments
* q2_files
    * `alignment.qza` : QIIME2 artifact of alignment
    * `dada2_seqs.qza` : QIIME2 artifact of DADA2 sequences
    * `dada2_stats.qza` : QIIME2 artifact of DADA2 stats
    * `dada2_table.qza` : QIIME2 artifact of DADA2 table
    * `filtered_dada2_seqs.qza` : QIIME2 artifact of chloroplast and mitochondria filtered sequences
    * `filtered_dada2_table.qza` : QIIME2 artifact of chloroplast and mitochondria filtered table
    * `imported.qza` : QIIME2 artifact of imported sequences
    * `masked_alignment.qza` : QIIME2 artifact of masked alignment
    * `rooted_tree.qza` : QIIME2 artifact of rooted tree
    * `taxonomy.qza` : QIIME2 artifact of taxonomic assignments
    * `tree.qza` : QIIME2 artifact of unrooted tree