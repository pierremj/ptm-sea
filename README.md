# A custom ptm-sea wrapper

A wrapper for running PTM-SEA from a phosphoproteomics differential abundance analysis summary table

## Introduction

The `ptmsea.R` script performs the following computations:

1.  Parse differential abundance analysis (DA) results from a CSV file

2.  Calculate signed logp values from fold-changes and t-test p-values

3.  Extract a motif surrounding the phosphosite by extracting +/-7 flanking amino acids

4.  Run ptm-sea by calling the `ssGSEA2.0` function

## Usage

A command line interphase is provided through `run_ptmsea.R`. This requires the following inputs, with defaults provided for internally used FASTA file and GMT gene set file.

    Options:
            -d CHARACTER, --dea=CHARACTER
                    CSV file with differential abundance results

            -o CHARACTER, --out=CHARACTER
                    Prefix to be added to output folder and files

            --gene_set=CHARACTER
                    PTM-SigDB gene set database 

            --fasta=CHARACTER
                    FASTA file used in proteome search 

            --out_dir=CHARACTER
                    Output directory[default= ./output]

            --ssgsea_dir=CHARACTER
                    Folder containing ssGSEA2.0 source files [default= ./src/ssGSEA2.0/]

            -h, --help
                    Show this help message and exit

Example:

    Rscript run_ptmsea.R -d ./sample_data/dea_PXD030674.csv -o test --gene_set ./sample_data/ptm.sig.db.all.flanking.mouse.v1.9.0.gmt --fasta ./sample_data/mouse.UP000000589_10090_with_isoforms.fasta --out_dir ./output --ssgsea_dir ./src/ssGSEA2.0/


## Requirements of *dea* Input File
The input file should be a CSV table with the following columns:
* `feature_label` The Uniprot Accession ID and the modified site in Proteome Discoverer format (PD column Modifications in Master Proteins). For example `Q80UK8 1xPhospho [S970]` for singly modified and `Q8VDM4 2xPhospho [S789; S802]` for double modified
* `contrast` Name of the comparison performed during differential analysis
* `logfc` The fold-change in log scale
* `pval` The nominal p-value of the t-test 


## Output

A folder named `<output_prefix>_gsea` will be created in the output directory. The folder will contain a GCT file with the signed log(p-values) extracted from the differential abundance analysis results. This file is used to run ssGSEA2.0. For more information about ssGSEA2.0 output, please consult the [ssGSEA2.0 package documentation](https://github.com/broadinstitute/ssGSEA2.0).

## References

The ssGSEA2.0 library is used as the main functionality of this project. If useful, please cite the original publication:

Krug, K., Mertins, P., Zhang, B., Hornbeck, P., Raju, R., Ahmad, R., . Szucs, M., Mundt, F., Forestier, D., Jane-Valbuena, J., Keshishian, H., Gillette, M. A., Tamayo, P., Mesirov, J. P., Jaffe, J. D., Carr, S. A., Mani, D. R. (2019). **A curated resource for phosphosite-specific signature analysis**, Molecular & Cellular Proteomics (in Press). <http://doi.org/10.1074/mcp.TIR118.000943>
