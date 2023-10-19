#!/usr/bin/env Rscript
library("optparse")
option_list = list(
  make_option(c("--dea","-d"), type="character", default=NULL, 
              help="CSV file with differential abundance results", metavar="character"),
  make_option(c("--out","-o"), type="character", default=NULL,  
              help="Prefix to be added to output folder and files", metavar="character"),
  make_option(c("--gene_set"), type="character", default="./data/ptm.sig.db.all.flanking.human.v1.9.0.gmt", 
              help="PTM-SigDB gene set database [default= %default]", metavar="character"),
  make_option(c("--fasta"), type="character", default="./data/Homo sapiens (SwissProt TaxID=9606_and_subtaxonomies) (1).fasta", 
              help="FASTA file used in proteome search [default= %default]", metavar="character"),
  make_option(c("--out_dir"), type="character", default="./output", 
              help="Output directory[default= %default]", metavar="character"),
  make_option(c("--ssgsea_dir"), type="character", default="./src/ssGSEA2.0/", 
              help="Folder containing ssGSEA2.0 source files [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if(is.null(opt$dea) | is.null(opt$out)){
  print_help(opt_parser)
  stop("Both input dea file and output name need to be provided\n
       Example:\n\tRscript run_ptmsea.R -d ./data/VAV1_Phos_PythonSummary_043021.xlsx -o VAV1")
}

script.dir <- opt$ssgsea_dir
source("./src/ptmsea.R") 



out_message <- run_ptmsea(input_dea = opt$dea, 
           fasta_input = opt$fasta, 
           out_path = opt$out_dir, 
           output_prefix = opt$out, 
           ssgsea2_dir = opt$ssgsea_dir, 
           gene_set_database = opt$gene_set)

print(out_message)



