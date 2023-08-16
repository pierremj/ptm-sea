#Libraries
library(tibble)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(cmapR) #Bioconductor
library(limma) #Bioconductor
library(Biostrings) #Bioconductor
library(devtools)
library(readxl)
library(readr)
library(pacman)

#' Run PTM-SEA
#' 
#' Perform PTM-SEA analysis using the ssGSEA2.0 algorithm developed
#' by Krug et al.
#'
#' @param input_dea 
#' @param fasta_input 
#' @param out_path 
#' @param output_prefix 
#' @param ssgsea2_dir 
#' @param gene_set_database 
#'
#' @return
#' @export
#'
#' @examples
#' run_ptmsea(input_dea = "./data/VAV1_Phos_PythonSummary_043021.xlsx",
#'  fasta_input <- "./data/Homo sapiens (SwissProt TaxID=9606_and_subtaxonomies) (1).fasta"
#'  out_path <- "./output/VAV1_gsea/"
#'  output_prefix <- "VAV1"
#'  ssgsea2_dir <- "./src/ssGSEA2.0/"
#'  gene_set_database <- "./data/ptm.sig.db.all.flanking.human.v1.9.0.gmt"
#'  )
run_ptmsea <- function(
    input_dea, 
    fasta_input, 
    out_path, 
    output_prefix, 
    ssgsea2_dir, 
    gene_set_database,
    pavito){
  
  #1. Extract p-val, logfc, and signed logpvals ---------------------------------
  if(!pavito){
    cat("Reading input excel file...\n")
    dea <- 
      read_xlsx(input_dea) %>%
      
      #Rename columns for convenience
      mutate(pep_id = `...1`,
             prot_id = FullAccession,
             seq_start = `Positions in Master Proteins`) %>%
      
      #Extract phosphorylation position
      mutate(ptm_site = Modifications) %>%
      separate_rows(ptm_site,sep = "; \\d") %>%
      filter(str_detect(ptm_site,"Phospho")) %>%
      mutate(ptm_site = str_match(ptm_site, "\\[(.+)\\]")[,2]) %>%
      separate_rows(ptm_site,sep = "; ") %>%
      mutate(ptm_site = str_remove(ptm_site,"\\(.+\\)")) %>%
      
      #Remove non-localized phosphosites
      filter(str_detect(ptm_site,"\\d")) %>%
      
      #Extract sequence start site
      mutate(seq_start = str_split_fixed(seq_start,";",2)[,1]) %>%
      mutate(seq_start = as.numeric(str_match(seq_start,"\\[(\\d+)-")[,2])) %>%
      
      #Offset phosphorylation position
      mutate(ptm_site = as.numeric(str_extract(ptm_site,"\\d+")) + seq_start - 1) %>% 
      
      #Select columns of interest
      select(prot_id,ptm_site,pep_id,starts_with("Ttest"),contains("/")) %>%
      select(-starts_with("Phos")) %>%
      select(-contains("_Phos")) %>%
      
      #Make into longer format for easy manipulation
      rename_with(function(x){
        grps <- str_match(x,".+(Grp\\d+).+(Grp\\d+)")
        str_c("fc;",grps[,2],"_v_",grps[,3])
      },starts_with("NormPhos")) %>%
      rename_with(function(x){
        grps <- str_match(x,".+(Grp\\d+).+(Grp\\d+)")
        str_c("pval;",grps[,3],"_v_",grps[,2])
      },starts_with("Ttest")) %>%
      pivot_longer(c(-prot_id,-ptm_site,-pep_id),
                   names_to = c(".value","contrast"),
                   names_sep = ";") %>%
      
      #Calculate signed log p-value 
      mutate(sign_logp = sign(log(fc))*-log10(pval))
    
    cat("Excel file parsing complete!\n")
  } else {
    cat("Reading input csv file...\n")
    dea <- 
      read_csv(input_dea) %>%
      
      mutate(
        prot_id = str_extract(feature_name, "^[[:alnum:]]+"),
        seq_start = str_match(feature_name,"^[[:alnum:]]+__(\\d+)")[,2],
        phosphosites = str_extract(feature_name, "\\dxphospho(__[y|t|s]\\d*)+")
      ) %>%
      mutate(phosphosites = str_extract_all(phosphosites,"_[s|t|y]\\d*")) %>%
      unnest(phosphosites) %>%
      filter(str_detect(phosphosites,"\\d")) %>%
      drop_na(logfc) %>%
      mutate(phosphosites = as.numeric(str_extract(phosphosites,"\\d+"))+ as.numeric(seq_start) - 1) %>%
      transmute(
        prot_id = toupper(prot_id),
        ptm_site = phosphosites,
        pep_id = feature_name,
        contrast = contrast,
        pval,fc = 2^logfc,
        sign_logp = -log10(pval)*sign(logfc)
      )
      
    
    cat("CSV file parsing complete!\n")
  }
  
  # 2. Add motif sequence ----------------------------------------------------
  
  #' Extract flanking sequences from a protein and phosphosite pair
  #'
  #' @param prot_id - Protein identifier as in names(fasta)
  #' @param ptm_site - Modification site as amino acid followed by number
  #' @param fasta - An AAStringSet created by Biostring
  #'
  #' @return +/-7 motif sequence
  extract_sequence <- function(prot_id,ptm_site,fasta){
    
    curr_fasta <- fasta[[prot_id]]
    seq_length <- length(curr_fasta)
    position <- ptm_site
    
    if(position > seq_length){
      return(NA)
    }
    
    start <- if(position < 8) 1 else position - 7 
    end <- if((seq_length - position) < 7) seq_length else position+7
    x <- as.character(curr_fasta[start:end])
    
    if(position + 7 > seq_length){
      x <- str_c(c(x,rep("_",7-(seq_length-position))),collapse="")
    }
    if(position <= 7){
      x <- str_c(c(rep("_",8-position),x),collapse="")
    }
    
    if(!(str_sub(x,8,8) %in% c("S","T","Y"))){
      return(NA)
    }
    
    x <- str_c(x,"-p")
    
    return(x)
  }
  cat("Extracting phosphosite motifs...\n")
  #Read fasta file
  fasta <- readAAStringSet(fasta_input)
  names(fasta) <- str_match(names(fasta),
                            "\\|([\\d\\w-]+)")[,2]
  
  #Extract motifs
  dea <- 
    dea %>%
    mutate(motif = map2_chr(prot_id,ptm_site,extract_sequence,fasta = fasta))
  cat("Motif extraction complete!\n")
  
  
  # 3. Create GCT file -----------------------------------------------------------
  
  cat("Creating GCT file for use as ssGSEA input...\n")
  #Extract signlogp values 
  #Collapse peptides to motifs, any motifs matched by more than a single peptide,
  #protein, or site will be concatenated with the separator '|'
  signlogp <- 
    dea %>%
    drop_na(motif) %>%
    group_by(motif,contrast) %>%
    nest() %>%
    mutate(data = map2(data,motif,function(x,y){
      # if(length(unique(x$prot_id))>1){
      #   warn(paste0("\nWarning: Multiple proteins match motif ", y))
      # }
      x <- x %>%
        slice_max(order_by = abs(sign_logp),
                  n=1)
      return(x)
    })) %>%
    unnest(data) %>%
    #Make into wide format
    select(-pval,-fc) %>%
    pivot_wider(id_cols = motif,
                names_from = contrast, 
                values_from = sign_logp,
                unused_fn = ~paste(unique(.),collapse="|",sep = "|"))
  
  mat <-  
    signlogp %>% 
    column_to_rownames(var = "motif") %>%
    select(-prot_id,-ptm_site,-pep_id) %>%
    as.matrix()
  #Fix vectors in first column
  
  rdesc <-  
    signlogp %>%
    select(motif,prot_id,ptm_site,pep_id)%>%
    mutate( id = motif) 
  
  cdesc <-
    data.frame(id = colnames(mat))
  
  full_out_path <- file.path(out_path,paste0(output_prefix,"_gsea"))
  if(!dir.exists(full_out_path)){
    dir.create(full_out_path)
  }
  
  gct_output <- 
    new('GCT',mat =mat,rdesc=rdesc,cdesc=cdesc)
  
  write_gct(gct_output,
            file.path(full_out_path ,"signlogp.gct"), appenddim = F)
  
  cat("GCT file created!\n")
  
  # 4. Run PTM-SEA ------------------------------------------------------------
  cat("Starting PTM-SEA algorithm...\n")
  script.dir <- path.expand("./src/ssGSEA2.0/")
  source(file.path(ssgsea2_dir,"src/ssGSEA2.0.R"))
  current_wd <- getwd() 
  
  #Copy gene set database to working directory
  file.copy(gene_set_database,full_out_path )
  
  setwd(full_out_path )
  #Setting sample.norm.type to none and correl.type to rank uses the actual signed log p-vals
  gsea.out <- ssGSEA2(input.ds = "signlogp.gct", output.prefix = output_prefix, 
                      gene.set.databases = basename(gene_set_database),
                      sample.norm.type = "rank", weight=0, 
                      correl.type="z.score", statistic="Kolmogorov-Smirnov",
                      output.score.type="NES", min.overlap=3, extended.output=T, 
                      global.fdr=F, nperm = 1000)
  setwd(current_wd)
  
  cat("PTM-SEA completed!\n")
  
  return(paste0("PTM-SEA analysis complete. See output in ",full_out_path))
}

