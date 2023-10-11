
##load data
data<- readRDS(".../Harmonized_data_CD68.rds")

### enrichment analysis
source("/Volumes/Bertolaz/DSP_RNAseq/DSP_data_Min_new/Github_codes/function_hallmark_enrichment.R")

directory = "path_to_...Hallmark_2022.txt"

require(readxl)
setwd(directory)
DEGs<- read_xlsx("DEGs_DZ_vs_LZ.xlsx", sheet = 2)

## enrichment UP genes - gsea
ris_enrich_LZ = funz_enrich_gsea( input_genes=DEGs$geni, correction="BH",
                                  dir= directory, pool = rownames(data) )

