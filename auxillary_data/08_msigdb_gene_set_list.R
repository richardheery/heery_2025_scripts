# Create a list with gene symbols for all MSigDB gene sets

# Load required packages
library(msigdb)

# Check latest available version of MsigDB. 2023.1 is most recent version available. 
getMsigdbVersions()

# Get GeneSetCollection for MSigDB for human with gene symbols. Took 2 minutes. 
system.time({msigdb_gsc = getMsigdb(org = 'hs', id = 'SYM', version = '2023.1')})

# Add KEGG gene sets to msigdb_gsc. Took 10 seconds
system.time({msigdb_gsc = appendKEGG(msigdb_gsc, version = "2023.1")})

# List all collections in MSigDB. Description of collections can be found at https://www.gsea-msigdb.org/gsea/msigdb/
# h is 50 hallmark pathways, C1 is genomic position gene sites (chromosome band), C2 is curated gene sets from pathway databases,
# C6 is oncology pathways and c8 is cell type signature gene sets
msigdb_collections = listCollections(msigdb_gsc)
msigdb_collections = setNames(msigdb_collections, msigdb_collections)

# Get gene symbols for all gene sets in collections
msigdb_collections_gene_set_list = lapply(msigdb_collections, function(x) 
  geneIds(subsetCollection(gsc = msigdb_gsc, collection = x)))

# List all available subcollections. Description of subsets can be found at https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#C1
msigdb_subcollections = listSubCollections(msigdb_gsc)
msigdb_subcollections = setNames(msigdb_subcollections, msigdb_subcollections)

# Get gene symbols for all gene sets in subcollections
msigdb_subcollections_gene_set_list = lapply(msigdb_subcollections, function(x) 
  geneIds(subsetCollection(gsc = msigdb_gsc, subcollection = x)))

# Put collections in alphabetical order
msigdb_subcollections_gene_set_list = msigdb_subcollections_gene_set_list[sort(names(msigdb_subcollections_gene_set_list))]

# Get names of all canonical pathways (CP) subcollections
cp_subcollections = grep("CP:", names(msigdb_subcollections_gene_set_list), value = T)

# Remove the name of the database from the pathway names from CP subcollections.
for(cp in cp_subcollections){
  names(msigdb_subcollections_gene_set_list[[cp]]) = 
    substring(names(msigdb_subcollections_gene_set_list[[cp]]), first = stringr::str_locate(names(msigdb_subcollections_gene_set_list[[cp]]), "_(.*?)$")[, 1] + 1)
}

# Get names of all GO subcollections
go_subcollections = grep("GO:", names(msigdb_subcollections_gene_set_list), value = T)

# Remove the name of the database from the pathway names from GO subcollections.
for(go in go_subcollections){
  names(msigdb_subcollections_gene_set_list[[go]]) = 
    substring(names(msigdb_subcollections_gene_set_list[[go]]), first = stringr::str_locate(names(msigdb_subcollections_gene_set_list[[go]]), "_(.*?)$")[, 1] + 1)
}

# Select some collections from msigdb_collections_gene_set_list and msigdb_subcollections_gene_set_list
selected_collections = c("h", "c1", "c6", "c8")
selected_subcollections = c("CP:BIOCARTA", "CP:KEGG", "CP:PID", "CP:REACTOME", "CP:WIKIPATHWAYS", "GO:BP", "GO:CC", "GO:MF", "TFT:GTRD")

# Add some collections from msigdb_collections_gene_set_list to msigdb_subcollections_gene_set_list
msigdb_complete_gene_set_list = c(msigdb_collections_gene_set_list[selected_collections], msigdb_subcollections_gene_set_list[selected_subcollections])

# Remove HALMMARK from start of Hallmark pathways
names(msigdb_complete_gene_set_list$h) = 
    substring(names(msigdb_complete_gene_set_list$h), first = stringr::str_locate(names(msigdb_complete_gene_set_list$h), "_(.*?)$")[, 1] + 1)

# Remove _TARGET_GENES from end of TFT:GTRD pathways
names(msigdb_complete_gene_set_list$`TFT:GTRD`) = gsub("_TARGET_GENES", "", names(msigdb_complete_gene_set_list$`TFT:GTRD`))

# Save msigdb_complete_gene_set_list
saveRDS(msigdb_complete_gene_set_list, "msigdb_complete_gene_set_list.rds")