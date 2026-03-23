# file name: GeNetME_Func
#
current_dir<-getwd() # current directory
set.seed(12345) # to make results reproducible 
#
#-------------------------------------------------------------------------------
# Packages and files
#-------------------------------------------------------------------------------
#
# Load libraries
#
library(clusterProfiler)
library(data.table)
library(DOSE)
library(enrichplot)
library(ggplot2)
library(GOplot)
library(httr)
library(igraph)
library(jsonlite)
library(magick)
library(Matrix)
library(org.Hs.eg.db)
library(pathview)
library(ReactomePA)
library(readxl)
library(stringr)
library(TissueEnrich)
#
#-------------------------------------------------------------------------------
# Parameters
#-------------------------------------------------------------------------------
#
STRING.co<-0.4 # cut-off for gene interaction in STRING API and STRING data
#
#-------------------------------------------------------------------------------
# Create output folder, if absent
#-------------------------------------------------------------------------------
#
folder_path<-file.path(current_dir,"Output")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
folder_path<-file.path(current_dir,"Output/ORA")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
} 
folder_path<-file.path(current_dir,"Output/ORA/KEGG")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
} 
folder_path<-file.path(current_dir,"Output/ORA/Reactome")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
folder_path<-file.path(current_dir,"Output/ORA/GO")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
#
#-------------------------------------------------------------------------------
# Build (if necessary) and load STRING database
#-------------------------------------------------------------------------------
#
folder_path<-file.path(current_dir,"Data")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
folder_path<-file.path(current_dir,"Data/STRING")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
#
# Load data for all species, then keep only H. sapiens
#
url<-paste0("https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz")
destfile<-"9606.protein.links.v12.0.txt.gz"
file_path<-file.path(current_dir,"Data/STRING/",destfile)
if(!file.exists(file_path)) {
  print("Downloading STRING database...")
  RETRY(
    verb = "GET",
    url = url,
    write_disk(file_path, overwrite = TRUE),
    times = 5,           # up to 5 attempts
    pause_min = 5,       # wait 5s between attempts
    terminate_on = c(404, 403) # don't retry on these errors
  )
}   
url<-paste0("https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz")
destfile<-"9606.protein.info.v12.0.txt.gz"
file_path<-file.path(current_dir,"Data/STRING/",destfile)
if(!file.exists(file_path)) {
  print("Downloading STRING database...")
  RETRY(
    verb = "GET",
    url = url,
    write_disk(file_path, overwrite = TRUE),
    times = 5,           # up to 5 attempts
    pause_min = 5,       # wait 5s between attempts
    terminate_on = c(404, 403) # don't retry on these errors
  )
}   
#
print("Loading STRING database...")
#
file.name<-paste0(current_dir,"/Data/STRING/9606.protein.links.v12.0.txt.gz")
STRING.matrix<-read.csv(file.name,sep=" ") # PPI scores
#
file.name<-paste0(current_dir,"/Data/STRING/9606.protein.info.v12.0.txt.gz")
STRING.names<-fread(file.name,sep="\t") # STRING preferred names
colnames(STRING.names)[1]<-"string_protein_id"
#
#-----------------------------------------------------------------------------
# Download ME/CFS module form Zhang S. 2025 
# We use supplementary table 2.
#-----------------------------------------------------------------------------
#
url<-paste0("https://www.medrxiv.org/content/medrxiv/early/2025/05/11/2025.04.15.25325899/DC9/embed/media-9.xlsx")
destfile<-"media-9.xlsx"
file_path<-file.path(current_dir,"Data",destfile)
if(!file.exists(file_path)) {
  print("Downloading data from Zhang S. et al. 2025 (https://pmc.ncbi.nlm.nih.gov/articles/PMC12047926/)...")
  RETRY(
    verb = "GET",
    url = url,
    write_disk(file_path, overwrite = TRUE),
    times = 5,           # up to 5 attempts
    pause_min = 5,       # wait 5s between attempts
    terminate_on = c(404, 403) # don't retry on these errors
  )
}  
#
#-----------------------------------------------------------------------------
# Download ME/CFS module form Sardell JM 2025
# We use supplementary table 3.
#-----------------------------------------------------------------------------
#
url<-paste0("https://www.medrxiv.org/content/medrxiv/early/2025/12/03/2025.12.01.25341362/DC2/embed/media-2.xlsx")
destfile<-"media-2.xlsx"
file_path<-file.path(current_dir,"Data",destfile)
if(!file.exists(file_path)) {
  print("Downloading data from Sardell JM. et al. 2025 (https://www.medrxiv.org/content/10.64898/2025.12.01.25341362v2)...")
  RETRY(
    verb = "GET",
    url = url,
    write_disk(file_path, overwrite = TRUE),
    times = 5,           # up to 5 attempts
    pause_min = 5,       # wait 5s between attempts
    terminate_on = c(404, 403) # don't retry on these errors
  )
}   
#
#-------------------------------------------------------------------------------
# Build (if necessary) and load NCBI database
#-------------------------------------------------------------------------------
#
folder_path<-file.path(current_dir,"Data/NCBI")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
} 
url <- "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz"
destfile <- "gene_info.gz"
file_path <- file.path(current_dir, "Data/NCBI/", destfile)
if (!file.exists(file_path) & !file.exists("Data/NCBI/human_gene_info.gz")) {
  print("Downloading NCBI database...")
  RETRY(
    verb = "GET",
    url = url,
    write_disk(file_path, overwrite = TRUE),
    times = 5,
    pause_min = 5,
    terminate_on = c(404, 403)
  )
}
#
file.name<-"Data/NCBI/human_gene_info.gz"
if (!file.exists(file.name)) {
  print("Editing NCBI database...")
  file.name<-"Data/NCBI/gene_info.gz"
  NCBI.names<-fread(file.name,sep="\t") 
  colnames(NCBI.names)[1]<-"tax_id"
  NCBI.names<-subset.data.frame(NCBI.names,tax_id==9606)
  NCBI.names<-subset.data.frame(NCBI.names,select=c("GeneID","Symbol","Synonyms","description"))
  colnames(NCBI.names)<-c("NCBI.id","name","Synonyms","description")
  file.name<-"Data/NCBI/human_gene_info.gz"
  fwrite(NCBI.names,file.name,sep="\t") 
  file.remove("Data/NCBI/gene_info.gz")
} else {
  print("Loading NCBI database...")
  NCBI.names<-fread(file.name,sep="\t")  
}
gc()
#
#-------------------------------------------------------------------------------
# Convert gene symbol to NCBI ID using local data base (not API)
#-------------------------------------------------------------------------------
#
Symbol2NCBI.db<-function(gene.symbol) {
  #
  # Search in "name" first
  #
  df<-subset.data.frame(NCBI.names,name==gene.symbol)
  if (nrow(df)>0) {
    return(as.character(df$NCBI.id[1]))
  } else { # Search among synonyms
    index<-which(grepl(paste0("\\b",gene.symbol,"\\b"),NCBI.names$Synonyms)) 
    if (length(index)==0) {
      return(NA)
    } else {
      return(as.character(NCBI.names$NCBI.id[index[1]]))
    }
  }
}
#
#-------------------------------------------------------------------------------
# This function asks STRING API its preferred name for a gene
#-------------------------------------------------------------------------------
#
STRING.name<-function(gene) {
  #
  species_id<-9606 # Homo sapiens
  base_url<-"https://string-db.org/api/json/get_string_ids"
  identifiers<-paste(gene,collapse="%0d")
  response<-httr::GET(base_url,query=list(identifiers=identifiers,species=species_id)) 
  #
  if (status_code(response)==200) {
    response_content<-rawToChar(response$content)
    Encoding(response_content)<-"UTF-8"
    data<-jsonlite::fromJSON(response_content) # parse
    if (is.data.frame(data)) {
      preferred.name<-data$preferredName
    } else {
      preferred.name<-NA
    }
  } else {
    stop("Failed to retrieve data. Please check the NCBI IDs and your internet connection.")
  }
  return(preferred.name)
}
#
#-------------------------------------------------------------------------------
# This function find genes that interact with the input gene, 
# using STRING database (not API)
#-------------------------------------------------------------------------------
#
STRING2<-function(gene) {
  #
  # Find all the interacting genes
  #
  df<-subset.data.frame(STRING.names,preferred_name==gene)
  gene<-df$string_protein_id
  df<-subset.data.frame(STRING.matrix,protein1==gene)
  df<-subset.data.frame(df,combined_score>=STRING.co*1000)
  #
  # edit the output
  #
  df<-df[,-1]
  colnames(df)<-c("name","score")
  #
  if (nrow(df)>0) {
    for (i in 1:nrow(df)) {
      temp.df<-subset.data.frame(STRING.names,string_protein_id==df$name[i])
      df$name[i]<-temp.df$preferred_name
    }
    df$score<-df$score/1000
    return(df)
  } else {
    return(NA)
  }
}
#
#-------------------------------------------------------------------------------
# For genes that appears more than once, this function generate a single row with
# all the information and with the highest score associated to each gene
#-------------------------------------------------------------------------------
#  
ListCollapse <- function(all.genes.zero) {
  unique.genes <- unique(all.genes.zero$NCBI.id)
  all.genes.unique <- all.genes.zero[1:length(unique.genes), ]
  for (i in 1:length(unique.genes)) {
    temp <- subset.data.frame(all.genes.zero, NCBI.id == unique.genes[i])
    all.genes.unique$NCBI.id[i]    <- temp$NCBI.id[1]
    all.genes.unique$name[i]       <- temp$name[1]
    all.genes.unique$list.count[i] <- length(unique(temp$list.name))
    all.genes.unique$list.name[i]  <- paste0(unique(temp$list.name), collapse="/")
  }
  return(all.genes.unique)
}
#
#-------------------------------------------------------------------------------
# This function build adjacency matrix from a list of genes using STRING database 
# (not API). 
#-------------------------------------------------------------------------------
#
GeneMatrix<-function(all.genes) {
  #
  # build a STRING database with only genes of interest and with score above STRING.co
  #
  df1<-STRING.names
  colnames(df1)[2]<-"name"
  df.names<-merge(all.genes,df1,by="name",all.y=F)
  df.names<-subset.data.frame(df.names,select=c("name","string_protein_id"))
  dt.names<-as.data.table(df.names) # data tables should be faster
  remove(df1,df.names)
  #
  df.matrix<-subset.data.frame(STRING.matrix,combined_score>=STRING.co*1000)
  df.matrix<-subset.data.frame(df.matrix,protein1%in%dt.names$string_protein_id)
  df.matrix<-subset.data.frame(df.matrix,protein2%in%dt.names$string_protein_id)
  dt.matrix<-as.data.table(df.matrix) # data tables should be faster
  remove(df.matrix)
  #
  # build gene.matrix
  #
  nodes<-dt.names$string_protein_id
  i<-match(dt.matrix$protein1,nodes)
  j<-match(dt.matrix$protein2,nodes)
  #
  Msp<-sparseMatrix(
    i = i, j = j, x = dt.matrix$combined_score/1000,
    dims = c(length(nodes), length(nodes)),
    dimnames = list(nodes, nodes)
  )
  #
  gene.matrix<-as.matrix(Msp)
  #
  id_to_name <- dt.names$name
  names(id_to_name) <- dt.names$string_protein_id
  rownames(gene.matrix) <- id_to_name[rownames(gene.matrix)]
  colnames(gene.matrix) <- id_to_name[colnames(gene.matrix)]
  #
  # Test symmetry and end
  #
  test<-isSymmetric(gene.matrix)
  if (!test) {
    stop("The adjacency matrix is not symmetric; there is an error!")
  } else {
    return(gene.matrix)  
  }
}
#
#-------------------------------------------------------------------------------
# This function performs Over-representation Analysis over several databases
#-------------------------------------------------------------------------------
#
ORA.fun<-function(all.genes,top.results,list.number) {
  #
  #-------------------------------------------------------------------------------
  # Over-representation Analysis with KEGG (Kyoto Encyclopedia of Genes and Genomes - 
  # biological pathways)
  #-------------------------------------------------------------------------------
  #
  KEGG_results<-enrichKEGG(gene=all.genes$NCBI.id,organism='hsa',universe=myuniverse)
  KEGG.ORA<-KEGG_results@result
  KEGG.ORA<-KEGG.ORA[,3:ncol(KEGG.ORA)]
  KEGG.ORA$method<-rep("KEGG",nrow(KEGG.ORA))
  #
  matrix.result<-KEGG.ORA
  matrix.result<-matrix.result[order(matrix.result$p.adjust),]
  if (nrow(matrix.result)>top.results) matrix.result<-matrix.result[1:top.results,]
  KEGG.ORA<-matrix.result
  #
  #-------------------------------------------------------------------------------
  # Over-representation Analysis with GO (Gene Ontology - functional role of genes within 
  # biological processes, molecular functions, and cellular components)
  #-------------------------------------------------------------------------------
  #
  GO_results<-enrichGO(gene=all.genes$NCBI.id,OrgDb=org.Hs.eg.db,universe=myuniverse,ont="CC")
  GO.ORA<-GO_results@result
  GO.ORA$method<-rep("GO",nrow(GO.ORA))
  #
  matrix.result<-GO.ORA
  matrix.result<-matrix.result[order(matrix.result$p.adjust),]
  if (nrow(matrix.result)>top.results) matrix.result<-matrix.result[1:top.results,]
  GO.ORA<-matrix.result
  #
  #-------------------------------------------------------------------------------
  # Over-representation Analysis using Reactome 
  #-------------------------------------------------------------------------------
  #
  pathway_results<-enrichPathway(gene=all.genes$NCBI.id,organism="human",universe=myuniverse)
  Rea.ORA<-pathway_results@result
  Rea.ORA$method<-rep("Reactome",nrow(Rea.ORA))
  #
  matrix.result<-Rea.ORA
  matrix.result<-matrix.result[order(matrix.result$p.adjust),]
  if (nrow(matrix.result)>top.results) matrix.result<-matrix.result[1:top.results,]
  Rea.ORA<-matrix.result
  #
  #-------------------------------------------------------------------------------
  # Over-representation Analysis using Disease Ontology (DO)
  #-------------------------------------------------------------------------------
  #
  DO.ORA<-enrichDO(gene=all.genes$NCBI.id,universe=myuniverse)
  DO.ORA<-DO.ORA@result
  DO.ORA$method<-rep("DO",nrow(DO.ORA)) 
  #
  matrix.result<-DO.ORA
  matrix.result<-matrix.result[order(matrix.result$p.adjust),]
  if (nrow(matrix.result)>top.results) matrix.result<-matrix.result[1:top.results,]
  DO.ORA<-matrix.result
  #
  #-------------------------------------------------------------------------------
  # Merge results from Over-representation Analysis and keep only significant ones
  #-------------------------------------------------------------------------------
  #
  if ("ONTOLOGY"%in%colnames(GO.ORA)) GO.ORA<-GO.ORA[,!names(GO.ORA) %in% "ONTOLOGY"]
  ORA<-rbind(KEGG.ORA,GO.ORA,Rea.ORA,DO.ORA)
  if (nrow(ORA)>0) {
    ORA<-subset.data.frame(ORA,p.adjust<=0.05)
  }
  if (nrow(ORA)>0) {
    ORA<-subset.data.frame(ORA,select=c("ID","Description","GeneRatio",
                                        "pvalue","p.adjust","geneID","method"))
    ORA$gene.name<-rep(NA,nrow(ORA))
    ORA$list.name<-rep(NA,nrow(ORA))
    ORA$list.count<-rep(0,nrow(ORA))
    for (i in 1:nrow(ORA)) {
      vector1<-stringr::str_split(ORA$geneID[i],"/")[[1]]
      vector2<-vector1
      vector3<-vector1
      for (j in 1:length(vector1)) {
        df<-subset.data.frame(all.genes,NCBI.id==as.numeric(vector1[j]))
        vector2[j]<-df$name
        vector3[j]<-df$list.name
      }
      ORA$gene.name[i]<-paste0(unique(vector2),collapse="/")
      ORA$list.name[i]<-paste0(unique(vector3),collapse="/")
      ORA$list.count[i]<-length(vector3)
    }
    for (i in 1:nrow(ORA)) {
      vector1<-unique(stringr::str_split(ORA$list.name[i],"/")[[1]])
      ORA$list.name[i]<-paste0(vector1,collapse="/")
      ORA$list.count[i]<-length(vector1)
    }
  } 
  #
  # Keep only terms with at least two controbuting gene lists
  #
  ORA<-subset.data.frame(ORA,list.count>=list.number)
  #
  # Write gene symbols in alphabetical order
  #
  if (nrow(ORA)>0) {
    for (i in 1:nrow(ORA)) {
      vec<-sort(str_split(ORA$gene.name[i],"/")[[1]])
      vec<-paste0(vec,collapse="/")
      ORA$gene.name[i]<-vec
    }
  } 
  #
  #-------------------------------------------------------------------------------
  # Plots for GO
  #-------------------------------------------------------------------------------
  #
  if (nrow(GO.ORA)>0) {
    index<-which(GO_results@result$ID%in%ORA$ID)
    GO_results@result<-GO_results@result[index,]
    options(ggrepel.max.overlaps = 1000)
    tiff("Output/ORA/GO/GO_ORA_goplot.tiff",width=20,height=8,units="in",res=600,compression="lzw")
    p<-goplot(GO_results,showCategory=10)
    print(p)
    dev.off()
    #
    # Save a jpeg version too
    #
    img<-image_read("Output/ORA/GO/GO_ORA_goplot.tiff")
    image_write(img,path="Output/ORA/GO/GO_ORA_goplot.jpeg",format="jpeg")
    #
    tiff("Output/ORA/GO/GO_ORA_cnetplot.tiff",width=20,height=20,units="in",res=600,compression="lzw")
    GO_results@result$geneID <- ORA$gene.name[match(GO_results@result$ID, ORA$ID)] # use gene symbols for this plot
    p<-cnetplot(GO_results,showCategory=10,layout="fr") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom") +
      ggtitle("Top Enriched GO Terms and Genes")
    print(p)
    dev.off()  
    #
    # Save a jpeg version too
    #
    img<-image_read("Output/ORA/GO/GO_ORA_cnetplot.tiff")
    image_write(img,path="Output/ORA/GO/GO_ORA_cnetplot.jpeg",format="jpeg")
  }
  #
  return(ORA)
}
#
#-------------------------------------------------------------------------------
# perform ORA with respect to tissue expression
#-------------------------------------------------------------------------------
# 
Tissue.ORA<-function(all.genes) {
  gs<-GeneSet(geneIds=all.genes$name,organism='Homo Sapiens',geneIdType=SymbolIdentifier())
  bk<-GeneSet(geneIds=STRING.names$preferred_name,organism='Homo Sapiens',geneIdType=SymbolIdentifier())
  output<-teEnrichment(gs,rnaSeqDataset=1,backgroundGenes=bk) # 1 for HPA, 2 for GTEx
  seEnrichmentOutput<-output[[1]]
  enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput),row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])
  enrichmentOutput$Tissue<-row.names(enrichmentOutput)
  enrichmentOutput$p.adjust<-10^-enrichmentOutput$Log10PValue # it is adjusted!
  #
  # plot
  #
  tiff("Output/ORA/Tissue_ORA.tiff",width=10,height=6,units="in",res=600,compression="lzw")
  p<-ggplot(enrichmentOutput,aes(x=reorder(Tissue,-Log10PValue),y=Log10PValue,
                              label = Tissue.Specific.Genes,fill = Tissue))+
    geom_bar(stat = 'identity')+
    geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +
    labs(x='', y = '-LOG10(P-Value)')+
    theme_bw()+
    theme(legend.position='none')+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title =
            element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          panel.grid.major= element_blank(),panel.grid.minor = element_blank())
  print(p)
  dev.off()
  #
  # Save a jpeg version too
  #
  img<-image_read("Output/ORA/Tissue_ORA.tiff")
  image_write(img,path="Output/ORA/Tissue_ORA.jpeg",format="jpeg")
  #
  return(enrichmentOutput)
}
#
#-------------------------------------------------------------------------------
# Prepare background genes
#-------------------------------------------------------------------------------
#
file.name<-"Data/My_Universe.tsv"
if (!file.exists(file.name)) {
  print("Building background for Over-representation analysis. This may take a while...")
  STRING.names$NCBI.id<-rep(NA,nrow(STRING.names))
  for (i in 1:nrow(STRING.names)) {
    STRING.names$NCBI.id[i]<-Symbol2NCBI.db(STRING.names$preferred_name[i])
  }
  myuniverse<-as.character(STRING.names$NCBI.id)  
  myuniverse<-na.omit(myuniverse)
  write.table(myuniverse,file=file.name,quote=F,row.names=F,col.names=T,sep="\t")
} else {
  myuniverse<-fread(file.name)
  myuniverse<-as.character(myuniverse$x)
}