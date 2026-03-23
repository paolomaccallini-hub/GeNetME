# file name: GeNetME_Main
#
source("GeNetME_Func.R")
#
#-----------------------------------------------------------------------------
# Read ME/CFS module form Zhang S. 2025 (https://pmc.ncbi.nlm.nih.gov/articles/PMC12047926/)
# We use supplementary table 2.
#-----------------------------------------------------------------------------
#
Module<-read_xlsx("Data/media-9.xlsx") 
Module<-subset.data.frame(Module,q_value<0.02) # 115 genes associated with ME/CFS  
#
# Use STRING preferred name for genes
#
gene.wo.STRING<-0
index<-c()
Module$name<-rep(NA,nrow(Module))
for (i in 1:nrow(Module)) {
  new.name<-STRING.name(Module$Gene[i])
  if (!is.na(new.name)) {
    Module$name[i]<-new.name
  } else if (is.na(new.name)) {
    print(paste(Module$Gene[i],"is not present in STRING's database"))
    gene.wo.STRING<-gene.wo.STRING+1
    index<-c(index,i)
  }
}
if (length(index)>0) Module<-Module[-index,] # only genes in STRING
#
# Edit
#
Module$NCBI.id<-rep(NA,nrow(Module))
Module$list.name<-rep("Zhang",nrow(Module))
#
# Add NCBI id
#
for (i in 1:nrow(Module)) {
  Module$NCBI.id[i]<-Symbol2NCBI.db(Module$name[i])
}
Zhang_module<-as.data.frame(Module[,c("name","NCBI.id","list.name")])
#
#-----------------------------------------------------------------------------
# Read ME/CFS module form Sardell JM 2025 (https://www.medrxiv.org/content/10.64898/2025.12.01.25341362v2)
# We use supplementary table 3.
#-----------------------------------------------------------------------------
#
Module<-read_xlsx("Data/media-2.xlsx",sheet=3,skip=2) 
Module<-data.frame(Gene=Module$`Gene name`)
#
# Use STRING preferred name for genes
#
gene.wo.STRING<-0
index<-c()
Module$name<-rep(NA,nrow(Module))
for (i in 1:nrow(Module)) {
  new.name<-STRING.name(Module$Gene[i])
  if (!is.na(new.name)) {
    Module$name[i]<-new.name
  } else if (is.na(new.name)) {
    print(paste(Module$Gene[i],"is not present in STRING's database"))
    gene.wo.STRING<-gene.wo.STRING+1
    index<-c(index,i)
  }
}
if (length(index)>0) Module<-Module[-index,] # only genes in STRING
#
# Edit
#
Module$NCBI.id<-rep(NA,nrow(Module))
Module$list.name<-rep("PL",nrow(Module))
#
# Add NCBI id
#
for (i in 1:nrow(Module)) {
  Module$NCBI.id[i]<-Symbol2NCBI.db(Module$name[i])
}
PL_module<-as.data.frame(Module[,c("name","NCBI.id","list.name")])
#
#-------------------------------------------------------------------------------
# Merge 
#-------------------------------------------------------------------------------
#
all.genes.zero<-rbind(Zhang_module,PL_module)
#
#-------------------------------------------------------------------------------
# Build a list with a single row for each gene
#-------------------------------------------------------------------------------
#
all.genes<-ListCollapse(all.genes.zero)
all.genes<-all.genes[,c("name","NCBI.id","list.name","list.count")]
#
#-------------------------------------------------------------------------------
# Build the adjacency matrix associated with the merged gene list
#-------------------------------------------------------------------------------
#
gene.matrix<-GeneMatrix(all.genes)
#
#-------------------------------------------------------------------------------
# Add components, degree, and scaled degree (score)
#-------------------------------------------------------------------------------
#
graph<-graph_from_adjacency_matrix(gene.matrix,mode="undirected",weighted=TRUE)
all.genes<-all.genes[match(rownames(gene.matrix),all.genes$name), ] 
#
com<-igraph::components(graph)
all.genes$comp<-com$membership
#
deg<-igraph::degree(graph)
all.genes$degree<-deg
#
strength <- igraph::strength(graph)
all.genes$strength <- strength
#
all.genes$score<-(all.genes$strength-min(all.genes$strength))/(max(all.genes$strength)-min(all.genes$strength))
#
#-------------------------------------------------------------------------------
# Save relevant results as RSD file
#-------------------------------------------------------------------------------
#
Disease_results_list<-list(all.genes,gene.matrix)
saveRDS(Disease_results_list,file="Output/Disease_results_list.rds")
write.table(all.genes,file="All_genes.tsv",sep="\t",row.names=FALSE,quote=FALSE)
#
#-------------------------------------------------------------------------------
# Read relevant results as RSD file
#-------------------------------------------------------------------------------
#
Disease_results_list<-readRDS(file="Output/Disease_results_list.rds")
all.genes<-Disease_results_list[[1]]
gene.matrix<-Disease_results_list[[2]]
#
#-------------------------------------------------------------------------------
# Build the graph associated with the Merged Gene List (MGL) and plot it
#-------------------------------------------------------------------------------
#
graph<-graph_from_adjacency_matrix(gene.matrix,mode="undirected",weighted=TRUE)
#
# Color the genes according to the corresponding Expanded Gene List (EGL)
#
colors<-hcl.colors(length(unique(all.genes$list.name)),palette="Dark 3",alpha=1)
vertex_colors<-c()
for (i in 1:nrow(all.genes)) {
  index<-which(unique(all.genes$list.name)==all.genes$list.name[i])
  vertex_colors[i]<-colors[index]
}
#
# Set the color of node labels
#
node.col<-rep("black",nrow(all.genes))
#
# Plot the image
#
tiff("Output/All_genes_graph.tiff",width=10,height=10,units="in",res=600,compression="lzw")
plot(graph,
     layout=layout_with_fr(graph,niter=30000,grid="nogrid",dim=2),
     vertex.size=3,
     vertex.label.cex=0.2,
     vertex.frame.color="black",
     vertex.color=vertex_colors,
     vertex.label.color=node.col,
     asp=1,
)
legend("topright",
       legend = unique(all.genes$list.name),
       col = "black",
       pt.bg = colors,
       pch = 21,
       pt.cex = 1.5,
       bty = "n")
dev.off()
#
# Save a jpeg version too
#
img<-image_read("Output/All_genes_graph.tiff")
image_write(img, path = "Output/All_genes_graph.jpeg", format = "jpeg")
#
# Save graph for cytoscape
#
edges<-as.data.frame(as_edgelist(graph))
weights<-E(graph)$weight
edges_df<-data.frame(source=edges[,1],target=edges[,2],interaction="interacts_with",
                     weight=weights)
file_name="Output/All_genes_cytoscape.tsv"
write.table(edges_df,file=file_name,sep="\t",row.names=FALSE,quote=FALSE)
#
#-------------------------------------------------------------------------------
# Select Disease module 
#-------------------------------------------------------------------------------
#
Disease_module<-all.genes
# Disease_module<-subset.data.frame(all.genes,list.name=="Zhang"|list.name="Zhang/PL")
# Disease_module<-subset.data.frame(all.genes,list.name=="PL"|list.name="Zhang/PL")
#
#-------------------------------------------------------------------------------
# Over-representation analysis (ORA)  
#-------------------------------------------------------------------------------
#
top.results<-100 # max number of terms included in the final result
list.number<-1 # the number of gene lists required for a term to be included
ORA<-ORA.fun(Disease_module,top.results,list.number)
file.name<-"Output/ORA/ORA_Disease_Module.tsv"
write.table(ORA,file=file.name,quote=F,row.names=F,col.names=T,sep="\t")
#
#-------------------------------------------------------------------------------
# ORA visualization
#-------------------------------------------------------------------------------
#
top.n<-5 # how many top terms?
for (ORA.type in c("hsa","GO","R-HSA","DO")) {
  myORA<-subset.data.frame(ORA,grepl(ORA.type,ORA$ID))
  if (nrow(myORA) < top.n) next  # add this
  terms<-myORA[,c("ID","Description","gene.name","p.adjust")]
  terms$Category<-rep("BP",nrow(terms))
  colnames(terms)<-c("ID","Term","Genes","adj_pval","Category")
  terms<-subset.data.frame(terms,select=c("Category","ID","Term","Genes","adj_pval"))
  #
  for (i in 1:nrow(terms)) {
    terms$Genes[i]<-paste0(str_split(terms$Genes[i],"/")[[1]],collapse=",")
  }
  all_names<-unique(strsplit(paste0(terms$Genes,collapse=","),",")[[1]])
  index<-which(Disease_module$name%in%all_names)
  gene_expression<-Disease_module[index,]
  gene_expression<-subset.data.frame(gene_expression,select=c("name","score"))
  colnames(gene_expression)<-c("ID","logFC")
  gene_expression$logFC<-as.numeric(gene_expression$logFC)
  circ<-circle_dat(terms,gene_expression)
  chord<-chord_dat(circ,gene_expression,terms$Term[1:top.n])
  if (ORA.type=="hsa") {
    title.str<-"Output/ORA/KEGG enrichment"
  } else if (ORA.type=="GO") {
    title.str<-"Output/ORA/GO enrichment"
  } else if (ORA.type=="R-HSA") {
    title.str<-"Output/ORA/Reactome enrichment"
  } else {
    title.str<-"Output/ORA/DO enrichment"
  }
  tiff(paste0(title.str,".tiff"),width=35,height=35,units="cm",res=600,compression="lzw")
  p<-GOChord(chord,space=0.001,title=title.str, 
             gene.order="logFC",
             gene.size = 5,
             process.label = 10,
             lfc.col=c('red','orange','yellow'))
  print(p)
  dev.off()
  #
  # Save a jpeg version too
  #
  img<-image_read(paste0(title.str,".tiff"))
  image_write(img,path=paste0(title.str,".jpeg"),format="jpeg")
}
#
#-------------------------------------------------------------------------------
# ORA Pathway visualization (KEGG)
#-------------------------------------------------------------------------------
#
pathway <- subset.data.frame(ORA, grepl("hsa", ORA$ID))$ID
pathway <- pathway[1:min(10, length(pathway))]
#
for (i in 1:length(pathway)) {
  p<-try(pathview(gene.data=as.character(Disease_module$NCBI.id),
              pathway.id=pathway[i],
              species="hsa",
              out.suffix="Disease_ORA",
              gene.idtype="entrez",
              kegg.native=T,
              node.sum="max",
              bins=list(gene=2),
              res=300,
              width=3000,
              height=2500))
  if (class(p)!="try-error") {
    p
  } else {
    try(dev.off())
  }
}
#
# Move to desired directory
#
files<-list.files(pattern="^hsa.*",full.names=T)
for (i in 1:length(files)) {
  if (grepl("Disease_ORA",files[i])) {
    file.rename(files[i],file.path("Output/ORA/KEGG",files[i]))   
  } else {
    file.remove(files[i])
  }
}
#
#-------------------------------------------------------------------------------
# ORA Pathway visualization (Reactome)
#-------------------------------------------------------------------------------
#
pathway<-subset.data.frame(ORA,grepl("R-HSA",ORA$ID))$ID
for (i in 1:length(pathway)) {
  entrez_ids<-subset.data.frame(ORA,ID==pathway[i])$geneID
  entrez_ids<-str_split(entrez_ids,"/")[[1]]
  writeLines(entrez_ids,paste0("Output/ORA/Reactome/",pathway[i],".txt"))
}
#
# submit this file at https://reactome.org/PathwayBrowser/#TOOL=AT
#
#-------------------------------------------------------------------------------
# ORA for tissue-specific gene expression
#-------------------------------------------------------------------------------
#
TORA<-Tissue.ORA(Disease_module)
file.name<-"Output/ORA/TORA_Disease_Module.tsv"
write.table(TORA,file=file.name,quote=F,row.names=F,col.names=T,sep="\t")
