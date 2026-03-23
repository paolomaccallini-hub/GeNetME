# file name: GeNetME_Main
#
source("GeNetME_Func.R")
#
#-------------------------------------------------------------------------------
# Parameters
#-------------------------------------------------------------------------------
#
STRING.co<-0.4 # cut-off for gene interaction in STRING API and STRING data
#
#-----------------------------------------------------------------------------
# Read ME/CFS module form Zhang S. 2025 (https://pmc.ncbi.nlm.nih.gov/articles/PMC12047926/)
# We use supplementary table 2.
#-----------------------------------------------------------------------------
#
Module<-read_xlsx("Data/Zhang_2025.xlsx") 
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
Module<-read_xlsx("Data/PL_2025.xlsx",sheet=3,skip=2) 
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
#-----------------------------------------------------------------------------
# Read ME/CFS module form Sardell JM 2025 (https://www.medrxiv.org/content/10.64898/2025.12.01.25341362v2)
# We use supplementary table 3.
#-----------------------------------------------------------------------------
#
Proteo_module<-read_xlsx("Data/Proteomics.xlsx",sheet=1) # total effect
#
# Use STRING preferred name for genes and remove genes not in STRING
#
index<-c()
gene.wo.STRING<-0
for (i in 1:nrow(Proteo_module)) {
  new.name<-STRING.name(Proteo_module$name[i])
  if (!is.na(new.name)) {
    Proteo_module$name[i]<-new.name
  } else if (is.na(new.name)) {
    print(paste(Proteo_module$name[i],"is not present in STRING's database"))
    gene.wo.STRING<-gene.wo.STRING+1
    index<-c(index,i)
  }
}
if (length(index)>0) {
  Results<-Proteo_module[-index,] # remove genes not in STRING
} else {
  Results<-Proteo_module
}
#
write_xlsx(Results,"Data/UKBB_proteomics.xlsx",col_names=T)
#
#-------------------------------------------------------------------------------
# Merge 
#-------------------------------------------------------------------------------
#
all.genes.zero<-rbind(Zhang_module,PL_module)
Exper_list<-all.genes.zero
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
# Build a graph with the expanded gene lists 
#-------------------------------------------------------------------------------
#
list.matrix<-ListMatrix(all.genes)
#
#-------------------------------------------------------------------------------
# Add components, degree, and scaled degree (score)
#-------------------------------------------------------------------------------
#
graph<-graph_from_adjacency_matrix(gene.matrix,mode="undirected",weighted=TRUE)
#
com<-igraph::components(graph)
all.genes$comp<-com$membership
#
deg<-igraph::degree(graph)
all.genes$degree<-deg
#
all.genes$score<-(all.genes$degree-min(all.genes$degree))/(max(all.genes$degree)-min(all.genes$degree))
#
#-------------------------------------------------------------------------------
# Save relevant results as RSD file
#-------------------------------------------------------------------------------
#
Disease_results_list<-list(Exper_list,all.genes,gene.matrix,list.matrix)
saveRDS(Disease_results_list,file="PF_output/Disease_results_list.rds")
file.name<-"PF_output/all_genes.tsv"
write.table(all.genes,file=file.name,quote=F,row.names=F,col.names=T,sep="\t")
#
#-------------------------------------------------------------------------------
# Read relevant results as RSD file
#-------------------------------------------------------------------------------
#
Disease_results_list<-readRDS(file="PF_output/Disease_results_list.rds")
Exper_list<-Disease_results_list[[1]]
all.genes<-Disease_results_list[[2]]
gene.matrix<-Disease_results_list[[3]]
list.matrix<-Disease_results_list[[4]]
#
#-------------------------------------------------------------------------------
# Build the graph associated with the Merged Gene List (MGL) and plot it
#-------------------------------------------------------------------------------
#
graph<-graph_from_adjacency_matrix(gene.matrix,mode="undirected",weighted=TRUE)
#
# Color the genes according to the corresponding Expanded Gene List (EGL)
#
all.genes<-all.genes[match(rownames(gene.matrix),all.genes$name), ] # correct the order!
colors<-hcl.colors(nrow(list.matrix),palette="Dark 3",alpha=1)
vertex_colors<-c()
for (i in 1:nrow(all.genes)) {
  index<-which(rownames(list.matrix)==all.genes$list.name[i])
  if (length(index)>0) {
    vertex_colors[i]<-colors[index]
  }
}
index<-which(is.na(vertex_colors))
vertex_colors[index]<-"white"
#
# Set the color of nodel labels
#
node.col<-rep("black",nrow(all.genes))
#
# Plot the image
#
tiff("PF_output/All_genes_graph.tiff",width=10,height=10,units="in",res=600,compression="lzw")
plot(graph,
     layout=layout_with_fr(graph,niter=30000,grid="nogrid",dim=2),
     vertex.size=3,
     vertex.label.cex=0.2,
     vertex.frame.color="black",
     vertex.color=vertex_colors,
     vertex.label.color=node.col,
     asp=1,
)
legend("topright",                       # or "bottomleft", etc.  
       legend = rownames(list.matrix),   # group names
       col = colors,                     # corresponding colors
       pch = 21,                         # filled circle
       pt.bg = colors,                   # fill color
       pt.cex = 1.5,                     # size of points
       bty = "n")                        # no box around legend
dev.off()
#
# Save a jpeg version too
#
img<-image_read("PF_output/All_genes_graph.tiff")
image_write(img, path = "PF_output/All_genes_graph.jpeg", format = "jpeg")
#
# Save graph for cytoscape
#
edges<-as.data.frame(as_edgelist(graph))
weights<-E(graph)$weight
edges_df<-data.frame(source=edges[,1],target=edges[,2],interaction="interacts_with",
                     weight=weights)
file_name="PF_output/All_genes_cytoscape.tsv"
write.table(edges_df,file=file_name,sep="\t",row.names=FALSE,quote=FALSE)
#
#-------------------------------------------------------------------------------
# Plot the graph associated with the expanded gene lists (EGLs)
#-------------------------------------------------------------------------------
#
# For each overlap different from zero, calculate a p value by hypergeometric test
#
NL<-nrow(list.matrix)
NB<-nrow(STRING.names) # total background (white and black balls: m+n)
list.matrix.pvalue<-list.matrix
for (i in 1:(NL-1)) {
  for (j in (i+1):NL) {
    list.matrix[i,j]<-as.numeric(list.matrix[i,j])
    if (list.matrix[i,j]==0) {
      list.matrix.pvalue[i,j]<-1
    } else {
      N1<-list.matrix[i,i] # size of EGL_i (number of white balls, m)
      N2<-list.matrix[j,j] # size of EGL_j (number of extracted balls, k)
      NR<-list.matrix[i,j] # size of overlap (number of white balls extracted, q) 
      list.matrix.pvalue[i,j]<-as.numeric(phyper(q=NR-1,m=N1,n=NB-N1,k=N2,lower.tail=F)) # P(X>NR-1)=P(X>=NR)
      list.matrix.pvalue[i,j]<-as.numeric(formatC(list.matrix.pvalue[i,j],format="e",digits=2))
    }
    list.matrix.pvalue[j,i]<-list.matrix.pvalue[i,j]
  }
}
for (i in 1:NL) {
  list.matrix.pvalue[i,i]<-NA
}
#
# Generate graph from list.matrix
#
temp.matrix<-list.matrix
diag(temp.matrix)<-0 # this is necessary to avoid loops
graph<-graph_from_adjacency_matrix(temp.matrix,mode="undirected",weighted=T)
#
# Generate graph from list.matrix.pvalue
#
temp.matrix<-list.matrix.pvalue
diag(temp.matrix)<-0 # this is necessary to avoid loops
graph.pvalue<-graph_from_adjacency_matrix(temp.matrix,mode="undirected",weighted=T)
#
# Prepare colors for the EGLs
#
colors<-hcl.colors(nrow(list.matrix),palette="Dark 3",alpha=1)
#
# Prepare width for edges
#
edge.w<-E(graph)$weight
edge.w.scaled<-(edge.w-min(edge.w))/(max(edge.w)-min(edge.w))*8+5
#
# Plot it
#
tiff("PF_output/All_genes_list_matrix.tiff",width=10,height=10,units="in",res=600,compression="lzw")
plot(graph,
     vertex.size=30,
     vertex.label.cex=2.,
     vertex.label=paste0(rownames(list.matrix),"\n(size = ",diag(list.matrix),")"),
     vertex.frame.color="black",
     vertex.color=colors,
     edge.label=paste0("overlap: ",E(graph)$weight,"\n p = ",E(graph.pvalue)$weight),
     edge.label.cex=1.5,
     edge.width=edge.w.scaled,
     asp=1,
)
dev.off()
#
#-------------------------------------------------------------------------------
# Select Disease module by components
#-------------------------------------------------------------------------------
#
# Disease_module<-subset.data.frame(all.genes,comp==1) # select only the main component
Disease_module<-all.genes # select all genes
#
#-------------------------------------------------------------------------------
# Over-representation analysis (ORA)  
#-------------------------------------------------------------------------------
#
top.results<-100 # max number of terms included in the final result
list.number<-1 # the number of gene lists required for a term to be included
ORA<-ORA.fun(Disease_module,top.results,list.number)
file.name<-"PF_output/ORA/ORA_Disease_Module.tsv"
write.table(ORA,file=file.name,quote=F,row.names=F,col.names=T,sep="\t")
#
#-------------------------------------------------------------------------------
# ORA visualization
#-------------------------------------------------------------------------------
#
top.n<-10 # how many top terms?
for (ORA.type in c("hsa","GO","R-HSA","DO")) {
  myORA<-subset.data.frame(ORA,grepl(ORA.type,ORA$ID))
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
    title.str<-"PF_output/ORA/KEGG enrichment"
  } else if (ORA.type=="GO") {
    title.str<-"PF_output/ORA/GO enrichment"
  } else if (ORA.type=="R-HSA") {
    title.str<-"PF_output/ORA/Reactome enrichment"
  } else {
    title.str<-"PF_output/ORA/DO enrichment"
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
pathway<-subset.data.frame(ORA,grepl("hsa",ORA$ID))$ID[1:10]
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
    file.rename(files[i],file.path("PF_output/ORA/KEGG",files[i]))   
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
  writeLines(entrez_ids,paste0("PF_output/ORA/Reactome/",pathway[i],".txt"))
}
#
# submit this file at https://reactome.org/PathwayBrowser/#TOOL=AT
#
#-------------------------------------------------------------------------------
# ORA for tissue-specific gene expression
#-------------------------------------------------------------------------------
#
Tissue.ORA(Disease_module)
#
#-------------------------------------------------------------------------------
# Hypergeometric tests for comparison with Proteomics 
#-------------------------------------------------------------------------------
#
# Prepare background universe
#
universe<-read_xlsx("Data/UKBB_proteomics.xlsx",sheet=1) # total effect 
NB<-nrow(universe)
#
# Prepare setA 
#
setA<-merge(Disease_module,universe,by="name") # keep only genes in universe
N1<-nrow(setA)
#
# Prepare setB (proteins significantly associated with ME/CFS)
#
setB<-subset.data.frame(universe,BH_pval_combined_te<=0.05) 
N2<-nrow(setB)
#
# Prepare intersection
#
intersection<-merge(setA,setB,by="name")
NR<-nrow(intersection)
#
# Perform hypergeometric test
#
myhypertest<-data.frame(setA="ME Module",setB="Proteomics")
pvalue<-as.numeric(phyper(q=NR-1,m=N1,n=NB-N1,k=N2,lower.tail=F))
myhypertest$sizeA<-N1
myhypertest$sizeB<-N2
myhypertest$overlap<-NR
myhypertest$background<-NB
myhypertest$pvalue<-formatC(pvalue,format="e",digits=2)
myhypertest$genes<-paste0(intersection$name,collapse="/")
#
file.name<-"PF_output/replication.tsv"
write.table(myhypertest,file=file.name,quote=F,row.names=F,col.names=T,sep="\t")
#
#-------------------------------------------------------------------------------
# Create a probability density for candidate genes
#-------------------------------------------------------------------------------
#
my.proband<-sample(STRING.names$preferred_name,1000)
mygenes<-data.frame(name=my.proband)
#
# For each gene calculate a score from interaction with genes of Disease Module
#
mygenes<-ScoreME(mygenes)
#
# Read previous results and add if different genes were stored
#
file.name<-"Data/Interaction_distribution.tsv"
if (file.exists(file.name)) {
  mygenes.old<-fread(file.name)
  index<-which(!mygenes.old$name%in%mygenes$name)
  mygenes<-rbind(mygenes.old[index,],mygenes)
}
hist(mygenes$score,25)
#
# Save result
#
write.table(mygenes,file=file.name,quote=F,row.names=F,col.names=T,sep="\t")
#
#-------------------------------------------------------------------------------
# Prioritize candidate genes of proband
#-------------------------------------------------------------------------------
#
my.proband<-c("LBR","COL4A4","ITGA9","CNTN5","PRKAB1","GLIPR1","RECQL","TNS2",
              "CHD9","EPX","ALOX12","GALR3","TFDP2","CNTN6","PCGF6","TYRO3",
              "CDH13","PTPRM","HEPH")
mygenes<-data.frame(Gene=my.proband)
#
# Use STRING preferred name for genes
#
gene.wo.STRING<-0
index<-c()
for (i in 1:nrow(mygenes)) {
  new.name<-STRING.name(mygenes$Gene[i])
  if (!is.na(new.name)) {
    mygenes$name[i]<-new.name
  } else if (is.na(new.name)) {
    print(paste(mygenes$Gene[i],"is not present in STRING's database"))
    gene.wo.STRING<-gene.wo.STRING+1
    index<-c(index,i)
  }
}
if (length(index)>0) mygenes<-mygenes[-index,] # only genes in STRING
#
# For each gene calculate a score from interaction with genes of Disease Module
#
mygenes<-ScoreME(mygenes)
#
# Compute a p-value for the scores
#
mygenes$pvalue<-rep(NA,nrow(mygenes))
file.name<-"Data/Interaction_distribution.tsv"
Distribution<-fread(file.name)
N<-nrow(Distribution)
for (i in 1:nrow(mygenes)) {
  index<-which(Distribution$score>=mygenes$score[i])
  mygenes$pvalue[i]<-length(index)/N
}
#
# Compute adjusted p-value for the scores
#
mygenes$p.adjust<-p.adjust(mygenes$pvalue,method="BH")
#
# Add NCBI id
#
mygenes$NCBI.id<-rep(NA,nrow(mygenes))
for (i in 1:nrow(mygenes)) {
  mygenes$NCBI.id[i]<-Symbol2NCBI.db(mygenes$name[i])
}
#
file.name<-"Proband.tsv"
write.table(mygenes,file=file.name,quote=F,row.names=F,col.names=T,sep="\t")

