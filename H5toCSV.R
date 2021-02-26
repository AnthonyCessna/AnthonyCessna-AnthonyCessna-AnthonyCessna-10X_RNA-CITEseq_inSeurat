
if(!require(Seurat)){
  install.packages("Seurat")
  library(Seurat)
}
if(!require(plotscale)){
  install.packages("plotscale")
  library(plotscale)
}
if(!require(tidyverse)){
  install.packages("tidyverse")
  library(tidyverse)
}

if(!require(cowplot)){
  install.packages("cowplot")
  library(cowplot)
}

if(!require(SeuratData)){
  devtools::install_github('satijalab/seurat-data')
  library(SeuratData)
}
 
if(!require(gridExtra)){
  install.packages("gridExtra")
  library(gridExtra)
}
if(!require(patchwork)){
  install.packages("patchwork")
  library(patchwork)
}
# change to directory files are in, or delete line and set directory through RStudio
setwd("~/Desktop/Seurat/Seurat")

######################################## PARAMETERS ###########################################

Generate_UMAPS <- TRUE
          WNN_UMAP <- TRUE
          RNA_Protein_UMAP <- FALSE
         


###############################################################################################

#lists h5 files in current directory
H5_files <- list.files(pattern = ".h5")

#reads in 10X h5 files
for(i in H5_files){

 assign(i, Read10X_h5(i)) 
  
}
# removes file extension from the name of the object
H5_names <- gsub(".h5", "", H5_files)
H5_names_count <- 1
for(i in H5_files){
  
  assign(H5_names[H5_names_count], get(i))
  H5_names_count <- H5_names_count + 1
  
}

# Assigns Gene expression data to variable for every h5 file
Gene_Expression_names <- character()
Gene_Expression_names_count <- 1
for(i in H5_names){
  
  H5 <- get(i)
  name <- paste0("Gene_Expression_", i)
  assign(name, H5$`Gene Expression`)
  Gene_Expression_names[Gene_Expression_names_count] <- name
  Gene_Expression_names_count <- Gene_Expression_names_count + 1
  
}

#gets rid of TCR genes
for(i in Gene_Expression_names){
  
  Gene_Expression <- get(i)
  vdj_names <- grepl("^TR[ABGD][VDJC]", rownames(Gene_Expression))
  Gene_Expression <- Gene_Expression[!vdj_names, ]
  assign(i, Gene_Expression)
}

#Adds -Gene suffix to gene names
for(i in Gene_Expression_names){
  
  Gene_Expression <- get(i)
  Gene_names <- rownames(Gene_Expression)
  Gene_names <- paste0(Gene_names, "-Gene")
  rownames(Gene_Expression) <- Gene_names
  assign(i, Gene_Expression)
}

# Assigns Antibody capture data to variable for every h5 file
Antibody_capture_names <- character()
Antibody_capture_names_count <- 1
for(i in H5_names){
  
  H5 <- get(i)
  name <- paste0("Antibody_capture_", i)
  assign(name, H5$`Antibody Capture`)
  Antibody_capture_names[Antibody_capture_names_count] <- name
  Antibody_capture_names_count <- Antibody_capture_names_count + 1
  
}

# Adds -Protein suffix to every Antibody capture parameter
for(i in Antibody_capture_names){
  
  Antibody_capture <- get(i)
  Protein_names <- rownames(Antibody_capture)
  Protein_names <- paste0(Protein_names, "-Protein")
  rownames(Antibody_capture) <- Protein_names
  assign(i, Antibody_capture)
}

# Creates Seurat object for every h5 file
Gene_Expression_names_count <- 1
for(i in H5_names){
  
  Seurat_Object <- CreateSeuratObject(counts = get(Gene_Expression_names[Gene_Expression_names_count]),
                                      project = i, min.cells = 10)
  assign(i, Seurat_Object)
  Gene_Expression_names_count <- Gene_Expression_names_count + 1
  
}

# creates percent mitochondrial DNA feature set for all Seurat Objects
for(i in H5_names){
  
  Seurat_Object <- get(i)
  Seurat_Object[["percent.mt"]] <- PercentageFeatureSet(Seurat_Object, pattern = "^MT-")
  assign(i, Seurat_Object)
}

# adds Assay object of Antibody capture to all Seurat Objects
Antibody_capture_names_count <- 1
for(i in H5_names){
  
  Seurat_Object <- get(i)
  Seurat_Object[["CITE"]] <- CreateAssayObject(counts = get(Antibody_capture_names[Antibody_capture_names_count]))
  assign(i, Seurat_Object)
  Antibody_capture_names_count <- Antibody_capture_names_count + 1
}

# nFeature_RNA is number of genes detected in each cell
# nCount_RNA is the total number of molecules detected in each cell
ViolinPlot_Names <- character()
ViolinPlot_Names_Count <- 1
for(i in H5_names){
  name <- paste0("ViolinPlot_", i)
  Violin <- VlnPlot(get(i), features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  assign(name, Violin)
  ViolinPlot_Names[ViolinPlot_Names_Count] <- name
  ViolinPlot_Names_Count <- ViolinPlot_Names_Count + 1
  
}

# save Violin plots to PDF in directory
for(i in ViolinPlot_Names){
  
  name <- paste0(i, ".pdf")
  ggsave(name, plot = get(i), path = "Plots")
  
}

# subsets data to filter out mitochondrial contamination and cells with really high or really low reads
for(i in H5_names){
  
  Seurat_Object <- subset(get(i), subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
  assign(i, Seurat_Object)
}


#Normalize data for Gene and Protein data
for(i in H5_names){
  
  Seurat_Object <- NormalizeData(get(i), assay = "CITE", normalization.method = "CLR", margin = 2)
  Seurat_Object <- NormalizeData(get(i))
  assign(i, Seurat_Object)
}

# Pick 500 most variable genes
for(i in H5_names){
  
  Seurat_Object <- FindVariableFeatures(get(i), selection.method = "vst", nfeatures = 504)
  assign(i, Seurat_Object)
  
}

# set names of variable features for antibody capture
for(i in H5_names){
  Seurat_Object <- get(i)
  VariableFeatures <- rownames(Seurat_Object[["CITE"]])
  VariableFeatures(Seurat_Object, assay = "CITE") <- VariableFeatures
  assign(i, Seurat_Object)

}

# Scale gene data and protein data
for(i in H5_names){
  
  Seurat_Object <- get(i)
  Seurat_Object <- ScaleData(Seurat_Object)
  Seurat_Object <- ScaleData(Seurat_Object, assay = "CITE")
  assign(i, Seurat_Object)
  
}
if(Generate_UMAPS == TRUE){
# Run PCA for Genes
for(i in H5_names){
  
  Seurat_Object <- get(i)
  Seurat_Object <- RunPCA(Seurat_Object, approx=FALSE)
  assign(i, Seurat_Object)
}




# Run PCA for Antibody capture
for(i in H5_names){
  
  Seurat_Object <- get(i)
  Seurat_Object <- RunPCA(Seurat_Object, reduction.name = 'apca', assay = "CITE", approx=FALSE)
  assign(i, Seurat_Object)
}


#create Elbow plot to pick number of dimensions of PCA for gene data
for(i in H5_names){
  
  Seurat_Object <- get(i)
  name <- paste0(i, "_Elbow.pdf")
  Plot <- ElbowPlot(Seurat_Object, reduction = "pca")
  ggsave(name, plot = Plot, path = "Plots")
  
  }
}

if(Generate_UMAPS == TRUE){
# Calculate closest neighbors
for(i in H5_names){
  
  Seurat_Object <- get(i)
  Seurat_Object <-  FindMultiModalNeighbors(
    Seurat_Object, reduction.list = list("pca", "apca"), 
    dims.list = list(1:15, 1:15), modality.weight.name = "RNA.weight",
    weighted.nn.name = "weighted.nn"
  )
  
  assign(i, Seurat_Object)
}

if(WNN_UMAP == TRUE){
# Run UMAP combining Gene and Protein data
for(i in H5_names){
  
  Seurat_Object <- get(i)
  Seurat_Object <- RunUMAP(Seurat_Object, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  Seurat_Object <- FindClusters(Seurat_Object, graph.name = "wsnn", algorithm = 3, resolution = 1, verbose = FALSE)
  assign(i, Seurat_Object)
}

#Plots WNN UMAP
for(i in H5_names){
  
  Seurat_Object <- get(i)
  p1 <- DimPlot(Seurat_Object, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
  name <- paste0(i, "_wnnUMAP.pdf")
  ggsave(name, plot = p1, path = "Plots")
  }
}

# Run UMAP with protein data and gene data individually
if(RNA_Protein_UMAP == TRUE){
for(i in H5_names){
  
  Seurat_Object <- get(i)
  Seurat_Object <- RunUMAP(Seurat_Object, reduction = 'pca', dims = 1:15, assay = 'RNA', 
                           reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
  
  Seurat_Object <- RunUMAP(Seurat_Object, reduction = 'apca', dims = 1:15, assay = 'CITE', 
                           reduction.name = 'cite.umap', reduction.key = 'citeUMAP_')
  assign(i, Seurat_Object)
  
}

# plots and prints individual UMAPS
for(i in H5_names){
  
  Seurat_Object <- get(i)
  p1 <- DimPlot(Seurat_Object, reduction = 'rna.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
  p2 <- DimPlot(Seurat_Object, reduction = 'cite.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
  name1 <- paste0(i, "_rnaUMAP.pdf")
  name2 <- paste0(i, "_citeUMAP.pdf")
  ggsave(name1, plot = p1, path = "Plots")
  ggsave(name2, plot = p2, path = "Plots")
  
  }
}
  
}

#extract vector of feature names for RNA assay
RNA_Feature_names <- character()
RNA_Feature_names_count <- 1
for(i in H5_names){
  
  Seurat_Object <- get(i)
  RNA_Feature_Vector <- Seurat_Object@assays[["RNA"]]@var.features
  name <- paste0("RNA_names_", i)
  assign(name, RNA_Feature_Vector)
  RNA_Feature_names[RNA_Feature_names_count] <- name
  RNA_Feature_names_count <- RNA_Feature_names_count + 1
  
}

#extract vector of feature names for Antibody capture assay
AntibodyCapture_Feature_names <- character()
AntibodyCapture_Feature_names_count <- 1
for(i in H5_names){
  
  Seurat_Object <- get(i)
  AntibodyCapture_Feature_vector <- Seurat_Object@assays[["CITE"]]@var.features
  name <- paste0("AntibodyCapture_names_", i)
  assign(name, AntibodyCapture_Feature_vector)
  AntibodyCapture_Feature_names[AntibodyCapture_Feature_names_count] <- name
  AntibodyCapture_Feature_names_count <- AntibodyCapture_Feature_names_count + 1
}

if(WNN_UMAP == TRUE){
# creates feature plots for Antibody capture data at specified width
AntibodyCapture_Feature_names_count <- 1
PlotWidth <- 3
AntibodyCapture_FeaturePlot_names_vector <- character()
for(i in H5_names){
  AntibodyCapture_FeaturePlots <- character()
  AntibodyCapture_FeaturePlot_begin <- 1
  AntibodyCapture_FeaturePlot_end <- PlotWidth
  PlotCount <- 1
  Seurat_Object <- get(i)
  DefaultAssay(Seurat_Object) <- 'CITE'
  AntibodyCapture_Feature_vector <- get(AntibodyCapture_Feature_names[AntibodyCapture_Feature_names_count])
  
  
  Vector_remainder <- length(AntibodyCapture_Feature_vector) %% PlotWidth
  
  if(is.nan(Vector_remainder/Vector_remainder)){
  while_end <- ((length(AntibodyCapture_Feature_vector) - Vector_remainder) / PlotWidth)
  }else{
    
    while_end <- ((length(AntibodyCapture_Feature_vector) - Vector_remainder) / PlotWidth) + 1
  }
  while_count <- 1
  while(while_count <= while_end){
  
    if(Vector_remainder == 0){
      
      if(AntibodyCapture_FeaturePlot_end <= length(AntibodyCapture_Feature_vector) - Vector_remainder){
        
        FeatureVector <- AntibodyCapture_Feature_vector[AntibodyCapture_FeaturePlot_begin:AntibodyCapture_FeaturePlot_end]
        Plot <- FeaturePlot(Seurat_Object, features = FeatureVector, reduction = 'wnn.umap', max.cutoff = 10, 
                            cols = c("lightgray","red"), ncol = PlotWidth, keep.scale = "all")
        Plotname <- paste0("AntibodyCapturePlot_", PlotCount, "_", i)
        assign(Plotname, Plot)
        AntibodyCapture_FeaturePlots[PlotCount] <- Plotname
        PlotCount <- PlotCount + 1
        AntibodyCapture_FeaturePlot_begin <- AntibodyCapture_FeaturePlot_begin + PlotWidth
        AntibodyCapture_FeaturePlot_end <- AntibodyCapture_FeaturePlot_end + PlotWidth 
        while_count <- while_count + 1
      
      
      }
    } else{
      
      if(AntibodyCapture_FeaturePlot_end <= length(AntibodyCapture_Feature_vector) - Vector_remainder){
        
        FeatureVector <- AntibodyCapture_Feature_vector[AntibodyCapture_FeaturePlot_begin:AntibodyCapture_FeaturePlot_end]
        Plot <- FeaturePlot(Seurat_Object, features = FeatureVector, reduction = 'wnn.umap', max.cutoff = 10, 
                            cols = c("lightgray","red"), ncol = PlotWidth, keep.scale = "all")
        Plotname <- paste0("AntibodyCapturePlot_", PlotCount, "_", i)
        assign(Plotname, Plot)
        AntibodyCapture_FeaturePlots[PlotCount] <- Plotname
        PlotCount <- PlotCount + 1
        AntibodyCapture_FeaturePlot_begin <- AntibodyCapture_FeaturePlot_begin + PlotWidth
        AntibodyCapture_FeaturePlot_end <- AntibodyCapture_FeaturePlot_end + PlotWidth 
        while_count <- while_count + 1
      }else{
        FeatureVector <- AntibodyCapture_Feature_vector[AntibodyCapture_FeaturePlot_begin : length(AntibodyCapture_Feature_vector)]
        Plot <- FeaturePlot(Seurat_Object, features = FeatureVector, reduction = 'wnn.umap', max.cutoff = 10, 
                            cols = c("lightgray","red"), ncol = PlotWidth, keep.scale = "all" )
        
        if(Vector_remainder == 1){
          Plot <- Plot + plot_spacer() + plot_spacer()
        }else if(Vector_remainder == 2){
          
          Plot <- Plot + plot_spacer()
        }
        Plotname <- paste0("AntibodyCapturePlot_", PlotCount, "_", i)
        assign(Plotname, Plot)
        AntibodyCapture_FeaturePlots[PlotCount] <- Plotname
        PlotCount <- PlotCount + 1
        
        while_count <- while_count + 1
      }
      
      
      
    }
  
  }
  AntibodyCaptureName <- paste0("AntibodyCapture_FeaturePlots_",i)
  assign(AntibodyCaptureName, AntibodyCapture_FeaturePlots)
  AntibodyCapture_FeaturePlot_names_vector[AntibodyCapture_Feature_names_count] <- AntibodyCaptureName
  AntibodyCapture_Feature_names_count <- AntibodyCapture_Feature_names_count + 1
  
}


# saves Antibody FeaturePlots as PDF's 
for(i in AntibodyCapture_FeaturePlot_names_vector){
  
  name <- i
  AntibodyCapture_Feature_names_count <- 1
  pdf(paste0("Plots/", name, ".pdf"), width = 11, height = 7)
  AntibodyCapture_FeaturePlots <- get(i)
  
  if(length(AntibodyCapture_FeaturePlots) %% 2 == 0){
    
    while(AntibodyCapture_Feature_names_count < length(AntibodyCapture_FeaturePlots)){
      
      tempPlot <- get(AntibodyCapture_FeaturePlots[AntibodyCapture_Feature_names_count])/
              get(AntibodyCapture_FeaturePlots[AntibodyCapture_Feature_names_count + 1])
      
      plot(tempPlot)
      AntibodyCapture_Feature_names_count <- AntibodyCapture_Feature_names_count + 2
      
    }
    
  }else{
    
    while(AntibodyCapture_Feature_names_count < length(AntibodyCapture_FeaturePlots)){
      
      tempPlot <- get(AntibodyCapture_FeaturePlots[AntibodyCapture_Feature_names_count])/
        get(AntibodyCapture_FeaturePlots[AntibodyCapture_Feature_names_count + 1])
      
      plot(tempPlot)
      AntibodyCapture_Feature_names_count <- AntibodyCapture_Feature_names_count + 2
    }
  
    DummyRow <- plot_spacer() + plot_spacer() + plot_spacer() 
    tempPlot <- get(AntibodyCapture_FeaturePlots[AntibodyCapture_Feature_names_count])/
                DummyRow
    plot(tempPlot)
  }
  dev.off()
  }
  
}

# creates feature plots for RNA data
RNA_Feature_names_count <- 1
RNA_FeaturePlot_names_vector <- character()
for(i in H5_names){
  
  RNA_FeaturePlots <- character()
  RNA_FeaturePlot_Begin <- 1
  RNA_FeaturePlot_end <- PlotWidth
  
  Seurat_Object <- get(i)
  DefaultAssay(Seurat_Object) <- 'RNA'
  RNA_Feature_Vector <- get(RNA_Feature_names[RNA_Feature_names_count])
  
  for(j in 1:((length(RNA_Feature_Vector)) / PlotWidth)){
    
    FeatureVector <- RNA_Feature_Vector[RNA_FeaturePlot_Begin:RNA_FeaturePlot_end]
    Plot <- FeaturePlot(Seurat_Object, features = FeatureVector, reduction = 'wnn.umap', max.cutoff = 10, 
                        cols = c("lightgray","red"), ncol = PlotWidth, keep.scale = "all")
    Plotname <- paste0("RNAPlot_", j, "_", i)
    assign(Plotname, Plot)
    RNA_FeaturePlots[j] <- Plotname
    RNA_FeaturePlot_Begin <- RNA_FeaturePlot_Begin + PlotWidth
    RNA_FeaturePlot_end <- RNA_FeaturePlot_end + PlotWidth
  }
  
  RNACaptureName <- paste0("RNA_FeaturePlots_", i)
  assign(RNACaptureName, RNA_FeaturePlots)
  RNA_FeaturePlot_names_vector[RNA_Feature_names_count] <- RNACaptureName
  RNA_Feature_names_count <- RNA_Feature_names_count + 1
}

# saves RNA feature plot as a pdf
for(i in RNA_FeaturePlot_names_vector){
  
  name <- i
  RNA_Feature_names_count <- 1
  pdf(paste0("Plots/", name, ".pdf"), width = 11, height = 7)
  RNA_FeaturePlots <- get(i)
  
  while(RNA_Feature_names_count < length(RNA_FeaturePlots)){
    
    tempPlot <- get(RNA_FeaturePlots[RNA_Feature_names_count])/
                get(RNA_FeaturePlots[RNA_Feature_names_count + 1])
    
    plot(tempPlot)
    RNA_Feature_names_count <- RNA_Feature_names_count + 2
    
  }
  dev.off()
}

# saves cleaned SeuratObjects as files
for(i in H5_names){
  
 Seurat_Object <- get(i)
 name <- paste0(i, ".rds")
 saveRDS(Seurat_Object, name)
}


# Extract gene counts
Gene_Expression_names_count <- 1
for(i in H5_names){
  
  Seurat_Object <- get(i)
  DefaultAssay(Seurat_Object) <- 'RNA'
  name <- Gene_Expression_names[Gene_Expression_names_count]
  CountsVariableGenes <- t(as.matrix(Seurat_Object@assays$RNA@scale.data[Seurat_Object@assays$RNA@var.features,]))
  Barcodes <- rownames(CountsVariableGenes)
  CountsVariableGenes <- as_tibble(CountsVariableGenes)
  CountsVariableGenes$Barcode <- Barcodes
  CountsVariableGenes <- CountsVariableGenes %>% relocate(Barcode, .before = 1)
  assign(name, CountsVariableGenes)
  Gene_Expression_names_count <- Gene_Expression_names_count + 1
}

# Extracts Antibody capture counts
Antibody_capture_names_count <- 1
for(i in H5_names){
  
  Seurat_Object <- get(i)
  name <- Antibody_capture_names[Antibody_capture_names_count]
  CountsAntibodyCapture <- t(as.matrix(Seurat_Object@assays$CITE@scale.data))
  assign(name, CountsAntibodyCapture)
  Antibody_capture_names_count <- Antibody_capture_names_count + 1
  
}

#cleans and makes combined cleaned protein and gene data table
Matrix_names <- character()
Matrix_names_count <- 1
Antibody_capture_names_count <- 1
H5_names_count <- 1
for(i in Gene_Expression_names){
  
  Total_table <- cbind2(get(i), get(Antibody_capture_names[Antibody_capture_names_count]))
  Total_table <- as_tibble(Total_table)
  #Total_table <- Total_table %>% mutate(Cell_ID = row_number()) %>% relocate(Cell_ID)
  
  name <- paste0("CountMatrix_", H5_names[H5_names_count])
  Matrix_names[Matrix_names_count] <- name
  assign(name, Total_table)
  
  Matrix_names_count <- Matrix_names_count + 1
  Antibody_capture_names_count <- Antibody_capture_names_count + 1
  H5_names_count <- H5_names_count + 1
}

# Creates DataTables of only Antibody capture data
for(i in Antibody_capture_names){
  
  AntibodyTable <- get(i)
  AntibodyTable <- as_tibble(AntibodyTable)
  assign(i, AntibodyTable)
  
}

# gets metadata tables from all Seurat Objects
Matrix_names_count <- 1
for(i in H5_names){
  
  Seurat_Object <- get(i)
  MetaData <- Seurat_Object@meta.data
  MetaData <- MetaData %>% rownames_to_column(var = "Barcode")
  CountMatrix <- get(Matrix_names[Matrix_names_count])
  CountMatrix <- CountMatrix %>% full_join(MetaData, by = "Barcode")
  assign(Matrix_names[Matrix_names_count], CountMatrix)
  
  Matrix_names_count <- Matrix_names_count + 1
}

# extracts UMAP coordinates and adds them to Adata table
Matrix_names_count <- 1

for(i in H5_names){
  
  Seurat_Object <- get(i)
  UMAP_coords <- Seurat_Object@reductions[["wnn.umap"]]@cell.embeddings
  UMAP_coords <- as.data.frame(UMAP_coords)
  UMAP_coords <- UMAP_coords %>% rownames_to_column(var = "Barcode")
  CountMatrix <- get(Matrix_names[Matrix_names_count])
  CountMatrix <- CountMatrix %>% full_join(UMAP_coords, by = "Barcode")
  assign(Matrix_names[Matrix_names_count], CountMatrix)
  Matrix_names_count <- Matrix_names_count + 1
}

# Write CSV of Protein data only
for(i in Antibody_capture_names){
  
   AntibodyTable <- get(i)
   name <- paste0(i, ".csv")
   write_csv(AntibodyTable, name)
  
}

# Write CSV of protein and gene counts
for(i in Matrix_names){
  
  CountsTable <- get(i)
  name <- paste0(i, ".csv")
  write_csv(CountsTable, name)
  
}
                                                             





