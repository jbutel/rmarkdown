# Pipeline for processing Affymetrix Microarray Data

---

**Date:** August 9th, 2022  
**Revision:** -  
**Document Number:** - 

**Submitted by:**  
Julia Butelet (GL4U Intern 2022)

**Approved by:**  

---

# Table of contents  

- [Pipeline for processing Affymetrix Microarray Data](#pipeline-for-processing-affymetrix-microarray-data)
- [Table of contents](#table-of-contents)
- [Software used](#software-used)
  - [1. Set Up Environment](#1-set-up-environment)
  - [2. Import Raw Data](#2-import-raw-data)
  - [3. Quality Analysis of Raw Data](#3-quality-analysis-of-raw-data)
  - [4. Normalization](#4-normalization)
  - [5. Quality Analysis Visualization for Normalized Data](#5-quality-analysis-visualization-for-normalized-data)
  - [6. Differential Expression (DE) Analysis](#6-differential-expression-de-analysis)
  - [7. Gene Annotations](#7-gene-annotations)
  - [8. Creating Final Table](#8-creating-final-table)


---

# Software used  

|Program|Version|Relevant Links|
|:------|:------:|:-------------|
|R|4.2.1|[https://www.r-project.org/](https://www.r-project.org/)|
|Bioconductor|3.15.0|[https://bioconductor.org](https://bioconductor.org)|
|oligo|1.58.|[https://www.bioconductor.org/packages/release/bioc/html/oligo.html](https://www.bioconductor.org/packages/release/bioc/html/oligo.html)|
|biomaRt|2.50.3|[https://bioconductor.org/packages/release/bioc/html/biomaRt.html](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)|
|dplyr|1.0.9|[https://dplyr.tidyverse.org/](https://dplyr.tidyverse.org/)|
|AnnotationDbi|1.56.2|[https://bioconductor.org/packages/release/bioc/html/AnnotationDbi.html](https://bioconductor.org/packages/release/bioc/html/AnnotationDbi.html)
|tibble|3.1.8|[https://tibble.tidyverse.org/](https://tibble.tidyverse.org/)|
|Biobase|2.54.0|[https://www.bioconductor.org/packages/release/bioc/html/Biobase.html](https://www.bioconductor.org/packages/release/bioc/html/Biobase.html)|
|limma|3.52.2|[https://www.bioconductor.org/packages/release/bioc/html/limma.html](https://www.bioconductor.org/packages/release/bioc/html/limma.html)


---

## 1. Set Up Environment

```R
## Install required R packages if not already installed ##

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("oligo")
BiocManager::install("limma")
BiocManager::install("biomaRt")
BiocManager::install("AnnotationDbi")
BiocManager::install("Biobase")

## Import libraries ##

library(oligo)
library(limma)
library(biomaRt)
library(dplyr)
library(tibble)
library(AnnotationDbi)
library(Biobase)

```

<br>

---

## 2. Import Raw Data

```R
    #Reading in runsheet with metadata

csv_file <- read.table(file = r"(C:\Users\linde\rmarkdown\affy_processing\data_runsheets\GLDS_6_runsheet.csv)", 
                       sep=",",
                       header = TRUE,
                       check.names = FALSE
                       )

    #Accessing arrray data files located in column name 'Path to Raw Data File'

array_data_files <- csv_file["Path to Raw Data File"]

    #Importing and read *.CEL files using oligo function 'read.celfiles()'

raw_data <- oligo::read.celfiles(array_data_files) 

    #'raw_data' type - ExpressionFeatureSet Class

#   Parameter Definitions:

#file – path of runsheet file to be accessed
#sep - setting delimiter, default is white space
#header - determines if first row in dataframe is column names
#check.names - if TRUE, variable names are checked and corrected for synaticially validity
```

<br>

---

## 3. Quality Analysis of Raw Data

```R
##Assessing quality of raw data through visualization of plots

# First plot: density plot
# hist() function plots the density estimates for each sample
# Purpose: Identify technical variation between samples
    
density_plot <- hist(raw_data, transfo=log2, which=c("pm", "mm", "bg", "both", "all"),
                    nsample=10000, target = "core", main = "Density of raw intensities for multiple arrays")

# Second Plot: Pseudoimage
# image() function produces a pseudoimage for each array
# Purpose: Identify spatial abnormalities

pseudoimage <- image(raw_data, transfo=log2)

# Third Plot: MA plot
# MAplot() function creates MA plots using a reference array
# Purpose: Comparison of intensity values from one sample to another

MA_plot <- MAplot(raw_data, ylim=c(-2, 4))

#   Parameter Definitions:

# transfo - used to scale the data
# which - defines specific probe types
# nsample - sample size used to produce plot
# target - specifies group of meta-probeset
# main - main title 
# ylim - scale for y axis


##Assessing quality through statistics

# IQRray

# The IQRray statistic is obtained by ranking all the probe intensities from
# a given array and by computing the average rank for each probe set. The
# interquartile range (IQR) of the probe sets average ranks serves then as
# quality score. 
	#data - ExpressionFeatureSet object obtained after reading in celfiles into R 
	
	#obtaining intensity values for perfect match (pm) probes

pm_data<-pm(raw_data)
	
	#ranking probe intensities for every array

rank_data<-apply(pm_data,2,rank)
	
	#obtaining names of probeset for every probe
    
probeNames<-probeNames(raw_data)
	
	#function computing IQR of mean probe ranks in probesets 

get_IQR<-function(rank_data,probeNames){round(IQR(sapply(split(rank_data,probeNames),mean)),digits=2)}
	
	#computing arIQR score

IQRray_score<-apply(rank_data,2,get_IQR,probeNames=probeNames)	
```

<br>

---

## 4. Normalization

```R
##Probe Level normalization

    #Converting raw data from ExpressionFeatureSet object to matrix 

raw_table <- exprs(raw_data)

    #Performing background correction

backgroundCorrection <- preprocessCore::rma.background.correct(raw_table)

    #Performing Normalization using Quantile Normalization

probelevel_table <- preprocessCore::normalize.quantiles(backgroundCorrection, keep.names=TRUE)

    #Writing probe level normalization to file

write.csv(probelevel_table, r"(C:\Users\linde\rmarkdown\notebook\probelevel_table.csv)")

##Probeset Level Normalization

    #RMA (Robust Multi-array Average) is most commonly used method for preprocessing normalization of 
    #Affymetrix microarray data. There are three steps in the RMA algorithm: Convolution Background 
    #Correction, Quantile Normalization, and Summarization using Median Polish.

    #Performing RMA

normRMA <- rma(raw_data)

    #Converting normalized data to a matrix for easier accessibility 

norm_probeset_table <- exprs(normRMA)  

    #Writing probeset level normalization to file

write.csv(norm_probeset_table, r"(C:\Users\linde\rmarkdown\notebook\norm_probeset_table.csv)")
```

<br>

---

## 5. Quality Analysis Visualization for Normalized Data

```R
##QA to ensure normalization worked properly
##Samples should appear more similar to each other after normalization
##Same types of plots as QA of raw data for easy comparison

    #Normalized MA plot

MAplot_norm <- MAplot(normRMA, ylim=c(-2, 4))

    #Normlized Density plot

densityPlot_norm <- hist(normRMA, main = "Density of raw intensities for multiple arrays")
```

<br>

---

## 6. Differential Expression (DE) Analysis

```R
##Getting Factor Values

    #Creating a list of Factor Values

fv <- data.frame(csv_file["Factor Value"])

all_factor_values <- fv[,1]

    #Obtaining the name of each factor value

factor_values <- unique(all_factor_values)

    #Making the names for both lists syntactically valid

all.factor.values <- make.names(all_factor_values)

factor.values <- make.names(factor_values)

##Creating Design Matrix

    #Adding factor values to phenoData

ph = raw_data@phenoData

ph@data[ ,2] = c(all.factor.values)

colnames(ph@data)[2]="source"

    #Grouping factor values

groups = ph@data$source

    #Transforming group names into statistical factors

f = factor(groups,levels = factor.values)

    #Creating matrix of values of the grouping variable

design = model.matrix(~ 0 + f) 

colnames(design) <- factor.values

##Forming Contrast Matrix

    #Fit design matrix to a linear model using lmFit()

fit = limma::lmFit(norm_probeset_table,design)

fit.groups <- colnames(fit$design) [which(fit$assign == 1)]

fit.index <- which(levels(levels) %in% fit.groups)

    #Creating comparisons for contrast matrix

fit.group.names <- gsub(" ", "_", sub(",", "_", unique(csv_file$'Factor Value')))

combos <- combn(fit.groups, 2)

combos.names <- combn(fit.group.names, 2)

contrasts <- c(paste(combos[1,], combos[2,], sep = "-", paste(combos[2,], combos[1,], sep="-")))

contrast.names <- c(paste(combos.names[1,], combos.names[2,], sep = "v"), paste(combos.names[2,], 
                    combos.names[1,], sep = "v"))

    #Creating Contrast Matrix

cont.matrix <- limma::makeContrasts(contrasts = contrasts, levels=design)

contrast.fit <- limma::contrasts.fit(fit, cont.matrix)

##Testing and Accessing Results (Two Factor Values)

    #Performing statistical t-test using eBayes()

contrast.fit.eb <- limma::eBayes(contrast.fit)

    #Extracting log fold change values

log_fold_change <- contrast.fit.eb$coefficients

    #Extracting p-values

p_value <- contrast.fit.eb$p.value

    #Changing names

log.fold.change <- c(log_fold_change) 

p.value <- c(p_value)

    #Adjusting the p-values using the Benjamini-Hochberg methods

adjusted_pvalue <- p.adjust(p_value,method="BH")

#Combining normalized expressions table with test results

dge_dataframe <- cbind(norm_probeset_table, log.fold.change, p.value, adjusted_pvalue)

## If dataset has three factor values:

    #GLDS-25 (difficult to automate with more than 2 factor values)s

    #Making a dataframe of log fold change and p-values for easier accessibility

log_fold_change_df <- data.frame(log_fold_change)
p_value_df <- data.frame(p_value)
lfcdf <- rename(log_fold_change_df, "Log Fold Change GC-SF"="Ground.Control..Space.Flight..Space.Flight..Ground.Control.",
    "Log Fold Change GC-VC" = "Ground.Control..Vivarium.Control..Vivarium.Control..Ground.Control.",
    "Log Fold Change SF-VC" = "Space.Flight..Vivarium.Control..Vivarium.Control..Space.Flight.")
 
pvdf <- rename(p_value_df, "P-value GC-SF"="Ground.Control..Space.Flight..Space.Flight..Ground.Control.",
    "P-value GC-VC" = "Ground.Control..Vivarium.Control..Vivarium.Control..Ground.Control.",
    "P-value SF-VC" = "Space.Flight..Vivarium.Control..Vivarium.Control..Space.Flight.")

    #Extracting each set of p-values 

pvalue_gcsf <- c(p_value_df[,1])
pvalue_gcvc <- c(p_value_df[,2])
pvalue_sfvc <- c(p_value_df[,3])

    #Adjusting each set of p-values

adjusted_pvalue_gcsf <- p.adjust(p_value_df[,1],method="BH")
adjusted_pvalue_gcvc <- p.adjust(p_value_df[,2],method="BH")
adjusted_pvalue_sfvc <- p.adjust(p_value_df[,3],method="BH")

    #Adding to the dataframe

de_dataframe <- cbind(norm_probeset_table, lfcdf, pvdf, adjusted_pvalue_gcsf,
                 adjusted_pvalue_gcvc, adjusted_pvalue_sfvc)

```

<br>

---

## 7. Gene Annotations

```R
##Using specific attributes from example dataset GLDS-6

    #Searching Ensembl 
    
listEnsembl()

    #Selecting mart object pointing to the database

ensembl <- useEnsembl(biomart = "genes")

    #Searching datasets in this database

searchDatasets(mart = ensembl, pattern = "Rat") 

    #Selecting dataset

ensembl <- useDataset(dataset = "rnorvegicus_gene_ensembl", mart = ensembl)

    #Finding the microarray probe ID attribute
    
searchAttributes(mart = ensembl, pattern = "affy_")

    #Generate mapped dataframe using getBM()

my_affy_ids <- rownames(norm_probeset_table)
gene_annotations <- getBM(
    attributes = c(
        "affy_rat230_2",
        "ensembl_gene_id",
        "uniprot_gn_symbol",
        "go_id"
        ), 
    filters = "affy_rat230_2", 
    values = my_affy_ids,  
    mart = ensembl)

#   Parameter Definitions:
#attributes - desired attributes extracted from specified database
#filters - filters to be used in the query
#values - query IDs
#mart - object connected to desird database and dataset

    #Converting multimappings from multiple rows into list fields

gene_annotations <- gene_annotations %>% 

        #grouping by probe ID

         group_by(affy_rat230_2) %>% 

         summarise( 

            #Convert uniprot_gn_symbol to a string of unique values in the group

            SYMBOL = toString(unique(uniprot_gn_symbol)), 

            #Convert Ensembl IDs to a string of unique IDs in the group

            ENSEMBL = toString(unique(ensembl_gene_id)), 

            #Convert go_ids to a string of unique IDs in the group
            
            GO_Ids = toString(unique(go_id))
            )
```

<br>

---

## 8. Creating Final Table

```R
##Joining tables to map to items in DE dataframe using left.join()

    #Creating a column with probe IDS in DE dataframe

final_dataframe <- rownames_to_column(de_dataframe, var="probesets") 
        
    #Joining using probe IDS as common factor

        %>% left_join(gene_annotations, by = c("probesets" = "affy_rat230_2"))

    #Writing final dataframe to file

write.csv(final_dataframe, r"(C:\Users\linde\rmarkdown\notebook\GLDS_6_DGE.csv)")
```

<br>

---

**Pipeline Input data:**

- Array data files (*.CEL)
- Metadata runsheet containing sample names, factor values, and paths to array data files (*.csv)

**Pipeline Output data:**

- Normalized probe level expression values (*.csv)
- Normalized probeset level expression values (*.csv)
- Final table with statistical results and gene annotations (*.csv)