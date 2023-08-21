# Load required packages

library(TCGAbiolinks)

# Set the working directory where you want to save the data
setwd("/path/to/your/directory")

# Specify the cancer type and data type
cancer_type <- "SKCM"  # Skin Cutaneous Melanoma
data_type <- "Gene expression"

# Query TCGA data
query <- GDCquery(project = "TCGA-SKCM", 
                  data.category = "Transcriptome Profiling", 
                  data.type = data_type)

# Download the data
GDCdownload(query)

# Load the downloaded data
data_path <- query$output$filename
expression_data <- read.table(data_path, header = TRUE, row.names = 1, sep = "\t")

# Get the sample details
sample_details <- GDCprepare(query)

# Print a summary of the loaded data and sample details
cat("Expression data dimensions:", dim(expression_data), "\n")
cat("Sample details dimensions:", dim(sample_details), "\n")

# You can now work with the expression data and sample details as needed
