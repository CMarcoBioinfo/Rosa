#!/usr/bin/env Rscript

check_install_load_package <- function(package_name, source = "CRAN") {
    if (!requireNamespace(package_name, quietly = TRUE)) {
        if (source == "Bioconductor") {
            BiocManager::install(package_name)
        } else {
            install.packages(package_name)
        }
    }
    library(package_name, character.only = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
}

args <- commandArgs(trailingOnly = TRUE)

check_install_load_package("stringr")

extension <- if (any(args %in% "--text")) ".txt" else ".xlsx"
bool_xlsx <- if (any(args %in% "--text")) FALSE else TRUE
input_file <- NULL
lenght_file <- NULL
graphics <- FALSE
StatAnalysis <- "YES"
thr=1


for (i in seq_along(args)) {
    if (args[i] == "--input"){
        input_file <- args[i + 1]
    } else if (args[i] == "--outputSignificant") {
        output_file_significant <- args[i + 1]
    } else if (args[i] == "--outputUnique") {
        output_file_unique <- args[i + 1]
    } else if (args[i] == "--Graphics") {
        graphics <- TRUE
    } else if (args[i] == "--directory") {
        chem_dessin <- args[i + 1]
    } else if (args[i] == "--length") {
        length <- args[i + 1]
        length <- as.integer(length)
    }
}

write_file <- function(df, output, xlsx, append_bool = FALSE, colNames = TRUE, startRow=1, cols_to_highlight=NULL) {
    if(xlsx) {
        if (append_bool && file.exists(output)) {
            wb <- loadWorkbook(output)
            writeData(wb, sheet = "Results", df, startCol = 1, startRow = startRow, rowNames = FALSE, colNames = colNames)
            saveWorkbook(wb, file = output, overwrite = TRUE)        
        } else {
             write.xlsx(df, file = output, rowNames = FALSE, colNames = colNames, sheetName = "Results", encoding = "UTF-8", na = "", startRow = startRow)
        }

        if (!is.null(cols_to_highlight) && file.exists(output)) {
            wb <- loadWorkbook(output)
            for (col in cols_to_highlight) {
                addStyle(wb, "Results", style = createStyle(fgFill = "#FFFF00"), cols = col, rows = 1:(nrow(df)+1), gridExpand = TRUE)
            }
            saveWorkbook(wb, file = output, overwrite = TRUE)
        }
    } else {
        write.table(df, file = output, row.names=FALSE, col.names=colNames, dec=".",sep="\t",quote=FALSE , append=append_bool)
    }
}

read_file <- function(file, xlsx) {
    if(xlsx) {
        data_genes <- read_excel(file)
        data_genes <- as.data.frame(data_genes)
        return(data_genes)
    } else {
        data_genes <- read.csv(file, sep = "\t" , header = TRUE, check.names = FALSE)
        data_genes <- as.data.frame(data_genes)
        return(data_genes)
    }
}

data_genes <- NULL
if (bool_xlsx) {
    check_install_load_package("readxl")
    check_install_load_package("openxlsx")
}

data_genes <- read_file(input_file, bool_xlsx)
significant_genes <- data_genes[!is.na(data_genes$Significative) & data_genes$Significative != "No" & !is.na(data_genes$filterInterpretation) & data_genes$event_type != "NoData" & data_genes$event_type != "Physio", ]
write_file(significant_genes, output_file_significant, bool_xlsx)
unique_genes <- data_genes[!is.na(data_genes$filterInterpretation) & data_genes$event_type != "NoData" & data_genes$event_type != "Physio" & data_genes$filterInterpretation == "Unique junction", ]
write_file(unique_genes, output_file_unique, bool_xlsx)

length_input <- 8 + length
EchName <- names(data_genes)[c(9:length_input)] 
SampleOutput <- names(data_genes)[c((length_input+1):(length_input+length))] 
SampleInput <- EchName


message("######################")
message("#Filter significant genes...")
message("######################")

new_columns <- data.frame(nbSignificantSamples = NA, p_value = NA, SignificanceLevel = "")    
significant_genes <- cbind(significant_genes, new_columns)
significant_genes$nbSignificantSamples <- str_count(significant_genes$Significative,"p-value")

for(j in 1:length(EchName)){
    message(paste("   Filter significant genes for:",EchName[j]))
    output_EchName = paste(chem_dessin,EchName[j], "/", EchName[j],".significant_junctions",extension,sep="")
    P_EchName <- paste("P",EchName[j],sep="_") 
    cols_to_highlight <- which(names(significant_genes) %in% c(P_EchName, EchName[j],"nbSignificantSamples"))
    write_file(significant_genes[0,], output_EchName, bool_xlsx)
    k <- 2
    for(i in 1:nrow(significant_genes))
        if(str_detect(significant_genes$Significative[i], EchName[j])){
            pattern <- paste0(EchName[j], ", p-value = ([0-9.eE+-]+)")
            p_value_sample <- as.numeric(str_match(significant_genes$Significative[i], pattern)[, 2])
            SignificanceLevel <- ''
            SignificanceLevel[p_value_sample<0.05& p_value_sample>0.01] <- '*'
            SignificanceLevel[p_value_sample<0.01& p_value_sample>0.001] <- '**'
            SignificanceLevel[p_value_sample<0.001] <- '***'            
            significant_genes$p_value[i] <- p_value_sample
            significant_genes$SignificanceLevel <- SignificanceLevel
            write_file(significant_genes[i,], output_EchName, bool_xlsx, append_bool=TRUE, colNames=FALSE, startRow=k)
            k <- k + 1
    }
    output <- read_file(output_EchName,bool_xlsx)
    write_file(output, output_EchName, bool_xlsx, append_bool=TRUE, cols_to_highlight=cols_to_highlight)
}


message("######################")
message("#Filter unique genes...")
message("######################")

unique_genes <- read_file(output_file_unique, bool_xlsx)
for(j in 1:length(EchName)){
    message(paste("   Filter unique genes for:",EchName[j]))
    P_EchName <- paste("P",EchName[j],sep="_") 
    output_EchName = paste(chem_dessin,EchName[j], "/", EchName[j],".unique_junctions",extension,sep="")
    unique_filter <- unique_genes[unique_genes[[P_EchName]] != 0,]

    cols_to_highlight <- which(names(unique_filter) %in% c(P_EchName, EchName[j], "nbSamp"))
    write_file(unique_filter, output_EchName, bool_xlsx, append_bool=TRUE, cols_to_highlight=cols_to_highlight)
}


