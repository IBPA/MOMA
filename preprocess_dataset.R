getExponentialData <- function(data) {
   meta    <- read.table("Dataset/Meta.txt",header=T,sep="\t")
   indices <- NULL
   for (i in 1:nrow(data)) {
      gphase <- as.character(meta[match(data[i,"ID"],meta$ID),"Growth.Phase"])
      if (!is.na(gphase) && (gphase == "exponential" || gphase == "mid-exponential" || gphase == "exponential (predicted)"))
         indices <- append(indices,i)
   }
   return(data[indices,])
}

averageSameConditionProfiles <- function(data) {
   conditions <- unique(data[,"Cond"])
   data.avg   <- NULL
   for (i in 1:length(conditions)) {
      data.avg <- rbind(data.avg,colMeans(data[data[,"Cond"]==conditions[i], grep("m.b",colnames(data))]))
   }
   rownames(data.avg) <- conditions
   return(data.avg)
}

#### option {factor or nofactor} : if there are factors, then we min-max normalize only non-factor variables.
normalizeOmics <- function(data, option) {
   indices_factor       <- NULL
   data_factors         <- NULL
   data_without_factors <- NULL
   if (option == "factor") {
      for (i in 1:ncol(data)) {
         if (is.factor(data[,i]))
            indices_factor <- append(indices_factor,i)
      }
      data_factors         <- data[,indices_factor]
      data_without_factors <- data[,-c(indices_factor)]
   } else {
      data_without_factors <- data
   }
   data.norm           <- apply(data_without_factors,2,function(x) (x-min(x))/(max(x)-min(x)))
   colnames(data.norm) <- colnames(data_without_factors)
   data.norm           <- data.norm[,!is.nan(data.norm[1,])]
   if (option == "factor") {
      data.norm           <- cbind(data_factors, data.norm)
      colnames(data.norm) <- colnames(data)
   }
   return(data.norm)
}

expandMetadata <- function(d) {
    meta.medium                   <- read.table("Dataset/Meta.Medium.txt",header=T,sep="\t")
    meta.medium.nonbioinfo        <- c("ID","Base.Medium", "Description", "Link", "PMID", "Defined", "Description.1", "X")
    meta.strain                   <- read.table("Dataset/Meta.Strain.txt",header=T,sep="\t", quote="")
    meta.strain.nonbioinfo        <- c("Strain.Name", "Source", "PMID", "Alternate.Names", "Plasmid", "Sex", "Lambda", "Comments")

    d.with_meta_info <- addMetadata(d,                "medium", meta.medium, meta.medium.nonbioinfo, "ID",          2, 0)
    d.with_meta_info <- addMetadata(d.with_meta_info, "strain", meta.strain, meta.strain.nonbioinfo, "Strain.Name", 1, "no")
    d.with_meta_info <- addBinaryVariables(d.with_meta_info, "stress", 3, c("(Intercept)", "none"))
    d.with_meta_info <- addBinaryVariables(d.with_meta_info, "gp", 4, c("(Intercept)", "na_WT"))

    return(d.with_meta_info)
}

getBinaryVariables <- function(d, meta_idx, columns_to_remove) {
   factor_var <- NULL
   for (row_idx in 1:nrow(d)) {
      cond          <- as.character(rownames(d)[row_idx])
      cond.vector   <- unlist(strsplit(cond,"[.]"))
      factor_var    <- append(factor_var, cond.vector[meta_idx])
   }
   factor_var        <- factor(factor_var)
   binary_matrix     <- model.matrix(~factor_var)
   colnames(binary_matrix) <- gsub("factor_var","", colnames(binary_matrix))
   binary_matrix     <- binary_matrix[, -which(colnames(binary_matrix) %in% columns_to_remove)]
   return(binary_matrix)
}

addBinaryVariables <- function(d, meta_category, meta_idx, columns_to_remove) {
   binary_matrix            <- getBinaryVariables(d, meta_idx, columns_to_remove)
   multi_perturbation_names <- colnames(binary_matrix)[grep(";", colnames(binary_matrix))]
   for (multi_perturbation_name in multi_perturbation_names) {
      perturbations <- unlist(strsplit(multi_perturbation_name,"[;]"))
      for (perturbation in perturbations) {
         if (perturbation %in% colnames(binary_matrix))
            binary_matrix[,perturbation] <- binary_matrix[,perturbation] + binary_matrix[,multi_perturbation_name]
         else {
            binary_matrix <- cbind(binary_matrix, binary_matrix[,multi_perturbation_name])
            colnames(binary_matrix)[ncol(binary_matrix)] <- perturbation
         }
      }
      binary_matrix <- binary_matrix[,-which(colnames(binary_matrix) %in% multi_perturbation_name)]
   }
   colnames(binary_matrix) <- paste0(meta_category, ".", colnames(binary_matrix))
   new_d                   <- cbind(binary_matrix, d)
   rownames(new_d)         <- rownames(d)
   return(new_d)
}

addMetadata <- function(d, meta_category, meta, meta.nonbioinfo, id, meta_idx, empty_value) {
   d.with_meta_info     <- NULL
   for (row_idx in 1:nrow(d)) {
      cond                <- as.character(rownames(d)[row_idx])
      cond.vector         <- unlist(strsplit(cond,"[.]"))
      meta.bioinfo        <- as.numeric(meta[which(cond.vector[meta_idx]==meta[,id]),-which(colnames(meta) %in% meta.nonbioinfo)]!=empty_value)
      names(meta.bioinfo) <- paste0(meta_category,".", colnames(meta)[-which(colnames(meta) %in% meta.nonbioinfo)])
      d.with_meta_info    <- rbind(d.with_meta_info, c(meta.bioinfo, d[row_idx,]))
   }
   rownames(d.with_meta_info) <- rownames(d)
   return(d.with_meta_info)
}

args                           <- commandArgs(TRUE)
transcriptome_data             <- read.table(args[1],header=T,sep="\t")
transcriptome_data.exp         <- getExponentialData(transcriptome_data)               # extract exponential data only
transcriptome_data.avg         <- averageSameConditionProfiles(transcriptome_data.exp) # average expression levels for each condition
transcriptome_data.norm        <- normalizeOmics(transcriptome_data.avg,"nofactor")    # apply min-max normalization
transcriptome_data.with_meta   <- expandMetadata(transcriptome_data.norm)              # expand metadata features from "Cond" variable
transcriptome_data.with_meta   <- cbind(rownames(transcriptome_data.with_meta), transcriptome_data.with_meta)
colnames(transcriptome_data.with_meta)[1] <- "Cond"

write.table(transcriptome_data.with_meta, args[2], sep="\t", quote=F, row.names=F)
