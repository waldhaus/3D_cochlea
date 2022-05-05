## Functions to identify gradually expressed genes (GEGs) along the tonotopic axis of the cochlea

## Step: 1: Identify the differentially expressed genes and further determine the GEGs using KL divergence 

#' regression() function fits the gene expression and PC1 using LOWESS method
#' @param pca_vector A 1D vector represents the x in the regression. Here we used PC1 as x. 
#' @param expression_vector A 1D vector represents the y in the regression. Here we used normalized gene expression
#' vector as y. The length of y has to be the same as the length of x.
#' @return The function returns two elements where each element is a 1D vector with the same length as the input vector.
#' The first element is the sorted smoothed x vector. The second element is the sorted smoothed y vector correspnding to
#' gene expression vector.
regression = function(pca_vector,expression_vector){
  library(stats)
  # x is PCA coordinates, and y is gene expression vector
  lowess_results = lowess(x = pca_vector,y = expression_vector, f = 0.5)
  lowess_x = lowess_results$x
  lowess_y = lowess_results$y
  
  return(list(lowess_x,lowess_y))
}


#' find_de_genes() function identifies the differentially expressed genes, determines the smoothed regression
#' line for gene experssion vector using LOWESS regression, and calculates the KL divergence between theoretical line
#' and regression line to determine the GEGs.
#' @param data A 2D cell-by-gene gene expression matrix
#' @param pca_df A data frame with rows represents the cells corresponding to the rows in data. The columns are PCs 
#' from PCA analysis. The data frame has at least 1 PC. 
#' @param groups A meta data information indicates the two groups we want to identify the DE genes. Here, we compared
#' two different clusters from Seurat.
#' @param thres A threshold to define the minimized difference between max and min of regressed gene expression vector.
#' By default, `thres` is 1. 
#' @param kl_thres KL divergence threshold to identify genes satifying the linear model. The smaller `kl_thres`, the
#' more stringent the results will be. By default, `kl_thres` is 0.3.
#' @return The function returns 5 elements as a list. The first element is gradually expressed gene list where the genes
#' are highly expressed in apex. The second element is gradually expressed gene list where the genes are highly expressed 
#' in base. The third element is a data frame containing the statistics analysis from this function. Specifically, it 
#' includes p-values from t-tests, BH adjusted p-values, FDR adjusted p-values, Bonferroni adjusted p-values, KL divergence
#' values, and a direction flag (0 represents the gene is highly expression in apex, vice versa.). The fourth element is 
#' a data frame containing the LOWESS regressed PC results. The fifth element is a data frame containing the LOWESS 
#' regressed gene expression results. The fourth and fifth elements have the same format with rows represent genes and 
#' columns represent cells.
find_de_genes = function(data, pca_df, groups, thres = 1, kl_thres = 0.3){
  library(philentropy)
  # Differential analysis using t-tests
  data_group1 = data[which(groups$seurat_cluster == unique(groups$seurat_cluster)[1]),]
  data_group2 = data[which(groups$seurat_cluster == unique(groups$seurat_cluster)[2]),]
  
  p_list = vector()
  de_gene_list = vector()
  kl_div_list = vector()
  regression_pca = data.frame(matrix(data = NA,nrow = ncol(data),ncol = nrow(data)))
  regression_expression = data.frame(matrix(data = NA,nrow = ncol(data),ncol = nrow(data)))
  direction_index = vector()
  
  # Loop for each gene
  for(i in 1:ncol(data)){
    if(i %% 1000 == 0){
      print(i)
    }
    # t-test
    p_list[i] = t.test(data_group1[,i],data_group2[,i], paired = FALSE, var.equal = FALSE)$p.value
    
    # LOWESS regression
    x = regression(pca_vector = pca_df[,1],expression_vector = data[,i])[[1]]
    y = regression(pca_vector = pca_df[,1],expression_vector = data[,i])[[2]]
    
    regression_pca[i,] = x
    regression_expression[i,] = y
    
    # Entropy
    if(p_list[i] < 0.05 & (max(y) - min(y) > thres)){
      index_y = which.max(y)
      
      # Generate therotical line with positive slope or negative slope
      if(index_y > (length(y)/2)){
        yy = seq(min(y),max(y),length.out = length(y))
        direction_index[i] = 0
      }
      else{
        yy = seq(max(y),min(y),length.out = length(y))
        direction_index[i] = 1
      }
      
      # Calculate the entropy between expression vector and theoretical vector
      kl_matrix = rbind(yy,y)
      kl_div_list[i] = KL(kl_matrix,unit = "log",est.prob = "empirical")
    }
    else{
      kl_div_list[i] = NA
      direction_index[i] = NA
    }
  }
  
  # Prepare dataset for return
  final_df = data.frame(p_value = p_list)
  final_df$p_value_BH = p.adjust(final_df$p_value,method = "BH")
  final_df$p_value_FDR = p.adjust(final_df$p_value,method = "fdr")
  final_df$p_value_BF = p.adjust(final_df$p_value,method = "bonferroni")
  final_df$kl_div = kl_div_list
  final_df$gradient_direction = direction_index
  final_df$gradient_direction = factor(final_df$gradient_direction,levels = c(0,1),labels = c("Apex","Base"))
  rownames(final_df) = colnames(data)
  
  rownames(regression_pca) = colnames(data)
  colnames(regression_pca) = rownames(data)
  
  rownames(regression_expression) = colnames(data)
  colnames(regression_expression) = rownames(data)
  
  # Setup a KL threshold for DE genes selection
  de_gene_list_apex = rownames(final_df)[which(final_df$kl_div <= kl_thres & final_df$gradient_direction == "Apex")]
  de_gene_list_base = rownames(final_df)[which(final_df$kl_div <= kl_thres & final_df$gradient_direction == "Base")]
  
  return(list(de_gene_list_apex,de_gene_list_base,final_df,regression_pca,regression_expression))
}


## Step 2: Identify the overlapped GEGs between two time points

#' filter_de_genes() function is to identity overlapped gradually expressed genes from two datasets (e.g. E12 and E14)
#' @param R1 A data frame from group 1 containing smoothed regression expression values for each gene and each cell. 
#' The `R1` data frame with genes on the rows and cells on the columns.
#' @param R2 A data frame from group 2 containing smoothed regression expression values for each gene and each cell. 
#' The `R2` data frame with genes on the rows and cells on the columns.
#' @param E1 A 1D vector containing KL divergence values from group 1. The lenght of E1 should match the number of rows
#' in data fram `R1`.
#' @param E2 A 1D vector containing KL divergence values from group 2. The lenght of E1 should match the number of rows
#' in data fram `R2`.
#' @param threshold A threshold to filter the genes based on KL divergence. By default, `threshold` is 0.9
#' @param diff A cutoff to filter the genes whose expression difference between maximum value and minimum value is 
#' larger than the cutoff By default, `diff` is 1.
#' @return The function returns 1 element, which is the overlapped gene list where the genes are gradually expressed
#' along the tonotopic axis.
filter_de_genes = function(R1, R2, E1, E2, threshold = 0.9, diff = 1){
  intersected_de_genes = c()
  
  # For loop each gene
  for(i in 1:nrow(R1)){
    if(i %% 1000 == 0){
      print(i)
    }
    
    # Thresholds
    if(is.na(E1[i]) != TRUE & is.na(E2[i]) != TRUE){
      # KL divergence threshold
      if(E1[i] < threshold & E2[i] < threshold){
        # Regression threshold
        if((max(R1[i,]) - min(R1[i,])) > diff & (max(R2[i,]) - min(R2[i,])) > diff){
          intersected_de_genes = append(intersected_de_genes,rownames(R1)[i])
        }
      }
    }
  }
  # Return a DE gene list for both datasets
  return(intersected_de_genes)
}


#' identify_direction() function is to identify the direction of the gradients
#' @param E1_results A data frame from `find_de_genes()` function including the direction information with `Apex` indicating
#' the gene is highly expressed in apex, vice versa. 
#' @param E2_results A second data frame from `find_de_genes()` function including the direction information with `Apex` 
#' indicating the gene is highly expressed in apex, vice versa. 
#' @param gene_list A vector containing genes of interest from `filter_de_genes()` function. The gene list is where users want
#' to identify the directions. 
#' @return The function returns 2 element as a list, which includes a gene list with genes highly expressed in apex and a gene 
#' list with genes highly expressed in base. 
identify_direction = function(E1_results, E2_results, gene_list){
  # Identify the direction of highly expressed genes
  E1_apex = rownames(E1_results)[which(E1_results$gradient_direction == "Apex")]
  E1_base = rownames(E1_results)[which(E1_results$gradient_direction == "Base")]
  
  E2_apex = rownames(E2_results)[which(E2_results$gradient_direction == "Apex")]
  E2_base = rownames(E2_results)[which(E2_results$gradient_direction == "Base")]
  
  # Identify intersected genes
  intersected_apex_list = intersect(E1_apex,E2_apex)
  intersected_base_list = intersect(E1_base,E2_base)
  
  # Separate gene list into two directions
  gene_list_apex = gene_list[which(gene_list %in% intersected_apex_list)]
  gene_list_base = gene_list[which(gene_list %in% intersected_base_list)]
  
  return(list(gene_list_apex, gene_list_base))
}


