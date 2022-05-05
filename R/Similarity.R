## Functions to test the two hypothesis model about cochlea extension using similarity method

#' cosine_distance() function is to calculate cosine distance between two vectors
#' @param X1 Matrix1 with n obersavation on rows where columns represent features
#' @param X2 Matrix2 with n obersavation on rows where columns represent features
#' @param index Indicator values to locate the number of rows for both `X1` and `X2`. The `index`
#' is a vector containing two values.
#' @return The function returns 1 element, which is a similarity matrix calculated by 
#' 1-cosine_distance. The higher value, the more similarity between two cells. 
cosine_distance = function(X1, X2, index){
  if(index[1] == 1){
    print(index[2])
  }
  
  A = X1[index[1],]
  B = X2[index[2],]
  
  # cosine distance is 1 - similarity
  cos_dis = 1 - sum(A*B)/sqrt(sum(A^2)*sum(B^2))
  
  return(cos_dis)
}

#' calculate_cell_similarity() function is to calculate pairwise cell similarity
#' @param data1 Data frame 1 or matrix1 with n obersavation on rows where columns represent features
#' @param data2 Data frame 2 or matrix2 with n obersavation on rows where columns represent features
#' @param method The method to calculate pairwise cell similarity. There are a few options to choose
#' `cosine` represents cosine distance and `Euclidean` indicates Euclidean distance. 
#' @return The function returns 1 element, which is a matrix of pairwise cell similarity. 
calculate_cell_similarity = function(data1,data2,method = c("cosine","Euclidean")){
  library(pdist)
  
  c1 = nrow(data1)
  c2 = nrow(data2)
  
  # Calculate cosine distance
  if(method == "cosine"){
    cmb = expand.grid(i = 1:c1,j = 1:c2)
    result_matrix = matrix(apply(cmb,1, function(x) cosine_distance(data1, data2, x)), c1, c2)
  }
  # Calculate Eucledian distance
  if(method == "Euclidean"){
    euclidean_matrix = pdist::pdist(data1, data2)
    result_matrix = as.matrix(euclidean_matrix)
  }
  
  return(result_matrix)
}


#' distance_to_similarity() function is to convert
#' @param distance_matrix Matrix calculated from calculate_cell_similarity() function. Each element is the distance_
#' matrix is the distance between two datasets
#' @param method The method to convert distance matrix to similarity matrix. `inverse` and `gaussian_kernel` are for
#' Euclidean distance. `inverse` take the inverse of distance with 1 to avoid divided by 0. `gaussian_kernel` has 
#' an additional parameter `gaussian_sigma` to control the band width. `cosine` similarity if for cosine distance, 
#' which is basically 1-cosine_distance.
#' @param gaussian_sigma The band width of gaussian kernel. By default, `gaussian_sigma` is 1
#' @return The function returns 1 element as matrix format, which is the similarity matrix between two datasets
distance_to_similarity = function(distance_matrix, method = c("inverse","cosine","gaussian_kernel"),gaussian_sigma = 1){
  if(method == "cosine"){
    similarity_matrix = 1 - distance_matrix
  }
  else if(method == "inverse"){
    similarity_matrix = 1/(1 + distance_matrix)
  }
  else if(method == "gaussian_kernel"){
    similarity_matrix = exp(-distance_matrix^2/(2*gaussian_sigma^2))
  }
  
  return(similarity_matrix)
}

#' generate_meta_cell() is to generate meta cells to denoise the dataset
#' @param data A dataframe or a matrix with n obersavation on rows where columns represent features
#' @param order_index A vector of order index of cells to group the cells
#' @param chunk The number of chunk to divide for meta cells. By default, `chunk` is NULL. 
#' @param ncell The number of cells to generate meta cell. By default, `ncells` is 10.
#' @return The function returns 1 element, which is an average ordered matrix/data frame.
generate_meta_cell = function(data,order_index,chunk = NULL,ncells = 10){
  ordered_data = data[order_index,]
  nr = nrow(ordered_data)
  
  # Using the same number of cell for each chunk
  if(is.null(chunk) == TRUE){
    split_list = split(ordered_data, rep(1:ceiling(nr/ncells), each=ncells, length.out = nr))
    
    # For loop to calculate mean gene expression value for each meta cell
    ave_ordered_data = c()
    for(i in 1:length(split_list)){
      chunk_df = split_list[[i]]
      ave_chunk_df = base::colMeans(chunk_df)
      ave_ordered_data = rbind(ave_ordered_data,ave_chunk_df)
    }
  }
  # Use the a fixed chunk to generate meta cells
  else{
    ncells = nr %/% (chunk-1)
    split_list = split(ordered_data, rep(1:chunk, each = ncells,length.out = nr))
    
    # For loop to calculate mean gene expression value for each meta cell
    ave_ordered_data = c()
    for(i in 1:length(split_list)){
      chunk_df = split_list[[i]]
      ave_chunk_df = base::colMeans(chunk_df)
      ave_ordered_data = rbind(ave_ordered_data,ave_chunk_df)
    }
  }
  
  return(ave_ordered_data)
}