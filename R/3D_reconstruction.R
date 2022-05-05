## Functions for 3D reconstruction of the developing cochlea

#' find_pca() function run PCA analysis on cell-by-gene matrix
#' @param data A cell-by-gene matrix with cells on the rows and genes on the columns. 
#' @param rotate A Boolean variable to indicate whether we want to rotate the PC axis for better visualization
#' @param pca_scale A logical value indicating whether the variables should be scaled to have unit 
#' variance before the analysis takes place. The default is FALSE 
find_pca = function(data,rotate = T, pca_scale = FALSE){
  # Run PCA
  pc = prcomp(data,center = T,scale. = pca_scale,rank. = 20)
  pc = as.data.frame(pc$x)[,1:6]
  
  if(rotate == TRUE){
    # Flip the order of PCs
    pc = pc[,c(2,1,3,4,5,6)]
    pc[,1] = -pc[,1]
    colnames(pc) = c("PC1","PC2","PC3","PC4","PC5","PC6")
  }
  else{
    colnames(pc) = c("PC1","PC2","PC3","PC4","PC5","PC6")
  }
  
  return(pc)
}

#' find_angle() function will re-arrange the data points on y axis for the scaled 2D PCA plot. 
#' @param x A number which you want to rotate. For example, x could be PC2 from PCA results
#' @param y A number which you want to rotate. For example, x could be PC3 from PCA results
find_angle = function(x,y){
  ratio = atan(y/x)
  
  if(x<0){
    ratio =  ratio + pi
  }
  if(ratio < 1.5){
    ratio = ratio + 2*pi
  }
  
  return(ratio)
}

#' get_2d() function will  re-arrange the data points based on PCs for scaled 2D PCA plot.
#' @param position position is a 2-column data frame containing x and y for find_angle() function. 
get_2d = function(position){
  result = vector()
  for(i in 1:nrow(position)){
    result[i] = find_angle(position[i,1],position[i,2])
  }
  
  return(result)
}

#' plot_2D() function projects PCA results with rotation on a scaled PCA plot to show the cell type distribution
#' @param pca_df A data frame from PCA results with cells on rows and PCs on columns
#' @param meta_df A meta information which includes the celltrails states, meta-cluster states, Seurat cluster, etc.
#' @param color Column name of meta_df indicates the color code in the pairwise PCA plots. It includes `celltrails_state`,
#' `meta_cluster`, `seurat_cluster`, and `seurat_cluster_anno`.
#' @param color_list Manually chosen colors to represent different clusters. By default, it is `color_list_meta_cluster`. 
#' @param gene The name of gene wants to project on the plot. By default, gene is NULL, which plot meta information.
#' @param seurat_object A Seurat object which includes the gene expression data for projection
#' @param ratio Aspect ratio of the plot. By default, `ratio` is 1.5. 
#' @param legend_position Legend position inherited from ggplot package. By default, `legend_position` is on the right side 
#' of the plot.
#' @param apex_flip A Boolean variable to indicate whether you want to flip the apex-to-base axis for better visualization
#' @param medial_flip A Boolean variable to indicate whether you want to flip the medial-lateral axis for better visualization
#' @return The function return 1 element which is a ggplot object. 
plot_2D = function(pca_df,meta_df,color,color_list = color_list_meta_cluster,gene = NULL,seurat_object,ratio = 1.5,legend_position="right",apex_flip=T, medial_flip = F){
  library(viridis)
  library(ggplot2)
  
  # Rotate PC2 and PC3 to get new x axis
  temp_y = get_2d(pca_df[,2:3])
  
  # Visualization
  if(apex_flip == TRUE){
    visual_df = data.frame(V1 = temp_y, V2 = -pca_df[,1])
  }
  else{
    visual_df = data.frame(V1 = temp_y, V2 = pca_df[,1])
  }
  if(medial_flip == TRUE){
    visual_df$V1 = abs(visual_df$V1 - 3*pi)
  }
  
  visual_df = cbind(visual_df,meta_df)
  
  if(is.null(gene)){
    g1 = ggplot(visual_df,aes(x = V1,y = V2, color = visual_df[,color])) + geom_point(size = 3) +
      scale_color_manual(values = color_list) + 
      theme_classic() +
      theme(axis.title = element_text(size = 30),axis.text = element_text(size = 25),legend.title = element_text(size = 30),legend.text = element_text(size = 25),legend.position = legend_position) +
      theme(axis.line.x  = element_line(size = 0.8),axis.line.y  = element_line(size = 0.8),axis.ticks = element_line(size = 0.8)) + 
      labs(x = "",y = "",color = "Cell Types") + xlim(0,3*pi) + 
      theme(aspect.ratio = ratio)
  }
  else{
    expression_vec = seurat_object@assays$RNA@data[which(rownames(seurat_object@assays$RNA@data) %in% gene),]
    visual_df$gene = expression_vec
    visual_df$gene[which(visual_df$gene == 0)] = NA
    
    color_map = viridis(15, option = "C")
    
    g1 = ggplot(visual_df,aes(x = V1,y = V2, color = gene)) + geom_point(size = 3) +
      scale_color_gradientn(colors = color_map,na.value = "grey80") + 
      theme_classic() +
      theme(axis.title = element_text(size = 30), axis.text = element_text(size = 25),legend.title = element_text(size = 30),legend.text = element_text(size = 25),legend.position = legend_position) +
      theme(axis.line.x  = element_line(size = 0.8),axis.line.y  = element_line(size = 0.8),axis.ticks = element_line(size = 0.8)) + 
      labs(x = "",y = "",color = gene) + xlim(0,3*pi) + 
      theme(aspect.ratio = ratio)
  }
  
  return(g1)
}

#' scale_2D() function normalizes the datasets with unit scale on a circle
#' @param df A 2-column data frame containing PC2 and PC3 for normalization
#' @return This function return a data frame after normalization
scale_2D = function(df){
  length_norm = sqrt(df[,1]^2+df[,2]^2)
  df[,1] = df[,1]/length_norm
  df[,2] = df[,2]/length_norm
  
  return(df)
} 

#' plot_cylinder() function projects 2D scaled PCA plots on a cylinder surface to represent 3D structure of the organ of Corti
#' @param pca_df A data frame from PCA results with cells on rows and PCs on columns
#' @param meta_df A meta information which includes the celltrails states, meta-cluster states, Seurat cluster, etc.
#' @param color Column name of meta_df indicates the color code in the pairwise PCA plots. It includes `celltrails_state`,
#' `meta_cluster`, `seurat_cluster`, and `seurat_cluster_anno`.
#' @param color_list Manually chosen colors to represent different clusters. By default, it is `color_list_meta_cluster`. 
#' @param gene The name of gene wants to project on the plot. By default, gene is NULL, which plot meta information.
#' @param seurat_object A Seurat object which includes the gene expression data for projection
#' @param ratio A scaling to adjust the height of the cylinder to simulate the cochlear duct extension. By default, `ratio`
#' is 0.6. 
#' @param x_flip A Boolean variable to indicate whether you want to flip the x-axis for better visualization. By default, 
#' `x_flip` is FALSE. 
#' @param x_flip A Boolean variable to indicate whether you want to flip the y-axis for better visualization. By default, 
#' `y_flip` is FALSE. 
#' @param x_flip A Boolean variable to indicate whether you want to flip the z-axis for better visualization. By default, 
#' `z_flip` is FALSE. 
#' @return This function saves a 3D cylinder reconstruction plot in your default folder 
plot_cylinder = function(pca_df,meta_df,color,color_list = color_list_meta_cluster,gene=NULL,seurat_object,ratio = 0.6,x_flip = F, y_flip = F, z_flip = F){
  library(plotly)
  
  # Run scale_2D for normalization
  pca_df[,2:3] = scale_2D(pca_df[,2:3])
  
  # Visualization
  visual_df = data.frame(V1 = pca_df[,3],V2 = pca_df[,2],V3 = pca_df[,1])
  visual_df = cbind(visual_df,meta_df)
  visual_df$V3 = ratio * visual_df$V3
  
  # Adjust the axis for better visualization
  if(z_flip == TRUE){
    visual_df$V3 = -visual_df$V3
  }
  if(x_flip == TRUE){
    visual_df$V1 = -visual_df$V1
  }
  if(y_flip == TRUE){
    visual_df$V2 = -visual_df$V2
  }
  
  # Prepare the cylinder structure
  angels = seq(0,2 * pi, length.out = nrow(visual_df))
  h_grid = seq(min(visual_df[,3]),max(visual_df[,3]),length.out = 10)
  
  cylinder_list = list()
  for(i in 1:length(h_grid)){
    cylinder_list[[i]] = data.frame(x = sin(angels),y = cos(angels),z = h_grid[i])
  }
  cylinder_df = cylinder_list[[1]]
  for(i in 2:length(cylinder_list)){
    cylinder_df = rbind(cylinder_df,cylinder_list[[i]])
  }
  
  # Plot meta info or gene expression values
  if(is.null(gene)){
    p1 = plot_ly(width = 10 * 100, height = 8 * 100) %>% 
      add_trace(visual_df, x = visual_df$V1, y = visual_df$V2, z = visual_df$V3, color = visual_df[,color], colors = color_list, type="scatter3d", mode="markers",marker = list(size = 20)) %>%
      add_trace(data = cylinder_df, x = cylinder_df$x, y = cylinder_df$y, z = cylinder_df$z,mode = "lines",type = "scatter3d",line=list(color='#CCCCCC',width = 10),showlegend = FALSE) %>% 
      layout(scene = list(xaxis = list(title = '',zeroline = FALSE,showline = FALSE,showticklabels = FALSE,showgrid = FALSE), 
                          yaxis = list(title = '',zeroline = FALSE,showline = FALSE,showticklabels = FALSE,showgrid = FALSE),
                          zaxis = list(title = '',zeroline = FALSE,showline = FALSE,showticklabels = FALSE,showgrid = FALSE,range = c(min(visual_df[,3]), 10)))) %>% 
      layout(showlegend = FALSE) 
    orca(p1, "cylinder_project.svg")
  }
  else{
    expression_vec = seurat_object@assays$RNA@data[which(rownames(seurat_object@assays$RNA@data) %in% gene),]
    visual_df$gene = expression_vec
    visual_df$gene[which(visual_df$gene == 0)] = NA
    
    p1 = plot_ly() %>% 
      add_trace(visual_df, x = visual_df$V1, y = visual_df$V2, z = visual_df$V3, color = visual_df$gene, type="scatter3d", mode="markers",marker=list(colorscale="Viridis")) %>% 
      add_trace(data = cylinder_df, x = cylinder_df$x, y = cylinder_df$y, z = cylinder_df$z,mode = "lines",type = "scatter3d",line=list(color='#CCCCCC',width = 6),showlegend = FALSE) %>% 
      layout(scene = list(xaxis = list(title = ''), yaxis = list(title = ''),zaxis = list(title = '')))
    orca(p1, paste(gene,"cylinder_project.svg",sep = "_"))
  }
  
  return(p1)
}


#' plot_donut() function projects scaled PCA plots on a donut shaped plot to show the organ of Corti structure
#' @param pca_df A data frame from PCA results with cells on rows and PCs on columns
#' @param meta_df A meta information which includes the celltrails states, meta-cluster states, Seurat cluster, etc.
#' @param color Column name of meta_df indicates the color code in the pairwise PCA plots. It includes `celltrails_state`,
#' `meta_cluster`, `seurat_cluster`, and `seurat_cluster_anno`.
#' @param color_list Manually chosen colors to represent different clusters. By default, it is `color_list_meta_cluster`. 
#' @param gene The name of gene wants to project on the plot. By default, gene is NULL, which plot meta information.
#' @param seurat_object A Seurat object which includes the gene expression data for projection
#' @param top_flip A Boolean element to indicate whether top and bottom cells should be flipped. We want to setup roof cells
#' at the top of the donut plot. By default, `top_flip` is False.
#' @param apex_flip A Boolean element to indicate whether apex and base should be flipped. We want to setup apical cells at
#' the inner cicle and basal cells at the outer circle. By default, `apex_flip` is False.
#' @param left_flip A Boolean element to indicate whether medial and lateral cells should be flipped. We want to setup 
#' medial cells at left and lateral cells at right. By default, `left_flip` is False.
#' @param enrich The name of GO term wants to project on the plot. By default, enrich is NULL, which plot gene information.
#' @param legend_position Legend position inherited from ggplot package. By default, `legend_position` is on the right side 
#' of the plot.
#' @return The function returns 1 element as a plot, which is a donut plot with color code of either meta information, or 
#' gene expression/GO term enrichment scores on single cell resolution.
plot_donut = function(pca_df,meta_df,color,color_list = color_list_meta_cluster,gene=NULL,seurat_object,top_flip=F,apex_flip=F,left_flip=F,enrich=NULL,legend_position="right"){
  library(plotly)
  library(viridis)
  
  # Run scale_2D for normalization
  if(apex_flip == T){
    pca_df[,1] = -pca_df[,1]
  }
  pca_df[,2:3] = scale_2D(pca_df[,2:3])
  
  # Rescale the data and project into a donunt-shape plot
  axis_min = min(pca_df[,1]) * 1.1
  axis_max = max(pca_df[,1]) * 1.1
  
  # Flip top and bottom
  if(top_flip == F){
    h = (-pca_df[,1] + axis_min)/(axis_max-axis_min) + 1.5
  }
  else{
    h = (pca_df[,1] - axis_min)/(axis_max-axis_min) -1.5
  }
  
  # Visualization
  visual_df = data.frame(V1 = pca_df[,3]*h,V2 = pca_df[,2]*h)
  if(left_flip == T){
    visual_df$V1 = -visual_df$V1 
  }
  visual_df = cbind(visual_df,meta_df)
  
  # Prepare the donut structure
  angels = seq(0,2*pi,length.out = 360)
  donut_df1 = data.frame(V1 = 0.5 * sin(angels),V2 = 0.5 * cos(angels))
  donut_df2 = data.frame(V1 = 1.5 * sin(angels),V2 = 1.5 * cos(angels))
  donut_df = rbind(donut_df1,donut_df2)
  
  # Plot meta info or gene expression values
  if(is.null(gene)){
    g1 = ggplot() +
      geom_point(data = visual_df,aes(x = V1, y = V2, color = visual_df[,color]),size=5)+
      scale_color_manual(values = color_list) +
      annotate("path",x=0.5*cos(seq(0,2*pi,length.out=360)),y=0.5*sin(seq(0,2*pi,length.out=360)),size = 1.5,color="#CCCCCC") +
      annotate("path",x=1.5*cos(seq(0,2*pi,length.out=360)),y=1.5*sin(seq(0,2*pi,length.out=360)),size = 1.5,color="#CCCCCC") +
      theme_classic() +
      labs(x = "",y = "",color="Cell Type")+
      theme(axis.text = element_text(size = 25),legend.title = element_text(size = 30),legend.text = element_text(size = 25),legend.position = legend_position) +
      theme(axis.line.x  = element_line(size = 0.8),axis.line.y  = element_line(size = 0.8),axis.ticks = element_line(size = 0.8)) + 
      labs(x = "",y = "",color = "Cell Types") + 
      theme(aspect.ratio = 1) +
      theme(panel.grid.major.x = element_blank(),axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.line = element_blank())
  }
  else{
    if(is.null(enrich) == F){
      expression_vec = meta_df[,enrich]
    }
    else{
      expression_vec = seurat_object@assays$RNA@data[which(rownames(seurat_object@assays$RNA@data) %in% gene),]
    }
    
    visual_df$gene = expression_vec
    visual_df$gene[which(visual_df$gene == 0)] = NA
    mid = mean(visual_df$gene,na.rm=T)
    
    color_map = viridis(15, option = "C")
    
    g1 = ggplot() +
      geom_point(data = visual_df,aes(x = V1, y = V2, color = gene),size=5)+
      scale_color_gradientn(colors = color_map,na.value = "grey80") +
      annotate("path",x=0.5*cos(seq(0,2*pi,length.out=360)),y=0.5*sin(seq(0,2*pi,length.out=360)),size = 1.5,color="#CCCCCC") +
      annotate("path",x=1.5*cos(seq(0,2*pi,length.out=360)),y=1.5*sin(seq(0,2*pi,length.out=360)),size = 1.5,color="#CCCCCC") +
      theme_classic() +
      theme(axis.text = element_text(size = 25),legend.title = element_text(size = 30),legend.text = element_text(size = 25),legend.position = legend_position) +
      theme(axis.line.x  = element_line(size = 0.8),axis.line.y  = element_line(size = 0.8),axis.ticks = element_line(size = 0.8)) + 
      labs(x = "",y = "",color = gene) + 
      theme(aspect.ratio = 1) +
      theme(panel.grid.major.x = element_blank(),axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.line = element_blank())
  }
  
  return(g1)
}





