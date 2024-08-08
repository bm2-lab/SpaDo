
DrawCluster<-function (data, label = NULL, point_size = 1, method = c("tsne", 
    "umap"), draw_cluster_text = TRUE, calculated = TRUE, pca = TRUE, 
    perplexity = 100, plot = TRUE, legend_plot = TRUE, seed = 1) 
{
    require(ggplot2)
    require(dplyr)
    set.seed((seed))
    method = method[1]
    if (calculated) {
        if (method == "tsne") {
            require(Rtsne)
            tsneresult2 <- Rtsne(t(data), perplexity = perplexity, 
                pca = pca)
            X <- as.data.frame(tsneresult2$Y)
        }
        else if (method == "umap") {
            require(umap)
            umapresult1 <- umap(t(data))
            X <- as.data.frame(umapresult1$layout)
        }
        else {
            print("method must be tsne or umap.")
            break
        }
    }
    else {
        X = data
    }
    if (length(label) == 0) {
        label <- array(1, dim(X)[1])
        labelname = c(1)
    }
    label<-label[colnames(data)]
    labelname <- names(table(label))
    p <- ggplot(X, aes(x = X[, 1], y = X[, 2]))
    cell_group = factor(label)
    if (method == "tsne") {
        p <- p + geom_point(aes(color = cell_group), size = point_size) + 
            xlab("tSNE1") + ylab("tSNE2")
    }
    else {
        p <- p + geom_point(aes(color = cell_group), size = point_size) + 
            xlab("umap1") + ylab("umap2")
    }
    if (draw_cluster_text) {
        Label_cal <- X
        Label_cal$cluster <- label
        cluster_x_y <- Label_cal %>% group_by(cluster) %>% summarise(x_median = median(V1), 
            y_median = median(V2))
        p <- p + annotate("text", x = cluster_x_y$x_median, y = cluster_x_y$y_median, 
            label = cluster_x_y$cluster)
    }
    if (plot) {
        if (legend_plot) {
            mytheme <- theme_bw() + theme(plot.title = element_text(size = rel(1.5), 
                hjust = 0.5), axis.title = element_text(size = rel(1)), 
                axis.text = element_text(size = rel(1)), panel.grid.major = element_line(color = "white"), 
                panel.grid.minor = element_line(color = "white"), 
                legend.text = element_text(size = 10), legend.title = element_text(size = 15))
            p <- p + mytheme + guides(colour = guide_legend(override.aes = list(size = 4)))
            print(p)
        }
        else {
            mytheme <- theme_bw() + theme(plot.title = element_text(size = rel(1.5), 
                hjust = 0.5), axis.title = element_text(size = rel(1)), 
                axis.text = element_text(size = rel(1)), panel.grid.major = element_line(color = "white"), 
                panel.grid.minor = element_line(color = "white"), 
                legend.text = element_text(size = 10), legend.title = element_text(size = 15), 
                legend.position = "none")
            p <- p + mytheme + guides(colour = guide_legend(override.aes = list(size = 4)))
            print(p)
        }
    }
    return(list(p = p, x = X, cell_group = cell_group))
}




Correlation<-function(matrix,method=c("pearson","spearman","cosine","euclidean"),cpu_num=4){
  cosdist <- function(x1,x2){
    n1 <- sqrt(sum(x1^2))
    n2 <- sqrt(sum(x2^2))
    d <- as.numeric(x1 %*% x2) / n1 / n2
    d
  }
  euclidean<-function(x1,x2){
    return(sqrt(t(x1-x2) %*% (x1-x2)))
  }
  method<-method[1]
  if(!(method %in% c("pearson","spearman","cosine","euclidean"))){
    print("method must be 'pearson','spearman','cosine' or 'euclidean'.")
    break
  }
  require(parallel)
  cpu_num_set <- makeCluster(getOption("cluster.cores", cpu_num))
  sample_num<-ncol(matrix)
  cor_matrix<-matrix(rep(1,sample_num^2),sample_num)
  colnames(cor_matrix)<-colnames(matrix)
  row.names(cor_matrix)<-colnames(matrix)
  simi_cor<-function(vec,matrix,method){
    if(method=="cosine"){
      cor_matrix<-apply(matrix,2,cosdist,x2=vec)
    }else if(method=="euclidean"){
      cor_matrix<-apply(matrix,2,euclidean,x2=vec)
    }else{
      cor_matrix<-apply(matrix,2,cor,y=vec,method=method)
    }
    return(cor_matrix)
  }
  cor_matrix<-parApply(cpu_num_set,matrix,2,simi_cor,matrix=matrix,method=method)
  stopCluster(cpu_num_set)
  return(cor_matrix)
}

SpatialNormalize<-function (expression_profile, ST_method = c("osmFISH", "merFISH", 
    "STARmap", "Others")) 
{
    ST_method <- ST_method[1]
    if (ST_method %in% c("osmFISH", "merFISH", "STARmap")) {
        normalize_scale <- function(x) {
            return(x/sum(x))
        }
        expression_profile[is.na(expression_profile)] <- 0
        expression_profile <- expression_profile[apply(abs(expression_profile), 
            1, sum) > 0, ]
        expression_profile <- apply(expression_profile, 2, normalize_scale)
        expression_profile[is.na(expression_profile)] <- 0
        expression_profile[is.nan(expression_profile)] <- 0
        expression_profile <- log(expression_profile + 1)
    }
    else {
        normalize_scale_factor <- function(x) {
            return(10000 * (x/sum(x)))
        }
        expression_profile[is.na(expression_profile)] <- 0
        expression_profile <- expression_profile[apply(abs(expression_profile), 
            1, sum) > 0, ]
        expression_profile <- apply(expression_profile, 2, normalize_scale_factor)
        expression_profile[is.na(expression_profile)] <- 0
        expression_profile[is.nan(expression_profile)] <- 0
        expression_profile <- log(expression_profile + 1)
        return(expression_profile)
    }
}


InitialClustering<-function (expression_profile, user_offered = FALSE, sample_information_user_offered = NULL, 
    nfeatures = 2000, resolution = (ifelse(user_offered == FALSE, 
        2, FALSE))) 
{
    initial_clustering_result <- list()
    initial_clustering_result$expression_profile <- expression_profile
    require(Seurat)
    sce_seurat <- CreateSeuratObject(expression_profile)
    sce_seurat<-NormalizeData(sce_seurat,normalization.method = NULL,verbose=FALSE)
    sce_seurat <- FindVariableFeatures(sce_seurat, nfeatures = nfeatures,verbose=FALSE)
    sce_seurat <- ScaleData(sce_seurat,verbose=FALSE)
    sce_seurat <- RunPCA(sce_seurat, features = VariableFeatures(object = sce_seurat), 
        ndims.print = 1, nfeatures.print = 1,verbose=FALSE)
    if (user_offered) {
        initial_clustering_result$sample_information <- sample_information_user_offered
    }
    else {
        sce_seurat <- FindNeighbors(sce_seurat, reduction = "pca")
        sce_seurat <- FindClusters(sce_seurat, resolution = resolution)
        sample_cluster <- Idents(sce_seurat)
        names(sample_cluster) <- colnames(expression_profile)
        initial_clustering_result$sample_information <- sample_cluster
    }
    initial_clustering_result$high_var_genes <- VariableFeatures(object = sce_seurat)
    initial_clustering_result$sce_seurat <- sce_seurat
    return(initial_clustering_result)
}

SpatialPlot<-function (initial_clustering_result, sample_information_coordinate) 
{
    sce_seurat <- initial_clustering_result$sce_seurat
    Idents(sce_seurat) <- initial_clustering_result$sample_information
    require(ggplot2)
    sce_seurat <- RunUMAP(sce_seurat, dims = 1:15,verbose=FALSE)
    sample_information_coordinate$cluster <- Idents(sce_seurat)
    cluster_plot <- DimPlot(sce_seurat, reduction = "umap")
    coordinate_plot <- ggplot(sample_information_coordinate, 
        aes(x = X, y = Y, colour = cluster)) + geom_point(size = 1.2)
    plot(cluster_plot)
    plot(coordinate_plot)
}

SpatialCellTypeDistribution<-SpatialCellTypeDistribution<-function (sample_information_coordinate, sequence_resolution = c("single_cell", 
    "spot"), sample_information_cellType = NULL, sample_information_decon = NULL, 
    neighbour_search_method = (ifelse(sequence_resolution == 
        "single_cell", "KNN", "radius")), k = (ifelse(neighbour_search_method == 
        "KNN", 30, FALSE)), r = (ifelse(neighbour_search_method == "radius", 2, FALSE))) 
{
    SpatialKNN <- function(sample_information_coordinate) {
        require(Seurat)
        test_coordinate <- sample_information_coordinate
        test_coordinate_expand <- test_coordinate
        for (i in 1:5) {
            test_coordinate_expand <- cbind(test_coordinate, 
                test_coordinate_expand)
        }
        colnames(test_coordinate_expand) <- paste(colnames(test_coordinate_expand), 
            1:ncol(test_coordinate_expand))
        sce_seurat <- CreateSeuratObject(t(test_coordinate_expand))
        sce_seurat<-NormalizeData(sce_seurat,normalization.method = NULL,verbose=FALSE)
        sce_seurat <- ScaleData(sce_seurat,verbose=FALSE)
        sce_seurat <- RunPCA(sce_seurat, features = rownames(sce_seurat), 
            ndims.print = 1, nfeatures.print = 1, npcs = 5,verbose=FALSE)
        sce_seurat@reductions$pca@cell.embeddings <- as.matrix(test_coordinate)
        sce_seurat <- FindNeighbors(sce_seurat, reduction = "pca", 
            dims = 1:2, return.neighbor = T, k.param = 50,verbose=FALSE)
        knn_sample <- sce_seurat@neighbors$RNA.nn@nn.idx
        knn_value <- sce_seurat@neighbors$RNA.nn@nn.dist
        cell_names <- sce_seurat@neighbors$RNA.nn@cell.names
        row.names(knn_sample) <- cell_names
        row.names(knn_value) <- cell_names
        knn_sample <- t(apply(knn_sample, 1, function(x, y) {
            y <- cell_names
            return(y[x])
        }))
        return(list(knn_value = knn_value, knn_sample = knn_sample))
    }
    SpatialKNN_result <- SpatialKNN(sample_information_coordinate)
    require(reshape2)
    neighbour_search_method = neighbour_search_method[1]
    knn_sample <- SpatialKNN_result$knn_sample
    knn_value <- SpatialKNN_result$knn_value
    knn_choose <- list()
    if (neighbour_search_method == "KNN") {
        if (k < 2) {
            print("k must be bigger than 1")
            break
        }
        for (i in 1:nrow(knn_sample)) {
            knn_choose[[i]] <- knn_sample[i, 1:k]
        }
    }
    else if (neighbour_search_method == "radius") {
        for (i in 1:nrow(knn_sample)) {
            knn_choose[[i]] <- knn_sample[i, knn_value[i, ] <= 
                r]
        }
    }
    else {
        print("Parameter 'neighbour_search_method' must be 'KNN' or 'radius'!")
        break
    }
    names(knn_choose) <- row.names(knn_sample)
    knn_choose_cellType <- knn_choose
    sequence_resolution <- sequence_resolution[1]
    if (sequence_resolution == "single_cell") {
        names(knn_choose_cellType) <- sample_information_cellType[names(knn_choose)]
        for (i in 1:length(knn_choose_cellType)) {
            knn_choose_cellType[[i]] <- sample_information_cellType[knn_choose[[i]]]
        }
        cellType_names <- names(table(sample_information_cellType))
        knn_matrix_list <- list()
        knn_matrix_cellType_list <- list()
        for (i in 1:length(cellType_names)) {
            knn_matrix_cellType <- knn_choose_cellType[names(knn_choose_cellType) == 
                cellType_names[i]]
            knn_matrix <- knn_choose[names(knn_choose_cellType) == 
                cellType_names[i]]
            knn_matrix_list[[i]] <- knn_matrix
            knn_matrix_cellType_list[[i]] <- knn_matrix_cellType
        }
        names(knn_matrix_list) <- cellType_names
        names(knn_matrix_cellType_list) <- cellType_names
        knn_matrix <- list(knn_matrix_cellType_list = knn_matrix_cellType_list, 
            knn_matrix_list = knn_matrix_list)
        cellType_names <- names(table(sample_information_cellType))
        cell_type_distribution_each_list <- list()
        cell_type_distribution <- matrix(NA, 0, length(cellType_names))
        colnames(cell_type_distribution) <- cellType_names
        cell_type_distribution_rownames <- c()
        for (k in 1:length(cellType_names)) {
            test_cellType_matrix <- knn_matrix[[1]][[k]]
            test_matrix <- knn_matrix[[2]][[k]]
            cell_type_distribution_rownames <- c(cell_type_distribution_rownames, 
                names(test_matrix))
            cell_type_distribution_each <- matrix(0, length(test_cellType_matrix), 
                length(cellType_names))
            colnames(cell_type_distribution_each) <- cellType_names
            for (i in 1:length(test_cellType_matrix)) {
                cellType_table <- table(factor(test_cellType_matrix[[i]], 
                  levels = cellType_names))
                if (sum(cellType_table) == 0) {
                  cell_type_distribution_each[i, ] <- cellType_table
                }
                else {
                  cell_type_distribution_each[i, ] <- cellType_table/sum(cellType_table)
                }
            }
            row.names(cell_type_distribution_each) <- names(test_matrix)
            cell_type_distribution_each_list[[k]] <- cell_type_distribution_each
            cell_type_distribution <- rbind(cell_type_distribution, 
                cell_type_distribution_each)
        }
        names(cell_type_distribution_each_list) <- cellType_names
        row.names(cell_type_distribution) <- cell_type_distribution_rownames
        return(cell_type_distribution)
    }
    else if (sequence_resolution == "spot") {
        sample_information_cellType <- rep(c("a", "b", "c","d","e"), 
            length.out = nrow(sample_information_decon))
        names(sample_information_cellType) <- row.names(sample_information_decon)
        names(knn_choose_cellType) <- sample_information_cellType[names(knn_choose)]
        for (i in 1:length(knn_choose_cellType)) {
            knn_choose_cellType[[i]] <- sample_information_cellType[knn_choose[[i]]]
        }
        cellType_names <- names(table(sample_information_cellType))
        knn_matrix_list <- list()
        knn_matrix_cellType_list <- list()
        for (i in 1:length(cellType_names)) {
            knn_matrix_cellType <- knn_choose_cellType[names(knn_choose_cellType) == 
                cellType_names[i]]
            knn_matrix <- knn_choose[names(knn_choose_cellType) == 
                cellType_names[i]]
            knn_matrix_list[[i]] <- knn_matrix
            knn_matrix_cellType_list[[i]] <- knn_matrix_cellType
        }
        names(knn_matrix_list) <- cellType_names
        names(knn_matrix_cellType_list) <- cellType_names
        knn_matrix <- list(knn_matrix_cellType_list = knn_matrix_cellType_list, 
            knn_matrix_list = knn_matrix_list)
        cellType_names <- colnames(sample_information_decon)
        cell_type_distribution <- matrix(NA, 0, length(cellType_names))
        colnames(cell_type_distribution) <- cellType_names
        cell_type_distribution_rownames <- c()
        for (k in 1:length(knn_matrix[[2]])) {
            test_matrix <- knn_matrix[[2]][[k]]
            cell_type_distribution_rownames <- c(cell_type_distribution_rownames, 
                names(test_matrix))
            for (i in 1:length(test_matrix)) {
                if (length(test_matrix[[i]]) == 1) {
                  cell_type_distribution <- rbind(cell_type_distribution, 
                    sample_information_decon[test_matrix[[i]], 
                      ])
                }
                else {
                  cell_type_distribution <- rbind(cell_type_distribution, 
                    apply(sample_information_decon[test_matrix[[i]], 
                      ], 2, mean))
                }
            }
        }
        row.names(cell_type_distribution) <- cell_type_distribution_rownames
        return(cell_type_distribution)
    }
    else {
        print("sequence_resolution must be one of 'single_cell' and 'spot'!")
        break
    }
}

DistributionDistance <- function(cell_type_distribution,distance = c("JSD", "manhattan"),no_cores=1) {
    method_choose <- distance[1]
    data_matrix<-t(cell_type_distribution)
    if (method_choose == "manhattan") {
        propor_dis <- dist(x = cell_type_distribution, method = "manhattan")
        propor_dis <- as.matrix(propor_dis)
        return(propor_dis)
    }else if(method_choose=="JSD"){
        if(no_cores<1){
            print("n_core cannot be less than 1!")
            break
        }
        if(no_cores==1){
            require(philentropy)
            propor_dis <- philentropy::distance(cell_type_distribution, method = method_choose)
            row.names(propor_dis) <- row.names(cell_type_distribution)
            colnames(propor_dis) <- row.names(cell_type_distribution)
            return(propor_dis)   
        }else{
            require(parallel)
            max_cores <- detectCores()
            if(no_cores>max_cores){
                print("no_core cannot bigger than the max cores in your computer!")
                break
            }else{
                require(lsa)
                kl_divergence <- function(p, q) {
                    p <- p + .Machine$double.eps
                    q <- q + .Machine$double.eps
                    return(sum(p * log(p / q)))
                }
                jsd <- function(p, q) {
                    m <- (p + q) / 2
                    return((kl_divergence(p, m) + kl_divergence(q, m)) / 2)
                }
                calculate_jsd_matrix <- function(data_matrix,no_cores) {
                    n_samples <- ncol(data_matrix)
                    jsd_matrix <- matrix(0, n_samples, n_samples)
                    cl <- makeCluster(no_cores)
                    clusterExport(cl, list("data_matrix", "jsd", "kl_divergence"))
                    jsd_matrix <- parLapply(cl, 1:n_samples, function(i) {
                        sapply(1:n_samples, function(j) {
                            if (i == j) {
                                return(0)
                            } else {
                                    return(jsd(data_matrix[,i], data_matrix[,j]))
                            }
                        })
                    })  
                  stopCluster(cl)
                  jsd_matrix <- do.call(cbind, jsd_matrix)
                  return(jsd_matrix)
                }
                distances <- calculate_jsd_matrix(data_matrix, no_cores)
                return(distances)   
            }           
        }
    }
}

DomainHclust<-function (distribution_distance, autoselection = TRUE, auto_resolution = c(0,1,2,3,4), domain_num = 10) 
{
    ## fastcluster::hclust has faster speed than stats::hclust. 
    ## Besides, fastcluster::hclust is without the limitation that size cannot be NA nor exceed 65536
    library(fastcluster)
    result_hc <- fastcluster::hclust(d = as.dist(distribution_distance), method = "ward.D2")
    if (autoselection) {
        require(dynamicTreeCut)
        auto_resolution <- auto_resolution[1]
        cluster_hc_dy <- cutreeDynamic(result_hc, distM = distribution_distance, deepSplit = auto_resolution)
        cluster_hc_dy_num <- length(table(cluster_hc_dy))
    }else {
        cluster_hc_dy_num = domain_num
    }
    cluster_hc_all = cutree(result_hc, k = 1:cluster_hc_dy_num)
    cluster_hc_all_order <- cluster_hc_all[result_hc$order, ]
    cluster_hc_all_order_rename <- as.data.frame(cluster_hc_all_order)
    colnames(cluster_hc_all_order_rename) <- paste("Domain_level", 
        1:cluster_hc_dy_num, sep = "_")
    cluster_hc_all_order_rename_c <- apply(cluster_hc_all_order_rename, 
        2, function(x) {
            return(as.character(x))
        })
    row.names(cluster_hc_all_order_rename_c) <- row.names(cluster_hc_all_order_rename)
    cluster_hc_all_order_rename_c <- as.data.frame(cluster_hc_all_order_rename_c)
    for (i in 1:ncol(cluster_hc_all_order_rename_c)) {
        cluster_hc_all_order_rename_c[, i] <- paste("Domain", i, cluster_hc_all_order_rename_c[, i], sep = "_")
    }
    return(list(hclust_result_df = cluster_hc_all_order_rename_c, 
        hclust_result_model = result_hc, hclust_result_order = cluster_hc_all_order))
}

DomainPlot<-function (domain_hclust, distribution_distance, sample_information_coordinate = NULL, 
    k = ncol(domain_hclust[[3]]), size = 1, shape = 19, aspect_ratio = 1) 
{
    require(ggplot2)
    if (k > ncol(domain_hclust[[3]])) {
        print(paste("The lagerest domain number is : ", ncol(domain_hclust[[3]]), 
            "k cannot be larger!", sep = ","))
        break
    }
    cluster_hc_all_order <- domain_hclust[[3]]
    cluster_hc = domain_hclust[[1]][, k]
    result_hc <- domain_hclust[[2]]
    names(cluster_hc) <- row.names(domain_hclust[[1]])
    a <- DrawCluster(t(distribution_distance), method = "umap", 
        label = cluster_hc[row.names(distribution_distance)])
    colour_use <- ggplot_build(a$p)$data[[1]]$colour
    colour_use_unique <- unique(colour_use)
    names(colour_use_unique) <- 1:length(colour_use_unique)
    colour_use_unique_ok <- colour_use_unique[unique(cluster_hc_all_order[, 
        k])]
    source("http://addictedtor.free.fr/packages/A2R/lastVersion/R/code.R")
    b <- A2Rplot(result_hc, k = k, boxes = FALSE, col.up = "gray50", 
        col.down = colour_use_unique_ok, show.labels = F, legend = T, 
        only.tree = T, type = c("triangle"))
    sample_information_coordinate$cluster_hc <- domain_hclust[[1]][row.names(sample_information_coordinate), 
        k]
    c <- ggplot(sample_information_coordinate, aes(x = X, y = Y, 
        colour = factor(cluster_hc))) + geom_point(size = size, 
        shape = shape) + theme(axis.text.x = element_text(vjust = 0.5, 
        hjust = 0.5, angle = 45), aspect.ratio = aspect_ratio)
    print(c)
}

DomainCharacterization<-function (domain_hclust, cell_type_distribution, k = ncol(domain_hclust[[1]]), 
    single_cell = TRUE, sample_information_cellType = NULL, plot = TRUE) 
{
    cluster_hc = domain_hclust[[1]][, k]
    names(cluster_hc) <- row.names(domain_hclust[[1]])
    if (!single_cell) {
        selected_samples <- names(cluster_hc[cluster_hc == names(table(cluster_hc))[1]])
        selected_propor <- cell_type_distribution[selected_samples, 
            ]
        propor_mean <- apply(selected_propor, 2, mean)
        df <- data.frame(propor_mean)
        for (i in 2:length(table(cluster_hc))) {
            selected_samples <- names(cluster_hc[cluster_hc == 
                names(table(cluster_hc))[i]])
            selected_propor <- cell_type_distribution[selected_samples, 
                ]
            propor_mean <- apply(selected_propor, 2, mean)
            df[, i] <- propor_mean
        }
        colnames(df) <- names(table(cluster_hc))
        df_plot <- melt(df)
        df_plot$cellType <- rep(row.names(df), ncol(df))
        colnames(df_plot) <- c("Domain", "Proportion", "Cell_type")
        plot_c <- ggplot(df_plot, aes(Domain, Proportion, fill = Cell_type)) + 
            geom_bar(stat = "identity", color = "white", size = 0.2) + 
            guides(fill = guide_legend(reverse = F)) + scale_y_continuous(expand = c(0, 
            0)) + theme(axis.text.x = element_text(vjust = 0.5, 
            hjust = 0.5, angle = 45), aspect.ratio = 1)
        if (plot) {
            plot(plot_c)
        }
        return(list(domain_average_distribution = df_plot))
    }
    else {
        sample_information_test <- sample_information_cellType
        cell_type_distribution_df <- as.data.frame(cell_type_distribution)
        cell_type_distribution_df$domain <- cluster_hc[row.names(cell_type_distribution_df)]
        selected_samples <- names(sample_information_test[sample_information_test == 
            names(table(sample_information_test))[1]])
        selected_propor <- table(cluster_hc[selected_samples])/length(cluster_hc[selected_samples])
        df_cellType <- as.data.frame(matrix(0, length(table(sample_information_test)), 
            length(table(cluster_hc))))
        row.names(df_cellType) <- names(table(sample_information_test))
        colnames(df_cellType) <- names(table(cluster_hc))
        df_cellType[1, names(selected_propor)] <- selected_propor
        for (i in 2:length(table(sample_information_test))) {
            selected_samples <- names(sample_information_test[sample_information_test == 
                names(table(sample_information_test))[i]])
            selected_propor <- table(cluster_hc[selected_samples])/length(cluster_hc[selected_samples])
            df_cellType[i, names(selected_propor)] <- selected_propor
        }
        occurence_cellType <- Correlation(t(df_cellType))
        selected_samples <- names(cluster_hc[cluster_hc == names(table(cluster_hc))[1]])
        selected_propor <- cell_type_distribution[selected_samples, 
            ]
        propor_mean <- apply(selected_propor, 2, mean)
        df <- data.frame(propor_mean)
        for (i in 2:length(table(cluster_hc))) {
            selected_samples <- names(cluster_hc[cluster_hc == 
                names(table(cluster_hc))[i]])
            selected_propor <- cell_type_distribution[selected_samples, 
                ]
            propor_mean <- apply(selected_propor, 2, mean)
            df[, i] <- propor_mean
        }
        colnames(df) <- names(table(cluster_hc))
        df_plot <- melt(df)
        df_plot$cellType <- rep(row.names(df), ncol(df))
        colnames(df_plot) <- c("Domain", "Proportion", "Cell_type")
        df_cellType_plot <- melt(t(df_cellType))
        colnames(df_cellType_plot) <- c("Domain", "Cell_type", 
            "Proportion")
        if (plot) {
            library(pheatmap)
            library(ggplot2)
            plot_b <- ggplot(df_cellType_plot, aes(Cell_type, 
                Proportion, fill = Domain)) + geom_bar(stat = "identity", 
                color = "white", size = 0.2) + guides(fill = guide_legend(reverse = F)) + 
                scale_y_continuous(expand = c(0, 0)) + theme(axis.text.x = element_text(vjust = 0.5, 
                hjust = 0.5, angle = 45), aspect.ratio = 1)
            plot_c <- ggplot(df_plot, aes(Domain, Proportion, 
                fill = Cell_type)) + geom_bar(stat = "identity", 
                color = "white", size = 0.2) + guides(fill = guide_legend(reverse = F)) + 
                scale_y_continuous(expand = c(0, 0)) + theme(axis.text.x = element_text(vjust = 0.5, 
                hjust = 0.5, angle = 45), aspect.ratio = 1)
            pheatmap(occurence_cellType)
            plot(plot_b)
            plot(plot_c)
        }
        return(list(occurence_rate = occurence_cellType, cell_type_average_distribution = df_cellType_plot, 
            domain_average_distribution = df_plot))
    }
}


SpatialCellTypeDistribution_multiple<-function (sample_information_coordinate_list, sequence_resolution = c("single_cell", 
    "spot"), sample_information_cellType_list = NULL, sample_information_decon_list = NULL,r=2,k=30) 
{
    sequence_resolution <- sequence_resolution[1]
    cell_type_distribution_list <- list()
    if (sequence_resolution == "single_cell") {
        if (is.null(sample_information_cellType_list)) {
            print("'sample_information_cellType_list' shouldn't be NULL!")
        }
        else {
            for (i in 1:length(sample_information_coordinate_list)) {
                cell_type_distribution_list[[i]] <- SpatialCellTypeDistribution(sample_information_coordinate = sample_information_coordinate_list[[i]], 
                  sequence_resolution = "single_cell", sample_information_cellType = sample_information_cellType_list[[i]],k=k)
            }
        }
    }
    else if (sequence_resolution == "spot") {
        if (is.null(sample_information_decon_list)) {
            print("'sample_information_decon_list' shouldn't be NULL!")
        }
        else {
            for (i in 1:length(sample_information_coordinate_list)) {
                cell_type_distribution_list[[i]] <- SpatialCellTypeDistribution(sample_information_coordinate = sample_information_coordinate_list[[i]], 
                  sequence_resolution = "spot", sample_information_decon = sample_information_decon_list[[i]],r=r)
            }
        }
    }
    names(cell_type_distribution_list) <- names(sample_information_coordinate_list)
    cellType_union <- colnames(cell_type_distribution_list[[1]])
    for (i in 2:length(cell_type_distribution_list)) {
        cellType_union <- union(cellType_union, colnames(cell_type_distribution_list[[i]]))
    }
    for (j in 1:length(cell_type_distribution_list)) {
        if (length(cellType_union) - length(colnames(cell_type_distribution_list[[j]])) > 
            0) {
            missing_data_1 <- matrix(0, nrow(cell_type_distribution_list[[j]]), 
                length(cellType_union) - length(colnames(cell_type_distribution_list[[j]])))
            colnames(missing_data_1) <- setdiff(cellType_union, 
                colnames(cell_type_distribution_list[[j]]))
            cell_type_distribution_list[[j]] <- cbind(cell_type_distribution_list[[j]], 
                missing_data_1)
            cell_type_distribution_list[[j]] <- cell_type_distribution_list[[j]][, 
                cellType_union]
        }
        row.names(cell_type_distribution_list[[j]]) <- paste(names(cell_type_distribution_list)[j], 
            row.names(cell_type_distribution_list[[j]]), sep = "_")
    }
    cell_type_distribution_combine <- cell_type_distribution_list[[1]]
    datasets_name <- names(cell_type_distribution_list)
    datasets_lable <- c(rep(datasets_name[1], nrow(cell_type_distribution_list[[1]])))
    names(datasets_lable) <- row.names(cell_type_distribution_list[[1]])
    for (k in 2:length(cell_type_distribution_list)) {
        cell_type_distribution_combine <- rbind(cell_type_distribution_combine, 
            cell_type_distribution_list[[k]])
        datasets_lable2 <- c(rep(datasets_name[k], nrow(cell_type_distribution_list[[k]])))
        names(datasets_lable2) <- row.names(cell_type_distribution_list[[k]])
        datasets_lable <- c(datasets_lable, datasets_lable2)
    }
    return(list(cell_type_distribution_combine = cell_type_distribution_combine, 
        datasets_lable = factor(datasets_lable)))
}

DomainPlot_multiple<-function (domain_hclust, distribution_distance, k = ncol(domain_hclust[[3]]), 
    size = 1, shape = 19, aspect_ratio = 1) 
{
    require(ggplot2)
    if (k > ncol(domain_hclust[[3]])) {
        print(paste("The lagerest domain number is : ", ncol(domain_hclust[[3]]), 
            "k cannot be larger!", sep = ""))
        break
    }
    cluster_hc_all_order <- domain_hclust[[3]]
    cluster_hc = domain_hclust[[1]][, k]
    result_hc <- domain_hclust[[2]]
    names(cluster_hc) <- row.names(domain_hclust[[1]])
    a <- DrawCluster(t(distribution_distance), method = "umap", 
        label = cluster_hc[row.names(distribution_distance)])
    colour_use <- ggplot_build(a$p)$data[[1]]$colour
    colour_use_unique <- unique(colour_use)
    names(colour_use_unique) <- 1:length(colour_use_unique)
    colour_use_unique_ok <- colour_use_unique[unique(cluster_hc_all_order[, 
        k])]
    source("http://addictedtor.free.fr/packages/A2R/lastVersion/R/code.R")
    b <- A2Rplot(result_hc, k = k, boxes = FALSE, col.up = "gray50", 
        col.down = colour_use_unique_ok, show.labels = F, legend = T, 
        only.tree = T, type = c("triangle"))
}

DomainDistribution<-function (hclust_result_df, datasets_lable_combine, domain_num = ncol(hclust_result_df)) 
{
    domain_result_each_list <- list()
    SO_tree <- hclust_result_df[, domain_num]
    names(SO_tree) <- row.names(hclust_result_df)
    domain_names <- names(table(SO_tree))
    datasets_cellNum <- table(datasets_lable_combine)
    datasets_num <- length(datasets_cellNum)
    datasets_name <- names(datasets_cellNum)
    proportion_dataset_domain <- matrix(0, domain_num, datasets_num)
    Term_table <- table(SO_tree)
    row.names(proportion_dataset_domain) <- names(Term_table)
    colnames(proportion_dataset_domain) <- names(datasets_cellNum)
    for (j in 1:length(Term_table)) {
        Datasets_term_cellNum <- table(datasets_lable_combine[names(SO_tree[SO_tree == 
            names(Term_table)[j]])])
        proportion_dataset_domain[j, ] <- Datasets_term_cellNum/datasets_cellNum
    }
    for (i in 1:datasets_num) {
        domain_result_each_list[[i]] <- hclust_result_df[names(datasets_lable_combine[datasets_lable_combine == 
            datasets_name[i]]), ]
    }
    names(domain_result_each_list) <- datasets_name
    return(list(domain_distribution = proportion_dataset_domain, 
        domain_result_each_list = domain_result_each_list))
}

SpatialReference<-function (cellType_distri, sample_information_region) 
{
    Feature_cluster <- function(expression_profile, sample_information) {
        sample_information <- sort(sample_information)
        expression_profile <- expression_profile[, names(sample_information)]
        num_each_class <- table(as.character(sample_information))
        feature_matrix <- matrix(0, length(num_each_class), nrow(expression_profile))
        row.names(feature_matrix) <- names(num_each_class)
        a = 1
        for (i in 1:(length(num_each_class))) {
            expression_profile_choose <- as.data.frame(expression_profile[, 
                a:(a + num_each_class[i] - 1)])
            a = a + num_each_class[i]
            feature_matrix[i, ] <- apply(expression_profile_choose, 
                1, mean)
        }
        colnames(feature_matrix)<-row.names(expression_profile)
        return(feature_matrix)
    }
    term_feature <- Feature_cluster(expression_profile = t(cellType_distri), 
        sample_information = sample_information_region)
   return(list(cellType_distri= cellType_distri, 
        term_feature = t(term_feature),  
        cellType = row.names(t(cellType_distri))))
}

SpatialQuery<-function (spatialReference_result, cellType_distri_query) 
{
    Get_query_hvg <- function(expression_profile, high_varGene_names) {
        missing_num <- length(high_varGene_names) - length(intersect(row.names(expression_profile), 
            high_varGene_names))
        missing_features <- setdiff(high_varGene_names, intersect(row.names(expression_profile), 
            high_varGene_names))
        missing_rate <- missing_num/length(high_varGene_names)
        if (missing_num > 0) {
            print(paste("The missing features in the query data is ", 
                missing_features, seq = ""))
            print(paste("The number of missing features in the query data is ", 
                missing_num, seq = ""))
            print(paste("The rate of missing features in the query data is ", 
                missing_rate, seq = ""))
            missing_data <- matrix(0, missing_num, ncol(expression_profile))
            row.names(missing_data) <- missing_features
            expression_profile <- rbind(expression_profile, missing_data)
        }
        expression_profile_hvg <- expression_profile[high_varGene_names, 
            ]
        return(expression_profile_hvg)
    }
    Query_result <- function(expression_profile_query_hvg, 
        feature_matrix) {
        require(philentropy)
        result <- data.frame(cluster_label = 1:ncol(expression_profile_query_hvg), cluster_distance = rep(0, ncol(expression_profile_query_hvg)))
        row.names(result) <- colnames(expression_profile_query_hvg)
        cor_result <- rep(0, ncol(feature_matrix))
        for(i in 1:nrow(result)){
        cor_result<-apply(feature_matrix,2,function(x){jensen_shannon(P=x,Q = expression_profile_query_hvg[,i],unit = 'log',testNA=T)})
        names(cor_result) <- colnames(feature_matrix)
        result[i, 1] <- names(sort(cor_result,decreasing=F))[1]
        result[i, 2] <- cor_result[result[i, 1]]
        }
        return(result)
    }
    cellType_distri_query_hvg <- Get_query_hvg(t(cellType_distri_query), 
        spatialReference_result$cellType)
    query_result <- Query_result(cellType_distri_query_hvg, 
        spatialReference_result$term_feature)
    query_result$Query_cell_id <- row.names(query_result)
    query_result$Predict_terms <- as.character(query_result$cluster_label)
    query_result$Distance <- query_result$cluster_distance
    query_result <- query_result[, c("Query_cell_id", 
        "Predict_terms", "Distance")]
    return(query_result)
}


