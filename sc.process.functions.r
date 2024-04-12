suppressMessages(library(ROGUE))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))
suppressMessages(library(Seurat))
suppressMessages(library(SeuratObject))
suppressMessages(library(dplyr))
suppressMessages(library(doParallel))
suppressMessages(library(qvalue))
suppressMessages(library(DESeq2))
suppressMessages(library(MASS))

{
  
  calculate_rogue <- function(seu.exp, seu.meta, resolution,sample = 'orig.ident') {
    resolutions <- paste('RNA_snn_res.', resolution, sep = '')
    rogue.res <- rogue(seu.exp, labels = seu.meta[[resolutions]], samples = seu.meta[[sample]], platform = "UMI", span = 0.6)
    rogue.res[is.na(rogue.res)] = 0
    cluster.mean <- rowMeans(rogue.res) %>% as.data.frame() %>% t() %>% as.data.frame()
    colnames(cluster.mean) = 'Mean'
    cluster.mean$res <- resolutions
    return(cluster.mean)
  }
  
  find_best_res <- function(seu.exp, seu.meta, resolutions,sample = 'orig.ident', mc.cores=4) {
    time = Sys.time()
    resolutions.list <- mclapply(resolutions, function(i) calculate_rogue(seu.exp, seu.meta, i,sample), mc.cores = mc.cores)  # 设置核心数量
    resolutions_df <- do.call(rbind.data.frame, resolutions.list)
    max_values <- resolutions_df[which.max(resolutions_df$Mean), "res"]
    new.time = Sys.time()-time
    print(paste('Best res: ', max_values))
    print(paste('Time cost: ', new.time))
    return(resolutions_df)
  }
  
  train_model = function(srt, sce, cluster = 'majorCluster'){
    
    que <- as.SingleCellExperiment(srt)
    assay(sce, "logcounts") <- assay(sce,"exprs")
    assay(sce, "counts") <- assay(sce,"exprs")
    seu <- as.Seurat(sce)
    seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 5000)
    
    sample <- c()
    for (i in unique(seu@meta.data[[cluster]])) {
      if(length(which(seu@meta.data[[cluster]]==i))<=5000){
        sample <- c(sample, which(seu@meta.data[[cluster]] == i))
      }else{
        sample <- c(sample,sample(which(seu@meta.data[[cluster]] == i),5000))
      }
    }
    
    sce <- sce[,sample]
    exp <- assay(sce,'exprs')
    rowFilter <- rowSums(exp)
    colFilter <- colSums(exp)
    exp <- exp[rowFilter > 100,]
    exp <- t(exp) %>% as.matrix() %>%　as.data.frame()
    #exp <- apply(exp,2,as.numeric)
    exp <- exp[,colnames(exp) %in% rownames(que)]
    exp <- exp[,colnames(exp) %in% VariableFeatures(seu)]
    exp <- exp[,!(colnames(exp) %in% grep("^MT-|^RPL|^RPS|^IG|^TR", colnames(exp), value = T))]
    colnames(exp) <- gsub('-','_',colnames(exp))
    exp$cluster <- factor(colData(sce)[rownames(exp),cluster])
    
    #rf <- randomForest::randomForest(majorCluster~., data = exp)
    #saveRDS(rf,file= 'rf_rf.bycellType.rds')
    ranger <- ranger::ranger(cluster~., data = exp)
    
    return(ranger)
  }
  
  process_cluster <- function(srt, r = 'RNA_snn_res.0.2', ranger_model, mark,
                              prefix = 'Sub', cutoff = (1 / length(levels(ranger_model$predictions)))) {
    
    Idents(srt)<-srt@meta.data[[r]]
    if(!(is.null(mark))){
      
      mark = mark
      
    }else{ mark<-FindAllMarkers(srt, only.pos=T, logfc.threshold=0.5) }
    top <- mark %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
    
    res <- c()
    anot <- c()
    ratio = list()
    for(i in unique(srt@meta.data[[r]])){
      srt1 <- srt[, srt@meta.data[[r]] == i]
      que <- as.SingleCellExperiment(srt1)
      exp_pred <- assay(que, "logcounts") %>% t() %>% as.data.frame()
      colnames(exp_pred) <- gsub('-', '_', colnames(exp_pred))
      pred <- predict(ranger_model, exp_pred)
      
      #pt = max(table(pred$predictions))/sum(table(pred$predictions))
      pt = as.data.frame.array(table(pred$predictions)/sum(table(pred$predictions)))
      colnames(pt) = 'ratio'
      pt$sub = names(table(pred$predictions))
      pt$res = i
      print(paste(i,": ", names(which.max(table(pred$predictions)))," ",max(pt$ratio), sep = ''))
      ratio[[i]] = pt
      #cutoff = (1 / length(levels(ranger_model$predictions)))
      
      # if(pt > cutoff){
      #   clu = names(which.max(table(pred$predictions)))
      # }else{
      #   gene = top[top$cluster == i,][['gene']]
      #   clu = paste(prefix,'.','C',i,'.',gene,sep = '')
      # }
      # 
      # res = c(res, i)
      # anot = c(anot, clu)
      
    }
    
    #clusters = data.frame(res = res, predict = anot)
    return(ratio)
  }
  
}


{
  
  diff_NB = function(data,group='NGroup') {
    res1 = data.frame(matrix(nrow = 0, ncol = 6))
    cells = unique(data$sub)
    data$group = data[[group]]
    for (j in 1:length(cells)){
      dat.test = data %>% filter(sub==cells[j]) %>% 
        dplyr::select(sub_num,group,sf.DESeq)
      # only 1 response group, not comparable
      if (length(unique(dat.test$group)) == 1){
        next
      } else {
        # mod.out <- glm.nb(sub_num~NGroup+offset(sf.DESeq),
        #                             data=dat.test)
        mod.out <- tryCatch({glm.nb(sub_num~group+offset(sf.DESeq),
                                    data=dat.test)},
                            error = function(e) {
                              glm(sub_num~group+offset(sf.DESeq),
                                  data=dat.test,family=poisson)})
        mod.summary <- summary(mod.out)
        coeff.name <- grep("^group",rownames(mod.summary$coefficients),
                           value=T)
        p.value <- mod.summary$coefficients[coeff.name,"Pr(>|z|)"]
        logFC <- -(mod.summary$coefficients[coeff.name, "Estimate"]/log(2))
        res.val = c(cells[j], logFC, p.value, 
                    nrow(dat.test[dat.test$group==sort(
                      unique(dat.test$group))[1],]),
                    nrow(dat.test[dat.test$group==sort(
                      unique(dat.test$group))[2],]))
        res1=rbind(res1, res.val)
      }}
    
    g1 = paste('n',sort(unique(dat.test$group))[1],sep = '')
    g2 = paste('n',sort(unique(dat.test$group))[2],sep = '')
    colnames(res1) =c('sub','LFC','p',g1,g2)
    res1$sub = as.character(res1$sub)
    #res1$dataset = as.character(res1$dataset)
    res1$p = as.numeric(res1$p)
    res1$LFC = as.numeric(res1$LFC)
    res1[[g1]] = as.numeric(res1[[g1]])
    res1[[g2]] = as.numeric(res1[[g2]])
    p = res1$p
    qobj = qvalue(p=p)
    qvalues = qobj$qvalues
    res1$qvalues = qvalues
    print(sort(unique(dat.test$group)))
    
    #res1$weight = -qnorm(res1$qvalues/2) * res1$LFC 
    #print(res1$weight)
    #print(res1)
    return(res1)
  }
  
}