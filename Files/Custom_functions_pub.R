## Make labels for general plots
makeLab = function(x,pc) {
  paste0("PC",pc,": ",x,"% variance")
}


## Silence print/cat statements to console (from Hadley Wickham)
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
}


## Save pheatmap images similar to ggsave
save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


PCAplot <- function(data, genotypes, conditions){
  ## Calculating PC components
  pcs = prcomp(t(data), center = TRUE)
  percentVar = round(((pcs$sdev) ^ 2 / sum((pcs$sdev) ^ 2)* 100), 2)
  if(missing(conditions)){
    Lib.Id <- genotypes
    print(ggplot(as.data.frame(pcs$x), aes(PC1,PC2, color = Lib.Id, shape = Lib.Id), environment = environment()) +
            xlab(makeLab(percentVar[1],1)) + ylab(makeLab(percentVar[2],2)) + geom_point(size = 8) + theme_grey() +
            theme(legend.text = element_text(size = 16, face = "bold"),
                  legend.title = element_text(size = 16, colour = "black", face = "bold"),
                  plot.title = element_blank(),
                  axis.title = element_text(size = 18, face = "bold"),
                  axis.text.x = element_text(size = 16, face = "bold", color = "black"),
                  axis.text.y = element_text(size = 16, face = "bold", color = "black"),
                  plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")))
  }
  else{
    print(ggplot(as.data.frame(pcs$x), aes(PC1,PC2, color = genotypes, shape = conditions), environment = environment()) +
            xlab(makeLab(percentVar[1],1)) + ylab(makeLab(percentVar[2],2)) + geom_point(size = 8) + theme_grey() +
            theme(legend.text = element_text(size = 16, face = "bold"),
                  legend.title = element_text(size = 16, colour = "black", face = "bold"),
                  plot.title = element_blank(),
                  axis.title = element_text(size = 18, face = "bold"),
                  axis.text.x = element_text(size = 16, face = "bold", color = "black"),
                  axis.text.y = element_text(size = 16, face = "bold", color = "black"),
                  plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")))
  }
}


## Runs ChromXR as a single function
#   Useful `filename` = paste("~/test/Figures/State_",i,"_pheatmap.png", sep = "")
doChromXR <- function(os_type = "mac", chromhmm, dir = paste(chromhmm, "/MYOUTPUT/", sep = ""), statenum,
                      samplenames, conditions, hmapcol, mdscol = rev(hmapcol)) {
  orig.wd <- getwd()
  setwd(chromhmm)
  cat(paste("RUNNING ON:", os_type))
  cat("\nFirst generate the ChromState .txt files")
  if(os_type == "mac"){
    change_num <- paste0("sed -i \"\" ", "-r", " 's/num_states=", "[0-9]+","/num_states=", statenum, "/' chromXRloop.sh")
  }
  if(os_type == "linux"){
    change_num <- paste0("sed -i -r 's/num_states=", "[0-9]+","/num_states=", statenum, "/' chromXRloop.sh")
  }
  if(os_type == "win10"){
    change_num <- paste0("bash -c \"sed -i -r 's/num_states=", "[0-9]+","/num_states=", statenum, "/' chromXRloop.sh\"")
  }
  system(change_num)
  if (os_type == "win10"){
    cat("\nWorking... please wait a couple of minutes.")
    quiet(system('bash -c "./chromXRloop.sh"', intern = TRUE))
  } else{
    cat("\nWorking... please wait a couple of minutes.")
    quiet(system("./chromXRloop.sh", intern = TRUE))
  }
  cat("\nFinished segmenting genome!")
  setwd(dir)
  on.exit(expr = setwd(orig.wd))
  cat("\nNow we make the plots... this will take some time.")
  for (i in 1:statenum){
    i <<- i
    cat(paste0("\nState number = ", i))
    tryCatch({
      dat <- read.table(paste("CombinedMatrix-",i,"-10000bps.txt",sep =""),header = FALSE,sep = "\t",quote = "",
                        row.names = 1, na.strings = FALSE, stringsAsFactors = FALSE)
      colnames(dat) <- samplenames
      ## sample type labels for plot
      condition.type <- factor(conditions)
      ## filtering low variable regions
      ## PCA plot
      ## data is near binary, so probably not best to use PCA here
      dat.state1 <- varFilter(as.matrix(dat), var.cutoff = 0.5)
      p1  <- PCAplot(data = dat.state1, genotypes = condition.type)
      ggsave(filename = paste("State_",i,"_PCA.png", sep = ""), dpi = 600, width = 8, height = 6, path = getwd())
      ## MDS Plots ##
      d <- dist(t(dat.state1), method = "binary")
      mds.plot <- cmdscale(d = d, eig = T, k = 2)
      mdsDist <- data.frame(genotypes = condition.type, x = mds.plot$points[,1], y = mds.plot$points[,2])
      ggplot(mdsDist, aes(x = x, y = y, color = genotypes)) + geom_point(size = 8) +
        scale_color_manual(values = mdscol) +
        ylab("MDS Coordinate 2") + xlab("MDS Coordinate 1") + theme_grey() +
        theme(legend.text = element_text(size = 18, face = "bold"),
              legend.title = element_text(size = 18, colour = "black", face = "bold"),
              axis.title = element_text(size = 18, face = "bold"),
              axis.text.x = element_text(size = 18, face = "bold", color = "black"),
              axis.text.y = element_text(size = 18, face = "bold", color = "black"),
              plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
      ggsave(filename = paste("State_",i,"_MDS.png", sep = ""), dpi = 600, width = 8, height = 6, path = getwd())  
      ## Annotation Data frame
      df1 <- data.frame(Condition = condition.type)
      ##Heatmap
      # pheatmap
      annotation <- df1
      rownames(annotation) <- colnames(dat.state1)
      Condition <- hmapcol
      names(Condition) <- unique(conditions)
      anno_colors <- list(Condition = Condition)
      heatmap <- pheatmap(mat = log2(dat.state1+1), cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE,
                          show_colnames = TRUE, clustering_distance_rows = "correlation", clustering_distance_cols = "euclidean",
                          clustering_method = "complete", annotation_col = annotation, annotation_colors = anno_colors,
                          fontsize = 12, fontface="bold",
                          border_color="white")
      save_pheatmap_png(heatmap, filename = paste("State_",i,"_pheatmap.png", sep = ""))
    }, error = function(e){cat("\nERROR :", conditionMessage(e))})
  }
  cat(paste0("\nAll done! Figures are in ", dir))
}