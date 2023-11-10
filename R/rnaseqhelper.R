## Mehmet Tekman, 2023

#' @importFrom tidyverse group_by


library(tidyverse)
library(RColorBrewer)

gois_list <- list(
  bra.down = c("Msgn1", "Osr1", "Rspo3", "Fgf8", "Wnt3a"),
  eo.down.repr = c("Dmdx1", "Dpf3", "Foxa1", "Hey1", "Hhex", "Tcf7l2", "Tle2"),
  eo.down = c("Mesp1", "Foxa2", "Sox17", "Lhx1", "Cer1"),
  pp.up = c("Lefty1", "Lefty2", "Nodal", "Wnt8a", "Fgf5", "Otx2", "Cldn6"),
  episc.up = c("Nanog", "Pou5f1", "Sox2", "L1td1", "Utf1"),
  ne.up = c("Sox1", "Sox3", "Olig3", "Gbx2", "Pou3f1", "Msx3", "Slc7a3",
            "Zic1", "Zic2", "Nkx1-2", "Epha2", "Efnb1", "Nkx6-2"),
  me.down = c("Pdgfra", "Foxc1", "Foxc2", "Isl1", "Kdr",
              "Mixl1", "Hand1", "Myh6"),
  hox.genes = c("Hoxa1", "Hoxa10", "Hoxa11", "Hoxa11os", "Hoxa13", "Hoxa2",
                "Hoxa3", "Hoxa4", "Hoxa5", "Hoxa6", "Hoxa7", "Hoxa9", "Hoxaas2",
                "Hoxaas3", "Hoxb1", "Hoxb13", "Hoxb2", "Hoxb3", "Hoxb3os",
                "Hoxb4", "Hoxb5", "Hoxb5os", "Hoxb6", "Hoxb7", "Hoxb8", "Hoxb9",
                "Hoxc10", "Hoxc11", "Hoxc12", "Hoxc13", "Hoxc4", "Hoxc5",
                "Hoxc6", "Hoxc8", "Hoxc9", "Hoxd1", "Hoxd10", "Hoxd11",
                "Hoxd12", "Hoxd13", "Hoxd3", "Hoxd3os1", "Hoxd4", "Hoxd8",
                "Hoxd9"),
  main.genes = c("T", "Eomes", "Wnt3", "Wnt3a", "Nodal", "Fgf8", "Mixl1",
                 "EGFP"),
  ac.1 = c("Fgf8", "Eomes", "Foxa2", "Meis1", "Otx2", "Lhx1",
           "Mixl1", "Cyp26a1", "T", "Cdx1", "Cdx2", "Axin2", "Rspo3", "Wnt3a"),
  ac.2 = c("Axin2", "Cdx1", "Cdx2", "Cyp26a1", "Dkk1", "Meis1", "Mesp1"),
  other1 = c("Foxc2", "Gsc", "Mixl1", "Foxc1", "Pdgfra", "Kdr", "Tbx1",
             "Amot", "Eya1", "Prrx1", "Lhx1", "Tnnt1", "Isl1", "Tbx20",
             "Myh7", "Myh6", "Eya2", "Tbx18", "Snai1", "Snai2", "Hand1",
             "Hand2", "Mesp1", "Pax2", "Mesp2", "Col2a1", "Myocd", "Pax3",
             "Wt1", "Dll3", "Prickle1", "Msgn1", "Rspo3", "Osr1", "Twist1",
             "Twist2", "Tbx6", "Meox1", "Gata6", "Gata4", "Foxa2", "Sox17",
             "Lama1", "Tgfa", "Fn1", "Cer1", "Tgfb2", "Bmper", "Chrd", "Lgr5",
             "Bmp2", "Fzd7", "Wnt3a", "Bmp7", "Cfc1", "Dkk1", "Fgf1",
             "Tdgf1", "Wnt3", "Tgfb1", "Nodal", "Dact1", "Bmp4", "Dll1",
             "Hey1", "Hhex", "Id1", "Dmbx1", "Dpf3", "Sall1", "Foxa1",
             "Tcf7l2", "Hesx1", "Tcf7l1", "Hdac7", "Otx2", "Fbn2", "Otx1"),
  other2 <- c("Six3", "Fabp7", "Ntrk2", "Zic1", "Mab21l2", "Nefl", "Olig3",
              "Pou4f2", "Nkx1-2", "Sox3", "Sema3c", "Epha2", "Zic3", "Efnb1",
              "Radil", "Syt11", "Slc7a3", "Nes", "Zic5", "Snph", "Nfasc",
              "Pou4f1", "Elavl3", "Nrcam", "Grin1", "Gfra3", "Phox2a", "Nkx6-2",
              "Pax6", "Nsg2", "Msx3", "Ephb1", "Ncan", "Nova2", "Zic2", "Grik3",
              "Epha1", "Bcl11a", "Hoxa2", "Tubb3", "Sox1", "Neurod1", "Neurog1",
              "Stmn2", "Atxn1", "Cntn2", "Neurl1a", "Sema3c", "Gap43", "Fgf5",
              "Tbx3", "Cldn6", "Pou3f1", "Lefty2", "Gbx2", "Nanog", "L1td1",
              "Wnt8a", "Pou5f1", "Lefty1", "Utf1", "Esrrb", "Sox2", "Lin28a",
              "Dnmt3a", "Cbx7", "Kdm2b", "Atf7ip2")
)


generate_labelling_table <- function(tab, tabn) {
  dat.labs <- tabn %>% group_by(cluster) %>%
    summarise(n=n_distinct(gene)) %>%
    mutate(ftext = paste0("cluster ", cluster, ", ", n, " genes"))
  
  ## Sanity check that all factors are in summary
  miss = which(!(levels(tabn$cluster) %in% dat.labs$cluster))
  if (length(miss)>0) {
    temp = cbind(cluster=miss, n=0, ftext=paste0("cluster ", miss, ", 0 genes"))
    dat.labs <- rbind(dat.labs, temp)
  }
  dat.labs <- dat.labs %>% arrange(cluster)
  dat.labs
  return(dat.labs)
}


clusterGenePlots <- function(tab, score.thresh=0, out.dir="gene_cluster") {
  tabn = tab %>% filter(score >= score.thresh)
  ## Keep all clusters even if they're empty
  tabn$cluster = factor(tabn$cluster, levels=unique(sort(tab$cluster)))
  
  ## yes
  dat.labs = generate_labelling_table(tab, tabn)

  time.breaks = sort(unique(sort(tab$time)))

  
  ## labelling function for facet
  lab.fun <- function(s) { (dat.labs %>% filter(cluster == s))$ftext }
  
  p1 = tabn %>% ggplot(aes(x=time, y=value, fill=condition, 
                           colour=condition, group=condition)) +  
    stat_summary(fun.data = 'mean_sdl', geom = 'ribbon', alpha=0.1, colour=NA) +
    stat_summary(fun = mean, geom = 'line') + 
    facet_wrap("cluster", labeller=labeller(cluster=lab.fun), drop = FALSE) +  
    ylab("Scaled Expression") +
    ggtitle("Gene Trends by cluster", 
            subtitle = paste0("Genes (", length(unique(tabn$gene)), 
                              ") with cluster affinity scores >= ", 
                              score.thresh)) +
    scale_x_continuous(breaks=time.breaks) +
    theme_bw()

  newoutprefix = file.path(out.dir, paste0("gene_plots-k", numClusters(tab), "-score", score.thresh))
  ggplot2::ggsave(plot=p1, filename=paste0(newoutprefix, ".svg"),
                  dpi=800, width = 10, height = 10, unit = "in")

  tabn %>% write_tsv(paste0(newoutprefix, ".tsv"))
  return(p1)
}

numClusters <- function(tab) {
  return(length(unique(tab$cluster)))
}

clusterGenePlotsWithVaryingScore <- function(
                                             tab, score_list=c(0, 0.5, 0.9, 0.99), out.dir="gene_cluster")
{
  library(patchwork)

  if (!dir.exists(out.dir)) {
    dir.create(out.dir)
  }

  res = lapply(score_list, function(x) {
    return(clusterGenePlots(tab, score.thresh=x, out.dir=out.dir))
  })

  newoutprefix = file.path(out.dir, 
                           paste0("gene_plots-k", numClusters(tab)))

  rplot = wrap_plots(res, ncol=2, guides = "collect") & 
    theme(plot.title = element_blank(), axis.title = element_blank())
  ggplot2::ggsave(plot=rplot, filename=paste0(newoutprefix, "_montage.svg"),
                  dpi=800, width = 15, height = 15, unit = "in")

  return(rplot)
} 

highQualityGenes <- function(sam.mat,
                             smallestGroupSize = 3,
                             detectability = 10) {
  ## Generate a list of high quality genes from
  ## an average sample wide-format matrix
  
  res = rowSums(
    sam.mat >= detectability
  ) >= smallestGroupSize
  keep.genes = res[res == TRUE]
  drop.genes = res[res == FALSE]
  write_tsv(as.data.frame(res), 
            paste0("smallestGroup",
                   smallestGroupSize, "-detected",
                   detectability, "_keep.genes"))
  message(paste0("Dropped: ", length(drop.genes),
                 " genes (", 
                 as.integer(100 * length(drop.genes)/ length(res)), "%)"))
  return(keep.genes)
}

GenePlotsByGene <- function(norm.long, scale.long, gois_list, 
                            outprefix="gene.lists", out.dir="gene_lists") {
  ## For a normalised long matrix and a scaled long matrix, following column pattern of:
  ## - gene	condition	time	replicate	value
  ## As well as a list of genes in the pattern on:
  ## - list(mygeneset1 = c("Eomes", "T"),
  ##        mygeneset2 = c("Lefty2", etc, etc.)

  if (!dir.exists(out.dir)) {
    dir.create(out.dir)
  }
  
  plotDims <- function(n) {
    ## Calculates plot and width height for a given N plots
    w = ceiling(sqrt(n))
    h = ceiling(n / w)
    return(list(w=w,h=h))
  }

  time.breaks = sort(unique(sort(norm.long$time)))

  pgene.list <- lapply(names(gois_list), function(glist_name) {
    glist <- gois_list[[glist_name]]
    genes.found = unique((norm.long %>% filter(gene %in% glist))$gene)

    if (length(genes.found) < 1) {
      message("no genes found for: ", glist_name)
      return(NULL)
    }
    
    pgene.norm <- norm.long %>% filter(gene %in% genes.found) %>%
      ggplot(aes(x=time, y=value, 
                 fill=condition, colour=condition, group=condition)) +
      stat_summary(fun.data = 'mean_sdl', geom = 'ribbon', alpha=0.1, colour=NA) +
      stat_summary(fun = mean, geom = 'line') +
      geom_jitter(width=0.1) +
      scale_y_log10() +
      facet_wrap("gene", scales="free_y") +
      ylab("Log10 Normalised Expression") +
      ggtitle("Normalised Expression Time Plots, grouped by Genes of Interest") +
      scale_x_continuous(breaks=time.breaks) +
      theme_bw()

    pdims = plotDims(length(genes.found))
    ggsave(plot=pgene.norm, 
           filename=file.path(out.dir, paste0(outprefix, ".", glist_name, ".normalised.svg")),
           dpi=800, width = pdims$w * 2, height = pdims$h * 1.5, unit = "in")

    pgene.scale <- scale.long %>% filter(gene %in% genes.found) %>%
      ggplot(aes(x=time, y=value, 
                 fill=condition, colour=condition, group=condition)) +
      stat_summary(fun.data = 'mean_sdl', geom = 'ribbon', alpha=0.1, colour=NA) +
      stat_summary(fun = mean, geom = 'line') +
      geom_jitter(width=0.1) +
      facet_wrap("gene", scales="free_y") +
      ylab("Scaled Expression") +
      ggtitle("Scaled Expression Time Plots, grouped by Genes of Interest") +
      scale_x_continuous(breaks=time.breaks) +
      theme_bw()

    ggsave(plot=pgene.scale, 
           filename=file.path(out.dir, paste0(outprefix, ".", glist_name, ".scaled.svg")),
           dpi=800, width = pdims$w * 2, height = pdims$h * 1.5, unit = "in")

    return(list(norm=pgene.norm, scale=pgene.scale))
  })
}


## We try Josephus's clustering
doKMeans <- function(data, k, plot_list, out.dir="heatmaps_k", 
                     highlight.genes=NULL, ...) {
  ## Perform KMeans clustering
  ## data:  a normalised (not-scaled) wide-format matrix.
  ## k: number of partitions
  ## plot_list: a list of column names for specifying the order of columns
  condition_list = plot_list[[1]]

  out.dir = paste0(out.dir, k)

  if (!dir.exists(out.dir)) {
    dir.create(out.dir)
  }
  
  output_svg_prefix = file.path(out.dir, paste0("genes.k", k, "."))
  
  scaledata <- t(scale(t(data)))
  scaledata <- scaledata[complete.cases(scaledata), ] ## Remove non-zero
  ##gc()
  kClust <- kmeans(scaledata, centers = k, nstart = 1000, iter.max = 30)
  kClusters <- kClust$cluster

  clust.centroid <- function(i, dat, clusters) {
    ##message("-[    clust.centroid ]-")
    ind = (clusters == i)
    colMeans(dat[ind, ])
  }
  kClustcentroids <- sapply(levels(factor(kClusters)),
                            clust.centroid, scaledata, kClusters) ## is a matrix

###Kmolten <- reshape::melt(kClustcentroids)
  Kmolten = as.data.frame(as.table(kClustcentroids))
  colnames(Kmolten) <- c("sample", "cluster", "value")
  ## ensure correct factorizing
  Kmolten$sample <- factor(Kmolten$sample)

  clusters_unique <- sort(unique(Kmolten$cluster))
  cores <- list()
  Kmolten.list <- list()
  scores <- list()
  corescores <- list()

  for (i in clusters_unique) {
    core_i <- Kmolten[Kmolten$cluster == i, ]
    core_i$sample <- factor(core_i$sample, levels = condition_list)
    K_i <- scaledata[kClusters == i, ]
    corescore_i <- function(x) cor(x, core_i$value)
    score_i <- apply(K_i, 1, corescore_i)
    ##K_i_molten <- reshape::melt(K_i)
    K_i_molten <- as.data.frame(as.table(K_i))
    colnames(K_i_molten) <- c("gene", "sample", "value")
    K_i_molten <- merge(K_i_molten, score_i, by.x = "gene", by.y = "row.names", all.x = T)
    colnames(K_i_molten) <- c('gene', 'sample', 'value', 'score')
    ##K_i_molten <- K_i_molten[order(K_i_molten$score), ]
    K_i_molten <- K_i_molten[order(K_i_molten$value), ]
    Kmolten.list[[i]] <- K_i_molten
    scores[[i]] <- score_i
    corescores[[i]] <- corescore_i
    cores[[i]] <- core_i
  }
  ## name collected values
  names(Kmolten.list) <- paste("K", 1:k, "molten", sep = '')
  names(cores) <- paste("core", 1:k, sep = '')
  names(scores) <- paste("score", 1:k, sep = '')

  scaledata.k <-cbind(scaledata,
                      cluster = kClust$cluster,
                      score = Reduce(`c`, scores)[rownames(scaledata)])
  scaledata.k <- scaledata.k[order(scaledata.k[, "cluster"], scaledata.k[, "score"]), ]

  select_df <- function(df, val, col.selector) {
    df[ df[, col.selector] == val, ]
  }
  df2dflist <- function(df, col.selector) { # actually it is split()
    message("-[df2dflist]-")
    col.vals <- unique(df[, col.selector])
    dfl <- lapply(seq(col.vals),
                  function(i) select_df(df,
                                        val = col.vals[i],
                                        col.selector))
    names(dfl) <- col.vals
    dfl
  }
  scaledata_list <- df2dflist(scaledata.k, "cluster")
  names(scaledata_list) <- paste("cluster", 1:k, sep = "")

  gaps.idxs <- cumsum(table(scaledata.k[, "cluster"])) # scaledata.k is one matrix! only column 'cluster'

  show.genes = rownames(scaledata)
  if (!is.null(highlight.genes)) {
    ## Blank out any non-listed genes
    show.genes[!(show.genes %in% highlight.genes)] = ""
  }
  
  singleHeatmap <- function(ordered_cols, outfix="") {
    library(pheatmap)

    main=paste0(nrow(scaledata), " genes, clustered into k=", 
                k, " ordered by ", outfix)
    
    colors.kr <- colorRampPalette(c("black", "#bb0000"))(100)
    ##colors.kr <- brewer.pal(10,"PuOr")
    svg(tempfile(fileext = ".svg"))
    ph <- pheatmap::pheatmap(
                      scaledata.k[, ordered_cols],
                      cluster_rows = F,
                      cluster_cols = F,
                      show_rownames = T,
                      labels_row=show.genes,
                      cellwidth = 40,
                      col = colors.kr,
                      fontsize_row = 8,
                      border_color = NA,
                      gaps_row = gaps.idxs,          # gap after each block
                      main = main,
                      ...
                    )
    dev.off()
    
    svg(paste0(output_svg_prefix, outfix, ".svg"))
    graphics::plot.new()
    ##print(ph)
    print(better_pheatmap(ph))
    dev.off()
  }

  for (name in names(plot_list)) {
    message("-[singleHeatmap]- ", name)
    singleHeatmap(plot_list[[name]], name)
  }
  return(scaledata.k)
}

getFillD <- function(rh) {
  getatts = function (x) attributes(x)$names

  ## If it has fill and dimension, we modify it 
  if ("gp" %in% getatts(rh)) {
    if ("fill" %in% getatts(rh$gp)) {
      hasdim = dim(rh$gp$fill)
      if (!is.null(hasdim))
      {
        message("Changed: ", rh, " ", hasdim)
        rh$gp$col <- rh$gp$fill
        rh$gp$lwd <- 3
      }
    }
  }

  if ("children" %in% getatts(rh)) {
    for (child in names(rh$children)) {
      rh$children[[child]] = getFillD(rh$children[[child]])
    }
  }
}


better_heatmap <- function(ph)
{
  lwdd = 2
  message("-[better_pheatmap]-")
  ## taken from
  ##https://stackoverflow.com/questions/44318690/no-border-color-in-pheatmap-in-r
  library(purrr)

  if (class(ph) == "Heatmap") {
    ## ComplexHeatmap
    ph = grid.grabExpr(draw(ph))
    hmaps = ph$children
    for (hm in names(hmaps)) {
      ## Get all Rects
      ##rects = grep("rect", names(ph$children[[hm]]$children), value=T)
      rects = names(ph$children[[hm]]$children)

      if (length(rects) > 0) {
        ##rects = names(ph$children[[hm]]$children)
        for (rr in rects) {
          ## Check for rects with tabular data
          hasdim = dim(ph$children[[hm]]$children[[rr]]$gp$fill)
          if (!is.null(hasdim))
          {
            message("--Here5")
            ph$children[[hm]]$children[[rr]]$gp$col <- ph$children[[hm]]$children[[rr]]$gp$fill
            ph$children[[hm]]$children[[rr]]$gp$lwd <- lwdd
          }
        }
      }
      else 
      {
        if ("gp" %in% names(ph$children[[hm]])) {
          ## Check for rects with tabular data
          hasdim = dim(ph$children[[hm]]$gp$fill)
          if (!is.null(hasdim))
          {
            message("Setting Cols for: ", hm)
            ph$children[[hm]]$gp$col <- ph$children[[hm]]$gp$fill
            ph$children[[hm]]$gp$lwd <- lwdd
          }
        }
      }
    }
    return(ph)
  }
  else {
    message("Here5")
    ## PHEATMAP
    ## Extract the right grob
    grob_classes <- purrr::map(ph$gtable$grobs, class)
    idx_grob <- which(purrr::map_lgl(grob_classes, function(cl) 'gTree' %in% cl))[1]
    grob_names <- names(ph$gtable$grobs[[idx_grob]]$children)
    idx_rect <- grob_names[grep('rect', grob_names)][1]
    
    ## Remove borders around cells
    ph$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$col <- ph$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$fill
    ph$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$lwd <- 3    
    return(ph)
  }
}


runDESeq <- function(tab, phenotype_data, min.detect=50, min.samples=8)
{
  library(tidyverse)
  library(DESeq2)
  
  sub_as = (tab[rowSums(tab > min.detect) > min.samples,])
  print(dim(sub_as))

  ddsObj = DESeqDataSetFromMatrix(
    countData = sub_as,
    colData = phenotype_data,
    design = ~Condition
  )
  ddsObj = DESeq(ddsObj)

  vsd1 <- vst(ddsObj, blind=FALSE)
  p1 = plotPCA(vsd1, intgroup=c("Condition"),returnData =F )
  vsd2 <- vst(ddsObj, blind=TRUE)
  p2 = plotPCA(vsd2, intgroup=c("Condition"),returnData =F )
  return(list(tab = tab, phenotype = phenotype_data,
              bFalse = p1, bTrue = p2, 
              vFalse = vsd1, vTrue = vsd2,
              ddsObj=ddsObj, sub=sub_as))
}

pcaAndMatrices <- function(res, out.dir="deseq2")
{
  ## Initial input matrices
  dir.create(out.dir, recursive=T, showWarnings=FALSE)

  res$tab %>% as.data.frame %>% 
    rownames_to_column("gene") %>% write_tsv(file.path(out.dir, "input_matrix.tsv"))
  res$phenotype %>% write_tsv(file.path(out.dir, "phenotype_data.tsv"))
  
  saveRDS(res, file.path(out.dir, "deseq2obj.rds"))
  
  res$sub %>% as.data.frame %>% 
    rownames_to_column("gene") %>%    
    write_tsv(file.path(out.dir, "input_matrix.filt.tsv"))
  
  counts(res$ddsObj, normalized=TRUE) %>% as.data.frame %>% rownames_to_column("gene") %>%
    write_tsv(file.path(out.dir, "input_matrix.filt.normalized.tsv"))
  
  assay(res$vFalse) %>% as.data.frame %>% rownames_to_column("gene") %>%
    write_tsv(file.path(out.dir, "input_matrix.filt.normalized.vst_corrected.tsv"))
  
  pdf(file.path(out.dir, "input_matrix.filt.normalized.vst_corrected_PCA.pdf"), 
      width=7, height=7)
  print(res$bFalse)
  print(res$bTrue)
  dev.off()
}

volcanoPlot <- function(dsqres, degenes, title, 
                        curve=list(sd=0.15, sc=10,offset=1), 
                        curve.scale=1, curve.show=FALSE, ylim=NULL)
{
  library(patchwork)
  library(ggrepel)
  options(repr.plot.height = 12, repr.plot.width = 12) 

  ## Extract relevant info from DESeq results
  ana = dsqres %>% 
    select(c(gene, log2FoldChange, padj)) %>% 
    mutate(mLog10Padj = -log10(padj)) %>%
    arrange(desc(mLog10Padj))

  ## The main function -- do not edit?
  cust.fun <- function(x, sd=0.15, sc=10, offset=1) {
    offset + dnorm(x, mean = 0, sd=sd) * sc
  }

  ## The slope function -- highlight these genes
  cfun = function(x) {
    return(cust.fun(x,
                    sd=curve$sd * curve.scale,
                    sc=curve$sc * curve.scale,
                    offset=curve$offset))
  }

  ## We highlight genes in the zoomed zone fitting the curve, but the 
  ## main DE genes are shown as shapes and higlighted
  red = ana %>%  
    mutate(isTopN = gene %in% degenes) %>%
    mutate(highlight = isTopN | (mLog10Padj > cust.fun(log2FoldChange)))
  ##mutate(isTopN = mLog10Padj > cfun(log2FoldChange))

  max_x = max(abs(ana$log2FoldChange)) + 0.05 ##symmetry
  
  plot1 = red %>% 
    ggplot(aes(x=log2FoldChange, 
               y=mLog10Padj, 
               colour=highlight,
               shape=isTopN,
               label=gene)) + 
    geom_point() + 
    scale_colour_manual(values = c("TRUE"  = "red",
                                   "FALSE" = "grey")) +
    scale_shape_manual(values =  c("TRUE"  =  5,
                                   "FALSE" = 19)) + 
    scale_x_continuous(lim=c(-max_x,max_x), 
                       breaks = waiver(),
                       n.breaks = 10) +
    geom_label_repel(
      data=red %>% filter(highlight == TRUE) %>% head(15),
      box.padding = 0.5,
      max.overlaps = 30,
      ##max.overlaps = 20,
      colour = "black") +
    ggtitle(title)

  if (curve.show) {
    plot1 = plot1 + geom_function(fun = cfun, n = 100, colour = "blue")
  }
  if (!is.null(ylim)) {
    plot1 = plot1 + scale_y_continuous(lim=ylim)
  }    
  return(plot1)
}



pairwiseHMapAndVolcano <- function(ddsObj,
                                   transformed.counts = NULL, 
                                   numer = "FRT", denom = "TAF2", 
                                   top.ngenes.tocluster = 2000,
                                   top.ngenes.tohighlight = 50,
                                   lFC.zoom = 1.5, pAdj.zoom = 20,
                                   kmeans = 2, out.dirprefix=".")
{
  ntitle = paste0(numer, " vs ", denom)
  outdir = file.path(out.dirprefix,
                     tolower(gsub("[^A-Za-z0-9]", "_", ntitle)))

  dir.create(outdir, recursive=T, showWarnings=FALSE)

  ## 1. Perform RNA-seq contrast between numerator and denominator
  dsqres <- DESeq2::results(
                      ddsObj,
                      contrast = c("Condition", numer, denom),
                      cooksCutoff = Inf,
                      independentFiltering = FALSE) %>% 
    as.data.frame %>% rownames_to_column("gene") %>%
    mutate(mLog10Padj = -log10(padj))

  norm.counts = counts(ddsObj,normalized=TRUE)

  top.genes.tocluster = (dsqres %>% arrange(desc(mLog10Padj)) %>% 
                         head(top.ngenes.tocluster))$gene
  top.genes.tohighlight = (dsqres %>% arrange(desc(mLog10Padj)) %>% 
                           head(top.ngenes.tohighlight))$gene

  write_tsv(data.frame(top.genes.tocluster=top.genes.tocluster), 
            file.path(outdir, paste0("clustered_genes.top",
                                     top.ngenes.tocluster, ".tsv")))
  write_tsv(data.frame(top.genes.tohighlight = top.genes.tohighlight), 
            file.path(outdir, paste0("volcano_genes.tophighlight",
                                     top.ngenes.tohighlight, ".tsv")))
  
  sample.columns = c(grep(paste0("^", numer), colnames(norm.counts), value = T),
                     grep(paste0("^", denom), colnames(norm.counts), value = T))
  
  ## Volcano Plots
  p1 = volcanoPlot(dsqres,
                   top.genes.tohighlight,
                   ntitle, 
                   curve=list(sd=0.3, sc=60, offset=10), 
                   curve.show=F)
  ## Volcano Plots zoomed in
  p2 = volcanoPlot(dsqres %>% filter(abs(log2FoldChange) < lFC.zoom & 
                                     mLog10Padj < pAdj.zoom), 
                   top.genes.tohighlight,
                   paste0(ntitle, " (zoomed)"),
                   curve=list(sd=0.25, sc=5, offset=8), ylim=c(0,22), 
                   curve.show=F)
  
  hmaps = wrap_plots(list(p1,p2), ncol=1, guides = "collect") & 
    theme_bw() + 
    theme(legend.position="none")

  volcano.svg = file.path(outdir, "volcano_pairwise.svg")
  deseq2.out = file.path(outdir, "deseq2_results.tsv")
  
  svg(volcano.svg, width=8, height=9)
  print(hmaps)
  dev.off()
  message("Saved Volcano: ", volcano.svg)

  dsqres %>% 
    mutate(isTopN.gene = gene %in% top.genes.tohighlight) %>% 
    write_tsv(deseq2.out)
  message("Saved DESeq2 Results: ", deseq2.out)

  for (kmk in kmeans) {
    message("Calculating k=", kmk)

    ## Pairwise Heatmap of Normalised and Transformed Counts
    niceKMeansHeatmap_NormAndTrans(
      norm.counts, transformed.counts, 
      genes.tocluster = top.genes.tocluster, 
      sample.columns = sample.columns,
      genes.tohighlight = top.genes.tohighlight,
      dsqres, kmk,
      out.dir = file.path(outdir, paste0("kmeans", kmk)), 
      heatprefix = "heatmap_pairwise",
      ntitle = "Heatmap Pairwise")
    
    ## Global Heatmap of Normalised and Transformed Counts
    niceKMeansHeatmap_NormAndTrans(
      norm.counts, transformed.counts, 
      genes.tocluster = top.genes.tocluster, 
      sample.columns = NULL,
      genes.tohighlight = top.genes.tohighlight,
      dsqres, kmk,
      out.dir = file.path(outdir, paste0("kmeans", kmk)), 
      heatprefix = "heatmap_all",
      ntitle = "Heatmap All")
  }
}


niceKMeansHeatmap_NormAndTrans <- function(norms, trans, 
                                           genes.tocluster=NULL, 
                                           sample.columns=NULL, 
                                           genes.tohighlight=NULL, 
                                           dsqres, 
                                           kmeans, out.dir, 
                                           heatprefix, ntitle) {
  library(tidyverse)

  if (is.null(genes.tocluster)) {
    message("using all genes in normalised matrix for clustering")
    genes.tocluster = rownames(norms)
  }
  if (is.null(sample.columns)) {
    message("using all samples in normalised matrix for clustering")
    sample.columns = colnames(norms)
  }
  
  ## Heatmaps
  res.dsqres = niceKMeansHeatmap(
    norms[genes.tocluster, sample.columns], 
    k = kmeans, 
    out.dir = out.dir,
    heatprefix = heatprefix,
    prefix.title = paste0(ntitle, " :"),
    highlight.genes = genes.tohighlight
  )
  ## Merge Norm cluster
  ##message("SAVING")
  ##saveRDS(list(dsqres, res.dsqres), 
  ##        paste0(format(Sys.time(), "%H-%M-%S--"), "BUMM.rds"))
  
  dsq.dsq = left_join(dsqres, 
                      res.dsqres$clusters,
                      by = c("gene" = "gene")) %>% 
    mutate(norm_cluster = case_when(
             is.na(cluster) ~ "not in norm heatmap",
             TRUE ~ cluster)) %>% select(-cluster)

  ## Heatmaps using Corrected Normalized Counts
  if (!is.null(trans)) {
    message("Using transformed counts too")

    ## Heatmaps
    res.dsqres.trans = niceKMeansHeatmap(
      trans[genes.tocluster, sample.columns], 
      k = kmeans, 
      out.dir = out.dir,
      heatprefix = paste0(heatprefix, ".vst_corrected"),
      prefix.title = paste0(ntitle, " (vst corrected) :"),
      highlight.genes = genes.tohighlight
    )
    ## Merge Trans cluster
    dsq.dsq = left_join(dsq.dsq, 
                        res.dsqres.trans$clusters,
                        by = c("gene" = "gene")) %>% 
      mutate(trans_cluster = case_when(
               is.na(cluster) ~ "not in trans heatmap",
               TRUE ~ cluster)) %>% select(-cluster)
  }

  save.cluster = file.path(out.dir, 
                           paste0("deseq2.results.cluster.k", kmeans, ".tsv"))
  write_tsv(dsq.dsq, save.cluster)
  message("Saved DESeq k", kmeans, "fill:rgb(33.72549%,0%,0%);: ", save.cluster)
}


niceKMeansHeatmap <- function(
                              norm_counts, k, out.dir="heatmaps_k",
                              heatprefix = "heatmap", prefix.title = "", 
                              highlight.genes=NULL,
                              width.in = 6, height.in = 6)
{
  ## normalised counts should be in the correct sample order already
  suppressPackageStartupMessages(library(ComplexHeatmap))
  options(repr.plot.height = height.in, repr.plot.width = width.in)

  if (!dir.exists(out.dir)) {
    dir.create(out.dir)
  }
  
  if (is.null(highlight.genes)) {
    ## If no genes given, show top N
    top.genes = head(names(sort(rowMeans(norm_counts), 
                                decreasing=T)), 30)
    top.title = paste0(nrow(norm_counts), " DE genes, top ",
                       length(top.genes), " highlighted")
  } else {
    top.genes = highlight.genes
    top.title = paste0(nrow(norm_counts), " DE genes, ",
                       length(top.genes), " highlighted")
  }
  
  scaledata <- t(scale(t(norm_counts)))
  scaledata <- scaledata[complete.cases(scaledata), ] ## Remove non-zero

  ha = rowAnnotation(foo = anno_mark(
                       at = which(rownames(scaledata) %in% top.genes), 
                       labels = top.genes))

  if (k > 5) {
    rtitle = "%s"
  } else {
    rtitle = "Clust %s"
  }

  hm_now = Heatmap(scaledata, 
                   column_title = paste0(prefix.title, "  ", top.title),
                   row_km = k,
                   ##border = TRUE,
                   cluster_row_slices = FALSE, ## Arrange the K? NO
                   cluster_rows = TRUE, ## Prevents annotation bunching if TRUE
                   show_row_dend = FALSE,
                   row_gap = unit(3, "mm"),
                   name = "mat", 
                   row_title = rtitle,
                   col = colorRampPalette(c("black", "#bb0000"))(100),
                   cluster_columns = FALSE,
                   show_row_names = FALSE,
                   right_annotation = ha,
                   column_names_rot = 45,
                   row_names_gp = gpar(fontsize = 4))

  hm_now_drawn = draw(hm_now)

  ## Get cluster assignments
  cluster.assignments = (function() {
    cluster.list = lapply(row_order(hm_now_drawn), 
                          function(x) rownames(scaledata)[x])
    cluster.table = do.call(rbind, lapply(names(cluster.list), 
                                          function(n) {
                                            data.frame(gene = cluster.list[[n]], cluster=n)
                                          }))
    return(cluster.table)
  })()

  ## Unite cluster assignments with Norm, Scaled, and Pvalue
  output_prefix = file.path(out.dir, paste0(heatprefix, ".k", k, "."))

  save.svg = paste0(output_prefix, "svg")
  save.svgold = paste0(output_prefix, "_old.svg")
  save.pdf = paste0(output_prefix, "pdf")
  save.png = paste0(output_prefix, "png")
  
  save.scale = paste0(output_prefix, "scale.tsv")
  save.norm = paste0(output_prefix, "norm.tsv")
  save.cluster = paste0(output_prefix, "clusters.tsv")

  ## Hack to remove the weird grey lines
  ## TODO: Fix, scaling issue...
###bhm = better_heatmap(hm_now)
  
  ##svg(save.svgold, width = width.in, height = height.in)
  ##print(hm_now)
  ##dev.off()

  svg(save.svg, width = width.in, height = height.in)
  ##grid.draw(bhm)
  print(hm_now)
  dev.off()

  pdf(save.pdf, width = width.in, height = height.in)
  ##grid.draw(bhm)
  ##grid.newpage()
  print(hm_now)
  dev.off()
  
  write_tsv(as.data.frame(norm_counts) %>% rownames_to_column("gene"),
            save.norm)
  message("Saved Norm: ", save.norm)

  write_tsv(as.data.frame(scaledata) %>% rownames_to_column("gene"), 
            save.scale)
  message("Saved Scale: ", save.norm)

  return(list(plot = hm_now, clusters=cluster.assignments, scaled=scaledata))
}
