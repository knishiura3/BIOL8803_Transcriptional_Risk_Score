library(ggplot2)
library(grid)
library(gridExtra) 
library(gtable)

# ########################################################################
# ## Regional association plots                                         ##
# ##                                                                    ##
# ## James Staley                                                       ##
# ## University of Bristol                                              ## 
# ## Email: james.staley@bristol.ac.uk                                  ##
# ##                                                                    ##
# ## 19/02/19                                                           ##
# ########################################################################


##########################################################
##### Legend #####
##########################################################

#' g_legend
#'
#' g_legend 
#' @param gplot a ggplot
#' @import ggplot2 grid gridExtra gtable
#' @author James R Staley <james.staley@bristol.ac.uk>
#' @export
g_legend<-function(gplot){
  tmp <- ggplot_gtable(ggplot_build(gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


##########################################################
##### Add legend to regional association plots #####
##########################################################

#' add_g_legend
#'
#' add_g_legend
#' @param g a ggplot
#' @param legend legend
#' @import ggplot2 grid gridExtra gtable
#' @author James R Staley <james.staley@bristol.ac.uk>
#' @export
add_g_legend <- function(g, legend){
  lheight <- sum(legend$height)*1.5
  g <- grid.arrange(g, legend, ncol = 1, heights = unit.c(unit(1, "npc") - lheight, lheight))
  return(g)
}

##########################################################
##### Stacked regional association plot #####
##########################################################

#' stack_assoc_plot
#'
#' stack_assoc_plot plots stacked regional association plots
#' @param markers data.frame of markers with markername (marker), chromosome (chr) and position (pos) 
#' @param z matrix of Z-scores or probabilities with one column for each trait
#' @param corr matrix of correlation statistics between markers
#' @param corr.top correlation statistics between the top marker and the rest of the markers
#' @param traits trait names
#' @param ylab the y-axis label
#' @param type the type of the plot either log10p or probabilities
#' @param x.min start of region
#' @param x.max end of region
#' @param top.marker the top associated marker, i.e. the marker with the largest -log10p or probability
#' @param legend add r2 legend
#' @import ggplot2 grid gridExtra gtable
#' @author James R Staley <james.staley@bristol.ac.uk>
#' @export
stack_assoc_plot_custom <- function(markers, z, corr=NULL, corr.top=NULL, traits, ylab=NULL, type="log10p", x.min=NULL, x.max=NULL, top.marker=NULL, legend=TRUE){
  
  # Error messages
  if(!(type=="log10p" | type=="prob")) stop("the type of plot has to be either log10p or prob")
  if(length(traits)!=ncol(z)) stop("the number of traits is not the same as the number of columns for the Z-scores")
  if(nrow(markers)!=nrow(z)) stop("the number of markers is not the same as the number of rows for the Z-scores")
  if(!is.null(corr)){if(ncol(corr)!=nrow(markers) | nrow(corr)!=nrow(markers)) stop("corr has to have the same dimensions as the number of rows in the markers dataset")}
  if(!is.null(corr.top)){if(length(corr.top)!=nrow(markers)) stop("corr.top has to have the same length as the number of rows in the markers dataset")}
  # if(any(rownames(corr)!=markers$marker)) stop("corr has to have the same markers in the same order as the markers dataset")
  if(any(names(markers)!=c("marker", "chr", "pos"))) stop("dataset needs to include marker, chr and pos columns in that order")
  if(length(unique(markers$chr))>1) stop("there should only be markers from one chromosome in the markers dataset")   
  if(!(markers$chr[1] %in% 1:22)) stop("the plotting tool is only for autosomal chromosomes") 
  if(any(is.na(markers))) stop("there are missing markers in your marker dataset") 
  # if(any(is.na(z))) stop("there are missing values in the Z-score matrix")
  if(class(markers$pos)!="integer") stop("the pos variable has to be an integer")
  if(is.null(corr) & !is.null(corr.top) & is.null(top.marker)) stop("top.marker must be defined if corr.top is provided")
  if(is.null(corr) & !is.null(corr.top)){if(length(corr.top)!=nrow(markers)) stop("corr.top has to have the same length as the number of rows in the markers dataset")}
  if(!is.null(top.marker) & length(which(top.marker==markers$marker))==0) stop("top.marker is not contained in the markers dataset")
  if(!is.null(top.marker) & length(which(top.marker==markers$marker))>1) stop("top.marker maps to multiple markers in the markers dataset")
  
  # Coerce data
  markers$marker <- as.character(markers$marker)
  chr <- as.integer(markers$chr[1])
  r2_legend <- legend
  if(is.null(x.min)){x.min <- min(as.integer(markers$pos))}
  if(is.null(x.max)){x.max <- max(as.integer(markers$pos))}
  if((x.max - x.min)>10000000) stop("the plotting tool can plot a maximum of 10MB")
  
  # mlog10p
  if(type=="log10p"){
    mlog10p <- suppressWarnings(apply(z, 2, function(x){-(log(2) + pnorm(-abs(x), log.p=T))/log(10)}))
    # Kenji Nishiura:  don't impose a -log10(pval) ceiling 
    # mlog10p[mlog10p>1000 & !is.na(mlog10p)] <- 1000
  }
  if(type=="log10p"){ylab <- expression("-log"["10"]*paste("(",italic("p"),")"))}else{if(is.null(ylab)){ylab <- "Probability"}}
 
  # Genes
  gene.region <- gassocplot::genes[gassocplot::genes$chr==chr & !(gassocplot::genes$end<x.min) & !(gassocplot::genes$start>x.max),]
  gene.region$start[gene.region$start<x.min] <- x.min
  gene.region$end[gene.region$end>x.max] <- x.max
  gene.region <- gene.region[with(gene.region, order(start)), ]
  ngenes <- nrow(gene.region)

  # Max and min
  x.min <- x.min - 0.02*(x.max - x.min)
  x.max <- x.max + 0.02*(x.max - x.min)
  
  # Correlation matrix
  if(is.null(corr) & is.null(corr.top)){r2_legend <- FALSE; corr <- matrix(NA, nrow=nrow(markers), ncol=nrow(markers))}

  # Recombination plot
  recombination.plot <- plot_recombination_rate_stack(chr, x.min, x.max)

  # Gene plot
  if(ngenes==0){gene.plot <- plot_gene_zero(chr, x.min, x.max, stack=TRUE)}
  if(ngenes>0 & ngenes<=5){gene.plot <- plot_gene_two(gene.region, chr, x.min, x.max, stack=TRUE)}
  if(ngenes>5 & ngenes<=10){gene.plot <- plot_gene_five(gene.region, chr, x.min, x.max, stack=TRUE)}
  if(ngenes>10 & ngenes<=25){gene.plot <- plot_gene_ten(gene.region, chr, x.min, x.max, stack=TRUE)}
  if(ngenes>25){gene.plot <- plot_gene_fifteen(gene.region, chr, x.min, x.max, stack=TRUE)}

  # Top marker
  if(length(top.marker)!=0){if(is.na(top.marker)){top.marker <- NULL}}
    
  # Association plot
  for(i in length(traits):1){
    if(type=="log10p"){
      data <- data.frame(marker=markers$marker, chr=as.integer(markers$chr), pos=as.integer(markers$pos), stats=mlog10p[,i], stringsAsFactors=F)
    }else{
      data <- data.frame(marker=markers$marker, chr=as.integer(markers$chr), pos=as.integer(markers$pos), stats=z[,i], stringsAsFactors=F)    
    }
    marker.plot <- plot_assoc_stack(data, corr, corr.top, x.min, x.max, top.marker, ylab, type)
    legend <- g_legend(marker.plot)
    if(i==length(traits)){g <- plot_regional_gene_assoc(recombination.plot, marker.plot, gene.plot, traits[i], ngenes)}
    if(i<length(traits)){
      g1 <- plot_regional_assoc(recombination.plot, marker.plot, traits[i])
      g <- gtable:::rbind_gtable(g1, g, "last")
      panels <- g$layout$t[grep("panel", g$layout$name)]
      g$heights[panels[1]] <- unit(3,"null") 
    } 
  }

  # Combined plot
  if(r2_legend==T){
    combined.plot <- add_g_legend(g, legend)
  }else{
    combined.plot <- g
  }

  return(combined.plot)
}