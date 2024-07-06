#library
library("GenomicRanges")
library("IRanges")
library("S4Vectors")
library("BiocGenerics")
library("GenomeInfoDb")
library(Rgb)
library(GenomicRanges)
library(pbmcapply)
library(outliers)
library(dplyr)
library("stringr")
library("data.table")
#############################################################################
########################## Poly(A)Clust #####################################
#############################################################################
###### Cluster based on gene region
polyAClust<-function(polyA,cutoff=24){
  ################ forward process
  forward.PolyA<-polyA[polyA$strand=="+",]
  forward.PolyAClus<-data.frame()
  for (i in 1:length(unique(forward.PolyA$seqnames))){
    Forward.PolyA=forward.PolyA[forward.PolyA$seqnames==unique(forward.PolyA$seqnames)[i],]
    forward.geneid<-data.frame()
    for (i in 1:length(unique(Forward.PolyA$gene_id))){
      Forward.geneid=Forward.PolyA[Forward.PolyA$gene_id==unique(Forward.PolyA$gene_id)[i],]
      Forward.filter <- data.frame()

      # Loop until there are no rows left with the maximum score
      while (nrow(Forward.geneid) > 0) {
        Forward.idx <- Forward.geneid[which(Forward.geneid$score == max(Forward.geneid$score)), ][1,]
        Forward.split <- Forward.geneid[Forward.geneid$coord > Forward.idx$coord -cutoff & Forward.geneid$coord < Forward.idx$coord +cutoff, ]
        Forward.split1 <- data.frame(
          seqnames = as.character(unique(Forward.split$seqnames)),
          start = min(Forward.split$coord),
          end = max(Forward.split$coord),
          coord = Forward.idx$coord,
          coordRead = Forward.idx$score,
          PAnum = length(Forward.split$coord),
          tot_read = sum(Forward.split$score),
          strand = unique(Forward.idx$strand),
          geneid= as.character(unique(Forward.split$gene_id))
        )

        Forward.geneid <- Forward.geneid[!Forward.geneid$coord %in% Forward.split$coord, ]
        Forward.filter <- rbind(Forward.filter, Forward.split1)
        Forward.filter<-Forward.filter[order(Forward.filter$start),]
      }
      forward.geneid<-rbind(forward.geneid,Forward.filter)
    }
    forward.PolyAClus<-rbind(forward.PolyAClus,forward.geneid)
  }

  ########### Reverse process
  reverse.PolyA<-polyA[polyA$strand=="-",]
  reverse.PolyAClus<-data.frame()
  for (i in 1:length(unique(reverse.PolyA$seqnames))){
    Reverse.PolyA=reverse.PolyA[reverse.PolyA$seqnames==unique(forward.PolyA$seqnames)[i],]
    reverse.geneid<-data.frame()
    for (i in 1:length(unique(Reverse.PolyA$gene_id))){
      Reverse.geneid=Reverse.PolyA[Reverse.PolyA$gene_id==unique(Reverse.PolyA$gene_id)[i],]
      Reverse.filter <- data.frame()

      # Loop until there are no rows left with the maximum score
      while (nrow(Reverse.geneid) > 0) {
        Reverse.idx <- Reverse.geneid[which(Reverse.geneid$score == max(Reverse.geneid$score)), ][1,]
        Reverse.split <- Reverse.geneid[Reverse.geneid$coord > Reverse.idx$coord -cutoff & Reverse.geneid$coord < Reverse.idx$coord +cutoff, ]
        Reverse.split1 <- data.frame(
          seqnames = as.character(unique(Reverse.split$seqnames)),
          start = min(Reverse.split$coord),
          end = max(Reverse.split$coord),
          coord = Reverse.idx$coord,
          coordRead = Reverse.idx$score,
          PAnum = length(Reverse.split$coord),
          tot_read=sum(Reverse.split$score),
          strand = unique(Reverse.idx$strand),
          geneid= as.character(unique(Reverse.split$gene_id))
        )

        Reverse.geneid <- Reverse.geneid[!Reverse.geneid$coord %in% Reverse.split$coord, ]
        Reverse.filter <- rbind(Reverse.filter, Reverse.split1)
        Reverse.filter<- Reverse.filter[order( Reverse.filter$start),]
      }
      reverse.geneid<-rbind(reverse.geneid,Reverse.filter)
    }
    reverse.PolyAClus<-rbind(reverse.PolyAClus,reverse.geneid)
  }
  ########################### combine ######################################################
  polyAClus<-rbind(forward.PolyAClus,reverse.PolyAClus)
  return(polyAClus)
}
#########
#### Cluster including intergenic region
polyAClust_genome<-function(polyA,cutoff=24){
  ################ forward process
  forward.PolyA<-polyA[polyA$strand=="+",]
  forward.PolyAClus<-data.frame()
  for (i in 1:length(unique(forward.PolyA$seqnames))){
    Forward.PolyA=forward.PolyA[forward.PolyA$seqnames==unique(forward.PolyA$seqnames)[i],]
    Forward.filter <- data.frame()
    # Loop until there are no rows left with the maximum score
    while (nrow(Forward.PolyA) > 0) {
      Forward.idx <- Forward.PolyA[which(Forward.PolyA$score == max(Forward.PolyA$score)), ][1,]
      Forward.split <- Forward.PolyA[Forward.PolyA$coord > Forward.idx$coord -cutoff & Forward.PolyA$coord < Forward.idx$coord +cutoff, ]
      Forward.split1 <- data.frame(
        seqnames = as.character(unique(Forward.split$seqnames)),
        start = min(Forward.split$coord),
        end = max(Forward.split$coord),
        coord = Forward.idx$coord,
        coordRead = Forward.idx$score,
        PAnum = length(Forward.split$coord),
        tot_read = sum(Forward.split$score),
        strand = unique(Forward.idx$strand),
        geneid= as.character(unique(Forward.split$gene_id))
      )

      Forward.PolyA <- Forward.PolyA[!Forward.PolyA$coord %in% Forward.split$coord, ]
      Forward.filter <- rbind(Forward.filter, Forward.split1)
      Forward.filter<-Forward.filter[order(Forward.filter$start),]
    }
    forward.PolyAClus<-rbind(forward.PolyAClus,Forward.filter)
  }

  ########### Reverse process
  reverse.PolyA<-polyA[polyA$strand=="-",]
  reverse.PolyAClus<-data.frame()
  for (i in 1:length(unique(reverse.PolyA$seqnames))){
    Reverse.PolyA=reverse.PolyA[reverse.PolyA$seqnames==unique(forward.PolyA$seqnames)[i],]
    Reverse.filter <- data.frame()
    # Loop until there are no rows left with the maximum score
    while (nrow(Reverse.PolyA) > 0) {
      Reverse.idx <- Reverse.PolyA[which(Reverse.PolyA$score == max(Reverse.PolyA$score)), ][1,]
      Reverse.split <- Reverse.PolyA[Reverse.PolyA$coord > Reverse.idx$coord -cutoff & Reverse.PolyA$coord < Reverse.idx$coord +cutoff, ]
      Reverse.split1 <- data.frame(
        seqnames = as.character(unique(Reverse.split$seqnames)),
        start = min(Reverse.split$coord),
        end = max(Reverse.split$coord),
        coord = Reverse.idx$coord,
        coordRead = Reverse.idx$score,
        PAnum = length(Reverse.split$coord),
        tot_read=sum(Reverse.split$score),
        strand = unique(Reverse.idx$strand),
        geneid= as.character(unique(Reverse.split$gene_id))
      )

      Reverse.PolyA <-Reverse.PolyA[!Reverse.PolyA$coord %in% Reverse.split$coord, ]
      Reverse.filter <- rbind(Reverse.filter, Reverse.split1)
      Reverse.filter<- Reverse.filter[order( Reverse.filter$start),]
    }
    reverse.PolyAClus<-rbind(reverse.PolyAClus,Reverse.filter)
  }
  ########################### combine ######################################################
  polyAClus<-rbind(forward.PolyAClus,reverse.PolyAClus)
  return(polyAClus)
}
######################################################################################
###################### Final cluster process (Proposed)
polyAclust<-function(polyA,cutoff){
  removeIn.data<-polyA[polyA$type != "intergenic",]
  polyA_gr<-removeIn.data
  gene_clust<-polyAClust(polyA_gr,cutoff =8)
  polyA_genr<-polyA[!rownames(polyA)%in%rownames(polyA_gr),]
  interg_clust<-polyAClust_genome(polyA_genr,cutoff=8)
  #pro_genome$c1<-rep(1,nrow(pro_genome))
  final_clust<-rbind(gene_clust,interg_clust)
  return(final_clust)
}

#############################################################################
buildGenomicRanges <- function(seqname,position,score,strand = '*'){
  points.gr = GRanges(seqnames = seqname,ranges = IRanges(start = position,width = 1),
                      strand = strand,score = score)
  points.gr = sort(points.gr)
}
########
##############################################
#       simple clustering by distance        #
##############################################

simpleCluster <- function(polyA,max.gapwidth=24){
  points.gr= buildGenomicRanges(seqname = polyA$seqnames, position = polyA$coord,
                                strand = polyA$strand,score = polyA$score)
  # Cluster points by distance
  range.gr = reduce(points.gr,min.gapwidth=max.gapwidth,with.revmap=T,ignore.strand=FALSE)

  # Sum the score
  range.gr$tot_score = sum(extractList(points.gr$score,range.gr$revmap))

  # Select the center
  idx = BiocGenerics::which.max(extractList(points.gr$score,range.gr$revmap))
  range.gr$center = start(points.gr)[idx+c(0,cumsum(lengths(range.gr$revmap))[1:(length(range.gr$revmap)-1)])]
  range.gr$score = points.gr$score[idx+c(0,cumsum(lengths(range.gr$revmap))[1:(length(range.gr$revmap)-1)])]
  # Return result
  return(range.gr)
}

##############################################
#             findPeaks function             #
##############################################

findPeaks <- function(sub_pos,sub_wts,min_delta=24){

  # Declare
  group = pos = wts = NULL

  ND <- length(sub_pos)

  if(ND<=2){
    res <- tibble(group=1,start=min(sub_pos),end=max(sub_pos),sum.wts=sum(sub_wts),center=sub_pos[which.max(sub_wts)])
    return(res)
  }

  dc <- 0.5
  #paste('Computing Rho with gaussian kernel of radius: ',dc,collapse = '')

  rho = rep(0,ND)
  #
  # Gaussian kernel
  #
  for(i in 1:(ND-1)){
    for(j in (i+1):ND){
      tmp_dist = abs(sub_pos[i]-sub_pos[j])
      tmp_wts = sub_wts[i]*sub_wts[j]
      rho[i]=rho[i]+exp(-(tmp_dist/dc)^2)*sub_wts[j];
      rho[j]=rho[j]+exp(-(tmp_dist/dc)^2)*sub_wts[i];
    }
  }
  for(i in 1:ND){
    rho[i]=rho[i]+sub_wts[i];
  }

  maxd <- abs(sub_pos[1] - sub_pos[ND])
  rho_sorted = sort(rho,decreasing = T,index.return=T)
  ordrho <- rho_sorted$ix
  rho_sorted <- rho_sorted$x

  delta <- rep(-1.,ND)
  nneigh <- rep(0,ND)

  for(ii in 2:ND){
    delta[ordrho[ii]] = maxd
    for(jj in 1:(ii-1)){
      tmp_dist = abs(sub_pos[ordrho[ii]]-sub_pos[ordrho[jj]])
      if(tmp_dist<=delta[ordrho[ii]]){
        delta[ordrho[ii]] = tmp_dist
        nneigh[ordrho[ii]] = ordrho[jj]
      }
    }
  }
  delta[ordrho[1]]=max(delta)

  decision.data <- data.frame(rho=rho,delta=delta,rho.delta=rho*delta)
  outlier <- rho>max(rho)/2 & delta>max(delta)/2 & delta>min_delta

  tmp <- decision.data$rho.delta
  tmp[which(outlier)] <- 0
  outlier <- (outlier | scores(tmp, type="chisq", prob=0.99))& delta>min_delta

  if(!sum(outlier)){
    res <- tibble(group=1,start=min(sub_pos),end=max(sub_pos),sum.wts=sum(sub_wts),center=sub_pos[which.max(sub_wts)])
    return(res)
  }

  NCLUST = 0
  cl = rep(-1,ND)
  icl = c()
  for(i in 1:ND){
    if(outlier[i]){
      NCLUST = NCLUST + 1
      cl[i] = NCLUST
      icl[NCLUST] = i
    }
  }
  #paste('NUMBER OF CLUSTERS: ',NCLUST,collapse = '')

  # Assignation
  for(i in 1:ND){
    if (cl[ordrho[i]]==-1){
      cl[ordrho[i]] = cl[nneigh[ordrho[i]]];
    }
  }

  center <- rep(1,ND)
  center[icl] <- 2
  res <- data.frame(pos= sub_pos,wts = sub_wts,group = cl,center=center)
  res <- res %>% group_by(group) %>% summarise(start=min(pos),end=max(pos),sum.wts = sum(wts),center = pos[center==2],.groups = 'drop')
}

####################
##############
###################################################################################
################################## QuantifyR
Cluster.PolyA <- function(polyA, max.gapwidth = 24, mc.cores =1){
  points.gr= buildGenomicRanges(seqname = polyA$seqnames, position = polyA$coord,
                                strand = polyA$strand,score = polyA$score)
  # Check parameters.
  simple.clusters = simpleCluster(polyA, max.gapwidth = max.gapwidth)
  simple.clusters.df = as.data.frame(simple.clusters)
  simple.clusters.df$split_label = NA

  # Re-clustering of polyA sites with large width using weighted density peak calling algorithm.
  print('Re-clustering by weighted density peak clustering algorithm!')
  idx = which(simple.clusters.df$width > max.gapwidth)
  pos = extractList(points.gr@ranges@start,simple.clusters$revmap[idx])
  wts = extractList(points.gr$score,simple.clusters$revmap[idx])

  split.clusters = pbmcmapply(findPeaks,pos,wts,SIMPLIFY = F, mc.cores = mc.cores)

  lens = sapply(split.clusters,nrow)
  split.clusters.df = bind_rows(split.clusters[lens!=1], .id = "split_label")

  split.clusters.df$seqnames = rep(simple.clusters.df$seqnames[idx[lens!=1]],times = lens[lens!=1])
  split.clusters.df$strand = rep(simple.clusters.df$strand[idx[lens!=1]],times = lens[lens!=1])
  split.clusters.df$width = split.clusters.df$end - split.clusters.df$start + 1
  split.clusters.df = dplyr::rename(split.clusters.df,score=sum.wts)

  # Generate final clusters.
  polyA = simple.clusters.df[-idx[lens!=1],c("seqnames","start","end","width","strand","score","center")]
  polyA = rbind(polyA,split.clusters.df[,c("seqnames","start","end","width","strand","score","center")])
  return(polyA)
}
###################################################################################
################################ Peaklu ###########################################
###################################################################################

clusterByPeak <- function(tss.dt, peakDistance=100, localThreshold=0.02, extensionDistance=30) {
  # create copy for reference later
  copied.dt <- copy(tss.dt)
  setkey(tss.dt, pos)
  ##define variable as a NULL value
  pos = peak = ID = forward = reverse = V1 = V2 = chr = NULL
  # get peakID
  # TODO could potentially by optimized more
  peakID <- vapply(seq_len(tss.dt[,.N]), function(x) {
    id <- 0
    temp <- tss.dt[x,pos]
    ##ZL
    if(tss.dt[x,pos] == tss.dt[pos>temp-peakDistance & pos<temp+peakDistance,][which(tags == max(tags)),pos][1]){
      id <- x
    }else{id <- 0}
    return(id)
  }, numeric(1))

  # manipulate data.table to collapse clustered rows
  tss.dt[, peak := peakID]
  tss.dt[, ID := .I]
  ###############################################################################
  ##local filtering
  ###############################################################################
  unique_strands <- unique(tss.dt$strand)

  localF <- list()  # Initialize an empty list to store results

  for (strand in unique_strands) {
    if (strand == "+") {
      localF[[strand]] <- sapply(peakID[peakID > 0], function(i) {
        temp <- tss.dt[pos >= tss.dt$pos[i] & pos <= tss.dt$pos[i] + peakDistance & tss.dt$strand == strand, ]
        temp$ID[which(temp$tag < tss.dt$tags[i] * localThreshold)]
      })
    } else {
      localF[[strand]] <- sapply(peakID[peakID > 0], function(i) {
        temp <- tss.dt[pos >= tss.dt$pos[i] - peakDistance & pos <= tss.dt$pos[i] & tss.dt$strand == strand, ]
        temp$ID[which(temp$tag < tss.dt$tags[i] * localThreshold)]
      })
    }
  }

  # Combine results if needed
  localF<- c(localF$`+`, localF$`-`)
  # Filter indices to remove only those within the range of tss.dt
  indices_to_remove <- unlist(localF)
  indices_to_remove <- indices_to_remove[indices_to_remove <= nrow(tss.dt) & indices_to_remove > 0]
  indices_to_remove<-unique(indices_to_remove)
  # Filter out-of-range indices
  valid_indices <- indices_to_remove[indices_to_remove <= nrow(tss.dt) & indices_to_remove > 0]

  # Remove rows from tss.dt
  if(length(valid_indices) > 0) {
    tss.dt <- tss.dt[-valid_indices,]
  }
  ###############################################################################
  ###############################################################################
  tss.dt[, forward := ifelse(data.table::shift(pos,1,type="lead") < pos + extensionDistance, 1, 0)] #
  tss.dt[, reverse := ifelse(data.table::shift(pos,1,type="lag") > pos - extensionDistance, 1, 0)]
  tss.dt <- tss.dt[,list(peak=max(peak),start=min(pos),end=max(pos),tags=sum(tags)),by=list(rleid(peak, forward, reverse))]##ZL?

  # get start and end boundaries for clusters
  # TODO revisit this code for better optimization
  clusters <- lapply(as.list(tss.dt[peak>0,rleid]), function(x) {
    start <- tss.dt[x,start]
    end <- tss.dt[x,end]

    if (x-1>0 && tss.dt[x-1,!peak>0] && tss.dt[x-1,end] > start - extensionDistance) {
      start <- tss.dt[x-1,start]
      if (x-2>0 && tss.dt[x-2,!peak>0] && tss.dt[x-2,end] > start - extensionDistance) {
        start <- tss.dt[x-2,start]
      }
    }
    if (x+1<tss.dt[,.N] && tss.dt[x+1,!peak>0] && tss.dt[x+1,start] < end + extensionDistance) {
      end <- tss.dt[x+1,end]
      if (x+2<tss.dt[,.N] && tss.dt[x+2,!peak>0] && tss.dt[x+2,start] < end + extensionDistance) {
        end <- tss.dt[x+2,end]
      }
    }
    list(start, end)
  })

  clusters <- rbindlist(clusters)

  # deal with overlapping clusters here
  # TODO this section needs some more optimization/work

  rowVec <- which(clusters$V2 >= data.table::shift(clusters$V1,1,type="lead"))
  if (length(rowVec)>0) {
    ###############################################################################
    ###############################################################################
    for(i in seq_len(length(rowVec))){clusters$V1[rowVec[i]+1] = clusters$V1[rowVec[i]]}
    clusters <- clusters[-rowVec,]
  }##


  #  clusters <- unique(clusters)

  # get full clustering data
  # core promoter boundaries are calculated here (i.e. cumsum distribution)
  tss_clusters <- lapply(as.list(seq_len(clusters[,.N])), function(i) {
    start <- clusters[i,V1]
    end <- clusters[i,V2]
    #copied.dt[, ID := .I]##NEW April18
    cluster.data <- copied.dt[pos >= start & pos <= end, ]
    tags.sum <- cluster.data[,sum(tags)]##NEW Sep25, tags -> tags.sum
    q1 <- cluster.data[which(cumsum(tags) > 0.1*tags.sum),min(pos)]
    q9 <- cluster.data[order(-pos)][which(cumsum(tags) > 0.1*tags.sum),max(pos)]
    list(i
         ,cluster.data[,chr[[1]]]
         ,start
         ,end
         ,cluster.data[,strand[[1]]]
         #,cluster.data[which(ID %in% peakID), pos]  # NEW - use id column to find the intersection for the peakID vector (should hopefully only be 1!)
         ,cluster.data[which.max(tags),pos]
         ,tags.sum
         ,cluster.data[,max(tags)]
         ,q1
         ,q9
         ,q9 - q1 + 1)
  })

  # set names
  tss_clusters <- rbindlist(tss_clusters)
  setnames(tss_clusters, c( "cluster"
                            , "chr", "start", "end", "strand"
                            , "dominant_tss", "tags", "tags.dominant_tss"
                            , "q_0.1", "q_0.9", "interquantile_width" ))
  return(tss_clusters)
}

