wkdir <- "~/cey/tDNA/NLR_T-DNASeq_distr"
setwd(wkdir)

gene_path <- "gene_TAIR10_T-DNASeq.bed"
exon_path <- "exon_mergeByGene_TAIR10_T-DNASeq.bed"
cds_path <- "CDS_mergeByGene_TAIR10_T-DNASeq.bed"
tdna_path <- "../data/T-DNASeq.bed"

gene_dat <- read.table(gene_path,sep="\t")
exon_dat <- read.table(exon_path,sep="\t")
cds_dat <- read.table(cds_path,sep="\t")
tdna_dat <- read.table(tdna_path,sep="\t")

gene_dat <- gene_dat[,c(1:3,5,7,9,11,12,14)]
exon_dat <- cbind(exon_dat[,c(1:3)],rep("exon",nrow(exon_dat)),exon_dat[,c(4,6,7,9)])
cds_dat <- cbind(cds_dat[,c(1:3)],rep("CDS",nrow(cds_dat)),cds_dat[,c(4,6,7,9)])
tdna_dat <- tdna_dat[,c(1:3,5)]

colnames(gene_dat) <- c("chrom","chromStart","chromEnd","feature","strand","geneID","tDNAStart","tDNAEnd","tDNA")
colnames(exon_dat) <- c("chrom","chromStart","chromEnd","feature","geneID","tDNAStart","tDNAEnd","tDNA")
colnames(cds_dat) <- c("chrom","chromStart","chromEnd","feature","geneID","tDNAStart","tDNAEnd","tDNA")
colnames(tdna_dat) <- c("chrom","chromStart","chromEnd","tDNA")

library(Gviz)
library(dplyr)

nlrs_geneID <- c(n="../NLRs/NBS.txt",
                 nl="../NLRs/NBS-LRR.txt",
                 cn="../NLRs/CC-NBS.txt",
                 cnl="../NLRs/CC-NBS-LRR.txt",
                 tn="../NLRs/TIR-NBS.txt",
                 tnl="../NLRs/TIR-NBS-LRR.txt")
nlrs_raw <- list()
for (name in names(nlrs_geneID)){
    nlrs_raw[[name]] <- toupper(as.vector(read.table(nlrs_geneID[[name]])[,1]))
}

filter_genes <- function(data,geneIDs){
    filtered <- list()

    filter_inner <- function(feature){
        filter_temp <- data.frame(matrix(ncol=ncol(data[[feature]])+1,nrow=0))
        for (name in names(geneIDs)){
            temp <- filter(data[[feature]],geneID %in% geneIDs[[name]])
            temp <- cbind(temp,group=rep(name,nrow(temp)))
            filter_temp <- rbind(filter_temp,temp)
        }
        return(filter_temp)
    }

    filtered[["gene"]] <- filter_inner("gene")
    filtered[["exon"]] <- filter_inner("exon")
    filtered[["cds"]] <- filter_inner("cds")
    return(filtered)
}

filtered_nlrs <- filter_genes(list("gene"=gene_dat,"exon"=exon_dat,"cds"=cds_dat),
                              nlrs_raw)
filtered_all <- list("gene"=cbind(gene_dat,group=rep(".",nrow(gene_dat))),
                     "exon"=cbind(exon_dat,group=rep(".",nrow(exon_dat))),
                     "cds"=cbind(cds_dat,group=rep(".",nrow(cds_dat))))

data <- filtered_all

plot_tDNA_distr(filtered_nlrs,"cds",nlrs_raw[["tnl"]])
plot_tDNA_distr(filtered_nlrs,"exon",nlrs_raw[["cnl"]])
plot_tDNA_distr(filtered_nlrs,"cds",c("AT3G44400","AT3G44480","AT3G44630","AT3G44670")) ## DM2

plot_cluster(filtered_nlrs,"chr4",9488400,9565600,"DM8")
plot_cluster(filtered_all,"chr4",9488400,9565600,"DM8")
## https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005990
plot_cluster(filtered_nlrs,"chr3",16044000,16097000,"DM2-1",showTDNA=F)
plot_cluster(filtered_all,"chr3",16044000,16097000,"DM2-1",showTDNA=F)
plot_cluster(filtered_nlrs,"chr3",16164900,16235000,"DM2-2")
plot_cluster(filtered_all,"chr3",16164900,16235000,"DM2-2")




plot_tDNA_distr <- function(data,feature,genes){
    ## get entries for genes in "genes"
    sel_raw <- list()
    sel_raw[["gene"]] <- filter(data[["gene"]],geneID %in% genes)
    sel_raw[["exon"]] <- inner_join(filter(data[["exon"]],geneID %in% genes),
                                    select(sel_raw[["gene"]],geneID,strand),by="geneID")
    sel_raw[["cds"]] <- inner_join(filter(data[["cds"]],geneID %in% genes),
                                   select(sel_raw[["gene"]],geneID,strand),by="geneID")
    sel_raw[["tDNA"]] <- filter(sel_raw[["gene"]],tDNA != ".")

    ## reformat chromStart and chromEnd so they're relative to the start of CDS
    sel <- list()
    sel[["feat"]] <- data.frame(matrix(ncol=ncol(sel_raw[[feature]]),nrow=0))
    sel[["cds"]] <- data.frame(matrix(ncol=ncol(sel_raw[["cds"]]),nrow=0))
    sel[["tDNA"]] <- data.frame(matrix(ncol=ncol(sel_raw[["tDNA"]]),nrow=0))
    for (currID in unique(sel_raw[[feature]]$geneID)){
        temp_feat <- filter(sel_raw[[feature]],geneID==currID)
        temp_cds <- filter(sel_raw[["cds"]],geneID==currID)
        temp_tDNA <- filter(sel_raw[["tDNA"]],geneID==currID)
        strand <- temp_feat$strand[[1]]
        if (strand == "+") {
            start <- min(temp_cds$chromStart)
            temp_feat <- mutate(temp_feat,chromStart=chromStart-start,
                                chromEnd=chromEnd-start)
            temp_cds <- mutate(temp_cds,chromStart=chromStart-start,
                               chromEnd=chromEnd-start)
            temp_tDNA <- mutate(temp_tDNA,chromStart=tDNAStart-start,
                                chromEnd=tDNAEnd-start)
        } else {
            start <- max(temp_cds$chromEnd)
            cstart <- temp_feat$chromStart
            cend <- temp_feat$chromEnd
            temp_feat <- mutate(temp_feat,chromStart=-(cend-start),
                                chromEnd=-(cstart-start),
                                strand=rep("+",nrow(temp_feat)))
            cstart <- temp_cds$chromStart
            cend <- temp_cds$chromEnd
            temp_cds <- mutate(temp_cds,chromStart=-(cend-start),
                               chromEnd=-(cstart-start),
                               strand=rep("+",nrow(temp_cds)))
            temp_tDNA <- mutate(temp_tDNA,chromStart=-c(tDNAEnd-start),
                                chromEnd=-c(tDNAStart-start),
                                strand=rep("+",nrow(temp_tDNA)))
        }
        sel[["feat"]] <- rbind(sel[["feat"]],temp_feat)
        sel[["cds"]] <- rbind(sel[["cds"]],temp_cds)
        sel[["tDNA"]] <- rbind(sel[["tDNA"]],temp_tDNA)
    }

    ## move all features to the right so that everything is >=0
    subtract_min <- function(df_ls){
        min_coord <- min(c(df_ls[["feat"]]$chromStart,
                           df_ls[["cds"]]$chromStart,df_ls[["tDNA"]]$chromStart))
        out <- list()
        out[["feat"]] <- mutate(df_ls[["feat"]],chromStart=chromStart-min_coord,
                                chromEnd=chromEnd-min_coord)
        out[["cds"]] <- mutate(df_ls[["cds"]],chromStart=chromStart-min_coord,
                               chromEnd=chromEnd-min_coord)
        out[["tDNA"]] <- mutate(df_ls[["tDNA"]],chromStart=chromStart-min_coord,
                                chromEnd=chromEnd-min_coord)
        return(list(shift=min_coord,data=out))
    }
    ## sel2 <- subtract_min(tnl)
    ## sel_5 <- as.character(unique(sel_raw[["gene"]]$geneID)[1:5])
    sel_all <- as.character(unique(sel_raw[["gene"]]$geneID))
    ## sel_sub <- lapply(sel2,function(x) {filter(x,geneID %in% sel_5)})
    sel_sub_stats <- subtract_min(sel)
    sel_sub <- sel_sub_stats[["data"]]

    gtrack <- GenomeAxisTrack()

    tracks <- list()
    for (currSel in sel_all){
        gene_feat <- filter(sel_sub[["feat"]],geneID==currSel)
        gene_cds <- filter(sel_sub[["cds"]],geneID==currSel)
        gene_tDNA <- filter(sel_sub[["tDNA"]],geneID==currSel)
        track1 <- AnnotationTrack(start=gene_feat$chromStart,end=gene_feat$chromEnd,
                                  chromosome="chr0",group=gene_feat$geneID,
                                  strand=gene_feat$strand,shape="smallArrow",
                                  fill="grey",col="dark grey",showId=T,name=currSel)
        trackcds <- AnnotationTrack(start=gene_cds$chromStart,end=gene_cds$chromEnd,
                                  chromosome="chr0",group=gene_cds$geneID,
                                  strand=gene_cds$strand,shape="box",
                                  fill="orange",col="dark orange",showId=F,name=currSel)
        if (nrow(gene_tDNA) > 0){
            track2 <- AnnotationTrack(start=gene_tDNA$chromStart,end=gene_tDNA$chromEnd,
                                      chromosome="chr0",name=currSel,group=gene_tDNA$tDNA,shape="box",
                                      col="red",fill="red",showId=F,stacking="dense",alpha=0.8)
            tracks[[currSel]] <- OverlayTrack(trackList = list(track1,trackcds,track2),name=currSel)
        } else {
            tracks[[currSel]] <- OverlayTrack(trackList = list(track1,trackcds),name=currSel)
        }
    }
    htrack <- HighlightTrack(trackList=gtrack,start=-sel_sub_stats[["shift"]],width=3)
    plotTracks(append(tracks,htrack),groupAnnotation="group",name=feature)
}

## Example
plot_tDNA_distr(filtered_nlrs,"cds",nlrs_raw[["tnl"]])
plot_tDNA_distr(filtered_nlrs,"cds",c("AT3G44400","AT3G44480","AT3G44630","AT3G44670")) ## DM2



plot_cluster <- function(data,chr,start,end,cluster_name="",showTDNA=F){
    ## extract NLRs in DM8 cluster
    cluster <- list()
    cluster[["gene"]] <- filter(data[["gene"]],start<=chromStart & chromStart<=end & chrom==chr)
    cluster[["exon"]] <- inner_join(filter(data[["exon"]],start<=chromStart & chromStart<=end & chrom==chr),
                                    select(cluster[["gene"]],geneID,strand),by="geneID")
    cluster[["cds"]] <- inner_join(filter(data[["cds"]],start<=chromStart & chromStart<=end & chrom==chr),
                                   select(cluster[["gene"]],geneID,strand),by="geneID")
    cluster[["tDNA"]] <- filter(tdna_dat,start<=chromStart & chromStart<=end & chrom==chr)

    gtrack <- GenomeAxisTrack()

    gene <- data.frame(distinct(cluster[["gene"]][,c(1:3,5,6)]))
    gene <- arrange(mutate(gene,strand=as.character(strand)),chromStart)
    track0 <- AnnotationTrack(start=gene$chromStart,end=gene$chromEnd,
                              chromosome=gene$chrom,strand=gene$strand,
                              group=gene$geneID,name=cluster_name,
                              fill="grey",col="dark grey",showId=T,cex.group=1)

    gene_ex <- data.frame(distinct(cluster[["exon"]][,c(1:3,5,10)]))
    ## gene_ex <- format_arrows(gene_ex)
    track1 <- AnnotationTrack(start=gene_ex$chromStart,end=gene_ex$chromEnd,
                              chromosome=gene_ex$chrom,strand=gene_ex$strand,
                              group=gene_ex$geneID,name="Exon",shape="smallArrow",
                              fill="dark green",col="dark green",showId=F)

    gene_cds <- data.frame(distinct(cluster[["cds"]][,c(1:3,5,10)]))
    ## gene_cds <- format_arrows(gene_cds)
    track2 <- AnnotationTrack(start=gene_cds$chromStart,end=gene_cds$chromEnd,
                              chromosome=gene_cds$chrom,group=gene_cds$geneID,
                              strand=gene_cds$strand,name="CDS",shape="box",
                              fill="orange",col="dark orange",showId=F)
    plotTracks(list(gtrack,track0,track1,track2),groupAnnotation="group")

    gene_tDNA <- distinct(select(filter(cluster[["gene"]],tDNA!="."),chrom,tDNAStart,tDNAEnd,tDNA))
    if (showTDNA) {
        track3 <- AnnotationTrack(start=cluster[["tDNA"]]$chromStart,end=cluster[["tDNA"]]$chromEnd,
                                  chromosome=cluster[["tDNA"]]$chrom,name="T-DNA",
                                  group=cluster[["tDNA"]]$tDNA,col="red",fill="red",
                                  showId=T,stackHeight=0.5)
    } else {
        track3 <- AnnotationTrack(start=cluster[["tDNA"]]$chromStart,end=cluster[["tDNA"]]$chromEnd,
                              chromosome=cluster[["tDNA"]]$chrom,name="T-DNA",group=cluster[["tDNA"]]$tDNA,
                              col="red",fill="red",showId=showTDNA,stacking="dense",alpha=0.8)
    }
    ot <- OverlayTrack(trackList = list(track1,track2),name="Exon,CDS")
    plotTracks(list(gtrack,track0,ot,track3),groupAnnotation="group",fontsize=15)
}

## Example
plot_cluster(filtered_nlrs,"chr4",9488400,9565600,"DM8")
## https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005990
plot_cluster(filtered_nlrs,"chr3",16044000,16097000,"DM2a",showTDNA=F)
plot_cluster(filtered_nlrs,"chr3",16164900,16235000,"DM2b")



## figure out arrows
format_arrows <- function(df){
    df <- arrange(mutate(df,strand=as.character(strand)),chromStart)
    for (row in seq(nrow(df))){
        if (((df[row,"strand"] == "+") &
             (row == nrow(df) | !identical(df[row+1,"geneID"],df[row,"geneID"]))) |
            ((df[row,]$strand == "-") &
             (row == 1 | !identical(df[row-1,"geneID"],df[row,"geneID"])))){
            next
        } else {
            df[row,"strand"] <- "*"
        }
    }
    return(df)
}

## ## get all TNLs
## tnl_raw <- list()
## tnl_raw[["gene"]] <- filter(nlrs_gene,geneID %in% nlrs_raw[["tnl"]])
## tnl_raw[["exon"]] <- inner_join(filter(nlrs_exon,geneID %in% nlrs_raw[["tnl"]]),
##                                 select(tnl_raw[["gene"]],geneID,strand),by="geneID")
## tnl_raw[["cds"]] <- inner_join(filter(nlrs_exon,geneID %in% nlrs_raw[["tnl"]]),
##                                select(tnl_raw[["gene"]],geneID,strand),by="geneID")
## tnl_raw[["tDNA"]] <- filter(tnl_raw[["gene"]],tDNA != ".")

## ## reformat chromStart and chromEnd so they're relative to the start of CDS
## tnl <- list()
## tnl[["gene"]] <- data.frame(matrix(ncol=ncol(tnl_raw[["gene"]]),nrow=0))
## tnl[["exon"]] <- data.frame(matrix(ncol=ncol(tnl_raw[["exon"]]),nrow=0))
## tnl[["cds"]] <- data.frame(matrix(ncol=ncol(tnl_raw[["cds"]]),nrow=0))
## tnl[["tDNA"]] <- data.frame(matrix(ncol=ncol(tnl_raw[["tDNA"]]),nrow=0))
## for (currID in unique(tnl_raw[["gene"]]$geneID)){
##     temp_gene <- filter(tnl_raw[["gene"]],geneID==currID)
##     temp_exon <- filter(tnl_raw[["exon"]],geneID==currID)
##     temp_cds <- filter(tnl_raw[["cds"]],geneID==currID)
##     temp_tDNA <- filter(tnl_raw[["tDNA"]],geneID==currID)
##     strand <- temp_gene$strand[[1]]
##     if (strand == "+") {
##         start <- min(temp_cds$chromStart)
##         temp_gene <- mutate(temp_gene,chromStart=chromStart-start,
##                             chromEnd=chromEnd-start)
##         temp_exon <- mutate(temp_exon,chromStart=chromStart-start,
##                             chromEnd=chromEnd-start)
##         temp_cds <- mutate(temp_cds,chromStart=chromStart-start,
##                            chromEnd=chromEnd-start)
##         temp_tDNA <- mutate(temp_tDNA,chromStart=tDNAStart-start,
##                             chromEnd=tDNAEnd-start)
##     } else {
##         start <- max(temp_cds$chromEnd)
##         cstart <- temp_gene$chromStart
##         cend <- temp_gene$chromEnd
##         temp_gene <- mutate(temp_gene,chromStart=-(cend-start),
##                             chromEnd=-(cstart-start),
##                             strand=rep("+",nrow(temp_gene)))
##         cstart <- temp_exon$chromStart
##         cend <- temp_exon$chromEnd
##         temp_exon <- mutate(temp_exon,chromStart=-(cend-start),
##                             chromEnd=-(cstart-start),
##                             strand=rep("+",nrow(temp_exon)))
##         cstart <- temp_cds$chromStart
##         cend <- temp_cds$chromEnd
##         temp_cds <- mutate(temp_cds,chromStart=-(cend-start),
##                             chromEnd=-(cstart-start),
##                            strand=rep("+",nrow(temp_cds)))
##         temp_tDNA <- mutate(temp_tDNA,chromStart=-c(tDNAEnd-start),
##                             chromEnd=-c(tDNAStart-start),
##                             strand=rep("+",nrow(temp_tDNA)))
##     }
##     tnl[["gene"]] <- rbind(tnl[["gene"]],temp_gene)
##     tnl[["exon"]] <- rbind(tnl[["exon"]],temp_exon)
##     tnl[["cds"]] <- rbind(tnl[["cds"]],temp_cds)
##     tnl[["tDNA"]] <- rbind(tnl[["tDNA"]],temp_tDNA)
## }

## subtract_min <- function(df_ls){
##     min_coord <- min(c(df_ls[["gene"]]$chromStart,df_ls[["exon"]]$chromStart,
##                        df_ls[["cds"]]$chromStart,df_ls[["tDNA"]]$chromStart))
##     out <- list()
##     out[["gene"]] <- mutate(df_ls[["gene"]],chromStart=chromStart-min_coord,
##                             chromEnd=chromEnd-min_coord)
##     out[["exon"]] <- mutate(df_ls[["exon"]],chromStart=chromStart-min_coord,
##                             chromEnd=chromEnd-min_coord)
##     out[["cds"]] <- mutate(df_ls[["cds"]],chromStart=chromStart-min_coord,
##                            chromEnd=chromEnd-min_coord)
##     out[["tDNA"]] <- mutate(df_ls[["tDNA"]],chromStart=chromStart-min_coord,
##                             chromEnd=chromEnd-min_coord)
##     return(out)
## }
## tnl2 <- subtract_min(tnl)
## tnl_5 <- as.character(unique(tnl_raw[["gene"]]$geneID)[1:5])
## tnl_all <- as.character(unique(tnl_raw[["gene"]]$geneID))
## tnl_sub <- lapply(tnl2,function(x) {filter(x,geneID %in% tnl_5)})
## tnl_sub <- tnl2

## gtrack <- GenomeAxisTrack()

## tracks <- list()
## for (tnl in tnl_all){
##     gene_cds <- filter(tnl_sub[["cds"]],geneID==tnl)
##     gene_tDNA <- filter(tnl_sub[["tDNA"]],geneID==tnl)
##     track1 <- AnnotationTrack(start=gene_cds$chromStart,end=gene_cds$chromEnd,
##                               chromosome="chr0",group=gene_cds$geneID,
##                               strand=gene_cds$strand,name=tnl,shape="smallArrow",
##                               fill="grey",col="dark grey",showId=T)
##     if (nrow(gene_tDNA) > 0){
##         track2 <- AnnotationTrack(start=gene_tDNA$chromStart,end=gene_tDNA$chromEnd,
##                                   chromosome="chr0",name=tnl,group=gene_tDNA$tDNA,shape="box",
##                                   col="red",fill="red",showId=F,stacking="dense",alpha=0.8)
##         tracks[[tnl]] <- OverlayTrack(trackList = list(track1,track2),name=tnl)
##     } else {
##         tracks[[tnl]] <- track1
##     }
## }
## plotTracks(append(gtrack,tracks),groupAnnotation="group",name="TNL")


## ## extract NLRs in DM8 cluster
## dm8 <- list()
## dm8[["gene"]] <- filter(nlrs_gene,9488400<=chromStart & chromStart<=9565600 & chrom=="chr4")
## dm8[["exon"]] <- inner_join(filter(nlrs_exon,9488400<=chromStart & chromStart<=9565600 & chrom=="chr4"),
##                             select(dm8[["gene"]],geneID,strand),by="geneID")
## dm8[["cds"]] <- inner_join(filter(nlrs_cds,9488400<=chromStart & chromStart<=9565600 & chrom=="chr4"),
##                            select(dm8[["gene"]],geneID,strand),by="geneID")
## dm8[["tDNA"]] <- filter(tdna_dat,9488400<=chromStart & chromStart<=9565600 & chrom=="chr4")

## ## figure out arrows
## format_arrows <- function(df){
##     df <- arrange(mutate(df,strand=as.character(strand)),chromStart)
##     for (row in seq(nrow(df))){
##         if (((df[row,"strand"] == "+") &
##              (row == nrow(df) | !identical(df[row+1,"geneID"],df[row,"geneID"]))) |
##             ((df[row,]$strand == "-") &
##              (row == 1 | !identical(df[row-1,"geneID"],df[row,"geneID"])))){
##             next
##         } else {
##             df[row,"strand"] <- "*"
##         }
##     }
##     return(df)
## }

## gtrack <- GenomeAxisTrack()

## gene <- data.frame(distinct(dm8[["gene"]][,c(1:3,5,6)]))
## gene <- arrange(mutate(gene,strand=as.character(strand)),chromStart)
## track0 <- AnnotationTrack(start=gene$chromStart,end=gene$chromEnd,
##                           chromosome=gene$chrom,strand=gene$strand,
##                           group=gene$geneID,name="DM8: NLR Genes",
##                           fill="grey",col="dark grey",showId=T,cex.group=1)

## gene_ex <- data.frame(distinct(dm8[["exon"]][,c(1:3,5,10)]))
## gene_ex <- format_arrows(gene_ex)
## track1 <- AnnotationTrack(start=gene_ex$chromStart,end=gene_ex$chromEnd,
##                           chromosome=gene_ex$chrom,strand=gene_ex$strand,
##                           group=gene_ex$geneID,name="Exon",
##                           fill="dark green",col="dark green",showId=F)

## gene_cds <- data.frame(distinct(dm8[["cds"]][,c(1:3,5,10)]))
## gene_cds <- format_arrows(gene_cds)
## track2 <- AnnotationTrack(start=gene_cds$chromStart,end=gene_cds$chromEnd,
##                          chromosome=gene_cds$chrom,group=gene_cds$geneID,
##                          strand=gene_cds$strand,name="CDS",
##                          fill="orange",col="dark orange",showId=F)
## plotTracks(list(gtrack,track0,track1,track2),groupAnnotation="group")

## gene_tDNA <- distinct(select(filter(dm8[["gene"]],tDNA!="."),chrom,tDNAStart,tDNAEnd,tDNA))
## track3 <- AnnotationTrack(start=dm8[["tDNA"]]$chromStart,end=dm8[["tDNA"]]$chromEnd,
##                           chromosome=dm8[["tDNA"]]$chrom,name="T-DNA",group=dm8[["tDNA"]]$tDNA,
##                           col="red",fill="red",showId=F,stacking="dense")
## ## ot <- OverlayTrack(trackList = list(track2,track3))
## plotTracks(list(gtrack,track0,track1,track2,track3),groupAnnotation="group",fontsize=15)

library(GenomicRanges)
gene_tDNA <- select(filter(dm8[["gene"]],tDNA!="."),chrom,tDNAStart,tDNAEnd)
gene_tDNA <- mutate(gene_tDNA,yaxis=rep(0,nrow(gene_tDNA)))
colnames(gene_tDNA) <- c("chromosome","start","end","T-DNA")
gene_tDNA <- makeGRangesFromDataFrame(gene_tDNA,start.field="tDNAStart",end.field="tDNAEnd",
                                      keep.extra.columns=T)
track_points <- DataTrack(data=gene_tDNA,name="T-DNA",genome="TAIR10")
plot(track_points)
ot <- OverlayTrack(trackList = list(track,track_points))
plotTracks(ot,groupAnnotation="group",stackHeight=0.3)
