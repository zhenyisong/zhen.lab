# Title     : the R script to atriumization on heart regeneration
# Objective : to find the initiation factor leading to the start of heart regeneration
# Created by: Yisong Zhen
# Created on: 11/6/2019
#----
pkgs <- c( 'tidyverse','Rsubread','org.Mm.eg.db','edgeR',
           'limma', 'DESeq2','DESeq','genefilter','grid', 
           'openxlsx','pheatmap','gridExtra','ggrepel',
           'annotate','clusterProfiler', 'cowplot',
           'RColorBrewer','affy','ComplexHeatmap',
           'cluster','factoextra','RDAVIDWebService',
           'Rsamtools', 'devtools','mouse4302.db',
           'GenomicAlignments', 'magrittr', 'network',
           'reshape2', 'stringr', 'rtracklayer','mgu74a.db',
           'TxDb.Mmusculus.UCSC.mm10.knownGene', 'GGally',
           'Mus.musculus', 'org.Mm.eg.db', 'FactoMineR',
           'BSgenome.Mmusculus.UCSC.mm10')



load.lib           <- lapply(pkgs, require, character.only = TRUE)


#---
# define data and file deposition layout
# this maybe the constant
# absolute path
# <annotation>
#    ----genome : the raw genome fasta seqeunces
#    ----index  : rsubread indexed genome data
#    ----annot  : all the genome annotation data
# <project>
#    ----<project name>
#            ----rawdata  : deposition SRA data or/and fastq.gz data
#            ----rdata    : save the processed .Rdata from the R procedure
#            ----rsubread : save the bam aligned RNA sequencing data
#---



data.root.dir          <- '/home/zhenyisong/cloudatlas/projects'
annotation.root.dir    <- '/home/zhenyisong/cloudatlas/annotation'
genome.annotation.dir  <- paste0(annotation.root.dir,'/annot') %>% file.path()
rsubread.index.dir     <- paste0(annotation.root.dir, '/index') %>% file.path()

genome.fasta.file      <- paste0(annotation.root.dir,
                                '/genome/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.fa')
raw.data.dir           <- paste0(data.root.dir, '/GSE64403')
#change the file tag index and gtf files downloading
rsubread.index.file    <- paste0( rsubread.index.dir, '/mm10.00.b.tab')
mRNA.GRCm38.GTF.file   <- paste0( genome.annotation.dir,'/igenomes.mm10.ensembl.gtf') %>%file.path()


# data output
onedrive.path          <- paste0('C:\\Users\\zheny\\OneDrive\\fuwaihospital\\',
                                 'publication\\2020papers\\atriumization')


# microarray data path
#

setwd('C:\\Users\\zheny\\Documents\\Rdata_projects')
load('GSE64403.Rdata')
GSE1479.raw.path     <- 'C:\\Users\\zheny\\Documents\\Rdata_projects\\GSE1479'
GSE775.raw.path      <- 'C:\\Users\\zheny\\Documents\\Rdata_projects\\GSE775'

# expand the data deposition layout
# define the dir structure or file names
#----

create.a.project.layout <- function(raw.path) {
    project.raw.data.dir    <- paste0(raw.data.dir, '/data') %>% file.path()
    project.output.dir      <- paste0(raw.data.dir,'/rsubread') %>% file.path()
    project.Rdata.dir       <- paste0(raw.data.dir, '/rdata') %>% file.path()
    if (! dir.exists(project.Rdata.dir)) dir.create( project.Rdata.dir,
                                                     showWarnings = FALSE, recursive = TRUE)
    project.layout <- list(raw.data.dir = project.raw.data.dir,
                           output.dir   = project.output.dir,
                           rdata.dir    = project.Rdata.dir)
    return(project.layout)
}

generate.output.filename    <- function(double.ends.rawdata.dir, output.dir) {
    project.files       <- double.ends.rawdata.dir %>%
                           list.files( pattern   = 'SRR.*.fastq.gz$',
                                       all.files = FALSE, full.names  = TRUE,
                                       recursive = FALSE, ignore.case = FALSE, include.dirs = F)
    read.1.files       <- project.files[grep('_1',project.files)]
    read.2.files       <- project.files[grep('_2',project.files)]
    project.output.filenames  <- basename(read.1.files) %>%
                                 sub(pattern = '_1.fastq.gz', replacement = '') %>%
                                 paste0(output.dir,'/', . ,'.bam')
    if (! dir.exists(output.dir)) dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)
    files.layout <- list(output = project.output.filenames, readfile1 = read.1.files, readfile2 = read.2.files )
    return(files.layout)
}

create.rsubread.index <- function(subread.index.dir, subread.file, base.string = 'mm10') {
    if (! dir.exists(subread.index.dir)) dir.create(subread.index.dir,
                                                    showWarnings = FALSE, recursive = TRUE)
    setwd(subread.index.dir)
    if (! file.exists(subread.file)) buildindex( basename   = base.string,
                                                 indexSplit = TRUE,
                                                 memory     = 4000,
                                                 reference  = genome.fasta.file )
    return(base.string)
}


'
project.layout      <- create.a.project.layout(raw.data.dir)
filenames.layout    <- generate.output.filename(project.layout$raw,project.layout$output)
rsubread.basename   <- create.rsubread.index(rsubread.index.dir, rsubread.index.file)
'

align.raw.data      <- function(raw.data.layout, index.dir, index.file) {
    project.layout      <- create.a.project.layout(raw.data.layout)
    filenames.layout    <- generate.output.filename(project.layout$raw,project.layout$output)
    rsubread.basename   <- create.rsubread.index(index.dir, index.file)
    align( index          = rsubread.basename,
           readfile1      = filenames.layout$readfile1,
           readfile2      = filenames.layout$readfile2,
           input_format   = 'gzFASTQ',
           type           = 'rna',
           output_file    = filenames.layout$output,
           output_format  = 'BAM',
           PE_orientation = 'fr',
           nthreads       = 2,
           indels         = 1,
           maxMismatches  = 3,
           phredOffset    = 33,
           unique         = T )
}

get.featureCount.result <- function(raw.data.layout,GTF.file) {
    project.layout          <- create.a.project.layout(raw.data.layout)
    filenames.layout        <- generate.output.filename(project.layout$raw,project.layout$output)
    featureCount.result     <- featureCounts( filenames.layout$output, useMetaFeatures = TRUE,
                                              countMultiMappingReads = FALSE,
                                              strandSpecific         = 0,
                                              isPairedEnd            = TRUE,
                                              requireBothEndsMapped  = TRUE,
                                              autosort               = TRUE,
                                              nthreads               = 2,
                                              annot.ext              = GTF.file,
                                              isGTFAnnotationFile    = TRUE,
                                              GTF.featureType        = 'exon',
                                              GTF.attrType           = 'gene_id',
                                              allowMultiOverlap      = TRUE)
    return(featureCount.result)
    
}


#---
# Figure 1
#---

#----
# get the atrium specific genes from microarray data
# original data processing code is in the script 'GSE1479_Microarray.R'
# at the fuwai dir published_paper?
#----

get.microarray.exprs <- function(data.path) {
    exprs.data  <- data.path %>% 
                   ReadAffy(celfile.path = .) %>% 
                   rma() %>% exprs()
    return(exprs.data)
}

secure.GSE1479.PCA <- function(GSE1479.exprs) {
    complete.data.pca  <- prcomp(t(GSE1479.exprs))
    return(complete.data.pca)
}

get.PCA.plot <- function(GSE1479.exprs) {
    pca.result        <-  secure.GSE1479.PCA(GSE1479.exprs)
    whole.av.attr     <- factor( c( rep(1,6),rep(2:3,5, each = 3)), 
                                 levels = 1:3, labels = c('heart','ventricle','atrium'))
    GSE1479.colnames  <- c('E10.5_1','E10.5_2','E10.5_3','E11.5_1','E11.5_2',
                                  'E11.5_3','E12.5_1','E12.5_2','E12.5_3','E12.5_1A',
                                  'E12.5_2A','E12.5_3A','E13.5_1','E13.5_2','E13.5_3',
                                  'E13.5_1A','E13.5_2A','E13.5_3A','E14.5_1','E14.5_2',
                                  'E14.5_3','E14.5_1A','E14.5_2A','E14.5_3A','E16.5_1',
                                  'E16.5_2','E16.5_3','E16.5_1A','E16.5_2A','E16.5_3A',
                                  'E18.5_1','E18.5_2','E18.5_3','E18.5_1A','E18.5_2A','E18.5_1A')
    GSE1479_pca_ggplot <- data.frame( PC1 = pca.result$x[,1], 
                                      PC2 = pca.result$x[,2], 
                                      AV  = whole.av.attr) %>% 
                          ggplot() + 
                          geom_point( aes(x = PC1, y = PC2, color = AV), size = 2.5) + 
                          geom_text_repel( aes(x = PC1, y = PC2, 
                                           label = GSE1479.colnames),
                                           size = 2.0) +
                          scale_colour_manual( name   = 'tissue types',
                                               values = c( '#F5D547','#DB3069','#1446A0'),
                                               labels = c( 'heart','ventricle','atrium')) +   
                          theme_classic() +
                          theme( aspect.ratio    = 1, 
                                 legend.position = c(0.9, 0.5),
                                 legend.title    = element_text(size = 12),
                                 legend.text     = element_text(size = 10.5))
    return(GSE1479_pca_ggplot)
}

'
get.PCA.plot(GSE1479.exprs)

'

get.av.genes.from.PCA <- function(pca.result) {
    PC1 <- pca.result$rotation[,'PC1'] %>% unlist
    PC2 <- pca.result$rotation[,'PC2'] %>% unlist
    gene.symbols <- rownames(pca.result$rotation) %>%
                    mget(mouse4302SYMBOL,ifnotfound = NA) %>%
                    make.names(unique = T)%>% unlist()
    
    names(PC1) <- gene.symbols 
    names(PC2) <- gene.symbols 
    PC1 <- sort(PC1,decreasing = F)
    PC2 <- sort(PC2,decreasing = F)
    
}

get.av.genes.from.DGE <- function(GSE1479.exprs) {

    chamber.norm.exprs <- GSE1479.exprs[,7:36]
    av.fac <- rep(rep(1:0,each = 3),5) %>% 
              factor(levels = 0:1, labels = c('atrium','ventricle'))
    atrium.median      <- apply(chamber.norm.exprs[,av.fac == 'atrium'], 1, median)
    ventricle.median   <- apply(chamber.norm.exprs[,av.fac == 'ventricle'], 1, median)
    chamber.contrast   <- ventricle.median - atrium.median
    chamber.exprs      <- chamber.norm.exprs[chamber.contrast >= 1 | chamber.contrast <= -1,]
    
    pt   <- apply(chamber.exprs, 1, function(x) t.test(x ~ av.fac)$p.value);
    pw   <- apply(chamber.exprs, 1, function(x) wilcox.test(x ~ av.fac)$p.value)
    selected.chamber.genes <- chamber.exprs[pt < 0.001 & pw < 0.001,]
    genes.name   <- rownames(selected.chamber.genes) %>%
                    getSYMBOL('mouse4302')
    result <- hclust( as.dist( 1 - cor(t( selected.chamber.genes ))),
                      method = 'complete') %>%
              cutree(k = 2, h = 5)
    atrium.probe.ID       <- names(result)[result == 2]
    ventricle.probe.ID    <- names(result)[result == 1]
    atrium.gene.names     <- mget(atrium.probe.ID,mouse4302SYMBOL,ifnotfound = NA) %>%
                             unique() %>% na.omit() %>% unlist()
    atrium.EntrezID       <- mapIds(mouse4302.db, 
                                    keys = atrium.probe.ID, 
                                    c('ENTREZID'), keytype = 'PROBEID') %>%
                             unique() %>% na.omit() %>% unlist
    atrium.gene.names     <- atrium.gene.names[!is.na(atrium.gene.names)]
    ventricle.gene.names  <- mget(ventricle.probe.ID,mouse4302SYMBOL,ifnotfound = NA) %>%
                             unique() %>% na.omit() %>% unlist()
    ventricle.EntrezID    <- mapIds(mouse4302.db, 
                                    keys = ventricle.probe.ID, 
                                    c('ENTREZID'), keytype = 'PROBEID') %>%
                             unique() %>% na.omit() %>% unlist
    # remove NA
    ventricle.gene.names  <- ventricle.gene.names[!is.na(ventricle.gene.names)]
    av.groups             <- list( atrium       = atrium.gene.names , 
                                   ventricle    = ventricle.gene.names,
                                   atrium.ID    = atrium.EntrezID,
                                   ventricle.ID = ventricle.EntrezID,
                                   atrium.probeID    = atrium.probe.ID,
                                   ventricle.probeID = ventricle.probe.ID)
    return(av.groups)
}


get.chamber.scree.plot <- function(GSE1479.exprs) {
    chamber.list           <- get.av.genes.from.DGE(GSE1479.exprs)
    probe.ID               <- rownames(GSE1479.exprs)
    chamber.probeID        <- c(chamber.list$atrium.probeID, chamber.list$ventricle.probeID)
    chamber.genes.exprs    <- GSE1479.exprs[probe.ID %in% chamber.probeID, ]
    GSE1479_pca            <- PCA(t(chamber.genes.exprs), scale.unit = FALSE, graph = FALSE)
    GSE1479_scree_ggplot   <- fviz_eig(GSE1479_pca, addlabels = TRUE, main = '') + 
                              theme(aspect.ratio = 1) + theme_minimal()
    return(GSE1479_scree_ggplot)
    
}

#---
# https://david.ncifcrf.gov/content.jsp?file=DAVID_API.html#approved_list
# https://stackoverflow.com/questions/31480579/r-david-webservice-sudden-transport-error-301-error-moved-permanently
# https://github.com/YuLab-SMU/clusterProfiler/issues/35
#---

get.atrium.DAVID.plot <- function(GSE1479.exprs) {
    av.genes      <- get.av.genes.from.DGE(GSE1479.exprs)
    av.ids        <- paste(av.genes$atrium.ID, collapse = ',')
    david.service <- DAVIDWebService$new(email = 'zhenyisong@fuwaihospital.org',
                                         url   = 'https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/')
    david.param   <- addList(david.service, av.ids, 
                             idType   = 'ENTREZ_GENE_ID',
                             listName = 'atrium_gene_set', 
                             listType = 'Gene')
    #setCurrentSpecies(david.service,10090)
    setAnnotationCategories(david.service, c('UP_KEYWORDS'))
    result.david <- getClusterReport(david.service)
    data.final   <- members(result.david)[[1]] %>% as.data.frame %$%
                    {data.frame( termlabels  = Term,
                     LogPvalue = - log(Bonferroni) ) } %>%
                    ggplot( aes( x = reorder(termlabels, LogPvalue), 
                                 y = LogPvalue)) + 
                    geom_bar( stat = 'identity', width = 0.6, 
                              position = position_dodge(width = 0.05), size = 15) +
                    theme(axis.text.x = element_text(angle = 60,hjust = 1, size = 6),
                          axis.text.y = element_text(hjust = 1, size = 10),
                          axis.title  = element_text(size = 12)) +
                    ylab('-log(p.adjust)') + 
                    xlab('Enrichment of SP Terms') + 
                    coord_flip() +
                    theme(aspect.ratio = 1)
    
    
    return(data.final)
    
}

"
get.atrium.DAVID.plot(GSE1479.exprs)
"

get.atrium.GO.BF.plot  <- function(GSE1479.exprs) {
    av.genes     <- get.av.genes.from.DGE(GSE1479.exprs)
    go.datafarme <- enrichGO( gene  = av.genes$atrium.ID,
                              OrgDb =  org.Mm.eg.db,
                              ont   = 'BP',
                              pAdjustMethod = 'BH',
                              qvalueCutoff  = 0.05)
    atrium.GO.plot <- go.datafarme %>% as.data.frame() %$%
                      {data.frame( Golabels  = paste(ID, Description, sep = ' '),
                                   LogPvalue = - log(p.adjust) ) } %>% .[1:9,] %>%
                      ggplot( aes( x = reorder(Golabels, LogPvalue), 
                                   y = LogPvalue)) + 
                      geom_bar( stat = 'identity', width = 0.7, 
                                position = position_dodge(width = 0.1), size = 20) +
                      theme(axis.text.x = element_text(angle = 60,hjust = 1, size = 7),
                            axis.text.y = element_text(hjust = 1, size = 7),
                            axis.title  =  element_text(size = 11) ) +
                      ylab('-log(p.adjust)') + 
                      xlab('Enrichment of GO BP') + 
                      theme(aspect.ratio = 1)
    return(atrium.GO.plot)
}

get.figure.1 <- function(GSE1479.exprs) {
    figure.1.A <- get.PCA.plot(GSE1479.exprs)
    figure.1.B <- get.chamber.scree.plot(GSE1479.exprs)
    figure.1.C <- get.atrium.GO.BF.plot(GSE1479.exprs)
    figure.1.D <- get.atrium.DAVID.plot(GSE1479.exprs)
    figure.1.final <- ggdraw() +
                      draw_plot( figure.1.A, 
                                 x =.00, y = .5,  width = .5, height = .5) +
                      draw_plot( figure.1.B,   
                                x =  0.5, y =  0.5, width = .5, height = .5) +
                      draw_plot( figure.1.C, 
                                x = 0.0, y = 0.0,  width = .5, height = .5) +
                      draw_plot( figure.1.D,  
                                 x = .5, y =  0,  width = .5, height = .5) +
                      draw_plot_label( label = LETTERS[1:4],
                                       x     = c(0,.5,0,.5),
                                       y     = c(1,1,.5,.5),
                                       size  = 10)
    return(figure.1.final)
}




#---
# Figure 2.
#---

.get.GSE775.groups <- function() {
    groups <- list( h.1  = c(1:3,24:26,42:44),
                    h.4  = c(13:15,36:38,54:56),
                    h.24 = c(7:9,30:32,48:50),
                    h.48 = c(10:12,33:35,51:53),
                    w.1  = c(4:6,27:29,45:47),
                    w.8  = c(16:18,39:41,57:59) )
    return(groups)
}


secure.tidy.data <- function(GSE775.exprs, gene.name = 'Sln') {
    infarction.norm.exprs <- GSE775.exprs
    groups                <- .get.GSE775.groups()
    gene.probe.id         <- mget(gene.name,revmap(mgu74aSYMBOL)) %>% unlist() %>% .[1]
    gene.norm.exprs       <- infarction.norm.exprs[gene.probe.id, groups %$% 
                                                   c(h.1,h.4,
                                                     h.24,h.48,w.1,w.8)]
    gene.tidy.data  <- data.frame( time  = factor(rep(1:6,each = 9),
                                                  levels = 1:6,
                                                  labels = c('h.1','h.4','h.24','h.48','w.1','w.8')),
                                  exprs  = infarction.norm.exprs[ gene.probe.id,
                                                                  groups %$% c(h.1,h.4,h.24,h.48,w.1,w.8)],
                                  region = factor( rep(rep(1:3,each = 3),6),
                                                   levels = 1:3,
                                                   labels = c('lv','ilv','nilv')))
    return(gene.tidy.data)
}
    
# TODO
# the figure show two legends. why?
# https://stackoverflow.com/questions/35108443/ggplot2-make-legend-key-fill-transparent
#---
get.nppa.regression.plot <- function(data.exprs, gene.name = 'Nppa') {
    gene.tidy.data <- secure.tidy.data(data.exprs,gene.name)
    gene.figure    <- ggplot(data = gene.tidy.data, aes(x = time, y = exprs, group = region) ) +
                      geom_point(aes(color = region, shape = region), show.legend = F) + 
                      geom_smooth(aes(color = region), method = lm, se = T, 
                                  show.legend = F, alpha = 0.1) + 
                      geom_smooth(aes(color = region), method = lm, se = F, alpha = 0.1) +
                      theme(legend.position = c(0.2, 0.8)) +
                      theme(axis.title.x = element_blank()) +
                      scale_color_manual( values = c('#540D6E','#EE4266','#FFD23F'), 
                                          name   = 'regions') +
                      labs(y = 'microarray exprs value') + 
                      theme(legend.title.align = 0.5, legend.title = element_text(size = 13),
                            legend.key = element_rect(colour = 'transparent', fill = 'transparent')) +
                      theme( aspect.ratio = 1)
    return( gene.figure )
}

"
get.nppa.regression.plot(data.affy)
"

get.acta2.boxplot <- function(data.exprs, gene.name = 'Acta2') {
    gene.tidy.data <- secure.tidy.data(data.exprs,gene.name)
    gene.figure    <- ggplot(data = gene.tidy.data, aes( x     = time, 
                                                         y     = exprs) ) +
                      geom_boxplot(aes(fill = region),
                                       position = position_dodge((0.5))) +
                      stat_summary(fun.y = mean, geom = 'point',
                                   shape = 18, size = 2, color = 'grey') +
                      scale_fill_manual( values = c('#540D6E','#EE4266','#FFD23F'), 
                                         name   = 'regions') +
                      theme(legend.position = c(0.1, 0.8), 
                            legend.title.align = 0.5,
                            legend.title = element_text(size = 13)) +
                      theme( aspect.ratio = 1)
    return(gene.figure)
}

"
get.acta2.boxplot(GSE775.exprs)
"

# Pheatmap error: 'gpar' element 'fill' must not be length 0 , 
# what to do about that in an assymmetrical matrix
# https://www.biostars.org/p/240404/

get.atrium.marker.heatmap <- function( data.affy) {
    atrium.markers     <- c('Sln','Fgf12','Myl7','Nppa')
    groups             <- .get.GSE775.groups() %$% c(h.1,h.4,h.24,h.48,w.1,w.8)
    gene.probe.ids     <- mget(atrium.markers,revmap(mgu74aSYMBOL)) %>% unlist()
    gene.norm.exprs    <- data.affy[gene.probe.ids, groups]
    column.annotation  <- factor( rep(rep(1:3,each = 3),6),
                                levels = 1:3,
                                labels = c('lv','ilv','nilv')) %>%
                         data.frame(regions = .) %>%
                         data.frame(time  = factor(rep(1:6,each = 9),
                                                   levels = 1:6,
                                                   labels = c('h.1','h.4','h.24','h.48','w.1','w.8')),.) 

    rownames(column.annotation)  <- colnames(gene.norm.exprs)
    column.colors                 <- list(time = c(h.1 = '#320D6D', h.4 = '#FFBFB7',h.24 = 'white',
                                                   h.48= '#FFD447', w.1 = '#700353',w.8 = '#4C1C00'),
                                          regions = c(lv = '#540D6E',ilv = '#EE4266',nilv = '#FFD23F'))
    
    heat.map <- pheatmap(gene.norm.exprs, scale = 'none', cluster_cols = F, cluster_rows = F,
                         color = colorRampPalette(c('#0B3954', '#BFD7EA', '#FF6663'))(50),
                         cellheight = 23, cellwidth = 6, border_color = NA,
                         labels_row     = atrium.markers, 
                         show_colnames  = F, annotation_col = column.annotation,
                         annotation_colors = column.colors, silent = T)
    return(heat.map)
}



atriumization.GSEA.analysis <- function( GSE775.exprs, GSE1479.exprs, 
                                         group = 'w.1', contrast = 'ilv - lv', chamber = 'atrium') {
    av.genes    <- get.av.genes.from.DGE(GSE1479.exprs)
    chamber.set <- av.genes$atrium
    if (chamber == 'ventricle') chamber.set <- av.genes$ventricle
   
    limma.result     <- get.GSE775.limma.result(GSE775.exprs,group, contrast) 
    chamber.index    <- rownames(limma.result) %>% getSYMBOL('mgu74a') %in% chamber.set
    chamber.probe.id <- rownames(limma.result)[chamber.index]
    t.statistic      <- limma.result$t[chamber.index] %>% 
                        sum()/sum(chamber.index) * sqrt(sum(chamber.index))
    logFC.data            <- limma.result$logFC
    names(logFC.data)     <- rownames(limma.result)
    logFC.data.sorted     <- sort(logFC.data, decreasing = TRUE)
    logFC.data.sorted.id  <- seq(1:length(logFC.data.sorted))
    chamber.pseudo.binary <- names(logFC.data.sorted) %>% 
                             getSYMBOL('mgu74a') %in% chamber.set
    chamber.pseudo.id        <- logFC.data.sorted.id[chamber.pseudo.binary]
    names(logFC.data.sorted) <- logFC.data.sorted.id
    
    chamber.gsea.set         <- data.frame( diseaseId = rep(chamber,length(chamber.pseudo.id)), 
                                            geneId    = chamber.pseudo.id, check.names = TRUE)
    chamber.GSEA             <- GSEA(logFC.data.sorted, 
                                     TERM2GENE = chamber.gsea.set,
                                     minGSSize = 20, pvalueCutoff = 1)
    #(1 - pnorm(t.statistic))
    atrium.GSEA.plot         <- gseaplot(chamber.GSEA, geneSetID = chamber,
                                         by = 'runningScore', pvalue_table = T,
                                         color.line = '#20BF55', color.vline = '#0B4F6C' ) +
                                theme(axis.title = element_text(size = 12),
                                      axis.text  = element_text(size = 8)) + 
                                labs(caption = '(data source: GSE775)') +
                                theme( aspect.ratio = 1)
    result <- list( class.t      =  1 - pnorm(t.statistic), 
                    chamber.gsea = chamber.GSEA, 
                    gsea.plot    = atrium.GSEA.plot)
    return(result)
}

"

gsea.result <- atriumization.GSEA.analysis(GSE775.exprs,GSE1479.exprs)

"

get.figure.2 <- function(GSE775.exprs,GSE1479.exprs) {
    figure.2.A <- get.nppa.regression.plot(GSE775.exprs,'Nppa')
    figure.2.B <- get.acta2.boxplot(GSE775.exprs)
    figure.2.C <- get.atrium.marker.heatmap(GSE775.exprs)$gtable
    figure.2.D <- atriumization.GSEA.analysis(GSE775.exprs,GSE1479.exprs) %$% 
                  gsea.plot
    figure.2.final <- ggdraw() +
                      draw_plot( figure.2.A, 
                                 x =.00, y = .5,  width = .5, height = .5) +
                      draw_plot( figure.2.B,   
                                 x =  0.5, y =  0.5, width = .5, height = .5) +
                      draw_plot( figure.2.C, 
                                 x = 0.0, y = 0.0,  width = .5, height = .5) +
                      draw_plot( figure.2.D,  
                                 x = .5, y =  0,  width = .5, height = .5) +
                      draw_plot_label( label = LETTERS[1:4],
                                 x     = c(0,.5,0,.5),
                                 y     = c(1,1,.5,.5),
                                 size  = 10)
    return(figure.2.final)
    
}


#---
# Figure 3
#---


#---
# SLN expression pattern for 4
# test hypothesis, the anterior-poserior boundary will disappeear
# ventricle
# see experimental design
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1570682
# Boyer data
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64403
# iP0  : 0 day postnatal myocardium, isolated CMs
# iDS7 : 7 days post-sham surgery myocardium, isolated CMs
# iDR7 : 7 days post-resection myocardium, isolated CMs
# vP0  : 0 day postnatal ventricular myocardium
# ex0hr: Explanted, dissociated adult cardiomyocytes, 0hr
# vAd  : Adult ventricular myocardium
# vD7R : 7 days post-resection surgery ventricular myocardium
# vD7S : 7 days post-sham surgery ventricular myocardium
#---

.get.rsubread.exprs <- function(rsubread.result) {
    rsubread.gene.symbols    <- mapIds( org.Mm.eg.db, keys = rsubread.result$annotation$GeneID,
                                        keytype = 'ENSEMBL', column = 'SYMBOL') %>% 
                                make.names(unique = T)
    rsubread.exprs           <- rsubread.result$counts
    rownames(rsubread.exprs) <- rsubread.gene.symbols
    return(rsubread.exprs)
}

.get.isolated.cell.deseq2 <- function(rsubread = GSE64403.results) {
    GSE64403.isolated.cells <- .get.rsubread.exprs(rsubread)[,1:12]
    GSE64403.cardiac.cells.condition   <- factor( rep(1:4,each = 3),levels = 1:4,
                                                  labels = c('iP0','iP4','iD7S','iD7R'))
    GSE64403.atrium.pValue             <- newCountDataSet( GSE64403.isolated.cells, 
                                                           GSE64403.cardiac.cells.condition ) %>%
        estimateSizeFactors() %>% estimateDispersions() %>%
        nbinomTest('iD7S', 'iD7R')
    rsubread.gene.symbols              <- mapIds( org.Mm.eg.db, keys = rsubread$annotation$GeneID,
                                                  keytype = 'ENSEMBL', column = 'SYMBOL') %>% 
        make.names(unique = T)
    rownames(GSE64403.atrium.pValue)   <- rsubread.gene.symbols
    #atrium.gene.P.value                 <- GSE64403.atrium.pValue[GSE64403.atrium.pValue$id == 'Nppa',]
    return(GSE64403.atrium.pValue)
}

.get.heart.tissue.deseq2    <- function(rsubread = GSE64403.results) {
    GSE64403.heart.tissue              <- .get.rsubread.exprs(rsubread)[,29:36]
    GSE64403.cardiac.tissue.condition  <- factor( rep(1:4,each = 2),levels = 1:4,
                                                  labels = c('vD1S','vD1R','vD7R','vD7S'))
    GSE64403.tissue.pValue             <- newCountDataSet( GSE64403.heart.tissue, 
                                                           GSE64403.cardiac.tissue.condition ) %>%
                                          estimateSizeFactors() %>% estimateDispersions() %>%
                                          nbinomTest('vD7S', 'vD7R')
    rsubread.gene.symbols              <- mapIds( org.Mm.eg.db, keys = rsubread$annotation$GeneID,
                                                  keytype = 'ENSEMBL', column = 'SYMBOL') %>% 
                                          make.names(unique = T)
    rownames(GSE64403.tissue.pValue)   <- rsubread.gene.symbols
    #atrium.gene.P.value                 <- GSE64403.tissue.pValue[GSE64403.tissue.pValue$id == 'Sln',]
    return(GSE64403.tissue.pValue)
}

secure.GSE64403.GSEA <- function(GSE1479.exprs, chamber = 'atrium') {
    tissue.deseq2        <- .get.heart.tissue.deseq2()
    av.genes             <- get.av.genes.from.DGE(GSE1479.exprs)
    chamber.set          <- av.genes$atrium
    if (chamber == 'ventricle') chamber.set <- av.genes$ventricle
    logFC.data               <- tissue.deseq2$log2FoldChange
    names(logFC.data)        <- rownames(tissue.deseq2)
    logFC.data               <- logFC.data[!is.na(logFC.data) & !is.infinite(logFC.data)]
    logFC.data.sorted        <- sort(logFC.data, decreasing = TRUE)
    logFC.data.sorted.id     <- seq(1:length(logFC.data.sorted))
    chamber.pseudo.binary    <- names(logFC.data.sorted) %in% chamber.set
    chamber.pseudo.id        <- logFC.data.sorted.id[chamber.pseudo.binary]
    names(logFC.data.sorted) <- logFC.data.sorted.id
    chamber.gsea.set         <- data.frame( diseaseId = rep(chamber,length(chamber.pseudo.id)), 
                                            geneId    = chamber.pseudo.id, check.names = TRUE)
    chamber.GSEA             <- GSEA(logFC.data.sorted, 
                                     TERM2GENE = chamber.gsea.set,
                                     minGSSize = 20, pvalueCutoff = 1)
    atrium.GSEA.plot         <- gseaplot(chamber.GSEA, geneSetID = chamber, 
                                         by = 'runningScore', 
                                         color.line = '#20BF55', color.vline = '#0B4F6C' ) +
                                theme(axis.title = element_text(size = 12),
                                      axis.text  = element_text(size = 8)) + 
                                labs(caption = '(data source: GSE64403)') +
                                theme( aspect.ratio = 1)
    return(atrium.GSEA.plot)
    
}

"

foo <- secure.GSE64403.GSEA(GSE1479.exprs)
"

secure.Sln.from.GSE64403 <- function(GSE64403.rsubread = GSE64403.results) {
    GSE64403.atrium.pValue             <- .get.isolated.cell.deseq2(GSE64403.rsubread)
    GSE64403.cardiac.cells             <- .get.rsubread.exprs(GSE64403.rsubread)[,7:12]
    GSE64403.cardiac.cells.condition   <- factor(rep(1:2,each = 3),
                                                 levels = 1:2,
                                                 labels = c('iDS','iDR'))
    atrium.tidy.data                   <- data.frame(exprs = GSE64403.cardiac.cells['Sln',], 
                                                     treat = GSE64403.cardiac.cells.condition )
    atrium.p.adj.Value                 <- paste0( 'italic(\'P\')[italic(\'adi\')]', 
                                                  ' ==',GSE64403.atrium.pValue['Sln',]$padj)
    prelim.figure.A <- ggplot(atrium.tidy.data, aes(x = treat, y = log(exprs))) + 
                       geom_jitter( aes(color = treat), 
                                    width = 0.05, height = 0.5, size = 3.0) +
                       scale_x_discrete(labels = c('sham','surgical')) +
                       scale_colour_manual( values = c('#1C3144','#D00000'), 
                             name   = 'treatment',labels = c('sham','surgical')) +
                       scale_y_continuous(breaks = seq(-5,10,2),limits = c(-2,10)) +
                       geom_segment(aes(x = 1.0, y = 9.0, xend = 1.0, yend = 9.5)) + 
                       geom_segment(aes(x = 1.0, y = 9.5, xend = 2.0, yend = 9.5)) + 
                       geom_segment(aes(x = 2.0, y = 9.0, xend = 2.0, yend = 9.5)) +
                       annotate( geom  = 'text', x = 1.5, y = 9.9, 
                                 label = atrium.p.adj.Value, parse = T) +
                       theme_classic()       +
                       theme(legend.position = c(0.8, 0.18))            +
                       theme(axis.title.x = element_blank(), legend.title.align = 0.5) +
                       labs(y = 'log 2(RPKM exprs value)') + 
                       theme( aspect.ratio = 1)
    
    return(prelim.figure.A)
}

get.AV.postSurgery.heatmap <- function(GSE64403.rsubread = GSE64403.results ) {
    av.genes               <- c('Sln','Tbx5','Nr2f2','Myl7','Myl4',
                                'Hey2','Irx4','Kcne1','Myh7','Pln')
    sample.info            <- data.frame( treat  = c('vAd','vAd','vD1S','vD1S','vD1R','vD1R',
                                                     'vD7R','vD7R','vD7S','vD7S'))
    rusubread.vsd.exprs    <- DESeqDataSetFromMatrix( countData = GSE64403.rsubread$counts[,27:36],
                                                      colData   = sample.info,
                                                      design    = ~ treat) %>% 
                              DESeq2::varianceStabilizingTransformation() %>%
                              assay()
    rsubread.gene.symbols  <- mapIds( org.Mm.eg.db, keys = GSE64403.rsubread$annotation$GeneID,
                                      keytype = 'ENSEMBL', column = 'SYMBOL') %>% 
                              make.names(unique = T)
    rownames(rusubread.vsd.exprs) <- rsubread.gene.symbols
    gene.norm.exprs               <- rusubread.vsd.exprs[av.genes,]
    column.annotation             <- data.frame(samples = factor(rep(1:5,each = 2), 
                                                                labels = c('vAd','vD1S',
                                                                           'vD1R','vD7R','vD7S')) )
    column.colors                 <- list(samples = c(vAd  = '#1C3114', vD1S = '#D00000',vD1R = '#FFBA08',
                                                      vD7R = '#A2AEBB', vD7S = '#3F88C5'))
    rownames(column.annotation)  <- colnames(gene.norm.exprs)
    heat.map <- pheatmap(gene.norm.exprs, scale = 'none', cluster_cols = F, cluster_rows = F,
                         color = colorRampPalette(c('#0B3954', '#BFD7EA', '#FF6663'))(50),
                         cellheight = 23, cellwidth = 6, labels_row = av.genes, 
                         border_color = NA, annotation_col = column.annotation,
                         show_colnames  = F, annotation_colors = column.colors, silent = T)
    return(heat.map)
    
}

get.AV.isolated.heatmap <- function(GSE64403.rsubread = GSE64403.results) {
    av.genes               <- c('Sln','Tbx5','Nr2f2','Myl7','Myl4',
                                'Hey2','Irx4','Kcne1','Myh7','Pln')
    sample.info            <- data.frame( treat  = c(rep('iP0',3),rep('iP4',3),rep('iP7S',3),
                                                     rep('iD7R',3)))
    rusubread.vsd.exprs    <- DESeqDataSetFromMatrix( countData = GSE64403.rsubread$counts[,1:12],
                                                      colData   = sample.info,
                                                      design    = ~ treat) %>% 
                              DESeq2::varianceStabilizingTransformation() %>%
                              assay()
    rsubread.gene.symbols  <- mapIds( org.Mm.eg.db, keys = GSE64403.rsubread$annotation$GeneID,
                                      keytype = 'ENSEMBL', column = 'SYMBOL') %>% 
                              make.names(unique = T)
    rownames(rusubread.vsd.exprs) <- rsubread.gene.symbols
    gene.norm.exprs               <- rusubread.vsd.exprs[av.genes,]
    column.annotation             <- data.frame(samples = factor(rep(1:4,each = 3),
                                                                 labels = c('iP0','iP4','iP7S','iP7R')))
    rownames(column.annotation)   <- colnames(gene.norm.exprs)
    column.colors                 <- list(samples = c(iP0  = '#D00000',iP4  = '#FFBA08',
                                                      iP7S = '#A2AEBB',iP7R = '#3F88C5'))
    heat.map <- pheatmap(gene.norm.exprs, scale = 'none', cluster_cols = F, cluster_rows = F,
                         cellheight = 23, cellwidth = 6, border_color = NA,
                         color = colorRampPalette(c('#0B3954', '#BFD7EA', '#FF6663'))(50),
                         labels_row     = av.genes, annotation_col = column.annotation,
                         show_colnames  = F, annotation_colors = column.colors, silent = T)
    return(heat.map)
}


get.figure.3 <- function(GSE64403.rsubread,GSE1479.exprs) {
    figure.3.A <- get.AV.postSurgery.heatmap(GSE64403.rsubread) %$% gtable
    figure.3.B <- secure.GSE64403.GSEA(GSE1479.exprs)
    figure.3.C <- get.AV.isolated.heatmap(GSE64403.rsubread) %$% gtable
    figure.3.D <- secure.Sln.from.GSE64403(GSE64403.rsubread) 

    figure.3.final <- ggdraw() +
                      draw_plot( figure.3.A, 
                                 x =.00, y = .5,  width = .5, height = .5) +
                      draw_plot( figure.3.B,   
                                 x =  0.5, y =  0.5, width = .5, height = .5) +
                      draw_plot( figure.3.C, 
                                 x = 0.0, y = 0.0,  width = .5, height = .5) +
                      draw_plot( figure.3.D,  
                                 x = .5, y =  0,  width = .5, height = .5) +
                      draw_plot_label( label = LETTERS[1:4],
                                       x     = c(0,.5,0,.5),
                                       y     = c(1,1,.5,.5),
                                       size  = 10)
    return(figure.3.final)
}




# TODO

#---
# Figure 4
#---


#----
# the Original data 
# from Broad institute
# http://software.broadinstitute.org/gsea/msigdb/cards/PID_RETINOIC_ACID_PATHWAY.html
#----
.get.probe.id <- function(gene.name) {
    #---
    # Sln mouse time expression pattern
    # mget('1448394_at',mouse4302GENENAME)
    # mget('17888',revmap(mouse4302ENTREZID))
    #---
    probe.ids   <- mget(gene.name,revmap(mgu74aSYMBOL)) %>% unlist() %>% .[1]
    return(probe.ids)
}

get.GSE775.limma.result <- function (GSE775.exprs,
                                     group = 'w.1', contrast = 'ilv - lv') {
    
    infarction.exprs <- GSE775.exprs
    group.selection  <- .get.GSE775.groups()[c(group)] %>% unlist()
    
    
    
    f.1                 <- factor( rep(1:3,each = 3),
                                   levels = 1:3,labels = c('lv','ilv','nilv'))
    design.1            <- model.matrix(~ 0 + f.1)
    colnames(design.1)  <- levels(f.1)
    contrast.matrix.1   <- makeContrasts(contrasts = contrast,levels = design.1)
    limma.result <- infarction.exprs[,group.selection] %>% 
        lmFit(design.1) %>% 
        contrasts.fit(contrast.matrix.1) %>%
        eBayes() %>% 
        topTable(adjust = 'BH', number = Inf)
    return(limma.result)
}

.secure.RA.signaling.pathway <- function(spiece = 'mouse', manual = T) {
    
    RA.components.human <- list(HDAC3 = 8841,	PRKCG = 5582, CDK1  = 983,  RBP1  = 5947, RARA = 5914,
                                RARB  = 5915, VDR	  = 7421, RARG  = 5916, PRKCA = 5578, PRKACA = 5566,
                                RXRA  = 6256, MAPK1 = 5594, RXRB  = 6257, AKT1  = 207 , MAPK8  = 5599,
                                RXRG  = 6258, NRIP1 = 8204, CDK7  = 1022, CCNH  = 902,  MNAT1  = 4331,
                                EP300 = 2033, HDAC1 = 3065, NCOA2 = 10499, NCOA1 = 8648, MAPK14 = 1432,
                                CREBBP = 1387, KAT2B = 8850, NCOR2 = 9612, NCOA3 = 8202 )
    RA.components.mouse <- list(Hdac3 = 15183, Prkcg = 18752, Cdk1 = 12534, Rbp1 = 19659, Rara = 19401,
                                Rarb  = 218772, Vdr = 22337, Rarg = 19411, Prkca = 18750, Prkaca = 18747,
                                Rxra  = 20181, Mapk1 = 26413, Rxrb = 20182, Akt1 = 11651, Mapk8 = 26419,
                                Rxrg  = 20183, Nrip1 = 268903, Cdk7 = 12572, Ccnh = 66671, Mnat1 = 17420,
                                Ep300 = 328572, Hdac1 = 433759, Ncoa2 = 17978, Ncoa1 = 17977, Mapk14 = 26416,
                                Crebbp = 12914, Kat2b = 18519, Ncor2 = 20602, Ncoa3 = 17979)
    RA.target.genes     <- list(Nr2c1 = 22025,Rara = 19401, Rarb  = 218772,Rarg = 19411,Hoxa1 = 15394,
                                Hoxb1 = 15407,Hoxa4 = 15401, Hoxb4 = 15412, Hoxd4 = 15436, Hoxb5 = 15413,
                                Cdx1 = 12590, Pax6 = 18508, Pdx1 = 18609, Pitx2 = 18741, Tlx2 = 21909,
                                Cebpe = 110794, Egr1 = 13653, Ets1 = 23871, Foxa1 = 15375,Hnf1a = 21405,
                                Neurog2 = 11924, Olig2 = 50913, Pou1f1 = 18736, Pou5f1 = 18999, Stat1 = 20846,
                                Adh7 = 11529, Aldh1a1 = 11668, Crabp2 = 12904, Cyp26a1 = 13082, Rbp1 = 19659,
                                Icam1 = 15894, Itgb3 = 16416, Lamb1 = 16777, Epo = 13856, Fgf8 = 14179,
                                Gh = 14599, Gnrh1 = 14714, Tshb = 22094, Lefty1 = 13590, Mdk = 17242,
                                Nodal = 18119, Oxt = 18429)
    if(manual == T) return(RA.target.genes)
    if (spiece == 'mouse') return(RA.components.mouse)
    if (spiece == 'human') return(RA.components.human)
    
}

.get.manual.curation.RA.genes <- function(spiece = 'mouse') {
    ra.mouse.genes <- list(Nr2c1 = 22025, Hoxa1 = 15394, Hoxb5 = 15413,Aldh1a1 = 11668,Aldh1a2 = 19378)
    ra.human.genes <- list(NR2C1 = 7181, HOXA1 = 3198,HOXB5 = 3215,ALDH1A1 = 216,ALDH1A2 = 395844)
    if (spiece == 'mouse') return(ra.mouse.genes)
    if (spiece == 'human') return(ra.human.genes)
}

.get.RA.infarction.tidydata <- function(data.affy, gene.name = 'Aldh1a2') {
    ra.gene.1       <- secure.tidy.data(data.affy, gene.name)
    ra.gene.1.probe <- .get.probe.id(gene.name)
    h.24.pvalue.1   <- get.GSE775.limma.result(data.affy, 
                                               group    = 'h.24', 
                                               contrast = 'ilv - lv') %>%
                       {.[ra.gene.1.probe,]} %$% adj.P.Val %>%
                       paste0( 'italic(\'P\')[italic(\'adi\')]', ' == ',.)
    h.24.pvalue.2   <- get.GSE775.limma.result(data.affy, 
                                               group    = 'h.24', 
                                               contrast = 'nilv - lv') %>%
                       {.[ra.gene.1.probe,]} %$% adj.P.Val %>%
                       paste0( 'italic(\'P\')[italic(\'adi\')]', ' == ',.)
    
    h.48.pvalue.1   <- get.GSE775.limma.result(data.affy, 
                                               group    = 'h.48', 
                                               contrast = 'ilv - lv') %>%
                       {.[ra.gene.1.probe,]} %$% adj.P.Val %>%
                       paste0( 'italic(\'P\')[italic(\'adi\')]', ' == ',.)
    h.48.pvalue.2   <- get.GSE775.limma.result(data.affy, 
                                               group    = 'h.48', 
                                               contrast = 'nilv - lv') %>%
                       {.[ra.gene.1.probe,]} %$% adj.P.Val %>%
                       paste0( 'italic(\'P\')[italic(\'adi\')]', ' == ',.)
    tidy.data       <- list( cleaned.data  = ra.gene.1,
                             h.24.pvalue.1 = h.24.pvalue.1,
                             h.24.pvalue.2 = h.24.pvalue.2,
                             h.48.pvalue.1 = h.48.pvalue.1,
                             h.48.pvalue.2 = h.48.pvalue.2)
    return(tidy.data)
    
}

get.RA.infarction.plot.1 <- function(data.affy, gene.name = 'Aldh1a2') {
    cleaned.data    <- .get.RA.infarction.tidydata(data.affy,gene.name)
    ra.gene.1       <- cleaned.data$cleaned.data
    h.24.pvalue.1   <- cleaned.data$h.24.pvalue.1
    h.24.pvalue.2   <- cleaned.data$h.24.pvalue.2
    h.48.pvalue.1   <- cleaned.data$h.48.pvalue.1
    h.48.pvalue.2   <- cleaned.data$h.48.pvalue.2
    ra.gene.1.plot  <- ggplot(ra.gene.1, aes(x = time, y = exprs)) +
                       scale_y_continuous(breaks = seq(5,10,1),limits = c(5,10)) +
                       geom_boxplot(aes(fill = region), 
                                    position = position_dodge(0.5)) +
                       scale_fill_manual( name   = 'regions',
                                          values = c( '#F5D547','#DB3069','#1446A0'),
                                          labels = c( 'lv','ilv','nilv')) +
                       geom_segment(aes(x = 2.8, y = 6.8, xend = 2.8, yend = 9.0)) + 
                       geom_segment(aes(x = 2.8, y = 9.0, xend = 3.0, yend = 9.0)) + 
                       geom_segment(aes(x = 3.0, y = 9.0, xend = 3.0, yend = 8.0)) +
                       annotate( geom  = 'text', x = 2.8, y = 9.1, 
                                 label =  h.24.pvalue.1, size = 3.5, parse = T) +
                       geom_segment(aes(x = 2.8, y = 5.5, xend = 2.8, yend = 6.1)) +
                       geom_segment(aes(x = 2.8, y = 5.5, xend = 3.2, yend = 5.5)) +
                       geom_segment(aes(x = 3.2, y = 5.5, xend = 3.2, yend = 6.2)) +
                       annotate( geom  = 'text', x = 3.0, y = 5.3, 
                                 label =  h.24.pvalue.2, size = 3.5, parse = T) +
                       geom_segment(aes(x = 3.8, y = 6.7, xend = 3.8, yend = 8.4)) + 
                       geom_segment(aes(x = 3.8, y = 8.4, xend = 4.0, yend = 8.4)) + 
                       geom_segment(aes(x = 4.0, y = 8.2, xend = 4.0, yend = 8.4)) +
                       annotate( geom  = 'text', x = 4.0, y = 8.6, 
                                 label =  h.48.pvalue.1, size = 3.5, parse = T) +
                       geom_segment(aes(x = 3.8, y = 5.8, xend = 3.8, yend = 6.0)) + 
                       geom_segment(aes(x = 3.8, y = 5.8, xend = 4.2, yend = 5.8)) + 
                       geom_segment(aes(x = 4.2, y = 5.8, xend = 4.2, yend = 6.6)) +
                       annotate( geom  = 'text', x = 4.2, y = 5.6, 
                                 label =  h.48.pvalue.2, size = 3.5, parse = T) +
                       theme_classic()       +
                       theme(legend.position = c(0.85, 0.7), legend.title.align = 0.5)            +
                       labs(y = 'log2 gene exprs') + 
                       theme( aspect.ratio = 1)
    return(ra.gene.1.plot)
                       
}

"
get.RA.infarction.plot.1(data.affy)
"

get.RA.infarction.plot.2 <- function(data.affy, gene.name = 'Pbx1') {
    cleaned.data    <- .get.RA.infarction.tidydata(data.affy,gene.name)
    ra.gene.1       <- cleaned.data$cleaned.data
    h.24.pvalue.1   <- cleaned.data$h.24.pvalue.1
    h.24.pvalue.2   <- cleaned.data$h.24.pvalue.2
    h.48.pvalue.1   <- cleaned.data$h.48.pvalue.1
    h.48.pvalue.2   <- cleaned.data$h.48.pvalue.2
    ra.gene.1.plot  <- ggplot(ra.gene.1, aes(x = time, y = exprs)) +
                       scale_y_continuous(breaks = seq(6,8,1),limits = c(6,8)) +
                       geom_boxplot(aes(fill = region), position = position_dodge(0.5)) +
                       geom_segment(aes(x = 2.8, y = 6.9, xend = 2.8, yend = 7.8)) + 
                       geom_segment(aes(x = 2.8, y = 7.8, xend = 3.1, yend = 7.8)) + 
                       geom_segment(aes(x = 3.1, y = 6.8, xend = 3.1, yend = 7.8)) +
                       annotate( geom  = 'text', x = 2.9, y = 7.9, 
                                 label =  h.24.pvalue.1, parse = T, size = 3.0) +
                       geom_segment(aes(x = 2.8, y = 6.2, xend = 2.8, yend = 6.3)) +
                       geom_segment(aes(x = 2.8, y = 6.2, xend = 3.3, yend = 6.2)) +
                       geom_segment(aes(x = 3.3, y = 6.2, xend = 3.3, yend = 6.4)) +
                       annotate( geom  = 'text', x = 3.0, y = 6.1, 
                                 label =  h.24.pvalue.2, parse = T, size = 3.0) +
                       geom_segment(aes(x = 3.8, y = 6.9, xend = 3.8, yend = 7.4)) + 
                       geom_segment(aes(x = 3.8, y = 7.4, xend = 4.0, yend = 7.4)) + 
                       geom_segment(aes(x = 4.0, y = 6.8, xend = 4.0, yend = 7.4)) +
                       annotate( geom  = 'text', x = 3.8, y = 7.5, 
                                 label =  h.48.pvalue.1, parse = T, size = 3.0) +
                       geom_segment(aes(x = 3.8, y = 6.2, xend = 3.8, yend = 6.6)) + 
                       geom_segment(aes(x = 3.8, y = 6.2, xend = 4.2, yend = 6.2)) + 
                       geom_segment(aes(x = 4.2, y = 6.2, xend = 4.2, yend = 6.5)) +
                       annotate( geom  = 'text', x = 3.9, y = 6.1, 
                                 label =  h.48.pvalue.2, parse = T, size = 3.0) +
                       theme_classic()       +
                       labs(y = 'log2 gene exprs') + 
                       theme( aspect.ratio = 1)
    return(ra.gene.1.plot)
    
}


"
get.RA.infarction.plot.1(data.affy)
.get.RA.infarction.tidydata(data.affy,'Hoxa1')
"

.get.RA.regeneration.tidy.data <- function(data.rsubread, gene.name = 'Aldh1a2') {
    GSE64403.heart.tissue              <- .get.rsubread.exprs(data.rsubread)[,27:36]
    GSE64403.cardiac.tissue.condition  <- factor( rep(1:5,each = 2),levels = 1:5,
                                                  labels = c('vAd','vD1S','vD1R','vD7R','vD7S'))
    GSE64403.tissue.exprs              <- newCountDataSet( GSE64403.heart.tissue, 
                                                           GSE64403.cardiac.tissue.condition ) %>%
                                          estimateSizeFactors() %>% estimateDispersions()
    GSE64403.tissue.pValue.1           <- GSE64403.tissue.exprs %>% nbinomTest('vD1S', 'vAd')
    GSE64403.tissue.pValue.2           <- GSE64403.tissue.exprs %>% nbinomTest('vD1R', 'vD1S')
    GSE64403.tissue.pValue.3           <- GSE64403.tissue.exprs %>% nbinomTest('vD1R', 'vAd')
    rsubread.gene.symbols              <- mapIds( org.Mm.eg.db, keys = data.rsubread$annotation$GeneID,
                                                  keytype = 'ENSEMBL', column = 'SYMBOL') %>% 
                                          make.names(unique = T)
    rownames(GSE64403.tissue.pValue.1)   <- rsubread.gene.symbols
    rownames(GSE64403.tissue.pValue.2)   <- rsubread.gene.symbols
    rownames(GSE64403.tissue.pValue.3)   <- rsubread.gene.symbols
    GSE64403.tidy.data   <- data.frame(groups = c('vAd','vD1S','vD1R'),
                                       exprs  = c(GSE64403.tissue.pValue.1[gene.name,]$baseMeanB,
                                                  GSE64403.tissue.pValue.1[gene.name,]$baseMeanA,
                                                  GSE64403.tissue.pValue.2[gene.name,]$baseMeanA))
    vD1R.vAd.pValue      <- GSE64403.tissue.pValue.3[gene.name,] %$% padj %>%
                            paste0( 'italic(\'P\')[italic(\'adi\')]', ' == ',.)
    vD1S.vAd.pValue      <- GSE64403.tissue.pValue.1[gene.name,] %$% padj %>%
                            paste0( 'italic(\'P\')[italic(\'adi\')]', ' == ',.)
    vD1R.vD1S.pValue     <- GSE64403.tissue.pValue.2[gene.name,] %$% padj %>%
                            paste0( 'italic(\'P\')[italic(\'adi\')]', ' == ',.)
    tidy.data            <- list( cleaned.data     = GSE64403.tidy.data,
                                  vD1R.vAd.pValue  = vD1R.vAd.pValue,
                                  vD1S.vAd.pValue  = vD1S.vAd.pValue,
                                  vD1R.vD1S.pValue = vD1R.vD1S.pValue)
    return(tidy.data)
    
}

get.RA.regeneration.plot.1 <- function(data.rsubread, gene.name = 'Aldh1a2' ) {
    cleaned.data        <- .get.RA.regeneration.tidy.data(data.rsubread, gene.name)
    GSE64403.tidy.data  <- cleaned.data$cleaned.data
    vD1R.vAd.pValue     <- cleaned.data$vD1R.vAd.pValue
    vD1S.vAd.pValue     <- cleaned.data$vD1S.vAd.pValue
    vD1R.vD1S.pValue    <- cleaned.data$vD1R.vD1S.pValue
    
    ra.gene.1.plot      <- ggplot(GSE64403.tidy.data, aes(x = groups, y = exprs)) +
                           geom_bar(aes(fill = groups),stat = 'identity', width = 0.7) + 
                           scale_fill_manual( name   = 'samples',
                                              values = c( '#F5D547','#DB3069','#1446A0'),
                                              labels = c( 'vAd','vD1R','vD1S')) +
                           scale_y_continuous(breaks = seq(0,1200,200),limits = c(0,1200)) +
                           geom_segment(aes(x = 1.1, y = 200, xend = 1.1, yend = 800)) + 
                           geom_segment(aes(x = 1.1, y = 800, xend = 2.0, yend = 800)) + 
                           geom_segment(aes(x = 2.0, y = 700, xend = 2.0, yend = 800)) +
                           annotate( geom  = 'text', x = 1.5, y = 840, 
                                     label =  vD1R.vAd.pValue, size = 3.5,  parse = T) +
                           geom_segment(aes(x = 2.2, y = 700, xend = 2.2, yend = 900)) + 
                           geom_segment(aes(x = 2.2, y = 900, xend = 2.9, yend = 900)) + 
                           geom_segment(aes(x = 2.9, y = 840, xend = 2.9, yend = 900)) +
                           annotate( geom  = 'text', x = 2.5, y = 940, 
                                     label =  vD1R.vD1S.pValue, size = 3.5, parse = T) +
                           geom_segment(aes(x = 0.9, y = 200, xend = 0.9, yend = 1100)) + 
                           geom_segment(aes(x = 0.9, y = 1100, xend = 3.1, yend = 1100)) + 
                           geom_segment(aes(x = 3.1, y = 840, xend = 3.1, yend = 1100)) +
                           annotate( geom  = 'text', x = 1.8, y = 1160, 
                                     label = vD1S.vAd.pValue, size = 3.5, parse = T) +
                           theme_classic()       +
                           theme(legend.position = c(0.95, 0.9), 
                                 legend.title.align = 0.5,
                                 axis.title.x = element_blank(),
                                 legend.title = element_text(size = 13),
                                 axis.title   = element_text(size = 12)) +
                           labs(y = 'expression value of RPKM') +
                           theme( aspect.ratio = 1)
    return(ra.gene.1.plot)

}

"
get.RA.regeneration.plot.1(data.rsubread)
"

get.RA.regeneration.plot.2 <- function(data.rsubread, gene.name = 'Pbx1' ) {
    cleaned.data        <- .get.RA.regeneration.tidy.data(data.rsubread, gene.name)
    GSE64403.tidy.data  <- cleaned.data$cleaned.data
    vD1R.vAd.pValue     <- cleaned.data$vD1R.vAd.pValue
    vD1S.vAd.pValue     <- cleaned.data$vD1S.vAd.pValue
    vD1R.vD1S.pValue    <- cleaned.data$vD1R.vD1S.pValue
    
    ra.gene.2.plot      <- ggplot(GSE64403.tidy.data, aes(x = groups, y = exprs)) +
        geom_bar(aes(fill = groups),stat = 'identity') + 
        scale_y_continuous(breaks = seq(0,1200,200),limits = c(0,1200)) +
        geom_segment(aes(x = 1.1, y = 200, xend = 1.1, yend = 800)) + 
        geom_segment(aes(x = 1.1, y = 800, xend = 2.0, yend = 800)) + 
        geom_segment(aes(x = 2.0, y = 700, xend = 2.0, yend = 800)) +
        annotate( geom  = 'text', x = 1.5, y = 840, label =  vD1R.vAd.pValue, parse = T) +
        geom_segment(aes(x = 2.2, y = 700, xend = 2.2, yend = 900)) + 
        geom_segment(aes(x = 2.2, y = 900, xend = 2.9, yend = 900)) + 
        geom_segment(aes(x = 2.9, y = 840, xend = 2.9, yend = 900)) +
        annotate( geom  = 'text', x = 2.5, y = 940, label =  vD1R.vD1S.pValue, parse = T) +
        geom_segment(aes(x = 0.9, y = 200, xend = 0.9, yend = 1100)) + 
        geom_segment(aes(x = 0.9, y = 1100, xend = 3.1, yend = 1100)) + 
        geom_segment(aes(x = 3.1, y = 840, xend = 3.1, yend = 1100)) +
        annotate( geom  = 'text', x = 1.8, y = 1160, label = vD1S.vAd.pValue, parse = T) +
        theme( aspect.ratio = 1)
    return(ra.gene.2.plot)
    
}


"
get.RA.regeneration.plot.1(data.rsubread)
.get.RA.regeneration.tidy.data(data.rsubread,'Rarg')

"



# TODO
RA.signaling.GSEA.Regeneration <- function(GSE64403.rsubread) {
    GSE64403.heart.tissue              <- .get.rsubread.exprs(GSE64403.rsubread)[,27:36]
    GSE64403.cardiac.tissue.condition  <- factor( rep(1:5,each = 2),levels = 1:5,
                                                  labels = c('vAd','vD1S','vD1R','vD7R','vD7S'))
    GSE64403.tissue.exprs              <- newCountDataSet( GSE64403.heart.tissue, 
                                                           GSE64403.cardiac.tissue.condition ) %>%
                                          estimateSizeFactors() %>% estimateDispersions()
    GSE64403.tissue.pValue.1           <- GSE64403.tissue.exprs %>% nbinomTest('vD1S', 'vAd')
    GSE64403.tissue.pValue.2           <- GSE64403.tissue.exprs %>% nbinomTest('vD1R', 'vD1S')
    GSE64403.tissue.pValue.3           <- GSE64403.tissue.exprs %>% nbinomTest('vD1R', 'vAd')
    rsubread.gene.symbols              <- mapIds( org.Mm.eg.db, keys = GSE64403.rsubread$annotation$GeneID,
                                                  keytype = 'ENSEMBL', column = 'SYMBOL') %>% 
                                          make.names(unique = T)
    rownames(GSE64403.tissue.pValue.1)   <- rsubread.gene.symbols
    rownames(GSE64403.tissue.pValue.2)   <- rsubread.gene.symbols
    rownames(GSE64403.tissue.pValue.3)   <- rsubread.gene.symbols
    GSE64403.logFC.data           <- GSE64403.tissue.pValue.1$log2FoldChange
    names(GSE64403.logFC.data)    <- GSE64403.tissue.pValue.1$id
    GSE64403.logFC.data           <- GSE64403.logFC.data[is.finite(GSE64403.logFC.data)]
    GSE64403.logFC.sorted.data    <- sort(GSE64403.logFC.data, decreasing = T)
    RA.set                        <- .secure.RA.signaling.pathway() %>% names()
    RA.gsea.set                   <- data.frame( diseaseId = rep('RA',length(RA.set)), 
                                                 geneId    = RA.set, check.names = TRUE)
    RA.GSEA               <- GSEA(GSE64403.logFC.sorted.data, 
                                  TERM2GENE = RA.gsea.set, minGSSize = 15, 
                                  pvalueCutoff = 1)
}

"
RA.signaling.GSEA.Regeneration(GSE64403.results)

"

RA.signaling.GSEA <- function( GSE775.exprs, GSE1479.exprs, 
                               group = 'w.1', contrast = 'ilv - lv') {
    limma.result     <- get.GSE775.limma.result(GSE775.exprs, group, contrast) 
    RA.set           <- .secure.RA.signaling.pathway() %>% names()
    RA.index         <- rownames(limma.result) %>% getSYMBOL('mgu74a') %in% RA.set
    RA.probe.id      <- rownames(limma.result)[RA.index]
    t.statistic      <- limma.result$t[RA.index] %>% 
                        sum()/sum(RA.index) * sqrt(sum(RA.index))
    logFC.data            <- limma.result$logFC
    names(logFC.data)     <- rownames(limma.result)
    logFC.data.sorted     <- sort(logFC.data, decreasing = TRUE)
    logFC.data.sorted.id  <- seq(1:length(logFC.data.sorted))
    RA.pseudo.binary      <- names(logFC.data.sorted) %>% 
                             getSYMBOL('mgu74a') %in% RA.set
    RA.pseudo.id          <- logFC.data.sorted.id[RA.pseudo.binary]
    names(logFC.data.sorted) <- logFC.data.sorted.id
    RA.gsea.set           <- data.frame( diseaseId = rep('RA',length(RA.pseudo.id)), 
                                         geneId    = RA.pseudo.id, check.names = TRUE)
    RA.GSEA               <- GSEA(logFC.data.sorted, 
                                  TERM2GENE = RA.gsea.set, minGSSize = 15, 
                                  pvalueCutoff = 1)
    return(RA.GSEA)
}


#---
# https://stackoverflow.com/questions/53524503/can-ggplot2-produce-flowcharts
# https://stackoverflow.com/questions/52820756/roundrectgrob-equivalent-for-geom-rect-in-ggplot
# https://stackoverflow.com/questions/4973898/combining-paste-and-expression-functions-in-plot-labels
#---
get.hypothesis.schema.plot <- function() {
    # paste0( 'italic(\'P\')[italic(\'adi\')]', ' == ',.)
    node.a.label     <- 'immune switch'
    edge.a.label     <- 'Raldh2'
    node.b.label     <- 'patterning & morphogenesis'
    edge.b.label     <- 'a.v reporgramming'
    node.c.label.1   <- 'atriumization'
    node.c.label.2   <- '(increased expression of atrium genes)'
    schema.ggplot  <- ggplot() +
                      scale_x_continuous(limits = c(0,110)) +
                      scale_y_continuous(limits = c(0,3)) +
                      statebins:::geom_rrect(aes(xmin = 0,xmax = 16.0,
                                    ymin = 1.85,ymax = 2.1, fill = '#1A535C') ,
                                    inherit.aes = F) +
                      geom_text(aes(x = 8.0, y = 2.0, label = node.a.label,
                                    fontface = 'italic'), parse = F) +
                      geom_segment(aes(x = 17, y = 2, xend = 23, yend = 2),
                                   arrow = arrow(length = unit(0.1, 'inches')) ) +
                      geom_text(aes(x = 20, y = 1.9, label = edge.a.label,
                                     fontface = 'italic'), parse = F) +
                      geom_segment(aes(x = 24, y = 1.86, xend = 24, yend = 1.94),
                                   arrow = arrow(length = unit(0.08, 'inches')) ) +
                      statebins:::geom_rrect(aes(xmin = 26,xmax = 58,
                                    ymin = 1.85,ymax = 2.1, fill = 'grey')) +
                      geom_text(aes(x = 42, y = 2.0, label = node.b.label, 
                                    fontface = 'italic') ,parse = F) +
                      statebins:::geom_rrect(aes(xmin = 70,xmax = 102,
                                ymin = 1.85, ymax = 2.1,fill = 'grey' )) +
                      geom_segment(aes(x = 59, y = 2, xend = 69, yend = 2),
                                   arrow = arrow(length = unit(0.1, 'inches')) ) +
                      geom_text(aes(x = 85.5, y = 2.05, label = node.c.label.1,
                                    fontface = 'italic'),parse = F) +
                      geom_text(aes(x = 85.5, y = 1.95, label = node.c.label.2),
                                    fontface = 'italic', size = 3,parse = F) +
                      theme_nothing()

   return(schema.ggplot) 
}


get.figure.4 <- function(data.affy,data.rsubread) {
    figure.4.A <- get.RA.infarction.plot.1(data.affy,'Aldh1a2')
    figure.4.B <- get.RA.regeneration.plot.1(data.rsubread,'Aldh1a2')
    figure.4.C <- get.hypothesis.schema.plot()
    
    figure.4.final <- ggdraw() +
                      draw_plot( figure.4.A, 
                                 x =.0, y = .5,  width = .5, height = .5) +
                      draw_plot( figure.4.B,   
                                 x =  0.5, y =  0.5, width = .5, height = .5) +
                      draw_plot( figure.4.C, 
                                 x = 0.0, y = 0.0,  width = 1.0, height = .5) +
                      draw_plot_label( label = LETTERS[1:3],
                                 x     = c(0,.5,0),
                                 y     = c(1,1,.5),
                                 size  = 10)
    return(figure.4.final)
    
}


#---
# double check
#---

"
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE775

 	This dataset is a time series (1 hour [h], 4 hours, 24 hours, 48 hours, 1 week [w], and 8 weeks) 
 	intended to compare normal functioning left ventricles [lv + lv2] with infarcted [ilv] and 
 	non-infarcted left ventricles [nilv]. Ilv samples are taken from the region between the LAD 
 	artery and the apex on a mouse with myocardial infarction. Lv2 samples are from the same region 
 	in a sham operated mouse. Nilv samples are taken from the region above the infartion and the 
 	left ventricle [lv] samples mimic that region in a sham mouse. The lv and lv2 samples can be 
 	compared as both are from normal functioning hearts.
 	
"

.datset.GSE775.check <- function(data.path = GSE775.raw.path) {
    setwd(onedrive.path)
    GSE775.exprs  <- get.microarray.exprs(data.path)
    gene.symbols  <- rownames(GSE775.exprs) %>%
                     getSYMBOL('mgu74a') %>%
                     make.names(unique = T)
    group.names   <- c('lv_1h_1','lv_1h_2','lv_1h_3',
                       'lv_1w_1','lv_1w_2','lv_1w_3',
                       'lv_24h_1','lv_24h_2','lv_24h_3',
                       'lv_48h_1','lv_48h_2','lv_48h_3',
                       'lv_4h_1','lv_4h_2','lv_4h_3',
                       'lv_8w_1','lv_8w_2','lv_8w_3',
                       'lv2_1h','lv2_1w','lv2_24h','lv2_48h','lv2_4h',
                       'MI_ilv_1h_1','MI_ilv_1h_2','MI_ilv_1h_3',
                       'MI_ilv_1w-1','MI_ilv_1w_2','MI_ilv_1w_3',
                       'MI_ilv_24h','MI_ilv_24h_2','MI_ilv_24h_3',
                       'MI_ilv_48h_1','MI_ilv_48h_2','MI_ilv_48h_3',
                       'MI_ilv_4h_1','MI_ilv_4h_2','MI_ilv_4h_3',
                       'MI_ilv_8w_1','MI_ilv_8w_2','MI_ilv_8w_3',
                       'MI_nilv_1h_1','MI_nilv_1h_2','MI_nilv_1h_3',
                       'MI_nilv_1w_1','MI_nilv_1w_2','MI_nilv_1w_3',
                       'MI_nilv_24h_1','MI_nilv_24h_2','MI_nilv_24h_3',
                       'MI_nilv_48h_1','MI_nilv_48h_2','MI_nilv_48h_3',
                       'MI_nilv_4h_1','MI_nilv_4h_2','MI_nilv_4h_3',
                       'MI_nilv_8w_1','MI_nilv_8w_2','MI_nilv_8w_3')
    colnames(GSE775.exprs) <- group.names
    groups                 <- .get.GSE775.groups() %>% unlist()
    rownames(GSE775.exprs) <- gene.symbols
    GSE775.final.data      <- GSE775.exprs[,groups]
    GSE775.set             <- createWorkbook(title = 'GSE775.exprs.data') 
    addWorksheet(GSE775.set, 'microarray expression value') 
    writeData(GSE775.set, 1, GSE775.final.data, rowNames = T) 
    saveWorkbook(GSE775.set, file = 'GSE775.xlsx', overwrite = TRUE)
    
}

"
vAd  : Adult ventricular myocardium
vD7S : 7 days post-sham surgery ventricular myocardium
"

.datset.GSE64403.check <- function(data.rsubread = GSE64403.results) {
    setwd(onedrive.path)
    GSE64403.RPKM <- .get.rsubread.exprs(data.rsubread)
    group.names   <- c('iP0_1','iP0_2','iP0_3','iP4_1','iP4_2','iP4_3',
                       'iD7S_1','iD7S_2','iD7S_3','iD7R_1','iD7R_2','iD7R_3',
                       'ex0hr_1','ex0hr_2','ex24hr_1','ex24hr_2','ex48hr_1','ex48hr_2',
                       'ex72hr_1','ex72hr_2','vP0_1','vP0_2','vP4_1','vP4_2',
                       'vP7_1','vP7_2','vAd_1','vAd_2','vD1S_1','vD1S_2','vD1R_1','vD1R_2',
                       'vD7R_1','vD7R_2','vD7S_1','vD7S_1')
    colnames(GSE64403.RPKM)  <- group.names
    GSE64403.RPKM['Aldh1a2',]
    GSE64403.set             <- createWorkbook(title = 'GSE64403.RPKM.data') 
    addWorksheet(GSE64403.set , 'GSE64403 RPKM value') 
    writeData(GSE64403.set , 1, GSE64403.RPKM, rowNames = T) 
    saveWorkbook(GSE64403.set , file = 'GSE64403.RPKM.xlsx', overwrite = TRUE)
}

export_excel_av_genes <- function(GSE1479.exprs) {
    av.genes <- get.av.genes.from.DGE(GSE1479.exprs)
    setwd(onedrive.path)
    atrium.set.probeID      <- av.genes$atrium.probeID
    atrium.set.symbol       <- mget(atrium.set.probeID,mouse4302SYMBOL,ifnotfound = NA) %>% unlist
    atrium.set.EntrezID     <- mget(atrium.set.probeID,mouse4302ENTREZID,ifnotfound = NA) %>% unlist
    ventricle.set.probeID   <- av.genes$ventricle.probeID
    ventricle.set.symbol    <- mget(ventricle.set.probeID,mouse4302SYMBOL,ifnotfound = NA) %>% unlist
    ventricle.set.EntrezID  <- mget(ventricle.set.probeID,mouse4302ENTREZID,ifnotfound = NA) %>% unlist
    atrium.set              <- data.frame(probe.ID    = atrium.set.probeID,
                                          EntrezID    = atrium.set.EntrezID,
                                          gene.symbol = atrium.set.symbol)
    ventricle.set           <- data.frame(probe.ID    = ventricle.set.probeID,
                                          EntrezID    = ventricle.set.EntrezID,
                                          gene.symbol = ventricle.set.symbol)
    chamber.set <- createWorkbook(title = 'chamber.Specific.Gene') 
    addWorksheet(chamber.set, 'atriumGenes') 
    addWorksheet(chamber.set, 'ventricleGenes') 
    writeData(chamber.set, 1, atrium.set) 
    writeData(chamber.set, 2, ventricle.set) 
    saveWorkbook(chamber.set, file = 'avGenelist.xlsx', overwrite = TRUE)

}


# function already written and perform real analysis now!

# align the RNA-seq raw data
align.raw.data(raw.data.layout = raw.data.dir, 
               index.dir       = rsubread.index.dir, 
               index.file      = rsubread.index.file )
# get the featureCount result using Rsubread!
GSE64403.results <- get.featureCount.result(raw.data.layout = raw.data.dir, 
                                            GTF.file        = mRNA.GRCm38.GTF.file)

"
get the final figures!
"

GSE775.result    <- get.microarray.exprs(GSE775.raw.path)
GSE1479.result   <- get.microarray.exprs(GSE1479.raw.path)
figure.1         <- get.figure.1(GSE1479.result)
figure.2         <- get.figure.2(GSE775.result,GSE1479.result)
figure.3         <- get.figure.3(GSE64403.results,GSE1479.result)
figure.4         <- get.figure.4(GSE775.result,GSE64403.results)

"
save the figures!
"
ggsave(filename = 'Figure.1.png', plot = figure.1, dpi = 600, unit = 'mm',
       scale = 1, width = 300, height = 250, limitsize = F, path = onedrive.path)
ggsave(filename = 'Figure.2.png', plot = figure.2, dpi = 600, unit = 'mm',
       scale = 1, width = 300, height = 250, limitsize = F, path = onedrive.path)
ggsave(filename = 'Figure.3.png', plot = figure.3, dpi = 600, unit = 'mm',
       scale = 1, width = 300, height = 250, limitsize = F, path = onedrive.path)
ggsave(filename = 'Figure.4.png', plot = figure.4, dpi = 600, unit = 'mm',
       scale = 1, width = 300, height = 250, limitsize = F, path = onedrive.path)


"
export the original expression data
in order to double check the expression trend.

"
.datset.GSE775.check()
.datset.GSE64403.check()

"

data.affy          <- GSE775.result 
data.rsubread      <- GSE64403.results
GSE775.exprs       <- GSE775.result
GSE1479.exprs      <- GSE1479.result
GSE64403.rsubread  <- GSE64403.results



 

RA.signaling.GSEA(GSE775.result,GSE1479.result,'w.8')

#.get.probe.id('Sln')




setwd(project.layout$rdata)
save.image('GSE64403.Rdata')
quit('no')