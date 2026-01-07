
library(data.table)
library(ggplot2)



################################################################################
# load GSEA gene signatures
################################################################################

load.gene.signatures=function(sigs=NULL, recreate=F){
  
library(hypeR)
library(stringr)

    #updated chemical perturbations signatures
    simpleCache::simpleCache("all_gene_signatures" %>% addversion, {
        cp.sigs= qusage::read.gmt(getExternalFile("l1000_cp.gmt", "https://lincs-dcic.s3.amazonaws.com/LINCS-sigs-2021/gmt/l1000_cp.gmt"))
  
  lig.sigs= qusage::read.gmt(getExternalFile("l1000_lig.gmt","https://lincs-dcic.s3.amazonaws.com/LINCS-sigs-2021/gmt/l1000_lig.gmt"))

      
    gsetpars=list(
    CELLTYPES=list(set=msigdb_gsets(species="Homo sapiens", category="C8"), name="CELLTYPES"),
    PANGLAO=list(set=hypeR::enrichr_gsets(genesets = "PanglaoDB_Augmented_2021"), name="PANGLAODB"),
    DESCARTES=list(set=hypeR::enrichr_gsets(genesets = "Descartes_Cell_Types_and_Tissue_2021"), name="DESCARTES"),
    WIKIPATHWAYS=list(set=hypeR::enrichr_gsets(genesets ="WikiPathway_2021_Human"), name="WIKIPATHWAYS"),
    GOMF=list(set=msigdb_gsets(species="Homo sapiens", category="C5", subcategory="GO:MF")), 
    GPDOWN=list(set=hypeR::enrichr_gsets(genesets = "Gene_Perturbations_from_GEO_down"), name="GPDOWN"), 
    GPUP=list(set=hypeR::enrichr_gsets(genesets = "Gene_Perturbations_from_GEO_up"), name="GPUP"),
      TFPERT=list(set=hypeR::enrichr_gsets(genesets = "TF_Perturbations_Followed_by_Expression"), name="TFPERT"),
      #=list(set=hypeR::enrichr_gsets(genesets = "Kinase_Perturbations_from_GEO_down"), name="KINPERTDOWN"),
      #KINPERTUP=list(set=hypeR::enrichr_gsets(genesets = "Kinase_Perturbations_from_GEO_up"), name="KINPERTUP"),  
      CHEMPERTDOWN=list(set=hypeR::gsets$new(cp.sigs[grepl("down$", names(cp.sigs))], name="LINCS1000_chempert_down_downloaded", version="v2021"), name="CHEMPERTDOWN"),
      CHEMPERTUP=list(set=hypeR::gsets$new(cp.sigs[grepl("up$", names(cp.sigs))], name="LINCS1000_chempert_up_downloaded", version="v2021"), name="CHEMPERTUP"),
      #CHEMPERTDOWN=list(set=hypeR::enrichr_gsets(genesets = "LINCS_L1000_Chem_Pert_down"), name="CHEMPERTDOWN"),
      #CHEMPERTUP=list(set=hypeR::enrichr_gsets(genesets = "LINCS_L1000_Chem_Pert_up"), name="CHEMPERTUP"),
      LIGPERTDOWN=list(set=hypeR::enrichr_gsets(genesets = "Ligand_Perturbations_from_GEO_down"), name="LIGPERTDOWN"),
      LIGPERTUP=list(set=hypeR::enrichr_gsets(genesets = "Ligand_Perturbations_from_GEO_up"), name="LIGPERTUP"), 
      LIGPERTLINCSDN=list(set=hypeR::gsets$new(lig.sigs[grepl("down$", names(lig.sigs))], name="LINCS1000_ligandpert_down_downloaded", version="v2021"), name="LIGPERTLINCSDN"),
      LIGPERTLINCSUP=list(set=hypeR::gsets$new(lig.sigs[grepl("up$", names(lig.sigs))], name="LINCS1000_ligandpert_up_downloaded", version="v2021"), name="LIGPERTLINCSUP"),
      CGP=list(set=msigdb_gsets(species="Homo sapiens",category="C2", subcategory = "CGP" ), name="CGP")
    #ARCHS4CELL_LINES=list(set=ARCHCELLLINES, name="ARCHS4_Cell-Lines")
  )
    }, assignToVar="gsetpars", recreate=recreate)
  
    fcat("Existent gene sets are:", paste(names(gsetpars), collapse=","))
if(is.null(sigs)){
  gsetpars
  }else{
  gsetpars[sigs]
  } 

}


resultsDir <- function(...) {
  paste0(config$out_root, "/atac/", ...)
}

pipelineDir <- function(...) {
  paste0(config$out_root, "/pipe/", ...)
}

getBamPath <- function(pipe_name, genome_ver = basename(config$genome_chromatin)) {
  return(pipelineDir("out/", pipe_name, "/aligned_", genome_ver, "/", pipe_name, "_sort_dedup.bam"))
}

getPeaksPath <- function(pipe_name, genome_ver = basename(config$genome_chromatin)) {
  return(pipelineDir("out/", pipe_name, "/peak_calling_", genome_ver, "/", pipe_name, "_peaks.narrowPeak"))
}

dtToGr <- function(dt, chrCol="chrom", startCol="start", endCol="end", metaCols=c()) {
  library("GenomicRanges")
  
  argList <- list()
  for(n in metaCols) {
    argList[[n]] <- dt[,get(n)]
  }
  
  argList$ranges <- IRanges(dt[,as.numeric(get(startCol))],dt[,as.numeric(get(endCol))])
  argList$seqnames <- dt[,get(chrCol)]
  
  do.call(GRanges, args=argList)
}
grToDt <- function(gr, chrCol="chrom", startCol="start", endCol="end") {		
  dt <- data.table(chrom=as.character(seqnames(gr)), start=start(gr), end=end(gr))
  setnames(dt, c(chrCol, startCol, endCol))
  dt
}

liftOver <- function(gr, from="hg19", to=config$genome_build) {
  
  if(from!=to) {
    # download liftOver chain file:
    chnF <- resultsDir("downloaded_files/", from, "To", capFirst(to), ".over.chain")
    if(!file.exists(chnF)) {
      dir.create(resultsDir("downloaded_files"), showWarnings=FALSE)
      download.file(paste0("https://hgdownload-test.gi.ucsc.edu/goldenPath/",from,"/liftOver/", basename(chnF), ".gz"), destfile=paste0(chnF,".gz"))
      R.utils::gunzip(paste0(chnF,".gz"))
    }
    
    chn <- rtracklayer::import.chain(chnF)	
    gr <- rtracklayer::liftOver(gr, chn)
  }
  
  return(gr)
}
simpleMotifName <- function(x) {
  gsub("^(.+) .+$", "\\1",x)
}

# adapted from https://github.com/satijalab/seurat/blob/a1294c4d363780548dbf9cc4a4abb3a6078a6d64/R/utilities.R
# to make it such that it can be run on a matrix without dependecies on other Seurat features
# (also made it a lot faster)
calculateModuleScore <- function(
    assay.data,
    features,
    nbin = 24,
    ctrl = 100,
    name = 'Cluster',
    seed = 1,
    nthread = 8
) {
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  features.old <- features
  
  cluster.length <- length(x = features)
  
  data.avg <- Matrix::rowMeans(x = assay.data)
  data.avg <- data.avg[order(data.avg)]
  data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e30, n = nbin, labels = FALSE, right = FALSE)
  names(x = data.cut) <- names(x = data.avg)
  
  # merged three for-loops into one parallelized apply function for speed:
  feature.scores.allinone <- BiocParallel::bplapply(1:length(features),  function(i) {
    features.use <- features[[i]]
    # pick out `ctrl` random features at a similar expression level as controls for each feature in the query set:
    ctrl.use <- unique(unlist(lapply(features.use, function(f) {
      names(x = sample(
        x = data.cut[which(x = data.cut == data.cut[f])],
        size = ctrl,
        replace = FALSE
      ))
    })))
    
    ctrl.scores <- Matrix::colMeans(x = assay.data[ctrl.use, , drop = FALSE])
    features.scores <- Matrix::colMeans(x = assay.data[features.use, , drop = FALSE])
    
    return(features.scores - ctrl.scores)
  }, BPPARAM=BiocParallel::MulticoreParam(nthread))
  
  names(feature.scores.allinone) <- names(features)
  return(simplify2array(feature.scores.allinone))
}

options("RCACHE.DIR"=paste0(config$out_root, "/rcache/"))

CACHE_GENE_ANNOT <- "gene_annot"
CACHE_ATAC_PEAKS <- "atac_peaks"
CACHE_ATAC_PEAKS_ANNOTATED <- "atac_peaks_annot"
CACHE_ATAC_DDS <- "atac_dds"
CACHE_ATAC_DDS_RES <- "atac_dds_res"
CACHE_ATAC_META <- "atac_meta"
CACHE_ATAC_REGION_SETS <- "atac_roi"
CACHE_ATAC_COMPARISONS <- "atac_cmps"
CACHE_ATAC_MOTIF_LABELS <- "atac_mtf_lbl"
CACHE_ATAC_ENRICHDB_PREFIX <- "atac_enrich"
CACHE_ATAC_GENE_ASSIGNMENTS <- "atac_gene_assign"
CACHE_ATAC_ENRICHMENT_RESULTS <- "atac_enrich_res"
CACHE_ATAC_REGION_ORDER <- "atac_region_order"
#CACHE_RNA_META <- "rna_meta"
#CACHE_RNA_DDS <- "rna_dds"
#CACHE_RNA_NORM <- "rna_norm"
#CACHE_RNA_EXPR <- "rna_expr"
CACHE_SCRNA_MARKERS <- "scrna_markers"
CACHE_SCRNA_DATA <- "scrna_data"
CACHE_SCRNA_DATA_AGGREGATED <- "scrna_data_agg"
CACHE_SCRNA_EXPR <- "scrna_expr"
CACHE_SCRNA_UMAP <- "scrna_umap"
CACHE_SCRNA_MUT_GENES <- "scrna_mut_genes"
CACHE_SCRNA_MUT_SCORE <- "scrna_mut_score"

dir.create(resultsDir(), showWarnings=FALSE)

################################################################################
# full annotation of a peaks matrix in one go
################################################################################

annotate.peaks.global=function(peaks_dt_global, analysis.version){
region.annot.list=list()

msg("Meulemann et al. (2020)")
# https://www.meuleman.org/research/dhsindex/
n <- "dhs"

simpleCache(paste0("region_annotations_", n, analysis.version),{
#annot <- versionedCache(n, instruction={
	# download and read Regulatory Index annotations:
	dhsIndexF <- getExternalFile("meulemann_2020_DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz", "https://www.meuleman.org/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz")
	dhsIndex <- fread(dhsIndexF)
	# convert to GRanges and find overlaps with peaks:
	dhsIndex <- dtToGr(dhsIndex, chrCol="seqname", metaCols=c("component"))	
	o <- GenomicRanges::findOverlaps(peaks_dt_global, dhsIndex)
	# collapse multiple overlaps into one string:
	o2 <- as.data.table(o)[,.(dhs_type=unique(dhsIndex$component[subjectHits])),by=queryHits][order(dhs_type),.(annot=paste0(dhs_type, collapse=";")), by=.(rid=peaks_dt_global$Geneid[queryHits])]
	o2
	}, assignToVar="o2",reload=T )
	
region.annot.list[[n]]=o2

regionAnnotTypes<-  n

msg("Zerbino et al. (2015)")
# http://www.ensembl.org/info/genome/funcgen/regulatory_build.html
n <- "ensreg"

simpleCache(paste0("region_annotations_", n, analysis.version),{

	# download and read Regulatory Build:
  ensRegBuildF <- getExternalFile("zerbino_2015_homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20210107.gff.gz", "http://ftp.ensembl.org/pub/release-104/regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20210107.gff.gz")
	ensRegBuild <- rtracklayer::import(ensRegBuildF)
	seqlevels(ensRegBuild) <- paste0("chr",seqlevels(ensRegBuild))
	# find overlaps with peaks:
	o <- GenomicRanges::findOverlaps(peaks_dt_global, ensRegBuild)
	# collapse multiple overlaps into one string:
	o3 <- as.data.table(o)[,.(ensreg_type=unique(ensRegBuild$feature_type[subjectHits])),by=queryHits][order(ensreg_type),.(ensreg_type=paste0(ensreg_type, collapse=";")), by=.(rid=peaks_dt_global$Geneid[queryHits])]

	o3
}, assignToVar="o3",reload=T )
region.annot.list[[n]]=o3	

regionAnnotTypes <- c(regionAnnotTypes, n)

msg("Zhang et al. (2021)")
# http://catlas.org/humanenhancer/#!/
# https://www.cell.com/cell/fulltext/S0092-8674(21)01279-4#supplementaryMaterial
n <- "catlas"

simpleCache(paste0("region_annotations_", n, analysis.version),{
#annot <- versionedCache(n, instruction={
	# download and read 1.1mio cCRE's:
	catlasF <- getExternalFile("zhang_2021_cCRE_hg38.tsv.gz", "http://renlab.sdsc.edu/kai/Key_Processed_Data/cCRE_hg38.tsv.gz")
	catlas <- fread(catlasF)
	
	# get simple categorization of cCRE modules (the only place where I found this was the header of a supplementary table with GO enrichments):
	metaF <- getExternalFile("zhang_2021_cCRE_hg38_meta.xlsx", "https://data.mendeley.com/public-files/datasets/yv4fzv6cnm/files/40c64a61-4e73-42c4-9067-b28c507410d8/file_downloaded")
	catlasMeta <- readxl::read_xlsx(metaF)
	if(!identical(as.integer(catlasMeta[1,-1]),1:150)) stop("CRE cluster metadta not ordered as expected")
	catlasClusts <- tolower(gsub("(\\s|\\&)+","_",gsub("\\.+\\d+$","",colnames(catlasMeta)[-1]))) # strings are in order of cluster number
		
	catlas[ , cre_name := paste_("CRE", `CRE module`, catlasClusts[as.numeric(`CRE module`)])]
		
	# find overlaps with peaks:
	o <- findOverlaps(peaks_dt_global, dtToGr(catlas, "#Chromosome", "hg38_Start", "hg38_End"))
	# collapse multiple overlaps into one string:
	o4 <- as.data.table(o)[,.(catlas_type=unique(catlas$cre_name[subjectHits])),by=queryHits][order(catlas_type),.(catlas_type=paste0(catlas_type, collapse=";")), by=.(rid=peaks_dt_global$Geneid[queryHits])]

	o4	
}, assignToVar="o4",reload=T )
	
region.annot.list[[n]]=o4

#}, buildEnvir=list(peaks_dt=peaks_dt_global[, .(chrom, start, end, rid)]))
regionAnnotTypes <- c(regionAnnotTypes, n)

msg("Gao et al. (2019)")
# https://academic.oup.com/nar/article/48/D1/D58/5628925
# http://www.enhanceratlas.org/index.php
n <- "enhatlas"

simpleCache(paste0("region_annotations_", n, analysis.version),{
#annot <- versionedCache(n, instruction={
	# download and read 1.1mio cCRE's:
	tmp <- getExternalFile("enhanceratlas_2.0_Species_enhancer.RData", "http://www.enhanceratlas.org/data/download/Species_enhancer.RData")
	tmp <- load(tmp)
	rm(list=setdiff(tmp,"HS")) # keep only human

	dt <- data.table(HS[,1])
	dt[, c("chrom","start", "end") := tstrsplit(V1, "[:-]", fixed=FALSE)]
	dt[, i:=1:.N]
		
	gr <- unlist(liftOver(dtToGr(dt, metaCols="i"), from="hg19"))
		
	iToCelltypes <- as.data.table(reshape2::melt(sapply(colnames(HS)[-1], function(x) {
		 which(HS[,x]>=25) # arbitrary definition of what constitutes and active enhancer (also I don't really know what the numbers mean, as it's not described anywhere)
	}, simplify=FALSE)))
	
	iToCelltypes[, .N, by=L1]

	# find overlaps with peaks:
	o <- unique(as.data.table(findOverlaps(peaks_dt_global, gr))[, .(queryHits, i=gr$i[subjectHits])])
	
	# merge with cell type annotations:
	o <- merge(o, iToCelltypes, by.x="i", by.y="value", allow.cartesian=TRUE)
	
	# collapse multiple overlaps into one string:
	o5 <- o[order(L1),.(enhatlas_type=paste0(L1, collapse=";")), by=.(rid=peaks_dt_global$Geneid[queryHits])]
	
		o5	
}, assignToVar="o5",reload=T )
	
region.annot.list[[n]]=o5
#	return(o)
#}, buildEnvir=list(peaks_dt=peaksDt[, .(chrom, start, end, rid)]))
regionAnnotTypes <- c(regionAnnotTypes, n)

msg("Vu & Ernst (2020)")
# https://github.com/ernstlab/full_stack_ChromHMM_annotations
n <- "chromhmm"

simpleCache(paste0("region_annotations_", n, analysis.version),{
#annot <- versionedCache(n, instruction={
	# download and read Regulatory Index annotations:
	chromhmmF <- getExternalFile("hg38lift_genome_100_segments.bed.gz", "https://public.hoffman2.idre.ucla.edu/ernst/1G6UT/hg38lift_genome_100_segments.bed.gz")
	chromhmm <- fread(chromhmmF)
	setnames(chromhmm, c("chrom","start","end","chromhmm"))
	# convert to GRanges and find overlaps with peaks:
	chromhmm <- dtToGr(chromhmm, metaCols=c("chromhmm"))	
	o <- findOverlaps(peaks_dt_global, chromhmm)
	# collapse multiple overlaps into one string:
	o6 <- as.data.table(o)[,.(chromhmm_type=unique(chromhmm$chromhmm[subjectHits])),by=queryHits][order(chromhmm_type),.(annot=paste0(chromhmm_type, collapse=";")), by=.(rid=peaks_dt_global$Geneid[queryHits])]
	#return(o)
#}, buildEnvir=list(peaks_dt=peaksDt[, .(chrom, start, end, rid)]))

		o6
}, assignToVar="o6",reload=T )

region.annot.list[[n]]=o6


	regionAnnotTypes <- c(regionAnnotTypes, n)


names(region.annot.list) = regionAnnotTypes



for(j in 1:length(region.annot.list)){

region.annot.list[[j]]=region.annot.list[[j]] %>% givecolnames(., nms=c("Geneid", "annotation"))
  
}


region.annot.list
}





################################################################################
# gene peak and tfbs annotation routines
################################################################################

prepare.motif.annotation.materials=function(){
  fcat("importing gtf...")
geneAnnot=rtracklayer::import(dataset.info$gtfpath)

#get the annotation of only protein coding genes. perhaps a bit restrictive for this pipeline
selGeneAnnot <- geneAnnot %>% as.data.frame %>% dplyr::filter(type=="gene", gene_type=="protein_coding") %>% GRanges

#imported from ncnb2_code/atac/04,  doesn't work 
#selGeneAnnot <- geneAnnot[(as.character(seqnames(geneAnnot)) %in% (peaksDt[,unique(chrom)]) & geneAnnot$type=="gene" & geneAnnot$gene_type=="protein_coding"),]

## motif annotation files
  fcat("importing motifs frome HOCOMOCO v11...")
motifFile <- getExternalFile("HOCOMOCOv11_full_HUMAN_mono_jaspar_format.txt", "https://hocomoco11.autosome.org/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_HUMAN_mono_jaspar_format.txt") 
motifFileAnnot <- getExternalFile("HOCOMOCOv11_full_annotation_HUMAN_mono.tsv", "https://hocomoco11.autosome.org/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_annotation_HUMAN_mono.tsv") 
#### from 04 annotate_genes



#### from annotate motifs
  fcat("preparing motif scanning auziliar variables...")
simpleCache::simpleCache("hocomoco_v11", {
	hocomoco <- universalmotif::read_jaspar(motifFile)
	hocomoco <- XMatrixList(lapply(hocomoco, function(m) universalmotif::convert_motifs(m, class = "TFBSTools-PFMatrix")), use.names = F, type = "PFMatrixList",    matrixClass = "PFMatrix")
	names(hocomoco) <- lget(hocomoco, "name", slot=TRUE)
	return(hocomoco)
}, assignToVar="mtfSet")

hocoAnnot <- fread(motifFileAnnot)
setkey(hocoAnnot, Model)

print(mtfSet)
print(head(hocoAnnot))

# T should be TBXT, the others are to be dplyr::filtered out (both on Y chromosome)
hocoAnnot[`Transcription factor`=="T", `Transcription factor`:="TBXT"]

mtfSet <- mtfSet[hocoAnnot[names(mtfSet), `Transcription factor`] %in% geneAnnot$gene_name]

mtfNames <- hocoAnnot[names(mtfSet), `Transcription factor`] 
names(mtfNames) <- names(mtfSet)
mtfLabels <- sprintf("%s (%s)", mtfNames, names(mtfNames))
names(mtfLabels) <- names(mtfSet)

relevant.chromosomes= c(paste0("chr", 1:22), "chrX")


fcat("Variables ready.")
list(mtfSet=mtfSet, 
  mtfNames=mtfNames, 
  mtfLabels=mtfLabels, 
  selGeneAnnot=selGeneAnnot,
  relevant.chromosomes=relevant.chromosomes)
  
  
  
}




get.global.peak.annotations=function(recreate=F){
  
  library(TFBSTools)
library(motifmatchr)
  library(GenomicRanges)
  
  mtf.prep.list=prepare.motif.annotation.materials()

  list2env(mtf.prep.list, .GlobalEnv)
  ################################################################################
# arrange the peaks with the highest PC loadings into a GRanges object
################################################################################
  
  fcat("formatting peak annotation as GRanges...")
simpleCache("peaks_dt_global_annotations" %>% addversion, {
peaks_dt_global=norm.peak.vals %>% as.data.frame %>% dplyr::mutate(seqnames=Chr,start=Start, end=End) %>% dplyr::select(Geneid,seqnames, start, end) %>% filter(seqnames %in% !!relevant.chromosomes) %>% GRanges
#peaks_dt_global=norm.peak.vals %>% as.data.frame %>% dplyr::mutate(seqnames=Chr,start=Start, end=End) %>% dplyr::select(Geneid,seqnames, start, end) %>% GRanges



fcat("Annotating TSS and overlapping...")
peak.gene.annotations=list(
peaks_dt_global=peaks_dt_global ,  
closest_tss_global=geneAssignmentStrategies$closest_tss(peaks_dt_global, gene_annot=selGeneAnnot, window_size=1000000),
    
overlapping_promoters_global=geneAssignmentStrategies$promo(peaks_dt_global, gene_annot=selGeneAnnot)
  )
peak.gene.annotations
}, assignToVar="peak.gene.annotations", reload=T)

 list2env(peak.gene.annotations, .GlobalEnv)
 
simpleCache(paste0("matched_motifs_allpeaks_version", analysis.version), {
#mtfMat <- versionedCache(paste0("motif_counts_", dplyr::selected_pc), instruction={
    fcat("matching motifs...")
	if(!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly=T)){
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", update=FALSE, ask=FALSE)
	}
  
	  mtfs <- motifmatchr::matchMotifs(mtfSet, peaks_dt_global, genome=config$genome_build, out="scores") 
	
	    fcat("counting motifs...")
	  mtf_mat <- as.matrix(motifmatchr::motifCounts(mtfs))
	rownames(mtf_mat) <- peaks_dt_global$Geneid
	 
	 fcat("formatting output...")
	trimtfname= Vectorize(function(x) strsplit(x, split="_HUMAN")[[1]][1], USE.NAMES=F)
	
matrix.tfnames= trimtfname(colnames(mtf_mat))	
mtf_mat=mtf_mat %>% givecolnames(., nms=matrix.tfnames)
	
mtf_mat}, assignToVar="mtf_mat", recreate=recreate, reload=!recreate)
	#return(mtfMat)
#}, buildEnvir=list(
 # mtf_set=mtfSet, 
  #peaks_dt=peaksDt[,.(chrom,start,end,rid)],
  #genome_build=config$genome_build
#))

simpleCache("mtf_mat_gathered" %>% addversion, {
  
  unique.tfs=unique(colnames(mtf_mat))

mtf.mat.summed.list=lapply(unique.tfs, function(tf){
   fcat(tf)
  mtchs=colnames(mtf_mat)==tf
  if(sum(mtchs)==1){
    mtf_mat[,colnames(mtf_mat)==tf] %>% coerce.to.matrix
  }else{
     apply(mtf_mat[,colnames(mtf_mat)==tf], 1, sum) %>% coerce.to.matrix %>% giverownames(., nms=rownames(mtf_mat)) %>% givecolnames(.,nms= tf)
  }

}) 

mtf_mat_gathered=mtf.mat.summed.list %>% Reduce(cbind, .) %>% givecolnames(., nms=unique.tfs)
mtf_mat_gathered

}, assignToVar="mtf_mat_gathered", recreate=recreate, reload=!recreate)



simpleCache(paste_("motif_matrix_calculations_matrixid", digest::digest(mtf_mat_gathered)) %>% addversion, {
  
    list(mtf_props=apply(mtf_mat_gathered>0,2, function(x) sum(x)/length(x)),
mtf_sums=apply(mtf_mat_gathered>0,2, function(x) sum(x)),
nonmtf_sums=apply(mtf_mat_gathered==0,2, function(x) sum(x)))
  }, assignToVar="mtf.props.list", recreate=recreate, reload=!recreate)


  fcat("Done.")
list(peak.gene.annotations=peak.gene.annotations, 
  mtf_mat=mtf_mat, 
  mtf_mat_gathered=mtf_mat_gathered,
  mtf.props.list=mtf.props.list
  )


}


################################################################################
# make barplots of enhancer annotation
################################################################################

annotation.barplots=function(target.peaks, region.annot.list, charlimit=30, numtopenhancers=20, color="red", title=""){
  library(patchwork)
#target peaks is an array/data frame where each peak is a row and columns are Geneid(name of the peak) 
regions=target.peaks 

lapply(1:length(region.annot.list),function(j){



df=region.annot.list[[j]] %>% filter(Geneid %in% regions) %>% dplyr::mutate(annotation=substr(annotation, 1, charlimit)) %>% group_by(annotation) %>% summarise(counts=n()) %>% arrange(-counts) %>% head(numtopenhancers) 
ordenhancers= df %>% pull(annotation)

pl=ggplot(df, aes(x=factor(annotation, levels=ordenhancers), y=counts ))+
  geom_col(fill=color)+theme_classic()+rotatex(90)+ggtitle(paste0(title, ": ", names(region.annot.list)[[j]]))
pl
}) %>% Reduce('+',.)

}

################################################################################
################################################################################
#specific functions block
################################################################################
################################################################################

fixcolnames= function(x) {
 colnames(x)= make.names(colnames(x)) 
  x
}

editname=Vectorize(function(xp) strsplit(xp,split="_")[[1]][1], USE.NAMES=F)

splitsamplename= Vectorize(function(x)  strsplit(x, split="_R")[[1]][1], USE.NAMES=F )







plotPCA_vst <- function (object, metadata=NULL, ids=NULL, ntop = 500, assaynumber=length(assays(object))) {
    if (is.null(ids)){
    ids<-rownames(metadata) 
    }
    rv         <- rowVars(assay(object[, ids], assaynumber), useNames=T)
    select     <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca       <- prcomp(t(assay(object[, ids], assaynumber)[select, ]), center=TRUE, scale=FALSE)
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    #LUIS 20230405
    if(!is.null(metadata)){
    df         <- cbind( metadata[ids, ] , pca$x)
    }
    #Order points so extreme samples are more likely to get label
    ord        <- order(abs(rank(df$PC1)-median(df$PC1)), abs(rank(df$PC2)-median(df$PC2)))
    df         <- df[ord,]
    attr(df, "percentVar") <- data.frame(PC=seq(along=percentVar), percentVar=100*percentVar)
    return(list(df, pca))
}


sumduplicates= function(vec){
 reslt=list()
 for(k in unique(names(vec))){
 reslt[[k]]= sum(vec[names(vec)==k])
     
 }
 
 return(unlist(reslt))
}




make.gene.vector=Vectorize(function(x) x  %>% strsplit(., split=",") %>% Reduce(c, .), USE.NAMES=F)


summarisename= Vectorize(function(nm, total.print.length=30, first.words=4){
  #this function takes the first "first" words of the name, then if this exceeds total.print.length, 
  #we slice the middle section. 
  nm2=strsplit(nm, split=" ")[[1]]
  nm=nm2[1:(ifelse(first.words>length(nm2), length(nm2), first.words))] %>% paste(., collapse=" ")
  
  if(nchar(nm)>total.print.length){
    
    extend=round(total.print.length/2)
    a=substr(nm, 1, extend)
    b=substr(nm, nchar(nm)-extend, nchar(nm))
    out=paste0(a,"...", b)
    
  }else{
    out=nm
  }
  return(out)
}, USE.NAMES=F)


colortarget= Vectorize(function(xx) if(is.na(xx)){return("black")}else{
  if( xx==targetcelltype) return("red") else{if( xx=="test") return("blue") else return("grey")}}, USE.NAMES=F)

################################################################################
# genome plotting functions
################################################################################

toBin <- function(coord, offset, res) {
	round((as.numeric(coord) - as.numeric(offset)) / as.numeric(res))
}
extendAndFlipRegions <- function(regs, ext=1000) {
	i <- start(regs)>end(regs)
	tmp <- start(regs)
	start(regs)[i] <- end(regs)[i]
	end(regs)[i] <- tmp[i]

	start(regs) <- pmax(start(regs)-ext, 0)
	end(regs) <- end(regs)+ext
	regs
}
themeCovPlot <- function(p, breaks_width=10000, pData=NULL, colorlist=allcolors) {
	cols <- colorlist$sample_group
	if(!is.null(pData)) cols <- cols[pData[,unique(sample_group)]]
	p + xlab(sprintf("Linearized genome (chromosomal coordinate [mb], mean of %dbp bin)", GENOME_BIN_SIZE)) + ylab("Normalized read count") + scale_x_continuous(expand=c(0,0), breaks=function(x) seq(floor(min(x)/breaks_width)*breaks_width, ceiling(max(x)/breaks_width)*breaks_width, by=breaks_width) , labels=function(x) ifelse(order(x)==1 | order(x)==length(x), sprintf("%.1f", x/1000000),"")) + defTheme() + themeBoxed() + theme(legend.position="bottom") + scale_fill_manual(values=cols, guide="none") + scale_color_manual(values=cols, guide="none") # cowplot::theme_minimal_hgrid(line_size=0.25) + t
}
plotCoverageHM <- function(hm_data, gene_annot=NULL, sample_groups, id="x", row_annot=NA, annot_colors=NA, cell_width=NA, fn_col_agg = mean, ul=NULL, group_colors=COLOR_PALETTES$condition, no_numbers=FALSE, no_labels=FALSE) {	
	colAnnotList <- list()
	
	if(!is.null(gene_annot)) colAnnotList$genes <- gene_annot
	
	grpMeans <- t(sapply(unique(sample_groups), function(x) { apply(hm_data[sample_groups==x,,drop=F], 2, fn_col_agg, na.rm=T)
	} ))		
	if(is.null(ul)) ul <- ceiling(max(grpMeans, na.rm=T))	
	colAnnotList <- c(colAnnotList, sapply(unique(sample_groups), function(grp) {
		axisParam <- ComplexHeatmap::default_axis_param("column")
		if(no_numbers==TRUE) axisParam[["gp"]] <- grid::gpar(fontsize = 0)	
		ComplexHeatmap::anno_barplot(grpMeans[grp,], axis_param=axisParam, gp=grid::gpar(fill=group_colors[[grp]], col=NA), ylim=c(0,ul))
	}, simplify=F))
	
	gpx <- grid::gpar()
	if(no_labels==TRUE) {
		gpx <- grid::gpar(fontsize=0, col=NA)
	}
	ha <- ComplexHeatmap::HeatmapAnnotation(genes=colAnnotList$genes, psc=colAnnotList$psc, nmp=colAnnotList$nmp, nc=colAnnotList$nc, sap=colAnnotList$sap, sapm=colAnnotList$sapm, wt=colAnnotList$wt, chr17q=colAnnotList[["17q"]], chr17q1q=colAnnotList[["17q1q"]], chr17q1q_dox=colAnnotList[["17q1q_dox"]] , col=annot_colors, na_col="white", annotation_name_rot=90, annotation_name_gp=gpx)
	
	if(is.data.frame(row_annot)) {
		row_annot <- ComplexHeatmap::HeatmapAnnotation(df=row_annot, which = "row", col=annot_colors)
	}
	
	rowGaps <- cumsum(table(sample_groups)[unique(sample_groups)])
	
	p <- ComplexHeatmap::pheatmap(as.matrix(hm_data), cellwidth=cell_width, cellheight=10, gaps_row=rowGaps, show_colnames=F, column_title=id, name="Normalized read count", annotation_colors=annot_colors, right_annotation=row_annot,  top_annotation = ha,  border_color=NA, cluster_cols=F, cluster_rows=F, col=HEATMAP_COLORS)
	# optional -- fix width: 	cellwidth=8  
	return(p)
}


na2nothing= function(x) ifelse(is.na(x), "", paste0("-", x))
testcompress=function(x) ifelse(grepl("test",x), "test", x)



geneAssignmentStrategies <- list(
	# assign peak to closest gene by peak-to-TSS distance (within max. distance):
	closest_tss = function(gr, gene_annot=selGeneAnnot, window_size=MAX_DIST_GENE) {	
		tss <- promoters(gene_annot, 0, 1)
		assoc <- as.data.table(GenomicRanges::distanceToNearest(gr, tss, ignore.strand=T))
		assoc <- assoc[, .(gr_index=queryHits, gene_symbol=tss@elementMetadata$gene_name[subjectHits], gene_dist=distance)]
			assoc<- assoc %>% dplyr::mutate(Geneid=gr[assoc$gr_index,]$Geneid)
		return(assoc[gene_dist<=window_size,])	
	},
	# assign peaks to all genes with "significant correlation" (within max. distance):
	sigcor = function(gr, gene_annot=selGeneAnnot, pos_cor_only=TRUE, pthresh=GENE_ASSOC_SIGCOR_THRESH, window_size=MAX_DIST_GENE, shuffle_cor=dShuffleCors, merge_cor=dMergeCor) {
		fn <- abs
		if(pos_cor_only) fn <- identity
	
		shuffle_cor <- shuffle_cor[rid%in%gr$rid & gene_symbol%in%gene_annot@elementMetadata$gene_name,]
		merge_cor <- merge_cor[rid%in%gr$rid & gene_symbol%in%gene_annot@elementMetadata$gene_name,]
	
		minCor <- as.numeric(round(shuffle_cor[is.finite(pearson),quantile(fn(pearson), 1-pthresh, na.rm=T)],1))
		msg(minCor)
		
		return(merge_cor[fn(pearson)>=minCor, .(gr_index=structure(1:length(gr), names=gr$rid)[rid], gene_symbol, gene_dist)])
	},
	# assign peaks to genes with overlapping promoters:
	promo = function(gr, gene_annot=selGeneAnnot, prox_width=5000) {
		promos <- GenomicRanges::promoters(gene_annot)
		assoc <- as.data.table(GenomicRanges::findOverlaps(gr, promos, ignore.strand=T))
		print(head(assoc))
		assoc[, gene_dist:= GenomicRanges::distance(gr[queryHits], promos[subjectHits])]
		assoc <- assoc[, .(gr_index=queryHits, gene_symbol=promos@elementMetadata$gene_name[subjectHits], gene_dist)]
		assoc<- assoc %>% dplyr::mutate(Geneid=gr[assoc$gr_index,]$Geneid)
		
		return(assoc)	
	},
	# assign peaks to genes with "significant correlation" or with overlapping promoter:
	promo_or_sigcor = function(gr, gene_annot=selGeneAnnot, pos_cor_only=TRUE, pthresh=GENE_ASSOC_SIGCOR_THRESH, window_size=MAX_DIST_GENE, shuffle_cor=dShuffleCors, merge_cor=dMergeCor) {
		sigRes <- geneAssignmentStrategies$sigcor(gr, gene_annot=gene_annot, pos_cor_only=pos_cor_only, pthresh=pthresh, window_size=window_size, shuffle_cor=shuffle_cor, merge_cor=merge_cor)
		promoRes <- geneAssignmentStrategies$promo(gr, gene_annot=gene_annot)
		
		return(unique(rbind(sigRes, promoRes)[,.(gene_dist=min(gene_dist)),by=.(gr_index, gene_symbol)]))
	}
)


annotate.close.gene.vector= function(targetpeaks){
 full_join( overlapping_promoters_global %>% dplyr::mutate(overlapping.promoter=gene_symbol)  %>% dplyr::filter(Geneid %in% targetpeaks)  %>% dplyr::select(Geneid, overlapping.promoter) , closest_tss_global %>% dplyr::mutate(closest.tss=gene_symbol, gene.dist=gene_dist) %>% dplyr::filter(Geneid %in% targetpeaks) %>% dplyr::select(Geneid, closest.tss, gene.dist), by="Geneid")
  #%>% arrange(factor(Geneid, levels=targetpeaks))
   
}


## rough merge of tables, contains duplicate intervals in geneid.
annotate.close.genes=function(targetpeaks){
  
  
  
  #case vector
if(!is.list(targetpeaks) && is.vector(targetpeaks)){
  gene.annotated.peaks=annotate.close.gene.vector(targetpeaks)
}
 #case list 
 if(is.list(targetpeaks)){
 gene.annotated.peaks= lapply(targetpeaks, function(x) annotate.close.gene.vector(x)) %>% givename(., names(targetpeaks))
}
   
  if(is.data.frame(targetpeaks)){
gene.annotated.peaks=targetpeaks %>% left_join(., overlapping_promoters_global %>% dplyr::mutate(overlapping.promoter=gene_symbol) %>% dplyr::select(Geneid, overlapping.promoter), by="Geneid") %>% left_join(., closest_tss_global %>% dplyr::mutate(closest.tss=gene_symbol, gene.dist=gene_dist) %>% dplyr::select(Geneid, closest.tss, gene.dist), by="Geneid")

  }
  gene.annotated.peaks
}



 getCoverageInInterval <- function(reads, regs, res=50, nthreads=10) {
  	library(Rsubread)
  
    # define the bins in the regions of interest:
  	bins <- rblapply(1:length(regs), function(i) {
  	  #fcat("region",i)
  		gr <- regs[i]
  		s <- round(start(gr)/res)*res
  		e <- round(end(gr)/res)*res
  		
  		#fcat("start", s, "end", e)
  		if(s > e) {
  			x <- s; s <- e; e <- x
  		}
  		if(abs(s-e)<res){
  		  bins=NULL}else{
  		    

  		bins <- data.table(Chr=as.character(seqnames(gr)), Start=seq(s,e-res,by=res), End=seq(s+res,e,by=res), Strand=as.character(strand(gr)), gene_name=gr$gene_name, start.original=gr$start.original, end.original=gr$end.original)
  		bins[,bin_num:=1:.N]
  		bins[,locusID:=paste_(gene_name,i,bin_num)]
  		bins[,geneID:=paste_(i,bin_num)]
  		
  		}
  		bins
  	}, "reg_id")
  	setkey(bins, geneID)
  	
  	# then use featureCounts to calcuate the overlapping read counts:
  	dtCov <- Rsubread::featureCounts(files=reads, annot.ext=as.data.frame(bins), allowMultiOverlap=T, strandSpecific=0, nthreads=nthreads, isPairedEnd=T)
  	#colnames(dtCov$counts) <- names(reads)		
  	dtCov <- melt(as.data.table(dtCov$counts,keep.rownames="bin_id"), id.vars="bin_id", variable.name="library_name", value.name="count")
  	dtCov <- cbind(dtCov, bins[dtCov$bin_id,.(reg_id,bin_num,bin_start=Start,bin_end=End, start.original, end.original)]) %>% as.data.frame %>% dplyr::mutate(Experiment=editname(as.character(library_name)))
  	dtCov
 }
 
 get.all.embeddings= function(pc) pca.data %>% dplyr::arrange(!!sym(pc)) %>% dplyr::select(!!sym(pc), celltype)
get.sample.peakvals= function(peakid, samples=NULL){
  
  if(!is.null(samples)){ mt=norm.peak.vals[ peakid, samples]}else{ mt=norm.peak.vals[peakid, !colnames(norm.peak.vals) %in% c("Geneid", "Chr", "Start", "End", "Strand", "Length")]}
  
  mt %>% t
  
}

get.samples.from.class=function(x, targetlist) pca.data %>% filter(!!sym(targetlist[[x]]$classvar) %in%  targetlist[[x]]$targets) %>% pull(Experiment)

get.sample.class.metadata=function(x, targetlist) pca.data %>% filter(!!sym(targetlist[[x]]$classvar) %in%  targetlist[[x]]$targets) %>% select(all_of(c(dataset.info$reagents, dataset.info$treatment_id_cols)))
  


################################################################################
# function to test tfbs enrichment based on a fisher and hypergeometric test
################################################################################

test.tfbs.enrichment=function(targetpeaks, mtf_mat_gathered, simple.output=T, min.frequency.in.target=0, target.name="target", mtf.props.list=NULL, filtered=T, background.peaks=NULL){
fcat("loading precalculated features of motif matrix if any...")
  
  if(is.null(mtf.props.list)){

    allpeakannots=get.global.peak.annotations()
    list2Env(allpeakannotations, .GlobalEnv)

  }
  
  
  mtfsums.raw=colSums(mtf_mat_gathered)/nrow(mtf_mat_gathered)
  
 if(is.null(background.peaks)){
   
  background.peaks= rownames(mtf_mat_gathered)
 }else{
   
  background.peaks= intersect(rownames(mtf_mat_gathered), background.peaks) 
 }
  
  gc()
  ##mtf_mat_gathered is a matrix where all tfbs have been counted in all the peaks, and all the sites corresponding to the same tfbs have been summed up
  if(length(targetpeaks)<=1){
   fcat("warning: lenngth of target peaks is not more than 1. returning NULL")
    return(NULL) 
    
  }else{
    
    
fcat("Calculating initial frequencies...")
  obsmat=mtf_mat_gathered[intersect(rownames(mtf_mat_gathered), targetpeaks),]
  restmat=mtf_mat_gathered[setdiff(background.peaks, targetpeaks),]

targetsums=apply(obsmat>0,2, function(x) sum(x))
targetsize=nrow(obsmat)
targetsums.raw=colSums(obsmat)/targetsize
targetprops=apply(obsmat>0,2, function(x) sum(x)/length(x))


restsums=apply(restmat>0,2, function(x) sum(x))

restsize=nrow(restmat)
restsums.raw=colSums(restmat)/restsize
restprops=apply(restmat>0,2, function(x) sum(x)/length(x))

fisher.results=matrix(nrow=2, ncol=ncol(mtf_mat_gathered), data=NA)
colnames(fisher.results)=colnames(mtf_mat_gathered)
rownames(fisher.results)=c("pvalue", "odds.ratio")
contingency.mats=list()
fcat("Looping through TFs...")

if(!is.null(background.peaks)){
#background matrix (without the target peaks) for fisher test

#calculations combining target and rest matrices for hypergeometric
combined.mtfsums=targetsums+restsums
combined.mat.size=restsize+targetsize
combined.nonmtfsums=combined.mat.size-(restsums+targetsums)
combined.props=combined.mtfsums/(combined.mat.size)
}else{


combined.mtfsums=allpeakannots$mtf.props.list$mtf_sums
combined.mat.size=allpeakannots$mtf.props.list$mtf_size
combined.nonmtfsums=allpeakannots$mtf.props.list$nonmtf_sums
combined.props=allpeakannots$mtf.props.list$mtf_props
}


final.list=lapply(colnames(mtf_mat_gathered), function(j){
  contmat=matrix(nrow=2, ncol=2, data=NA) %>% givecolnames(., nms=c("mtf.present","mtf.absent")) %>% giverownames(., nms=c("target","rest"))
  

contmat["target","mtf.present"]=targetsums[j]
contmat["target","mtf.absent"]=targetsize-targetsums[j]
contmat["rest","mtf.present"]=restsums[j]
contmat["rest","mtf.absent"]=restsize-restsums[j]

ft=fisher.test(contmat)

list(fisher.pvalue=ft$p.value,
fisher.odds.ratio=ft$estimate,
contingency.matrix=contmat,
  global.tfbs.probability=allpeakannots$mtf.props.list$mtf_props[j],
  sample.tfbs.probability=targetprops[j],
  bg.tfbs.probability=restprops[j],
  target.rawsums=targetsums.raw[j],
  rest.rawsums=restsums.raw[j],
  all.rawsums=mtfsums.raw[j],
  fold.enrichment.allpeaks=targetprops[j]/allpeakannots$mtf.props.list$mtf_props[j],
  fold.enrichment.bg=targetprops[j]/restprops[j],
  hypergeom.pvalue.lower.allpeaks=phyper(targetsums[j],allpeakannots$mtf.props.list$mtf_sums[j],allpeakannots$mtf.props.list$nonmtf_sums[j], targetsize ),
  hypergeom.pvalue.higher.allpeaks=1-phyper(targetsums[j],allpeakannots$mtf.props.list$mtf_sums[j],allpeakannots$mtf.props.list$nonmtf_sums[j], targetsize ),
  hypergeom.pvalue.lower.bg=phyper(targetsums[j],combined.mtfsums[j],combined.nonmtfsums[j], targetsize ),
  hypergeom.pvalue.higher.bg= 1-phyper(targetsums[j],combined.mtfsums[j],combined.nonmtfsums[j], targetsize )
  )

}) %>% givename(., colnames(mtf_mat_gathered))

fcat("Preparing output...")
if(simple.output){
fe.allpeaks=lapply(final.list, function(x) x$fold.enrichment.allpeaks) %>% Reduce(c, .)
fe.bg=lapply(final.list, function(x) x$fold.enrichment.bg) %>% Reduce(c, .)
sums.target=lapply(final.list, function(x) x$fold.enrichment.bg) %>% Reduce(c, .)
phighers.all=lapply(final.list, function(x) x$hypergeom.pvalue.higher.allpeaks) %>% Reduce(c, .)
phighers.bg=lapply(final.list, function(x) x$hypergeom.pvalue.higher.bg) %>% Reduce(c, .)
sampleprob=lapply(final.list, function(x) x$sample.tfbs.probability) %>% Reduce(c, .)
rawsums.target=lapply(final.list, function(x) x$target.rawsums) %>% Reduce(c, .)
rawsums.rest=lapply(final.list, function(x) x$rest.rawsums) %>% Reduce(c, .)
rawsums.all=lapply(final.list, function(x) x$all.rawsums) %>% Reduce(c, .)
bgprob=lapply(final.list, function(x) x$bg.tfbs.probability) %>% Reduce(c, .)
globalprob=lapply(final.list, function(x) x$global.tfbs.probability) %>% Reduce(c, .)
fisher.pvalues=lapply(final.list, function(x) x$fisher.pvalue) %>% Reduce(c, .)
fisher.odds.ratios=lapply(final.list, function(x) x$fisher.odds.ratio) %>% Reduce(c, .)

tfdf=list(fold.enrichment=fe.allpeaks, fold.enrichment.bg=fe.bg, rawsums.all=rawsums.all, rawsums.target=rawsums.target,rawsums.rest=rawsums.rest, hypergeom.higher=phighers.all, hypergeom.higher.bg=phighers.bg, frequency.in.target=sampleprob, frequency.in.bg= bgprob,fisher.odds.ratio.bg=fisher.odds.ratios, frequency.global=globalprob, fisher.pvalue.contrast=fisher.pvalues) %>% as.data.frame 
tfdf=tfdf %>% dplyr::mutate(p.adj=p.adjust(tfdf$hypergeom.higher, method="BH" ), p.adj.bg=p.adjust(tfdf$hypergeom.higher.bg, method="BH" ),relevance=frequency.in.target*log2(fold.enrichment), target=!!target.name) %>% names2col(., "gene_symbol")

if(filtered){
atfdf=tfdf %>% dplyr::filter(p.adj<=0.05, abs(log2(fold.enrichment))>=1,frequency.in.target>=min.frequency.in.target) %>% arrange(-fold.enrichment)
}else{
atfdf=tfdf %>% arrange(-fold.enrichment) 
}
gc()
return(atfdf)

}
  }
 
}

################################################################################
# normalise columns in matrix
################################################################################
normalise.columns <- function(x) {
  y <- lapply(1:ncol(x), function(i) {
    x[, i, drop = FALSE] / max(x[, i])
  }) %>% Reduce(cbind, .)
  return(y)
}



################################################################################
# function as above but 
################################################################################



################################################################################
# function to do the above but for a list
################################################################################


test.tfbs.enrichment.list=function(targetpeaks.list, mtf_mat, simple.output=T, min.frequency.in.target=0, mtf.props.list=NULL){
  
  fcat("loading precalculated features of motif matrix if any...")
  
  if(is.null(mtf.props.list)){
  simpleCache(paste_("motif_matrix_calculations_matrixid", digest::digest(mtf_mat)) %>% addversion, {
  
    list(mtf_props=apply(mtf_mat>0,2, function(x) sum(x)/length(x)),
mtf_sums=apply(mtf_mat>0,2, function(x) sum(x)),
nonmtf_sums=apply(mtf_mat==0,2, function(x) sum(x)))
  }, assignToVar="mtf.props.list", reload=T)
  }
  
  if(is.null(names(targetpeaks.list))) {
    names(targetpeaks.list)= paste0("peakset_", 1:length(targetpeaks.list))
  }
    
  finalmat=lapply(1:length(targetpeaks.list), function(xxx){
    
  
    rr=test.tfbs.enrichment(targetpeaks.list[[xxx]], mtf_mat_gathered=mtf_mat, simple.output=T, min.frequency.in.target=0, target.name=names(targetpeaks.list)[xxx], mtf.props.list=mtf.props.list)
  }) %>% Reduce(rbind, .)
  
  finalmat
  
}




################################################################################
# misc functions to collect information about genes
################################################################################



get.locus.genes=Vectorize(function(x, type="tss", return.distance=F){
  #type can be "tss" or "promoter"
  if(type=="tss"){
  mat=closest_tss_global
  }else{
  mat=overlapping_promoters_global
    }
  
  ops= mat[Geneid==x,]
  if(nrow(ops)!=0){
    
    if(return.distance && type=="tss"){
       return(paste(ops %>% pull(gene_dist), collapse=","))
    }else{
    return(paste(ops %>% pull(gene_symbol), collapse=","))
    }
      }else{
     return(NA)
  }
}
  , USE.NAMES=F)


get.closest.tss=function(peaks, limit) closest_tss_global %>% dplyr::filter(gene_dist<=limit ,  Geneid %in% peaks ) %>% arrange(gene_dist) %>% pull(gene_symbol)


get.closest.tss.df=function(peaks, limit) closest_tss_global %>% dplyr::filter(gene_dist<=limit ,  Geneid %in% peaks ) %>% arrange(gene_dist) 
