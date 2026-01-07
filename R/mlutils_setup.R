################################################################################
# adding a path for installing libraries
################################################################################

dir.create( "~/mnt_out/containerlibrary/", recursive=T)
.libPaths(c(.libPaths(), "~/mnt_out/containerlibrary/"))

################################################################################
# installing and loading packages decision tree
################################################################################

install_and_load <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    # Try installing with BiocManager first
    tryCatch({
      BiocManager::install(pkg, dependencies = TRUE, update = FALSE, lib="~/mnt_out/containerlibrary/")
    }, error = function(e) {
      message(paste("BiocManager could not install", pkg, "- trying install.packages instead."))
      install.packages(pkg, dependencies = TRUE, lib="~/mnt_out/containerlibrary/")
    })
    
    # Load the package after installation
    library(pkg, character.only = TRUE)
  }
}


################################################################################
# load libraries and packages
################################################################################

packages <- c(
    "GSVA", #calculating signature scores
  "ggplot2",
  "viridis",
  "patchwork",
  "TFBSTools", 
  "ggwordcloud", #making word clouds
  "simpleCache", #caching intermediate data for reproducibility
  "ggVennDiagram",
  "ComplexHeatmap",
  "ggrepel", 
  "ggrastr", 
  "patchwork",
  "dplyr", 
  "data.table",
  "ggwordcloud", 
  "simpleCache",
  "tidyr",
  "DESeq2", 
  "pROC",
  "randomForest",
  "motifmatchr",
  "ggVennDiagram", 
  "ggsignif"#,
  #  "qusage"# problemativ at 4.1.3
)


# Apply the function to each package
sapply(packages, install_and_load)

#set cache dir

setCacheDir("~/mnt_out/rcache")


#requirement for hypeR installation 
#devtools::install_version("msigdbr", version = "7.4.1")
install_and_load("hypeR")


################################################################################
#formatting functions
################################################################################

remove.pipes=function(x) gsub("\\|","-",x)

givenames=
  function(x, nam) {
    names(x)= nam
    x
  } 

givecolnames=function(x, ind, nms){
  colnames(x)[ind]= nms
  x
  
}


givename=
  function(x, nam) {
    names(x)= nam
    x
  } 


################################################################################
# updating configs
################################################################################
          
update_config <- function(config_a_path, config_b_path, output_path = "./config_modified.yaml") {
  # Read files
  config_a <- readLines(config_a_path)
  config_b <- readLines(config_b_path)

  # Parse A into named list
  parse_line <- function(line) {trimws
    parts <- strsplit(line, ":", fixed = TRUE)[[1]]
    if (length(parts) >= 2) {
      key <- trimws(parts[1])
      value <- paste(parts[-1], collapse = ":") |> trimws()
      return(setNames(list(value), key))
    }
    NULL
  }

  # Build a named vector of overrides
  overrides <- lapply(config_a, parse_line)
  overrides <- overrides[!sapply(overrides, is.null)]
  overrides <- do.call(c, overrides)

  # Update config B lines
  updated_b <- sapply(config_b, function(line) {
    parts <- strsplit(line, ":", fixed = TRUE)[[1]]
    if (length(parts) >= 2) {
      key <- trimws(parts[1])
      if (key %in% names(overrides)) {
        return(paste0(key, ": ", overrides[[key]]))
      }
    }
    return(line)
  })

  # Write to file
  writeLines(updated_b, output_path)
}


update_config2 <- function(config_a_path, config_b_path, output_path = "./config_modified.yaml") {
  # Read files
  config_a <- readLines(config_a_path)
  config_b <- readLines(config_b_path)
  
  # Helper: extract YAML key from a line (before ":")
  extract_key <- function(line) {
    parts <- strsplit(line, ":", fixed = TRUE)[[1]]
    if (length(parts) >= 2) {
      return(trimws(parts[1]))
    }
    return(NA)
  }
  
  # Get keys from config B
  keys_b <- sapply(config_b, extract_key)
  keys_b <- keys_b[!is.na(keys_b) & nzchar(keys_b)]
  
  # Filter out from A any line whose key appears in B
  filtered_a <- config_a[!(sapply(config_a, extract_key) %in% keys_b)]
  
  # Append B to A
  combined <- c(filtered_a, config_b)
  
  # Write output
  writeLines(combined, output_path)
}


################################################################################
# ACCESSORY FUNCTIONS. quick access functions to important info
################################################################################

get.experiments=function(condition, meta) meta %>% filter(Condition==!!condition) %>% pull(Experiment)

get.dscategory= Vectorize(function(n) dataset.info$dscategory[n], USE.NAMES=F)

get.dscategory.fromname= Vectorize(function(nm) dataset.info$dscategory[dataset.info$dataset_name==nm], USE.NAMES=F)
################################################################################
# assembling workspace variables
################################################################################

prep.dir.for.out.root= function(x) gsub("/home/rstudio/mnt_out|~/mnt_out", config$out_root_host, x) 

prep.dir.external= function(x) gsub("/home/rstudio/mnt_out|~/mnt_out", config$out_root_host, x) %>%  gsub("/home/rstudio|~", config$project_root_host, .)


  prep.dir.for.in.root <- function(x) {
  gsub(config$out_root_host, config$out_root, x, fixed = TRUE)
}

prep.dir.for.project.root= function(x) gsub("~", config$project_root_host, x)



################################################################################
# constructing nfcore paths from within the pipeline
################################################################################

make.countpath=function(rootdir) paste0(rootdir, "bwa/merged_library/macs2/broad_peak/consensus/consensus_peaks.mLb.clN.featureCounts.txt")

make.bigwigpath=function(rootdir) paste0(rootdir, "bwa/merged_library/bigwig/")
make.profilepath= function(x) paste0(x, "/bwa/merged_library/deeptools/plotprofile/")

make.bampath=function(rootdir) paste0(rootdir, "/bwa/merged_library/")

make.allbampath=function(rootdir) paste0(rootdir, "bwa/merged_library/*.bam")

make.bampath.rep=function(rootdir) paste0(rootdir, "/bwa/merged_replicate/")


clean.extension.bam=function(x) gsub(".mLb.clN.sorted.bam", "", x)

get.nfcore.root=function(datasetnumber, dataset.info) dataset.info$dataset_paths_nfcore[datasetnumber] 

################################################################################
# get and assemble counts paths
################################################################################

get.counts.path.self=function(dataset, dataset.info){
 
if(is.na(dataset.info$dataset_paths_counts_self[dataset])||dataset.info$dataset_paths_counts_self[dataset]=="NA"){

  paste0(dataset.info$dataset_paths_nfcore[dataset], "bwa/merged_library/macs2/broad_peak/consensus/consensus_peaks.mLb.clN.featureCounts.txt")
  
}else{
  dataset.info$dataset_paths_counts_self[dataset]
}
}

get.counts.path.ref=function(dataset, dataset.info){
  
 
  
  if(is.na(dataset.info$dataset_paths_counts_ref[dataset] %>% as.logical)||dataset.info$dataset_paths_counts_self[dataset]=="NA"){
    #if it is not explicitly declared and it is a reference
     if(dataset.info$is_reference[dataset] & sum(dataset.info$is_reference)==1){
 out=get.counts.path.self(dataset, dataset.info)
}else{
    
  out=paste0("/home/rstudio/mnt_out/ml2cell_", dataset.info$dataset_name[dataset], "_mapto_", dataset.info$dataset_name[dataset.info$is_reference]) %>% paste0(., "/allmergedcounts_withreference.txt")
}}else{
  
  out=dataset.info$dataset_paths_counts_ref[dataset] %>% paste0(., "/allmergedcounts_withreference.txt")
}
  
out
}


################################################################################
# assemble the nextflow ml2cell command
################################################################################

assemble.ml2cell.nextflow=function(test.datasetnumber, reference.datasetnumber, dataset.info){
  
test.dataset=dataset.info$dataset_name[test.datasetnumber]
reference.dataset=dataset.info$dataset_name[reference.datasetnumber]

bamlocation= make.allbampath(get.nfcore.root(test.datasetnumber, dataset.info))
  annotation= make.countpath(get.nfcore.root(reference.datasetnumber, dataset.info))
  refdata=annotation
  outdir= get.mappedcounts.path(test.dataset, reference.dataset, dataset.info)
  
  cmd <- sprintf(
    "nextflow run /home/rstudio/ml2cell_code/nextflow/channeltest.nf --bamlocation \"%s\" -w %s --annotation %s --endparameter '\\-p' -with-report --refdata %s --outdir %s",
    bamlocation,
    dataset.info$nextflow_workdir,
    annotation,
    refdata,
    outdir
  )
  return(cmd)
}



################################################################################
# function to get colors when present or not
################################################################################

default.colors=list()
  default.colors[[2]]=c(`0` = "darkgoldenrod1", `1` = "darkgoldenrod",`FALSE` = "darkgoldenrod1", `TRUE` = "darkgoldenrod")

getcolors <- function(df, colname, default = 2, tf.names = TRUE) {
  # Ensure the global variable allcolors exists
  if (!exists("allcolors", inherits = TRUE)) {
    stop("Global variable 'allcolors' does not exist.")
  }

  # If colname is in allcolors, return the associated color vector
  if (colname %in% names(allcolors)) {
    return(allcolors[[colname]])
  }

  # If colname is a column in the dataframe
  if (colname %in% colnames(df)) {
    values <- unique(df[[colname]]) %>% removenas
    if(!is.null(default.colors[[length(values)]])){
      colors <- default.colors[[length(values)]]
    }else{
    colors <- grDevices::rainbow(length(values))
    }
    names(colors) <- values
    return(colors)
  }

  # If colname is not in allcolors and not in df, handle default
  if (default == 2) {
    colors <- grDevices::rainbow(2)
    if (tf.names) {
      names(colors) <- c("TRUE", "FALSE")
    }
    return(colors)
  }

  # Otherwise, return random unnamed colors of length `default`
  return(grDevices::rainbow(default))
}


################################################################################
#function to import all variables from the config to make them arrays when they are separated by commas
################################################################################


config.to.vectors=function(config){
dataset.info=lapply(names(config), function(x){
  #fcat(x)
  if(!is.null(config[[x]]) && is.character(config[[x]])){
  strsplit(config[[x]], split=",")[[1]]
  }else{
   config[[x]]
  }
}) %>% givename(., names(config))
  
dataset.info$is_reference=as.logical(dataset.info$is_reference) 
dataset.info$is_external=as.logical(dataset.info$is_external) 
dataset.info$is_test=as.logical(dataset.info$is_test) 
dataset.info$filter.frip=as.logical(dataset.info$filter.frip) 
dataset.info$filter.peaks=as.logical(dataset.info$filter.peaks)
dataset.info$frip_threshold=as.numeric(dataset.info$frip_threshold)
dataset.info$dataset_paths_counts_ref= lapply(1:length(dataset.info$dataset_paths_counts_ref), function(x) get.counts.path.ref(x, dataset.info)) %>% Reduce(c, .) 
dataset.info$dataset_paths_counts_self= lapply(1:length(dataset.info$dataset_paths_counts_self), function(x)get.counts.path.self(x, dataset.info)) %>% Reduce(c, .) 

dataset.info
}                               
#########################################


################################################################################
# incorporate csv data into an updated config file
################################################################################
fcat("Reading dataset_info.csv and incorporating data into config...")

# Load required packages
library(readr)

# Read the CSV
csv_data <- read_csv("~/dataset_info_example.csv")

# Create an output file
yaml_path <- "~/dataset.config.yaml"

if(!file.exists(yaml_path)){
fileConn <- file(yaml_path, open="a")
}else{
fileConn <- file(yaml_path, open="w")
#writeLines("", fileConn)
close(fileConn)
fileConn <- file(yaml_path, open="a")
}
# Loop through each column
for (cc in 1:length(names(csv_data))) {
 colname=names(csv_data)[cc]
  # Collapse all rows in the column into a comma-separated string
  values <- paste(csv_data[[colname]], collapse = ",")
  
  # Write the line to the YAML file
  writeLines(paste0(colname, ": ", values), con=fileConn )
}

# Close the file connection
close(fileConn)

#perform the updating

update_config2("~/config.yaml", "~/dataset.config.yaml", output_path = "~/config.yaml")

fcat("Dataset info has been prepared and written onto config.yaml.\n")


################################################################################
#actually import the config file and create vectors when needed, so all other 
#functions have access to the config
################################################################################
fcat("Importing config information...")

config <- yaml::read_yaml(file = '~/config.yaml')

config$is_reference=as.character(csv_data$dataset_type=="reference") %>% paste(., collapse=",")
config$is_test=as.character(csv_data$dataset_type=="test")%>% paste(., collapse=",")
config$is_external=as.character(csv_data$dataset_type=="external") %>% paste(., collapse=",")
config$dscategory=csv_data$dataset_type %>% paste(., collapse=",")

## transform config variables into vectors based on comma separation  
dataset.info=config.to.vectors(config)



################################################################################
# colors (some colors assigned programmatically)
################################################################################

write.color.file=function(colorlist, filename="colors"){
  
  path=paste0(config$project_root, "/", filename, ".csv")

library(purrr)
library(readr)
  
# convert to tidy data frame
df <- map_dfr(names(colorlist), function(varname) {
  tibble(
    variable = varname,
    value    = names(colorlist[[varname]]),
    color    = unname(colorlist[[varname]])
  )
})


# write to CSV
write_csv(df, path)

}


read.color.file <- function(filename="colors") {
  path=paste0(config$project_root, "/", filename, ".csv")
  
  df <- read_csv(path, show_col_types = FALSE)
  
  # check required columns
  stopifnot(all(c("variable", "value", "color") %in% names(df)))
  
  # split by variable and build named vectors
  nms=unique(df$variable)
  
  ilist=lapply(nms, function(x){
    mdf=df %>% dplyr::filter(variable==!!x)
  c(mdf$color) %>% givename(., mdf$value)
    }) %>% givename(., nms)
  
    #lapply(l, function(x) 
    #map(~ set_names(.x$color, .x$value)) %>%
    #set_names(map_chr(group_keys(df, variable), 1))
  ilist
}


allcolors=read.color.file()

################################################################################
#add current version to a character string
################################################################################

addversion= function(x) paste_(x, "version", analysis.version)

#this comes in handy when we have to recall static analyses from previous versions
addversion2=function(x) paste_(x, "version", analysis.version2)

################################################################################
# Labeling arrays and matrices
################################################################################

################################################################################
#force a named vector into a matrix
################################################################################
coerce.vector.to.matrix=function(vec)  matrix(vec, ncol = 1, dimnames = list(names(vec), NULL))

################################################################################
# global image parameters
################################################################################

w=250
pw=2.5
rr=600

################################################################################
#timestamp png or pdf plots
################################################################################

timestamp= function(nm, date.use=NULL)  ifelse(is.null(date.use), paste0(Sys.Date(),"_", nm), paste0(date.use,"_", nm))

################################################################################
# add prefix and suffix to a name vector
################################################################################

addprefix= (function(cc, prefix="C")  paste0(prefix, cc)) %>% Vectorize


################################################################################
# modifying the plot labels in ggplots in  a simple fashion
################################################################################

#resizing overall text in the plot
resizetext=function(x) theme_classic(base_size = x)

rotatex= function(x) theme(axis.text.x = element_text(angle = x, hjust=1, vjust=0.5))


################################################################################
# Trimming column names from a feature counts table in which the names have a long bam file suffix
################################################################################

cleannames <- function(a, sep="_"){
  ncn<- a %>%  lapply(., function(x) {strsplit(x, split=sep)[[1]][1]})  %>% unlist
  ncn
}


cleancolnames <- function(a){
  ncn<- a %>% colnames %>% lapply(., function(x) {strsplit(x, split="_")[[1]][1]})  %>% unlist
  colnames(a)<- ncn
  a
}


################################################################################
#FUNCTIONS TO EDIT METADATA FIELDS
################################################################################


assignhealth<- Vectorize(function(x){
  out=0
  if(any(grepl('Donor|4983|6792|CB102|BM-CD34', x))) {return("Healthy"); out=1} 
  if(any(grepl('SU', x))) {return("AML"); out=1} 
  if(out==0){
    return(x) 
  }
}, USE.NAMES=F)

assignhealth2<- Vectorize(function(x){
  out=0
  if(any(grepl('leukem', x))) {return("AML"); out=1} 
  else {return("Healthy"); out=1} 
  if(out==0){
    return(x) 
  }
}, USE.NAMES=F)


pcadf=function(sc){
  cbind(metadata(sc), Embeddings(sc, reduction="pca"))
}

getdate= function() format(Sys.time(), "%Y%b%d")


setlevels=function(so, varr, levls){
  Idents(so)=varr
  so@active.ident <- factor(x = so@active.ident, levels = levls)
  so}


col2names= function(dff, coll){
  rownames(dff)= dff[, coll]
  dff
}


names2col=function(dff, coln="gene") {
  dff[, coln]= rownames(dff) 
  dff
}

paste_=function(...) paste(..., sep="_")


showpalette=function(p, pname="generic"){
  sca=10
  rr=600
  w=250
  tpng(paste_("palette", pname), res=rr, he=sca*w, wi=sca*w*2)
  
  pl=ggplot(1:length(p) %>% as.data.frame %>% dplyr::mutate(h=1, colors=p, name=names(p)), aes(x=factor(.), y=factor(h), label=colors))+
    geom_col(aes(fill=factor(.)))+
    scale_fill_manual(values=namenums(p))+
    theme_classic()+
    coord_cartesian(ylim=c(0,1))+
    geom_text(y=0.5, color="black", angle=90)+
    #NoAxes()+
    NoLegend()+
    ggtitle(pname)
  
  if(is.null(names(p))){
    print(pl )
  }else{
    print(pl+geom_text(inherit.aes=T,aes(label=name), y=1.1, color="black", angle=90))  
  }
  dev.off()
  
  return(pl+geom_text(inherit.aes=T,aes(label=name), y=1.1, color="black", angle=90))
}



###########################
#Below: imports from scutils_lean
###########################



##########################################################################################################################
#general functions
##########################################################################################################################
###########################################################################################################################
cat("adding general functions...\n")
#permutate an array
shuffle=function(x) sample(x, length(x), replace=F)

#prepare gene DE table by removing NAs, filtering by cutoff p-value
filterorderDE<-function(x, thresh=0.05, lfthresh=1, netchange=FALSE, top=NA){
  
  if(netchange==TRUE){
    x=x[is.na(x[,pvaluecol])==FALSE & x[,pvaluecol ]<thresh & abs(x[,lfccol])>lfthresh,]
    ord=order(abs(x[,lfccol]), decreasing=T)
  }else{
    x=x[is.na(x[,pvaluecol])==FALSE & x[,pvaluecol ]<thresh & x[,lfccol]>lfthresh,]
    ord=order(x[, lfccol], decreasing=TRUE)
  }
  if(!is.na(top)){
    x[ord[1:top], ]
  }else{
    x[ord, ]
  }
}


reddf=function(sc, reduc="pca"){
  cbind(metadata(sc), Embeddings(sc, reduction=reduc))
}

harmonydf=function(sc){
  cbind(metadata(sc), Embeddings(sc, reduction="harmony"))
}

umapdf=function(sc, reductions=c("umap"), change.names=NULL){
  cbind(metadata(sc), lapply(1:length(reductions), function(x) Embeddings(sc, reduction=reductions[x])) %>% Reduce(cbind, .)   )
  
  
}

pcadf=function(sc){
  cbind(metadata(sc), Embeddings(sc, reduction="pca"))
}

getdate= function() format(Sys.time(), "%Y%b%d")


setlevels=function(so, varr, levls){
  Idents(so)=varr
  so@active.ident <- factor(x = so@active.ident, levels = levls)
  so}


col2names= function(dff, coll){
  rownames(dff)= dff[, coll]
  dff
}

condrbind=function(x,y){
  
  if(is.null(x)){
    return(y)}else{
      if(is.null(y)){return(x)}
      else{
        return(rbind(x,y))
      }}
  
}

condcbind=function(x,y){
  
  if(is.null(x)){
    return(y)}else{
      if(is.null(y)){return(x)}
      else{
        return(cbind(x,y))
      }}
  
}


triagecells.multi= function(so, ...){
  grps=list(...)# assemble groups into a list
  
  allgroups=lapply(grps, function(x) combineclusters2(so, x))
  mergername=function(x, refnames){
    if(x==0){
      return(refnames)
    }else{
      refnames[allgroups[[x]]]=LETTERS[x]
      return(mergername(x-1, refnames))
    }
  }
  out=as.character(so@meta.data[, "seurat_clusters"])
  final=mergername(length(allgroups), out)
  final= factor(final)
  final=factor(final, labels=1:length(levels(final)))
  final
}



###############
#integrate names into the dataframe
###############
namestocol=function(df){
  df[, "id"]=rownames(df)
  df
}

###############
#shortcut paste with "_" in between
###############

paste_=function(...) paste(..., sep="_")




tpng=function(nm,path=config$plotpath,  ...) png(paste0(path, timestamp(nm), ".png"),...)

tpdf=function(nm, path=config$plotpath, ...) pdf(paste0(path, timestamp(nm), ".pdf"),...)

#############################################################
#                   
#############################################################

replacepoint= function(num) paste(strsplit(as.character(num), split='\\.')[[1]], collapse="p")


#######################################################################################
#nearest neighbor operations.
#######################################################################################
#get a particular variable from nearest neighbors of a cell. 
#######################################################################################
getnnvar= function(so, cll, vr) metadata(so)[so@graphs$RNA_nn[cll, ]==1, vr]
getsnnvar= function(so, cll, vr) metadata(so)[so@graphs$RNA_snn[cll, ]==1, vr]

getnn= function(so, cll) metadata(so)[so@graphs$RNA_snn[cll, ]==1, ]



#######################################################################################
#get data slots from a seurat object
#######################################################################################

getcounts<- function(so) so@assays$RNA@counts



#######################################################################################
# Cluster merging operations
########################################################################################
####look for cells that belong to any of the cluster ids in the list and true them. specifically for scM
combineclusters= function(lst, outvar="seurat_clusters") Reduce('|', lapply(lst, function(el) metadata(scM)[,outvar]==el))
### more general to be applied to any seurat object
combineclusters2= function(so, lst, outvar="seurat_clusters") Reduce('|', lapply(lst, function(el) metadata(so)[,outvar]==el))
####################################
#merge clusters in list lst, in factor column
####################################
renamelabel=function(column, ..., renameto){
  grps=list(...)
  lapply(grps, function(lst){
    groupp=combineclusters(lst)
    sapply(1:length(scM@meta.data[,column]), function(x) ifelse(groupp[x],renameto, as.character(scM@meta.data[x, column])), USE.NAMES=F)
  })
}


triagecells.multi2=function(so, groupvar="seurat_clusters", ...){
  #... are the clusters in the group one after the other, ordered arbitrarily and grouped using c() when they are meant to be merged
  # warning: if some cluster groups overlap in any extent, cell indices will appear in two or more occasions and will be sequentially replaced!
  
  grps=list(...)
  grps2=grps %>% lapply(., function(x) paste0(addprefix(x, "^"), "$"))
  
  original=so %>% metadata %>% pull(get(groupvar))
  final.arr= rep(NA, length(original))
  indlist=lapply(grps2, function(gr) grepl( gr %>% paste0(., collapse="|"), original) %>% which )
  qlist=lapply(grps2, function(gr) gr %>% paste0(., collapse="|"))
  vector.stamp=function( vector, ind){
    if (ind<= length(grps)){ 
      vector[indlist[[ind]]]=ind
      return(vector.stamp( vector, ind+1))
    }else{
      vector[is.na(vector)]= paste_("original", original[is.na(vector)])
      return(vector) 
    }
  }
  
  out=vector.stamp(final.arr, 1)
}


###########################################
#replace each value in an array for another one.
###########################################

generalreplace=function(value, reff, target){
  
  target[which(reff==value)]
  
}

##############################################
#making the row names into a gene column
##############################################

names2col=function(dff, coln="gene") {
  dff[, coln]= rownames(dff) 
  dff
}

##########################################################
#mohamed's function to accomodate the slingshot data into a dataframe
##########################################################

tidy_ss_output <- function(ss_output){
  library(tidyr)
  library(dplyr)
  library(tibble)
  #cell embeddings
  cell_embedd <- ss_output@reducedDim
  #coordinates of the curve
  curve_coord <- lapply(names(ss_output@curves), function(x){
    `[[`(ss_output@curves[[x]],1) %>% 
      as.data.frame() %>% 
      rename_with(.,~paste0(x,"_",.))
  }) %>%
    bind_cols()
  #order of cell over infered curves
  cell_order <- lapply(names(ss_output@curves), function(x){
    `[[`(ss_output@curves[[x]],3) %>% 
      as.data.frame() %>% 
      rename_with(.,~paste0(x,"_lambda"))
  }) %>%
    bind_cols()
  #bind outputs in a single data frame
  cbind(cell_embedd,curve_coord,cell_order) %>% 
    as.data.frame() %>% 
    rownames_to_column("cell_id")
}


#############################################################################################################################################
#      PLOTTING FUNCTIONS
#############################################################################################################################################


#### modify hex colors up or down. 


dehash=function(color) substr(color, 2,7)
color2dec=function(color) hex2dec(substr(color, 2,7))
dec2color=function(dec) paste0('#', dec2hex(dec))
color2=function(dec) paste0('#', dec2hex(dec))

colormore=function(color, pts){
  
  newdec=color2dec(color)+pts
  
  if(newdec>hex2dec("FFFFFF")){newdec="FFFFFF"}
  if(newdec<hex2dec("000000")){newdec="000000"}
  
  return(dec2color(newdec))
}

colorless=function(color, pts){ #virtually the same as above but automatically subtracts the number so you don't have to worry about signs. 
  
  newdec=color2dec(color)-pts
  
  if(newdec>hex2dec("FFFFFF")){newdec="FFFFFF"}
  if(newdec<hex2dec("000000")){newdec="000000"}
  return(dec2color(newdec))
}


colormix=function(color1, color2){ #mixing two colors!
  getsubcolors= function(co) lapply(c(1,3,5), function(x) substr(dehash(co), x, x+1))
  
  
  newdec=color2dec(color)+colordec(color2)
  
  if(newdec>hex2dec("FFFFFF")){newdec="FFFFFF"}
  if(newdec<hex2dec("000000")){newdec="000000"}
  return(dec2color(newdec))
}

colorsubtract=function(color1, color2){ #subtracting two colors!
  
  newdec=abs(color2dec(color)-colordec(color2))
  
  if(newdec>hex2dec("FFFFFF")){newdec="FFFFFF"}
  if(newdec<hex2dec("000000")){newdec="000000"}
  return(dec2color(newdec))
}



colorbleach=function(color, pts){ ###ADDS PTS TO EACH OF THE r g b SECTIONS
  
  
  newdec=color2dec(color)+pts
  
  if(newdec>hex2dec("FFFFFF")){newdec="FFFFFF"}
  if(newdec<hex2dec("000000")){newdec="000000"}
  
  return(dec2color(newdec))
}

#########################################################
#ACCESSORY FUNCTIONS
#########################################################

metadata= function(x) x@meta.data


seuratprocessing=function(o){
  o=NormalizeData(o)
  o=FindVariableFeatures(o, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  o=ScaleData(o)
  o=RunPCA(o, dims=cfig$PCA_DIMS)
  o=ProjectDim(o)
  #o=RunUMAP(o, reduction="pca", dims=cfg$PCA_DIMS)
  o=FindNeighbors(o)
  o=FindClusters(o)
  o
}

SCTprocessing=function(o){
  o=NormalizeData(o)
  o=FindVariableFeatures(o, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  o=ScaleData(o)
  o=RunPCA(o, dims=cfig$PCA_DIMS)
  o=ProjectDim(o)
  #o=RunUMAP(o, reduction="pca", dims=cfg$PCA_DIMS)
  o=FindNeighbors(o)
  o=FindClusters(o)
  o
}

seuratprocessing_umap=function(o){
  o=NormalizeData(o)
  o=FindVariableFeatures(o, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  o=ScaleData(o)
  o=RunPCA(o, dims=cfig$PCA_DIMS)
  o=ProjectDim(o)
  o=RunUMAP(o, reduction="pca", dims=cfg$PCA_DIMS)
  o=FindNeighbors(o)
  o=FindClusters(o)
  o
}
getpname=function(x) unique(x@meta.data[, "orig.ident"])


seuratprocessing_markers1=function(so, verbose=F, umapdims=30, mincpm=3, minfc=1, fdr=0.05, selecttop=5, vfeatures=2000){
  so=SCTransform(so, verbose=verbose, variable.features.n=vfeatures)
  so=RunPCA(so, verbose=verbose)
  so=FindNeighbors(so, verbose=verbose)
  so=FindClusters(so, verbose=verbose)
  #so=RunUMAP(so, dims=1:di, verbose=verbose, return.model=T)
  
  
  cat(paste0("finding markers using the Riffle permutation test...\n"))  
  
  m=cb_pos_markers(counts=so@assays$RNA@counts, grouping=so %>% metadata %>% pull("seurat_clusters") )
  
  clusters=m %>% pull(group) %>% unique %>% as.numeric %>% sort
  ## arrange and filter best markers
  
  if(is.null(selecttop)){
    lt=lapply(clusters, function(x) m %>% filter(group==x,FDR<=fdr, logFC>=minfc, logCPM>=mincpm) %>% arrange(-logFC))
    topstr=""  
  }else{
    lt=lapply(clusters, function(x) m %>% filter(group==x,FDR<=fdr, logFC>=minfc, logCPM>=mincpm) %>% arrange(-logFC) %>% head(n=selecttop))
  }
  mp=lt %>% Reduce(rbind,.)
  
  cat("Number of filtered markers: ", nrow(mp), "(", (nrow(mp)/nrow(m))*100, "% markers post filtering)\n")
  so$RNA@misc$markers <- m
  so$RNA@misc$top_markers <- mp
  cat("Dataset", nm, ": calculating residuals for missing genes in scale data...\n")
  so <- GetResidual(so, features = so$RNA@misc$top_markers$feature, verbose = F)
  
  gc()
  so
}


seuratmarkers.delegate=function(so, method="deseq", group_column="seurat_clusters", replicate_column=NULL, padj=0.05, minfc=0.5,selecttop=NULL,minrate=0.5, fullreload=T, fullrecreate=F, getresidual=F, min.cells.per.group=3){
  fcat(paste0("finding markers using the DElegate package...\n") ) 
  #simpleCache(paste_("markers_DElegate_seuratproject",Project(so), "on_groups_of", group_column, "splitbyreps", replicate_column, "method", method), {
  
  #function (object, meta_data = NULL, group_column = NULL, replicate_column = NULL, 
  #method = "edger", min_rate = 0.05, min_fc = 1, lfc_shrinkage = NULL, 
  #verbosity = 1)   
  if(min.cells.per.group){
    
    countstab=so %>% metadata %>% group_by(!!sym(group_column)) %>% summarise(counts=n()) %>% as.data.frame
    
    goodcats=countstab[ countstab$counts>=min.cells.per.group, group_column]
    fcat("enough cells found for", paste(goodcats, collapse=","))
    
    so = so[, (metadata(so) %>% filter(!!sym(group_column) %in% goodcats) %>% rownames)]
    
  }
  
  
  m= DElegate::FindAllMarkers2(object=so, group_column=group_column, replicate_column=replicate_column, method=method, min_rate=minrate)
  m
  #m=cb_pos_markers(scwt@assays$RNA@counts, grouping=sma %>% metadata %>% pull() )
  #}, assignToVar="m", reload=fullreload, recreate=fullrecreate)
  
  #print(m)
  
  if(group_column %in% (allcolors %>% names)){
    clusters=allcolors[[group_column]][allcolors[[group_column]] %>% names %in% (so %>% metadata %>% pull(group_column) %>% unique)] %>% names
  }else{
    clusters=so %>% metadata %>% pull(group_column) %>% unique 
    
  }
  
  ## arrange and filter best markers.
  ## currently filtering such that the minimum rate for in group cells is at least minrate. 
  pa=padj
  mi=minfc
  if(is.null(selecttop)){
    lt=lapply(clusters, function(x) m %>% filter(group1==x,padj<=pa, log_fc>=mi, rate1>=minrate) %>% arrange(-log_fc))
    topstr=""  
  }else{
    lt=lapply(clusters, function(x) m %>% filter(group1==x,padj<=pa, log_fc>=mi, rate1>=minrate) %>% arrange(-log_fc) %>% head(n=selecttop))
  }
  mp=lt %>% Reduce(rbind,.)
  
  cat("Number of filtered markers: ", nrow(mp), "(", (nrow(mp)/nrow(m))*100, "% markers post filtering)\n")
  so[["RNA"]]@misc$markers <- m
  so[["RNA"]]@misc$top_markers <- mp
  if(getresidual){
    so <- GetResidual(so, features = so$RNA@misc$top_markers$feature, verbose = F)
  }
  so
}



################################################################################
# calculate markers using DE, each versus each other
################################################################################


seuratmarkers.delegateDE=function(so, method="deseq", group_column="seurat_clusters", compare="all_vs_all", replicate_column=NULL, padj=0.05, minfc=0.5,selecttop=NULL,minrate=0.5, fullreload=T, fullrecreate=F, getresidual=F, min.cells.per.group=3, lfc_shrinkage=NULL){
  #the other option for "compare" is each vs rest"
  fcat(paste0("finding markers using the DElegate package...\n") ) 
  simpleCache(paste_("markers_DElegate_seuratproject",Project(so), "on_groups_of", group_column, "splitbyreps", replicate_column, "method", method), {
    
    #function (object, meta_data = NULL, group_column = NULL, replicate_column = NULL, 
    #method = "edger", min_rate = 0.05, min_fc = 1, lfc_shrinkage = NULL, 
    #verbosity = 1)   
    if(min.cells.per.group){
      
      countstab=so %>% metadata %>% group_by(!!sym(group_column)) %>% summarise(counts=n()) %>% as.data.frame
      
      goodcats=countstab[ countstab$counts>=min.cells.per.group, group_column]
      fcat("enough cells found for", paste(goodcats, collapse=","))
      
      so = so[, (metadata(so) %>% filter(!!sym(group_column) %in% goodcats) %>% rownames)]
      
    }
    
    
    m= DElegate::findDE(object=so, group_column=group_column, compare=compare, replicate_column=replicate_column, method=method,lfc_shrinkage=NULL)
    m
    #m=cb_pos_markers(scwt@assays$RNA@counts, grouping=sma %>% metadata %>% pull() )
  }, assignToVar="m", reload=fullreload, recreate=fullrecreate)
  
  #print(m)
  
  if(group_column %in% (allcolors %>% names)){
    clusters=allcolors[[group_column]][allcolors[[group_column]] %>% names %in% (so %>% metadata %>% pull(group_column) %>% unique)] %>% names
  }else{
    clusters=so %>% metadata %>% pull(group_column) %>% unique 
    
  }
  
  ## arrange and filter best markers.
  ## currently filtering such that the minimum rate for in group cells is at least minrate. 
  pa=padj
  mi=minfc
  if(is.null(selecttop)){
    lt=lapply(clusters, function(x) m %>% filter(group1==x,padj<=pa, log_fc>=mi, rate1>=minrate) %>% arrange(-log_fc))
    topstr=""  
  }else{
    lt=lapply(clusters, function(x) m %>% filter(group1==x,padj<=pa, log_fc>=mi, rate1>=minrate) %>% arrange(-log_fc) %>% head(n=selecttop))
  }
  mp=lt %>% Reduce(rbind,.)
  
  cat("Number of filtered markers: ", nrow(mp), "(", (nrow(mp)/nrow(m))*100, "% markers post filtering)\n")
  so[["RNA"]]@misc$markers <- m
  so[["RNA"]]@misc$top_markers <- mp
  if(getresidual){
    so <- GetResidual(so, features = so$RNA@misc$top_markers$feature, verbose = F)
  }
  so
}


############################################
#function to  get the sample from a barcode.
############################################
getsamplegroupfrombc=function(bc){
  if (bc=="Doublet"){
    return("Doublet")
  }else{
    if (bc=="Negative"){
      return("Negative")
    }else{
      
      grp=paste(unname(as.matrix(dA[multiseq_sequence== bc][,.(sample_group, stage,treatment, experiment_group)][1])[1,c(1:4)]), collapse="_")
      grp
      
    }
  }
  
}
####### get only stage from barcode
getstagefrombc=function(bc){
  if (bc=="Doublet"){
    return("Doublet")
  }else{
    if (bc=="Negative"){
      return("Negative")
    }else{
      
      grp=paste(unname(as.matrix(dA[multiseq_sequence== bc][,.(stage)][1])[1,1]), collapse="_")
      grp
      
    }
  }
}


####### get only stage from barcode
getconditionfrombc=function(bc){
  if (bc=="Doublet"){
    return("Doublet")
  }else{
    if (bc=="Negative"){
      return("Negative")
    }else{
      
      grp=paste(unname(as.matrix(dA[multiseq_sequence== bc][,.(sample_group)][1])[1,1]), collapse="_")
      grp
      
    }
  }
}

####### get only induction from barcode
gettreatmentfrombc=function(bc){
  if (bc=="Doublet"){
    return("Doublet")
  }else{
    if (bc=="Negative"){
      return("Negative")
    }else{
      
      grp=paste(unname(as.matrix(dA[multiseq_sequence== bc][,.(treatment)][1])[1,1]), collapse="_")
      grp
      
    }
  }
}

####replicate
getreplicatefrombc=function(bc){
  if (bc=="Doublet"){
    return("Doublet")
  }else{
    if (bc=="Negative"){
      return("Negative")
    }else{
      
      grp=paste(unname(as.matrix(dA[multiseq_sequence== bc][,.(replicate)][1])[1,1]), collapse="_")
      grp
      
    }
  }
}


recreateAll <- F

#dA <- fread(metapath)
#dA <- merge(dA[bsf_name!="",], fread(barcodepath)[,.(multiseq_id, multiseq_sequence)], by="multiseq_id", all.x=T)
#setkey(dA, sample_name)

#################make roc curve
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}


############Quality Ccontrol QC plots function
logvarname= function(x) paste0("log10_")


qcplots=function(so, 
                 plotred=FALSE, 
                 plotsize=10, 
                 plotname="qcplots",
                 plot.title=NULL,
                 demux= NULL, 
                 topx=10000,
                 topf=20000, 
                 ww=1, 
                 hh=1,
                 colorby="seurat_clusters",
                 format="png",
                 nott=scales::comma,
                 path="./", 
                 cutoffs=NULL, 
                 cellnumbers=NULL, 
                 seriate.method=NULL,
                 plot.empty=NULL, 
                 groupvar="group1",
                 lfcvar="log_fc"){
  
  Idents(so)=colorby  
  
  markers= so@assays$RNA@misc$top_markers
  
  if(is.null(markers)){
    fcat("There are no markers! please add markers under @assays$RNA@misc$top_markers" )
    
    
    
    cells=NULL
  }else{
    if(is.null(seriate.method)){            
      cells <- WhichCells(so, downsample = 100)
    }else{
      cells= seriatecells(so, clusvar=colorby, meth=seriate.method, groupvarname=groupvar, lfcvarname = lfcvar)
    }
  }   
  
  
  if(is.null(plot.title)){
    plot.title=Project(so)
  }
  
  
  require(scales)
  tsz=7
  
  
  so$log10_nCount_RNA <- log10(so$nCount_RNA)
  so$log10_nFeature_RNA <- log10(so$nFeature_RNA)
  
  qc1=ggplot(so@meta.data, aes(x=nCount_RNA))+
    geom_histogram(binwidth=20, alpha=0.5)+
    labs(x = "Read depth",y = "Frequency")+
    guides(fill=guide_legend())+
    coord_cartesian(xlim=c(0, topf))+
    scale_y_continuous(labels = nott)+
    scale_x_continuous(labels = nott)+
    theme_classic()+
    labs(title="Reads per barcode")+theme_classic()+theme(axis.text=element_text(size=tsz), axis.title=element_text(size=tsz));
  
  #distribution of counts per barcode
  qc2=ggplot(so@meta.data, aes(x=nFeature_RNA))+
    labs(x = "Num. Genes",y = "Frequency")+
    geom_histogram(binwidth=20, alpha=0.5)+
    labs(title="Genes per barcode")+
    scale_y_continuous(labels = nott)+
    scale_x_continuous(labels = nott)+
    theme_classic()+
    coord_cartesian(xlim=c(0, topx))+theme_classic()+theme(axis.text=element_text(size=tsz), axis.title=element_text(size=tsz));
  
  
  if(!(colorby %in% names(allcolors)) || allcolors[[colorby]] %>% length != so %>% metadata %>% pull(!!sym(colorby)) %>% unique %>% length){
    
    cats=so %>% metadata %>% pull(!!sym(colorby)) %>% unique
    allcolors[[colorby]]=  randomcolors(cats %>% length)
    names(allcolors[[colorby]])= cats
  }
  
  qc3 <- FeatureScatter(so, group.by=colorby, feature1 = "nFeature_RNA", feature2 = "percent.mt", raster=T)+
    scale_color_manual(values=allcolors[[colorby]])+
    labs(x = "Num. Genes",y = "% Mitochondrial")+
    ggtitle("")+theme_classic()+
    scale_x_continuous(labels = scientific)+
    theme(axis.text.x=element_text(size=tsz), axis.text.y=element_text(size=tsz),axis.title=element_text(size=tsz))+NoLegend()
  qc4 <- FeatureScatter(so, group.by=colorby,feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster=T)+
    scale_color_manual(values=allcolors[[colorby]])+
    labs(x = "Read depth",y = "Num.Genes")+
    ggtitle("")+
    scale_y_continuous(labels = nott)+
    scale_x_continuous(labels = nott)+
    theme_classic()+
    theme(axis.text.x=element_text(size=tsz), axis.text.y=element_text(size=tsz), axis.title=element_text(size=tsz),legend.position="top",legend.title=element_blank())
  qc5 <- FeatureScatter(so, group.by=colorby,feature1 = "nCount_RNA", feature2 = "percent.mt", raster=T)+
    scale_color_manual(values=allcolors[[colorby]])+
    labs(x = "Read depth",y ="% Mitochondrial")+
    ggtitle("")+theme_classic()+
    scale_x_continuous(labels = scientific)+
    theme(axis.text.x=element_text(size=tsz), axis.text.y=element_text(size=tsz), axis.title=element_text(size=tsz))+NoLegend()  
  
  
  if(plotred!=FALSE){
    
    dms=10 
    umap1=DimPlot(so, group.by = colorby, reduction=plotred, label=T, label.size=6, label.color="black", cols= allcolors[[colorby]], raster=T)+
      ggtitle(paste0(plot.title, ": ", plotred, "\n", ncol(so), " cells"))+NoAxes()#+NoLegend()#+scale_color_manual( values=qccols);tsne1 
    
    umapmt=FeaturePlot(so, feature="percent.mt", reduction=plotred, raster=T)+
      ggtitle("% Mitochondrial")+NoAxes()#+NoLegend()#+scale_color_manual( values=qccols);tsne1 
    
    umapreads=FeaturePlot(so,  feature="nCount_RNA",reduction=plotred, raster=T)+
      ggtitle("Reads/BC")+NoAxes()#+NoLegend()#+scale_color_manual( values=qccols);tsne1 
    
    umapgenes=FeaturePlot(so,feature="nFeature_RNA",  reduction=plotred, raster=T)+
      ggtitle("Genes/BC")+NoAxes()#+NoLegend()#+scale_color_manual( values=qccols);tsne1 
    
    vlnreads= ggplot(so %>% metadata, aes(x=factor(!!sym(colorby)), fill=factor(!!sym(colorby)), y=nCount_RNA))+geom_violin()+scale_y_continuous(trans='log10')+scale_fill_manual(values=allcolors[[colorby]])+theme_classic()+NoLegend()
    vlngenes= ggplot(so %>% metadata, aes(x=factor(!!sym(colorby)), fill=factor(!!sym(colorby)), y=nFeature_RNA))+geom_violin()+scale_fill_manual(values=allcolors[[colorby]])+theme_classic()+NoLegend()
    vlnmt= ggplot(so %>% metadata, aes(x=factor(!!sym(colorby)), fill=factor(!!sym(colorby)), y=percent.mt))+scale_fill_manual(values=allcolors[[colorby]])+geom_violin()+theme_classic()+NoLegend()
    
    
    if( !is.null(cutoffs)){
      lcolor="red"
      vlnreads= vlnreads+geom_hline(yintercept= cutoffs$counts, color=lcolor)
      vlngenes= vlngenes+geom_hline(yintercept= cutoffs$features, color=lcolor)
      vlnmt= vlnmt+geom_hline(yintercept= cutoffs$mito, color=lcolor)
    }
    
    if( !is.null(cellnumbers)){
      ncolor="black"
      sumdata= so %>% metadata %>% group_by(seurat_clusters) %>% summarise(counts=n(), nCount_RNA=median(nCount_RNA), nFeature_RNA=median(nFeature_RNA), percent.mt=median(percent.mt))
      
      vlnreads= vlnreads+geom_text(data=sumdata, aes( label=paste0("n=",counts)), angle=90, col=ncolor)
      vlngenes= vlngenes+geom_text(data=sumdata, aes( label=paste0("n=",counts)), angle=90, col=ncolor)
      vlnmt= vlnmt+geom_text(data=sumdata, aes( label=paste0("n=",counts)), angle=90, col=ncolor)
    }
    
    if(!is.null(demux) && !is.null(so$RNA@misc$top_markers)){
      
      demuxplot=ggplot(so %>% metadata, aes(x=seurat_clusters))+geom_bar(aes(fill=!!sym(demux)), position="fill")+theme_classic()
      
      
      
      layout="AAAAMMMMM
             AAAALLLLL
                    AAAALLLLL
                    AAAALLLLL
                    AAAALLLLL
                    EEBBIIIII
                    EEBBIIIII
                    FFCCJJJJJ
                    FFCCJJJJJ
                    GGDDKKKKK
                    GGDDKKKKK
                    "
      
      
      
      fcat("Preparing heatmap...")
      
      #if(!is.null(plot.empty)){
      # emptydropsplot=ggplot(so %>% metadata[cells, ] %>% dplyr::mutate(ord=1:length(cells)), aes(x=ord, y=1, fill=isemptydroplet))+geom_bar(color=NA)
      #}else{
      #emptydropsplot=ggplot(so %>% metadata[cells, ] %>% dplyr::mutate(ord=1:length(cells)), aes(x=ord, y=1, fill=isemptydroplet))+geom_bar(color=NA)
      #}             
      
      
      
      ph <- DoHeatmap(so, features = so$RNA@misc$top_markers$feature, group.colors=allcolors[[colorby]], slot = "scale.data", cells = cells) + NoLegend()
      
      
      
      
      qc_assembly=wrap_plots(A = umap1, B = qc3, C = qc4,
                             D = qc5, E = umapmt, F=umapreads, G=umapgenes, I=vlnmt, J=vlnreads, K=vlngenes, L=ph,M=demuxplot,
                             design = layout)
    }else{
      
      
      qc_assembly=((umap1)/(qc3+qc4+qc5))|(umapmt+umapreads+umapgenes)|(vlnmt+vlnreads+vlngenes)#/(qc6+qc7+qc8)
    }
    ###################
    #plotting with tsne  
    ###################
    tpng(plotname, path=path, width=w*plotsize*ww, height=w*plotsize*hh, res=600)
    print(qc_assembly)
    dev.off()
    
  }else{
    
    #ppc1=DimPlot(object = subset(so, subset = percent.mt<10), cols= c("grey", "blue"), reduction = "pca");
    ###################
    #plotting without tsne  
    ###################
    
    layout <- "
####NN
AABBNN
AABBNN
CCDDEE
CCDDEE
"
    qc_assembly=wrap_plots(A = qc1, qc2 = qc2, C = qc3, D=qc4+NoLegend(), E=qc5, N=qc4,  design = layout)
    
    
    png.pdf.prop=12
    if(format=="png"){
      tpng(plotname, path=path, width=w*plotsize*ww/png.pdf.prop, height=w*plotsize*hh/png.pdf.prop, res=600)
      print(qc_assembly)
      dev.off()
    }
    if(format=="pdf"){
      tpdf(plotname, path=path, width=pw*plotsize*ww/png.pdf.prop, height=pw*plotsize*hh/png.pdf.prop)
      print(qc_assembly)
      dev.off()
    }
    
  }
  list(plot=qc_assembly, cells=cells)
}


os <- function() {
  require(dplyr)
  gc();
  message("Objects in MB:")
  objects = ls(envir=.GlobalEnv);
  classes = sapply(objects, function(object) { class(get(object))[1] });
  sizes = sapply(objects, function (object) { object.size(get(object)) } )
  a = data.frame(MB=sizes/1e6, class=classes)
  ord = order(sizes, decreasing=TRUE)
  a2 = a[ord,];
  a2 = a2[! a2$class == "function",]
  print(head(a2, 30))
  message("Sum: ", signif(sum(sizes/1e6),4), " MB (", signif(sum(sizes/1e9),4), " GB)")
  # fcat("total: ", a2 %>% summarise(total=sum(MB)/1000) %>% pull(total) %>% as.numeric , "GB" 
  a2
}
#os()



###########################################################################
# FUNCTION TO KEEP THE BARE BONES MINIMUM SEURAT DATA, REMOVING ALL UNNNECESSARY ASSAYS ETC.
###########################################################################

slimseurat=function(so){
  DietSeurat( so,counts = TRUE,data = TRUE,scale.data = TRUE,features = NULL,
              assays = c("RNA","refAssay"))#,dimreducs = c("pca", "umap"))
}


###data=FALSE curently not supported
slimerseurat=function(so){
  DietSeurat( so,counts = TRUE,scale.data = FALSE,features = NULL,
              assays = "RNA", dimreducs=NULL)#,dimreducs = c("pca", "umap"))
}


############################################################################
#  function to get a plotting table with some genes and some cells together with the metadata
############################################################################
join_meta_exp=function(so, genes, cells=NULL, normrow=FALSE, lognormrow=FALSE, metacols=NULL, assay="RNA"){
  
  if(is.null(cells)){cells= colnames(so)}
  
  m=as.data.frame(t(as.matrix(so[genes,cells ]@assays[[assay]]@data)))
  if(normrow){
    m=apply(m, 1, normvec)
  }else{
    if(lognormrow){
      m=apply(m, 1, lognormvec)
    }
    
  }
  
  if (is.null(metacols)){metacols=colnames(metadata(so))}  
  cbind( metadata(so)[cells, metacols], m)
  
}

#############################################################################
#function to attach metadata and gene reads.
#############################################################################


genesdf=function(so, genes) as.matrix((so@assays$RNA@data))[genes,  ]    


####### function to change the legend size inside guide
legendsize=function(x) guide_legend(override.aes = list(size = x) )

##############################################################################
#resizing overall text in the plot
resizetext=function(x) theme_classic(base_size = x)

rotatex= function(x) theme(axis.text.x = element_text(angle = x))


################################################################################
#Heatmap accessory functions
################################################################################

###function to find the names of the genes in each cluster
genesinclusters=function(phe, hvarmat,  top=NULL){
  
  lapply(row_order(phe), function(x) {
    
    if(!is.null(top)){
      rownames(hvarmat)[x][1:top]
    }else{
      rownames(hvarmat)[x]
    }
  }
  )
}



topnclustmembers=function(phe, top=3){
  
  Reduce(c, lapply(row_order(phe), function(x) {
    x[1:top]
    
  }
  ))
}
###############################################################

getheatmapcolors=function(ph, var, orientation='top_annotation'){
  base::eval(glue::glue("pha@{orientation}@anno_list${var}@color_mapping@colors"))
  colls}

#####################################################################

namenums=function(x, prefix=NULL, zero=F){
  if(is.null(prefix)){
    
    if(zero){n=1}else{n=0}
    names(x)= (1:length(x))-n; return(x)}else{
      names(x)= paste0(prefix, 1:length(x))
    }
  
  x
}

########################################################################

showpalette=function(p, pname="generic"){
  sca=10
  rr=600
  w=250
  tpng(paste_("palette", pname), res=rr, he=sca*w, wi=sca*w*2)
  
  pl=ggplot(1:length(p) %>% as.data.frame %>% dplyr::mutate(h=1, colors=p, name=names(p)), aes(x=factor(.), y=factor(h), label=colors))+
    geom_col(aes(fill=factor(.)))+
    scale_fill_manual(values=namenums(p))+
    theme_classic()+
    coord_cartesian(ylim=c(0,1))+
    geom_text(y=0.5, color="black", angle=90)+
    #NoAxes()+
    NoLegend()+
    ggtitle(pname)
  
  if(is.null(names(p))){
    print(pl )
  }else{
    print(pl+geom_text(inherit.aes=T,aes(label=name), y=1.1, color="black", angle=90))  
  }
  dev.off()
  
  return(pl+geom_text(inherit.aes=T,aes(label=name), y=1.1, color="black", angle=90))
}
#showpalette(p)

##################################################################################
#GGUMAP wrapper  and derivatives
##################################################################################

ggumap=function(so, umap.df=NULL, colorby="seurat_clusters", labelby="seurat_clusters", sz=0.02, ssca=300, colorlist=NULL, reductions=c("umap"), reduction.key="UMAP_", legend.size=3, color.labels.by="defaultcolor"){
  require(ggrepel)
  require(ggrastr)
  dim1=paste0(reduction.key, 1)
  dim2=paste0(reduction.key, 2)
  
  ##############################################################################
  # sourcing data
  ##############################################################################
  if(!is.null(umap.df)){
    cats=umap.df %>% pull(!!sym(colorby)) %>% unique
    dat= umap.df %>% dplyr::mutate(UMAP_1=!!sym(dim1), UMAP_2=!!sym(dim2)) %>% dplyr::mutate(none="")
  }else{
    cats=so %>% metadata %>% pull(!!sym(colorby)) %>% unique
    dat=umapdf(so, reductions=reductions) %>% dplyr::mutate(UMAP_1=!!sym(dim1), UMAP_2=!!sym(dim2)) %>% dplyr::mutate(none="")
  }
  ##############################################################################
  # preparing colorlist
  ##############################################################################
  
  # colorlist doesn't exist.
  if (is.null(colorlist[[colorby]])){
    colorlist[[colorby]]= rainbow(length(cats)) %>% givename(., cats)
    fcat("colors used:")
    
    colorlist[[colorby]] %>% dput
  }
  
  ### actual categories are smaller than total number of colors
  if (length(colorlist[[colorby]] %>% names) > length(cats) && any(cats %in% colorlist[[colorby]] %>% names) ){
    fcat("adjusting colorlist to present labels")
    
    colorlist[[colorby]]= allcolors[[colorby]][cats]
  }
  ### there are more categories than colors present
  if (length(colorlist[[colorby]] %>% names) < length(cats) && any(cats %in% colorlist[[colorby]] %>% names) ){
    fcat("adjusting colorlist to present labels")
    missinglabels= setdiff(cats, colorlist[[colorby]] %>% names)
    
    colorlist[[colorby]]= c(allcolors[[colorby]][cats], randomcolors(missinglabels) %>% givename(., missinglabels))
  }
  colorlist[["defaultcolor"]]=c(one="black")
  colorlist[[color.labels.by]]=c(colorlist[[color.labels.by]] )
  colorlist[[colorby]]=c(colorlist[[colorby]] )
  if(colorby=="seurat_clusters"){
    colorlist[["seurat_clusters"]]=adjust.seurat.colors(so)
  }
  ##############################################################################
  # constructing plot
  ##############################################################################
  
  
  
  dat2=dat %>% group_by(!!sym(labelby)) %>% summarise(UMAP_1=jitter(mean(UMAP_1)), UMAP_2=jitter(mean(UMAP_2)), defaultcolor="one")
  #catch some pesky history dependent bug with geom_point_rast
  
  gpuc<- tryCatch(
    {
      # First try with geom_point_rast
      print(ggplot(dat, aes(x=UMAP_1, y=UMAP_2, label=!!sym(labelby)))+
              geom_point_rast(size=sz, aes(color=!!sym(colorby)))+
              scale_color_manual(values=colorlist[[colorby]], guide = guide_legend(override.aes = list(size = legend.size, shape=15) ))+
              theme_classic()+
              geom_text(data= dat2,  color="black", size=sz*ssca))
      
      ggplot(dat, aes(x=UMAP_1, y=UMAP_2, label=!!sym(labelby)))+
        geom_point_rast(size=sz, aes(color=!!sym(colorby)))+
        scale_color_manual(values=colorlist[[colorby]], guide = guide_legend(override.aes = list(size = legend.size, shape=15) ))+
        theme_classic()+
        geom_text(data= dat2,  color="black", size=sz*ssca)
    },
    error = function(e) {
      if (grepl("Empty raster", e$message)) {
        message("geom_point_rast failed (Empty raster). Falling back to geom_point().")
        ggplot(dat, aes(x=UMAP_1, y=UMAP_2, label=!!sym(labelby)))+
          geom_point(size=sz, aes(color=!!sym(colorby)))+
          scale_color_manual(values=colorlist[[colorby]], guide = guide_legend(override.aes = list(size = legend.size, shape=15) ))+
          theme_classic()+
          geom_text(data= dat2,  color="black", size=sz*ssca)
      } else {
        stop(e) # rethrow unexpected errors
      }
    }
  )
  
  
  gpuc}


ggumap3=function(so, umap.df=NULL,
                 colorby="seurat_clusters",
                 labelby="seurat_clusters",
                 glassworkby="seurat_clusters",
                 glasswork.params=NULL,
                 sz=0.02,
                 ssca=300,
                 colorlist=allcolors,
                 reductions=c("umap"),
                 reduction.key="UMAP_", 
                 legend.size=3){
  require(ggrepel)
  require(concaveman)
  
  dim1=paste0(reduction.key, 1)
  dim2=paste0(reduction.key, 2)
  
  ##############################################################################
  # sourcing data
  ##############################################################################
  if(!is.null(umap.df)){
    cats=umap.df %>% pull(!!sym(colorby)) %>% unique
    dat= umap.df %>% dplyr::mutate(UMAP_1=!!sym(dim1), UMAP_2=!!sym(dim2)) %>% dplyr::mutate(none="")
  }else{
    cats=so %>% metadata %>% pull(!!sym(colorby)) %>% unique
    dat=umapdf(so, reductions=reductions) %>% dplyr::mutate(UMAP_1=!!sym(dim1), UMAP_2=!!sym(dim2)) %>% dplyr::mutate(none="")
  }
  ##############################################################################
  # preparing colorlist
  ##############################################################################
  
  # colorlist doesn't exist.
  if (is.null(colorlist[[colorby]])){
    colorlist[[colorby]]= randomcolors(length(cats)) %>% givename(., cats)
  }
  
  ### actual categories are smaller than total number of colors
  if (length(colorlist[[colorby]] %>% names) > length(cats) && any(cats %in% colorlist[[colorby]] %>% names) ){
    fcat("adjusting colorlist to present labels")
    
    colorlist[[colorby]]= allcolors[[colorby]][cats]
  }
  ### there are more categories than colors present
  if (length(colorlist[[colorby]] %>% names) < length(cats) && any(cats %in% colorlist[[colorby]] %>% names) ){
    fcat("adjusting colorlist to present labels")
    missinglabels= setdiff(cats, colorlist[[colorby]] %>% names)
    
    colorlist[[colorby]]= c(allcolors[[colorby]][cats], randomcolors(missinglabels) %>% givename(., missinglabels))
  }
  
  colorlist[["defaultcolor"]]=c(one="black")
  
  ##############################################################################
  # constructing plot
  ##############################################################################
  
  gpuc=ggplot(dat, aes(x=UMAP_1, y=UMAP_2, label=!!sym(labelby)))+
    geom_point_rast(size=sz, aes(color=!!sym(colorby)))+
    scale_color_manual(values=colorlist[[colorby]], guide = guide_legend(override.aes = list(size = legend.size) ))+
    theme_classic()+
    #NoLegend()+
    #geom_text_repel(data=dat %>% group_by(!!sym(labelby)) %>% summarise(UMAP_1=mean(UMAP_1), UMAP_2=mean(UMAP_2)), size=sz*ssca)
    geom_text(data=dat %>% group_by(!!sym(labelby)) %>% summarise(UMAP_1=mean(UMAP_1), UMAP_2=mean(UMAP_2), defaultcolor="one"), aes(color=color.labels.by), size=sz*ssca)
  
  if(glassworkby=="none" |is.na(glassworkby) |is.null(glassworkby)| glassworkby=="null" ){
    gpuc
  }else{
    
    if(!is.null(glasswork.params)){
      msg("making glasswork hulls")
      glasswork.params=list(
        size=1,
        fillvar="null",
        linecolor="#000000",
        concavity=2, 
        fillcols="#FFFFFF00",
        clus=dat %>% pull(!!sym(glassworkby)) %>% unique %>% as.character
      )
    }
    
    polys=clusterhull3(NULL, umap.df=dat, clusvar=glassworkby, clus=glasswork.params$clus, fillvar="null", linecolor=glasswork.params$linecolor, size=glasswork.params$size, reduction.key=reduction.key, fillcols=glasswork.params$fillcols, concavity=glasswork.params$concavity)
    gpuc+polys+NoLegend()
    
  }}




















ggumap2=function(so, colorby="seurat_clusters", labelby="seurat_clusters", sz=0.02, ssca=300, colorlist=allcolors, reductions=c("umap"), reduction.key="UMAP_"){
  require(ggrepel)
  cats=so %>% metadata %>% pull(!!sym(colorby)) %>% unique
  if (length(colorlist[[colorby]])!= length(cats)){
    colorlist[[colorby]]= randomcolors(length(cats))
  }
  #udimnames=c(paste0(reduction.key, 1), paste0(reduction.key, 2))
  dim1=paste0(reduction.key, 1)
  dim2=paste0(reduction.key, 2)
  dat=umapdf(so, reductions=reductions) %>% dplyr::mutate(UMAP_1=!!sym(dim1), UMAP_2=!!sym(dim2)) %>% dplyr::mutate(none="")
  gpuc=ggplot(dat, aes(x=!!sym(udimnames[1]), y=!!sym(udimnames[2]), label=!!sym(labelby)))+
    geom_point_rast(size=sz, aes(color=!!sym(colorby)))+
    scale_color_manual(values=colorlist[[colorby]])+
    theme_classic()+
    #NoLegend()+
    #geom_text_repel(data=dat %>% group_by(!!sym(labelby)) %>% summarise(UMAP_1=mean(UMAP_1), UMAP_2=mean(UMAP_2)), size=sz*ssca)
    
    geom_text(data=dat %>% group_by(!!sym(labelby)) %>% summarise(UMAP_1=mean(UMAP_1), UMAP_2=mean(UMAP_2)), size=sz*ssca)
  
  gpuc}


####return a list containng the umap and the labels separate. 
ggumapl=function(so, colorby="seurat_clusters", labelby="seurat_clusters", sz=0.02, ssca=300, colorlist=allcolors, reductions=c("umap"), reduction.key="UMAP_"){
  library(ggrepel)
  dim1=paste0(reduction.key, 1)
  dim2=paste0(reduction.key, 2)
  dat=umapdf(so, reductions=reductions) %>% rename(UMAP_1=!!sym(dim1), UMAP_2=!!sym(dim2)) %>% dplyr::mutate(none="")
  gpuc=ggplot(dat, aes(x=UMAP_1, y=UMAP_2, label=!!sym(labelby)))+
    geom_point(size=sz, aes(color=!!sym(colorby)))+
    scale_color_manual(values=colorlist[[colorby]])+
    theme_classic()
  #NoLegend()+
  #geom_text_repel(data=dat %>% group_by(!!sym(labelby)) %>% summarise(UMAP_1=mean(UMAP_1), UMAP_2=mean(UMAP_2)), size=sz*ssca)
  
  list(gpuc, 
       geom_text_repel(data=dat %>% group_by(!!sym(labelby)) %>% summarise(UMAP_1=mean(UMAP_1), UMAP_2=mean(UMAP_2)), size=sz*ssca)
  )
}

###beta
ggumapc=function(so, colorby="seurat_clusters", labelby="seurat_clusters", sz=0.02, ssca=300, colorlist=allcolors, reductions=c("umap"), reduction.key="UMAP_"){
  library(ggrepel)
  if (sapply(so %>% metadata %>% pull(!!sym(colorby)), is.double) %>% any){
    colorf=scale_color_manual(values=colorlist[[colorby]])
  }else{colorf=scale_color_gradient(values=colorlist[[colorby]])}
  
  dim1=paste0(reduction.key, 1)
  dim2=paste0(reduction.key, 2)
  dat=umapdf(so, reductions=reductions) %>% rename(UMAP_1=!!sym(dim1), UMAP_2=!!sym(dim2))
  gpuc=ggplot(dat, aes(x=UMAP_1, y=UMAP_2, label=!!sym(labelby)))+
    geom_point(size=sz, aes(color=!!sym(colorby)))+
    scale_color_manual(values=colorlist[[colorby]])+
    theme_classic()+
    NoLegend()+
    geom_text_repel(data=dat %>% group_by(!!sym(labelby)) %>% summarise(UMAP_1=mean(UMAP_1), UMAP_2=mean(UMAP_2)), size=sz*ssca)
  
  gpuc}



ggumaph=function(so, colorby="seurat_clusters", sz=0.02, ssca=300, colorlist=allcolors, reductions=c("umap"), reduction.key="UMAP_"){
  library(ggrepel)
  dim1=paste0(reduction.key, 1)
  dim2=paste0(reduction.key, 2)
  dat=umapdf(so, reductions=reductions) %>% rename(UMAP_1=!!sym(dim1), UMAP_2=!!sym(dim2)) %>% dplyr::mutate(none="", ccondition=factor(ccondition, levels=names(colorlist[["ccondition"]])))
  gpuc=ggplot(dat, aes(x=UMAP_1, y=UMAP_2))+
    geom_point(size=sz, aes(color=!!sym(colorby)))+
    gghighlight()+
    scale_color_manual(values=colorlist[[colorby]])+
    theme_classic()
  #NoLegend()
  
  gpuc}

ggumaps=function(so, colorby="seurat_clusters", sz=0.02, ssca=300, colorlist=allcolors, reductions=c("umap"), reduction.key="UMAP_"){
  library(ggrepel)
  dim1=paste0(reduction.key, 1)
  dim2=paste0(reduction.key, 2)
  dat=umapdf(so, reductions=reductions) %>% rename(UMAP_1=!!sym(dim1), UMAP_2=!!sym(dim2)) %>% dplyr::mutate(none="", ccondition=factor(stage, levels=names(colorlist[["stage"]])))
  gpuc=ggplot(dat, aes(x=UMAP_1, y=UMAP_2))+
    geom_point(size=sz, aes(color=!!sym(colorby)))+
    gghighlight()+
    scale_color_manual(values=colorlist[[colorby]])+
    theme_classic()
  
  gpuc}
#################################

addprefix= (function(cc, prefix="C")  paste0(prefix, cc)) %>% Vectorize(., USE.NAMES=F)



givecolnames=function(x, ind=NULL, nms){
  if (is.null(ind)){
    ind= 1:length(nms) 
  }
  colnames(x)[ind]= nms
  x
  
}

giverownames=function(x, nms){
  rownames(x)= nms
  x
  
}

#####

findgeneindex=function(so, x)   (so %>% rownames) %>% grepl(paste0("^", x,"$"), .) %>% which


########################

na2zero=function(x) ifelse(is.na(x),0, x) 



maptoref= function(query, refdataset, mapvar, pdms=1:50, return.merge.metadata=T, return.query.metadata=F, return.lean=T, querylabel="query", reflabel="reference", customcol=NULL, mappinglabel=paste0(querylabel,"_to_",reflabel)){
  ###if query or ref do not have an SCT refAssay, compute them
  nameflag=0
  if ((intersect(query %>% colnames, refdataset %>% colnames) %>% length)>0){
    colnames(query)= paste0("query_", colnames(query))
    
    colnames(refdataset)=paste0("ref_", colnames(refdataset))
    nameflag=1
  }             
  
  cat("Making reference...\n\n\n")
  DefaultAssay(query)="refAssay"
  DefaultAssay(refdataset)="refAssay"
  refdataset <- RunPCA(refdataset, features = VariableFeatures(refdataset), verbose=F)
  refdataset <- RunUMAP(refdataset, reduction = "pca", dims = pdms, , assay="refAssay", return.model=TRUE, verbose=F)
  ref <- AzimuthReference(
    object = refdataset,
    refUMAP = "umap",
    refDR = "pca",
    refAssay = "refAssay",
    metadata = mapvar,
    dims = pdms,
    k.param = 31,
    reference.version = "1.0.0"
  )
  cat("Finding transfer anchors...\n\n\n")
  anchors <- FindTransferAnchors(
    reference = ref,
    query = query,
    k.filter = NA,
    reference.neighbors = "refdr.annoy.neighbors",
    reference.assay = "refAssay",
    query.assay = "refAssay",
    reference.reduction = "refDR",
    normalization.method = "SCT",
    features = intersect(rownames(x = ref), VariableFeatures(object = query)),
    dims = pdms,
    n.trees = 20,
    mapping.score.k = 100
  )
  cat("Transfering data...\n\n\n")
  query <- TransferData(
    reference = ref,
    query = query,
    dims = pdms,
    anchorset = anchors,
    refdata = mapvar,
    n.trees = 20,
    store.weights = TRUE
  )
  
  
  cat("Integrating embeddings...\n\n\n")
  # Calculate the embeddings of the query data on the reference SPCA
  query <- IntegrateEmbeddings(
    anchorset = anchors,
    reference = ref,
    query = query,
    reductions = "pcaproject",
    reuse.weights.matrix = TRUE,
    dims.to.integrate=pdms
  )
  
  cat("Finding query's neighbors in ref...\n\n\n")
  query[["query_ref.nn"]] <- FindNeighbors(
    object = Embeddings(ref[["refDR"]]),
    query = Embeddings(query[["integrated_dr"]]),
    return.neighbor = TRUE,
    l2.norm = TRUE
  )
  
  # The reference used in the app is downsampled compared to the reference on which
  # the UMAP model was computed. This step, using the helper function NNTransform,
  # corrects the Neighbors to account for the downsampling.
  cat("Azimuth nntransform...\n\n\n")
  query <- Azimuth:::NNTransform(
    object = query,
    meta.data = ref[[]]
  )
  
  # Project the query to the reference UMAP.
  cat("projecting query into the reference umap...\n\n\n")
  query[["mapfun.umap"]] <- RunUMAP(
    object = query[["query_ref.nn"]],
    reduction.model = ref[["refUMAP"]],
    reduction.key = 'mapfun.umap_',
    dims=pdms
  )
  
  cat("adding mapping score...\n\n\n")
  query <- AddMetaData(
    object = query,
    metadata = MappingScore(anchors = anchors),
    col.name = "mapping.score"
  )
  
  cat("preparing the merge...\n\n\n")
  query[["origin"]]=querylabel
  refdataset[["origin"]]=reflabel
  
  refdataset=AddMetaData(refdataset, reddf(ref, "refUMAP")[, c("refumap_1","refumap_2")])
  refdataset=AddMetaData(refdataset, metadata(refdataset)[, c(mapvar)], col.name=paste0("mapfun_", mapvar))
  
  query=AddMetaData(query, reddf(query, "mapfun.umap")[, c("mapfunumap_1")], col.name="refumap_1")
  query=AddMetaData(query, reddf(query, "mapfun.umap")[, c("mapfunumap_2")], col.name="refumap_2")
  query=AddMetaData(query, metadata(query)[, c("predicted.id")], col.name= paste0("mapfun_",mapvar))
  query=AddMetaData(query, metadata(query)[, c("predicted.id.score")], col.name= paste0("mapfun_",mapvar, "_predicted.id.score"))
  query=AddMetaData(query, metadata(query)[, c("mapping.score")], col.name= paste0("mapfun_",mapvar, "_mapping.score"))
  
  nm1=paste0("mapfun_",mapvar)
  nm2=paste0("mapfun_",mapvar, "_predicted.id.score")
  nm3=paste0("mapfun_",mapvar, "_mapping.score")
  
  
  expandcols=function(x){
    x[,nm2]=NA
    x[,nm3]=NA
    
    x}
  
  
  if(!is.na(customcol)){
    getcols=function(x) x[, c("origin", "refumap_1", "refumap_2",nm1, nm2, nm3, customcol)]
  }else{
    getcols=function(x) x[, c("origin", "refumap_1", "refumap_2",nm1, nm2, nm3)  ]
  }
  if(return.merge.metadata){
    
    
    #B %>% left_join(A )
    rbind(query %>% metadata %>% getcols , refdataset %>% metadata %>% expandcols %>% getcols)
    
  }else{
    if (return.query.metadata){ 
      if (nameflag==1){
        colnames(query)= sapply(strsplit(colnames(query), split="query_"), function(x) x[2])
      }
      
      
      if(return.lean){
        
        return(
          
          query %>% metadata  %>% getcols
          
        )}else{ return(query %>% metadata)}
      
      
      
    }
  }
}



maptoref2= function(query, refdataset, ref.assay="refAssay", query.assay="refAssay",  mapvar, normalization.method="SCT", reference.pca="pca", reference.umap="umap", project.query=FALSE,n.trees=20, pdms=1:50, return.merge.metadata=T, return.query.metadata=F, return.lean=T, querylabel="query", reflabel="reference", customcol=NA, mappinglabel=paste0(querylabel,"_to_",reflabel), process.reference=F){
  ###if query or ref do not have an SCT refAssay, compute them
  nameflag=0
  if ((intersect(query %>% colnames, refdataset %>% colnames) %>% length)>0){
    colnames(query)= paste0("query_", colnames(query))
    
    colnames(refdataset)=paste0("ref_", colnames(refdataset))
    nameflag=1
  }             
  
  
  cat("Making reference...\n\n\n")
  DefaultAssay(query)=query.assay
  
  if(process.reference==T){ # this assumes that the assay exists
    
    if(!(ref.assay %in% refdataset@assays %>% names)){
      error(paste0("Reference assay", ref.assay, " does not exist. Make sure to include this assay before running."))
    }else{
      
      DefaultAssay(refdataset)=ref.assay  
      refdataset <- RunPCA(refdataset, features = VariableFeatures(refdataset), verbose=F, reduction.name="mappingpca")
      refdataset <- RunUMAP(refdataset, reduction = "mappingpca", dims = pdms, , assay=ref.assay, return.model=TRUE, verbose=F)
    }
  }
  
  cat("Finding transfer anchors...\n\n\n")
  
  
  anchors <- FindTransferAnchors(
    reference = refdataset,
    query = query,
    k.filter = NA,
    reduction="pcaproject", #
    reference.neighbors = NULL,
    reference.assay = ref.assay,
    query.assay = query.assay,
    reference.reduction = reference.reduction,
    project.query=project.query,
    normalization.method = normalization.method,
    features = intersect(rownames(x = refdataset), VariableFeatures(object = query)),
    dims = pdms,
    n.trees = n.trees,
    mapping.score.k = 100
  )
  cat("Transfering data...\n\n\n")
  query <- TransferData(
    reference = refdataset,
    query = query,
    dims = pdms,
    anchorset = anchors,
    refdata = mapvar,
    n.trees = n.trees,
    store.weights = TRUE
  )
  
  
  cat("Integrating embeddings...\n\n\n")
  # Calculate the embeddings of the query data on the reference SPCA
  query <- IntegrateEmbeddings(
    anchorset = anchors,
    reference = refdataset,
    query = query,
    reductions = "pcaproject",
    reuse.weights.matrix = TRUE,
    dims.to.integrate=pdms,
    new.reduction.name="integrated_dr"
  )
  
  cat("Finding query's neighbors in ref...\n\n\n")
  query[["query_ref.nn"]] <- FindNeighbors(
    object = Embeddings(refdataset[["mappingpca"]]),
    query = Embeddings(query[["integrated_dr"]]),
    return.neighbor = TRUE,
    l2.norm = TRUE
  )
  
  # The reference used in the app is downsampled compared to the reference on which
  # the UMAP model was computed. This step, using the helper function NNTransform,
  # corrects the Neighbors to account for the downsampling.
  cat("Azimuth nntransform...\n\n\n")
  query <- Azimuth:::NNTransform(
    object = query,
    meta.data = refdataset[[]]
  )
  
  # Project the query to the reference UMAP.
  cat("projecting query into the reference umap...\n\n\n")
  query[["mapfun.umap"]] <- RunUMAP(
    object = query[["query_ref.nn"]],
    reduction.model = refdataset[["refUMAP"]],
    reduction.key = 'mapfun.umap_',
    dims=pdms
  )
  
  cat("adding mapping score...\n\n\n")
  query <- AddMetaData(
    object = query,
    metadata = MappingScore(anchors = anchors),
    col.name = "mapping.score"
  )
  
  cat("preparing the merge...\n\n\n")
  query[["origin"]]=querylabel
  refdataset[["origin"]]=reflabel
  
  refdataset=AddMetaData(refdataset, reddf(ref, "refUMAP")[, c("refumap_1","refumap_2")])
  refdataset=AddMetaData(refdataset, metadata(refdataset)[, c(mapvar)], col.name=paste0("mapfun_", mapvar))
  
  query=AddMetaData(query, reddf(query, "mapfun.umap")[, c("mapfunumap_1")], col.name="refumap_1")
  query=AddMetaData(query, reddf(query, "mapfun.umap")[, c("mapfunumap_2")], col.name="refumap_2")
  query=AddMetaData(query, metadata(query)[, c("predicted.id")], col.name= paste0("mapfun_",mapvar))
  query=AddMetaData(query, metadata(query)[, c("predicted.id.score")], col.name= paste0("mapfun_",mapvar, "_predicted.id.score"))
  query=AddMetaData(query, metadata(query)[, c("mapping.score")], col.name= paste0("mapfun_",mapvar, "_mapping.score"))
  
  nm1=paste0("mapfun_",mapvar)
  nm2=paste0("mapfun_",mapvar, "_predicted.id.score")
  nm3=paste0("mapfun_",mapvar, "_mapping.score")
  
  
  expandcols=function(x){
    x[,nm2]=NA
    x[,nm3]=NA
    
    x}
  
  
  if(!is.na(customcol)){
    getcols=function(x) x[, c("origin", "refumap_1", "refumap_2",nm1, nm2, nm3, customcol)]
  }else{
    getcols=function(x) x[, c("origin", "refumap_1", "refumap_2",nm1, nm2, nm3)  ]
  }
  if(return.merge.metadata){
    
    
    #B %>% left_join(A )
    rbind(query %>% metadata %>% getcols , refdataset %>% metadata %>% expandcols %>% getcols)
    
  }else{
    if (return.query.metadata){ 
      if (nameflag==1){
        colnames(query)= sapply(strsplit(colnames(query), split="query_"), function(x) x[2])
      }
      
      
      if(return.lean){
        
        return(
          
          query %>% metadata  %>% getcols
          
        )}else{ return(query %>% metadata)}
      
      
      
    }
  }
}



maptoref.pt1= function(query, 
                       refdataset, 
                       ref.assay="refAssay",
                       query.assay="refAssay",
                       mapvar, 
                       normalization.method="SCT",
                       reference.pca="pca", 
                       reference.umap="umap",
                       reference.reduction=reference.pca,
                       reference.neighbors=NULL,
                       project.query=FALSE,
                       n.trees=20,
                       pdms=1:50,
                       return.merge.metadata=T, 
                       return.query.metadata=F,
                       return.lean=T,
                       querylabel="query", 
                       reflabel="reference", 
                       customcol=NULL, 
                       mappinglabel=paste0(querylabel,"_to_",reflabel),
                       process.reference=F){
  ###if query or ref do not have an SCT refAssay, compute them
  
  nameflag=0
  if ((intersect(colnames(query), colnames(refdataset)) %>% length )>0){
    msg("labelling cell ids according to dataset of origin...")
    query=RenameCells(query, new.name= paste0("query_", colnames(query)))
    
    refdataset=RenameCells(refdataset, new.name= paste0("ref_", colnames(refdataset)))
    nameflag=1
  }             
  fcat("finished labeling cell ids")
  
  cat("Making reference...\n\n\n")
  DefaultAssay(query)=query.assay
  
  if(process.reference==T){ # this assumes that the assay exists
    
    if(!(ref.assay %in% refdataset@assays %>% names)){
      error(paste0("Reference assay", ref.assay, " does not exist. Make sure to include this assay before reprocessing."))
    }else{
      
      DefaultAssay(refdataset)=ref.assay  
      refdataset <- RunPCA(refdataset, features = VariableFeatures(refdataset), verbose=F, reduction.name="mappingpca")
      refdataset <- RunUMAP(refdataset, reduction = "mappingpca", dims = pdms, , assay=ref.assay, return.model=TRUE, verbose=F)
    }
  }else{
    fcat("setting", reference.pca, "and", reference.umap, "as mapping reductions")
    refdataset@reductions$mappingpca= refdataset@reductions[[reference.pca]]
    refdataset@reductions$mappingumap= refdataset@reductions[[reference.umap]]
  }
  
  # if(!is.null(reference.neighbors){
  #  nnname=paste0(umapred, "_snn")
  #  fcat("Attempting to use", nnname, "as reference neighbors...")
  # reference.neighbors=nnname
  #}
  
  if(is.null(reference.reduction)){
    fcat("Setting", reference.pca, "as the reference reduction")
    reference.reduction=reference.pca
  }
  
  
  cat("Finding transfer anchors...\n\n\n")
  
  
  anchors <- FindTransferAnchors(
    reference = refdataset,
    query = query,
    k.filter = NA,
    reduction="pcaproject", # method of projection
    reference.neighbors = reference.neighbors,
    reference.assay = ref.assay,
    query.assay = query.assay,
    reference.reduction = reference.reduction,
    project.query=project.query,
    normalization.method = normalization.method,
    features = intersect(rownames(x = refdataset), VariableFeatures(object = query)),
    dims = pdms,
    n.trees = n.trees,
    mapping.score.k = 100
  )
  cat("Transfering data...\n\n\n")
  query <- TransferData(
    reference = refdataset,
    query = query,
    dims = pdms,
    anchorset = anchors,
    refdata = mapvar,
    n.trees = n.trees,
    store.weights = TRUE
  )
  list(query=query, anchors=anchors, ref=refdataset)
  
  
  
}





maptoref.pt2= function(query, refdataset, anchors, ref.assay="refAssay", query.assay="refAssay",  mapvar, normalization.method="SCT", reference.pca="pca", reference.umap="umap",reference.reduction=reference.pca, reference.neighbors=NULL, project.query=FALSE,n.trees=20, pdms=1:50, return.merge.metadata=T, return.query.metadata=F, return.lean=T, querylabel="query", reflabel="reference", customcol=NULL, mappinglabel=paste0(querylabel,"_to_",reflabel), process.reference=F){
  
  
  query=Seurat::MapQuery(
    anchor=anchors,
    query=query,
    reference=refdataset,
    refdata = mapvar,
    new.reduction.name = "mapfun.pca",
    reference.reduction = reference.reduction,
    reference.dims = pdms,
    query.dims = pdms,
    reduction.model = reference.umap,
    transferdata.args = list(),
    integrateembeddings.args = list(),
    projectumap.args = list(
      reduction.name=reference.umap,
      reduction.key= prepare.rk(reference.umap)
    ),
    verbose = TRUE
  )   
  
  
  
  cat("adding mapping score...\n\n\n")
  query <- AddMetaData(
    object = query,
    metadata = MappingScore(anchors = anchors),
    col.name = "mapping.score"
  )    
  
  
  cat("preparing the merge...\n\n\n")
  query[["origin"]]=querylabel
  refumap1=prepare.rk(reference.umap) %>% paste0(., "1")
  refumap2=prepare.rk(reference.umap) %>% paste0(., "2")
  
  #refdataset=AddMetaData(refdataset, reddf(ref, "refUMAP")[, c("refumap_1","refumap_2")])
  refdataset=AddMetaData(refdataset, metadata(refdataset)[, c(mapvar)], col.name=paste0("mapfun_", mapvar))
  
  query=AddMetaData(query, reddf(query, reference.umap)[, refumap1], col.name=refumap1)
  query=AddMetaData(query, reddf(query, reference.umap)[, refumap2], col.name=refumap2)
  query=AddMetaData(query, metadata(query)[, c("predicted.id")], col.name= paste0("mapfun_",mapvar))
  query=AddMetaData(query, metadata(query)[, c("predicted.id.score")], col.name= paste0("mapfun_",mapvar, "_predicted.id.score"))
  query=AddMetaData(query, metadata(query)[, c("mapping.score")], col.name= paste0("mapfun_",mapvar, "_mapping.score"))
  
  nm1=paste0("mapfun_",mapvar)
  nm2=paste0("mapfun_",mapvar, "_predicted.id.score")
  nm3=paste0("mapfun_",mapvar, "_mapping.score")
  
  
  expandcols=function(x){
    x[,nm2]=NA
    x[,nm3]=NA
    
    x}
  
  
  if(!is.null(customcol)){
    getcols=function(x) x[, c("origin", refumap1, refumap2,nm1, nm2, nm3, customcol)]
  }else{
    getcols=function(x) x[, c("origin", refumap1, refumap2,nm1, nm2, nm3)  ]
  }
  if(return.merge.metadata){
    
    
    #B %>% left_join(A )
    rbind(query %>% metadata %>% getcols , refdataset %>% metadata %>% expandcols %>% getcols)
    
  }else{
    if (return.query.metadata){ 
      #if (nameflag==1){
      # colnames(query)= sapply(strsplit(colnames(query), split="query_"), function(x) x[2])
      #}
      
      
      if(return.lean){
        
        return(
          
          query %>% metadata  %>% getcols
          
        )}else{ return(query %>% metadata)}
      
      
      
    }
  }   
  
  
  
  
}



maptoref.workflow=function(query, refdataset, ...){
  #function(query, refdataset, ref.assay="refAssay", query.assay="refAssay",  mapvar, normalization.method="SCT", reference.pca="pca", reference.umap="umap", project.query=FALSE,n.trees=20, pdms=1:50, return.merge.metadata=T, return.query.metadata=F, return.lean=T, querylabel="query", reflabel="reference", customcol=NA, mappinglabel=paste0(querylabel,"_to_",reflabel), process.reference=F)
  mappedcells=maptoref.pt1(query, refdataset, ...)
  mapmetadata=maptoref.pt2(query=mappedcells$query, anchors=mappedcells$anchors, refdataset=mappedcells$ref, ...)
  mapmetadata
  
}


labeltransfer= function(query, 
                        refdataset, 
                        ref.assay="refAssay",
                        query.assay="refAssay",
                        mapvar, 
                        normalization.method="SCT",
                        reference.pca="pca", 
                        reference.umap="umap",
                        reference.reduction=reference.pca,
                        reference.neighbors=NULL,
                        project.query=FALSE,
                        n.trees=20,
                        pdms=1:50,
                        return.lean=T,
                        querylabel="query", 
                        reflabel="reference", 
                        customcol=NULL, 
                        mappinglabel=paste0(querylabel,"_to_",reflabel),
                        process.reference=F){
  ###if query or ref do not have an SCT refAssay, compute them
  
  nameflag=0
  if ((intersect(colnames(query), colnames(refdataset)) %>% length )>0){
    msg("labelling cell ids according to dataset of origin...")
    query=RenameCells(query, new.name= paste0("query_", colnames(query)))
    
    refdataset=RenameCells(refdataset, new.name= paste0("ref_", colnames(refdataset)))
    nameflag=1
  }             
  fcat("finished labeling cell ids")
  
  cat("Making reference...\n\n\n")
  DefaultAssay(query)=query.assay
  
  if(process.reference==T){ # this assumes that the assay exists
    
    if(!(ref.assay %in% refdataset@assays %>% names)){
      error(paste0("Reference assay", ref.assay, " does not exist. Make sure to include this assay before reprocessing."))
    }else{
      
      DefaultAssay(refdataset)=ref.assay  
      refdataset <- RunPCA(refdataset, features = VariableFeatures(refdataset), verbose=F, reduction.name="mappingpca")
      refdataset <- RunUMAP(refdataset, reduction = "mappingpca", dims = pdms, , assay=ref.assay, return.model=TRUE, verbose=F)
    }
  }else{
    fcat("setting", reference.pca, "and", reference.umap, "as mapping reductions")
    refdataset@reductions$mappingpca= refdataset@reductions[[reference.pca]]
    refdataset@reductions$mappingumap= refdataset@reductions[[reference.umap]]
  }
  
  # if(!is.null(reference.neighbors){
  #  nnname=paste0(umapred, "_snn")
  #  fcat("Attempting to use", nnname, "as reference neighbors...")
  # reference.neighbors=nnname
  #}
  
  if(is.null(reference.reduction)){
    fcat("Setting", reference.pca, "as the reference reduction")
    reference.reduction=reference.pca
  }
  
  
  cat("Finding transfer anchors...\n\n\n")
  
  
  anchors <- FindTransferAnchors(
    reference = refdataset,
    query = query,
    k.filter = NA,
    reduction="pcaproject", # method of projection
    reference.neighbors = reference.neighbors,
    reference.assay = ref.assay,
    query.assay = query.assay,
    reference.reduction = reference.reduction,
    project.query=project.query,
    normalization.method = normalization.method,
    features = intersect(rownames(x = refdataset), VariableFeatures(object = query)),
    dims = pdms,
    n.trees = n.trees,
    mapping.score.k = 100
  )
  cat("Transfering data...\n\n\n")
  query <- TransferData(
    reference = refdataset,
    query = query,
    dims = pdms,
    anchorset = anchors,
    refdata = mapvar,
    n.trees = n.trees,
    store.weights = TRUE
  )
  #list(query=query, anchors=anchors, ref=refdataset)
  
  
  cat("adding mapping score...\n\n\n")
  query <- AddMetaData(
    object = query,
    metadata = MappingScore(anchors = anchors),
    col.name = "mapping.score"
  )    
  
  
  cat("preparing the merge...\n\n\n")
  query[["origin"]]=querylabel
  
  #refdataset=AddMetaData(refdataset, reddf(ref, "refUMAP")[, c("refumap_1","refumap_2")])
  refdataset=AddMetaData(refdataset, metadata(refdataset)[, c(mapvar)], col.name=paste0("mapfun_", mapvar))
  
  query=AddMetaData(query, metadata(query)[, c("predicted.id")], col.name= paste0("mapfun_",mapvar))
  query=AddMetaData(query, metadata(query)[, c("predicted.id.score")], col.name= paste0("mapfun_",mapvar, "_predicted.id.score"))
  query=AddMetaData(query, metadata(query)[, c("mapping.score")], col.name= paste0("mapfun_",mapvar, "_mapping.score"))
  
  nm1=paste0("mapfun_",mapvar)
  nm2=paste0("mapfun_",mapvar, "_predicted.id.score")
  nm3=paste0("mapfun_",mapvar, "_mapping.score")
  
  
  expandcols=function(x){
    x[,nm2]=NA
    x[,nm3]=NA
    
    x}
  
  
  if(!is.null(customcol)){
    getcols=function(x) x[, c("origin", nm1, nm2, nm3, customcol)]
  }else{
    getcols=function(x) x[, c("origin",nm1, nm2, nm3)  ]
  }
  
  
  if(return.lean){
    
    return(
      
      query %>% metadata  %>% getcols
      
    )}else{ return(query %>% metadata)}
  
  
  
}   












maptoref.pt3= function(query, refdataset, anchors, ref.assay="refAssay", query.assay="refAssay",  mapvar, normalization.method="SCT", reference.pca="pca", reference.umap="umap", project.query=FALSE,n.trees=20, pdms=1:50, return.merge.metadata=T, return.query.metadata=F, return.lean=T, querylabel="query", reflabel="reference", customcol=NA, mappinglabel=paste0(querylabel,"_to_",reflabel), process.reference=F){
  
  cat("Finding query's neighbors in ref...\n\n\n")
  query[["query_ref.nn"]] <- FindNeighbors(
    object = Embeddings(refdataset[[reference.pca]]),
    query = Embeddings(query[["integrated_dr"]]),
    return.neighbor = TRUE,
    l2.norm = TRUE
  )
  
  
  
  
  # The reference used in the app is downsampled compared to the reference on which
  # the UMAP model was computed. This step, using the helper function NNTransform,
  # corrects the Neighbors to account for the downsampling.
  cat("Azimuth nntransform...\n\n\n")
  query <- Azimuth:::NNTransform(
    object = query,
    meta.data = refdataset[[]]
  )
  
  # Project the query to the reference UMAP.
  cat("projecting query into the reference umap...\n\n\n")
  query[["mapfun.umap"]] <- RunUMAP(
    object = query[["query_ref.nn"]],
    reduction.model = refdataset[["refUMAP"]],
    reduction.key = 'mapfun.umap_',
    dims=pdms
  )
  
  cat("adding mapping score...\n\n\n")
  query <- AddMetaData(
    object = query,
    metadata = MappingScore(anchors = anchors),
    col.name = "mapping.score"
  )
  
  cat("preparing the merge...\n\n\n")
  query[["origin"]]=querylabel
  refdataset[["origin"]]=reflabel
  
  refdataset=AddMetaData(refdataset, reddf(ref, "refUMAP")[, c("refumap_1","refumap_2")])
  refdataset=AddMetaData(refdataset, metadata(refdataset)[, c(mapvar)], col.name=paste0("mapfun_", mapvar))
  
  query=AddMetaData(query, reddf(query, "mapfun.umap")[, c("mapfunumap_1")], col.name="refumap_1")
  query=AddMetaData(query, reddf(query, "mapfun.umap")[, c("mapfunumap_2")], col.name="refumap_2")
  query=AddMetaData(query, metadata(query)[, c("predicted.id")], col.name= paste0("mapfun_",mapvar))
  query=AddMetaData(query, metadata(query)[, c("predicted.id.score")], col.name= paste0("mapfun_",mapvar, "_predicted.id.score"))
  query=AddMetaData(query, metadata(query)[, c("mapping.score")], col.name= paste0("mapfun_",mapvar, "_mapping.score"))
  
  nm1=paste0("mapfun_",mapvar)
  nm2=paste0("mapfun_",mapvar, "_predicted.id.score")
  nm3=paste0("mapfun_",mapvar, "_mapping.score")
  
  
  expandcols=function(x){
    x[,nm2]=NA
    x[,nm3]=NA
    
    x}
  
  
  if(!is.null(customcol)){
    getcols=function(x) x[, c("origin",nm1, nm2, nm3, customcol)]
  }else{
    getcols=function(x) x[, c("origin",nm1, nm2, nm3)  ]
  }
  if(return.merge.metadata){
    
    
    #B %>% left_join(A )
    rbind(query %>% metadata %>% getcols , refdataset %>% metadata %>% expandcols %>% getcols)
    
  }else{
    if (return.query.metadata){ 
      if (nameflag==1){
        colnames(query)= sapply(strsplit(colnames(query), split="query_"), function(x) x[2])
      }
      
      
      if(return.lean){
        
        return(
          
          query %>% metadata  %>% getcols
          
        )}else{ return(query %>% metadata)}
      
      
      
    }
  }
}

#######################################################################################################
#Function to generate hulls from umap clusters that light up according to some metric
#######################################################################################################
clusterhullfill=function(so, umap.df= NULL, clusvar, clus=NULL, linecolor="white" ,fillvar=NULL, fillcols=NULL, size=0.05){
  # if squares...
  #maxx=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% pull(UMAP_1) %>% max
  #maxy=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>%pull(UMAP_2) %>% max
  #minx=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% pull(UMAP_1) %>% min
  #miny=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% pull(UMAP_2) %>% min
  
  if(is.null(umap.df)){
    umap.df= so %>% umapdf
    
  }
  
  if(is.null(clus)){
    clus=umap.df %>% pull(clusvar) %>% unique %>% sort
    
  }
  
  catpols=function(cc, pols=NULL){ 
    #if(is.null(fillvar)){
    cl=umap.df %>% filter(!!sym(clusvar)==clus[cc]) %>% select(UMAP_1, UMAP_2, !!sym(fillvar))
    #fillvar=NA
    #}else{
    #cl=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% select(UMAP_1, UMAP_2, fillvar)
    #}
    
    clhull=cl[chull(cl), ]
    
    if(is.null(pols)){
      pols<-list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size))
    }else{
      pols<-append(pols, geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size))
    }
    
    if(cc==length(clus)){
      return(pols)
    }else{
      return(catpols(cc+1, pols))}
    
    
  }
  
  catpols(1)
}



################################################################################
#
################################################################################

clusterhullfill2=function(so, umap.df= NULL, clusvar, clus=NULL, linecolor="white" ,fillvar=NULL, fillcols=NULL, size=0.05, concavity=2){
  # if squares...
  #maxx=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% pull(UMAP_1) %>% max
  #maxy=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>%pull(UMAP_2) %>% max
  #minx=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% pull(UMAP_1) %>% min
  #miny=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% pull(UMAP_2) %>% min
  
  if(is.null(umap.df)){
    umap.df= so %>% umapdf
    
  }
  
  if(is.null(clus)){
    clus=umap.df %>% pull(clusvar) %>% unique %>% sort
    
  }
  
  catpols=function(cc, pols=NULL){ 
    #if(is.null(fillvar)){
    cl=umap.df %>% filter(!!sym(clusvar)==clus[cc]) %>% select(UMAP_1, UMAP_2, !!sym(fillvar))
    #fillvar=NA
    #}else{
    #cl=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% select(UMAP_1, UMAP_2, fillvar)
    #}
    
    #clhull=cl[chull(cl), ]
    clhull=concaveman(cl, concavity=concavity)
    
    if(is.null(pols)){
      pols<-list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size))
    }else{
      pols<-append(pols, geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size))
    }
    
    if(cc==length(clus)){
      return(pols)
    }else{
      return(catpols(cc+1, pols))}
    
    
  }
  
  catpols(1)
}

#############################################################################
#Code for seurat utilities which are useful
#############################################################################
NoLegend=function (...)
{
  no.legend.theme <- theme(legend.position = "none", validate = TRUE,
                           ...)
  return(no.legend.theme)
}


NoAxes=function (..., keep.text = FALSE, keep.ticks = FALSE)
{
  blank <- element_blank()
  no.axes.theme <- theme(axis.line.x = blank, axis.line.y = blank,
                         validate = TRUE, ...)
  if (!keep.text) {
    no.axes.theme <- no.axes.theme + theme(axis.text.x = blank,
                                           axis.text.y = blank, axis.title.x = blank, axis.title.y = blank,
                                           validate = TRUE, ...)
  }
  if (!keep.ticks) {
    no.axes.theme <- no.axes.theme + theme(axis.ticks.x = blank,
                                           axis.ticks.y = blank, validate = TRUE, ...)
  }
  return(no.axes.theme)
}



#######################importing manually functions to run hull sample from mvgps

hull_sample=function (X, num_grid_pts = 500, grid_type = "regular", 
                      trim_hull = FALSE, trim_quantile = NULL) 
{
  X_rslt <- X_check(X)
  assign("X", X_rslt$X)
  assign("m", X_rslt$m)
  grid_type <- match.arg(grid_type, choices = c("regular", 
                                                "random", "hexagonal"))
  if (trim_hull == TRUE) {
    if (is.null(trim_quantile)) 
      stop("trim_hull set to TRUE but trim_quantile not specified.", 
           call. = FALSE)
    if (trim_quantile < 0.5 | trim_quantile > 1) 
      stop("trim_quantile must be between [0.5, 1]", 
           call. = FALSE)
    trim_upper <- apply(X, 2, quantile, trim_quantile)
    trim_lower <- apply(X, 2, quantile, 1 - trim_quantile)
    X_trim <- sapply(seq_len(m), function(x) {
      ifelse(X[, x] > trim_upper[x], NA, ifelse(X[, x] < 
                                                  trim_lower[x], NA, X[, x]))
    })
    colnames(X_trim) <- colnames(X)
    X <- na.omit(X_trim)
  }
  if (m == 2) {
    hpts <- chull(X)
    hpts <- c(hpts, hpts[1])
    hpts_vs <- as.matrix(X[hpts, ])
    m <- Polygon(hpts_vs)
    ps <- Polygons(list(m), 1)
    sps <- SpatialPolygons(list(ps))
    sp_grid_pts <- spsample(sps, n = num_grid_pts, type = grid_type)
    grid_pts <- coordinates(sp_grid_pts)
    colnames(grid_pts) <- colnames(X)
  }
  else {
    hpts <- geometry::convhulln(X)
    hpts_ind <- unique(c(hpts))
    hpts_vs <- X[hpts_ind, ]
    grid_pts <- NULL
  }
  return(list(hpts_vs = hpts_vs, grid_pts = grid_pts, X = X))
}



X_check=function (X) 
{
  X <- as.matrix(X)
  m <- ncol(X)
  if (!any(apply(X, 2, is.numeric))) 
    stop("X must be numeric", call. = FALSE)
  if (m < 2) 
    stop("Exposure is not multivariate", call. = FALSE)
  return(list(X = X, m = m))
}


####################################################################
#function to make cluster hulls filled with points that signal a correspondence of cells, with points jittered to minimise overplotting effects
####################################################################


clusterhullpoint=function(so, umap.df= NULL, query.df, clusvar, clus, linecolor="white" ,fillvar=NULL, colorpoints="stage", fillcols=NULL, size=0.05, pointsize=0.05, reduction.key="UMAP_", concavity=2, downsample=1, rescale.cells=NULL){
  require(sf)
  
  
  if(is.null(umap.df)){
    umap.df= so %>% umapdf %>% dplyr::mutate(UMAP_1= !!sym(paste0(reduction.key, "1")), UMAP_2=!!sym(paste0(reduction.key, "2")))
    
  }else{
    umap.df= umap.df  %>% dplyr::mutate(UMAP_1= !!sym(paste0(reduction.key, "1")), UMAP_2=!!sym(paste0(reduction.key, "2")))
  }
  
  
  destinyclusters=query.df %>% group_by( !!sym(mapfun(clusvar))) %>% summarise(counts=n()) %>% as.data.frame %>% arrange(-counts)
  
  ##collect all the destinies of each cell in each  together
  
  if(is.null(rescale.cells)){
    des=destinyclusters %>% dplyr::mutate(temptot= sum(counts), frac=counts/temptot, newcounts= ceiling(frac*temptot) ) %>% col2names(., mapfun(clusvar))
  }else{
    des=destinyclusters %>% dplyr::mutate(temptot= sum(counts), newtot=!!rescale.cells, frac=counts/temptot, newcounts= ceiling(frac*newtot) ) %>% col2names(., mapfun(clusvar))
    
  }
  
  
  fcat("printing normal des:")
  print(destinyclusters)
  if(downsample!=1){
    fcat("printing downsampled des:")
    des= des %>% dplyr::mutate(counts=ceiling(counts*downsample))
    print(des)
  }
  getclusterfrac=function(x) des[x, "frac"]
  umap.df[["clusterfrac"]]=NA
  umap.df[["clusterfrac"]]=sapply(umap.df %>% pull(clusvar), getclusterfrac) %>% unname
  
  catpols=function(cc, pols=NULL){ 
    #if(is.null(fillvar)){
    cl=umap.df %>% dplyr::mutate(null="one") %>% filter(!!sym(clusvar)==clus[cc]) %>% dplyr::select(UMAP_1, UMAP_2, !!sym(fillvar))
    #fillvar=NA
    #}else{
    #cl=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% select(UMAP_1, UMAP_2, fillvar)
    #}
    
    #clhull=cl[chull(cl), ] #obsolete but keep in case we need to revert to convex hulls
    
    
    
    clhull=concaveman(as.matrix(cl[, c("UMAP_1", "UMAP_2")]), concavity=concavity)
    
    rownames(clhull)= 1:nrow(clhull) %>% addprefix(., "hullpoint_")
    
    clhull=as.data.frame(clhull) %>% givecolnames(., nms=c("UMAP_1", "UMAP_2")) %>% dplyr::mutate(null="one")
    
    
    dists=dist(rbind(clhull[, c("UMAP_1", "UMAP_2")], cl[,c("UMAP_1", "UMAP_2")])) %>% as.matrix
    ld= nrow(dists)
    lc= nrow(clhull)
    inds=lapply(rownames(clhull), function(r){ names(dists[r, (lc+1):ld])[findmin(dists[r,(lc+1):ld])]}) %>% unlist
    clhull[[fillvar]]= cl[inds, ][[fillvar]]
    clhull[[mapfun(varr)]]= cl[inds, ][[mapfun(varr)]]
    clhull[[varr]]= cl[inds, ][[mapfun(varr)]]
    
    
    #print(clhull)
    ##getting number of points and sampling the hull space uniformly for the number of points found there
    fcat("number of cells in", clus[cc], ":", query.df %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% nrow())
    #at this point downsampling should already have been taken care of
    numcells=des %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% pull(newcounts)
    originalcells=des %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% pull(counts)
    if(length(numcells)==0){
      numcells=0
    }
    #numcells=query.df %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% nrow()
    subsample=(numcells*downsample) %>% round
    #fcat("numcells is" , numcells)
    if(!is.null(colorpoints) & !is.na(numcells) & numcells>2){
      #hsamples=hull_sample(X=clhull[, c(1,2)] %>% as.matrix, query.df %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% nrow() , grid_type="random")$grid_pts
      
      polygon = sf::st_polygon(list(clhull[, c("UMAP_1","UMAP_2")] %>% data.matrix))
      
      # Sample 50 random points within the polygon
      hsamples = sf::st_sample(polygon, size=numcells)
      hsamples= sf::st_coordinates(hsamples) %>% givecolnames(., nms=c("UMAP_1","UMAP_2"))
      
      
      fcat("samples computed...")
      fcat("number of cells in hull sample for", clus[cc], ":", hsamples %>% nrow())
      #print(hsamples)
      ###making matrix of mock coordinates for each cell
      
      if(subsample==originalcells){
        rplce=F}else{rplce=T}
      pointsdf=query.df %>% dplyr::mutate(null=NA) %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% sample_n(subsample, replace=rplce) %>% dplyr::mutate(UMAP_1=hsamples[,1], UMAP_2=hsamples[,2]) 
      
      pointsdf[[varr]] = pointsdf[[mapfun(varr)]]
      
      if(is.null(pols)){
        fcat("creating poygon list...")
        pols<-list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size),
                   geom_point_rast(data=pointsdf, inherit.aes=T, size=pointsize,aes(x=UMAP_1, y=UMAP_2, color=!!sym(colorpoints))))
      }else{
        fcat("appending poygon to list...")
        pols<-append(pols, list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size),
                                geom_point_rast(data=pointsdf,inherit.aes=T,size=pointsize, aes(x=UMAP_1, y=UMAP_2, color=!!sym(colorpoints)))))
      }
      
      if(cc==length(clus)){
        fcat("finishing up. Attaching color scheme...")
        return(append(pols, scale_color_manual(values=allcolors[[colorpoints]])))
        #return(pols)
      }else{
        return(catpols(cc+1, pols))}
      
    }else{###if there are no query cells in cluster
      
      
      if(is.null(pols)){
        fcat("creating poygon list...")
        pols<-list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size))
      }else{
        fcat("appending poygon to list...")
        pols<-append(pols, list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size)))
      }
      
      if(cc==length(clus)){
        #msg("finishing up. Attaching color scheme...")
        fcat("finishing up. Attaching color scheme...")
        return(append(pols, scale_color_manual(values=allcolors[[colorpoints]])))
      }else{
        return(catpols(cc+1, pols))}
      
      
      
      
      
      
    }
    
    
    
  }
  
  catpols(1)
}



####just make cluster hulls
clusterhull=function(so, umap.df= NULL, clusvar, clus=NULL, linecolor="white" ,fillvar="null", colorpoints="stage", fillcols=NULL, size=0.05, pointsize=0.05, reduction.key="UMAP_"){
  require(sp)
  
  
  if(is.null(umap.df)){
    umap.df= so %>% umapdf     
  }
  
  umap.df= umap.df %>% dplyr::mutate(null=NA, UMAP_1= !!sym(paste0(reduction.key, "1")), UMAP_2=!!sym(paste0(reduction.key, "2")))
  
  
  if(is.null(clus)){
    clus= umap.df %>% pull(clusvar) %>% unique %>% sort
    
  }
  
  if( is.null(allcolors[[fillvar]])){
    
    allcolors[[fillvar]]= randomcolors((umap.df %>% pull(fillvar) %>% unique %>% length) )
  }
  
  catpols=function(cc, pols=NULL){ 
    #if(is.null(fillvar)){
    cl=umap.df %>% dplyr::mutate(null="one") %>% filter(!!sym(clusvar)==clus[cc]) %>% select(UMAP_1, UMAP_2, !!sym(fillvar))
    #fillvar=NA
    #}else{
    #cl=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% select(UMAP_1, UMAP_2, fillvar)
    #}
    
    clhull=cl[chull(cl), ]
    
    #print(clhull)
    ##getting number of points and sampling the hull space uniformly for the number of points found there
    #fcat("number of cells in", clus[cc], ":", query.df %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% nrow())
    #numcells=0
    #numcells=query.df %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% nrow()
    
    
    #if(!is.null(colorpoints) && numcells>2){
    #hsamples=hull_sample(X=clhull[, c(1,2)] %>% as.matrix, query.df %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% nrow() , grid_type="random")$grid_pts
    
    #fcat("samples computed...")
    #fcat("number of cells in hull sample for", clus[cc], ":", hsamples %>% nrow())
    #print(hsamples)
    ###making matrix of mock coordinates for each cell
    #pointsdf=query.df %>% dplyr::mutate(null=NA) %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% dplyr::mutate(UMAP_1=hsamples[,1], UMAP_2=hsamples[,2])
    
    
    if(is.null(pols)){
      fcat("creating poygon list...")
      pols<-list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size))
    }else{
      fcat("appending poygon to list...")
      pols<-append(pols, list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size)))
    }
    
    if(cc==length(clus)){
      fcat("finishing up. Attaching color scheme...")
      return(list(pols
                  #, scale_color_manual(values=allcolors[[colorpoints]])
                  , scale_fill_manual(values=allcolors[[fillvar]]) 
      ) %>% Reduce(append, .)
      )
      #return(pols)
    }else{
      return(catpols(cc+1, pols))}
  }
  catpols(1)
}

################################################################################
#uses concaveman for cluster hulls
################################################################################

clusterhull2=function(so, umap.df= NULL, clusvar, clus=NULL, linecolor="white" ,fillvar="null", colorpoints="stage", fillcols=NULL, size=0.05, pointsize=0.05, reduction.key="UMAP_", concavity=2){
  require(sp)
  
  
  if(is.null(umap.df)){
    umap.df= so %>% umapdf     
  }
  
  umap.df= umap.df %>% dplyr::mutate(null=NA, UMAP_1= !!sym(paste0(reduction.key, "1")), UMAP_2=!!sym(paste0(reduction.key, "2")))
  
  
  if(is.null(clus)){
    clus= umap.df %>% pull(clusvar) %>% unique %>% sort
    
  }
  
  if( is.null(allcolors[[fillvar]])){
    
    allcolors[[fillvar]]= randomcolors((umap.df %>% pull(fillvar) %>% unique %>% length) )
  }
  
  catpols=function(cc, pols=NULL){ 
    #if(is.null(fillvar)){
    cl=umap.df %>% dplyr::mutate(null="one") %>% filter(!!sym(clusvar)==clus[cc]) %>% select(UMAP_1, UMAP_2, !!sym(fillvar))
    #fillvar=NA
    #}else{
    #cl=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% select(UMAP_1, UMAP_2, fillvar)
    #}
    
    #clhull=cl[chull(cl), ]
    
    clhull=concaveman(as.matrix(cl[, c("UMAP_1", "UMAP_2")]), concavity=concavity)
    
    rownames(clhull)= 1:nrow(clhull) %>% addprefix(., "hullpoint_")
    
    clhull=as.data.frame(clhull) %>% givecolnames(., nms=c("UMAP_1", "UMAP_2")) %>% dplyr::mutate(null="one")
    # fillvar on the concave matrix gets assiged the color of the closest point
    
    
    
    
    if(fillvar != "null"){
      dists=dist(rbind(clhull[, c("UMAP_1", "UMAP_2")], cl[,c("UMAP_1", "UMAP_2")])) %>% as.matrix
      ld= nrow(dists)
      lc= nrow(clhull)
      clhull[[fillvar]]= cl[lapply(rownames(clhull), function(r){
        names(dists[r, (lc+1):ld])[findmin(dists[r,(lc+1):ld])]  
        
      }) %>% unlist, ][[fillvar]]
      
    }
    
    
    
    #print(clhull)
    ##getting number of points and sampling the hull space uniformly for the number of points found there
    #fcat("number of cells in", clus[cc], ":", query.df %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% nrow())
    #numcells=0
    #numcells=query.df %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% nrow()
    
    
    #if(!is.null(colorpoints) && numcells>2){
    #hsamples=hull_sample(X=clhull[, c(1,2)] %>% as.matrix, query.df %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% nrow() , grid_type="random")$grid_pts
    
    #fcat("samples computed...")
    #fcat("number of cells in hull sample for", clus[cc], ":", hsamples %>% nrow())
    #print(hsamples)
    ###making matrix of mock coordinates for each cell
    #pointsdf=query.df %>% dplyr::mutate(null=NA) %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% dplyr::mutate(UMAP_1=hsamples[,1], UMAP_2=hsamples[,2])
    
    
    if(is.null(pols)){
      fcat("creating poygon list...")
      pols<-list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size))
    }else{
      fcat("appending poygon to list...")
      pols<-append(pols, list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size)))
    }
    
    if(cc==length(clus)){
      fcat("finishing up. Attaching color scheme...")
      return(list(pols
                  #, scale_color_manual(values=allcolors[[colorpoints]])
                  , scale_fill_manual(values=allcolors[[fillvar]]) 
      ) %>% Reduce(append, .)
      )
      #return(pols)
    }else{
      return(catpols(cc+1, pols))}
  }
  catpols(1)
}


################################################################################
# non recursive version of clusterhull2 to avoid weird R issues
################################################################################


clusterhull3=function(so, umap.df= NULL, clusvar, clus=NULL, linecolor="white" ,fillvar="null", colorpoints="stage", fillcols=NULL, size=0.05, pointsize=0.05, reduction.key="UMAP_", concavity=2){
  require(sp)
  
  
  if(is.null(umap.df)){
    umap.df= so %>% umapdf     
  }
  
  umap.df= umap.df %>% dplyr::mutate(null=NA, UMAP_1= !!sym(paste0(reduction.key, "1")), UMAP_2=!!sym(paste0(reduction.key, "2")))
  
  
  if(is.null(clus)){
    clus= umap.df %>% pull(clusvar) %>% unique %>% sort
    
  }
  
  if( is.null(allcolors[[fillvar]])){
    
    allcolors[[fillvar]]= randomcolors((umap.df %>% pull(fillvar) %>% unique %>% length) )
  }
  
  pols=NULL
  polygon=function (x, y = NULL, density = NULL, angle = 45, border = NULL, 
                    col = NA, lty = par("lty"), ..., fillOddEven = FALSE) 
  {
    ..debug.hatch <- FALSE
    xy <- xy.coords(x, y, setLab = FALSE)
    if (is.numeric(density) && all(is.na(density) | density < 
                                   0)) 
      density <- NULL
    if (!is.null(angle) && !is.null(density)) {
      polygon.onehatch <- function(x, y, x0, y0, xd, yd, ..debug.hatch = FALSE, 
                                   ...) {
        if (..debug.hatch) {
          points(x0, y0)
          arrows(x0, y0, x0 + xd, y0 + yd)
        }
        halfplane <- as.integer(xd * (y - y0) - yd * (x - 
                                                        x0) <= 0)
        cross <- halfplane[-1L] - halfplane[-length(halfplane)]
        does.cross <- cross != 0
        if (!any(does.cross)) 
          return()
        x1 <- x[-length(x)][does.cross]
        y1 <- y[-length(y)][does.cross]
        x2 <- x[-1L][does.cross]
        y2 <- y[-1L][does.cross]
        t <- (((x1 - x0) * (y2 - y1) - (y1 - y0) * (x2 - 
                                                      x1))/(xd * (y2 - y1) - yd * (x2 - x1)))
        o <- order(t)
        tsort <- t[o]
        crossings <- cumsum(cross[does.cross][o])
        if (fillOddEven) 
          crossings <- crossings%%2
        drawline <- crossings != 0
        lx <- x0 + xd * tsort
        ly <- y0 + yd * tsort
        lx1 <- lx[-length(lx)][drawline]
        ly1 <- ly[-length(ly)][drawline]
        lx2 <- lx[-1L][drawline]
        ly2 <- ly[-1L][drawline]
        segments(lx1, ly1, lx2, ly2, ...)
      }
      polygon.fullhatch <- function(x, y, density, angle, ..debug.hatch = FALSE, 
                                    ...) {
        x <- c(x, x[1L])
        y <- c(y, y[1L])
        angle <- angle%%180
        if (par("xlog") || par("ylog")) {
          warning("cannot hatch with logarithmic scale active")
          return()
        }
        usr <- par("usr")
        pin <- par("pin")
        upi <- c(usr[2L] - usr[1L], usr[4L] - usr[3L])/pin
        if (upi[1L] < 0) 
          angle <- 180 - angle
        if (upi[2L] < 0) 
          angle <- 180 - angle
        upi <- abs(upi)
        xd <- cos(angle/180 * pi) * upi[1L]
        yd <- sin(angle/180 * pi) * upi[2L]
        if (angle < 45 || angle > 135) {
          if (angle < 45) {
            first.x <- max(x)
            last.x <- min(x)
          }
          else {
            first.x <- min(x)
            last.x <- max(x)
          }
          y.shift <- upi[2L]/density/abs(cos(angle/180 * 
                                               pi))
          x0 <- 0
          y0 <- floor((min(y) - first.x * yd/xd)/y.shift) * 
            y.shift
          y.end <- max(y) - last.x * yd/xd
          while (y0 < y.end) {
            polygon.onehatch(x, y, x0, y0, xd, yd, ..debug.hatch = ..debug.hatch, 
                             ...)
            y0 <- y0 + y.shift
          }
        }
        else {
          if (angle < 90) {
            first.y <- max(y)
            last.y <- min(y)
          }
          else {
            first.y <- min(y)
            last.y <- max(y)
          }
          x.shift <- upi[1L]/density/abs(sin(angle/180 * 
                                               pi))
          x0 <- floor((min(x) - first.y * xd/yd)/x.shift) * 
            x.shift
          y0 <- 0
          x.end <- max(x) - last.y * xd/yd
          while (x0 < x.end) {
            polygon.onehatch(x, y, x0, y0, xd, yd, ..debug.hatch = ..debug.hatch, 
                             ...)
            x0 <- x0 + x.shift
          }
        }
      }
      if (missing(col) || is.null(col)) {
        col <- par("fg")
      }
      else if (any(is.na(col))) {
        col[is.na(col)] <- par("fg")
      }
      if (is.null(border)) 
        border <- col
      if (is.logical(border)) {
        if (!is.na(border) && border) 
          border <- col
        else border <- NA
      }
      start <- 1
      ends <- c(seq_along(xy$x)[is.na(xy$x) | is.na(xy$y)], 
                length(xy$x) + 1)
      num.polygons <- length(ends)
      col <- rep_len(col, num.polygons)
      if (length(border)) 
        border <- rep_len(border, num.polygons)
      if (length(lty)) 
        lty <- rep_len(lty, num.polygons)
      if (length(density)) 
        density <- rep_len(density, num.polygons)
      angle <- rep_len(angle, num.polygons)
      i <- 1L
      for (end in ends) {
        if (end > start) {
          if (is.null(density) || is.na(density[i]) || 
              density[i] < 0) 
            .External.graphics(C_polygon, xy$x[start:(end - 
                                                        1)], xy$y[start:(end - 1)], col[i], NA, lty[i], 
                               ...)
          else if (density[i] > 0) {
            polygon.fullhatch(xy$x[start:(end - 1)], xy$y[start:(end - 
                                                                   1)], col = col[i], lty = lty[i], density = density[i], 
                              angle = angle[i], ..debug.hatch = ..debug.hatch, 
                              ...)
          }
          i <- i + 1
        }
        start <- end + 1
      }
      .External.graphics(C_polygon, xy$x, xy$y, NA, border, 
                         lty, ...)
    }
    else {
      if (is.logical(border)) {
        if (!is.na(border) && border) 
          border <- par("fg")
        else border <- NA
      }
      .External.graphics(C_polygon, xy$x, xy$y, col, border, 
                         lty, ...)
    }
    invisible()
  }
  for(cc in 1:length(clus)){
    
    cl=umap.df %>% dplyr::mutate(null="one") %>% filter(!!sym(clusvar)==clus[cc]) %>% dplyr::select(UMAP_1, UMAP_2, !!sym(fillvar))
    
    
    clhull=concaveman(as.matrix(cl[, c("UMAP_1", "UMAP_2")]), concavity=concavity)
    
    rownames(clhull)= 1:nrow(clhull) %>% addprefix(., "hullpoint_")
    
    clhull=as.data.frame(clhull) %>% givecolnames(., nms=c("UMAP_1", "UMAP_2")) %>% dplyr::mutate(null="one")
    # fillvar on the concave matrix gets assiged the color of the closest point
    
    if(fillvar != "null"){
      dists=dist(rbind(clhull[, c("UMAP_1", "UMAP_2")], cl[,c("UMAP_1", "UMAP_2")])) %>% as.matrix
      ld= nrow(dists)
      lc= nrow(clhull)
      clhull[[fillvar]]= cl[lapply(rownames(clhull), function(r){
        names(dists[r, (lc+1):ld])[findmin(dists[r,(lc+1):ld])]  
        
      }) %>% unlist, ][[fillvar]]
      
    }
    
    if(is.null(pols)){
      fcat("creating poygon list...")
      pols<-list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size))
    }else{
      fcat("appending poygon to list...")
      pols<-append(pols, list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size)))
    }
    
    if(cc==length(clus)){
      fcat("finishing up. Attaching color scheme...")
      return(list(pols
                  #, scale_color_manual(values=allcolors[[colorpoints]])
                  , scale_fill_manual(values=allcolors[[fillvar]]) 
      ) %>% Reduce(append, .)
      )
      #return(pols)
    }
  } 
  return(pols) 
}




randomcolors= function(n)
  
  lapply(1:n, function(x) paste0("#",sample(c(0:9,LETTERS[1:6]), 1),
                                 sample(c(0:9,LETTERS[1:6]), 1),
                                 sample(c(0:9,LETTERS[1:6]), 1),
                                 sample(c(0:9,LETTERS[1:6]), 1),
                                 sample(c(0:9,LETTERS[1:6]), 1),
                                 sample(c(0:9,LETTERS[1:6]), 1))) %>% Reduce(c, .)




demapfun=function(x) strsplit(x, split="mapfun_")[[1]][2]
mapfun=function(x) x %>% addprefix(., "mapfun_")


#####################################

prefixids.smart= function(so){
  colnames(so@assays$RNA@counts)= paste_(Project(so),  colnames(so@assays$RNA@counts))
  so
}
prefixids= function(so, pref=""){
  if ( !("patientID" %in% (so %>% metadata %>% colnames)) ){
    if ( ("orig.ident" %in% (so %>% metadata %>% colnames)) ){
      RenameCells(so, new.names= so %>% metadata %>% dplyr::mutate(newcellid=paste0(pref,orig.ident,"_", so %>% colnames)) %>% pull(newcellid) )
    }else{
      RenameCells(so, new.names= so %>% metadata %>% dplyr::mutate(newcellid=paste0(pref, so %>% colnames)) %>% pull(newcellid) )
    }}else{
      RenameCells(so, new.names= so %>% metadata %>% dplyr::mutate(newcellid=paste0(pref, patientID, "_", so %>% colnames)) %>% pull(newcellid) )
    }
  
}

##assuming there is a patientID
prefixids= function(so, pref=""){
  RenameCells(so, new.names= so %>% metadata %>% dplyr::mutate(newcellid=paste0(pref, so %>% colnames)) %>% pull(newcellid) )
  
}


prefixids.matrix= function(mat, prefix){
  colnames(mat)= paste_(prefix,  colnames(mat))
  mat
}

prefixids.meta= function(mat, prefix){
  rownames(mat)= paste_(prefix,  rownames(mat))
  mat
}

alph=77
applyalpha=function(color) color %>% paste0(., alph) 


################################################################################
# SPECIFIC NCNB2 functions
################################################################################

# split the hto_demux name into different bits. 


extractname = Vectorize(function(x) {
  #cat(x, "\n")
  x=gsub("-", "_", x)
  if (!grepl("_", x)) {
    return(c(NA, NA, NA, NA) %>% givename(., c("cell.line", "condition", "stage", "replicate")))
  } else{
    #HTO-H7_17q_day3_rep1_CMO
    
    #small replacement on the condition section of the string
    #cat("1")
    if (!grepl("^H7", x)) {
      x0=paste0("H7_", x)
    }else{x0=x}
    
    x1 = gsub("H7WT", "H7_WT", x0)
    x2 = gsub("17q1qMYCN_", "17qx_", x1) #transform all the potential17q1qs into 17qx
    x3 = gsub("17qMYCN_", "17q1qMYCNdox_", x2) #transform all the potential mycns into mycndox
    x4 = gsub("MYCNdox_", "MYCN_", x3) #transform all dox into mycn
    x5 = gsub("17qx", "17q1q", x4) #tranform all qx into 17q1q
    
    
    x6 = gsub("_CMO", "", x5)
    x7= gsub("rep","R", x6)
    x8= gsub("day", "D", x7)
    
    xf=x8
    
    #make a list with all the parts thereafter
    xs1 = strsplit(xf, split = "_")[[1]]
    #print(x2)
    names(xs1) = c("cell.line", "condition", "stage", "replicate")
    xs2 = as.list(xs1)
    xs2$condition = paste0("c", xs2$condition)
    
    xs2 %>% unlist
  }
  
}, USE.NAMES=F)

################################################################################
# screen print shortcuts
################################################################################


#concatenate a message for the screen, add screen format

fcat= function(...) cat(paste(...,sep=" "), "\n")
msg=fcat
showdataset=function(so)
  fcat("Dataset", so$dsname %>% unique %>% paste(., collapse="-")) 


################################################################################
# dataset filtering on the go
################################################################################

filterds=function(so, dbcutoff=NULL, filter.empty.drops=NULL, demuxvar="sample_name", qc.pars=NULL){
  if(!is.null(qc.pars)){
    qcpars=qc.pars 
  }
  
  nm=Project(so)
  if(!(nm %in% names(qcpars))){
    st=paste(nm, "does not have an assigned QC filter")
    stop(st)
  }else{
    
    fcat("Updating filter metadata...")  
    so=AddMetaData(so, so %>% metadata %>% dplyr::mutate(counts.cutoff=qcpars[[nm]]$counts, features.cutoff= qcpars[[nm]]$features, mito.cutoff=qcpars[[nm]]$mito, qcc.pass=((nFeature_RNA>=features.cutoff & nCount_RNA >= counts.cutoff & percent.mt<= mito.cutoff )))  %>% dplyr::select(counts.cutoff, features.cutoff, mito.cutoff, qcc.pass))
  }
  if(demuxvar %in% (so %>% metadata %>% colnames)){
    so=so[, so %>% metadata %>% filter(!(!!sym(demuxvar) %in% c("Negative", "Doublet", "Multiplet", "Unassigned", "Blank", NA)), nCount_RNA>=qcpars[[nm]]$counts, nFeature_RNA>=qcpars[[nm]]$features, percent.mt<=qcpars[[nm]]$mito ) %>% rownames] 
    
  }else{
    so=so[, so %>% metadata %>% filter(nCount_RNA>=qcpars[[nm]]$counts, nFeature_RNA>=qcpars[[nm]]$features, percent.mt<=qcpars[[nm]]$mito ) %>% rownames] 
  }  
  
  if( !is.null(dbcutoff) & ("scdblscore" %in% (so %>% metadata %>% colnames))){
    so=so[, so %>% metadata %>% filter(scdblscore<=dbcutoff) %>% rownames]
    so=AddMetaData(so, metadata= so %>% metadata %>% dplyr::mutate(scDblFinder.cutoff=dbcutoff) %>% select(scDblFinder.cutoff))
  }
  if( !is.null(filter.empty.drops) & ("isemptydroplet" %in% (so %>% metadata %>% colnames))){
    so=so[, so %>% metadata %>% filter(isemptydroplet==FALSE) %>% rownames]
  }
  
  
  so   
}


filteremptydrops=function(so, filter.empty.drops=NULL){
  nm=Project(so)
  
  if( !is.null(filter.empty.drops) & ("isemptydroplet" %in% (so %>% metadata %>% colnames))){
    so=so[, so %>% metadata %>% filter(isemptydroplet==FALSE) %>% rownames]
  }
  
  
  so   
}



stagename= Vectorize(function(nm) stages[[nm]], USE.NAMES=F)


################################################################################
#reorder cells in clusters according to different methods (e.g., seriation or signature, i.e. how trong the markers are)
################################################################################

cellinfo.getcells=function(so, clusvar="seurat_clusters",ncells=100, clusters=NULL, fullrecreate=T, replicatevar=NULL, fullreload=F, ...){ #replicatevar=NULL,  clusters=NULL){
  
  Idents(so)=clusvar
  if(is.null(clusters)){
    clusids= so[[clusvar]] %>% pull(clusvar) %>% as.character %>% unique %>% sort
  }else{
    clusids=clusters 
  }
  
  fcat("Getting cells for each cluster...")
  cellist=lapply(as.character(clusids), function(x) {
    out=WhichCells(so, idents=x, downsample=ncells)
    #fcat("total number of cells:", length(out))
    out
  }
  
  ) %>% givename(., as.character(clusids))
  
  cellist
}


cellinfo.getmarkers=function(so, clusvar="seurat_clusters", selecttop=NULL){
  
  if(is.null(so$RNA@misc$top_markers)){  
    so=seuratmarkers.delegate(so, minfc=1, group_column=clusvar, replicate_column=replicatevar, fullrecreate=fullrecreate, fullreload=fullreload, ...)
  }
  
}





seriatecells=function(so, clusvar="seurat_clusters", clusters=NULL, meth="signature", markers=NULL, ncells=100, ass="SCT", groupvarname="group1",lfcvarname="log_fc", extended.output=T, deduped=T, nmarkers=NULL){
  
  Idents(so)=clusvar
  if(is.null(clusters)){
    clusids= so[[clusvar]] %>% pull(clusvar) %>% as.character %>% unique %>% sort
  }else{
    clusids=clusters 
  }
  
  
  markers= so@assays$RNA@misc$top_markers
  
  
  if(is.null(markers)){
    stop("There are no markers! please add markers under @assays$RNA@misc$top_markers\n" ) 
  }else{
    fcat("Splitting groups...")
    markerlist0=lapply(clusids, function(x) markers %>% filter(!!sym(groupvarname)==x)) %>% givename(., clusids)
    
    hasgenes= Vectorize(function(tabb) as.logical(nrow(tabb)>0), USE.NAMES=F)
    
    
    fcat("current length of markerlist0", length(markerlist0))
    which.have.genes=hasgenes(markerlist0)
    markerlist0=markerlist0[which.have.genes]
    fcat("current length of markerlist0 after filtering", length(markerlist0))
    fcat("these clusters have genes:", paste(markerlist0 %>% names, collapse=" ")) 
    
    fcat("current length of clusids", length(clusids))
    clusids=as.character(clusids)[which.have.genes]
    fcat("current length of clusids after filter", length(clusids))
    names(markerlist0)= as.character(clusids) 
    
    
    if(deduped==T){
      fullgns= lapply(markerlist0, function(x) x %>% pull(feature)) %>% Reduce(c, .)
      
      ###when duplicates, keep the gene as marker in the cluster where it shows highest log_fold change
      fcat("Assigning genes in more than one cluster to the highest expressing cluster...")    
      deduped=lapply(unique(fullgns), function(x){
        lapply(markerlist0, function(y) y %>% filter(feature==x)) %>% Reduce(condrbind, .) %>% arrange(-!!sym(lfcvarname)) %>% head(n=1)  %>% select(feature,!!sym(lfcvarname), !!sym(groupvarname)) 
        
      }) %>% Reduce(rbind, .)
      
      deduped2= deduped %>% dplyr::mutate(ncluster=sapply(deduped[[groupvarname]], function(x) generalreplace(x, clusids,1:length(clusids)), USE.NAMES=F )) %>% arrange(ncluster)
      
      #arrange genes from highest log fold change to lowest per cluster
      markers.revised= deduped2 %>% group_by(ncluster) %>% arrange(-!!sym(lfcvarname), .by_group=T) %>% ungroup()
      
      
    }else{
      markers.revised=markers # forwarding the original table
    }
    
    
    markergaps= markers.revised %>% pull(ncluster) %>% as.integer %>% diff %>% as.logical %>% which
    markervector=markers.revised %>% pull(feature) 
    
    if(is.null(nmarkers)){    
      markerlist=lapply(as.character(clusids), function(x) markers.revised %>% filter(!!sym(groupvarname)==x) %>% pull(feature))
    }else{
      markerlist=lapply(as.character(clusids), function(x) markers.revised %>% filter(!!sym(groupvarname)==x) %>% pull(feature) %>% head(nmarkers))
      
    }
    names(markerlist)=as.character(clusids)  
    markergaps=lapply(1:length(markerlist), function(x) rep(x, length(markerlist[[x]]))) %>% Reduce(c, .) %>% diff %>% as.logical %>% which
    
  }
  
  
  
  fcat("Getting cells for each cluster...")
  cellist=lapply(as.character(clusids), function(x) {
    out=WhichCells(so, idents=x, downsample=ncells)
    #fcat("total number of cells:", length(out))
    out
  }
  
  ) %>% givename(., as.character(clusids))
  
  
  ordcellist=tryCatch({
    
    if (meth=="rankcounts"){
      celltotals=list()
      ordcellist=lapply(1:length(cellist), function(x) {
        
        newcells= metadata(so)[cellist[[x]], ] %>% select(!!sym(clusvar), nCount_RNA) %>% arrange(-nCount_RNA) %>% slice_head(n=ncells.show) %>% rownames
        celltotals[[x]]<<- rep(x, length(newcells))
        newcells
      }
      )
      
      
    }else{
      if (meth=="signature"){
        fcat("Using signature method to rank cells...")
        celltotals=list()
        ordcellist=lapply(1:length(cellist), function(x) {
          
          print(markerlist[[x]])
          
          sumcells=apply(so[markerlist[[x]],cellist[[x]]]@assays[[ass]]@counts, 2, sum)
          
          
          #arranged from weakest to strongest  
          newcells= metadata(so)[cellist[[x]], ] %>% dplyr::mutate(markersum= sumcells, signature=markersum/nCount_RNA) %>% select(!!sym(clusvar), signature) %>% arrange(signature) %>% slice_head(n=ncells) %>% rownames
          celltotals[[x]]<<- rep(x, length(newcells))
          newcells
        }
        ) %>% givename(., clusids)
        
      }else{
        if (meth=="seriate"){
          library(seriation)
          celltotals=list()
          ordcellist=lapply(1:length(cellist), function(x) {
            
            
            ss=seriation::get_order(seriate(as.matrix(so[markerlist[[x]],cellist[[x]]]@assays$RNA@data), margin=2, method="BEA_TSP"))
            celltotals[[x]]<<- rep(x, length(ss))
            cellist[[x]][ss]
          }
          )  %>% givename(., clusids)
          
          
        }else{
          
          if(meth=="seurat"){
            fcat("Using Seurat Module Score method to rank cells...")
            celltotals=list()
            
            newnames=names(markerlist) %>% make.names
            so <- AddModuleScore(so, markerlist, name=newnames )
            
            
            ordcellist=lapply(1:length(clusids), function(x) {
              
              print(markerlist[[x]])
              
              modulename= paste0(newnames[x], x)
              
              #sumcells=apply(so[markerlist[[x]],so %>% metadata %>% pull(), 2, sum)
              
              newcells=so@meta.data[cellist[[x]], ] %>% filter(!!sym(clusvar)==clusids[x]) %>% select(modulename) %>% arrange(!!sym(modulename)) %>% rownames
              
              #arranged from weakest to strongest  
              #newcells= metadata(so)[cellist[[x]], ] %>% dplyr::mutate(markersum= sumcells, signature=markersum/nCount_RNA) %>% select(!!sym(clusvar), signature) %>% arrange(signature) %>% slice_head(n=ncells) %>% rownames
              celltotals[[x]]<<- rep(x, length(newcells))
              newcells
            }
            ) %>% givename(., clusids)
            
            
          }else{
            
            if(meth=="seurat2"){
              fcat("Using Seurat Module Score method to rank cells and retrieveng strongest", ncells, "cells in clusters...")
              celltotals=list()
              
              so <- AddModuleScore(so, markerlist, name=names(markerlist) )
              
              
              ordcellist=lapply(1:length(clusids), function(x) {
                
                print(markerlist[[x]])
                
                modulename= paste0(clusids[x], x)
                
                #sumcells=apply(so[markerlist[[x]],so %>% metadata %>% pull(), 2, sum)
                newcells=so@meta.data %>% filter(!!sym(clusvar)==clusids[x]) %>% select(modulename) %>% arrange(!!sym(modulename)) %>% rownames
                newcells=newcells[1:ncells]
                
                #arranged from weakest to strongest  
                #newcells= metadata(so)[cellist[[x]], ] %>% dplyr::mutate(markersum= sumcells, signature=markersum/nCount_RNA) %>% select(!!sym(clusvar), signature) %>% arrange(signature) %>% slice_head(n=ncells) %>% rownames
                celltotals[[x]]<<- rep(x, length(newcells))
                newcells
              }
              )
              
              
              
            }
            
          }
          
          
        }}}
    ordcellist} , error=function(e){fcat("There was trouble with seriating the cells. returning unseriated cells");
      print(e); 
      cellist
    })
  
  celltotals=lapply(1:length(ordcellist), function(x) rep(x, length(ordcellist[[x]])))
  
  fcat("making cell vector...")
  cellvector=Reduce(c, ordcellist)
  fcat("getting the gap positions...")
  cellgaps= celltotals %>% Reduce(c,.) %>% diff %>% as.logical %>% which  
  
  cell.reference=lapply(1:length(ordcellist), function(nn) rep(names(ordcellist)[nn], length(ordcellist[[nn]])   )) %>% Reduce(c, .)  
  marker.reference=lapply(1:length(markerlist), function(nn) rep(names(markerlist)[nn], length(markerlist[[nn]])   )) %>% Reduce(c, .)
  
  cell.annotation= as.data.frame(cell.reference) %>% giverownames(cellvector) %>% givecolnames(., 1, clusvar)
  marker.annotation= as.data.frame(marker.reference) %>% giverownames(markervector) %>% givecolnames(., 1, clusvar)
  
  
  
  
  if(extended.output==T){
    list(cells=cellvector,
         cell.list=ordcellist,
         gaps.cells=cellgaps,
         gaps.markers=markergaps,
         markers=markervector,
         markerlist=markerlist,
         cell.annotation=cell.annotation,
         marker.annotation=marker.annotation,
         clustering.variable=clusvar,
         input.cluster.ids=clusids,
         method=meth,
         cell.metadata=clusvar,
         marker.metadata=clusvar
    )}else{
      cellvector
    }
}



seriatecells2=function(so, clusvar="seurat_clusters", clusters=NULL, meth="signature", markers=NULL, ncells=100, ass="SCT", groupvarname="group",lfcvarname="logFC", extended.output=F, deduped=T, groups.match.markers=T){
  
  Idents(so)=clusvar
  if(is.null(clusters)){
    clusids= so[[clusvar]] %>% pull(clusvar) %>% as.character %>% unique %>% sort
  }else{
    clusids=clusters 
  }
  markers.are.list=F
  markers.are.table=F
  ## check if markers have been provided externally
  if(!is.null(markers) & class(markers)=="list"){
    markers.are.list=T
    markers.are.table=F
    fcat("External marker list has been provided.")
  }else{
    markers.are.list=F
    if(!is.null(markers) & (class(markers) %in% c("data.table", "matrix","data.frame"))){
      fcat("External marker table has been provided.")
      markers.are.table=T
    }
  }
  ### CASE WHEN EXTERNAL LIST OF MARKERS IS PROVIDED
  if(markers.are.list){    
    markerlist0=markers
    
    if(!is.null(names(markerlist0))){# if markerlist has names
      clusids=names(markerlist0)
      
      #all(clusters %in% (ref.cellinfo$markerlist %>% names))
      
    }else{ ### if markerlist has no names just make some
      fcat("external markers have no names. Automatically generating names...")
      clusids=(1:length(markerlist0)) %>% paste0("C", .)
      names(markerlist0)= clusids
      
    }}
  
  # assemble the markers into a simple dataframe to perform duplication operations
  if(deduped==T){
    fcat("deduplicating marker list...")
    fullgns=lapply(1:length(markerlist0), function(x) markerlist0[[x]]) %>% Reduce(c, .)
    fullgns.groups=lapply(1:length(markerlist0), function(x) rep(names(markerlist0)[x], length(markerlist0[[x]])) ) %>% Reduce(c, .)
    
    # in the absence of GE information, remove duplicates in the order they appear first time marker appears is considered the 
    #most relevant. 
    dupd= duplicated(fullgns)
    markers.df=as.data.frame(fullgns.groups[!dupd]) %>% giverownames(fullgns[!dupd]) %>% givecolnames(., 1, nms=clusvarname) 
    
    markerlist= lapply(as.character(clusids), function(x)  markers.df %>% filter(!!sym(clusvarname)==!!x) %>% rownames)
    names(markerlist)=as.character(clusids)
  }
  
  ### NO EXTERNAL MARKERS PROVIDED: MARKER TABLE EXPECTED WITHIN OBJECT
  if(!markers.are.list & !markers.are.table & is.null(markers)){
    #if external markers are not provided at all, start from the top markers table
    markers= so@assays$RNA@misc$top_markers
    if(is.null(markers)){
      markers= so@assays$RNA@misc$markers
      
      if(is.null(markers)){
        #if top_markers is not there, then throw an error
        stop("There are no markers! please add markers under @assays$RNA@misc$top_markers\n" ) 
      }else{
        markers.are.table=T 
      }
      
    }else{
      markers.are.table=T
    }
    
    if(markers.are.table){
      # if there is a table, collect the DE markers for each group. 
      fcat("Splitting groups and retrieving top markers from table...")
      markerlist0=lapply(clusids, function(x) markers %>% filter(!!sym(groupvarname)==x))
      
      hasgenes= Vectorize(function(tabb) as.logical(nrow(tabb)>0), USE.NAMES=F)
      fcat("current length of markerlist0", length(markerlist0))
      which.have.genes=hasgenes(markerlist0)
      markerlist0=markerlist0[which.have.genes]
      fcat("current length of markerlist0 after filtering", length(markerlist0))
      fcat("current length of clusids", length(clusids))
      clusids=as.character(clusids)[which.have.genes]
      fcat("current length of clusids after filter", length(clusids))
      names(markerlist0)= clusids 
      
      
      fcat("Markers found for the following groups:", paste(clusids, collapse=","))
      
    } 
    
    #removing duplicates when a marker DE table is provided (not a list)
    if(deduped==T & markers.are.table){
      fullgns= lapply(markerlist0, function(x) x %>% pull(feature)) %>% Reduce(c, .)
      
      ###when duplicates, keep the gene as marker in the cluster where it shows highest log_fold change
      fcat("Assigning genes in more than one cluster to the highest expressing cluster...")    
      deduped=lapply(unique(fullgns), function(x){
        lapply(markerlist0, function(y) y %>% filter(feature==x)) %>% Reduce(condrbind, .) %>% arrange(-!!sym(lfcvarname)) %>% head(n=1)  %>% select(feature,!!sym(lfcvarname), !!sym(groupvarname)) 
        
      }) %>% Reduce(rbind, .)
      
      deduped2= deduped %>% dplyr::mutate(ncluster=sapply(deduped[[groupvarname]], function(x) generalreplace(x, clusids,1:length(clusids)), USE.NAMES=F )) %>% arrange(ncluster)
      
      #arrange genes from highest log fold change to lowest per cluster
      markers.revised= deduped2 %>% group_by(ncluster) %>% arrange(-!!sym(lfcvarname), .by_group=T) %>% ungroup()
      
      
    }
    
    
    markergaps= markers.revised %>% pull(ncluster) %>% as.integer %>% diff %>% as.logical %>% which
    #markervector=markers.revised %>% pull(feature) 
    
    markerlist=lapply(as.character(clusids), function(x) markers.revised %>% filter(!!sym(groupvarname)==x) %>% pull(feature))
    
    names(markerlist)=as.character(clusids)  
    
  }
  
  ### COMMON DOWNSTREAM PROCESSING AFTER DEFINING THE MARKER LIST 
  
  markervector= markerlist %>% Reduce(c, .)
  markergaps=lapply(1:length(markerlist), function(x) rep(x, length(markerlist[[x]]))) %>% Reduce(c, .) %>% diff %>% as.logical %>% which
  
  
  fcat("Getting cells for each cluster...")
  cellclasses=so %>% metadata %>% pull(clusvar) %>% unique
  fcat("cellclasses are", cellclasses)
  cellist=lapply(cellclasses, function(x) {
    
    #  if(!(x %in% ( cellclasses ))  ){
    #  NULL}else{
    out=WhichCells(so, idents=x, downsample=ncells)
    #fcat("total number of cells:", length(out))
    out
    #}
  }
  
  ) %>% givename(., cellclasses)
  
  #fcat("removing clusters for which no cells were found")
  #nonnulls=lapply(cellist, function(x) is.null(x)) %>% unlist %>% unname
  #cellist= cellist[!nonnulls] %>% givename(., as.character(clusids[!nonnulls]))
  
  
  if (meth=="rankcounts"){
    celltotals=list()
    ordcellist=lapply(1:length(cellist), function(x) {
      
      newcells= metadata(so)[cellist[[x]], ] %>% select(!!sym(clusvar), nCount_RNA) %>% arrange(-nCount_RNA) %>% slice_head(n=ncells.show) %>% rownames
      celltotals[[x]]<<- rep(x, length(newcells))
      newcells
    }
    )
    
    
  }else{
    if (meth=="signature"){
      fcat("Using signature method to rank cells...")
      celltotals=list()
      ordcellist=lapply(1:length(cellist), function(x) {
        
        print(markerlist[[x]])
        
        sumcells=apply(so[markerlist[[x]],cellist[[x]]]@assays[[ass]]@counts, 2, sum)
        
        
        #arranged from weakest to strongest  
        newcells= metadata(so)[cellist[[x]], ] %>% dplyr::mutate(markersum= sumcells, signature=markersum/nCount_RNA) %>% select(!!sym(clusvar), signature) %>% arrange(signature) %>% slice_head(n=ncells) %>% rownames
        celltotals[[x]]<<- rep(x, length(newcells))
        newcells
      }
      )
      
      cellvector=Reduce(c, ordcellist)
      fcat("getting the gap positions...")
      cellgaps= celltotals %>% Reduce(c,.) %>% diff %>% as.logical %>% which
    }else{
      if (meth=="seriate"){
        library(seriation)
        celltotals=list()
        ordcellist=lapply(1:length(cellist), function(x) {
          
          
          ss=seriation::get_order(seriate(as.matrix(so[markerlist[[x]],cellist[[x]]]@assays$RNA@data), margin=2, method="BEA_TSP"))
          celltotals[[x]]<<- rep(x, length(ss))
          cellist[[x]][ss]
        }
        ) 
        
        
      }else{
        
        if(meth=="seurat"){
          fcat("Using Seurat Module Score method to rank cells...")
          celltotals=list()
          
          newnames=names(markerlist) %>% make.names
          #names(newnames)=names(markerlist)
          fcat("newnames is ", newnames)
          so=tryCatch({ AddModuleScore(so, markerlist, name=newnames)}, error=function(e){
            fcat("trying 12 bins because of numerical issues...")
            AddModuleScore(so, markerlist, nbin=12, name=newnames)})
          
          
          if(groups.match.markers){
            ordcellist=lapply(1:length(markerlist), function(x) {
              
              print(names(markerlist)[x])
              
              modulename= paste0(newnames[x], x)
              markername=names(markerlist)[x]
              #sumcells=apply(so[markerlist[[x]],so %>% metadata %>% pull(), 2, sum)
              newcells=so@meta.data[cellist[[markername]], ] %>% select(modulename) %>% arrange(!!sym(modulename)) %>% rownames
              if(!is.null(ncells) & ncells<= length(newcells)){
                #get the ncells strongest cells
                newcells=newcells[(length(newcells)-ncells):length(newcells)]
              }
              #arranged from weakest to strongest  
              #newcells= metadata(so)[cellist[[x]], ] %>% dplyr::mutate(markersum= sumcells, signature=markersum/nCount_RNA) %>% select(!!sym(clusvar), signature) %>% arrange(signature) %>% slice_head(n=ncells) %>% rownames
              celltotals[[x]]<<- rep(x, length(newcells))
              newcells
            }
            )
            
            cellvector=Reduce(c, ordcellist)
            fcat("getting the gap positions...")
            cellgaps= celltotals %>% Reduce(c,.) %>% diff %>% as.logical %>% which
          }else{ # if cell groups do not match markers, we just retrieve the ncells with the strongest signature
            fcat("cells and markers do not have corresponding ids. Getting cells with strongest signature")
            ordcellist=lapply(1:length(clusids), function(x) { #we evaluate the strength of all signatures regardless of the cell label
              
              print(names(markerlist)[x])
              
              modulename= paste0(newnames[clusids[x]], x)
              
              #arrange from biggest to smallest
              newcells=so@meta.data %>% select(modulename) %>% arrange(!!sym(modulename)) %>% rownames
              
              if(!is.null(ncells) & ncells<= length(newcells)){
                #get the ncells strongest cells
                newcells=newcells[(length(newcells)-ncells):length(newcells)]
              }
              #arranged from weakest to strongest  
              #newcells= metadata(so)[cellist[[x]], ] %>% dplyr::mutate(markersum= sumcells, signature=markersum/nCount_RNA) %>% select(!!sym(clusvar), signature) %>% arrange(signature) %>% slice_head(n=ncells) %>% rownames
              celltotals[[x]]<<- rep(x, length(newcells))
              newcells
            }
            )
            
            cellvector=Reduce(c, ordcellist)
            fcat("getting the gap positions...")
            cellgaps= celltotals %>% Reduce(c,.) %>% diff %>% as.logical %>% which
            
          } 
        }else{ #another method different from above
          
          if(meth=="seurat2"){
            fcat("Using Seurat Module Score method to rank cells and retrieveng strongest", ncells, "cells in clusters...")
            celltotals=list()
            
            so <- AddModuleScore(so, markerlist, name=names(markerlist) )
            
            
            ordcellist=lapply(1:length(clusids), function(x) {
              
              print(markerlist[[x]])
              
              modulename= paste0(clusids[x], x)
              
              #sumcells=apply(so[markerlist[[x]],so %>% metadata %>% pull(), 2, sum)
              newcells=so@meta.data %>% filter(!!sym(clusvar)==clusids[x]) %>% select(modulename) %>% arrange(!!sym(modulename)) %>% rownames
              newcells=newcells[1:ncells]
              
              #arranged from weakest to strongest  
              #newcells= metadata(so)[cellist[[x]], ] %>% dplyr::mutate(markersum= sumcells, signature=markersum/nCount_RNA) %>% select(!!sym(clusvar), signature) %>% arrange(signature) %>% slice_head(n=ncells) %>% rownames
              celltotals[[x]]<<- rep(x, length(newcells))
              newcells
            }
            )
            
            cellvector=Reduce(c, ordcellist)
            fcat("getting the gap positions...")
            cellgaps= celltotals %>% Reduce(c,.) %>% diff %>% as.logical %>% which
            
            
          }
          
        }
        
        
      }
      
    }
  }
  if(extended.output==T){
    list(cells=cellvector,cell.list=ordcellist , gaps.cells=cellgaps, gaps.markers=markergaps, markers=markervector, markerlist=markerlist)}else{
      cellvector
    }
}

################################################################################
# classify cells into singlet, doublet or negative based on the hto_demux
################################################################################

classifydemux= Vectorize(function(x){
  if(x!="Negative" & x!= "Doublet" & x!= "Unassigned"& x!= "Blank" & x!= "Multiplet"){
    return("Singlet")
  }else{
    return(x)
  }
}, USE.NAMES=F)


################################################################################
# retrieve a doublet scoring filter for a dataset
################################################################################

getdbcutoff= Vectorize(function(nm){
  dbcutoffs[[nm]]
  
}, USE.NAMES=F)


################################################################################
# paste two strings together with a comma in between (great with Reduce())
################################################################################


pastec=function(a,b) paste0(c(a,b), collapse=",")


################################################################################
# a more comprehensive condition taking into account the demux state, negative or doublet. 
################################################################################

get_condition_demux=Vectorize(function(condition, demux){
  
  if(is.na(condition)){
    
    if(!(is.na(demux))){
      return(demux)
    }else{
      return("Blank")
      
    }
    
  }else{
    return(condition)
  }
  
}, USE.NAMES=F)


################################################################################
# remove the last character from a string
################################################################################
chopstring=function(st) substr(st, 1, nchar(st)-1)


################################################################################
#fix bad names in 
################################################################################

fixnames= function(x){nms= x %>% colnames %>% make.names; colnames(x)=nms; fcat("new names are", nms); x}


################################################################################
#get the replicate
################################################################################
getrep= Vectorize(function(x) repmap[[x]], USE.NAMES=F)


################################################################################
# trim final one in the cell names
###############################################################################

trimone= Vectorize(function(x) strsplit(x, split="-")[[1]][1]) 


################################################################################
#
################################################################################

is.loaded= function(pkgname) !requireNamespace(pkgname, quietly = TRUE)

################################################################################
# make square legend
################################################################################

squarelegend = function(varr, size=1, shape=0, colorlist=allcolors){
  
  arglist=list()
  arglist[[varr]]= guide_legend(
    override.aes = list(
      shape = shape,  # 0: square, 1: circle
      size = size,   # Customize point sizes
      color = allcolors[[varr]],  # Customize point colors
      label = allcolors[[varr]] %>% names # Customize legend labels
    )
  )
  
  do.call(guides, arglist)
  
}

################################################################################
# compare category changes for cells across several category variables
################################################################################

transferplot= function(so=NULL, df=NULL,  labelvars, colorlist=allcolors, facetvar="all_cells", facet.ncol=NULL, facet.nrow=NULL, column.width=0.5, alpha=0.5, scale=5){
  
  
  sca=scale  
  freqvar=facetvar
  library(ggalluvial)
  library(ggrepel)
  
  #prepare colors if unavailable

  
  #lapply(labelvars, function(x) allcolors[[k]])    
  if(!is.null(df)){
    met=df %>% dplyr::mutate(all_cells="All cells")
  }else{
    if(!is.null(so)){  
      met= so %>% metadata %>% dplyr::mutate(all_cells="All cells")  
    }else{
      error("Please add either a Seurat object or a data frame as input!")
    }}
  
  
  fcat("preparing colors")
  for(k in 1:length(labelvars)){
    
    if(is.null(colorlist[[labelvars[k]]])){
      fcat(as.character(labelvars[k]))
      cats=met %>% pull(!!sym(as.character(labelvars[k]))) %>% unique
      fcat("unique categories:", paste(cats, collapse=" "))
      colorlist[[as.character(labelvars[k])]]<- randomcolors( cats %>% length) %>% givename(., cats)
      #allcolors[[labelvars[k]]]<<- colorlist[[labelvars[k]]]
      
    }
  }
  fcat("Colors used for plotting:")
  cat(code.for.list.of.vectors(colorlist[labelvars]))
  
  
   fcat("creating factors")
  
  # factors suspended for now
  for(k in 1:length(labelvars)){
    
    met = met %>% dplyr::mutate(!!labelvars[k]:=factor(!!sym(labelvars[k]), levels=colorlist[[labelvars[k]]] %>% names ) )
    
  }
 
    fcat("calculating frequencies t1")

  t1=met %>% dplyr::select(c(freqvar,labelvars)) %>% group_by(!!!syms(c(freqvar, labelvars))) %>% summarise( counts=n())  %>% group_by(!!sym(freqvar)) %>% dplyr::mutate(Freq=counts/sum(counts)) %>% ungroup() 
   fcat("calculating frequencies t2")
  t2= t1 %>% dplyr::mutate(ind= 1:nrow(t1)) %>% pivot_longer(., cols=labelvars, names_to="originalcolumn", values_to="cell.label")
  
   fcat("making plot")
  tpdf(path=config$plotpath, paste0("alluvialplot_labeltransfer", paste(labelvars, collapse="-")), width=pw*5/2*sca/10, height=pw*sca*3/14)
  ga=ggplot(t2,
            aes(y = Freq, x=factor(originalcolumn, levels=labelvars), stratum=cell.label, alluvium=ind, fill=cell.label))+ 
    geom_flow() +
    geom_stratum(width = column.width, alpha=alpha) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)))+#, min.segment.length=0.1) +
    scale_x_discrete(expand= c(.1, .1), labels = labelvars)+
    facet_wrap(~factor(get(freqvar), levels=(colorlist[[freqvar]] %>% names)), ncol=facet.ncol, nrow=facet.nrow)+
    #scale_fill_manual(values=colorlist[labelvars] %>% Reduce(c,.))+
    theme_classic()
  print(ga+NoLegend())
  dev.off()
  
  ga+NoLegend()
}

################################################################################
# get highest loadings from pca
################################################################################

gethighloadings=function(so, comps=1:5, nfeatures=10, reduction="pca", reduction.key="PC_"){
  rk=reduction.key
  compnms= paste0(rk, comps)
  flist= lapply(compnms, function(x){
    
    list(positive=
           so@reductions[[reduction]]@feature.loadings %>% as.data.frame %>% select(!!sym(x)) %>% arrange(-!!sym(x)) %>% head(nfeatures) %>% rownames
         , negative=
           so@reductions[[reduction]]@feature.loadings %>% as.data.frame %>% select(!!sym(x)) %>% arrange(!!sym(x)) %>% head(nfeatures) %>% rownames
    ) 
    
  }) %>% givename(compnms)
  flist
}


################################################################################

################################################################################
code.for.list.of.vectors <- function(list_of_vectors) {
  code0 <- "list(\n"
  codelist=list()
  nms=names(list_of_vectors)
  # Loop through each vector in the list
  for (i in seq_along(list_of_vectors)) {
    
    if(is.null(nms)){
      vector_name <- paste0("vector", i)
    }else{
      vector_name <- nms[i] 
    }
    code <-  paste0("  ", vector_name, " = c(\n")
    
    # Loop through each element in the vector
    for (element_name in names(list_of_vectors[[i]])) {
      element_value <- toString(list_of_vectors[[i]][element_name])
      code <- paste0(code, paste0("   ", element_name, " = '", element_value, "',\n"))
    }
    
    code <- paste0(code, ")\n")
    codelist[[i]]=code
  }
  code2= codelist %>% Reduce(c, .) %>% paste(., collapse=",")
  
  code3 <- paste0(code0, code2, ")\n") 
  
  return(code3 %>% gsub(",\n)", ")", .))
}

# Example usage:
#vector1 <- c("name1" = "value1", "name2" = "value2")
#vector2 <- c("name3" = "value3", "name4" = "value4")
#list_of_vectors <- list(vector1, vector2)

#code_string <- code.for.list.of.vectors(allcolors[labelvars])
#code_string
#cat(code_string)

################################################################################
# remove a field from a list
################################################################################

dropfield= function(lisst, x){ lisst[[x]]=NULL; lisst} 


################################################################################
# create an arbitrary colormap
################################################################################

breakColors = function(breaks, colors, center=0, tol=0.001)
{
  ## In case of explicit color definitions
  nbreaks = length(breaks)
  nclass  = nbreaks - 1
  if (!is.function(colors)) {
    ncolors = length(colors)
    if (ncolors > nclass) {
      warning("more colors than classes: ignoring ", ncolors-nclass, " last colors")
      colors = colors[1:nclass]
    } else if (nclass > ncolors) {
      stop(nclass-ncolors, " more classes than colors defined")
    }
  } else {
    ## Are the classes symmetric and of same lengths?
    clens = diff(breaks)
    aclen = mean(clens)
    if (aclen==0) stop("Dude, your breaks are seriously fucked up!")
    relerr = max((clens-aclen)/aclen)
    if ( (center %in% breaks) & (relerr < tol) ) { ## yes, symmetric
      ndxcen = which(breaks==center)
      kneg = ndxcen -1 
      kpos = nbreaks - ndxcen
      kmax = max(kneg, kpos)
      colors = colors(2*kmax)
      if (kneg < kpos) {
        colors = colors[ (kpos-kneg+1) : (2*kmax) ]
      } else if (kneg > kpos) {
        colors = colors[ 1 : (2*kmax - (kneg-kpos)) ]
      }
    } else {                                      ## no, not symmetric
      colors = colors(nclass)
    }
  }
  colors
}

################################################################################

removegene=function(genes,x) genes[genes!=x]

################################################################################
# format recuciton key based on an arbitrary id.
################################################################################


prepare.rk=function(x)  paste0(x, "_") %>% gsub("\\.", "", .)


################################################################################
#
################################################################################


sling_cell_df <- function(sling_ctr){
  raw_pseudotime <- slingPseudotime(sling_ctr)
  order_df <- apply(raw_pseudotime, 2, function(x){
    x[!is.na(x)] <- rank(x[!is.na(x)])
    x
  }) %>%
    as.data.frame()%>% 
    rename_with(.,~paste0("order","_",.)) %>% 
    dplyr::mutate(mean_order = apply(., 1, mean, na.rm = TRUE))
  
  
  pseudotime_df <- apply(raw_pseudotime, 2, function(x){
    x1 <- x[!is.na(x)]
    x[!is.na(x)] <- (x1-min(x1))/(max(x1)-min(x1))*100
    x
  })%>%
    as.data.frame()%>% 
    rename_with(.,~paste0("pseudotime","_",.))%>% 
    dplyr::mutate(mean_pseudotime = apply(., 1, mean, na.rm = TRUE))
  
  df <- bind_cols(raw_pseudotime,pseudotime_df, order_df) %>% 
    as.data.frame()
  
  return(df)
}


################################################################################
# ENTROPY FUNCTION
################################################################################

#https://rpubs.com/philjet/shannonentropy
#entropy function
entropy <- function(target) {
  freq <- table(target)/length(target)
  # vectorize
  vec <- as.data.frame(freq)[,2]
  #drop 0 to avoid NaN resulting from log2
  vec<-vec[vec>0]
  #compute entropy
  -sum(vec * log2(vec))
}



################################################################################
# findmin
################################################################################

findmin = function(arr){
  
  m= min(abs(arr))
  return(which(arr==m))
  
}


################################################################################
# get stage number from a "DX" format
################################################################################

get.stage.number= Vectorize(
  function(x) substr(x, 2, 10) %>% as.numeric
)

################################################################################
#
################################################################################

prepare.rk=function(x)  paste0(x, "_") %>% gsub("\\.", "", .)



################################################################################
# dplyr - ready get rows
################################################################################


getrows=function(x, rows) x %>% names2col(., "nm") %>% filter(nm %in% rows) %>% select(-nm)


################################################################################
# generalreplace, grep version
################################################################################

#value may be longer than elements in reff and so may contain reff
generalreplace.grepl=function(vec, refs, targets, return.sums=T){
  subss=list()
  finalvector=vec
  for(n in 1:length(refs)){
    reff=refs[n]
    which.match=grepl(reff, vec)
    subss[[n]]=rep(FALSE, length(vec))
    subss[[n]][which.match]=TRUE
    finalvector=gsub(reff, targets[n], finalvector )
    
  }
  allsums=subss %>% Reduce('+', .)  
  if(any(allsums)>1){
    warning("Same value substituted multiple times")
  }
  if(return.sums==T){
    list(final.vector=finalvector, number.replacements=allsums)
  }else{
    finalvector 
  }
  
}
################################################################################
#
################################################################################

reorderdf= function(df, varr, catorder, return.list=T){
  
  l=lapply(catorder, function(x){
    df %>% filter(!!sym(varr)==x)
  })
  if(return.list){
  }else{
    l=l %>% Reduce(rbind,.)
  }
  l
}





################################################################################
# wrapper to take 2 seurat objects with a mapfun variable and show the points of one in another's clusters (groups)
################################################################################

#project.to.glasswork=function(ref, query, )
################################################################################
# functions that emulate matlab zeros and ones
################################################################################
zeros <- function(rows, cols) {
  return(matrix(0, nrow = rows, ncol = cols))
}

ones <- function(rows, cols) {
  return(matrix(1, nrow = rows, ncol = cols))
}


remove.zerorows= function(matt){

  if(length(dim(matt))<=1){
    print("Argument not a matrix. returning as is")
    matt
    }else{
    matt[apply(matt, 1, sum)!=0,]}

}
  
  remove.zerocols=function(matt){
      if(length(dim(matt))<=1){
    print("Argument not a matrix. returning as is")
    matt
        }else{
    matt[,apply(matt, 2, sum)!=0]}}

get.variant.rows= function(matt){
  
  #non.zero.rows=apply(matt, 1, sum)!=0
  variant.rows=apply(matt, 1, sd)!=0
  #print(cbind(non.zero.rows, variant.rows))
  variant.rows
  
}

get.invariant.rows= function(matt){
  
  #non.zero.rows=apply(matt, 1, sum)!=0
  invariant.rows=apply(matt, 1, sd)==0
  #print(cbind(non.zero.rows, variant.rows))
  invariant.rows
  
}

get.nonzerorows= function(matt) (apply(matt, 1, sum)!=0) %>% unname  


fill.missing.markers.matrix=function(matt, markernames, fill.value=0){
  missing.ones=NULL 
  
  newlines=lapply(1:length(markernames), function(k){ 
    
    nl=NULL
    if(!(markernames[k] %in% rownames(matt))){
      
      
      
      #missing.ones=c(missing.ones, markernames[k])
      
      nc=length(colnames(matt))
      nl=matrix(data=0, nrow=1, ncol= nc, dimnames=list(markernames[k], colnames(matt)))
    }
    nl
  }) %>% Reduce(rbind, .)
  #fcat("number of columns of matt", length(colnames(matt)), "number of columns of newlines", ncol(newlines))
  matt=rbind(matt, newlines)
  #if(!is.null(missing.ones)){
  
  #fcat("missing genes filled:", paste(missing.ones, collapse=",")) 
  #}
  matt
}



################################################################################
#for plotting purposes, slightly jitter one position of a row
################################################################################


random.glitch.row=function(mat, rownumber=x, val=0.00001){
  
  s=sample(1: ncol(mat), 1)
  mat[  rownumber,  s]=mat[  rownumber,  s]+val
  
  mat
  
}

################################################################################
# randomly glitch all invariant rows for plotting purposes
################################################################################

glitch.invariants=function(matt, val=0.00001){
  
  invariant.rows= !get.variant.rows(matt)   
  
  
  for(k in 1:length(invariant.rows)){
    set.seed(invariant.rows[k])
    matt=random.glitch.row(matt, rownumber=invariant.rows[k], val=val)
    
    
  }
  set.seed(42)
  matt
}

################################################################################
#create a heatmap based directly from a cellinfo object
################################################################################


cellinfo.heatmap= function(so,
                           cellinfo,
                           assay="RNA",
                           name=NULL,
                           pheatmap.params=NULL, ## currently not used
                           genes.to.label=NULL, 
                           color.labels.by=NULL,
                           colorlist=NULL,
                           cmap=NULL,
                           colorbar.max=2,
                           colorbar.spacing=0.1,
                           show.legend=TRUE,
                           return.everything=F,
                           return.cellinfo=F){
  
  #function refresh cell info to update ubject
  
  
  
  # heatmap colors
  
  if(is.null(cmap)){
    cmap=c("#0B9988", "#f7f7f7", "#F9AD03") 
  }
  
  colorpalette=colorRampPalette(cmap, space="rgb")
  
  ####Heatmap color schemes
  #petrol gold
  
  br=seq(-colorbar.max,colorbar.max,colorbar.spacing)  ## numeric breaks of the color bins
  colls=colorpalette(length(br))   
  
  #cellinfo= refresh.cellinfo(cellinfo)
  
  
  
  
  
  fcat("getting data matrix...")
  mat=so[cellinfo$markers, cellinfo$cells ]@assays[[assay]]@data %>% as.matrix
  fcat("finding invariant and absent genes...")
  absent.genes= setdiff(cellinfo$markers, rownames(mat))
  variant.rows=get.variant.rows(mat)
  invariant.genes=mat[!variant.rows, ] %>% rownames 
  fcat("Warning: the following invariant genes have been removed:\n", invariant.genes %>% dput)
  cellinfo2=removemarkers(cellinfo, markers=union(invariant.genes, absent.genes))
  
  if(!is.null(cellinfo$cell.metadata)){
    fcat("adding cell annotation from variables in metadata...")
    cellinfo2=cellinfo2 %>% add.cell.annotation(so, ., cellinfo2$cell.metadata)
  }
  
  
  
  if(!is.null(genes.to.label)){
    fcat("finding genes to label...")
    # original markers
    gns=cellinfo2$markers
    
    
    geneguide=  1:length(gns)
    names(geneguide)=gns
    #mrkrs are genes of interest
    
    
    
    if(!is.null(color.labels.by) & !is.null(colorlist)){
      
      gene.order.info=geneguide[genes.to.label]
      
      
      getcol=Vectorize(function(x) 
        if(x %in% (allcolors[[color.labels.by]] %>% names)){
          allcolors[[color.labels.by]][x]}else{
            NA
          }, USE.NAMES=F)
      
      getposition=Vectorize(function(x) geneguide[x], USE.NAMES=F)
      
      colref= (cellinfo2$marker.annotation %>% names2col(., "gene") %>% filter(gene %in% genes.to.label) %>% dplyr::mutate(color=getcol(!!sym(color.labels.by)), position= getposition(gene) )  )
      genecols=colref %>% arrange(position) %>% pull(color) %>% unname
      
      
      linesgp= grid::gpar(col=genecols)
      labelsgp= grid::gpar(col=genecols)
    }else{
      linesgp=NULL
      labelsgp=NULL
    }
    ha = ComplexHeatmap::rowAnnotation(foo = ComplexHeatmap::anno_mark(at = geneguide[genes.to.label] %>% unname %>% as.numeric, 
                                                                       labels = genes.to.label,
                                                                       lines_gp= linesgp, 
                                                                       labels_gp= labelsgp))
    
    
    #geneguide[genes.to.label] %>% as.character %>% dput
    
  }else{ha=NULL}
  
  fcat("getting data matrix...")
  
  mat2=so[cellinfo2$markers, cellinfo2$cells ]@assays[[assay]]@data %>% as.matrix
  
  cats=cellinfo2$cell.metadata
  
  for(ct in cats){
    
    if(!is.null(allcolors[[ct]])){  
      cellinfo2$cell.annotation[[ct]]= factor(cellinfo2$cell.annotation[[ct]], levels=allcolors[[ct]] %>% names) 
    }
  }
  
  
  ph= ComplexHeatmap::pheatmap(mat2,
                               scale="row",
                               cluster_row=FALSE,
                               show_colnames=FALSE,
                               show_rownames=FALSE,
                               cluster_col=FALSE, 
                               breaks=br,
                               col=colls, 
                               annotation_col=cellinfo2$cell.annotation,
                               annotation_row=cellinfo2$marker.annotation,
                               gaps_col= cellinfo2$gaps.cells,
                               gaps_row= cellinfo2$gaps.markers,
                               annotation_colors=colorlist,
                               border_color=NA,
                               right_annotation=ha, 
                               fontsize=5,
                               main=name 
                               #show_heatmap_legend=show.legend
                               #column_gap = unit(.2, "mm")
  )
  
  
  if(return.everything){
    
    list(heatmap=ph, cellinfo=cellinfo2, matrix=mat2)
  }else{
    
    if(return.cellinfo){
      return(cellinfo2) 
    }else{
      
      return(ph)
    }
    
  }
  
}



################################################################################
#refresh cellinfo object
################################################################################


refresh.cellinfo=function(cellinfo, cell.group.label="cell.group", marker.group.label="marker.group"){
  fcat("readjusting")
  
  cellinfo$markers= cellinfo$markerlist %>% Reduce(c, .)
  cellinfo$cells= cellinfo$cell.list %>% Reduce(c, .)
  fcat("making references")
  cell.reference=lapply(1:length(cellinfo$cell.list), function(nn) rep(names(cellinfo$cell.list)[nn], length(cellinfo$cell.list[[nn]])   )) %>% Reduce(c, .)  
  marker.reference=lapply(1:length(cellinfo$markerlist), function(nn) rep(names(cellinfo$markerlist)[nn], length(cellinfo$markerlist[[nn]])   )) %>% Reduce(c, .)
  
  ##making sure there arent deduped genes
  
  tff=!duplicated(cellinfo$markers)
  cellinfo$markers=cellinfo$markers[tff]
  marker.reference=marker.reference[tff]
  
  
  tfc=!duplicated(cellinfo$cells)
  cellinfo$cells=cellinfo$cells[tfc]
  cell.reference=cell.reference[tfc]
  
  ## reassemble deduplicated lists
  markernames=cellinfo$markerlist %>% names
  cellnames=cellinfo$cell.list %>% names
  
  cellinfo$markerlist=lapply(names(cellinfo$markerlist), function(x) cellinfo$markers[marker.reference==x]) %>% givename(., markernames)
  cellinfo$cell.list=lapply(names(cellinfo$cell.list), function(x) cellinfo$cells[cell.reference==x]) %>% givename(., cellnames)
  
  
  cellinfo$gaps.markers=lapply(1:length(cellinfo$markerlist), function(x) rep(x, length(cellinfo$markerlist[[x]]))) %>% Reduce(c, .) %>% diff %>% as.logical %>% which
  cellinfo$gaps.cells= lapply(1:length(cellinfo$cell.list), function(x) rep(x, length(cellinfo$cell.list[[x]]))) %>% Reduce(c, .) %>% diff %>% as.logical %>% which
  
  fcat("making annotations")
  if(!is.null(cellinfo$clustering.variable)){
    cell.group.label=cellinfo$clustering.variable
    
    if(zsoverlap(names(cellinfo$markerlist), names(cellinfo$cell.list))==1 ){
      fcat("Labels of cells and markers are shared.\n Assuming", cellinfo$clustering.variable,"as common variable for markers and cells")
      marker.group.label=cellinfo$clustering.variable
    }
  }else{
    fcat("No clustering variable. Checking metadata labels")
    if(!is.null(cellinfo$cell.metadata)){
      fcat("Cell metadata found")
      cell.group.label=cellinfo$cell.metadata[1]
    }
    
    if(!is.null(cellinfo$marker.metadata)){
      fcat("Marker metadata found")
      marker.group.label=cellinfo$marker.metadata[1]
    }
  }
  
  #if(is.null(cellinfo$cell.annotation)){}
  previous.cell.annotation=cellinfo$cell.annotation
  fresh.cell.annotation=as.data.frame(cell.reference) %>% giverownames(., cellinfo$cells) %>% givecolnames(., 1, cell.group.label)
  
  
  
  if(is.null(previous.cell.annotation) || length(intersect(rownames(previous.cell.annotation), cellinfo$cells))==0){
    cellinfo$cell.annotation= fresh.cell.annotation
  }else{
    #bring in the previous cell information dat for the new cells.
    othercols=setdiff(colnames(previous.cell.annotation), cell.group.label) %>% unique
    
    added.cell.annotation=previous.cell.annotation %>% getrows(., rownames(fresh.cell.annotation)) %>% select(all_of(othercols))
    
    cellinfo$cell.annotation=cbind(fresh.cell.annotation, added.cell.annotation)
    
  }
  
  ## annotation of markers
  previous.marker.annotation=cellinfo$marker.annotation
  fresh.marker.annotation= as.data.frame(marker.reference) %>% giverownames(., cellinfo$markers) %>% givecolnames(., 1, marker.group.label)
  
  if(is.null(previous.marker.annotation)|| length(intersect(rownames(previous.marker.annotation), cellinfo$markers))==0){
    cellinfo$marker.annotation= fresh.marker.annotation
  }else{
    #bring in the previous marker information dat for the new cells, except of the first column
    othercols.markers=setdiff(colnames(previous.marker.annotation), marker.group.label) %>% unique
    
    added.marker.annotation=previous.marker.annotation %>% getrows(., rownames(fresh.marker.annotation)) %>% select(all_of(othercols.markers))
    cellinfo$marker.annotation=cbind(fresh.marker.annotation, added.marker.annotation)
    
  }
  
  
  cellinfo$cell.metadata= colnames(cellinfo$cell.annotation)
  cellinfo$marker.metadata= colnames(cellinfo$marker.annotation)
  
  cellinfo
}

################################################################################
#update clustering variable
################################################################################
update.group.vars=function(cellinfo, global=NULL, cells.group=NULL, markers.group=NULL){
  if(!is.null(global)){
    cellinfo$clustering.variable=global
    cellinfo$cell.metadata[1]=global
    cellinfo$marker.metadata[1]=global
    subst=colnames(cellinfo$cell.annotation)
    subst[1]=global
    colnames(cellinfo$cell.annotation)=subst
    substm=colnames(cellinfo$marker.annotation)
    substm[1]=global
    colnames(cellinfo$marker.annotation)=substm
    
  }else{
    
    if(!is.null(cells.group)){
      cellinfo$clustering.variable=cells.group
      cellinfo$cell.metadata[1]=cells.group
      
      subst=colnames(cellinfo$cell.annotation)
      subst[1]=cells.group
      colnames(cellinfo$cell.annotation)=subst
      
    }
    
    if(!is.null(markers.group)){
      
      
      cellinfo$marker.metadata[1]=markers.group
      
      substm=colnames(cellinfo$marker.annotation)
      substm[1]=markers.group
      colnames(cellinfo$marker.annotation)=substm
      
    }
    
  }
  cellinfo
  
}

################################################################################
# seriate genes given a cellinfo object
################################################################################

seriategenes=function(so, cellinfo, assay="SCT", seriate.groups=F, link.markers.cells=T){
  
  #mat=tryCatch({so[cellinfo$markers, ][[assay]]@data %>% as.matrix}, function(e) {
  # warning("problem extracting assay data matrix. extracting counts instead")
  #so[st, dscells][["RNA"]]@counts %>% as.matrix
  #})
  
  mat=so[cellinfo$markers, ][[assay]]@data %>% as.matrix
  absent.genes= setdiff(cellinfo$markers, rownames(mat))
  vr= get.variant.rows(mat) 
  
  cellinfo= removemarkers(cellinfo, union(mat[!vr,] %>% rownames, absent.genes))
  mat=mat[vr, ] 
  clusmeans=c()
  ct=1
  fcat("Seriating genes")
  ordered.gene.list=lapply(1: length(cellinfo$markerlist), function(nn){
    
    st= cellinfo$markerlist[[nn]]
    if(!link.markers.cells){
      dscells=colnames(mat)
    }else{
      dscells=cellinfo$cell.list[[nn]]
    }
    if(length(st)>1){
      ss=seriation::get_order(
        #assumes an ordered correspondence between cell markers and cells in cellinfo
        seriation::seriate(scale((mat[st, dscells]) %>% t) %>%t)
        , method="PCA")
      clusmeans[ct]<<- mean(colMeans(scale((mat[st, dscells]) %>% t) %>%t)*(sapply(1:length(dscells), function(x) x**2, USE.NAMES=F)))
      ct<<-ct+1
      st[ss %>% unname]}else{
        
        clusmeans[ct]<<- mean(scale((mat[st, dscells]) %>% as.numeric)*(sapply(1:length(dscells), function(x) x**2, USE.NAMES=F)))
        
        st
      }
  }) %>% givename(., cellinfo$markerlist %>% names)
  
  
  if(seriate.groups){
    names(clusmeans)= 1:nk
    neworder=clusmeans %>% sort %>% names
    newnames=names(cellinfo$markerlist)[neworder]
    cellinfo$markerlist=lapply(neworder, function(tt) ordered.gene.list[[tt]]) %>% givename(., newnames)
  }else{
    
    cellinfo$markerlist=ordered.gene.list 
  }
  
  cellinfo.final=refresh.cellinfo(cellinfo)
  
  cellinfo.final
}



################################################################################
# cellinfo related functions
################################################################################
adjust.markerlist= function(markerlist, n) lapply(markerlist, function(x) removenas(x[1:n])) %>% givename(., markerlist %>% names)


adjust.markerlist.tail= function(markerlist, n) lapply(markerlist, function(x) removenas(x[(length(x)-n):length(x)])) %>% givename(., markerlist %>% names)



removenas=function(x) x[!is.na(x)]
adjust.markers= function(cellinfo, n){
  cellinfo$markerlist= adjust.markerlist(cellinfo.markerlist, n)
  cellinfo=refresh.cellinfo(cellinfo)
  cellinfo
}


################################################################################
# adjust cells (get top)
################################################################################

adjust.cells= function(cellinfo, n){
  cellinfo$cell.list=lapply(cellinfo$cell.list, function(x){
    
    if(n<=length(x)){
      removenas(x[(length(x)-n):length(x)])
    }else{
      x
    }
    
    
    
  }) %>% givename(., cellinfo$cell.list %>% names) 
  cellinfo=refresh.cellinfo(cellinfo)
  cellinfo
}


################################################################################
# cellinfo collect values for a variable inside the cellinfo
################################################################################

cellinfo.get.values= function(cellinfo, variable, so){
  
  metacols= so %>% metadata  %>% colnames
  genenames= rownames(so)
  
  if(variable %in% metacols){
    
    for(xx in names(cellinfo$cell.list)){
      
      cellinfo$values[[variable]][[xx]]= (so %>% metadata)[cellinfo$cell.list[[xx]], variable] 
    }
    
  }
  
  
  if(variable %in% genenames){
    
    
    metamat=join_meta_exp(so, genes=variable, assay=DefaultAssay(so))
    for(xx in names(cellinfo$cell.list)){ 
      cellinfo$values[[variable]][[xx]]=  metamat[cellinfo$cell.list[[xx]],variable ]
      
      
    }
    
  }
  cellinfo
}


################################################################################
# cellinfo single dataset workflow
################################################################################


cellinfo.workflow.singleds=function(so, grouping.variable="seurat_clusters", assay="RNA", markers=NULL, recalculate.markers=F, marker.params=NULL, colorlist=NULL, ncells=300, return.seurat=F){
  
  if(assay=="RNA"){
    fcat("Warning: assay is set as RNA. consider using a SCT assay for better visualisation") 
  }
  # housekeeping.
  #check that there is an overlap between categories of clusvar and categories of group1 in markers table  
  ###############################################################################
  #Step A. locate or recalculate markers. 
  ###############################################################################
  
  if(is.null(so$RNA@misc$top_markers) & is.null(markers)){
    fcat("no markers found. Please calculate markers and pass them to argunment markers or to so$RNA@misc$top_markers") 
    
  }
  ##############
  # Step b.  generate cellinfo with marker signatures. 
  #############
  if(is.null(colorlist)){
    clusterss=NULL  
  }else{
    clusterss=names(colorlist[[grouping.variable]])
    
  }
  
  cellinfo= seriatecells(so, clusvar=grouping.variable, clusters= clusterss, meth="seurat", extended.output=T, deduped=T, ncells=ncells)
  
  outs= cellinfo.heatmap(so, cellinfo,return.everything=T, genes.to.label=adjust.markers(cellinfo, 2)$markers, colorlist=colorlist)
  if(return.seurat){
    outs[["seuratobject"]]=so 
  }
  
  outs
  
}  

################################################################################
# cellinfo- create ranked cell plot... single variable
################################################################################

cellinfo.ordered.cell.plot= function(cellinfo, variable, scale.x=T, clusters=names(cellinfo$cell.list), colorlist=NULL, groupvar=cellinfo$clustering.variable){
  
  if(is.null(colorlist)){
    
    colorlist=list()
    
    catnames=names(cellinfo$cell.list)
    colorlist[[groupvar]]= randomcolors(length(catnames)) %>% givename(., catnames)
  }
  
  # structure data in df for ggplot
  
  
  plotdf=lapply(clusters, function(gr){  
    numcells= length(cellinfo$cell.list[[gr]])
    
    if(scale.x){
      tot=numcells 
    }else{
      tot=1 
    }
    
    
    cellranks=(1:numcells)/numcells
    vals=cellinfo$values[[variable]][[gr]]
    
    mt1=cellranks %>% as.data.frame %>% givecolnames(., nms="relative.rank.in.group")
    mt1[, variable]=vals
    mt2=giverownames(mt1, cellinfo$cell.list[[gr]]) %>% dplyr::mutate(group=!!gr)
    
    mt2    
  }) %>% Reduce(rbind, .)
  
  ggplot(plotdf)+geom_path(aes(x=relative.rank.in.group, y=!!sym(variable), color=group))+theme_classic()+xlab("Cell's rank in group")+NoLegend()+scale_color_manual(values=colorlist[[groupvar]])#+ylab(!!sym(variable))
  
}




################################################################################
#get markers
################################################################################

get.markers=function(so){
  
  so$RNA@misc$markers 
}

get.top.markers=function(so){
  
  so$RNA@misc$top_markers 
}

get.significant.markers=function(so){
  
  so$RNA@misc$significant_markers
}


################################################################################
# format gene vector for enrichr
################################################################################

format.enrichr=function(x) x %>% paste(., collapse="\n") %>% cat


################################################################################
# function to subset a cellinfo to a few subgroups.
################################################################################
subset.cellinfo=function(cellinfo, labels){
  fcat("restricting to labels only present in list")
  tf=lapply(labels, function(x) x %in% (cellinfo$markerlist %>% names)) %>% Reduce(c, .)
  labels2=labels[tf]
  
  cellinfo$markerlist=cellinfo$markerlist[labels2]
  cellinfo$cell.list=cellinfo$cell.list[labels2]
  cellinfo$input.cluster.ids=labels2
  refresh.cellinfo(cellinfo)
}

fillmat=function(genemat, xx){ 
  if(is.null(genemat[[xx]]))
  {genemat[[xx]]=0}
  
  return(genemat)}

add.cell.annotation= function(so, cellinfo, vars, assay="SCT", overwrite=T){
  
  gene.vars=vars[vars %in% rownames(so)]
  non.gene.vars=vars[!(vars %in% rownames(so))]
  meta.vars= non.gene.vars[non.gene.vars %in% (metadata(so) %>% colnames)]
  non.meta.vars= non.gene.vars[!(non.gene.vars %in% (metadata(so) %>% colnames))]
  
  absent.vars=intersect(non.gene.vars, non.meta.vars)
  if(length(absent.vars)>0){
    fcat("Warning: the following variables are absent from the dataset: ", paste(absent.vars, collapse=" "))
  }
  
  clls=cellinfo$cell.annotation %>% rownames
  
  
  new.annotation=metadata(so)[clls,] %>% select(all_of(meta.vars))
  #%>% givecolnames(., nms=meta.vars)
  
  #if overwrite is false, we only incorporate new columns
  if(overwrite==F){
    addition.vars=setdiff(meta.vars, cellinfo$metadata)
  }else{
    #if overwrite is true then we allow replacement of old variables
    addition.vars=meta.vars 
  }
  
  for(mv in addition.vars){
    cellinfo$cell.annotation[[mv]]= new.annotation[[mv]]
  }
  
  
  if(length(gene.vars)>0){
    genemat=join_meta_exp(so, genes=gene.vars, cells=clls)
    for(gn in gene.vars){
      genemat=fillmat(genemat, gn)
    }
    #%>%  givecolnames(., nms=gene.vars)
    cellinfo$cell.annotation=cellinfo$cell.annotation %>% cbind(., genemat[clls,] %>% select(all_of(gene.vars)) )
  }
  
  cellinfo$cell.metadata=colnames(cellinfo$cell.annotation)
  cellinfo
  
  
}

remove.cell.annotation=function(cellinfo, variables){
  
  cellinfo$cell.annotation= cellinfo$cell.annotation %>% select(-all_of(variables)) 
  
  cellinfo$cell.metadata=  setdiff(colnames(cellinfo$cell.annotation), variables)
  cellinfo 
}

remove.marker.annotation=function(cellinfo, variables){
  
  cellinfo$marker.annotation= cellinfo$marker.annotation %>% select(-all_of(variables)) 
  
  cellinfo$marker.metadata=  setdiff(colnames(cellinfo$marker.annotation), variables)
  cellinfo 
}

################################################################################
# remove genes from a cell info object
################################################################################
removemarkers=function(cellinfo, markers){
  nms=names(cellinfo$markerlist)
  cellinfo$markerlist=lapply(1:length(cellinfo$markerlist), function(x){
    gns=cellinfo$markerlist[[x]]
    
    gns[!(gns %in% markers)]
    
  }) %>% givename(., nms)
  
  refresh.cellinfo(cellinfo)
}


################################################################################

#szymkiewicz_simpson_coefficient to calculate size of overlaps between vectors
zsoverlap <- Vectorize(function(vector1, vector2) {
  
  # Calculate the size of the intersection
  intersection_size <- length(intersect(vector1 , vector2))
  #fcat(intersection_size)
  # Calculate the minimum size of the two sets
  if(length(vector1)==0 || length(vector1)==0){
    fcat("Warning: some of the sets are zero length")
  }
  min_size <- min(length(vector1), length(vector2))
  #fcat(min_size)
  # Calculate the Szymkiewicz–Simpson coefficient
  szymkiewicz_simpson <- intersection_size / min_size
  
  return(szymkiewicz_simpson)
}, USE.NAMES=F)


################################################################################
# filter cells based on a prediciton score for a mapping.
################################################################################

filterps= function(so, varbl="seurat_clusters", th) so[, so %>% metadata %>% filter(!!sym(paste0(mapfun(varbl), "_predicted.id.score"))>=!!th ) %>% rownames]



################################################################################
# selectively remove heatmap legends from a ComlpexHeatmap object
################################################################################


removelegends= function(hm, legends, where="top"){
  
  annoslot=paste0(where, "_annotation")
  #for(ann in c("top_annotation", "right_annotation", "left_annotation", "bottom_annotation")){
  for(leg in legends){
    slot(hm, annoslot)@anno_list[[leg]]@show_legend=F
  }
  
  hm 
  #}
}


################################################################################
#   classify cells into tf, ligands and receptors
################################################################################


simpleCache(paste_("genes", "db", "human_tfs"), {
# m <- fread("http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv")
#m[,.N,by=`Is TF?`]
#Is TF?    N
# 1:    Yes 1639
# 2:     No 1126
#ccbr.tfs <- m[`Is TF?`=="Yes", `HGNC symbol`] #paste0(`Ensembl ID`,"-",`HGNC symbol`)]

## obtained from homococo v11, from all the TFs containing matrices for which DNA can be scanned

hocomoco.tfs= c("AHR", "AIRE", "ALX1", "ALX3", "ALX4", "ANDR", "ANDR", "ANDR", 
                "AP2A", "AP2B", "AP2C", "AP2D", "ARI3A", "ARI5B", "ARNT2", "ARNT", 
                "ARX", "ASCL1", "ASCL2", "ATF1", "ATF2", "ATF2", "ATF2", "ATF3", 
                "ATF4", "ATF6A", "ATF7", "ATOH1", "BACH1", "BACH2", "BARH1", 
                "BARH2", "BARX1", "BARX2", "BATF3", "BATF", "BATF", "BC11A", 
                "BCL6B", "BCL6", "BHA15", "BHE22", "BHE23", "BHE40", "BHE41", 
                "BPTF", "BRAC", "BRAC", "BRCA1", "BSH", "CDC5L", "CDX1", "CDX2", 
                "CEBPA", "CEBPB", "CEBPD", "CEBPE", "CEBPG", "CEBPZ", "CENPB", 
                "CLOCK", "COE1", "COT1", "COT1", "COT2", "COT2", "CPEB1", "CR3L1", 
                "CR3L2", "CREB1", "CREB3", "CREB5", "CREM", "CRX", "CTCFL", "CTCF", 
                "CUX1", "CUX2", "CXXC1", "DBP", "DDIT3", "DLX1", "DLX2", "DLX3", 
                "DLX4", "DLX5", "DLX6", "DMBX1", "DMRT1", "DPRX", "DRGX", "DUX4", 
                "DUXA", "E2F1", "E2F2", "E2F3", "E2F4", "E2F4", "E2F5", "E2F6", 
                "E2F7", "E2F8", "E4F1", "EGR1", "EGR2", "EGR2", "EGR3", "EGR4", 
                "EHF", "ELF1", "ELF2", "ELF3", "ELF5", "ELK1", "ELK3", "ELK4", 
                "EMX1", "EMX2", "EOMES", "EPAS1", "ERG", "ERR1", "ERR2", "ERR3", 
                "ESR1", "ESR1", "ESR2", "ESR2", "ESX1", "ETS1", "ETS2", "ETV1", 
                "ETV2", "ETV3", "ETV4", "ETV5", "ETV6", "ETV7", "EVI1", "FEV", 
                "FEZF1", "FIGLA", "FLI1", "FLI1", "FOSB", "FOSL1", "FOSL2", "FOS", 
                "FOXA1", "FOXA2", "FOXA3", "FOXB1", "FOXC1", "FOXC2", "FOXD1", 
                "FOXD2", "FOXD3", "FOXF1", "FOXF2", "FOXG1", "FOXH1", "FOXI1", 
                "FOXJ2", "FOXJ3", "FOXJ3", "FOXK1", "FOXL1", "FOXM1", "FOXO1", 
                "FOXO3", "FOXO4", "FOXO6", "FOXP1", "FOXP2", "FOXP3", "FOXQ1", 
                "FUBP1", "GABPA", "GATA1", "GATA1", "GATA2", "GATA2", "GATA3", 
                "GATA4", "GATA5", "GATA6", "GBX1", "GBX2", "GCM1", "GCM2", "GCR", 
                "GCR", "GFI1B", "GFI1", "GLI1", "GLI2", "GLI3", "GLIS1", "GLIS2", 
                "GLIS3", "GMEB2", "GRHL1", "GRHL2", "GSC2", "GSC", "GSX1", "GSX2", 
                "HAND1", "HAND1", "HBP1", "HEN1", "HES1", "HES5", "HES7", "HESX1", 
                "HEY1", "HEY2", "HIC1", "HIC2", "HIF1A", "HINFP", "HLF", "HLTF", 
                "HMBX1", "HME1", "HME2", "HMGA1", "HMGA2", "HMX1", "HMX2", "HMX3", 
                "HNF1A", "HNF1B", "HNF1B", "HNF4A", "HNF4G", "HNF6", "HOMEZ", 
                "HSF1", "HSF1", "HSF2", "HSF4", "HTF4", "HXA10", "HXA11", "HXA13", 
                "HXA1", "HXA2", "HXA5", "HXA7", "HXA9", "HXB13", "HXB1", "HXB2", 
                "HXB3", "HXB4", "HXB6", "HXB7", "HXB8", "HXC10", "HXC11", "HXC12", 
                "HXC13", "HXC6", "HXC8", "HXC9", "HXD10", "HXD11", "HXD12", "HXD13", 
                "HXD3", "HXD4", "HXD8", "HXD9", "ID4", "IKZF1", "INSM1", "IRF1", 
                "IRF2", "IRF3", "IRF4", "IRF5", "IRF7", "IRF8", "IRF9", "IRX2", 
                "IRX3", "ISL1", "ISL2", "ISX", "ITF2", "JDP2", "JUNB", "JUND", 
                "JUN", "KAISO", "KAISO", "KAISO", "KLF12", "KLF13", "KLF14", 
                "KLF15", "KLF16", "KLF1", "KLF3", "KLF4", "KLF5", "KLF6", "KLF8", 
                "KLF9", "LBX2", "LEF1", "LHX2", "LHX3", "LHX4", "LHX6", "LHX8", 
                "LHX9", "LMX1A", "LMX1B", "LYL1", "MAFA", "MAFB", "MAFF", "MAFF", 
                "MAFG", "MAFG", "MAFK", "MAFK", "MAF", "MAF", "MAX", "MAZ", "MAZ", 
                "MBD2", "MCR", "MECP2", "MEF2A", "MEF2B", "MEF2C", "MEF2D", "MEIS1", 
                "MEIS1", "MEIS2", "MEIS3", "MEOX1", "MEOX2", "MESP1", "MGAP", 
                "MITF", "MIXL1", "MLXPL", "MLX", "MNX1", "MSX1", "MSX2", "MTF1", 
                "MXI1", "MXI1", "MYBA", "MYBB", "MYB", "MYCN", "MYC", "MYF6", 
                "MYNN", "MYOD1", "MYOD1", "MYOG", "MZF1", "NANOG", "NANOG", "NDF1", 
                "NDF2", "NF2L1", "NF2L2", "NFAC1", "NFAC1", "NFAC2", "NFAC3", 
                "NFAC4", "NFAT5", "NFE2", "NFIA", "NFIA", "NFIB", "NFIC", "NFIC", 
                "NFIL3", "NFKB1", "NFKB2", "NFYA", "NFYB", "NFYC", "NGN2", "NKX21", 
                "NKX22", "NKX23", "NKX25", "NKX28", "NKX31", "NKX32", "NKX61", 
                "NKX61", "NKX62", "NOBOX", "NOTO", "NR0B1", "NR1D1", "NR1D1", 
                "NR1H2", "NR1H3", "NR1H3", "NR1H4", "NR1H4", "NR1I2", "NR1I2", 
                "NR1I3", "NR1I3", "NR2C1", "NR2C2", "NR2E1", "NR2E3", "NR2F6", 
                "NR4A1", "NR4A2", "NR4A3", "NR5A2", "NR6A1", "NRF1", "NRL", "OLIG1", 
                "OLIG2", "OLIG2", "OLIG3", "ONEC2", "ONEC3", "OSR2", "OTX1", 
                "OTX2", "OVOL1", "OVOL2", "OZF", "P53", "P53", "P5F1B", "P63", 
                "P63", "P73", "P73", "PATZ1", "PATZ1", "PAX1", "PAX2", "PAX3", 
                "PAX4", "PAX5", "PAX6", "PAX7", "PAX8", "PBX1", "PBX1", "PBX2", 
                "PBX3", "PBX3", "PDX1", "PDX1", "PEBB", "PHX2A", "PHX2B", "PIT1", 
                "PITX1", "PITX2", "PITX3", "PKNX1", "PLAG1", "PLAL1", "PO2F1", 
                "PO2F2", "PO2F3", "PO3F1", "PO3F2", "PO3F3", "PO3F4", "PO4F1", 
                "PO4F2", "PO4F3", "PO5F1", "PO5F1", "PO6F1", "PO6F2", "PPARA", 
                "PPARA", "PPARD", "PPARG", "PPARG", "PRD14", "PRDM1", "PRDM4", 
                "PRDM6", "PRGR", "PRGR", "PROP1", "PROX1", "PRRX1", "PRRX2", 
                "PTF1A", "PTF1A", "PURA", "RARA", "RARA", "RARA", "RARB", "RARG", 
                "RARG", "RARG", "RAX2", "RELB", "REL", "REST", "RFX1", "RFX1", 
                "RFX2", "RFX2", "RFX3", "RFX4", "RFX5", "RFX5", "RHXF1", "RORA", 
                "RORG", "RREB1", "RUNX1", "RUNX2", "RUNX3", "RXRA", "RXRA", "RXRB", 
                "RXRG", "RX", "SALL4", "SCRT1", "SCRT2", "SHOX2", "SHOX", "SIX1", 
                "SIX2", "SMAD1", "SMAD2", "SMAD3", "SMAD4", "SMCA1", "SMCA5", 
                "SNAI1", "SNAI2", "SOX10", "SOX10", "SOX11", "SOX13", "SOX15", 
                "SOX17", "SOX18", "SOX1", "SOX21", "SOX2", "SOX2", "SOX3", "SOX4", 
                "SOX5", "SOX7", "SOX8", "SOX9", "SOX9", "SP1", "SP1", "SP2", 
                "SP2", "SP3", "SP4", "SP4", "SPDEF", "SPI1", "SPIB", "SPIC", 
                "SPZ1", "SRBP1", "SRBP2", "SRF", "SRY", "STA5A", "STA5B", "STAT1", 
                "STAT1", "STAT2", "STAT3", "STAT4", "STAT6", "STF1", "SUH", "TAF1", 
                "TAL1", "TAL1", "TBP", "TBR1", "TBX15", "TBX19", "TBX1", "TBX20", 
                "TBX21", "TBX2", "TBX3", "TBX4", "TBX5", "TCF7", "TEAD1", "TEAD2", 
                "TEAD3", "TEAD4", "TEF", "TF2LX", "TF65", "TF7L1", "TF7L2", "TFAP4", 
                "TFCP2", "TFDP1", "TFE2", "TFE3", "TFEB", "TGIF1", "TGIF2", "THA11", 
                "THAP1", "THA", "THA", "THB", "THB", "TLX1", "TWST1", "TWST1", 
                "TYY1", "TYY2", "UBIP1", "UNC4", "USF1", "USF2", "VAX1", "VAX2", 
                "VDR", "VDR", "VENTX", "VEZF1", "VEZF1", "VSX1", "VSX2", "WT1", 
                "WT1", "XBP1", "Z324A", "Z354A", "ZBED1", "ZBT14", "ZBT17", "ZBT18", 
                "ZBT48", "ZBT49", "ZBT7A", "ZBT7B", "ZBTB4", "ZBTB4", "ZBTB6", 
                "ZEB1", "ZEP1", "ZEP2", "ZF64A", "ZFHX3", "ZFP28", "ZFP42", "ZFP82", 
                "ZFX", "ZFX", "ZIC1", "ZIC2", "ZIC3", "ZIC4", "ZIM3", "ZKSC1", 
                "ZKSC3", "ZN121", "ZN134", "ZN134", "ZN136", "ZN140", "ZN143", 
                "ZN148", "ZN214", "ZN219", "ZN232", "ZN250", "ZN257", "ZN260", 
                "ZN263", "ZN263", "ZN264", "ZN274", "ZN281", "ZN282", "ZN317", 
                "ZN320", "ZN322", "ZN329", "ZN331", "ZN333", "ZN335", "ZN335", 
                "ZN341", "ZN341", "ZN350", "ZN350", "ZN382", "ZN384", "ZN394", 
                "ZN394", "ZN410", "ZN418", "ZN418", "ZN423", "ZN436", "ZN449", 
                "ZN467", "ZN490", "ZN502", "ZN524", "ZN528", "ZN547", "ZN549", 
                "ZN554", "ZN554", "ZN563", "ZN563", "ZN582", "ZN586", "ZN589", 
                "ZN652", "ZN667", "ZN680", "ZN708", "ZN708", "ZN713", "ZN740", 
                "ZN768", "ZN770", "ZN770", "ZN784", "ZN816", "ZN816", "ZNF18", 
                "ZNF41", "ZNF41", "ZNF76", "ZNF85", "ZNF85", "ZNF8", "ZSC16", 
                "ZSC22", "ZSC31", "ZSCA4")
#union(ccbr.tfs, hocomoco.tfs)  
hocomoco.tfs
}, assignToVar="tfList",  recreate=T)
istf= function(x) x %in% tfList 
#lrtable=readRDS(file.path(params$resourcepath, "cell_cell_interactions/CellTalkDB/human_lr_pair.rds"))

#isligand=function(x) x %in%  (lrtable %>% pull(ligand_gene_symbol))
#isreceptor=function(x) x %in%  (lrtable %>% pull(receptor_gene_symbol))
#islrpair=function(x) x %in%  (lrtable %>% pull(lr_pair))


colorgene=function(x) {
  
  if(istf(x)){return("black")}else{
    if(isligand(x)){return("blue")}else{
      if(isreceptor(x) || (x %in% c("CRABP1", "RARRES1", "RARRES2", "CRABP2"))){return("red")}else{return("grey")}
    }
  }}

colorgene2=function(x) {
  
  if(istf(x)){return("black")}else{return("grey")}}





#alosteric protein database
#https://mdl.shsmu.edu.cn/ASD2023Common/static_file/archive_2023/ASD_Release_202309_AS.tar.gz

################################################################################
#
################################################################################

divide_legend_columns <- function(plot, num_columns, legend.position="bottom") {
  # Get the current legend
  current_legend <- plot + theme(legend.position = legend.position)  # Adjust the position if needed
  
  # Modify the legend to have x columns
  modified_legend <- current_legend +
    guides(color = guide_legend(nrow = 1, ncol = num_columns))
  
  return(modified_legend)
}



################################################################################
# No Legend title
################################################################################
NoLegendTitle <- function(plot) {
  # Modify the legend to remove the title
  guides(guide_legend(title = NULL))
  
}
################################################################################

remove_x_axis <- function() {
  # Modify the theme to remove x-axis
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
}

################################################################################
# scale text size
################################################################################

scale.text=function(scaling_factor){
  theme(
    text = element_text(size =  scaling_factor),
    axis.title = element_text(size =  scaling_factor),
    axis.text = element_text(size = scaling_factor),
    axis.title.x = element_text(size = scaling_factor),
    axis.text.x = element_text(size =  scaling_factor),
    axis.title.y = element_text(size =  scaling_factor),
    axis.text.y = element_text(size =  scaling_factor),
    legend.text = element_text(size =  scaling_factor),
    legend.title = element_text(size =  scaling_factor),
    plot.title = element_text(size =  scaling_factor),
    plot.subtitle = element_text(size =  scaling_factor),
    plot.caption = element_text(size =  scaling_factor)
  )
}



scale.text.general=function(plt, sca){
  plt$theme$text$size= plt$theme$text$size*sca
  plt
}


################################################################################
# masking types in kameneva et al cell types
################################################################################
masktype=Vectorize(function(x) if(x =="SCP" || x=="sympathoblasts" || x=="mesenchyme"){return(x)}else{return("other")}, USE.NAMES=F)

masktype3=Vectorize(function(x) if(x =="SCP" || x=="sympathoblasts" || x=="mesenchyme" || x=="chromaffin"){return(x)}else{return("other")}, USE.NAMES=F)

################################################################################
#replace markers within a seurat object
################################################################################

replacemarkers=function(so, newmarkers){
  if(!is.null(newmarkers)){
    so$RNA@misc$top_markers=newmarkers 
    return(so)
  }else{
    return(so)
  }}

################################################################################
# whether a gene is in chr 17 or chr1
################################################################################
#   
#   if(!exists("gtftable")){ 
# gtfpath=params$gtfpath
# 
# gtftable=rtracklayer::import(gtfpath)
# }
# genepositions=gtftable %>% as.data.frame %>% select(seqnames, start, gene_name, type)%>% filter(type=="gene")
# 
# 
# get.seqnames=Vectorize(function(x) (genepositions %>% filter(gene_name==x))$seqnames %>% as.character, USE.NAMES=F)
# 
# is.in.17=Vectorize(function(x){ if(!is.na(x) && grepl("^chr17$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
# is.in.1=Vectorize(function(x){ if(!is.na(x) && grepl("^chr1$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
# is.in.2=Vectorize(function(x){ if(!is.na(x) && grepl("^chr2$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
# is.in.3=Vectorize(function(x){ if(!is.na(x) && grepl("^chr3$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
# is.in.4=Vectorize(function(x){ if(!is.na(x) && grepl("^chr4$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
# is.in.5=Vectorize(function(x){ if(!is.na(x) && grepl("^chr5$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
# is.in.6=Vectorize(function(x){ if(!is.na(x) && grepl("^chr6$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
# is.in.7=Vectorize(function(x){ if(!is.na(x) && grepl("^chr7$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
# is.in.8=Vectorize(function(x){ if(!is.na(x) && grepl("^chr8$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
# is.in.9=Vectorize(function(x){ if(!is.na(x) && grepl("^chr9$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
# is.in.10=Vectorize(function(x){ if(!is.na(x) && grepl("^chr10$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
# is.in.11=Vectorize(function(x){ if(!is.na(x) && grepl("^chr11$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
# is.in.12=Vectorize(function(x){ if(!is.na(x) && grepl("^chr12$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
# is.in.13=Vectorize(function(x){ if(!is.na(x) && grepl("^chr13$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
# is.in.14=Vectorize(function(x){ if(!is.na(x) && grepl("^chr14$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
# is.in.15=Vectorize(function(x){ if(!is.na(x) && grepl("^chr15$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
# is.in.16=Vectorize(function(x){ if(!is.na(x) && grepl("^chr16$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
# is.in.17=Vectorize(function(x){ if(!is.na(x) && grepl("^chr17$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
# is.in.18=Vectorize(function(x){ if(!is.na(x) && grepl("^chr18$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
# is.in.19=Vectorize(function(x){ if(!is.na(x) && grepl("^chr19$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
# is.in.20=Vectorize(function(x){ if(!is.na(x) && grepl("^chr20$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
# is.in.21=Vectorize(function(x){ if(!is.na(x) && grepl("^chr21$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
# is.in.22=Vectorize(function(x){ if(!is.na(x) && grepl("^chr22$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)



################################################################################
#
################################################################################

rmquery=Vectorize(function(x){if(grepl("query",x)){
  strsplit(x, split="query_")[[1]][2]
}else{x} }, USE.NAMES=F)



rmquery.matrix= function(x){
  x %>% names2col(., "cellid") %>% dplyr::mutate(cellid=rmquery(cellid)) %>% col2names(., "cellid")
}




################################################################################
#  remove duplicates from a vector
################################################################################

deduplicate= function(x) x[!duplicated(x)]


################################################################################
# quick print a matrix's contents
################################################################################

quickprint=function(mat, wide=F, limit=5) if(wide){ mat[1:limit, ]}else{mat[1:limit, 1:limit]}


################################################################################
# fit hill function 
################################################################################

hillfunction <- function(x, k, n, beta) beta*((x**n)/((k**n)+(x**n)))

#model1 <- nls(fluorI ~ eDecay(t,myA,myT), data=ExpData, start=list(myA=10,myT=5))



################################################################################
# gglabel (auxiliary side bar plot to make a legend
################################################################################

manual.legend=function(colorlist, pname, size=7, trailing.spaces=15, bar.width=.6){
  
  
  p=colorlist[[pname]] %>% rev
  
  pl=ggplot(1:length(p) %>% as.data.frame %>% dplyr::mutate(h=.1, colors=p, name=names(p)), aes(x=factor(.), y=factor(h)))+
    geom_col(aes(fill=factor(.)))+
    scale_fill_manual(values=namenums(p))+
    coord_cartesian( ylim=c(0,10))+
    theme_classic()+
    #geom_text(y=0.5, color="black", angle=90)+
    NoAxes()+
    NoLegend()+
    ggtitle(paste0(paste(rep(" ", trailing.spaces), collapse=""), pname))
  
  return(pl+geom_text(inherit.aes=T,aes(label=name), y=1.05, size=size, color="black", hjust=0)+coord_flip(ylim=c(0,1/bar.width)))
}

################################################################################
# print an extract of a matrix without getting into huge columns or huge rows overprinting trouble
################################################################################

quickview=function(mat, wide=F, limit=5) if(wide){ mat[1:limit, ]}else{mat[1:limit, 1:limit]}

################################################################################
# roc analysis to quantify the signal of targets within selected peaks
################################################################################
make.roc.curves=function(targets, targetpeaks, pca.data, classvar="celltype"){


####make a function that takes the values for a peak for all samples and calculates a roc curve


is.target= function(x) ifelse(x %in% targets, TRUE, FALSE)

fcat("collecting embeddings...")
#ae= get.all.embeddings(selected.pc)
ae= pca.data %>% dplyr::mutate(class=is.target(!!sym(classvar))) %>% select(class, !!sym(classvar))
sp= get.sample.peakvals(peakid=targetpeaks[, "Geneid"], samples= rownames(ae)) %>% cbind(ae, .)


################################################################################
# generate roc.curves based on the individual promoter-embedded peaks to dplyr::select those that on its own do a great job of dplyr::selecting target cells. 
################################################################################
fcat("starting roc curve loop...")
#selected.peaks=peak.gene.results[[selected.pc]][[region.selection]]$Geneid[1:total.intervals] # dplyr::select peaks with overlapping promoters
#selected.names=peak.gene.results[[selected.pc]][[region.selection]]$gene_symbol[1:total.intervals] # dplyr::select peaks with overlapping promoters
fcat("Applying ROC to intervals")
nl=0
allrocs=lapply(1:nrow(targetpeaks), function(xx){
pp=targetpeaks[xx, "Geneid"]

#trying new function to make roc calculation more efficient
#fcat("starting to calculate roc matrix")
rocmat=get.roc.matrix(sp, signal.col=pp, class.col="class")
#fcat("done calculating roc")
#rocmat=roc.curve(sp[, c(pp, "class")],step=0.001, true.value=TRUE, false.value=FALSE, color="black", line.width=1, plot.symbol=20,make.plot=F )
#dev.off()
au=trapezium.auc2(rocmat)
#fcat(pp, nn, au)
if(nl==50){
cat("\n")
  nl<<-0
}else{
  nl<<-nl+1
cat(".")
}

list(
  #plot=ggplot(rocmat %>% as.data.frame)+geom_path(aes(x=FPs, y=TPs))+ggtitle(nn)+geom_abline(slope=1, intercept=0, color="red", linetype="dashed")+theme_classic()+NoAxes(),
     auc=au, 
     Geneid=pp,
     gene_symbol=nn)



}) 
allrocs
}

################################################################################
# roc analysis to quantify the signal of a collective list of targets within selected peaks
################################################################################

make.roc.curves.all=function(targetlist, targetpeaks, pca.data){


####make a function that takes the values for a peak for all samples and calculates a roc curve


is.target= function(x, targets) fifelse(x %in% targets, TRUE, FALSE)

fcat("collecting embeddings...")
#ae= get.all.embeddings(selected.pc)

sp= get.sample.peakvals(peakid=targetpeaks[, "Geneid"], samples= rownames(pca.data)) %>% cbind(., pca.data %>% select(!contains("PC")))



################################################################################
# generate roc.curves based on the individual promoter-embedded peaks to dplyr::select those that on its own do a great job of dplyr::selecting target cells. 
################################################################################
fcat("starting roc curve loop...")
#selected.peaks=peak.gene.results[[selected.pc]][[region.selection]]$Geneid[1:total.intervals] # dplyr::select peaks with overlapping promoters
#selected.names=peak.gene.results[[selected.pc]][[region.selection]]$gene_symbol[1:total.intervals] # dplyr::select peaks with overlapping promoters
fcat("Applying ROC to intervals")
nl=0
allrocs=lapply(1:nrow(targetpeaks), function(xx){
  if(nl==0){
    cat(xx, ".", targetpeaks[xx, "Geneid"])
  }
alltargets=lapply(1:length(targetlist), function(xy){

classvar=targetlist[[xy]]$classvar
spt =sp %>% dplyr::mutate(class=is.target(!!sym(classvar), targetlist[[xy]]$targets))   


pp=targetpeaks[xx, "Geneid"]
#fcat(pp)
#trying new function to make roc calculation more efficient
#fcat("starting to calculate roc matrix")
rocmat=get.roc.matrix(spt, signal.col=pp, class.col="class")
#fcat("done calculating roc")
#rocmat=roc.curve(sp[, c(pp, "class")],step=0.001, true.value=TRUE, false.value=FALSE, color="black", line.width=1, plot.symbol=20,make.plot=F )
#dev.off()
au=trapezium.auc2(rocmat)
#fcat(pp, nn, au)


list(
  #plot=ggplot(rocmat %>% as.data.frame)+geom_path(aes(x=FPs, y=TPs))+ggtitle(nn)+geom_abline(slope=1, intercept=0, color="red", linetype="dashed")+theme_classic()+NoAxes(),
     auc=au, 
     Geneid=pp,
     gene_symbol=nn, 
     target=names(targetlist)[xy]
     )

}) %>% bind_rows
  
if(nl==50){
cat("\n")
  nl<<-0
}else{
  nl<<-nl+1
cat(".")
}
  alltargets
}) %>% bind_rows
allrocs
}


################################################################################
# perform hypeR gene set enrichments and plots
################################################################################

get.hyper.enrichments= function(markerlist, gsetpars=load.gene.signatures(), test="hypergeometric", statt="JACCARD", markersid="",plotting.format="spaced5", colorlist=NULL, numbars=5, title="", make.barplot=F, make.heatmap=F, background=23467){
  library(hypeR)
    library(stringr)
    if(is.null(colorlist)){
  colorlist=randomcolors(markerlist %>% names %>% length) %>% givename(.,markerlist %>% names )
  fcat(dput(colorlist))
    }
allhyperdfs=lapply(gsetpars, function(gp){
    hyperwtdf=NULL
    gset=gp$set
    gsetname=gp$name
    
    get.gsgenes= Vectorize(function(genesetname) gset$genesets[[genesetname]], USE.NAMES=F)
    
    
    
    #tryCatch({
    hyperwtdf=lapply(1:length(markerlist), function(y){
                       cat(gsetname, "gene set",names(markerlist)[y], "\n")
                       
                       
                       x=markerlist[[y]] 
                       if(length(x)==0){
                         return()}
                       hyp=hypeR::hypeR(signature=x, genesets=gset, test=test, plotting=TRUE, background=background)
                       
                       if (statt=="JACCARD"){
                         df=hyp$data %>% dplyr::mutate(jaccard_index= overlap/(signature+geneset-overlap)) %>% dplyr::arrange(-jaccard_index) 
                       }else{
                         if(statt=="PVAL"){
                           df=hyp$data %>% dplyr::mutate(jaccard_index= overlap/(signature+geneset-overlap)) %>% dplyr::arrange(pval) 
                           
                         }    else {
                           
                           if(statt=="FDR"){
                             df=hyp$data %>% dplyr::mutate(jaccard_index= overlap/(signature+geneset-overlap)) %>% dplyr::arrange(fdr) 
                             
                           }
                         
                         }
                         
                       }
                       df =df %>% dplyr::mutate(geneset.size=geneset) 
                       df =df %>% dplyr::mutate(geneset=gsetname) 
                       df$clusterlabel=names(markerlist)[y]
                       #df$cluster=y+seq(0, 1-(1/numbars), 1/numbars)
                       df$guide=1
                       df
                     }) %>% Reduce(rbind,.)
    
    if(!is.null(hyperwtdf)){
      hyperwtdf[, "label"]= gsub("_", " ", hyperwtdf[, "label"]) %>% str_to_sentence
      hyperwtdf= hyperwtdf %>% names2col(., "original.geneset")
      pparsets=list(
        tight=list(barw=0.98, tsz=200, scf=0.3, sca=1.3, tallf=8,fatf=5, parstring="tight"),
        spaced=list(barw=0.95, tsz=180, scf=0.3, sca=1.3, tallf=12,fatf=5, parstring="spaced"),
        spaced2=list(barw=0.95, tsz=180, scf=0.3, sca=1.3, tallf=12,fatf=15, parstring="superwide"),
        spaced3=list(barw=0.90, tsz=160, scf=0.3, sca=1.3, tallf=18,fatf=15, parstring="superwide"), 
        spaced4=list(barw=0.90, tsz=80, scf=0.1, sca=8, tallf=25,fatf=15, parstring="fit"),
        spaced5=list(barw=0.90, tsz=20, scf=0.1, sca=8, tallf=25,fatf=15, parstring="fit")
      )
      parsel=pparsets[[dataset.info$gsea_plot_params_hyper]]
      fatf=parsel$fatf
      tallf=parsel$tallf
      tsz=parsel$tsz
      gtsz=parsel$tsz*parsel$scf
      barw=parsel$barw
      sca=parsel$sca
      parstring=parsel$parstring
      allstops=hyperwtdf$cluster
      
      ### add some formatting to the labels
      
      hyperwtdf=hyperwtdf %>% dplyr::mutate(original.label=label, short.label=summarisename(label, dataset.info$gsea_total_print_length), label= str_replace(label, paste0(str_to_sentence(gsetname)," "),""))
      
      ##we use pdf because sometimes the text in the categories is very long and goes beyond the edge of the artboards which we unclip in illustrator. 
      # all genes in marker list
    }
    
  hyperwtdf
})

  allhyperdfs=allhyperdfs %>% bind_rows %>% dplyr::filter(overlap>0)
  allhyperdfs=allhyperdfs %>% dplyr::mutate(pval.direct=dhyper(overlap,geneset.size, background-geneset.size, signature)) %>% group_by(clusterlabel, geneset) %>% dplyr::mutate(padj=p.adjust(pval, method="BH"), tlog.pval=ifelse(pval.direct>=0.05, 0, -log(pval.direct)), tlog.padj=ifelse(padj<=0.05,-log(padj), 0)) %>% ungroup() %>% as.data.frame %>% dplyr::select(short.label, geneset, clusterlabel, overlap, jaccard_index, background, geneset.size, signature, hits,pval,pval.direct, padj, fdr, tlog.pval, tlog.padj,original.label)
}


  




################################################################################
# extract colors from a ggplot, or any other parameter for that matter
################################################################################
# getcolors=function(plt, type){
# plot_data <- ggplot_build(plt)
#   
#   # Extract the colors from the data
#   colors <- unique(unlist(lapply(plot_data$data, function(layer) {
#     layer[[type]]
#   })))
#   colors
# }

################################################################################
# find if an array is empty
################################################################################
isempty= Vectorize(function(x) length(x)==0, USE.NAMES=F)

################################################################################
# coerce a one dimensional named vector into a data frame
################################################################################

vector.as.df=function(mt, nm="value"){
 mt2=list(peakid=unname(mt))
 mt2=mt2[[1]] %>% data.frame %>% giverownames(., names(mt)) %>% givecolnames(., nms=nm)
 mt2
}


################################################################################
# add a random set of colors to colorlist if nonexistent
################################################################################


 create.colors.if.empty=function(colorlist, colorvar, dff, seed=42){
   set.seed(seed)
if(is.null(colorlist[[colorvar]])){
  caats= dff %>% pull(colorvar) %>% unique
 colorlist[[colorvar]]= randomcolors(length(caats)) %>% givename(., caats)
}
  set.seed(42)
    return(colorlist)
 }



################################################################################
# generate full paths
################################################################################

#generate.full.paths=function(meta) meta %>% dplyr::mutate(bampath=paste0(bamroot,Experiment, "_REP", replicate,".mLb.clN.sorted.bam"), bigwigpath=paste0(bamroot,"/bigwig/",Experiment, "_REP", replicate,".mLb.clN.bigWig"))
robust.home.path= function(x) gsub("~", "/home/rstudio", x)

generate.full.paths=function(meta) meta %>% dplyr::mutate( bampath=paste0(robust.home.path(bamroot),Experiment, "_REP", replicate,".mLb.clN.sorted.bam"), bigwigpath=paste0(robust.home.path(bamroot),"/bigwig/",Experiment, "_REP", replicate,".mLb.clN.bigWig"),bampath.merged=paste0(robust.home.path(nfcoreroot),"/bwa/merged_replicate/",gsub("-r[0-9]","", Experiment),".mRp.clN.sorted.bam"),  bigwigpath.merged=paste0(robust.home.path(bamroot),"/bwa/merged_replicate/bigwig/",gsub("-r[0-9]","", Experiment),".mRp.clN.bigWig"))

remove.rn= function(x)  gsub("-r[0-9]+", "", x)

################################################################################
# import metadata from dataset number j as described in the config
################################################################################

get.metadata.from.dataset=function(j, ids=NULL, exclude= dataset.info$exclude_cell_type, celltype.var=dataset.info$cell_type_variable, dataset.info){  
  fcat("1")
metadata.main<- fread(dataset.info$dataset_paths_metadata[j])

fcat("1.1")
metadata.main= metadata.main %>% as.data.frame(.,  optional=T) %>% fixcolnames %>% dplyr::mutate(dsname=dataset.info$dataset_name[j], bamroot=make.bampath(get.nfcore.root(j, dataset.info)), nfcoreroot=get.nfcore.root(j, dataset.info), dscategory=dscategory, dsdatatype="ATAC", Experiment=!!sym(dataset.info$metadata_sample_name_col[j]), Condition=!!sym(dataset.info$cell_group_vars[j]), !!sym(dataset.info$cell_group_vars[j]):=as.character(!!sym(dataset.info$cell_group_vars[j]))) 
if(!("replicate" %in% colnames(metadata.main))){
  metadata.main$replicate=1
}

if(!any(grepl("_REP[0-9]$", metadata.main[[dataset.info$metadata_sample_name_col[j]]]))){

  metadata.main[[dataset.info$metadata_sample_name_col[j]]]= paste0(metadata.main[[dataset.info$metadata_sample_name_col[j]]], "_REP", metadata.main$replicate)
   
}

metadata.main= metadata.main %>% col2names(., dataset.info$metadata_sample_name_col[j])

metadata.main<- metadata.main %>% generate.full.paths




if(!dataset.info$is_reference[j]){
  fcat("2")
 metadata.main= metadata.main %>% dplyr::mutate( Condition=paste0(!!sym(dataset.info$cell_group_vars[j]), "_ds",j) )
}
  

if(is.null(ids)){
  fcat("3")
rownames(metadata.main)=metadata.main[, dataset.info$metadata_sample_name_col[j]]
ids=rownames(metadata.main)
}
#retrieve quality control files
fcat("4")
if(any(grepl("multiqc", list.files(dataset.info$dataset_paths_nfcore[j])))){
  
metadata.qc= fread(file=paste0( dataset.info$dataset_paths_nfcore[j], "/multiqc/", config$peak_width, "/multiqc_data/multiqc_samtools_stats_1.txt")) %>% as.data.frame 
rownames(metadata.qc) <- metadata.qc %>% pull(Sample) # %>% splitsamplename


metadata.frip=import.multiqc.x.x.info(j, filename="multiqc_mlib_frip_score-plot.txt", varname="frip.score")
metadata.peak.count=import.multiqc.x.x.info(j, filename="multiqc_mlib_peak_count-plot.txt", varname="peak.count")

#bind all matrices together

metadata.complete <- bind_cols(metadata.main[ids,],metadata.qc[ids,], metadata.frip[ids, ,drop=F],metadata.peak.count[ids, ,drop=F] ) %>% giverownames(., nms=ids)


}else{
 warning(paste0("dataset.", dataset.info$dataset_name[j], ": Please verify the path containing the multiqc output. \nProceeding without global quality control measures"))
  metadata.frip <- NULL
  metadata.qc <- NULL
  metadata.peak.count<-NULL
  metadata.complete=metadata.main
}


metalist=list(
  metadata.main=metadata.main,
  metadata.frip=metadata.frip,
  metadata.qc=metadata.qc,
  metadata.complete=metadata.complete,
  metadata.peak.count=metadata.peak.count)

if(!is.null(exclude)){
metalist=metadata.exclude.groups(metalist, is.reference=dataset.info$is_reference[j], celltype.var=celltype.var, exclude= exclude)  
  
}
metalist

}


################################################################################
# exclude dataset groups from metadata
################################################################################
metadata.exclude.groups=function(metalist, is.reference, celltype.var=dataset.info$cell_type_variable, exclude= dataset.info$exclude_cell_type){
  
excluded.ids=NULL
if(is.reference){
    excluded.ids<-metalist$metadata.complete %>% filter(!!sym(celltype.var) %in% exclude) %>% rownames
 metadata.filtered<- metalist$metadata.complete %>% filter(!(!!sym(celltype.var) %in% exclude))
 included.ids=setdiff(rownames(metalist$metadata.complete), excluded.ids)

}else{
 metadata.filtered<-metalist$metadata.complete
 included.ids=rownames(metalist$metadata.complete)
  
   
}

list(metadata.complete=metadata.filtered,
  metadata.main=metalist$metadata.main[included.ids,],
  metadata.qc=metalist$metadata.qc[included.ids,], 
  metadata.frip=metalist$metadata.frip[included.ids,,drop=F],
  metadata.peak.count=metalist$metadata.peak.count[included.ids,,drop=F],
  excluded.ids=excluded.ids, 
  included.ids=included.ids)
}

################################################################################
# exclude datasets below a specific frip score
################################################################################
metadata.filter.frip=function(metalist, frip.threshold=0.05){
  
excluded.ids=NULL

    excluded.ids<-metalist$metadata.complete %>% filter(frip.score<frip.threshold) %>% rownames
 metadata.filtered<- metalist$metadata.complete %>% filter(frip.score>=frip.threshold)
 included.ids=setdiff(rownames(metalist$metadata.complete), excluded.ids)



list(metadata.complete=metadata.filtered,
  metadata.main=metalist$metadata.main[included.ids,],
  metadata.qc=metalist$metadata.qc[included.ids,], 
  metadata.frip=metalist$metadata.frip[included.ids,,drop=F],
  metadata.peak.count=metalist$metadata.peak.count[included.ids, , drop=F],
  excluded.ids=metalist$excluded.ids,
  excluded.ids.frip=excluded.ids,
  included.ids=included.ids)
}

################################################################################
# expand colors based on an updated dataframe that presumably describes more categories than those present in the colorlist
################################################################################

expand.colors=function(colorlist, df, varr){ 
  if(varr %in% colnames(df) && length(unique(df[,varr])>0)){
  extracats=setdiff(as.character(unique(df[, varr])), names(colorlist[[varr]])) %>% removenas() ;
  if(length(extracats)>0){
  colorlist[[varr]]=c(colorlist[[varr]], randomcolors(length(extracats)) %>% givename(., extracats)); 
  }
  }
  colorlist 
  
  }

update.all.colors=function(colorlist, df){
for(j in names(colorlist)){
colorlist=expand.colors(colorlist, df, j)
}
  colorlist
}


################################################################################
# import multiqc x,x matrix information
################################################################################


import.multiqc.x.x.info= function(j, filename="multiqc_mlib_peak_count-plot.txt", varname="peak.count"){
  
  
  ff= fread(file=paste0( dataset.info$dataset_paths_nfcore[j], "/multiqc/", config$peak_width, "/multiqc_data/", filename)) %>% as.matrix 
ff <- ff %>% giverownames(., ff[, "Sample"]) 
ff<-ff[, -1]
colnames(ff) <- rownames(ff)
ff <- lapply(rownames(ff), function(x){l=list(); l[[x]]=as.numeric(ff[x,x]); l}) %>% unlist %>% as.data.frame %>% givecolnames(., nms=varname) 
ff

}
  
  
  
################################################################################
# import counts
################################################################################



  
   
  import.counts=function(j, dataset.info, recreate=F){ 
    simpleCache(paste0("count_list_dataset", dataset.info$dataset_name[j], "_path_", digest::digest(c(dataset.info$dataset_paths_counts_self[j],dataset.info$dataset_paths_counts_ref[j]) ) ) %>% addversion, {

counts.to.ref= fread(file= dataset.info$dataset_paths_counts_ref[j], skip=1) %>% as.data.frame %>% col2names(., "Geneid")  
counts.to.self= fread(file=dataset.info$dataset_paths_counts_self[j], skip=1) %>% as.data.frame %>% col2names(., "Geneid")


colnames(counts.to.self)= counts.to.self %>% colnames %>% clean.extension.bam
colnames(counts.to.ref)= counts.to.ref %>% colnames %>% clean.extension.bam

list(counts.to.self=counts.to.self, counts.to.ref=counts.to.ref)

    }, assignToVar="counts.list", recreate=recreate)
 counts.list     
      
  }



################################################################################
#
################################################################################




################################################################################
# functions to  obtain embeddings
################################################################################

get.all.embeddings= function(pc) pca.data %>% dplyr::arrange(!!sym(pc)) %>% dplyr::select(!!sym(pc), !!sym(classvar))
get.sample.peakvals= function(peakid, samples=NULL){
  
  if(!is.null(samples)){ mt=norm.peak.vals.t[ samples, peakid]}else{ mt=norm.peak.vals.t[ !colnames(norm.peak.vals) %in% c("Geneid", "Chr", "Start", "End", "Strand", "Length"),peakid ]}

if(length(peakid)==1){
  mt2=vector.as.df(mt, peakid)
}else{
rn=rownames(mt)  
mt2=apply(mt, 2, as.numeric)
rownames(mt2)=rn
}
rm(mt)
gc()
mt2
}


################################################################################
# functions to make roc curves
################################################################################

make.roc.curves=function(targets, targetpeaks, pca.data, classvar="celltype"){


####make a function that takes the values for a peak for all samples and calculates a roc curve


is.target= function(x) ifelse(x %in% targets, TRUE, FALSE)

fcat("collecting embeddings...")
#ae= get.all.embeddings(selected.pc)
ae= pca.data %>% dplyr::mutate(class=is.target(!!sym(classvar))) %>% select(class, !!sym(classvar))
sp= get.sample.peakvals(peakid=targetpeaks[, "Geneid"], samples= rownames(ae)) %>% cbind(ae, .)


################################################################################
# generate roc.curves based on the individual promoter-embedded peaks to dplyr::select those that on its own do a great job of dplyr::selecting target cells. 
################################################################################
fcat("starting roc curve loop...")
#selected.peaks=peak.gene.results[[selected.pc]][[region.selection]]$Geneid[1:total.intervals] # dplyr::select peaks with overlapping promoters
#selected.names=peak.gene.results[[selected.pc]][[region.selection]]$gene_symbol[1:total.intervals] # dplyr::select peaks with overlapping promoters
fcat("Applying ROC to intervals")
nl=0
allrocs=lapply(1:nrow(targetpeaks), function(xx){
pp=targetpeaks[xx, "Geneid"]

#trying new function to make roc calculation more efficient
#fcat("starting to calculate roc matrix")
rocmat=get.roc.matrix(sp, signal.col=pp, class.col="class")
#fcat("done calculating roc")
#rocmat=roc.curve(sp[, c(pp, "class")],step=0.001, true.value=TRUE, false.value=FALSE, color="black", line.width=1, plot.symbol=20,make.plot=F )
#dev.off()
au=trapezium.auc2(rocmat)
#fcat(pp, nn, au)
if(nl==50){
cat("\n")
  nl<<-0
}else{
  nl<<-nl+1
cat(".")
}

list(
  #plot=ggplot(rocmat %>% as.data.frame)+geom_path(aes(x=FPs, y=TPs))+ggtitle(nn)+geom_abline(slope=1, intercept=0, color="red", linetype="dashed")+theme_classic()+NoAxes(),
     auc=au, 
     Geneid=pp,
     gene_symbol=nn)



}) 
allrocs
}

################################################################################
# functions to make roc curves for all targets
################################################################################

make.roc.curves.all=function(targetlist, targetpeaks, pca.data){


####make a function that takes the values for a peak for all samples and calculates a roc curve

is.target= function(x, targets) fifelse(x %in% targets, TRUE, FALSE)

fcat("collecting embeddings...")
#ae= get.all.embeddings(selected.pc)

sp= get.sample.peakvals(peakid=targetpeaks[, "Geneid"], samples= rownames(pca.data)) %>% cbind(., pca.data %>% select(!contains("PC")))



################################################################################
# generate roc.curves based on the individual promoter-embedded peaks to dplyr::select those that on its own do a great job of dplyr::selecting target cells. 
################################################################################
fcat("starting roc curve loop...")
#selected.peaks=peak.gene.results[[selected.pc]][[region.selection]]$Geneid[1:total.intervals] # dplyr::select peaks with overlapping promoters
#selected.names=peak.gene.results[[selected.pc]][[region.selection]]$gene_symbol[1:total.intervals] # dplyr::select peaks with overlapping promoters
fcat("Applying ROC to intervals")
nl=0
allrocs=lapply(1:nrow(targetpeaks), function(xx){
  #if(nl==0){
  #  cat(xx, ".", targetpeaks[xx, "Geneid"])
  #}
alltargets=lapply(1:length(targetlist), function(xy){

classvar=targetlist[[xy]]$classvar
spt =sp %>% dplyr::mutate(class=is.target(!!sym(classvar), targetlist[[xy]]$targets))


pp=targetpeaks[xx, "Geneid"]
#fcat(pp)
#trying new function to make roc calculation more efficient
#fcat("starting to calculate roc matrix")
rocmat=get.roc.matrix(spt, signal.col=pp, class.col="class")
#fcat("done calculating roc")
#rocmat=roc.curve(sp[, c(pp, "class")],step=0.001, true.value=TRUE, false.value=FALSE, color="black", line.width=1, plot.symbol=20,make.plot=F )
#dev.off()
au=trapezium.auc2(rocmat)
#fcat(pp, nn, au)


list(
  #plot=ggplot(rocmat %>% as.data.frame)+geom_path(aes(x=FPs, y=TPs))+ggtitle(nn)+geom_abline(slope=1, intercept=0, color="red", linetype="dashed")+theme_classic()+NoAxes(),
     auc=au, 
     Geneid=pp,
     #gene_symbol=nn, 
     target=names(targetlist)[xy]
     )

}) %>% bind_rows
  
#if(nl==50){
#cat("\n")
#  nl<<-0
#}else{
#  nl<<-nl+1
#cat(".")
#}
  alltargets
}) %>% bind_rows
allrocs
}




make.roc.curves.fast=function(sample.peak.values, pca.data){


####make a function that takes the values for a peak for all samples and calculates a roc curve

is.target= function(x, targets) fifelse(x %in% targets, TRUE, FALSE)

fcat("collecting embeddings...")
#ae= get.all.embeddings(selected.pc)

sp= get.sample.peakvals(peakid=targetpeaks[, "Geneid"], samples= rownames(pca.data)) %>% cbind(., pca.data %>% select(!contains("PC")))



################################################################################
# generate roc.curves based on the individual promoter-embedded peaks to dplyr::select those that on its own do a great job of dplyr::selecting target cells. 
################################################################################
fcat("starting roc curve loop...")
#selected.peaks=peak.gene.results[[selected.pc]][[region.selection]]$Geneid[1:total.intervals] # dplyr::select peaks with overlapping promoters
#selected.names=peak.gene.results[[selected.pc]][[region.selection]]$gene_symbol[1:total.intervals] # dplyr::select peaks with overlapping promoters
fcat("Applying ROC to intervals")
nl=0
allrocs=lapply(1:nrow(targetpeaks), function(xx){
  #if(nl==0){
  #  cat(xx, ".", targetpeaks[xx, "Geneid"])
  #}
alltargets=lapply(1:length(targetlist), function(xy){

classvar=targetlist[[xy]]$classvar
spt =sp %>% dplyr::mutate(class=is.target(!!sym(classvar), targetlist[[xy]]$targets))


pp=targetpeaks[xx, "Geneid"]
#fcat(pp)
#trying new function to make roc calculation more efficient
#fcat("starting to calculate roc matrix")
rocmat=get.roc.matrix(spt, signal.col=pp, class.col="class")
#fcat("done calculating roc")
#rocmat=roc.curve(sp[, c(pp, "class")],step=0.001, true.value=TRUE, false.value=FALSE, color="black", line.width=1, plot.symbol=20,make.plot=F )
#dev.off()
au=trapezium.auc2(rocmat)
#fcat(pp, nn, au)


list(
  #plot=ggplot(rocmat %>% as.data.frame)+geom_path(aes(x=FPs, y=TPs))+ggtitle(nn)+geom_abline(slope=1, intercept=0, color="red", linetype="dashed")+theme_classic()+NoAxes(),
     auc=au, 
     Geneid=pp,
     #gene_symbol=nn, 
     target=names(targetlist)[xy]
     )

}) %>% bind_rows
  
#if(nl==50){
#cat("\n")
#  nl<<-0
#}else{
#  nl<<-nl+1
#cat(".")
#}
  alltargets
}) %>% bind_rows
allrocs
}


################################################################################
# coerce to matrix (similar to drop)
################################################################################


coerce.to.matrix=function(vec){
 if(is.vector(vec)){ 
  matrix(vec, ncol = 1, dimnames = list(names(vec), NULL))
 }else{
  vec 
 }
}

################################################################################

################################################################################
isempty.list= function(ll) lapply(ll , function(x) ifelse(length(x)<=1, T, F)) %>% unlist


################################################################################
# get gene set genes
################################################################################



################################################################################
# making enrichment heatmap
################################################################################
make.enrichment.heatmap=function(x, hyperdf, gsetpars=load.gene.signatures(), setid=NULL, cluster=F){

  
     if(is.null(setid)){
       
       setid=hyperdf %>% pull(clusterlabel) %>% unique %>% paste(., collapse="-")
     }
  
      fct=hyperdf$original.geneset %>% unique
      
      get.gsgenes= Vectorize(function(genesetname, gsetgroup) gsetpars[[gsetgroup]]$set$genesets[[genesetname]], USE.NAMES=F)
      
      get.genesetgroup= function(gsname, hyperdf=hyperdf) as.data.frame(hyperdf) %>% dplyr::filter(original.geneset==!!gsname) %>% pull(geneset) %>% unique
      
        if(length(x)>0 && length(fct)>0 ){
          #fcat("x has more than ozero elements")
        gsgenemat=matrix(nrow=length(x), ncol=length(fct), data=0)
      rownames(gsgenemat)=x
      colnames(gsgenemat)=fct
      
      for(fc in 1:nrow(hyperdf)){
        
        fcgenes=get.gsgenes(hyperdf[fc, "original.geneset"],hyperdf[fc, "geneset"] ) %>% as.vector
        
         if(!is.null(fcgenes)){
           #fcat("fcgenes is not null. making dataframe")
           odf=list(gene=x, geneset=fc, is.in.markerlist=(x %in% fcgenes)) %>% data.frame
           podf= odf %>% dplyr::filter(is.in.markerlist)
           #fcat("rows of podf", nrow(podf), "...")
           if(nrow(podf)>0){
           for(j in 1:nrow(podf)){
           gsgenemat[podf[j, "gene"], podf[j, "geneset"]]<-1
             
             
             }
           }
         }
      }
        
      
        }
           
            if(length(dim( remove.zerorows(gsgenemat))) <=1  || length(dim( remove.zerocols(gsgenemat))) <=1 ){
            warning("problem with final matrix dimensions (either too many 0's or wrong dimensions) Heatmap cannot be created") 
         
         hm=NULL 
         gsgenemat=NULL
         original.gset.names=NULL
       }else{
         
         
         
         
         gsgenemat=remove.zerocols(remove.zerorows(gsgenemat))
        
       }
        
      if(length(dim(gsgenemat)) >1){
        
      fcat("getting gset names")
      original.gset.names=colnames(gsgenemat)     
      colnames(gsgenemat)= lapply(colnames(gsgenemat), function(x) paste0(get.genesetgroup(x, hyperdf), "-", x)) %>% Reduce(c, .)
      
     # fcat("drawing heatmap...")
      colnames(gsgenemat)=summarisename(colnames(gsgenemat), dataset.info$gsea_total_print_length)
      #quickprint(gsgenemat) 
      hm=ComplexHeatmap::pheatmap(gsgenemat, cluster_row=cluster, cluster_col=cluster, border_color = NA, main=paste0("gene set ", setid, ": hit incidence in all gene sets"))
      #hm=NULL  #temporary patch
       
      }else{
      
         warning("either the markerlist or the gene set list is zero. Heatmap cannot be created.") 
        
      
         hm=NULL 
         gsgenemat=NULL
         original.gset.names=NULL
        }
      list(heatmap=hm, matrix=gsgenemat, original.gset.names=original.gset.names)
      
}


################################################################################
# scrutinise the gene sets for a keyword
################################################################################
get.relevant.genesets=function( keyword, gsetpars, exclude="placeholder"){
  
allgsetnames=lapply(gsetpars, function(x) names(x$set$genesets)) %>% Reduce(c, .) 
allgsets=lapply(gsetpars, function(x) x$set$genesets)  %>% Reduce(append, .) %>% givename(., allgsetnames)

relevant.gsetnames=allgsetnames[grepl(keyword, allgsetnames, ignore.case=T)]
relevant.gsetnames=relevant.gsetnames[!grepl(exclude, relevant.gsetnames)]

list(gene.sets=allgsets[relevant.gsetnames], relevant.geneset.names=relevant.gsetnames)
}



get.gene.gsea.relevance=function(genes,keyword,gsetpars,exclude="placeholder", label="",simple.output=T){

  geneset.search=get.relevant.genesets( keyword, gsetpars, exclude=exclude)


all.incidences=lapply(genes, function(gene) (lapply(geneset.search$gene.sets, function(x) gene %in% x) %>% unlist %>% sum)/ length(geneset.search$gene.sets)) %>% unlist() %>% givename(., genes)

if(simple.output){
all.incidences   
}else{
  list(incidence.percentage=all.incidences, relevant.genesets=geneset.search, df=coerce.vector.to.matrix(all.incidences) %>% givecolnames(., nms=paste0("incidence.in.", label, ".sets")) %>% as.data.frame )
}
  
}

################################################################################
#
################################################################################

################################################################################
#
################################################################################


coerce.vector.to.matrix=function(vec)  matrix(vec, ncol = 1, dimnames = list(names(vec), NULL))



################################################################################
# analysis of whether some peaks are in the range of the reference
################################################################################

range.analysis.workflow=function(target.peaks,target_cell_type, test.ids,vsd, pca.data , mtf_mat_gathered, sample.arrangement, keyword=config$gsea_celltype_searchterm, label=config$target_cell_type, exclude=config$gsea_celltype_exclude,enrichments=NULL,below.percentile=.1, pct.always.below=1){


target.ids=pca.data %>% filter(celltype==!!target_cell_type) %>% pull(Experiment)

test.vals=vsd[target.peaks , test.ids ]
ref.vals=vsd[target.peaks , target.ids]

################################################################################
# direct method
################################################################################
#rangetests=lapply( 1:length(target.peaks), function(j) {
#fcat(j, ":", peak.collection$roc[j])
 # rng=range(ref.vals[j, ])
  
  #isinrange= function(x, rng) ifelse(x>rng[1] && x<=rng[2], 0, isabove(x, rng))
 # isabove= function(x, rng) ifelse(x>rng[2], 1, isbelow(x, rng))
  #isbelow= function(x, rng) ifelse(x<=rng[1], -1, NA)
  #isexpressed=function(x) ifelse(x>=min.expression, -1, -2)
  
 # result=lapply(test.vals[target.peaks[j],], function(y) isinrange(y, rng)) %>% unlist()
 # result
  
#}) %>% Reduce(rbind, .) %>% givecolnames(.,nms= colnames(test.vals)) %>% giverownames(., target.peaks)

################################################################################
#percentile method
################################################################################
rangetests=lapply( 1:length(target.peaks), function(j) {

  ec=ecdf(ref.vals[j, ])
  result=ec(test.vals[target.peaks[j],test.ids])
  result=ifelse(result>=below.percentile, -1,0)
  
}) %>% Reduce(rbind,.) %>% givecolnames(.,nms= colnames(test.vals)) %>% giverownames(., target.peaks)

################################################################################
#  prepare the enrichments
################################################################################

howmany.peaks.in.range=apply(rangetests[target.peaks, test.ids], 2, function(x) sum(x==0))
howmany.peaks.below.range=apply(rangetests[target.peaks, test.ids], 2, function(x) sum(x==-1)) 
howmany.peaks.above.range=apply(rangetests[target.peaks, test.ids], 2, function(x) sum(x==1))

peak.range.count.df=list(in.range=howmany.peaks.in.range, below=howmany.peaks.below.range, above=howmany.peaks.above.range) %>% as.data.frame %>% cbind(., pca.data[names(howmany.peaks.in.range), ]) %>% pivot_longer(., cols=c("in.range", "below", "above"), names_to="peak.class", values_to="peak.class.count")

#ggplot(peak.range.count.df, aes(x=tmnt.full, y=peak.class.count, fill=peak.class))+geom_col(position="fill")

peak.freq.below.range=apply(rangetests[target.peaks, test.ids], 1, function(x) sum(x==-1)/length(x)) %>% sort(., decreasing=T)
peaks.always.below.range=peak.freq.below.range[peak.freq.below.range>= pct.always.below] %>% names

peak.freq.in.range=apply(rangetests[target.peaks,test.ids], 1, function(x) sum(x==0)/length(x)) %>% sort(., decreasing=T)
peaks.always.in.range=peak.freq.in.range[peak.freq.in.range==1] %>% names

peak.freq.above=apply(rangetests[target.peaks,test.ids], 1, function(x) sum(x==1)/length(x)) %>% sort(., decreasing=T)
peaks.always.above=peak.freq.in.range[peak.freq.above==1] %>% names

peak.freq.in.or.above=apply(rangetests[target.peaks,test.ids], 1, function(x) (sum(x==0)+sum(x==1))/length(x)) %>% sort(., decreasing=T)
peaks.in.or.above=peak.freq.in.or.above[peak.freq.in.or.above==1] %>% names

#best.sample.peaks.in.range=names(howmany.peaks.in.range)[1]
#which.peaks.in.range=rangetests[all.peaks.celltypes.auc[[target_cell_type]], best.sample.peaks.in.range]==0
#which.peaks.below.range=rangetests[all.peaks.celltypes.auc[[target_cell_type]], best.sample.peaks.in.range]==-1


#peaks.in.range= rownames(rangetests[which.peaks.in.range,])


################################################################################
# make heatmaps for the different classes of peaks.
################################################################################
#categorise peaks with many sorts of active and inactive peaks
fcat("organising peaks according to sensitivity")
above.partitions=lapply(seq(1, 0, -0.05), function(x)  peak.freq.above[peak.freq.above>x] %>% names) %>% Reduce(c, .)
in.range.partitions=lapply(seq(1, 0, -0.05), function(x)  peak.freq.in.range[peak.freq.in.range>x] %>% names) %>% Reduce(c, .)
below.partitions=lapply(seq(0, 1, 0.05), function(x)  peak.freq.below.range[peak.freq.below.range<x] %>% names) %>% Reduce(c, .)
pks=c(above.partitions, in.range.partitions, below.partitions) %>% unique

#define boundary between fractional peaks and eternally low peaks. 
gaps.peaks= length(pks)

pks.complete=c(pks, peaks.always.below.range)

#order arranged by number of peaks in range
smpls=names(howmany.peaks.in.range) %>% rev

# matrix subset to the relevant peaks being shown
range.mat=rangetests[pks.complete, ]

#contingency for lack of a pre-arrangement of samples
if(!is.null(sample.arrangement)){
  sample.arrangement=colnames(range.mat)
}

#seriating the matrix
range.mat.row.ord=seriation::get_order(seriation::seriate(range.mat, margin=1, method="PCA"))
range.mat.col.ord=seriation::get_order(seriation::seriate(range.mat, margin=2, method="PCA"))

fcat("plotting sensitive/insensitive peaks in different ways...")

br=seq(-1,1,1)  ## numeric breaks of the color bins
#colls=c("aquamarine4", "aquamarine", "darkgoldenrod1")  
colls=c("cyan4", "cyan3", "white")  

sca=3
tpdf(path=config$plotpath, paste0("heatmap_inrange_to_", target_cell_type, "rowseriation_pct.threshold_",below.percentile ) %>% addversion, wi=pw*sca*.6, he=pw*sca)
print(ComplexHeatmap::pheatmap(range.mat[range.mat.row.ord, sample.arrangement],
  cluster_row=F,
  cluster_col=F,
  annotation_col=pca.data[sample.arrangement %>% rev,dataset.info$annotvars],
  show_rownames=F,
  show_colnames=F,
  annotation_colors=allcolors,
  use_raster=T,
  breaks=br, 
  col=colls, 
  gaps_row=gaps.peaks 
  ) )
dev.off()

fcat("printing peaks to screen")
peaklist=list(inrange=pks, below.range=peaks.always.below.range) 
#fcat(peaklist %>% dput)

lapply(1:length(peaklist), function(pk){
  
  fwrite(norm.peak.vals[peaklist[[pk]], 2:4] %>% unname %>% giverownames(., NULL), file=paste0("~/metadata/", names(peaklist)[pk], ".bed"), sep="\t", row.names=F)
})


################################################################################
#  calculate the enrichments
################################################################################

fcat("calculating enrichments of sensitive/insensitive peaks...")
tfbs.enrichment.list=list(
#highly.matched.peaks=test.tfbs.enrichment(peak.freq.in.or.above[peak.freq.in.or.above>dataset.info$peak_sensitivity_threshold] %>% names, mtf_mat_gathered, target.name="highly.matched.peaks"),
below.range=test.tfbs.enrichment(peak.freq.below.range[peak.freq.below.range==1] %>% names, mtf_mat_gathered, target.name="always.below.range"),
#often.below.range=test.tfbs.enrichment(peak.freq.below.range[peak.freq.below.range>=dataset.info$peak_sensitivity_threshold] %>% names, mtf_mat_gathered, target.name="often.below.range"),
#high.treatment.sensitivity=test.tfbs.enrichment(peak.freq.below.range[peak.freq.below.range>=1-dataset.info$peak_sensitivity_threshold & peak.freq.below.range<=dataset.info$peak_sensitivity_threshold] %>% names, mtf_mat_gathered, target.name="high.treatment.sensitivity") 
  )




fcat("making data frame...")
tfbs.enrichment.df=tfbs.enrichment.list %>% bind_rows()
rownames(tfbs.enrichment.df)=NULL


################################################################################
# prepare contrasts for venn diagram
################################################################################

tfbs.enrichment.list2=tfbs.enrichment.list[lapply(tfbs.enrichment.list, function(x)  !is.null(x)) %>% Reduce(c, .)]
fcat("Making Venn Diagrams of TFs associated with each sensitivity group...")

mtfs.contrast=lapply(tfbs.enrichment.list2, function(x) x %>% pull(gene_symbol))  %>% givename(., names(tfbs.enrichment.list2))


gven=ggVennDiagram(mtfs.contrast, label_alpha=0 )
#ggplot2::scale_fill_gradient(low="blue",high = "pink")


fcat("extracting exclusive tfs for each sensitivity group...")


list.tfbs.exclusive=lapply(1:length(mtfs.contrast), function(x){
  setdiff(mtfs.contrast[[x]],    mtfs.contrast[1:length(mtfs.contrast)!=x] %>% Reduce(c, .)  ) 
}) %>% givename(., names(mtfs.contrast))


#list.tfbs.shared=list(shared.tfbs=intersect(mtfs.contrast$highly.matched.peaks, mtfs.contrast$below.range), 
 # all.tfbs=union(mtfs.contrast$highly.matched.peaks, mtfs.contrast$below.range))



fcat("calculating gene set enrichments for each sensitivity group's TFs over all genesets...")
if(is.null(enrichments)){
  
simpleCache(paste0("enrichments_targets_target_",target_cell_type, "_", digest::digest(target.peaks),"_list_", digest::digest(list.tfbs.exclusive), "_gsetpars_", digest::digest(gsetpars)) %>% addversion, {   
enrichments=lapply(1:length(gsetpars), function(x) get.hyper.enrichments(list.tfbs.exclusive, gsetpars = gsetpars[x], make.barplot=F, make.heatmap = F))  %>% givename(., names(gsetpars))

enrichments},assignToVar="enrichments", reload=T) 
}

fcat("starting to prepare perturbation heatmaps...")

  perts=c("LIGPERTUP|LIGPERTLINCSUP|CHEMPERTUP|KINPERTUP", "LIGPERTDOWN|LIGPERTLINCSDN|CHEMPERTDOWN|KINPERTDOWN", "CGP")
  
  
perturbation.heatmaps=lapply(names(list.tfbs.exclusive), function(cluslab){

fcat("sensitivity group", cluslab, "...")

lapply(perts, function(searchterms){
  
fcat("Looking for groups enrichments in ", searchterms, "...")
endf=lapply(enrichments[grepl(searchterms, names(enrichments))],  function(x) x$df) %>% Reduce(rbind, .) 

#we make the enrichment heatmap matrix whereby we see what gene sets share the same hits for biological insight.
#we keep the matrix to make a more tailored heatmap later. 

fcat("preparing the enrichment matrix...")

eh=make.enrichment.heatmap(x=list.tfbs.exclusive[[cluslab]], gsetpars=gsetpars, endf %>% dplyr::filter(clusterlabel==!!cluslab))

if( !is.null(eh$matrix) ){
mat=remove.zerocols(remove.zerorows(eh$matrix %>% t))


gset.annot=endf %>% dplyr::filter(clusterlabel==!!cluslab) %>% filter(original.geneset %in% eh$original.gset.names) %>% select(jaccard_index, pval)


fcat("finding the relative relevance of each gene in the genesets for cell type ", target_cell_type, " ...")

relev=get.gene.gsea.relevance(colnames(mat), gsetpars, simple.output=F, keyword=keyword, label=target_cell_type, exclude=exclude)
tf.annot=relev$df %>% givecolnames(., nms=paste0("relevance.for.", target_cell_type)) 

#
fcat("ranking the perturbations to put those hitting the most relevant TFs on top")
#find the profiles that most closely involve target-relevant genes
hammings=lapply(1:nrow(mat), function(j){
sum(!(tf.annot[,1]>0 & as.logical(mat[j,])))
}) %>% unlist
horder=order(hammings, decreasing=F)
col_fun <- colorRampPalette(c(0, 1), c( "white", "green"))

fcat("preparing annotations for the perturbation heatmap...")
column_ha <- HeatmapAnnotation(
  annotation = anno_simple(tf.annot[,1]),
  show_annotation_name = TRUE
)


row_ha <- ComplexHeatmap::rowAnnotation( jaccard.index = anno_barplot(gset.annot))
fcat("executing the perturbation heatmap...")
hm.perturbations=ComplexHeatmap::pheatmap(mat[horder, ], col=c("purple", "yellow"), cluster_row=F, main=paste0("potential perturbations to ", ifelse(grepl("UP", searchterms), "activate", "deactivate")," TFs at ",cluslab , " peaks"), annotation_col=tf.annot, right_annotation=row_ha)
tpdf(path=config$plotpath, paste0("heatmap_tfbs_",  "_target_", target_cell_type, cluslab,"_perturbation_signatures", remove.pipes(searchterms)) %>% addversion, wi=pw*sca*3, he=pw*sca*1)
print(hm.perturbations)
dev.off()

}else{
  mat=NULL
  hm.perturbations=NULL
  }
  
  list(matrix=mat, heatmap=hm.perturbations)

}) %>% givename(.,perts %>% remove.pipes())  
  
}) %>% givename(., names(list.tfbs.exclusive))

list(
  range.mat=range.mat,
  tfbs.enrichment.list=tfbs.enrichment.list, 
  tfbs.per.region.all=mtfs.contrast,
  tfbs.per.region.exclusive=list.tfbs.exclusive,
  perturbation.heatmaps=perturbation.heatmaps, 
  sensitive.peaks=pks, 
  insensitive.peaks=peaks.always.below.range, 
  peak.range.count.df,
  venn.diagram=gven, 
  enrichments=enrichments, 
  
  )

}


################################################################################
# function to get enrichments of a list and store it as a cache
################################################################################

get.enrichments.cache=function(list.tfbs, gsetpars, cores=10){
    
  fcat("calculating gene set enrichments for each sensitivity group's TFs over all genesets...")

  cachename=paste0("enrichments_list_", digest::digest(list.tfbs), "_gsetpars_", digest::digest(gsetpars)) %>% addversion
simpleCache(cachename, { 
  
enrichments=parallel::mclapply(1:length(gsetpars), function(x, ll, gp){

simpleCache(paste0("enrichment_", digest::digest(ll), "_gset_", digest::digest(gp[x])) %>% addversion, {   
enr= get.hyper.enrichments(ll, gsetpars = gp[x], make.barplot=F, make.heatmap = F)
enr
}, assignToVar="enr", reload=T)
enr  
},ll=list.tfbs,gp=gsetpars, mc.cores=cores, mc.cleanup=T)

lapply(enrichments, function(x) x$df) %>% Reduce(rbind, .)

},assignToVar="enrichments", recreate=T) 

list(df=enrichments, gsets= names(gsetpars), cachename=cachename)  
}
################################################################################
# module to make perturbation heatmaps based on a tf list
################################################################################

################################################################################
# prepare deeptools script programatically by selecting sample groups from the workspace
# and automatically adding a reference sampe to cross compare
################################################################################


prepare.deeptools.script= function(
      test.class.variable,
  test.class,
    peak.list,
  peaks.id="unnamed-peaks",
  norm.peak.vals,
  meta,
  dataset.info,
  test.filter=NULL,
  invert=F,
  max.test.samples=NULL,
    reftype=config$target_cell_type,
  ref.sample.number=1,
  ref.sample.id=NULL,
  outdir=  NULL, 
  outdir.internal=NULL,
  outfile.matrix=NULL,
outfile.script= NULL, 
  outfile.sortedregions=NULL, 
  num.clusters=4, 
  append.to.file=F, 
  cluster.using.samples=1,
  sort.using.samples=1,
  heatmap.height=28, 
  bin.size=25, 
  specific.samples=NULL,
  yMin=0, 
  yMax=1, 
  zMin=0,
  zMax=2, 
  num.processors=1,
  selected.bigwigs=NULL, #provide the exact path of the bigwigs and ignore everything else
  date.use=NULL, 
  replace=F){
  
  fcat("preparing paths")
  test.class.formatted=gsub(" ", "-", test.class)
  
  newdir=timestamp(paste(c("deeptools_ref",reftype, test.class.variable, gsub("\\|","-", test.class.formatted), peaks.id ), collapse="_"), date.use=date.use)

  
  if(is.null(outdir.internal)){
   coredir=paste0(dataset.info$out_root, "/deeptoolsrerun/")
   
      if(!dir.exists(coredir)){
  dir.create(coredir)
  }
        outdir.internal=paste0(coredir, newdir)
  }
  
  outdir.external=outdir.internal %>% prep.dir.for.out.root
  
  if(!dir.exists(outdir.internal)){
  dir.create(outdir.internal)
  }
  
  if(is.null(outfile.matrix)){
    outfile.matrix=paste0(outdir.external, "/integrated_matrix.gz")
  }else{
  outifle.matrix=outfile.matrix %>% prep.dir.for.out.root
  }
  
    if(is.null(outfile.sortedregions)){
    outfile.sortedregions=paste0(outdir.external, "/sortedregions.bed")
    }else{
    outifle.sortedregions=outfile.sortedregions %>% prep.dir.for.out.root
    }
  
    if(is.null(outfile.script)){
    outfile.script="~/bash/deeptools_to_refsample.sh" 
    }
  
    outfile.script.archive=paste0(outdir.internal, "/run_deeptools.sh")
    outfile.scripts.list=list(bash=outfile.script, archive=outfile.script.archive)
    
  
  fcat("preparing peak region files")
  
  peakfilenames=lapply(1:length(peak.list), function(y){
  peak.list.name.formatted=gsub(" ", "-", names(peak.list)[y])  
  peakfilename.internal=paste0(outdir.internal,"/", peak.list.name.formatted, ".bed")
  peakfilename.external=paste0(outdir.external,"/", peak.list.name.formatted, ".bed")
  fcat(y, "extracting bed data")
  matt=norm.peak.vals[peak.list[[y]], 2:4]
  rownames(matt)=NULL
  fcat(y, "writing bed file")
  fwrite(matt, file=peakfilename.internal, row.names=F, col.names=F, sep="\t")
    fcat("making list")
  list(internal=peakfilename.internal, external=peakfilename.external, name= paste(names(peak.list)[y], collapse="_") )
}) %>% givename(., names(peak.list))


  allpeakfilenames.external= lapply(peakfilenames, function(x) x$external)  
  
  fcat("finding reference sample") # pending: merge test samples into one file, but this depends on having mergebigwig
  
  refbigwigpath=dataset.info$dataset_paths_nfcore[as.logical(dataset.info$is_reference)] %>% make.bigwigpath
  pprofilepath=dataset.info$dataset_paths_nfcore[as.logical(dataset.info$is_reference)] %>% make.profilepath 

  all.bigwig.paths.reference=list.files(path=refbigwigpath, pattern=".bigWig") %>% paste0(refbigwigpath,.) %>% prep.dir.for.out.root
  all.dmatrix.paths.reference=list.files(path=pprofilepath, pattern=".scale_regions.computeMatrix.mat.gz") %>% paste0(pprofilepath, .)

  if(is.null(ref.sample.id)){
  refsample=(meta %>% filter(celltype==!!reftype) %>% pull(Experiment))[ref.sample.number]
}else{
   refsample=meta %>% filter(celltype==!!reftype) %>%  grepl(ref.sample.id, Experiment) %>% pull(Experiment)
}
  
bigwig.refsample=grep(refsample, all.bigwig.paths.reference, value=T)
dmatrix.refsample=grep(refsample, all.dmatrix.paths.reference, value=T)
  
  fcat("finding all bigwig files and assempling their paths assembling their paths")
  
  
  if(is.null(selected.bigwigs)){
 # assemble the source paths where to find the bigwigs 
    sourcepaths.bigwig.root=dataset.info$dataset_paths_nfcore %>% make.bigwigpath 
 
  # look for the collection of bigwig files and make their paths
  all.bigwigpaths.parts=lapply(1:length(sourcepaths.bigwig.root), function(x){
  corefiles=list.files(path=sourcepaths.bigwig.root[x], pattern=".bigWig")
  
 list(filenames=corefiles, filepaths=corefiles %>% paste0(sourcepaths.bigwig.root[x], .) )
})
  
  
bigwigpaths.assembled= lapply(all.bigwigpaths.parts, function(x) x$filepaths) %>% Reduce(c, .) %>% prep.dir.for.out.root

#with all the possible paths assemble, now retrieve the relevant ones based on filters

if(is.null(specific.samples)){
selected.bigwigs= grep( (meta %>% filter(grepl(!!test.class, !!sym(test.class.variable))) %>% pull(Experiment) %>% paste(., collapse="|")), bigwigpaths.assembled, value=T )
}else{

  
selected.bigwigs= lapply(specific.samples, function(x){
  
  grep( (meta %>% filter(Experiment==x) %>% pull(Experiment) %>% paste(., collapse="|")), bigwigpaths.assembled, value=T )
  

}) %>% Reduce(c,.)

}


if(!is.null(test.filter)){
  
    selected.bigwigs= grep(test.filter, selected.bigwigs, invert=invert, value=T)
     
  }

  fcat("checking to adjust the total number of test samples displayed")
  if(!is.null(max.test.samples)){
    if(max.test.samples<=length(selected.bigwigs)){
   selected.bigwigs=selected.bigwigs[1:max.test.samples] 
    }
  }
}else{
  
 bigwig.refsample="" 
}
################################################################################
# WRITING SRCIPT  
################################################################################
  fcat("writing script")
  
    outfile.template="~/mnt_out/matrix_B.gz"
  
################################################################################
# prepare computeMAtrix command
################################################################################

computematrix.command=c(
  "computeMatrix reference-point -a 3000 -b 3000 --sortRegions descend -S ", bigwig.refsample, 
  selected.bigwigs,  
  " -bs ", bin.size,
    " --missingDataAsZero",
  " --referencePoint center",
  " -p", num.processors,
  " -R ",
  allpeakfilenames.external,"-o",
   outfile.matrix,
  " --skipZeros --outFileSortedRegions",
  outfile.sortedregions
  )  %>% paste(., collapse=" ")  

################################################################################
# prepare plotHeatmap command
################################################################################


heatmap.command= paste0(
  "plotHeatmap -m ", 
  outfile.matrix, 
  " --sortUsingSamples ",sort.using.samples,
  " --zMin ", zMin, " --zMax ", zMax, " ", 
  " --yMin ", yMin," --yMax ", yMax, " ",
    "  --outFileSortedRegions ", outdir.external, "/sorted_regions_heatmap.bed ",
  "  --outFileNameMatrix ", outdir.external, "/heatmap_matrix.gz ",
  " --heatmapHeight ", heatmap.height,
  " -o ", 
  prep.dir.for.out.root(dataset.info$plotpath), 
  paste_(newdir,"integrated_heatmap.pdf")
  )

heatmap.command.cluster= paste0(
  "plotHeatmap -m ", 
  outfile.matrix, 
  " --sortUsingSamples ",sort.using.samples,
  " --zMin ", zMin, " --zMax ", zMax, " ", 
  " --yMin ", yMin," --yMax ", yMax, " ",
    " --hclust ", num.clusters, " --clusterUsingSamples ", cluster.using.samples,
  "  --outFileSortedRegions ", outdir.external, "/sorted_regions_heatmap_clustered.bed ",
  "  --outFileNameMatrix ", outdir.external, "/heatmap_matrix_clustered.gz ",
  " --heatmapHeight ", heatmap.height,
  " -o ", 
  prep.dir.for.out.root(dataset.info$plotpath), 
  paste_(newdir,"integrated_heatmap_clustered.pdf")
  )



fcat("creating output scripts...")
#there is a working script and an archive script going to bash and a script going to the deeptools output so it can be permantently reproduced later on

#conditions to create the script
# if there is no date specified,  if there is a date specified but replace is true, 
if((is.null(date.use)) || (!is.null(date.use) && replace) ){
lapply(outfile.scripts.list, function(of){
if(!append.to.file){
  fwrite(list("#!/usr/bin/bash"), file=of, append=F)  
}

fwrite(list(computematrix.command), file=of, append=append.to.file)  
fwrite(list(heatmap.command), file=of, append=T)    
fwrite(list(heatmap.command.cluster), file=of, append=T)    
 
})
}

list(
  selected.bigwigs=selected.bigwigs, 
  bigwigpaths.assembled=bigwigpaths.assembled,
  computematrix.command=computematrix.command,
  heatmap.command=heatmap.command, 
  heatmap.command.cluster=heatmap.command.cluster,
  path.matrix=paste0(outdir.internal, "/integrated_matrix.gz"),
  path.heatmap.matrix=paste0(outdir.internal, "/heatmap_matrix.gz"),
  path.heatmap.matrix.clustered=paste0(outdir.internal, "/heatmap_matrix_clustered.gz"),
  path.sortedregions=paste0(outdir.internal, "/sorted_regions_heatmap.bed"), 
  path.sortedregions.clustered=paste0(outdir.internal, "/sorted_regions_heatmap_clustered.bed"), 
  selected.bigwigs=selected.bigwigs,
  all.bigwigs=c(bigwig.refsample, selected.bigwigs),
  outdir.internal=outdir.internal, 
  outdir.external=outdir.external,
  outfile.scripts.list=outfile.scripts.list
  )

}



collect.sample.bigwigs= function(meta, test.class, test.class.variable, bigwigvar="bigwigpath"){
  
  selected.bigwigs=  lapply(1:length(test.class), function(tc){
  fcat("collecting", test.class[tc], " samples")
 meta %>% dplyr::filter(grepl(test.class[tc], !!sym(test.class.variable))) %>% pull(!!sym(bigwigvar))
}) %>% givename(., test.class)
  
fcat(class(selected.bigwigs))
  
list(list=selected.bigwigs,   vector=selected.bigwigs %>% Reduce(c, .) )  

}

prepare.deeptools.script2= function(
      test.class.variable,
  test.class,
    peak.list,
  peaks.id="unnamed-peaks",
  meta,
  dataset.info,
  test.filter=NULL,
  invert=F,
  max.test.samples=NULL,
    reftype=config$target_cell_type,
  refclass=config$cell_type_variable,
  ref.sample.number=1,
  ref.sample.id=NULL,
  outdir=  NULL, 
  outdir.internal=NULL,
  outfile.matrix=NULL,
outfile.script= NULL, 
  outfile.sortedregions=NULL, 
  num.clusters=4, 
  append.to.file=F, 
  cluster.using.samples=1,
  sort.using.samples=1,
  heatmap.height=28, 
  bin.size=25, 
  specific.samples=NULL,
  yMin=0, 
  yMax=1, 
  zMin=0,
  zMax=2, 
  num.processors=1,
  selected.bigwigs=NULL, #provide the exact path of the bigwigs and ignore everything else
  date.use=NULL, 
  replace=F, 
  usemerged=T, 
  run.external=F, one.per.class=T){
  
  prepdir= ifelse(run.external, prep.dir.for.out.root, identity)
  
  
  fcat("preparing paths")
  test.class.formatted=gsub(" ", "-", test.class)
  
  newdir=timestamp(paste(c("deeptools_ref",reftype, test.class.variable, gsub("\\|","-", test.class.formatted), peaks.id ), collapse="_"), date.use=date.use)

  
  if(is.null(outdir.internal)){
   coredir=paste0(dataset.info$out_root, "/deeptoolsrerun/")
   
      if(!dir.exists(coredir)){
  dir.create(coredir)
  }
        outdir.internal=paste0(coredir, newdir)
  }
  
  outdir.external=outdir.internal %>% prepdir
  
  if(!dir.exists(outdir.internal)){
  dir.create(outdir.internal)
  }
  
  if(is.null(outfile.matrix)){
    outfile.matrix=paste0(outdir.external, "/integrated_matrix.gz")
  }else{
  outifle.matrix=outfile.matrix %>% prepdir
  }
  
    if(is.null(outfile.sortedregions)){
    outfile.sortedregions=paste0(outdir.external, "/sortedregions.bed")
    }else{
    outifle.sortedregions=outfile.sortedregions %>% prepdir
    }
  
    if(is.null(outfile.script)){
    outfile.script=paste0("~/bash/deeptools", ifelse(run.external, "external","internal"),".sh" )
    }
  
    outfile.script.archive=paste0(outdir.internal, "/run_deeptools", ifelse(run.external, "external","internal"),".sh")
    outfile.scripts.list=list(bash=outfile.script, archive=outfile.script.archive)
    
  
  fcat("preparing peak region files")
  
  peakfilenames=lapply(1:length(peak.list), function(y){
  peak.list.name.formatted=gsub(" ", "-", names(peak.list)[y])  
  peakfilename.internal=paste0(outdir.internal,"/", peak.list.name.formatted, ".bed")
  peakfilename.external=paste0(outdir.external,"/", peak.list.name.formatted, ".bed")
  fcat(y, "extracting bed data")
  matt=norm.peak.vals[peak.list[[y]], 2:4]
  rownames(matt)=NULL
  fcat(y, "writing bed file")
  fwrite(matt, file=peakfilename.internal, row.names=F, col.names=F, sep="\t")
    fcat("making list")
  list(internal=peakfilename.internal, external=peakfilename.external, name= paste(names(peak.list)[y], collapse="_") )
}) %>% givename(., names(peak.list))


  allpeakfilenames.external= lapply(peakfilenames, function(x) x$external)  
  
 
  
  ##############################################################################
  # preparing reference
  ##############################################################################
  fcat("preparing reference")
  meta.ref= meta %>% dplyr::filter(dscategory=="reference", !!sym(refclass)==reftype)
  if(is.null(ref.sample.id)){
  meta.ref=meta.ref[ref.sample.id, ]
  }else{
  meta.ref=meta.ref[ref.sample.number, , drop=F]  
  }
  
  ##############################################################################
  # preparing and ordering bigwigs
  ##############################################################################
  fcat("preparing bigwigs")
  #decide whether to use paths from merged library or merged replicate
  bigwigvar= ifelse(usemerged, "bigwigpath.merged", "bigwigpath")
  
#dmatrix.refsample=grep(refsample, all.dmatrix.paths.reference, value=T)
  
  fcat("finding all bigwig files and assempling their paths assembling their paths")
  
  
  if(is.null(selected.bigwigs)){

  

if(is.null(specific.samples)){
  
    fcat("no specific samples requested")
  
  #no specific samples have been defined: then collect paths for all the test.class elements in test.variable
selected.bigwigs0=collect.sample.bigwigs(meta, test.class, test.class.variable, bigwigvar) 



}else{

  #specific samples have been chosen
  fcat("specific samples have been chosen")
  
  selected.bigwigs0= collect.sample.bigwigs(meta, specific.samples, "Experiment", bigwigvar) 


}
    
    if(one.per.class){
     selected.bigwigs=lapply(selected.bigwigs0$list, function(x) x[1]) %>% Reduce(c,. ) 
    }else{
     selected.bigwigs=selected.bigwigs0$vector 
    }
    
################################################################################
# join with reference
################################################################################
selected.bigwigs= c(meta.ref %>% pull(!!sym(bigwigvar)), selected.bigwigs)
selected.bigwigs= unique(selected.bigwigs) %>% prepdir

fcat(selected.bigwigs)
################################################################################
# remove samples with a specific pattern in the file name
################################################################################
fcat("filtering out unwanted samples")
if(!is.null(test.filter)){
  
    selected.bigwigs= grep(test.filter, selected.bigwigs, invert=invert, value=T)
     
  }

  fcat("checking to adjust the total number of test samples displayed")
  if(!is.null(max.test.samples)){
    if(max.test.samples<=length(selected.bigwigs)){
   selected.bigwigs=selected.bigwigs[1:max.test.samples] 
    }
  }
}
  
 bigwig.refsample="" 

################################################################################
# WRITING SRCIPT  
################################################################################
  fcat("writing script")
  
    outfile.template="~/mnt_out/matrix_B.gz"
  
################################################################################
# prepare computeMAtrix command
################################################################################

computematrix.command=c(
  "computeMatrix reference-point -a 3000 -b 3000 --sortRegions descend -S ", 
  selected.bigwigs,  
  " -bs ", bin.size,
    " --missingDataAsZero",
  " --referencePoint center",
  " -p", num.processors,
  " -R ",
  allpeakfilenames.external,"-o",
   outfile.matrix,
  " --skipZeros --outFileSortedRegions",
  outfile.sortedregions
  )  %>% paste(., collapse=" ")  %>% prepdir

################################################################################
# prepare plotHeatmap command
################################################################################


heatmap.command= paste0(
  "plotHeatmap -m ", 
  outfile.matrix, 
  " --sortUsingSamples ",sort.using.samples,
  " --zMin ", zMin, " --zMax ", zMax, " ", 
  " --yMin ", yMin," --yMax ", yMax, " ",
    "  --outFileSortedRegions ", outdir.external %>% prepdir , "/sorted_regions_heatmap.bed ",
  "  --outFileNameMatrix ", outdir.external %>% prepdir, "/heatmap_matrix.gz ",
  " --heatmapHeight ", heatmap.height,
  " -o ", 
  dataset.info$plotpath, 
  paste_(newdir,"integrated_heatmap.pdf") 
  ) %>% prepdir

heatmap.command.cluster= paste0(
  "plotHeatmap -m ", 
  outfile.matrix, 
  " --sortUsingSamples ",sort.using.samples,
  " --zMin ", zMin, " --zMax ", zMax, " ", 
  " --yMin ", yMin," --yMax ", yMax, " ",
    " --hclust ", num.clusters, " --clusterUsingSamples ", cluster.using.samples,
  "  --outFileSortedRegions ", outdir.external %>% prepdir , "/sorted_regions_heatmap_clustered.bed ",
  "  --outFileNameMatrix ", outdir.external %>% prepdir, "/heatmap_matrix_clustered.gz ",
  " --heatmapHeight ", heatmap.height,
  " -o ", 
  dataset.info$plotpath %>% prepdir, 
  paste_(newdir,"integrated_heatmap_clustered.pdf")
  ) %>% prepdir



fcat("creating output scripts...")
#there is a working script and an archive script going to bash and a script going to the deeptools output so it can be permantently reproduced later on

#conditions to create the script
# if there is no date specified,  if there is a date specified but replace is true, 
if((is.null(date.use)) || (!is.null(date.use) && replace) ){
lapply(outfile.scripts.list, function(of){
  
  fcat("writing to", of)
if(!append.to.file){
  fwrite(list("#!/usr/bin/bash"), file=of, append=F)  
}

fwrite(list(computematrix.command), file=of, append=append.to.file)  
fwrite(list(heatmap.command), file=of, append=T)    
fwrite(list(heatmap.command.cluster), file=of, append=T)    
 
})
}

list(
  selected.bigwigs=selected.bigwigs, 
  original.bigwigs.list=selected.bigwigs0$list,
  computematrix.command=computematrix.command,
  heatmap.command=heatmap.command, 
  heatmap.command.cluster=heatmap.command.cluster,
  path.matrix=paste0(outdir.internal, "/integrated_matrix.gz"),
  path.heatmap.matrix=paste0(outdir.internal, "/heatmap_matrix.gz"),
  path.heatmap.matrix.clustered=paste0(outdir.internal, "/heatmap_matrix_clustered.gz"),
  path.sortedregions=paste0(outdir.internal, "/sorted_regions_heatmap.bed"), 
  path.sortedregions.clustered=paste0(outdir.internal, "/sorted_regions_heatmap_clustered.bed"), 
  selected.bigwigs=selected.bigwigs,
  all.bigwigs= selected.bigwigs,
  outdir.internal=outdir.internal, 
  outdir.external=outdir.external,
  outfile.scripts.list=outfile.scripts.list
  )

}



################################################################################
#prepare bigwigmerge command programatically
################################################################################

prepare.bigwigmerge= function(
      test.class.variable,
  test.class,
    peak.list,
  peaks.id="unnamed-peaks",
  norm.peak.vals,
  meta,
  dataset.info,
  test.filter=NULL,
  invert=F,
  max.test.samples=NULL,
    reftype=config$target_cell_type,
  ref.sample.number=1,
  ref.sample.id=NULL,
  outdir=  NULL, 
  outdir.internal=NULL,
  outfile.matrix=NULL,
outfile.script= NULL, 
  outfile.sortedregions=NULL, 
  num.clusters=4, 
  append.to.file=F, 
  cluster.using.samples=1,
  sort.using.samples=1,
  heatmap.height=28, 
  bin.size=25, 
  specific.samples=NULL,
  yMin=0, 
  yMax=1, 
  zMin=0,
  zMax=2
  ){
  
  fcat("preparing paths")
  
  test.class.formatted=gsub(" ", "-", test.class)
  newdir=timestamp(paste(c("bigwigmerge", test.class.variable, gsub("\\|","-", test.class.formatted), peaks.id ), collapse="_"))

  
  if(is.null(outdir.internal)){
   coredir=paste0(dataset.info$out_root, "/bigwigmerge/")
   
      if(!dir.exists(coredir)){
  dir.create(coredir)
  }
        outdir.internal=paste0(coredir, newdir)
  }
  
  outdir.external=outdir.internal %>% prep.dir.for.out.root
  
  if(!dir.exists(outdir.internal)){
  dir.create(outdir.internal)
  }
  
  if(is.null(outfile.matrix)){
    outfile.matrix=paste0(outdir.external, "/integrated_matrix.gz")
  }else{
  outifle.matrix=outfile.matrix %>% prep.dir.for.out.root
  }
  
    if(is.null(outfile.sortedregions)){
    outfile.sortedregions=paste0(outdir.external, "/sortedregions.bed")
    }else{
    outifle.sortedregions=outfile.sortedregions %>% prep.dir.for.out.root
    }
  
    if(is.null(outfile.script)){
    outfile.script="~/bash/mergebigwig.sh"
    }
  
  chr.sizes.file=prep.dir.for.project.root("~/metadata/GRCh38_EBV.chrom.sizes.tsv")
  

  fcat("finding reference sample") # pending: merge test samples into one file, but this depends on having mergebigwig
  
  refbigwigpath=dataset.info$dataset_paths_nfcore[as.logical(dataset.info$is_reference)] %>% make.bigwigpath
  pprofilepath=dataset.info$dataset_paths_nfcore[as.logical(dataset.info$is_reference)] %>% make.profilepath 

  all.bigwig.paths.reference=list.files(path=refbigwigpath, pattern=".bigWig") %>% paste0(refbigwigpath,.) %>% prep.dir.for.out.root
  all.dmatrix.paths.reference=list.files(path=pprofilepath, pattern=".scale_regions.computeMatrix.mat.gz") %>% paste0(pprofilepath, .)

  if(is.null(ref.sample.id)){
  refsample=(meta %>% filter(celltype==!!reftype) %>% pull(Experiment))[ref.sample.number]
}else{
   refsample=meta %>% filter(celltype==!!reftype) %>%  grepl(ref.sample.id, Experiment) %>% pull(Experiment)
}
  
bigwig.refsample=grep(refsample, all.bigwig.paths.reference, value=T)
dmatrix.refsample=grep(refsample, all.dmatrix.paths.reference, value=T)
  
  fcat("finding all bigwig files and assempling their paths assembling their paths")
  
 # assemble the source paths where to find the bigwigs 
    sourcepaths.bigwig.root=dataset.info$dataset_paths_nfcore %>% make.bigwigpath 
 
  # look for the collection of bigwig files and make their paths
  all.bigwigpaths.parts=lapply(1:length(sourcepaths.bigwig.root), function(x){
  corefiles=list.files(path=sourcepaths.bigwig.root[x], pattern=".bigWig")
  
 list(filenames=corefiles, filepaths=corefiles %>% paste0(sourcepaths.bigwig.root[x], .) )
})
  
  
bigwigpaths.assembled= lapply(all.bigwigpaths.parts, function(x) x$filepaths) %>% Reduce(c, .) %>% prep.dir.for.out.root

#with all the possible paths assemble, now retrieve the relevant ones based on filters

if(is.null(specific.samples)){
selected.bigwigs= grep( (meta %>% filter(grepl(!!test.class, !!sym(test.class.variable))) %>% pull(Experiment) %>% paste(., collapse="|")), bigwigpaths.assembled, value=T )
}else{

  
selected.bigwigs= lapply(specific.samples, function(x){
  
  grep( (meta %>% filter(Experiment==x) %>% pull(Experiment) %>% paste(., collapse="|")), bigwigpaths.assembled, value=T )
  

}) %>% Reduce(c,.)

}


if(!is.null(test.filter)){
  
    selected.bigwigs= grep(test.filter, selected.bigwigs, invert=invert, value=T)
     
  }

  fcat("checking to adjust the total number of test samples displayed")
  if(!is.null(max.test.samples)){
    if(max.test.samples<=length(selected.bigwigs)){
   selected.bigwigs=selected.bigwigs[1:max.test.samples] 
    }
  }
  
  
  heatmap.command= paste0(
  "plotHeatmap -m ", 
  outfile.matrix, 
  " --sortUsingSamples ",sort.using.samples,
  " --zMin ", zMin, " --zMax ", zMax, " ", 
  " --yMin ", yMin," --yMax ", yMax, " ",
    "  --outFileSortedRegions ", outdir.external, "/merged_bigwigs_bedgraph.bg ",
  "  --outFileNameMatrix ", outdir.external, "/heatmap_matrix_clustered.gz ",
  " --heatmapHeight ", heatmap.height,
  " -o ", 
  prep.dir.for.out.root(dataset.info$plotpath), 
  paste_(newdir,"integrated_heatmap.pdf")
  )
  
  ##############################################################################
  # prepare escript
  ##############################################################################
  
  out.bedgraph.nonsorted=paste0(outdir.external, "/merged_bigwigs_bedgraph.bedGraph ")
  out.bedgraph.sorted=paste0(outdir.external, "/merged_bigwigs_bedgraph_sorted.bedGraph ")
  
  mbw.command= paste0("bigWigMerge ",
  paste(selected.bigwigs, collapse=" "),  " ", out.bedgraph.nonsorted )
  
sort.command=paste0("sort -k1,1 -k2,2n ", out.bedgraph.nonsorted, ">", out.bedgraph.sorted)
    
 final.bigwig=paste0(outdir.external, "/merged_bigwigs_", test.class.formatted, ".bw" )
btbb.command=paste0("bedGraphToBigWig ", out.bedgraph.sorted ," ", chr.sizes.file, " ", final.bigwig)




#fwrite(list( "#!/usr/bin/bash"), file=outfile.script, append=append.to.file)
fwrite(list( mbw.command), file=outfile.script, append=append.to.file)
fwrite(list( sort.command), file=outfile.script, append=T, quote=F)
fwrite(list( btbb.command), file=outfile.script, append=T, quote=F)


list(merge.bigwig.command=mbw.command, bedgraphtobigwig.command=btbb.command, sort.command=sort.command, final.bigwig=final.bigwig)



}



###original


de.matches.invitro.workflow.individual = function(
    envlist, recreate = F, ref = config$target_cell_type,
    lineage = NULL, refvar = "celltype", reload = T, cores = 10,
    comparison.id = paste0("vs", refvar, "_", ref), 
    cutoff.mean=0, min.reads=0, min.expts=0, altHypothesis=NULL, filter.normalized=F
){
  library(stringr)

  dds = envlist$dds
  meta2 = envlist$meta2
  dataset.info = envlist$dataset.info
  peaks = envlist$peaks
  meta3 = meta2

  fcat("entering loop")

  lapply(1:length(test.ids.list), function(k)  {

    refs = paste(c(ref, lineage), collapse = "|")

    allinvitro = colData(dds) %>%
      as.data.frame %>%
      dplyr::filter((Experiment %in% test.ids.list[[k]]) |
                      grepl(refs, !!sym(refvar))) %>%
      arrange(round) %>% pull(Experiment)

    positives = make.class.target.list(
      meta3[allinvitro, , drop = FALSE], "Condition"
    )

    # -------------------------
    # ORIGINAL SUBSET OF PEAKS
    # -------------------------
   dds.invitro= dds[Reduce(c, peaks) %>% unique ,allinvitro] #case 2: only include concordance peaks in the DE calculation

    
    
    rng = 1:length(positives)

    par.deseq2.invitro.eachvstarget = parallel::mclapply(
      rng,
      function(x, dd, allpos, pd) {

        lab = names(allpos)[x]
        if (lab == ref) {
          return(NULL)
        } else {
          classvar = "Condition"
          levs = c(lab, ref)

          is.target = function(x, targets)
            fifelse(x %in% targets, levs[1], levs[2])

          simpleCache::simpleCache(
            paste0(
              "deseq2_invitro_individual-", lab,
              "_vs_", ref,
              "_id_", comparison.id,
              "filtermean", cutoff.mean, 
              "minreads", min.reads, 
              "minexpts", min.expts
              #"_positivesid", digest::digest(allpos),
              #"_ddsid_", digest::digest(dd)
            ) %>% addversion,
            {
              dd[[paste0("target_vs_", ref)]] <-
                factor(is.target(dd[[classvar]], allpos[[x]]$targets),
                       levels = levs)

              design(dd) <- as.formula(paste0("~target_vs_", ref))
              
              if(cutoff.mean>0){
                 # -----------------------------------------------------------
              # NEW STEP: Remove peaks whose normalized baseMean <= 2
              # -----------------------------------------------------------
              dd <- estimateSizeFactors(dd)
              norm.mat <- counts(dd, normalized = TRUE)
              baseMean <- rowMeans(norm.mat)

              keep <- baseMean >= cutoff.mean   # Remove baseMean 
              dd <- dd[keep, ]
              rm(norm.mat)
              }
              
              if(any(c(min.expts, min.reads)>0)){
               
                matt <- counts(dd, normalized = filter.normalized)
                
                
                keep_by_condition <- sapply(levs, function(cond) {
                
                # Count how many samples per peak have >= y reads
                 rowSums(matt[, get.experiments(cond, meta3), drop = FALSE] >= min.reads) >= min.expts
                })

                # keep peaks that satisfy filtering in *all* conditions
                keep <- apply(keep_by_condition, 1, any) 
                dd <- dd[keep, ]
                 
              }

              dd <- DESeq(dd)

              res <- results(dd,
                             contrast = c(paste0("target_vs_", ref),
                                          levs[1], levs[2]), altHypothesis=altHypothesis)
            },
            assignToVar = "res",
            reload = reload, 
            recreate=!reload
          )

          return(res %>% as.data.frame %>% 
                   dplyr::mutate(target=lab, id = !!comparison.id) %>% names2col(., "Geneid"))
        }
      },
      allpos = positives,
      dd = dds.invitro,
      mc.cores = cores,
      mc.cleanup = TRUE
    ) %>% givename(., names(positives[rng]))

  }) %>% givename(., paste0("dataset", 1:length(test.ids.list)))
}








de.matches.invitro.workflow.sepbatches=function(envlist, recreate=F, ref=config$target_cell_type, refvar="celltype", reload=T, cores=1, comparison.id=paste0("vs", refvar, "_", ref) ){

  dds=envlist$dds
  meta2=envlist$meta2
  dataset.info=envlist$dataset.info
  selected.peaks=envlist$selected.peaks
  
  #meta3=conditional.suffix(meta2, search.var="dsname", target.var="Condition", pattern="Seruggia2025-3", suf="v")
  meta3=meta2
  colData(dds)= DataFrame(meta3)
  
  
## we can compare to all peaks or to selected peaks (concordance peaks)
  comp.types=c("all peaks", "selected peaks")
  
# go through each comparison type and eah batch  
alldes=lapply(comp.types, function(comp.type){  
 lapply(1:length(test.ids.list), function(k)  {
  
 #isolate the experiments of this one  batch together with the reference samples
allinvitro= colData(dds) %>% as.data.frame %>% dplyr::filter((Experiment %in%  test.ids.list[[k]])|grepl(ref,!!sym(refvar)) ) %>% arrange(round) %>% pull(Experiment)

# for the selected experiments, make a list of all existing conditions we want to loop through
positives=make.class.target.list(meta3[allinvitro, , drop=F], "Condition")


if(comp.type=="all peaks"){
dds.invitro= dds[,allinvitro] #case 1: include all peaks
}else{
dds.invitro= dds[selected.peaks ,allinvitro] #case 2: only include selected peaks in the DE calculation
}


rng=1:(length(positives))
par.deseq2.invitro.eachvstarget=parallel::mclapply(rng, function(x, dd, allpos, pd, ct){
  
  lab=names(allpos)[x]
  if(lab==ref){
  return(NULL)
  }else{
  classvar="Condition"
  levs=c(lab, ref)
  is.target= function(x, targets) fifelse(x %in% targets, levs[1], levs[2])
  
  
simpleCache::simpleCache(paste0("deseq2_invitro_batches_id_", comparison.id, "_",lab, "_vs_",ref, "_comparedpeaks", comp.type) %>% addversion, {
  

dd[[paste0("target_vs_",ref)]] <- factor(is.target(dd[[classvar]], allpos[[x]]$targets), levels=levs)


################################################################################
# custom filtering
################################################################################


# Re-run the DESeq pipeline with the new condition
design(dd) <- as.formula(paste0("~target_vs_", ref))
dd <- DESeq(dd)  
res <- results(dd, contrast = c(paste0("target_vs_", ref), levs[1], levs[2]))  

 }, assignToVar="res", recreate=T)

  return(res %>% as.data.frame %>% dplyr::mutate(id=!!comparison.id))
  }
  },allpos=positives,dd=dds.invitro,ct=comp.type, mc.cores=cores, mc.cleanup=T) %>% givename(., names(positives[rng]))

}) %>% givename(., paste0("dataset", 1:length(test.ids.list)) )
  }) %>% givename(., comp.types)
alldes  
}


de.matches.invitro.workflow.sepbatches2=function(envlist, recreate=F, ref=config$target_cell_type, refvar="celltype", reload=T, cores=1,comparison.id=paste0("vs", refvar, "_", ref) , selected.conditions="all"){

  dds=envlist$dds
  meta2=envlist$meta2
  dataset.info=envlist$dataset.info
  selected.peaks=envlist$selected.peaks
  
  #meta3=conditional.suffix(meta2, search.var="dsname", target.var="Condition", pattern="Seruggia2025-3", suf="v")
  meta3=meta2
  colData(dds)= DataFrame(meta3)
  
  
## we can compare to all peaks or to selected peaks (concordance peaks)
  #comp.types=c("all peaks", "selected peaks")
  
# go through each comparison type   
#alldes=lapply(comp.types, function(comp.type){  
 # go through each comparison type and each batch 
 
    if(selected.conditions!="all"){
      test.ids.list2= lapply(test.ids.list, function(x) x[x %in% selected.conditions])
    }else{
     test.ids.list2=test.ids.list 
    }
      
  
  alldes=lapply(1:length(test.ids.list2), function(k){
    
    if(!all(is.null(test.ids.list2[[k]]))){

 #isolate the experiments of this one  batch together with the reference samples
allinvitro= colData(dds) %>% as.data.frame %>% dplyr::filter((Experiment %in%  test.ids.list2[[k]])|grepl(ref,!!sym(refvar)) ) %>% arrange(round) %>% pull(Experiment)

# for the selected experiments, make a list of all existing conditions we want to loop through
positives=make.class.target.list(meta3[allinvitro, , drop=F], "Condition")


if(!is.null(selected.peaks)){

dds.invitro= dds[selected.peaks ,allinvitro] #case 2: only include selected peaks in the DE calculation
}

rng=1:(length(positives))

par.deseq2.invitro.eachvstarget=parallel::mclapply(rng, function(x, dd, allpos, rld,ct){
  
  
simpleCache::simpleCache(paste0("deseq2_invitro_sepbatches2_id_", comparison.id, "_vs_",ref, "_comparedpeaks", ct) %>% addversion, {

 dd[[refvar]] <- factor(dd[[refvar]]) 
  ## in this version the design is literally per condition in the batch
design(dd) <- as.formula(paste0("~", refvar))
dd <- DESeq(dd)  
dd
 }, assignToVar="dds.invitro", reload=rld)
  

  
  lab=names(allpos)[x]
  if(lab==ref){
  return(NULL)
  }else{
  classvar="Condition"
  levs=c(lab, ref)
  is.target= function(x, targets) fifelse(x %in% targets, levs[1], levs[2])

res <- results(dd, contrast = classvar, levs[1], levs[2])  

  return(res %>% as.data.frame %>% dplyr::mutate(id=!!comparison.id, target=lab, reference=ref))
  
  }},allpos=positives,dd=dds.invitro,rld, reload, ct=comp.type, mc.cores=cores, mc.cleanup=T) %>% givename(., names(positives[rng]))

par.deseq2.invitro.eachvstarget
   
}else{return(NULL)}


  }) %>% givename(., paste0("dataset", 1:length(test.ids.list2)) )
 

    
alldes  
}







de.matches.invitro.workflow.sepbatches2=function(envlist, recreate=F, ref=config$target_cell_type, refvar="celltype",
                                                 reload=T, cores=1, comparison.id=paste0("vs", refvar, "_", ref),
                                                 selected.conditions="all", test.ids.list){

  fcat("test: beginning\n")

  dds=envlist$dds
  meta2=envlist$meta2
  dataset.info=envlist$dataset.info
  selected.peaks=envlist$selected.peaks

  fcat("test: envlist loaded\n")

  meta3=meta2
  colData(dds)= DataFrame(meta3)

  fcat("test: colData assigned\n")

  # selected conditions filtering
  if(selected.conditions!="all"){
    fcat("test: filtering selected conditions\n")
    test.ids.list2= lapply(test.ids.list, function(x) x[x %in% selected.conditions])
  } else {
    test.ids.list2=test.ids.list 
  }

  fcat("test: test.ids.list2 prepared, length = ", length(test.ids.list2), "\n")

  alldes=lapply(1:length(test.ids.list2), function(k){

    fcat("test: starting dataset index k = ", k, "\n")



      # isolate invitro experiments
      allinvitro= meta3 %>% as.data.frame %>%
        dplyr::filter((Experiment %in%  test.ids.list2[[k]]) | grepl(ref, !!sym(refvar))) %>%
        arrange(round) %>% pull(Experiment)

      fcat("test: allinvitro computed", paste(allinvitro, collapse=","), "\n")
      
          if(length(allinvitro)>1){

      fcat("test: entering batch loop for k = ", k, "\n")
      
      

      # generate positives
      positives=make.class.target.list(meta3[allinvitro, , drop=F], "Condition")

      fcat("test: positives generated, length = ", length(positives), "\n")

      # subset peaks if provided
      if(!is.null(selected.peaks)){
        fcat("test: subsetting selected.peaks\n")
        dds.invitro= dds[selected.peaks ,allinvitro]
        fcat("test: dds.invitro subset done\n")
      }

      rng=1:(length(positives))
      fcat("test: rng defined, size = ", length(rng), "\n")

      par.deseq2.invitro.eachvstarget=parallel::mclapply(
        rng,
        function(x, dd, allpos, rld, ct){
          fcat("test: mclapply iteration x = ", x, "\n")

          simpleCache::simpleCache(
            paste0("deseq2_invitro_sepbatches2_id_", comparison.id, "_vs_",ref, "_comparedpeaks", ifelse(is.null(ct), "all", digest::digest(ct))) %>% addversion,
            {
              fcat("test: inside simpleCache block, x = ", x, "\n")

              dd[[refvar]] <- factor(dd[[refvar]])
              design(dd) <- as.formula(paste0("~", refvar))
              fcat("test: DESeq design set\n")

              dd <- DESeq(dd)
              fcat("test: DESeq completed for x = ", x, "\n")

              dd
            },
            assignToVar="dds.invitro", reload=rld
          )

          lab=names(allpos)[x]
          fcat("test: lab = ", lab, "\n")

          if(lab==ref){
            fcat("test: lab equals ref, returning NULL\n")
            return(NULL)
          } else {
            classvar="Condition"
            levs=c(lab, ref)

            res <- results(dd, contrast=classvar, levs[1], levs[2])
            fcat("test: results extracted for lab = ", lab, "\n")

            return(res %>% as.data.frame %>% dplyr::mutate(id=!!comparison.id, target=lab, reference=ref))
          }

        },
        allpos=positives,
        dd=dds.invitro,
        rld=reload,
        ct=selected.peaks,
        mc.cores=cores,
        mc.cleanup=T
      ) %>% givename(., names(positives[rng]))

      fcat("test: mclapply finished for k = ", k, "\n")

      par.deseq2.invitro.eachvstarget

    } else {
      fcat("test: NULL test.ids.list2[[k]]\n")
      return(NULL)
    }

  }) %>% givename(., paste0("dataset", 1:length(test.ids.list2)))

  fcat("test: alldes finished\n")

  alldes
}













de.matches.invitro.workflow.lineage=function(envlist, recreate=F, ref=config$target_cell_type, lineage=c("CMP","MEP"), refvar="celltype", reload=T, cores=10, comparison.id=paste0("vs", refvar, "_", ref)){

#this function includes  several relevant samples from the target lineage, those whose differentiall gene expression may be relevant to quantify variation in the discordant regions.   
  
  dds=envlist$dds
  meta2=envlist$meta2
  dataset.info=envlist$dataset.info
  
  meta3=conditional.suffix(meta2, search.var="dsname", target.var="Condition", pattern="Seruggia2025-3", suf="v")
  
  colData(dds)= DataFrame(meta3)
  
 lapply(1:length(test.ids.list), function(k)  {
  
  ### perform the DE , test dataset versus dataset 
refs= paste(c(ref, lineage), collapse="|")   

allinvitro= colData(dds) %>% as.data.frame %>% dplyr::filter((Experiment %in%  test.ids.list[[k]])|grepl(refs,!!sym(refvar)) ) %>% arrange(round) %>% pull(Experiment)


positives=make.class.target.list(meta3[allinvitro, , drop=F], "Condition")


#dds.invitro= dds[,allinvitro] #case 1: include all peaks
dds.invitro= dds[Reduce(c, concordance.peaks) %>% unique ,allinvitro] #case 2: only include concordance peaks in the DE calculation
rng=1:(length(positives))
par.deseq2.invitro.eachvstarget=parallel::mclapply(rng, function(x, dd, allpos, pd){
  
  lab=names(allpos)[x]
  if(lab==ref){
    
  return(NULL)
  }else{
  classvar="Condition"
  levs=c(lab, ref)
  is.target= function(x, targets) fifelse(x %in% targets, levs[1], levs[2])
  
  
simpleCache::simpleCache(paste0("deseq2_invitro_validate_target-",lab, "_vs_",ref, "_positivesid", digest::digest(allpos),"_ddsid_",digest::digest(dd),  "analysis", analysis.version), {
  

dd[[paste0("target_vs_",ref)]] <- factor(is.target(dd[[classvar]], allpos[[x]]$targets), levels=levs)

# Re-run the DESeq pipeline with the new condition
design(dd) <- as.formula(paste0("~target_vs_", ref))
dd <- DESeq(dd)  
res <- results(dd, contrast = c(paste0("target_vs_", ref), levs[1], levs[2]))  

 }, assignToVar="res", reload=reload)

  return(res %>% as.data.frame %>% dplyr::mutate(id=!!comparison.id))
  }
  },allpos=positives,dd=dds.invitro, mc.cores=cores, mc.cleanup=T) %>% givename(., names(positives[rng]))

}) %>% givename(., paste0("dataset", 1:length(test.ids.list)) )
}




targeted.changes.analysis=function(delist, concordance.peaks, exclude=NULL, pval=0.05, thr=0.2, filter=F, facet="treatment"){
  null.elements=lapply(delist, is.null) %>% unlist
  fcat("removing the following null elements: ", paste(names(null.elements[null.elements]), collapse=","))
  
delist=delist[ !null.elements %>% unname]

  
  
check.concordance=function(x){
ifelse(x %in% concordance.peaks.unique$concordant.open, "concordant.open", 
  ifelse(x %in% concordance.peaks.unique$discordant, "discordant", 
    ifelse(x %in% concordance.peaks.unique$concordant.closed, "concordant.closed", "other")))
}
  
  
demat.concord<-lapply(1:length(delist), function(x)  delist[[x]] %>% as.data.frame %>% dplyr::filter(abs(log2FoldChange)>=thr, pvalue<=pval) %>% names2col(., "Geneid") %>% dplyr::mutate(treatment=names(delist)[x], effect=ifelse(log2FoldChange>=1, "differentially open", "differentially closed"), is.discordant= Geneid %in% concordance.peaks$discordant, is.concordant.open=Geneid %in% concordance.peaks$concordant.open, is.concordant.closed=Geneid %in% concordance.peaks$concordant.closed) %>% arrange(log2FoldChange) %>% dplyr::mutate(concordance=check.concordance(Geneid))) %>% Reduce(rbind, .)


demat.depeaks= demat.concord %>% group_by(treatment, concordance) %>% summarise(counts=n())


if(!is.null(exclude)){
  demat.concord=demat.concord %>% dplyr::filter(!(treatment %in% exclude))
}



demat.summary= demat.concord %>% group_by(treatment, effect) %>% summarise(percent.discordant=sum(is.discordant)/length(concordance.peaks$discordant), percent.concordant.open=sum(is.concordant.open)/length(concordance.peaks$concordant.open), percent.concordant.closed=sum(is.concordant.closed)/length(concordance.peaks$concordant.closed))


demat.tidy=pivot_longer(demat.summary, contains("percent"), names_to="type", values_to="percentage") %>% as.data.frame %>% dplyr::mutate(percentage.=ifelse(effect=="differentially closed", -percentage, percentage), type=gsub("percent.", "", type))

plt=ggplot(demat.tidy, aes(y=factor( type, levels=c("discordant", "concordant.closed", "concordant.open") %>% rev), x=percentage., fill=effect))+geom_col()+facet_wrap(~treatment, ncol=1, scale="free_y")+scale_fill_manual(values=c("differentially closed"="red", "differentially open"="cyan"))+theme_classic()

plt=ggplot(demat.tidy, aes(y=treatment, x=percentage., fill=effect))+geom_col()+facet_wrap(as.formula(~factor(type, levels=c("concordant.open", "discordant", "concordant.closed") )), ncol=1, scale="free_y")+scale_fill_manual(values=c("differentially closed"="red", "differentially open"="cyan"))+theme_classic()

plt2=ggplot(demat.concord %>% dplyr::filter(!is.na(concordance)), aes( x=log2FoldChange, fill=concordance))+geom_density(alpha=0.5)+facet_wrap(as.formula(~factor(treatment)), ncol=1, scale="free_y")
#+scale_fill_manual(values=c("differentially closed"="red", "differentially open"="cyan"))+theme_classic()

list(df.tidy=demat.tidy, df=demat.concord, plot.percentages=plt, plot.effectsizes=plt2, df.depeaks=demat.depeaks)
}






umap.louvain.clustering <- function(umap_result, resolution=1) {
  
  # Step 1: Extract the nearest neighbors matrix from the UMAP result
  nn_matrix <- umap_result$nn$euclidean$idx
  
  # Step 2: Create an edge list from the nearest neighbors matrix
  # Create pairs of observations based on their nearest neighbors
  edges <- do.call(rbind, lapply(1:nrow(nn_matrix), function(i) {
    cbind(i, nn_matrix[i, ])
  }))
  
  # Step 3: Create a graph using the edge list
  g <- graph_from_edgelist(edges, directed = FALSE)
  
  # Step 4: Apply the Louvain clustering algorithm to the graph
  louvain_result <- cluster_louvain(g, resolution=resolution)
  
  # Step 5: Extract cluster membership
  cluster_membership <- membership(louvain_result)
  
  # Return the cluster membership vector
  return(cluster_membership)
}


################################################################################
# functions to process information from the deeptools output
################################################################################
get.sample.coords= function(k, extend.from.center=NULL, cols.per.sample){
 startpoint=cols.per.sample*(k-1)+1
 endpoint=cols.per.sample*k
  
 if(is.null(extend.from.center)){
   ctr=round((endpoint-startpoint)/2)
   
   (ctr-extend.from.center):(ctr+extend.from.center)
   
 }else{
 startpoint:endpoint
 }
}


deeptools.sample.umap=function(smpl, extend.from.center=30, md=.3, nn=20, ncomp=60, rs=1, matth=matth, matts=matts){ 
  submat=matth[,get.sample.coords(1, extend.from.center=extend.from.center)]
mattpca=prcomp(submat, scale=T)
umap.deeptools.object <- umap(mattpca$x[, 1:ncomp], init = "random",  min_dist = md,n_neighbors = nn, repulsion_strength=rs, ret_model=T)
umap.df= umap.deeptools.object$embedding %>% as.data.frame %>% givecolnames(., ind=c(1,2), nms=c("UMAP1", "UMAP2")) %>% cbind(.,matts) %>% dplyr::mutate(value.at.peak.center=submat[, round(extend.from.center/2)])
(ggplot(umap.df, aes(x=UMAP1, y=UMAP2, color=deepTools_group))+geom_point())+ggplot(umap.df, aes(x=UMAP1, y=UMAP2, color=value.at.peak.center))+geom_point()+scale_color_viridis(limits = c(0, 3), oob = scales::squish)
}


deeptools.cross.samples.umap=function(samples, extend.from.center=30, md=.3, nn=20, ncomp=60, rs=1,coll="median", cluster.resolution=1, matth=matth, matts=matts, cols.per.sample){ 
  submat.list=lapply(1:length(samples), function(xx){
  matth[,get.sample.coords(samples[xx], extend.from.center=extend.from.center, cols.per.sample=cols.per.sample)]
  }) 
  

  submat=submat.list %>% Reduce(cbind, .)
  
  submat.medians= apply(submat, 1, median)
  
  submat.centeredmedians=lapply(submat.list, function(mm) mm[, round(extend.from.center/2), drop=T]) %>% Reduce(cbind, .) %>% apply(., 1, median)
      
      peak.center.medians=lapply(1:length(samples), function(xx){
  matth[,get.sample.coords(samples[xx], extend.from.center=extend.from.center, cols.per.sample=cols.per.sample)] 
  }) %>% Reduce(cbind, .)
  
mattpca=prcomp(submat, scale=T)
umap.deeptools.object <- umap(mattpca$x[, 1:ncomp], init = "random",  min_dist = md,n_neighbors = nn, repulsion_strength=rs, ret_model=T, ret_nn=T)

louvain.result=umap.louvain.clustering(umap.deeptools.object, resolution=cluster.resolution)

hclust.result=hclust(dist(mattpca$x[, 1:ncomp]))$labels


umap.df= umap.deeptools.object$embedding %>% as.data.frame %>% givecolnames(., ind=c(1,2), nms=c("UMAP1", "UMAP2")) %>% cbind(.,matts) %>% dplyr::mutate( median= submat.medians, median.centers=submat.centeredmedians, louvain.cluster=factor(as.vector(louvain.result)), hcluster=hclust.result) 

clusternames=unique(umap.df$louvain.cluster)

peak.list=lapply(clusternames, function(ll){

  umap.df %>% filter(louvain.cluster==!!ll) %>% pull(Geneid)  
  
}) %>% givename(., clusternames)

plt=(ggplot(umap.df, aes(x=UMAP1, y=UMAP2, color=deepTools_group))+geom_point())+ggplot(umap.df, aes(x=UMAP1, y=UMAP2, color=!!sym(coll)))+geom_point()

list(umap.df=umap.df, plot=plt, peak.list=peak.list, submat=submat
  )

}


add.peak.annotation.to.mat=function(matt, start.ind=1){
  
  lapply(1:nrow(matt), function(y){
  
  annotation.peaks %>% dplyr::filter(Chr==matt[y, start.ind], Start==matt[y,start.ind+1], End==!!matt[y, start.ind+2])
  
  
}) %>% Reduce(rbind, .) %>% bind_cols(matt,.)

}


peak.clustering.workflow=function(matfile, matfile.heatmap, matfile.sortedregions, nsamples,which.samples.umap, umap.pars, annotation.peaks){

  
list2env(umap.pars)  
library(uwot)
library(igraph)
#library(leiden)
# import the matrix
fcat("importing matrices")
matfile=matfile %>% prep.dir.for.in.root
matfile.heatmap= matfile.heatmap %>% prep.dir.for.in.root
matfile.sortedregions= matfile.sortedregions %>% prep.dir.for.in.root

matt=fread(file=matfile, header=F)
matth=fread(file=matfile.heatmap, header=F) %>% as.data.frame
matts=fread(file=matfile.sortedregions, header=T) %>% as.data.frame


relevant.cols= 7:ncol(matth)
matth.annotation=matth[, 1:6]
matth=matth[rownames(matts), relevant.cols]
cols.per.sample=length(relevant.cols)/nsamples


fcat("annotating reordered peaks")

simpleCache(paste0("reordered_peaks_annotated_matrixid", digest::digest(matts)) %>% addversion, {
matts.annotated=add.peak.annotation.to.mat(matts)
}, assignToVar="matts.annotated", reload=T)
gaps.row=factor(matts[, "deepTools_group"]) %>% as.integer %>% diff %>% as.logical %>% which

 
#ComplexHeatmap::pheatmap(matth[,get.sample.coords(1, extend.from.center=100)], scale="none", cluster_row=F, cluster_col=F, breaks=seq(0,3,.02), gaps_row=gaps.row, annotation_row=matts[, "deepTools_group", drop=F])



#deeptools.sample.umap(1, 30, ncomp=10, md=1)

################################################################################
#make a umap now with each peak containing the information of multiple samples at once
################################################################################



fcat("processing cross sample umap")
peak.clustering.result=deeptools.cross.samples.umap(which.samples.umap, extend.from.center=umap.pars$extend.from.center, ncomp=umap.pars$ncomp, md=umap.pars$md, rs = umap.pars$rs,nn=umap.pars$nn, coll=umap.pars$coll, cluster.resolution=umap.pars$cluster.resolution, matth=matth,  matts=matts.annotated, cols.per.sample=cols.per.sample)
fcat("Done. retrieving result")

fcat("preparing heatmap")

hm=ComplexHeatmap::pheatmap(peak.clustering.result$submat %>% as.matrix, scale="row", cluster_row=F, cluster_col=F, breaks=seq(0,3,.02), gaps_row=gaps.row, annotation_row=peak.clustering.result$umap.df[, c("deepTools_group", "louvain.cluster"), drop=F])
gc()

list(peak.clustering.result=peak.clustering.result, 
  cross.sample.umap=deeptools.cross.samples.umap, 
  matth=matth, 
  matth.annotation=matth.annotation,
  matts=matts.annotated, 
  umap.pars=umap.pars, 
  heatmap=hm
  
)}


get.section=function(matt, section, num.sections=4){
 blocksize=round(nrow(matt)/num.sections) 
cc=1 #record of current matrix index
bc=1 #record of current section
blocks=list()
while(cc< nrow(matt)){ # go on while section is not reached
 if((cc+blocksize)>=nrow(matt)){
   blocks[[bc]]=cc:nrow(matt)
 }else{
   blocks[[bc]]=cc:(cc+blocksize)
 }
 cc=cc+blocksize+1
 bc=bc+1
}  
    sectionmatt= matt[blocks[section] %>% Reduce(c, .), ]
    sectionmatt
 }
  

read.deeptools=function(dtp, annotation.peaks, masks=list()){
  
  matfile=dtp$path.matrix 
  matfile.heatmap=dtp$path.heatmap.matrix
  matfile.sortedregions=dtp$path.sortedregions
  nsamples=length(dtp$all.bigwigs)

    

library(uwot)
library(igraph)
#library(leiden)
# import the matrix
fcat("importing matrices")
matfile=matfile %>% prep.dir.for.in.root
matfile.heatmap= matfile.heatmap %>% prep.dir.for.in.root
matfile.sortedregions= matfile.sortedregions %>% prep.dir.for.in.root

matt=fread(file=matfile, header=F)
matth=fread(file=matfile.heatmap, header=F) %>% as.data.frame
matts=fread(file=matfile.sortedregions, header=T) %>% as.data.frame


simpleCache(paste0("reordered_peaks_annotated_matrixid", digest::digest(matts)) %>% addversion, {
matts.annotated=add.peak.annotation.to.mat(matts)
}, assignToVar="matts.annotated", reload=T)

relevant.cols= 7:ncol(matth)
matth.annotation=matth[, 1:6]
matth=matth[rownames(matts), relevant.cols]


list(
  matth=matth, 
  matth.annotation=matth.annotation,
  matts=matts.annotated
)}
 
#


#deeptools.sample.umap(1, 30, ncomp=10, md=1)

################################################################################
#make a umap now with each peak containing the information of multiple samples at once
################################################################################





################################################################################
# parse deeptools script to find all paths inside
################################################################################

### function in progress.
# Function to parse the script and retrieve paths after specific flags
parse_script_paths <- function(script_path) {
  # Read the script from the provided path
  script_lines <- readLines(script_path)
  
  # Initialize the list to store the paths
  final.list <- list()
  
  # Define the flags you are looking for
  flag_patterns <- list(
    "-S" = "-S\\s+([^\\s+]+)",  # Single path after -S
    "--outFileSortedRegions" = "--outFileSortedRegions\\s+([^\\s]+)",  # Single path after --outFileSortedRegions
    "-R" = "-R\\s+([^\\s]+(\\s+[^-\\s][^\\s]+)*)"  # Multiple paths after -R (before next flag)
  )
  
  # Loop through each line in the script
  for (line in script_lines) {
    # For each flag, find matches using regular expressions
    for (flag in names(flag_patterns)) {
      pattern <- flag_patterns[[flag]]
      
      # Match the pattern to extract the paths
      match <- regmatches(line, regexec(pattern, line))
      
      # If a match is found
      if (length(match[[1]]) > 1) {
        # If the flag is -R (multiple paths), split the result by space
        if (flag == "-R") {
          paths <- unlist(strsplit(match[[1]][2], "\\s+"))
        } else {
          paths <- match[[1]][2]  # Single path for other flags
        }
        
        # Store the paths in the list using the flag as the name
        final.list[[gsub("-", "", flag)]] <- paths
      }
    }
  }
  
  # Return the final list of paths
  return(final.list)
}

################################################################################
#two level clustering workflow
################################################################################

two.level.clustering.workflow= function(matt, cor.threshold, clusk, clusk2, label="", allcolors=NULL){  
corr.tfs=cor(matt)
high.cor.tfs=colnames(corr.tfs)[apply(abs(corr.tfs)>=cor.threshold, 2, function(x) sum(x)>10) ] 
corr.tfs2=corr.tfs[high.cor.tfs, high.cor.tfs]

hm.corr=ComplexHeatmap::pheatmap(corr.tfs2, show_rownames=F, show_colnames=F) %>% draw
################################################################################
# define cluster numbers that will define the tf groupings
################################################################################

tfannot.corr=as.data.frame(cutree(as.hclust(column_dend(hm.corr)), k=clusk)) %>% givecolnames(., nms="tfcorr.cluster") %>% dplyr::mutate(tfcorr.cluster=factor(tfcorr.cluster))

clusters.use.corr=tfannot.corr %>% group_by(tfcorr.cluster) %>% summarise(counts=n()) %>% dplyr::filter(counts>2) %>% pull(tfcorr.cluster) %>% as.vector

hm.corr2=ComplexHeatmap::pheatmap(corr.tfs2, annotation_col=tfannot.corr[colnames(corr.tfs2),,drop=F],
  annotation_row=tfannot.corr[colnames(corr.tfs2),,drop=F], show_rownames=F, show_colnames=F) %>% draw()

################################################################################
# define (and fix) cluster colors
################################################################################


simpleCache(paste_("tf_corr_cluster_colors_numclusters", clusk,"clustersid", digest::digest(clusters.use.corr)) %>% addversion, {
 tfcols=hm.corr2@ht_list[[1]]@top_annotation@anno_list$tfcorr.cluster@color_mapping@colors 
 tfcols
}, assignToVar="tfcols",reload=T)
allcolors$TF.cluster=tfcols
names(allcolors$TF.cluster)=paste0("TF.cluster.", names(tfcols))
allcolors$tfcorr.cluster=tfcols


hm.corr3=ComplexHeatmap::pheatmap(corr.tfs2, annotation_col=tfannot.corr[colnames(corr.tfs2),,drop=F],
  annotation_row=tfannot.corr[colnames(corr.tfs2),,drop=F], show_rownames=F, show_colnames=F, annotation_colors=allcolors, use_raster=T) %>% draw()
################################################################################
# hide one of the legends
################################################################################
hm.corr3@ht_list[[1]]@left_annotation@anno_list$tfcorr.cluster@show_legend=F

################################################################################
# plot clusters at the first level 
################################################################################

################################################################################
#further divide samples into two more clusters
################################################################################

tfannot.corr.2=lapply(tfannot.corr %>% pull(tfcorr.cluster) %>% unique, function(xx){
  subsection=tfannot.corr %>% dplyr::filter(tfcorr.cluster==xx) 
  tfs=subsection %>% rownames
  
  # clustering via the correlation
  tr=hclust(dist(corr.tfs2[tfs, tfs] %>% t))
  # clustering via euclidean distance of the tfbs incidences across blocks
  tr2=hclust(dist(matt[, tfs] %>% t))
  
  tfannot.level2.corr=cutree(tr, k=clusk2) %>% coerce.to.matrix %>% givecolnames(., nms="tfcluster.level2.corr") 
  tfannot.level2.euclidean=cutree(tr2, k=clusk2) %>% coerce.to.matrix %>% givecolnames(., nms="tfcluster.level2.euclidean")
  
  tf.annot.level3= bind_cols(subsection, tfannot.level2.corr, tfannot.level2.euclidean) %>% arrange(tfcorr.cluster, tfcluster.level2.corr, tfcluster.level2.euclidean) %>% dplyr::mutate(TF.cluster=paste_(tfcorr.cluster, tfcluster.level2.corr, tfcluster.level2.euclidean))
}) %>% Reduce(rbind,.)


hm.corr4=ComplexHeatmap::pheatmap(corr.tfs2, annotation_col=tfannot.corr.2[colnames(corr.tfs2),,drop=F],
  annotation_row=tfannot.corr.2[colnames(corr.tfs2),,drop=F], show_rownames=F, show_colnames=F) %>% draw()

################################################################################
# hide one of the legends
################################################################################
hm.corr4@ht_list[[1]]@left_annotation@anno_list$tfcorr.cluster@show_legend=F



################################################################################
# prepare heatmap resources
################################################################################

#tfannotation correlation rearranged by correlation cluster
# to use the original correlation clusters, use tfcorr.cluster as a variable
tfacr=tfannot.corr.2 %>% arrange(TF.cluster) %>% select(TF.cluster) %>% dplyr::mutate(TF.cluster=paste0("TF.cluster.", TF.cluster))
tfacr.names=tfannot.corr.2 %>% arrange(TF.cluster) %>% rownames

clusters.use.corr=tfacr %>% group_by(TF.cluster) %>% summarise(counts=n()) %>% dplyr::filter(counts>2) %>% pull(TF.cluster) %>% as.vector





simpleCache(paste_("tf_corr_cluster_colors_numclusters_twolevels", clusk,clusk2, "clustersid", digest::digest(clusters.use.corr)) %>% addversion, {
 tfcols=hm.corr4@ht_list[[1]]@top_annotation@anno_list$TF.cluster@color_mapping@colors 
 tfcols
}, assignToVar="tfcols",reload=T)
allcolors$TF.cluster=tfcols
names(allcolors$TF.cluster)=paste0("TF.cluster.", names(tfcols))
allcolors$tfcorr.cluster=tfcols



hm.corr5=ComplexHeatmap::pheatmap(corr.tfs2, annotation_col=tfacr[colnames(corr.tfs2),,drop=F],
  annotation_row=tfacr[colnames(corr.tfs2),,drop=F], annotation_colors=allcolors, show_rownames=F, show_colnames=F) %>% draw()


sca=2
tpdf(path=config$plotpath, paste0("heatmap_high_correlation_tfs_numclusters_twolevels_",label, clusk, clusk2, "mincor", cor.threshold) %>% addversion, wi=pw*sca, he=pw*sca)
print(hm.corr5)
dev.off()


################################################################################
#define tf groups with more than 2 tfs to find perturbations for
################################################################################

tfbs.assoc.list.corr=lapply(clusters.use.corr, function(x){
  
tff=tfacr %>% filter(TF.cluster==!!x) %>% rownames

#seriate

tff[seriation::get_order(seriation::seriate(matt[, tff], margin=2, method="PCA"))]

}) %>% givename(., clusters.use.corr)

newtff=tfbs.assoc.list.corr %>% Reduce(c, .)

hmtf=ComplexHeatmap::pheatmap(matt[,newtff ]  , cluster_row=F, cluster_col=F, scale="none", annotation_col=tfacr[newtff,,drop=F], show_rownames=F, show_colnames=F, annotation_row=1:nrow(matt) %>% coerce.to.matrix %>% as.data.frame %>% givecolnames(., nms="block") %>% dplyr::mutate(megablock=factor(cut(block,breaks=seq(0,200, 30), labels=F))), annotation_colors=allcolors, use_raster=T) %>% draw()

sca=2
tpdf(path=config$plotpath, paste_("heatmap_chromatinblocks_vs_tfclusters_two_levels",label, clusk,clusk2, "mincor", cor.threshold) %>% addversion, wi=pw*sca, he=pw*sca)
print(hmtf)
dev.off()


list(heatmap.correlations=hm.corr5,
  heatmap.counts=hmtf,
  tfbs.assoc.list.corr=tfbs.assoc.list.corr
  

  )

}


################################################################################
#
################################################################################
prepare.name.annotation=function(matt, genes.to.label, row.where="right", column.where="bottom", linecolor=NULL, labelcolor=NULL){
  
     linesgp= grid::gpar(col=linecolor)
      labelsgp= grid::gpar(col=labelcolor)
  
  fcat("preparing guide vector list")
  guidelist=lapply(c(1,2), function(namedim){
  
    if(namedim==1){
      geneguide=1:dim(matt)[1]
      }else{
        
      geneguide= 1:dim(matt)[2]}
    
    fcat("geneguide is of length", length(geneguide))
    fcat("dimnames(matt)[[",namedim, "]] is of length", length(dimnames(matt)[[namedim]]))
    
    if(is.null(dimnames(matt)[[namedim]])){
      fcat("dimnames not found. filling in")
      dimnames(matt)[[namedim]]=1:dim(matt)[namedim]
    }
    
    names(geneguide)=dimnames(matt)[[namedim]]
    
    geneguide
  })
    #mrkrs are genes of interest
 #dput(guidelist)
  
      getposition.row=Vectorize(function(x) guidelist[[1]][x], USE.NAMES=F)
      getposition.col=Vectorize(function(x) guidelist[[2]][x], USE.NAMES=F)
      
      fcat("preparing annotation objects")
      selected.rows=intersect(names(guidelist[[1]]), genes.to.label)
      selected.cols=intersect(names(guidelist[[2]]), genes.to.label)
      
      fcat("selected rows")
      print(selected.rows)
      
      
       fcat("selected columns")
      print(selected.cols)
      
      
      
    if(length(selected.rows)>0){
      
    row.annomark=  ComplexHeatmap::anno_mark(at = guidelist[[1]][selected.rows] %>% unname %>% as.numeric, 
                                                                       labels = selected.rows,
                                                                       lines_gp= linesgp, 
                                                                       labels_gp= labelsgp)
      
      
    rowannot = ComplexHeatmap::HeatmapAnnotation(which="row", foo = row.annomark)
    }else{
     rowannot=NULL 
     row.annomark=NULL
    }
    
    
     if(length(selected.cols)>0){
       column.annomark=ComplexHeatmap::anno_mark(at = guidelist[[2]][selected.cols] %>% unname %>% as.numeric, 
                                                                       labels = selected.cols,
                                                                       lines_gp= linesgp,
                                                                       labels_gp= labelsgp, 
                                                                                  side="bottom")
      columnannot=ComplexHeatmap::HeatmapAnnotation(which="column", collabels = column.annomark ) 
      }else{
        columnannot=NULL
        column.annomark=NULL
        }

        
list(row=rowannot, row.annomark=row.annomark, column=columnannot, column.annomark=column.annomark)
}
    




getlast= function(x, n) x[(length(x)-(n-1)): length(x)]



zero2na= function(x) ifelse(x==0, NA, x)


################################################################################
# transform a factor into a logical matrix of true or false
################################################################################

create_logical_matrix <- function(factor_vector) {
  # Get the unique levels (categories) in the factor
  categories <- levels(factor_vector)
  
  # Create a matrix where each entry is TRUE if it matches the category, FALSE otherwise
  logical_matrix <- outer(factor_vector, categories, FUN = "==")
  
  # Set column names to be the category names
  colnames(logical_matrix) <- categories
  
  return(logical_matrix)
}



create_colored_list <- function(factor_vector) {
  # Get the unique levels (categories) of the factor
  categories <- levels(factor_vector)
  
  # Create a named list where each element has the vector c("TRUE"="firebrick", "FALSE"="white")
  colored_list <- setNames(
    replicate(length(categories), c("TRUE" = "firebrick", "FALSE" = "white"), simplify = FALSE),
    categories
  )
  
  return(colored_list)
}

################################################################################
#package variables in environment into a list
################################################################################
package.vars <- function(...) {
  # Capture the input arguments as a list
  args <- list(...)
  # Use names of arguments from their deparse expressions
  names(args) <- sapply(substitute(list(...))[-1], deparse)
  return(args)
}



################################################################################
# pca random forest workflow
################################################################################

random.forest.workflow= function(envlist, dataset.info, recreate=F){

  
config=dataset.info  
test.datasets=dataset.info$dataset_name[as.logical(dataset.info$is_test)]
allcolors=envlist$allcolors
pca.data=envlist$pca.data
refids=envlist$refids

  
gini.threshold=dataset.info$randomforest_gini_threshold

set.seed(dataset.info$randomseed)

generic.pc.list=paste0("PC", 1:length(envlist$refids))


 
# Train the Random Forest classifier
traininglevels=intersect(names(allcolors[[dataset.info$cell_type_variable]]), unique(pca.data[refids,"Condition" ]))   
fcat(traininglevels)

simpleCache(paste_("randomforest_model_seed",config$randomseed, "initialcomponents", config$randomforest_initial_components, "target",config$cell_type_variable, "testds_",paste(test.datasets, collapse="-"),"version", analysis.version), { 

rf_model <- randomForest(x = pca.data[refids, generic.pc.list], y = factor(pca.data[refids, config$cell_type_variable], levels=traininglevels ), ntree = config$randomforest_ntrees)

rf_model
}, assignToVar="rf_model", reload=!recreate, recreate=recreate)

# Print the trained model
print(rf_model)

rf.err.rate.generic=rf_model$err.rate[nrow(rf_model$err.rate), "OOB"]

confusion.training=rf_model$confusion %>% as.data.frame %>% dplyr::select(-class.error)
confusion.training2=lapply(1:nrow(confusion.training), function(x) confusion.training[x, ]/sum(confusion.training[x,])) %>% Reduce(rbind,. )
fcat("checkpoint A")

sca=1.5
rr=600
hmpgen=ComplexHeatmap::pheatmap(confusion.training2[traininglevels, traininglevels], cluster_col=F, cluster_row=F, main=paste0("Training set. All PCs. Error:", rf.err.rate.generic))
tpdf(path=config$plotpath, paste0("heatmap_randomForest_trainingset_confusionmatrix_analysisversion", analysis.version), wi=pw*sca, he=pw*sca)
print(hmpgen)
dev.off()

# Get variable importance
var_importance <- importance(rf_model)

# Print variable importance
print(var_importance)

# Plot variable importance
fcat("checkpoint B")
sca=1.5
var_imp_plot <- varImpPlot(rf_model)
arr.varimp <- var_imp_plot %>% as.data.frame %>% names2col(., "PC")%>% dplyr::mutate(rank=1:n()) %>% dplyr::filter(MeanDecreaseGini>=1, rank<=20) %>% arrange(MeanDecreaseGini)
ginilevels.generic=arr.varimp %>% pull(PC) 
tpdf(path=config$plotpath, paste_("var.importance.plot.generic", analysis.version), wi=pw*sca, he=pw*sca*1)

print(
  
 ggplot(arr.varimp, aes(x=MeanDecreaseGini, y=factor(PC, levels=ginilevels.generic)))+geom_point(color="black")+cowplot::theme_cowplot()+cowplot::background_grid()+
    geom_vline(xintercept=gini.threshold, color="red")
)  
dev.off()
#var_imp_plot
####prediction of test samples
fcat("checkpoint C")
test_data=pca.data %>% filter(dsname %in% test.datasets) %>% dplyr::select( all_of(generic.pc.list))

testids=pca.data %>% filter(dsname %in% test.datasets)  %>% arrange(tmnt.full) %>% rownames

# Make predictions on the test set
predictions.prob <- predict(rf_model, test_data, type="prob")[testids,]
predictions <- predict(rf_model, test_data)
# Create the confusion matrix
conf_matrix <- table(Actual = rownames(test_data), Predicted = predictions)

color_palette <- colorRampPalette(c("blue", "white", "red"))(100)

# Define the range for the color bar
# For example, we want the colors to range from -2 to 2
#breaks <- seq(0, 1, length.out=101)
sca=1.5
rr=600

hmpprobs.generic=ComplexHeatmap::pheatmap(predictions.prob, cluster_col=F, cluster_row=F, main=paste0("test prediction probabilities. error:",rf.err.rate.generic))
tpdf(path=config$plotpath, paste0("heatmap_randomForest_test_generic_predictionprobs_analysisversion", analysis.version), wi=pw*sca*1.1, he=pw*sca)
print(hmpprobs.generic)
dev.off()


hmpprobs.generic

fcat("checkpoint D")
################################################################################
# selecting pcs based on their Mean decrease Gini
################################################################################
var.importance.arranged=var_importance %>% as.data.frame %>% arrange(-MeanDecreaseGini)

#for versions before 20250810, use onlly rf_model inside digest
simpleCache(paste0("informative.pcs_rfmodel", digest::digest(list(rf_model, dataset.info)), "analysisversion_", analysis.version), {
informative.pcs=rownames(var.importance.arranged)[var.importance.arranged>=gini.threshold]
informative.pcs}, assignToVar="informative.pcs",reload=!recreate) 

fcat("checkpoint E")
################################################################################
# use the informative  pcs togenerate a second cleaner model with sleected features
################################################################################


simpleCache(paste_("randomforest_model_informativepcs_seed",config$randomseed, "components", config$randomforest_initial_components,"gini", dataset.info$randomforest_gini_threshold, digest::digest(informative.pcs), "target",config$cell_type_variable, "testds_",paste(test.datasets, collapse="-"),"version", analysis.version), { 

rf_model_informative <- randomForest(x = pca.data[refids, informative.pcs], y = factor(pca.data[refids, config$cell_type_variable], levels=traininglevels), ntree = config$randomforest_ntrees)
rf_model_informative
}, assignToVar="rf_model_informative", reload=!recreate)



rf.err.rate.informative=rf_model_informative$err.rate[nrow(rf_model_informative$err.rate), "OOB"]

print(rf_model_informative)

confusion.trainingi=rf_model_informative$confusion %>% as.data.frame %>% dplyr::select(-class.error)
confusion.training2i=lapply(1:nrow(confusion.trainingi), function(x) confusion.trainingi[x, ]/sum(confusion.trainingi[x,])) %>% Reduce(rbind,. )

fcat("checkpoint 1")
sca=1.5
rr=600
hm.training.informative=ComplexHeatmap::pheatmap(confusion.training2i[traininglevels, traininglevels], cluster_col=F, cluster_row=F, main=paste0("Training set (Informative PCs). Error:", rf.err.rate.informative))
tpdf(path=config$plotpath, "heatmap_randomForest_trainingset_confusionmatrix_informativepcs_analysisversion" %>% addversion, wi=pw*sca, he=pw*sca)
print(hm.training.informative)
dev.off()
print(hm.training.informative)
fcat("checkpoint 2")
# Get variable importance
var_importancei <- importance(rf_model_informative)

# Print variable importance
print(var_importancei)

# Plot variable importance
var_imp_ploti <- varImpPlot(rf_model_informative)
print(var_imp_ploti)



# Make predictions on the test set
predictions.probi <- predict(rf_model_informative, test_data, type="prob")[testids,]
predictionsi <- predict(rf_model_informative, test_data)
# Create the confusion matrix
conf_matrixi <- table(Actual = rownames(test_data), Predicted = predictionsi)


# Define the range for the color bar
# For example, we want the colors to range from -2 to 2
#breaks <- seq(0, 1, length.out=101)
edm2.samples=pca.data %>% dplyr::filter(dscategory=="test", tmnt2!="nt") %>% rownames
edm1.samples=pca.data %>% dplyr::filter(dscategory=="test",tmnt2=="nt", tmnt.full!="nt_nt") %>% rownames


ery.order.edm1=predictions.probi[edm1.samples, "Ery"] %>% sort(., decreasing=T) %>% names
ery.order.edm2=predictions.probi[edm2.samples, "Ery"] %>% sort(., decreasing=T) %>% names

sca=1.5
rr=600
hm.inf1=ComplexHeatmap::pheatmap(predictions.probi[ery.order.edm1, traininglevels], cluster_col=F, cluster_row=F, main="test prediction probabilities", annotation_row=pca.data[ery.order.edm1, dataset.info$annotvars ]  %>% dplyr::mutate(rbc.round1=factor(rbc.round1)),annotation_colors=allcolors, show_rownames=F, breaks=seq(0,.3,0.02))
hm.inf2=ComplexHeatmap::pheatmap(predictions.probi[ery.order.edm2,traininglevels ], cluster_col=F, cluster_row=F, main="test prediction probabilities", annotation_row=pca.data[ery.order.edm2, dataset.info$annotvars ] %>% dplyr::mutate(rbc.round1=factor(rbc.round1)), annotation_colors = allcolors, show_rownames = F, breaks=seq(0,.3,.02))

#round(edm1.samples/nrow(predictions.probi))

fcat("checkpoint 3")
tpdf(path=config$plotpath, paste0("heatmap_randomForest_test_EDM1_predictionprobs_informativepcs_analysisversion", analysis.version), wi=pw*sca*1.5, he=pw*sca*1.5)
print(hm.inf1)
dev.off()

fcat("checkpoint 4")
tpdf(path=config$plotpath, paste0("heatmap_randomForest_test_EDM2_predictionprobs_informativepcs_analysisversion", analysis.version), wi=pw*sca*1.5, he=pw*sca*1.5)
print(hm.inf2)
dev.off()


  result_list=list(
    var_importance= var_importance,
    var_importancei=    var_importancei,
    rf_model=    rf_model,
    rf_model_informative=    rf_model_informative,
    rf.err.rate.generic=    rf.err.rate.generic,
confusion.training2=confusion.training2,
predictions.prob=predictions.prob,
hmpprobs.generic=hmpprobs.generic,
var.importance.arranged=var.importance.arranged,
var.importance.arranged=var.importance.arranged,
confusion.training2i=confusion.training2i,
rf.err.rate.informative=rf.err.rate.informative,
predictions.probi=predictions.probi[, traininglevels],
    conf_matrix=conf_matrix,
conf_matrixi=conf_matrixi, 
    informative.pcs=informative.pcs
    )


  
  
  
}


################################################################################
# peak roc analysis workflow
################################################################################


compute.aucs.workflow <- function(input_list) {
  
library(parallel)
library(digest)
library(pROC)
library(simpleCache)
library(data.table)
  
  
  # Extract variables from the input list
  positives.celltypes <- input_list$positives.celltypes
  peaks <- input_list$peaks
  meta <- input_list$meta
  refids <- input_list$refids
  analysis.version <- input_list$analysis.version
  max.peaks.pull <- input_list$max.peaks.pull
  norm.peak.vals <- input_list$norm.peak.vals
  mc.cores<- input_list$mc.cores
  
  
  get.sample.peakvals=function(peakid, samples=NULL){
  
  if(!is.null(samples)){ mt=norm.peak.vals[ peakid, samples]}else{ mt=norm.peak.vals[peakid, !colnames(norm.peak.vals) %in% c("Geneid", "Chr", "Start", "End", "Strand", "Length")]}
  
  mt %>% t
  
}
  
  
  # Helper functions
  is.target <- function(x, targets) fifelse(x %in% targets, 1, 0)
  
  make.target.vector <- function(x) is.target(meta[refids, x$classvar], x$targets)
  
  # Precompute target vectors
  target.vector.list <- lapply(positives.celltypes, make.target.vector)
  
  # Precompute peak matrix
  sp <- get.sample.peakvals(peakid = peaks, samples = refids) %>% 
    cbind(., meta[refids, ])

  
  # Compute collective AUCs
  simpleCache(paste0("collective_aucs_celltypes_all_peaksid_", digest::digest(peaks)) %>% addversion, {
    
    parallel.auc.df <- mclapply(names(positives.celltypes), function(x, spv, pks, tvl, av) {
      
      simpleCache(paste0("collective_aucs_celltypes_top", max.peaks.pull, "perPC_target", x, 
                         "_peaks_id", digest::digest(pks), "analysis", av), {
        aucs <- lapply(pks, function(peak) {
          roc(response = tvl[[x]], predictor = spv[, peak], direction = "<", quiet = TRUE)$auc
        }) %>% unlist()
        
        data.frame(list(Geneid = pks, auc = aucs, target = x))
        
      }, assignToVar = "chunk.aucs", reload = TRUE)
      chunk.aucs
    }, spv = sp, pks = peaks, tvl = target.vector.list, mc.cores = mc.cores, av = analysis.version, mc.cleanup = TRUE)
    
    collective.aucs.celltypes <- Reduce(rbind, parallel.auc.df)
    
  }, assignToVar = "collective.aucs.celltypes", reload = TRUE)
  
  # Return the final object
  return(collective.aucs.celltypes)
}






################################################################################
# function to ensure that some packages are loaded, else install them, else fail
################################################################################

ensure.packages <- function(packages) {
  # Ensure BiocManager is installed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  # Function to load or install a package
  load_or_install <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(sprintf("Package '%s' is not installed. Attempting to install...", pkg))
      tryCatch(
        {
          BiocManager::install(pkg, ask = FALSE, update = FALSE)
          if (!requireNamespace(pkg, quietly = TRUE)) {
            install.packages(pkg)
          }
        },
        error = function(e) {
          stop(sprintf("Failed to install package '%s': %s", pkg, e$message), call. = FALSE)
        }
      )
    }
    # Try loading the package
    if (!require(pkg, character.only = TRUE)) {
      stop(sprintf("Failed to load package '%s' after installation attempts.", pkg), call. = FALSE)
    }
  }
  
  # Apply the function to each package
  invisible(lapply(packages, load_or_install))
}




get.perturbation.heatmaps=function(list.tfbs, gsetpars=load.gene.signatures(), enrichments=NULL, keyword=config$gsea_celltype_searchterm, label=config$target_cell_type, exclude=config$gsea_celltype_exclude, tf.annot.external=NULL, colorlist=NULL, cluster=F){
if(is.null(enrichments)){
enrichments=get.enrichments.cache(list.tfbs, gsetpars)
}
perts=c("LIGPERTUP|LIGPERTLINCSUP", "LIGPERTDOWN|LIGPERTLINCSDN", "CGP")
perturbation.heatmaps=lapply(names(list.tfbs), function(cluslab){
fcat("sensitivity group", cluslab, "...")
lapply(perts, function(searchterms){
fcat("Looking for groups enrichments in ", searchterms, "...")
endf= enrichments$df %>% filter(grepl(!!searchterms, geneset))
#we make the enrichment heatmap matrix whereby we see what gene sets share the same hits for biological insight.
#we keep the matrix to make a more tailored heatmap later.
fcat("preparing the enrichment matrix...")
eh=make.enrichment.heatmap(x=list.tfbs[[cluslab]], gsetpars=gsetpars, endf %>% dplyr::filter(clusterlabel==!!cluslab))
if( !is.null(eh$matrix) ){
mat=remove.zerocols(remove.zerorows(eh$matrix %>% t))
gset.annot=endf %>% dplyr::filter(clusterlabel==!!cluslab) %>% filter(original.geneset %in% eh$original.gset.names) %>% select(jaccard_index, pval)
fcat("finding the relative relevance of each gene in the genesets for cell type ", target_cell_type, " ...")
relev=get.gene.gsea.relevance(colnames(mat), gsetpars, simple.output=F, keyword=keyword, label=target_cell_type, exclude=exclude)
tf.annot=relev$df %>% givecolnames(., nms=paste0("relevance.for.", target_cell_type))
if(!is.null(tf.annot.external)){
tf.annot=tf.annot %>% bind_cols(., tf.annot.external[rownames(tf.annot),,drop=F ])
}
#
fcat("ranking the perturbations to put those hitting the most relevant TFs on top")
#find the profiles that most closely involve target-relevant genes
hammings=lapply(1:nrow(mat), function(j){
sum(!(tf.annot[,1]>0 & as.logical(mat[j,])))
}) %>% unlist
horder=order(hammings, decreasing=F)
#col_fun <- colorRamp2(c(0, 1), c( "white", "green"))
fcat("preparing annotations for the perturbation heatmap...")
column_ha <- HeatmapAnnotation(
annotation = anno_simple(tf.annot[,1]),
show_annotation_name = TRUE
)
row_ha <- ComplexHeatmap::rowAnnotation( jaccard.index = anno_barplot(gset.annot))
fcat("executing the perturbation heatmap...")
hm.perturbations=ComplexHeatmap::pheatmap(mat[horder, ], col=c("darkmagenta", "yellow"), cluster_row=F, cluster_col=cluster, main=paste0("potential perturbations to ", ifelse(grepl("UP", searchterms), "activate", "deactivate")," TFs at ",cluslab , " peaks"), annotation_col=tf.annot, right_annotation=row_ha, annotation_colors=colorlist)
tpdf(path=config$plotpath, paste0("heatmap_tfbs_",  "_target_", target_cell_type, label, cluslab,"_perturbation_signatures", remove.pipes(searchterms)) %>% addversion, wi=pw*sca*3, he=pw*sca*1)
print(hm.perturbations)
dev.off()
}else{
mat=NULL
hm.perturbations=NULL
}
list(matrix=mat, heatmap=hm.perturbations)
}) %>% givename(.,perts %>% remove.pipes())
}) %>% givename(., names(list.tfbs))
perturbation.heatmaps}





###############################################################################
# get the ath submatrix of a larger matrix equally divided into n submatrices
###############################################################################


get.submatrix <- function(mat, a, n) {
  # Check if inputs are valid
  if (!is.matrix(mat)) stop("Input must be a matrix.")
  if (a <= 0 || a > n) stop("Argument 'a' must be between 1 and 'n'.")
  if (n <= 0) stop("Argument 'n' must be a positive integer.")
  
  # Determine column ranges for splitting
  total_cols <- ncol(mat)
  col_ranges <- split(1:total_cols, cut(1:total_cols, n, labels = FALSE))
  
  # Extract the requested submatrix
  submatrix <- mat[, col_ranges[[a]], drop = FALSE]
  return(submatrix)
}

# get an interval of -span:span of the middle of a matrix
get.half.interval=function(matt, span=10){ md=round(ncol(matt)/2); matt[,(md-span):(md+span)]}

get.center.sums=function(matt, nsubmats=3, span=10){
  lapply(1:nsubmats, function(x)
get.half.interval(get.submatrix(matt %>% as.matrix, x, nsubmats), span) %>% rowSums
)
}

simple.heatmap=function(matt, ...){ ComplexHeatmap::pheatmap(matt, cluster_row=F, cluster_col=F, scale="none", ...)}
row.heatmap=function(matt, ...){ ComplexHeatmap::pheatmap(matt, cluster_row=F, cluster_col=F, scale="row", ...)}
cluster.heatmap=function(matt, ...){ ComplexHeatmap::pheatmap(matt, cluster_row=T, cluster_col=T, scale="none", ...)}





ratio.to.first <- function(input_matrix) {
  # Ensure the input is a matrix
  if (!is.matrix(input_matrix)) {
    stop("Input must be a matrix.")
  }
  
  # Get the number of columns in the matrix
  num_cols <- ncol(input_matrix)
  
  # Ensure the matrix has at least two columns
  if (num_cols < 2) {
    stop("Matrix must have at least two columns.")
  }
  
  # Extract the first column
  first_col <- input_matrix[, 1]
  
  # Compute the average of the second-to-last column
  second_to_last_col <- input_matrix[, num_cols - 1]
  avg_second_to_last <- mean(second_to_last_col)
  
  # Divide the average of the second-to-last column by the first column
  result <- avg_second_to_last / first_col
  
  return(result)
}



dump.table=function(x, ...){
  variable_name <- deparse(substitute(x))
  fwrite(x, file=paste0(config$out_root, "/metadata/", variable_name, ".csv"), ...)

}

###############################################################################
# seriate columns of a matrix
##############################################################################

seriatecols=function(x) x[,seriation::get_order(seriation::seriate(x, margin=2, method="PCA"))]

seriaterows=function(x, reverse.order=F){ 
  
  ord=seriation::get_order(seriation::seriate(x, margin=1, method="PCA"))
  if(reverse.order){
    ord=rev(ord)
  }
  x[ord,]
 
   
}

cumulative.lengths <- function(lst) {
  # Calculate the lengths of each vector
  lengths <- sapply(lst, length)
  
  # Compute the cumulative sum of lengths
  cumulative <- cumsum(lengths) %>% unname
  
  # Return the result as a named list

  
  
  return(cumulative)
}


################################################################################
# generate letter aliases for values in a column (to use as proxy for treatment)
################################################################################

make_alias <- function(df, colname, na.values="NA", na.replaceby=NA, replace=T) {
  # Ensure colname is treated as a string
  col <- df[[colname]]
  
  # Get unique sorted values

  uniq_vals <- unique(col)
   uniq_vals <- uniq_vals[!(uniq_vals %in% na.values)] 
     uniq_vals <- sort(uniq_vals)
 
  # Generate letter aliases (handles >26 values too)
  n <- length(uniq_vals)
  aliases <- sapply(seq_len(n), function(i) {
    catalog=c(LETTERS, paste0(LETTERS,"II"),paste0(LETTERS,"III") )

    excel_style <- catalog[i]
    excel_style
  })
  
  
  # Map values to aliases
  alias_map <- setNames(aliases, uniq_vals)
  
  
 alias_map <- c(alias_map, rep(na.replaceby, length(na.values)) %>% givename(., na.values))
 
 fcat(dput(alias_map))
  
 
 if(replace){
   #we store the original column with an extension
     new_colname <- paste0(colname, "_original")
  df[[new_colname]] <-df[[colname]]
   df[[colname]] <-alias_map[as.character(col)]
 }else{
  # Create new alias column
  new_colname <- paste0(colname, "_alias")
  df[[new_colname]] <- alias_map[as.character(col)]
 }
 
  return(df)
}


################################################################################
# function to generate treatment columns
################################################################################

generate_treatment_columns <- function(mat) {
  # Convert matrix to a data frame for easier handling
  df <- as.data.frame(mat)
  
  # Get the suffixes from column names
  suffixes <- unique(sub(".*\\.", "", colnames(df)))
  
  # Process each suffix
  for (suf in suffixes) {
    # Find columns with the current suffix
    cols_with_suffix <- grep(paste0("\\.", suf, "$"), colnames(df), value = TRUE)
    
    # Sort columns alphabetically and paste row values together
    df[[paste0("tmnt.", suf)]] <- apply(df[cols_with_suffix], 1, function(row) {
      # Replace NAs with "NA", then paste
      paste0(ifelse(is.na(row), "NA", row), collapse = "")
    })
  }
  
  # Create the tmnt.full column by combining all "tmnt.*" columns
  treatment_columns <- grep("^tmnt\\.\\d+$", colnames(df), value = TRUE)
  df$tmnt.full <- apply(df[treatment_columns], 1, function(row) {
    paste(row, collapse = "_")
  })
  
  # Return the updated matrix
  return(df)
}



#############################################
# second version generates round column
#############################################
generate_treatment_columns2 <- function(data) {
  # Identify the unique suffixes in the column names
  suffixes <- unique(gsub(".*\\.", "", grep("\\.", colnames(data), value = TRUE)))
  
  # Loop through the suffixes to create the treatment columns
  for (suffix in suffixes) {
    # Identify the columns with the current suffix
    cols_with_suffix <- grep(paste0("\\.", suffix, "$"), colnames(data), value = TRUE)
    # Sort columns alphabetically and concatenate their values
    data[[paste0("tmnt.", suffix)]] <- apply(data[cols_with_suffix], 1, function(row) {
      paste0(ifelse(is.na(row), "NA", row), collapse = "")
    })
  }
  
  # Add the "tmnt.full" column
  tmnt_cols <- grep("^tmnt\\.\\d+$", colnames(data), value = TRUE)
  data$tmnt.full <- apply(data[tmnt_cols], 1, function(row) {
    paste(row, collapse = "_")
  })
  
  # Add the "round" column
  data$round <- apply(data[tmnt_cols], 1, function(row) {
    # Determine the last non-"NANANA" column
    last_non_na_round <- max(which(row != "NANANANANA"), na.rm = TRUE)
    if (is.finite(last_non_na_round)) {
      paste0("Round", last_non_na_round)
    } else {
      NA # Case where all columns are "NANANA"
    }
  })
  
  return(data)
}



################################################################################
#
################################################################################


generate_treatment_columns3<- function(pd, dataset.info){
  pca.data <-pd
all.reagent.cols=dataset.info$reagents
treatment.id.cols=dataset.info$treatment_id_cols
#infer the number of steps in the experiment
treatment.steps= strsplit(all.reagent.cols, split="\\.") %>% lapply(., function(x) x[length(x)]) %>% unlist %>% unique



## collect all treatment reagents per step
treatments.list=lapply(treatment.steps, function(x){
  l=list(
  reagents=all.reagent.cols[grepl(x, all.reagent.cols)]
  )
  l
}) %>% givename(., treatment.steps)

#input is the number of the treatment round
get.treatment.na= function(tn){
rep(NA, length(treatments.list[[tn]]$reagents)) %>% paste(.,collapse="" )  
}



newcols=list()
#newcols.aliases=list()
for(k in 1:length(treatments.list)){
fcat("loop ", k)
  coln=paste0("tmnt", k)
  pca.data=pca.data %>% dplyr::mutate(!!sym(coln) := apply(select(., all_of(treatments.list[[k]]$reagents)), 1, function(x) paste(x, collapse="")))
  
  #make human readable aliases for the  treatment
  pca.data=pca.data %>% make_alias(., coln, na.values=get.treatment.na(k), na.replaceby="nt", replace=T) #make an alias but store the original
newcols[[k]]=coln
#newcols.aliases[[k]]=paste0(coln, "_alias")
  }


newcols= newcols %>% Reduce(c, .)

treatments.list= treatments.list %>% givename(., newcols)

################################################################################
# make round column
################################################################################
pca.data$round= lapply(newcols, function(k){
  ifelse(pca.data[[k]]==get.treatment.na(k),0,1) 
}) %>% bind_cols %>% as.data.frame %>% givecolnames(.,nms=newcols) %>% rowSums

make.round.text=function(rnd) ifelse(is.na(rnd), NA, paste0("Round", rnd))
pca.data$round= make.round.text(zero2na(pca.data$round))



pca.data= pca.data %>% dplyr::mutate(tmnt.full=apply(select(., all_of(newcols)), 1, paste, collapse = "_"), celltype=ifelse(dscategory=="reference",Condition, as.character(dsname)))

pca.data

}



get.last.element=function(x) x[length(x)]

get.element= function(a,rr, cc ) x[rr, cc]
get.row= function(a, rr) x[,cc]
get.column= function(a, cc) x[,cc ]

get.nonzero= function(x) x[x!=0]



################################################################################
# sum rows with the same id
################################################################################
sum.rows.by.id <- function(mat, ids) {
  fcat("enter sum.rows.by.id")
  # Convert to data frame if rownames are provided as IDs
  if (!missing(ids)) {
    mat <- cbind(ids=ids, as.data.frame(mat))
  }
  
  # Group by the 'ids' column and sum all numeric columns
  aggregated.matrix <- aggregate(. ~ ids, data =mat, FUN = sum)
  
  # Convert back to a matrix if necessary
  rownames(aggregated.matrix) <- aggregated.matrix$ids
  aggregated.matrix <- as.matrix(aggregated.matrix %>% select(-ids))# Remove the 'ids' column
  fcat("exit sum.rows.by.id")
  return(aggregated.matrix)
}


################################################################################
# search files matching strings in a vector inside a file
################################################################################


search.files.in.folder <- function(strings, folder_path) {
  # Ensure the folder path exists
  if (!dir.exists(folder_path)) {
    stop("The specified folder does not exist.")
  }
  
  # Get the list of files in the folder
  all_files <- list.files(folder_path, full.names = FALSE)
  
  # Create a function to find matching files for each string
  find_matches <- function(string) {
    matches <- grep(string, all_files, value = TRUE)
    if (length(matches) == 0) {
      return(NA)  # If no match, return NA
    } else {
      return(paste(matches, collapse = ","))  # If multiple matches, concatenate with commas
    }
  }
  
  # Apply the function to each string
  result <- sapply(strings, find_matches)
  
  # Return the vector of results
  return(result)
}



################################################################################
#  get cis peaks
################################################################################


get.cis.peaks=function(args){
 #find the closest peaks
  
  present.peaks=intersect(rownames(args$mtf_mat_gathered), args$peaks)
  
      if(length(present.peaks)< length(args$peaks)){
   fcat("warning: not all peaks found in motif matrix") 
  }
  closest.genes=args$closest_tss_global %>% filter(Geneid %in% present.peaks) %>% arrange(factor(Geneid, levels=present.peaks)) %>% select(Geneid, gene_symbol) 

  
  relevant.motif.freqs= apply(args$mtf_mat_gathered[present.peaks, , drop=F]>0 , 2, sum) %>% get.nonzero %>% sort(., decreasing=T)
  
  relevant.motif.names= apply(args$mtf_mat_gathered[present.peaks, , drop=F]>0 , 2, sum) %>% get.nonzero %>% sort(., decreasing=T) %>% names
  
  
  ## sum motifs from regions that are nearest to the same tf  
sum.mat=sum.rows.by.id(mtf_mat_gathered[closest.genes$Geneid, relevant.motif.names], closest.genes$gene_symbol) 

sum.mat=sum.mat[(rownames(sum.mat) %in% relevant.motif.names),]  
  ### which tfs regulate their own regions
  
  which.are.cis=intersect(closest.genes %>% pull(gene_symbol), relevant.motif.names)
  list(
      gene.motif.mat=sum.mat,
    closest.genes=closest.genes,
    relevant.motif.names=relevant.motif.names,
    which.are.cis=which.are.cis,
  percentage.cis.regions= sum(closest.genes$gene_symbol %in% which.are.cis)/nrow(closest.genes),
  percentage.cis.tfs= sum(relevant.motif.names %in% which.are.cis)/length(relevant.motif.names)
  
  )
  

  
}


get.tfs=function(args){
 #find the closest peaks
  
  closest.genes=args$closest_tss_global %>% filter(Geneid %in% args$peaks) %>% arrange(factor(Geneid, levels=args$peaks)) %>% select(Geneid, gene_symbol) %>% dplyr::mutate(is.tf=istf(gene_symbol)) 

list(
  peak.info=closest.genes,
 tf.info= closest.genes %>% filter(is.tf),
 counts= closest.genes %>% filter(is.tf) %>% group_by(gene_symbol) %>% summarise(counts=n()) %>% arrange(-counts),
 tfs=closest.genes %>% filter(is.tf) %>% group_by(gene_symbol) %>% summarise(counts=n()) %>% arrange(-counts) %>% pull(gene_symbol) %>% unique
  )
  
}
  
 
stdplot=function(x, title, w=1, h=1, sca=1.5){
 tpdf(path=config$plotpath, title, wi=pw*sca*w, he=pw*sca*h)
  print(x)
  dev.off()
}
  
  
capvalues=function(x, lim=5) ifelse(x>=lim, lim, x)



get.targets= function(tf, pks){
  
submat=mtf_mat_gathered[intersect(rownames(mtf_mat_gathered), pks), tf, drop=F]  

  targets=rownames(submat[submat[, 1]>0,, drop=F])
}


get.targets.mtfs= function(tf, pks){
  
submat=mtf_mat_gathered[intersect(rownames(mtf_mat_gathered), pks), tf, drop=F]  

  mtfs=submat[submat[, 1]>0,1]
  list(
  numbers=mtfs, 
  mean=mean(mtfs),
  std=sd(mtfs),
  median=median(mtfs),
  p25=quantile(mtfs, .25),
   p75=quantile(mtfs, .75), 
  )
  
}

shared.regulation= function(tf1, tf2, pks, mtf_mat_gathered){
  
submat=mtf_mat_gathered[intersect(rownames(mtf_mat_gathered), pks), c(tf1, tf2 ), drop=F] %>% as.data.frame %>% dplyr::mutate(both=(!!sym(tf1)>0 & !!sym(tf2)>0), either=(!!sym(tf1)>0 | !!sym(tf2)>0))

pcts=submat %>% summarise( both.pct= sum(both)/n(), either.pct=sum(either)/n())

both.peaks=submat %>% dplyr::filter(both) %>% rownames
either.peaks=submat %>% dplyr::filter(either) %>% rownames
both= length(both.peaks) 
either=length(either.peaks)

      either.pct=pcts$either.pct
    both.pct=pcts$both.pct
  package.vars(
    tf1,
    tf2,
    both.peaks,
    either.peaks,
    both,
    both.pct,
    either,
either.pct
  )
  
}


filter.matrix <- function(cor_matrix, threshold) {
  # Check if the input is a matrix
  if (!is.matrix(cor_matrix)) {
    stop("The input must be a correlation matrix.")
  }
  
  # Ensure the diagonal is excluded
  diag(cor_matrix) <- 0
  
  # Identify the rows and columns to keep based on the threshold
  keep_indices <- apply(cor_matrix, 1, function(row) any(abs(row) >= threshold))
  
  # Filter the matrix by the identified indices
  filtered_matrix <- cor_matrix[keep_indices, keep_indices]
  diag(filtered_matrix) <- 1
  return(filtered_matrix)
}



list.vals.in.cols= function(x){
 
  vals=lapply(colnames(x), function(xx){
  x[, xx] %>% unname  
  }) %>% Reduce(c, .) %>% as.vector() %>% unique 
  
  vals
   
}



make.cache.code <- function(var, tt, varstring)  {
  # Construct the string for simpleCache code
  code <- sprintf(
    "simpleCache('%s', { %s <- %s }, assignToVariable = TRUE, reload = TRUE)", 
    tt, varstring, var
  )
  
  # Return the generated code as a string
  return(code)
}




make.hyper.barplot=function(hdf, genesets=c("LIGPERTDOWN", "LIGPERTLINCSDN"), numbars=10, color="violet",  x="jaccard_index", colorvar="clusterlabel"){
prep.df=hdf %>% dplyr::filter(tlog.pval>0, geneset %in% !!genesets) %>% group_split(geneset, clusterlabel) %>% lapply(., function(x) x %>% arrange(-jaccard_index) %>% group_by(geneset) %>%  dplyr::mutate(rank=1:n(), cpval=categorize.p.value(pval)) %>% head(numbars)) %>% Reduce(rbind, .) %>% as.data.frame

 bp= ggplot(prep.df, aes(x=rank,  y=!!sym(x)) )+geom_col( alpha=0.7, aes(fill=!!sym(colorvar)))+geom_text(inherit.aes=F, aes(label=original.label, x=rank),hjust=0, y=0 )+
   #geom_text(inherit.aes=F, color="white",aes(label=hits, x=rank,hjust=0, y=0 ,vjust=-.5))+
  rotatex(90)+
   #scale_fill_viridis(option="C")+
  # scale_fill_gradient(low="#DDDDDD", high="violet")+
  #scale_fill_manual(values=allcolors$motif.clusters)+
  #scale_color_manual(values=allcolors$motif.clusters)+
  facet_grid(rows=vars(factor(clusterlabel)),cols=vars(factor(geneset, levels=genesets)), scales="free_x")+
  coord_flip()+theme_minimal()+scale_x_reverse()
 bp
  
}

make.hyper.barplot.padj=function(hdf, genesets=c("LIGPERTDOWN", "LIGPERTLINCSDN"), numbars=10, color="violet"){
prep.df=hdf %>% dplyr::filter(tlog.pval>0, geneset %in% !!genesets) %>% group_split(geneset, clusterlabel) %>% lapply(., function(x) x %>% arrange(-jaccard_index) %>% group_by(geneset) %>%  dplyr::mutate(rank=1:n(), cpval=categorize.p.value(pval)) %>% head(numbars)) %>% Reduce(rbind, .) %>% as.data.frame

 bp= ggplot(prep.df, aes(x=rank,  y=jaccard_index) )+geom_col(aes(fill=tlog.padj), alpha=0.7)+geom_text(inherit.aes=F, aes(label=short.label, x=rank,hjust=0, y=0 ))+
   #geom_text(inherit.aes=F, color="white",aes(label=hits, x=rank,hjust=0, y=0 ,vjust=-.5))+
  rotatex(90)+
   #scale_fill_viridis(option="C")+
   scale_fill_gradient(low="#DDDDDD", high=color)+
  #scale_fill_manual(values=allcolors$motif.clusters)+
  #scale_color_manual(values=allcolors$motif.clusters)+
  facet_grid(rows=vars(factor(clusterlabel)),cols=vars(factor(geneset, levels=genesets)), scales="free_x")+
  coord_flip()+theme_minimal()+scale_x_reverse()
 bp
  
}

make.hyper.barplot.pval=function(hdf, genesets=c("LIGPERTDOWN", "LIGPERTLINCSDN"), numbars=10, color="violet"){
prep.df=hdf %>% dplyr::filter(tlog.pval>0, geneset %in% !!genesets) %>% group_split(geneset, clusterlabel) %>% lapply(., function(x) x %>% arrange(-jaccard_index) %>% group_by(geneset) %>%  dplyr::mutate(rank=1:n(), cpval=categorize.p.value(pval)) %>% head(numbars)) %>% Reduce(rbind, .) %>% as.data.frame

 bp= ggplot(prep.df, aes(x=rank,  y=jaccard_index) )+geom_col(aes(fill=tlog.pval), alpha=0.7)+geom_text(inherit.aes=F, aes(label=short.label, x=rank,hjust=0, y=0 ))+
   #geom_text(inherit.aes=F, color="white",aes(label=hits, x=rank,hjust=0, y=0 ,vjust=-.5))+
  rotatex(90)+
   #scale_fill_viridis(option="C")+
   scale_fill_gradient(low="#DDDDDD", high=color)+
  #scale_fill_manual(values=allcolors$motif.clusters)+
  #scale_color_manual(values=allcolors$motif.clusters)+
  facet_grid(rows=vars(factor(clusterlabel)),cols=vars(factor(geneset, levels=genesets)), scales="free_x")+
  coord_flip()+theme_minimal()+scale_x_reverse()
 bp
  
}

make.hyper.barplot.cluster=function(hdf, genesets=c("LIGPERTDOWN", "LIGPERTLINCSDN"), numbars=10){
prep.df=hdf %>% dplyr::filter(tlog.pval>0, geneset %in% !!genesets) %>% group_split(geneset, clusterlabel) %>% lapply(., function(x) x %>% arrange(-jaccard_index) %>% group_by(geneset) %>%  dplyr::mutate(rank=1:n(), cpval=categorize.p.value(pval)) %>% head(numbars)) %>% Reduce(rbind, .) %>% as.data.frame

 bp= ggplot(prep.df, aes(x=rank,  y=jaccard_index) )+geom_col(aes(fill=clusterlabel), alpha=0.7)+geom_text(inherit.aes=F, aes(label=short.label, x=rank,hjust=0, y=0 ))+
   #geom_text(inherit.aes=F, color="white",aes(label=hits, x=rank,hjust=0, y=0 ,vjust=-.5))+
  rotatex(90)+
   #scale_fill_viridis(option="C")+
  # scale_fill_gradient(low="#DDDDDD", high="violet")+
  scale_fill_manual(values=allcolors$motif.clusters)+
  #scale_color_manual(values=allcolors$motif.clusters)+
  facet_grid(rows=vars(factor(clusterlabel)),cols=vars(factor(geneset, levels=genesets)), scales="free_x")+
  coord_flip()+theme_minimal()+scale_x_reverse()
 bp
  
}




shift.last.first <- function(vec) {
  if (length(vec) <= 1) {
    return(vec)  # No change needed for empty or single-element vectors
  }
  c(tail(vec, 1), head(vec, -1))
}

pairs.from.vectors <- function(vec1, vec2) {
  if (length(vec1) != length(vec2)) {
    stop("Both vectors must have the same length.")
  }
 lapply(1:length(vec1), function(x)
   
   c(vec1[x], vec2[x])
   
   
   )
}

# Function to check, install if necessary, and load packages
install_and_load <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    # Try installing with BiocManager first
    tryCatch({
      BiocManager::install(pkg, dependencies = TRUE, update = FALSE)
    }, error = function(e) {
      message(paste("BiocManager could not install", pkg, "- trying install.packages instead."))
      install.packages(pkg, dependencies = TRUE)
    })
    
    # Load the package after installation
    library(pkg, character.only = TRUE)
  }
}



################################################################################
# get peaks from bed in a file
################################################################################

read.bed.peaks <- function(path) {
  # Find all .bed files in the path (recursive)
  bed_files <- list.files(path, pattern = "\\.bed$", full.names = TRUE, recursive = TRUE)
  
  # Read each file into a data frame and name the list elements
  peak_list <- lapply(bed_files, function(f) read.table(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE))
  
  # Set names: just the file name without path or extension
  names(peak_list) <- tools::file_path_sans_ext(basename(bed_files))
  
  return(peak_list)
}



combine_first_non_na <- function(df, ..., new_col = "combined") {
  vars <- tidyselect::eval_select(rlang::expr(c(...)), df)

  df[[new_col]] <- apply(df[vars], 1, function(row) {
    first_non_na <- row[!is.na(row)]
    if (length(first_non_na) > 0) first_non_na[1] else NA
  })

  return(df)
}



#for(tl in 1:length(treatments.list)){
#metadata.seruggia2= metadata.seruggia2 %>% dplyr::mutate(!!(paste0("tmnt.",tl)):=paste(!!!syms(treatments.list[[tl]]$reagents), sep="-"))
#}

#metadata.seruggia2=metadata.seruggia2 %>% dplyr::mutate(tmnt.full= paste_(!!!syms(paste0("tmnt.", treatment.steps))))

#fwrite(metadata.seruggia2, file="~/metadata/samplesheet_seruggia_2024_2_conditions_annotated.csv")

################################################################################
#annotating our own sample sheet with the reagent info
################################################################################

extractrbc1 <- Vectorize(function(x){
  ifelse(grepl("RBC-step-1", x), substr(x, 12, 12), extractrbc12(x))
}, USE.NAMES = FALSE)

# Round 1, type 2
extractrbc12 <- Vectorize(function(x){
  ifelse(grepl("RBC-step1", x), substr(x, 11, 11), extractrbc21(x))
}, USE.NAMES = FALSE)

# Round 1 for second rounds
extractrbc21 <- Vectorize(function(x){
  ifelse(grepl("RBC-step-2--", x), substr(x, 13, 13), extractvemain(x, "1"))
}, USE.NAMES = FALSE)

# Round 2 for second rounds
extractrbc22 <- Vectorize(function(x){
  ifelse(grepl("RBC-step-2--", x), substr(x, 14, 14), extractvemain(x, "2"))
}, USE.NAMES = FALSE)

# Combo round
extractrbc.combo <- Vectorize(function(x){
  ifelse(grepl("RBC-step-2--", x), substr(x, 13, 14), extractvemain(x, "combo"))
}, USE.NAMES = FALSE)

# Main extractor
extractvemain <- Vectorize(function(x, mode){
  if (grepl("RBC-Ve2-", x)) {
    extr <- strsplit(x, split = '-')[[1]]
  } else {
    extr <- NA
  }

  if (is.na(extr)) {
    return(NA)
  }

  if (grepl("TGF|IL4", extr[3]) && !(mode %in% c("1"))) {
    return(paste0(extr[3:4], collapse = "-"))
  } else {
    return(extractvesub(extr[3], mode))
  }
}, vectorize.args = "x", USE.NAMES = FALSE)

# Sub extractor
extractvesub <- Vectorize(function(x, mode){
  ifelse(mode == "1", substr(x, 1, 1), 
         ifelse(mode == "2", substr(x, 2, 2), 
                x))
}, USE.NAMES = FALSE)




get.treatment=Vectorize(function(x,tmnt ){
 treatments=list(
  "1"= list(epo=1, insulin.heparin=0, hydrocortisone=0) ,
  "2"=list(epo=1, insulin.heparin=0, hydrocortisone=1) ,
  "3"=list(epo=1, insulin.heparin=1, hydrocortisone=0), 
  "4"=list(epo=1, insulin.heparin=1, hydrocortisone=1), 
  "5"=list(epo=3, insulin.heparin=0, hydrocortisone=0), 
  "6"=list(epo=3, insulin.heparin=0, hydrocortisone=1), 
  "7"=list(epo=3, insulin.heparin=1, hydrocortisone=0), 
  "8"=list(epo=3, insulin.heparin=1, hydrocortisone=1) 
 )
 code=c(A=1,B=2, C=3, D=4, E=5, F=6, G=7, H=8)
  
if(x %in% names(code)){
  out=treatments[[code[x]]][[tmnt]]
}else{
 if(x %in% names(treatments)){
 out=treatments[[x]][[tmnt]]
   
   }else{
    out=NA 
   }
  
}
  
}, USE.NAMES=F)

annotate.rbc=function(df) df %>% dplyr::mutate(rbc.round1=as.integer(extractrbc1(Experiment)), rbc.round2=extractrbc22(Experiment), rbc.combo=extractrbc.combo(Experiment))




################################################################################
#
################################################################################


vector.list.to.df= function(lst){
  
  bdf=lapply(1:length(lst), function(x){
    
  df=as.data.frame(lst[[x]]) %>% dplyr::mutate(class=names(lst)[[x]]) %>% givecolnames(.,1, nms="value")
  df$duplicated.inner=duplicated(df$value)
  df
  }) %>% bind_rows

bdf$duplicated.outer=duplicated(bdf$value)
  

    #bdf= bdf %>% dplyr::select(-duplicated.inner) %>% dplyr::select(-duplicated.outer)
    bdf.nodups=bdf %>% dplyr::filter(!duplicated.inner, !duplicated.outer)
    bdf.dups = bdf %>% dplyr::filter(duplicated.inner, duplicated.outer)
    

  
list(dups=bdf.dups, nodups=bdf.nodups %>% dplyr::select(-duplicated.inner) %>% dplyr::select(-duplicated.outer) %>% col2names(., "value"))


}

################################################################################
# subset and process a deseq2 object
################################################################################




subset.and.process.dds <- function(dds, peaks, samples = NULL, vst_name = "vst", varr="rbc.combo", meta) {
  # Remove duplicated peaks
  unique_peaks <- unique(peaks)

    # Subset metadata
  if(is.null(meta)){
  meta <- as.data.frame(colData(dds))
   
  
  }
  
  meta.f <- meta[samples,, drop=F] %>% as.data.frame %>% dplyr::filter(!is.na(!!sym(varr)))
  
  fcat("removed", nrow(meta)-nrow(meta.f), "samples with NA category")
  
  removed= setdiff(rownames(meta), rownames(meta.f))
 
  fcat("removed", paste(removed, collapse=","))
    
fcat("nrow(meta):", nrow(meta.f), "ncol(meta):",ncol(meta.f))
  
  
  
    # Subset counts
  counts_filtered <- counts(dds)[rownames(dds) %in% unique_peaks, rownames(meta.f), drop=F]
  
fcat("nrow(counts):", nrow(counts_filtered), "ncol(counts):",ncol(counts_filtered))


  # Build new DESeqDataSet
  dds_new <- DESeqDataSetFromMatrix(
    countData = counts_filtered,
    colData   = meta.f,
    design    = ~dsname + rbc.combo
  )
  
  # Estimate size factors
  dds_new <- estimateSizeFactors(dds_new)
  
  # Apply VST
  rld <- varianceStabilizingTransformation(dds_new, blind = FALSE)
  
  # Add vst assay
  assay(dds_new, vst_name) <- assay(rld)
  
  return(dds_new)
}





################################################################################
# sample pair point plot, ranking the samples based on the residuals of a linear model
################################################################################


sample.point.plot = function(s1, s2, peaks = NULL, lm_fit=NULL, sz=.2) {
  s1=make.names(s1)
  s2=make.names(s2)
  fcat(1)
  df = fixcolnames(counts)[peaks,c(s1, s2), drop = F] %>% as.data.frame 
head(df)
  # it linear model
  fcat(2)
  if(is.null(lm_fit)){
  lm_fit <- lm(as.formula(paste(s2, "~", s1)), data = df)
  }
  df$expected.linear=predict(lm_fit, data=fixcolnames(vsd)[pks, s1])
  # Calculate residuals
  fcat(3)
  df$residual <- resid(lm_fit)

  # Rank residuals by absolute value (can also use signed residuals if desired)
  fcat(4)
  df$res_rank <- rank(df$residual)

  # Normalize rank to [0,1] for color gradient
  fcat(5)
  df$res_rank_scaled <- (df$res_rank - mean(df$res_rank)) /(sd(df$res_rank)/df$expected.linear)

  fcat(6)
  plt = ggplot(df, aes(x = log(!!sym(s1)), y = log(!!sym(s2)))) +
    geom_point(aes(color = res_rank_scaled), size=sz) +
    scale_color_gradient2(low = "blue", mid = "grey", high = "red" )+
    coord_cartesian(xlim=c(0,10), ylim=c(0,10))

  fcat(7)
  print(plt)

  fcat(8)
  list(plot = plt, fit = lm_fit, df=df)
}



################################################################################
# 
################################################################################

make.class.target.list <- function(df, colname) {
  unique_vals <- unique(df[[colname]])
  result <- setNames(
    lapply(unique_vals, function(val) {
      list(classvar = colname, targets = as.character(val))
    }),
    unique_vals
  )
  return(result)
}




conditional.suffix <- function(df, search.var, target.var, pattern, suf) {
  
library(dplyr)

 
  
  df %>%
    dplyr::mutate(
      !!sym(target.var) := ifelse(
        str_detect(!!sym(search.var), pattern),
        paste0(!!sym(target.var), suf),
        !!sym(target.var)
      )
    )
}


################################################################################
# common prefix
################################################################################


common.prefix <- function(strings) {
  if (length(strings) == 0) return("")
  
  # Split all strings into characters
  split_strings <- strsplit(strings, "")
  
  # Find the minimum length among the strings
  min_length <- min(sapply(split_strings, length))
  
  prefix <- character()
  
  for (i in seq_len(min_length)) {
    # Take the ith character from all strings
    chars_at_i <- sapply(split_strings, `[[`, i)
    
    # If all characters are the same, add to prefix
    if (length(unique(chars_at_i)) == 1) {
      prefix <- c(prefix, chars_at_i[1])
    } else {
      break
    }
  }
  
  paste(prefix, collapse = "")
}

find.class.prefixes= function(mat, class.variable, prefix.variable, new.column.name){
  
classes=mat %>% pull(!!sym(class.variable)) %>% unique 

prefix.list=lapply(classes, function(x) {
  
  mat %>% dplyr::filter(class.variable ==x ) %>% pull(prefix.variable) %>% common.prefix
}) %>% givename(., classes)
fcat(prefix.list)

get.class.common.prefix=function(clss, prefix.list) prefix.list[[clss]]

mat %>% dplyr::mutate(!!sym(new.column.name):= get.class.common.prefix(!!sym(class.variable), !!prefix.list))

  
}


################################################################################
# transform vector into percentiles
################################################################################


transform.to.percentiles <- function(x) {
  # Ensure input is numeric
  if (!is.numeric(x)) stop("Input must be a numeric vector.")
  
  # Compute percentile rank for each element
  ranks <- rank(x, ties.method = "average")  # Handle ties by averaging ranks
  percentiles <- (ranks - 1) / (length(x) - 1)
  
  return(percentiles)
}


################################################################################
# theme clean minimal without lines inside the plot
################################################################################

theme_clean_minimal <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.major = element_blank(),  # remove major grid lines
      panel.grid.minor = element_blank(),  # remove minor grid lines
      panel.border = element_blank(),      # no border
      axis.line = element_blank()          # remove axis lines
    )
}



################################################################################
# english number
################################################################################

english.number<- function(number) {
  if (!requireNamespace("english", quietly = TRUE)) {
    install.packages("english")
  }
  library(english)
  
  result <- as.character(english::as.english(number))
  return(result)
}



################################################################################
# move cache with specific string directly to the current version
################################################################################
transfer.cache.direct=function(cachestring, obj=NULL, previous.version=NULL){
if(!is.null(previous.version)){  
new.string=paste_(cachestring, "version", previous.version)
  
simpleCache(new.string, assignToVar="obj", reload=T)
}  
simpleCache(cachestring %>% addversion, {obj}, assignToVar="obj", recreate=T)
obj  
}




pasteall=function(x, collapse="_"){
 paste(x, collapse=collapse) 
}

################################################################################
# generate reagent panel
################################################################################

plot_treatment_config <- function(df) {
  
  # Identify treatment columns (letters + dot + number)
  treat_cols <- grep("^[A-Za-z]+\\.[0-9]+$", names(df), value = TRUE)
  
  # Pivot longer
  df_long <- df %>%
    pivot_longer(cols = all_of(as.numeric(treat_cols)), names_to = "treatment", values_to = "value") %>%
    dplyr::mutate(value = replace_na(value, 0)) %>%   # treat NA as 0
    filter(value != 0)                         # drop zeros/NA
     
  # Plot
  ggplot(df_long, aes(x = tmnt.full, y = treatment)) +
    geom_point(aes(size = value), color = "black") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Configuration", y = "Treatment", size = "Value")
}
  

################################################################################
# sort plot paths
################################################################################


sort.plots.by.date <- function(path=config$plotpath, move = TRUE) {
  # Ensure the path exists
  if (!dir.exists(path)) {
    stop("The specified path does not exist.")
  }
  
  # Get all files in the path
  files <- list.files(path, full.names = TRUE)
  
  # Filter only files that start with YYYY-MM-DD
  pattern <- "^\\d{4}-\\d{2}-\\d{2}"
  valid_files <- files[grepl(pattern, basename(files))]
  
  if (length(valid_files) == 0) {
    message("No files with YYYY-MM-DD timestamp found.")
    return(invisible(NULL))
  }
  
  for (file in valid_files) {
    filename <- basename(file)
    
    # Extract YYYY-MM-DD
    timestamp <- sub(sprintf("(%s).*", pattern), "\\1", filename)
    
    # Extract YYYY-MM
    year_month <- substr(timestamp, 1, 7)
    
    # Create subfolder
    target_dir <- file.path(path, year_month)
    if (!dir.exists(target_dir)) {
      dir.create(target_dir, recursive = TRUE)
    }
    
    # Define target path
    target_file <- file.path(target_dir, filename)
    
    # Move or copy file
    if (move) {
      file.rename(file, target_file)
    } else {
      file.copy(file, target_file, overwrite = TRUE)
    }
  }
  
  message("Files organized successfully.")
}

################################################################################
# sea4rch patterns in the cache directory
################################################################################

grep.cache <- function(pattern, path = "~/mnt_out/rcache", ignore_case = TRUE, recursive = TRUE) {
  # Expand the path (so "~" works)
  path <- path.expand(path)
  
  # Check if the directory exists
  if (!dir.exists(path)) {
    stop("Directory does not exist: ", path)
  }
  
  # List files matching the pattern
  files <- list.files(
    path = path,
    pattern = pattern,
    ignore.case = ignore_case,
    recursive = recursive,
    full.names = TRUE
  )
  
  return(files)
}
  
  
################################################################################
# purify list from undesired intersections
################################################################################
make.unique.list=function(lst){
  attr(gplots::venn(lst), "intersections")[names(lst)]
}


find.membership <- Vectorize( function(x, lst) {
  # Check input types
  if (!is.character(x) || length(x) != 1) {
    stop("x must be a single string.")
  }
  if (is.null(names(lst))) {
    stop("lst must be a named list.")
  }

  # Find which list elements contain x
  hits <- names(lst)[sapply(lst, function(v) x %in% v)]

  # Return results
  if (length(hits) == 0) {
    return(NA_character_)
  } else if (length(hits) == 1) {
    return(hits)
  } else {
    return(paste(sort(hits), collapse = "&"))
  }
},vectorize.args="x", USE.NAMES=F)


xyplot=function(tab, x, y, colorby="MYCN_status", legend=T, legend.position=ifelse(legend, "bottom", "none"), sz=1, title=colorby){
# corr <- cor(analysis.list$gpcdf[[x]], analysis.list$gpcdf[[y]])
(ggplot(tab, aes(x=!!sym(x),y=!!sym(y) ))+
geom_point( aes(color=!!sym(colorby)), size=sz)+
scale_color_manual(values=allcolors[[colorby]])+
theme_classic()+theme(legend.position = legend.position)+
guides(color = guide_legend(nrow = 1))+
ggtitle(title)) %>% ggMarginal(., type="density", groupFill=T, groupColour=T)
}



################################################################################
# transform each column of a matrix into percentiles
################################################################################

make.pctile.columns <- function(vsd2) {
  vsdp <- lapply(1:ncol(vsd2), function(x) {
    transform.to.percentiles(vsd2[, x, drop = FALSE])
  }) %>%
    Reduce(cbind, .) %>%
    as.data.frame() %>%
    giverownames(nms = rownames(vsd2)) %>%
      givecolnames(nms=colnames(vsd2)) 
  
  return(vsdp)
}
  ##############################################################################
  # filter a value or higher
  ##############################################################################

filter.matrix.higher=function(vsdp, pthresh, negative.value=0){
vsdp %>% apply(., c(1,2), function(x) ifelse(x>=pthresh, x, negative.value))
}


################################################################################
# effect size plot on several peak groups, and treatments facetted by treatment
################################################################################

effect.size.plot=function(delist, concordance.peaks, thr=0, pval=0.05, exclude.cat="concordance", exclude="other", colorlist=allcolors){
  
  check.concordance=function(x){
ifelse(x %in% concordance.peaks$concordant.open, "concordant.open", 
  ifelse(x %in% concordance.peaks$discordant, "discordant", 
    ifelse(x %in% concordance.peaks$concordant.closed, "concordant.closed", "other")))
}
  
  
demat.concord<-lapply(1:length(delist), function(x){ 
  fcat(x)
  delist[[x]] %>% as.data.frame %>% 
    #dplyr::filter(abs(log2FoldChange)>=thr, pvalue<=pval) %>%
    names2col(., "Geneid") %>% dplyr::mutate(treatment=names(delist)[x], effect=ifelse(log2FoldChange>=1, "differentially open", "differentially closed")) %>% arrange(log2FoldChange) %>% dplyr::mutate(concordance=ifelse(abs(log2FoldChange)>=thr&pvalue<=pval, check.concordance(Geneid), "no DA")) %>% dplyr::filter(!!sym(exclude.cat)!=!!exclude)  })  %>% Reduce(rbind, .)


plt2=ggplot(demat.concord %>% dplyr::filter(!is.na(concordance)), aes( x=log2FoldChange, fill=concordance))+geom_density(alpha=0.5)+facet_wrap(as.formula(~factor(treatment)), ncol=1, scale="free_y")+theme_minimal()+scale_fill_manual(values=colorlist[["concordance"]])
plt2
}



################################################################################
#  effect size plot on several peak groups and treatments facetted by peak group
################################################################################
effect.size.plot2=function(delist, concordance.peaks, thr=0, pval=0.05, padj=Inf,exclude.cat="concordance", exclude="other", colorlist=allcolors,scale=NULL){
  
  check.concordance=function(x){
ifelse(x %in% concordance.peaks$concordant.open, "concordant.open", 
  ifelse(x %in% concordance.peaks$discordant, "discordant", 
    ifelse(x %in% concordance.peaks$concordant.closed, "concordant.closed", "other")))
}
  
  
demat.concord<-lapply(1:length(delist), function(x){ 
  fcat(x)
  delist[[x]] %>% as.data.frame %>% 
    #dplyr::filter(abs(log2FoldChange)>=thr, pvalue<=pval) %>%
    names2col(., "Geneid") %>% dplyr::mutate(treatment=names(delist)[x], effect=ifelse(log2FoldChange>=1, "differentially open", "differentially closed")) %>% arrange(log2FoldChange) %>% dplyr::mutate(concordance=ifelse(abs(log2FoldChange)>=thr&pvalue<=pval&padj<=!!padj, check.concordance(Geneid), "no DA")) %>% dplyr::filter(!!sym(exclude.cat)!=!!exclude)  })  %>% Reduce(rbind, .)


plt2=ggplot(demat.concord %>% dplyr::filter(!is.na(concordance)), aes( x=log2FoldChange, fill=treatment))+geom_density(alpha=0.5)+geom_vline(xintercept = 0)+facet_wrap(as.formula(~factor(concordance)), ncol=1, scale=scale)+theme_minimal()
plt2
}

################################################################################
#. 
###############################################################################

oneless= function(x) x[1:(length(x)-1)]


################################################################################
# annotate peak table with geneid
###############################################################################


annotate.peak.table=function(tabb, annotation.peaks, annot.columns=1:4, by="Geneid"){
  pks=tabb$Geneid
  
  left_join(annotation.peaks[pks, , drop=F] %>% select(annot.columns), tabb, by=by)
  
}





################################################################################
#
################################################################################


create.color.ramp=function(values, colors) {
  
  # create a color ramp function based on the given colors
  ramp <- colorRamp(colors)
  
  # store min/max of the values
  vmin <- min(values, na.rm = TRUE)
  vmax <- max(values, na.rm = TRUE)
  
  # return a function that converts numbers → interpolated hex colors
  function(x) {
    # rescale input to [0,1]
    t <- (x - vmin) / (vmax - vmin)
    t[t < 0] <- 0
    t[t > 1] <- 1
    
    # convert interpolated RGB to hex
    rgb(ramp(t), maxColorValue = 255)
    
    
  }
}




interpolate.colors <- function(values, colors) {
  # ensure numeric
  vals <- as.numeric(values)
  
  # get range
  rng <- range(vals, na.rm = TRUE)
  
  # interpolation function
  palette_fun <- colorRampPalette(colors)
  
  # map each unique value into 0–1 range
  scaled <- (unique(vals) - rng[1]) / diff(rng)
  scaled[scaled < 0] <- 0
  scaled[scaled > 1] <- 1

  # get interpolated colors
  interpolated <- rgb(colorRamp(colors)(scaled)/255)

  # return a named vector for convenience
  names(interpolated) <- unique(values)
  
  return(interpolated)
}



################################################################################
# PCA creation and projection functions
################################################################################


# PCA pre-processeor
#
# Generate all the necessary information to plot PCA from a DESeq2 object
# in which an assay containing a variance-stabilised matrix of counts is
# stored. Copied from DESeq2::plotPCA, but with additional ability to
# say which assay to run the PCA on.

# author Gavin Kelly
#

projectsamples=function(dds, refids=NULL, testids=NULL, assay="vst", n_top_var=Inf, label=NULL, return.reference=T){
################################################################################
# create pca of reference a
################################################################################  
metadata=as.data.frame(colData(dds))  
fcat("generating pca...")  
pca.list      <- plotPCA_vst(dds, ids=refids, metadata=metadata[refids, ], assay=vst_name, ntop=n_top_var)
refpca.df<-pca.list[[1]] 
refpca.obj<-pca.list[[2]]
fcat("running pca projection...")
#project samples onto the PCA
projected.df= scale(assay(dds[, testids], vst_name) %>% t, refpca.obj$center, refpca.obj$scale) %*% refpca.obj$rotation

#combine test and reference dataset
fcat("binding projection results...")
projected.df<- cbind(metadata[testids, ], projected.df[testids, ] )
if(return.reference){
pca.data<- bind_rows(refpca.df, projected.df)
}else{
  pca.data=projected.df
}

if(!is.null(label)){
pca.data[testids, "celltype"]=label
}

pca.data}




simplePCA=function(dds, assay="vst", n_top_var=Inf, label=NULL, return.full=F){
################################################################################
# create pca of reference a
################################################################################  
metadata=as.data.frame(colData(dds))  
fcat("generating pca...")  
pca.list      <- plotPCA_vst(dds, ids=rownames(metadata), metadata=metadata, assay=vst_name, ntop=n_top_var)
refpca.df<-pca.list[[1]] 
refpca.obj<-pca.list[[2]]


refpca.df[, "pca_id"]=label


if(return.full){
  list(df=refpca.df, pca.obj=refpca.obj)
}else{
refpca.df
}
  }




