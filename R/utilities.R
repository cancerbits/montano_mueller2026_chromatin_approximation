# example of an R file with a helper function

# pretty print a time difference between two proc.time() calls
time_diff <- function(start_time, end_time = NULL) {
  if (is.null(end_time)) {
    end_time <- proc.time()
  }
  dt_cpu <- lubridate::make_difftime(num = sum(end_time[c('user.self', 'sys.self')] - start_time[c('user.self', 'sys.self')]))
  dt_elapsed <- lubridate::make_difftime(num = end_time['elapsed'] - start_time['elapsed'])
  
  sprintf('Elapsed time: %1.2f %s; CPU time: %1.2f %s', 
          dt_elapsed, attr(x = dt_elapsed, which = 'units'),
          dt_cpu, attr(x = dt_cpu, which = 'units'))
}


metaDir <- function(...) {
  paste0(config$project_root, "/metadata/", ...)
}
dataDir <- function(...) {
  paste0(config$out_root, "/results/", ...)
}
resourceDir <- function(...) {
  paste0(config$resource_root, "/", ...)
}

resultsDir <- function(...) {
  paste0(config$out_root, "/", ...)
}
externalDataDir <- function(...) {
  paste0(config$project_root, "/external/", ...)
}
figuresDir <- function(...) {
  paste0(config$out_root, "/figures/", ...)
}
getExternalFile <- function(name, url) {
  #dir.create(resultsDir("downloaded_files"), showWarnings=FALSE)
  localFile <- resultsDir("/atac/downloaded_files/", name)
  if(!file.exists(localFile)) {
    if(RCurl::url.exists(url)) {
      download.file(url, destfile=localFile)
    } else {
      msgF("Warning: external resource at URL %s does not exist (anymore), using locally cached file", url)
      file.copy(externalDataDir("backup_downloaded_files/", name), localFile)
    }
  }
  return(localFile)
}
dtToDf <- function(dt, rownameCol=1) {
  df <- as.data.frame(dt)
  rownames(df) <- df[,rownameCol]
  if(is.numeric(rownameCol)) {
    df <- df[,-rownameCol,drop=FALSE] 
  }
  else {
    df <- df[,setdiff(colnames(df),rownameCol),drop=FALSE] 
  }
  df
}

gg <- function(p, name, width=NA, height=NA, scale=1, res=150, type="png", noWarn=TRUE, expand0=FALSE, addBase=FALSE) {	
  fname <- resultsDir(name, ".", gsub("cairo-","",type))		
  
  if(addBase) p <- ggAddPanelBase(p)
  if(expand0) p <- p + scale_y_continuous(expand=c(0,0))
  if(noWarn) {
    suppressWarnings(
      ggsave(filename=fname, device = type, plot=p, dpi=res, scale=scale, width=width, height=height, units="in")
    )
  }
  else {
    ggsave(filename=fname, device = type, plot=p, dpi=res, scale=scale, width=width, height=height, units="in")
  }
}
ggAddPanelBase <- function(p, expand0 = FALSE) {
  p <- p +  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)
  if(expand0) p <- p + scale_y_continuous(expand=c(0,0))
  return(p)
}	
pToSig = function(p, ns="", threshs=c(0.05, 0.01, 0.001)) {
  sapply(p, function(pval) ifelse(pval>threshs[1], ns, ifelse(pval<=threshs[3], "***", ifelse(pval<=threshs[2], "**", "*"))) )
}
repl <- function(x, a, b, by_name=F) {
  if(by_name) {
    x[names(x)==a] <- b
  } else {
    x[x==a] <- b
  }
  x
}
getColors <- function(categories, pal=NULL) {	
  n <- length(categories)
  cols <- NA
  if(!is.null(pal)) {
    cols <- RColorBrewer::brewer.pal(max(n,3),pal)
  }
  else if(n == 1) cols <- "black"
  else if(n <= 8) {
    cols <- RColorBrewer::brewer.pal(max(n,3),"Set2")
    if(n < 3) cols <- cols[1:n]
  }
  else if(n <= 9) {
    cols <- RColorBrewer::brewer.pal(n,"Set1")
  }
  else if(n <= 12) {
    cols <- RColorBrewer::brewer.pal(n,"Set3")
  }
  else cols <- rainbow(n)
  return(structure(cols,names=as.character(categories)))
}

msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M"), "\t", ...)
}
msgF <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M"), "\t", sprintf(...))
}
paste_ <- function(...) {
  paste(..., sep="_")
}
forEach <- function(...) {
  void <- sapply(...)
}
prettyIntegers <- function(x) {
  formatC(x, big.mark = ",", format = "d")
}
abscap <- function(x, cap=0.99) {
  thresh <- quantile(abs(x), cap, na.rm=T)
  i <- !is.na(x) & abs(x)>=thresh
  x[i] <- sign(x[i]) * thresh
  x
}
absmax <- function(x) {
  x[which.max(abs(x))]
}
absmin <- function(x) {
  x[which.min(abs(x))]
}

capFirst <- function(str) paste0(toupper(substr(str, 1, 1)), substr(str, 2, nchar(str)))

rblapply <- function(args, fun, id="id", ..., cores=1) {
  require(data.table)
  if(cores>1) {
    require(parallel)
    res <- parallel::mclapply(X=args, FUN=fun, ..., mc.cores=cores)
    names(res) <- names(args)
    res <- rbindlist(res, idcol=id, fill=T)
    res[,paste(id):=args[get(id)]] # args have been converted to numbers --> convert them back			
  } else {
    res <- rbindlist(sapply(X=args, FUN=fun, simplify=FALSE, ...), idcol=id, fill=T)
  }
  return(res)
} 
rblapplyDT <- function(dt, fun, idcol) {
  res <- apply(dt, 1, function(x) {
    res <- data.table(fun(as.list(x)))
    res[, paste(idcol):=x[idcol]]
    res
  })
  if(is.list(res)) res <- rbindlist(res)
  return(as.data.table(res))
}
defTheme <- function(p,flipX = FALSE, boldTitles=FALSE, withAxisLines=FALSE, withGrid = FALSE, topLegend=FALSE, noLegendTitle=FALSE, baseSize = 12, baseFamily = "") {
  library(ggplot2)
  
  ggTheme <- theme_bw(base_size = baseSize, base_family = baseFamily) %+replace% 
    theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), 
          panel.grid.major = element_blank(), axis.text=element_text(size=baseSize-1), strip.text = element_text(face=ifelse(boldTitles,"bold","plain"),size=max(baseSize-2,baseSize*0.9)), axis.title = element_text(face=ifelse(boldTitles,"bold","plain"),size=max(baseSize-2,baseSize*0.9)),
          strip.background = element_blank(), legend.key = element_blank())
  if(!withGrid) {
    ggTheme <- ggTheme %+replace% theme(panel.grid.minor = element_blank())
  }
  if(flipX) {
    ggTheme <- ggTheme %+replace% theme(axis.text.x = element_text(angle=90,hjust=1))
  }
  if(withAxisLines) {
    ggTheme <- ggTheme %+replace% theme(axis.line=element_line(), axis.line.x=element_line(), axis.line.y=element_line())
  }
  
  if(topLegend) {
    ggTheme <- ggTheme %+replace% theme(legend.position="top")
  }
  if(noLegendTitle) {
    ggTheme <- ggTheme %+replace% theme(legend.title=element_blank())
  }
  
  return(ggTheme)
}
themeBoxed <- function(...) {
  theme(panel.border=element_rect(colour="black", fill=NA), axis.line=element_blank(), ...)
}	

# extract one property or slot from each element in the list:
lget <- function(l, n, slot=F, ...) {
  sapply(l, function(x) {
    if(is.null(x)) return(NULL)
    ifelse(slot, slot(x,n), x[[n]])
  }, ...)
}

#' Create or retrieve a version-controlled simpleCache
#'
#' Extends simpleCache by definining dependencies (everything in buildEnvir).
#' A hashtag for the contents of buildEnvir will be automatically added to the 
#' cacheName. This allows for rebuilding to happen automatically if (and only
#' if) the inputs have changed.
#' 
#' Note, unlike simpleCache, this function returns the value of the cache 
#' directly instead of assigning it to a variable in the environment.
#' 
#' @param cacheName an arbitrary name for the cache 
#' @param buildEnvir a list of all variables that should be available for cache building and which will be tracked for consistency
#' @param loadEnvir see simpleCache::simpleCache (doesn't usually need to changed)
#'
#' @return the value of the serialized (or newly built) cache 
versionedCache <- function(cacheName, buildEnvir=NULL, loadEnvir=parent.frame(), ...) {
  cacheVer <- digest::digest(buildEnvir)
  
  msgF("cache %s, version=%s", cacheName, cacheVer)
  simpleCache::simpleCache(cacheName=paste0(cacheName, "___v_", cacheVer), buildEnvir=buildEnvir, assignToVariable="ret_tmp", loadEnvir=loadEnvir, ...)
  cDir <- simpleCache::getCacheDir()
  if(hasArg(cacheDir)) cDir <- cacheDir
  if(hasArg(instruction)) system(sprintf("touch %s/%s.RData", cDir, paste0(cacheName, "___v_", cacheVer)))
  
  return(get("ret_tmp", envir=loadEnvir))
}

#' Get latest version-controlled simpleCache
#'
#' Retrieve the latest version of a version-controlled simpleCache object.
#' 
#' @param cacheName an arbitrary name for the cache 
#' @param loadEnvir see simpleCache::simpleCache (doesn't usually need to changed)
#'
#' @return the value of the serialized cache 
latestCache <- function(cacheName, loadEnvir=parent.frame(), ...) {
  versionedCachePrefix <- paste0(cacheName, "___v_")
  
  availCaches <- grep(versionedCachePrefix, simpleCache::listCaches(), value=TRUE)
  availCacheFiles <- file.path(simpleCache::getCacheDir(), availCaches)
  modTimes <- sapply(availCacheFiles, file.mtime)
  
  availCaches <- availCaches[rev(order(modTimes))]
  msgF("cache %s, version=%s (latest of %d cache[s])", cacheName, gsub(versionedCachePrefix, "", availCaches[1]), length(availCaches))
  
  simpleCache::simpleCache(gsub(".RData", "", availCaches[1]), assignToVariable="ret_tmp", loadEnvir=loadEnvir, ...)
  
  return(get("ret_tmp", envir=loadEnvir))
}



