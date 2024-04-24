Collapse barcodes to final per-RBD/mutant phenotype scores
================
Tyler Starr
06/02/2023

- <a href="#setup" id="toc-setup">Setup</a>
- <a href="#calculate-per-variant-mean-scores-within-replicates"
  id="toc-calculate-per-variant-mean-scores-within-replicates">Calculate
  per-variant mean scores within replicates</a>
- <a href="#calculate-per-mutant-score-across-libraries"
  id="toc-calculate-per-mutant-score-across-libraries">Calculate
  per-mutant score across libraries</a>
- <a
  href="#correlations-among-backgrounds-and-to-prior-wuhan-hu-1-dms-data"
  id="toc-correlations-among-backgrounds-and-to-prior-wuhan-hu-1-dms-data">Correlations
  among backgrounds and to prior Wuhan-Hu-1 DMS data</a>
- <a href="#heatmaps" id="toc-heatmaps">Heatmaps!</a>

This notebook reads in the per-barcode titration Kds and expression
measurements from the `compute_binding_Kd` and
`compute_expression_meanF` scripts. It synthesizes these two sets of
results and calculates the final ‘mean’ phenotypes for each variant, and
generates some coverage and QC analyses.

``` r
#list of packages to install/load
packages = c("yaml","data.table","tidyverse","gridExtra","grid","bio3d","ggridges","ggrepel","GGally")
#install any packages not already installed
installed_packages <- packages %in% rownames(installed.packages())
if(any(installed_packages == F)){
  install.packages(packages[!installed_packages],
                   lib=c(paste("/uufs/chpc.utah.edu/common/home/",Sys.getenv("USER"),"/RLibs/",Sys.getenv("R_VERSION"),sep="")),
                   repos=c("http://cran.us.r-project.org"))
}
#load packages
invisible(lapply(packages, library, character.only=T))

knitr::opts_chunk$set(echo = T)
knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))


#read in config file
config <- read_yaml("config.yaml")

#make output directory
if(!file.exists(config$final_variant_scores_dir)){
  dir.create(file.path(config$final_variant_scores_dir))
}

#make pdb output directory
if(!file.exists(paste(config$final_variant_scores_dir,"/pdbs/",sep=""))){
  dir.create(file.path(paste(config$final_variant_scores_dir,"/pdbs/",sep="")))
}
```

Session info for reproducing environment:

``` r
sessionInfo()
```

    ## R version 4.1.3 (2022-03-10)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Rocky Linux 8.8 (Green Obsidian)
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /uufs/chpc.utah.edu/sys/spack/linux-rocky8-nehalem/gcc-8.5.0/intel-oneapi-mkl-2021.4.0-h43nkmwzvaltaa6ii5l7n6e7ruvjbmnv/mkl/2021.4.0/lib/intel64/libmkl_rt.so.1
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] grid      stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] GGally_2.1.2      ggrepel_0.9.1     ggridges_0.5.4    bio3d_2.4-4      
    ##  [5] gridExtra_2.3     forcats_0.5.1     stringr_1.4.0     dplyr_1.0.8      
    ##  [9] purrr_0.3.4       readr_2.1.2       tidyr_1.2.0       tibble_3.1.6     
    ## [13] ggplot2_3.4.1     tidyverse_1.3.1   data.table_1.14.2 yaml_2.3.5       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.11        lubridate_1.8.0    assertthat_0.2.1   digest_0.6.29     
    ##  [5] utf8_1.2.2         R6_2.5.1           cellranger_1.1.0   plyr_1.8.6        
    ##  [9] backports_1.4.1    reprex_2.0.1       evaluate_0.15      httr_1.4.7        
    ## [13] pillar_1.7.0       rlang_1.0.6        readxl_1.3.1       rstudioapi_0.13   
    ## [17] rmarkdown_2.13     munsell_0.5.0      broom_0.7.12       compiler_4.1.3    
    ## [21] modelr_0.1.8       xfun_0.30          pkgconfig_2.0.3    htmltools_0.5.2   
    ## [25] tidyselect_1.1.2   reshape_0.8.8      fansi_1.0.2        crayon_1.5.0      
    ## [29] tzdb_0.2.0         dbplyr_2.1.1       withr_2.5.0        jsonlite_1.8.7    
    ## [33] gtable_0.3.0       lifecycle_1.0.3    DBI_1.1.2          magrittr_2.0.2    
    ## [37] scales_1.2.1       cli_3.6.0          stringi_1.7.6      fs_1.5.2          
    ## [41] xml2_1.3.3         ellipsis_0.3.2     generics_0.1.2     vctrs_0.5.2       
    ## [45] RColorBrewer_1.1-2 tools_4.1.3        glue_1.6.2         hms_1.1.1         
    ## [49] parallel_4.1.3     fastmap_1.1.0      colorspace_2.0-3   rvest_1.0.2       
    ## [53] knitr_1.37         haven_2.4.3

## Setup

Read in tables of per-barcode expression and binding Kd measurements and
combine.

``` r
dt_bind <- data.table(read.csv(config$Titeseq_Kds_file),stringsAsFactors=F)
dt_expr <- data.table(read.csv(config$expression_sortseq_file),stringsAsFactors=F)
```

## Calculate per-variant mean scores within replicates

Calculate the mean binding and expression score collapsed by genotype.
Also output the number of barcodes across which a variant score was
determined in each library.

``` r
dt_bind[is.na(log10Ka),TiteSeq_avgcount:=NA]
dt_expr[is.na(expression),expr_count:=NA]

dt_bind[,mean_bind:=mean(log10Ka,na.rm=T),by=c("library","target","variant_class","aa_substitutions")]
dt_bind[,sd_bind:=sd(log10Ka,na.rm=T),by=c("library","target","variant_class","aa_substitutions")]
dt_bind[,n_bc_bind:=sum(!is.na(log10Ka)),by=c("library","target","variant_class","aa_substitutions")]
dt_bind[,avg_count_bind:=mean(TiteSeq_avgcount,na.rm=T),by=c("library","target","variant_class","aa_substitutions")]

dt_bind <- unique(dt_bind[,.(library,target,variant_class,aa_substitutions,n_aa_substitutions,mean_bind,sd_bind,n_bc_bind,avg_count_bind)])

dt_expr[,mean_expr:=mean(expression,na.rm=T),by=c("library","target","variant_class","aa_substitutions")]
dt_expr[,sd_expr:=sd(expression,na.rm=T),by=c("library","target","variant_class","aa_substitutions")]
dt_expr[,n_bc_expr:=sum(!is.na(expression)),by=c("library","target","variant_class","aa_substitutions")]
dt_expr[,avg_count_expr:=mean(expr_count,na.rm=T),by=c("library","target","variant_class","aa_substitutions")]

dt_expr <- unique(dt_expr[,.(library,target,variant_class,aa_substitutions,n_aa_substitutions,mean_expr,sd_expr,n_bc_expr,avg_count_expr)])
```

Some QC plots. First, look at distribution of number barcodes for
binding and expression measurements for single mutant detemrinations.
These are ‘left-justified’ histograms, so the leftmost bar represents
the number of genotypes for which no barcodes were collapsed to final
measurement in a pool.

``` r
par(mfrow=c(2,2))
hist(dt_bind[library=="pool1" & variant_class=="1 nonsynonymous",n_bc_bind],main="pool1, bind",right=F,breaks=max(dt_bind[library=="pool1" & variant_class=="1 nonsynonymous",n_bc_bind],na.rm=T),xlab="")
hist(dt_bind[library=="pool2" & variant_class=="1 nonsynonymous",n_bc_bind],main="pool2, bind",right=F,breaks=max(dt_bind[library=="pool2" & variant_class=="1 nonsynonymous",n_bc_bind],na.rm=T),xlab="")
hist(dt_expr[library=="pool1" & variant_class=="1 nonsynonymous",n_bc_expr],main="pool1, expr",right=F,breaks=max(dt_expr[library=="pool1" & variant_class=="1 nonsynonymous",n_bc_expr],na.rm=T),xlab="number barcodes collapsed")
hist(dt_expr[library=="pool2" & variant_class=="1 nonsynonymous",n_bc_expr],main="pool2, expr",right=F,breaks=max(dt_expr[library=="pool2" & variant_class=="1 nonsynonymous",n_bc_expr],na.rm=T),xlab="number barcodes collapsed")
```

<img src="collapse_scores_files/figure-gfm/hist_n_bc_per_mutant-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/histogram_n_bc_per_geno_sep-libs.pdf",sep=""),useDingbats=F))
```

What about how SEM tracks with number of barcodes collapsed? This could
help for choosing a minimum number of barcodes to use.

``` r
par(mfrow=c(2,2))
plot(dt_bind[library=="pool1" & variant_class=="1 nonsynonymous",n_bc_bind],
     dt_bind[library=="pool1" & variant_class=="1 nonsynonymous",sd_bind/sqrt(n_bc_bind)],
     pch=16,col="#00000005",main="pool1, bind",ylab="SEM",xlab="number barcodes collapsed")
plot(dt_bind[library=="pool2" & variant_class=="1 nonsynonymous",n_bc_bind],
     dt_bind[library=="pool2" & variant_class=="1 nonsynonymous",sd_bind/sqrt(n_bc_bind)],
     pch=16,col="#00000005",main="pool2, bind",ylab="SEM",xlab="number barcodes collapsed")
plot(dt_expr[library=="pool1" & variant_class=="1 nonsynonymous",n_bc_expr],
     dt_expr[library=="pool1" & variant_class=="1 nonsynonymous",sd_expr/sqrt(n_bc_expr)],
     pch=16,col="#00000005",main="pool1, expr",ylab="SEM",xlab="number barcodes collapsed")
plot(dt_expr[library=="pool2" & variant_class=="1 nonsynonymous",n_bc_expr],
     dt_expr[library=="pool2" & variant_class=="1 nonsynonymous",sd_expr/sqrt(n_bc_expr)],
     pch=16,col="#00000005",main="pool2, expr",ylab="SEM",xlab="number barcodes collapsed")
```

<img src="collapse_scores_files/figure-gfm/sem_v_n-bc-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/sem_v_n-bc.pdf",sep=""),useDingbats=F))
```

Format into a ‘mutation lookup table’, where we focus just on the single
mutants (and wildtype), breakup the string of mutations, and fill in the
table to also include any missing mutants.

``` r
dt_mutant_bind <- dt_bind[variant_class %in% "1 nonsynonymous",]

#split mutation string
#define function to apply
split_mut <- function(x){
  split <- strsplit(x,split="")[[1]]
  return(list(split[1],as.numeric(paste(split[2:(length(split)-1)],collapse="")),split[length(split)]))
}
dt_mutant_bind[,c("wildtype","position","mutant"):=split_mut(as.character(aa_substitutions)),by=aa_substitutions]

dt_mutant_bind <- dt_mutant_bind[,.(library,target,wildtype,position,mutant,mean_bind,sd_bind,n_bc_bind,avg_count_bind)]

aas <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","-")
#fill out missing values in table with a hideous loop, so the table is complete for all mutaitons (including those that are missing). If you are somebody who is reading this code, I apologize.
for(lib in c("pool1","pool2")){
  for(bg in as.character(unique(dt_mutant_bind$target))){
    for(pos in 1:max(dt_mutant_bind$position)){
      for(aa in aas){
        if(!(aa %in% as.character(dt_mutant_bind[library==lib & target==bg & position==pos,mutant]))){
          dt_mutant_bind <- rbind(dt_mutant_bind,list(lib, bg, dt_mutant_bind[target==bg & position==pos,wildtype][1],pos,aa),fill=T) #note this will leave NA for wildtype if a position is completely missing in both libraries
        }
      }
    }
  }
}
setkey(dt_mutant_bind,library,target,position,mutant)

#reset index for BA.2.86 which has site 483 deletion, or site 153 in 1-indexed numbering
dt_mutant_bind[target=="BA286" & position>=153,position:=position+1]
dt_mutant_bind[target=="BA286" & position==202,position:=153]
dt_mutant_bind[target=="BA286" & position==153,wildtype:="-"]
setkey(dt_mutant_bind,library,target,position,mutant)


#fill in wildtype values -- should vectorize in data table but being so stupid so just going to write for loop
for(bg in c("BA286","EG5","FLip")){
  for(lib in c("pool1","pool2")){
    dt_mutant_bind[library==lib & target==bg & wildtype==mutant, c("mean_bind","sd_bind","n_bc_bind","avg_count_bind"):=dt_bind[library==lib & target==bg & variant_class=="wildtype",.(mean_bind,sd_bind,n_bc_bind,avg_count_bind)]]
  }
}

#add delta bind measures
for(bg in c("BA286","EG5","FLip")){
  for(lib in c("pool1","pool2")){
    ref_bind <- dt_bind[library==lib & target==bg & variant_class=="wildtype",mean_bind]
    dt_mutant_bind[library==lib & target==bg,delta_bind := mean_bind - ref_bind]
  }
}

#repeat for expr
dt_mutant_expr <- dt_expr[variant_class %in% "1 nonsynonymous",]

#split mutation string
#define function to apply
split_mut <- function(x){
  split <- strsplit(x,split="")[[1]]
  return(list(split[1],as.numeric(paste(split[2:(length(split)-1)],collapse="")),split[length(split)]))
}
dt_mutant_expr[,c("wildtype","position","mutant"):=split_mut(as.character(aa_substitutions)),by=aa_substitutions]

dt_mutant_expr <- dt_mutant_expr[,.(library,target,wildtype,position,mutant,mean_expr,sd_expr,n_bc_expr,avg_count_expr)]

aas <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","-")
#fill out missing values in table with a hideous loop, so the table is complete for all mutaitons (including those that are missing). If you are somebody who is reading this code, I apologize.
for(lib in c("pool1","pool2")){
  for(bg in as.character(unique(dt_mutant_expr$target))){
    for(pos in 1:max(dt_mutant_expr$position)){
      for(aa in aas){
        if(!(aa %in% as.character(dt_mutant_expr[library==lib & target==bg & position==pos,mutant]))){
          dt_mutant_expr <- rbind(dt_mutant_expr,list(lib, bg, dt_mutant_expr[target==bg & position==pos,wildtype][1],pos,aa),fill=T)  #note this will leave NA for wildtype if a position is completely missing in both libraries
        }
      }
    }
  }
}
setkey(dt_mutant_expr,library,target,position,mutant)

#reset index for BA.2.86 which has site 483 deletion, or site 153 in 1-indexed numbering
dt_mutant_expr[target=="BA286" & position>=153,position:=position+1]
dt_mutant_expr[target=="BA286" & position==202,position:=153]
dt_mutant_expr[target=="BA286" & position==153,wildtype:="-"]
setkey(dt_mutant_expr,library,target,position,mutant)

#fill in wildtype values -- should vectorize in data table but being so stupid so just going to write for loop
for(bg in c("BA286","EG5","FLip")){
  for(lib in c("pool1","pool2")){
    dt_mutant_expr[library==lib & target==bg & wildtype==mutant, c("mean_expr","sd_expr","n_bc_expr","avg_count_expr"):=dt_expr[library==lib & target==bg & variant_class=="wildtype",.(mean_expr,sd_expr,n_bc_expr,avg_count_expr)]]
  }
}

#add delta expr measures
for(bg in c("BA286","EG5","FLip")){
  for(lib in c("pool1","pool2")){
    ref_expr <- dt_expr[library==lib & target==bg & variant_class=="wildtype",mean_expr]
    dt_mutant_expr[library==lib & target==bg,delta_expr := mean_expr - ref_expr]
  }
}
```

We have duplicates for each measurement. Let’s look at correlations!

``` r
par(mfrow=c(1,2))
x <- dt_mutant_expr[library=="pool1" & wildtype!=mutant,mean_expr]; y <- dt_mutant_expr[library=="pool2" & wildtype!=mutant,mean_expr]; plot(x,y,pch=16,col="#00000020",xlab="replicate 1",ylab="replicate 2",main="expression");model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

x <- dt_mutant_bind[library=="pool1" & wildtype!=mutant,mean_bind]; y <- dt_mutant_bind[library=="pool2" & wildtype!=mutant,mean_bind]; plot(x,y,pch=16,col="#00000020",xlab="replicate 1",ylab="replicate 2",main="binding affinity");model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")
```

<img src="collapse_scores_files/figure-gfm/plot_correlations-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/replicate_correlations.pdf",sep=""),useDingbats=F))
```

Binding correlations, separate per background:

``` r
par(mfrow=c(1,3))

x <- dt_mutant_bind[library=="pool1" & wildtype!=mutant & target=="BA286",mean_bind]; y <- dt_mutant_bind[library=="pool2" & wildtype!=mutant & target=="BA286",mean_bind]; plot(x,y,pch=16,col="#00000020",xlab="replicate 1",ylab="replicate 2",main="binding affinity, BA286");model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

x <- dt_mutant_bind[library=="pool1" & wildtype!=mutant & target=="EG5",mean_bind]; y <- dt_mutant_bind[library=="pool2" & wildtype!=mutant & target=="EG5",mean_bind]; plot(x,y,pch=16,col="#00000020",xlab="replicate 1",ylab="replicate 2",main="binding affinity, EG5");model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

x <- dt_mutant_bind[library=="pool1" & wildtype!=mutant & target=="FLip",mean_bind]; y <- dt_mutant_bind[library=="pool2" & wildtype!=mutant & target=="FLip",mean_bind]; plot(x,y,pch=16,col="#00000020",xlab="replicate 1",ylab="replicate 2",main="binding affinity, FLip");model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")
```

<img src="collapse_scores_files/figure-gfm/plot_correlations_separate-backgrounds-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/replicate_correlations_bind_separate-backgrounds.pdf",sep=""),useDingbats=F))
```

## Calculate per-mutant score across libraries

Collapse down to mean from both replicates, and total n barcodes between
the two/three replicates. Also record the number of the replicates the
variant was quantified within. Note, we are currently keeping a value
even if it’s determined from a single bc fit in a single pool. Later on,
we may want to require some combination of minimum number of bcs within
or between libraries for retention.

``` r
dt_final_bind <- copy(dt_mutant_bind)

dt_final_bind[ ,bind_tot:=mean(mean_bind,na.rm=T),by=c("target","position","mutant")]
dt_final_bind[ ,delta_bind_tot:=mean(delta_bind,na.rm=T),by=c("target","position","mutant")]
dt_final_bind[ ,n_bc_bind_tot:=sum(n_bc_bind,na.rm=T),by=c("target","position","mutant")]
dt_final_bind[ ,n_libs_bind_tot:=sum(!is.na(mean_bind)),by=c("target","position","mutant")]

#switch to spike indexing of postitions
dt_final_bind$position <- dt_final_bind$position + config$site_number_offset

#add single mutation string
dt_final_bind[,mutation:=paste(wildtype,position,mutant,sep=""),by=c("wildtype","position","mutant")]

dt_final_bind <- unique(dt_final_bind[,.(target,wildtype,position,mutant,mutation,bind_tot,delta_bind_tot,n_bc_bind_tot,n_libs_bind_tot)])

#repeat expr
dt_final_expr <- copy(dt_mutant_expr)

dt_final_expr[ ,expr_tot:=mean(mean_expr,na.rm=T),by=c("target","position","mutant")]
dt_final_expr[ ,delta_expr_tot:=mean(delta_expr,na.rm=T),by=c("target","position","mutant")]
dt_final_expr[ ,n_bc_expr_tot:=sum(n_bc_expr,na.rm=T),by=c("target","position","mutant")]
dt_final_expr[ ,n_libs_expr_tot:=sum(!is.na(mean_expr)),by=c("target","position","mutant")]

#switch to spike indexing of postitions
dt_final_expr$position <- dt_final_expr$position + config$site_number_offset

#add single mutation string
dt_final_expr[,mutation:=paste(wildtype,position,mutant,sep=""),by=c("wildtype","position","mutant")]

dt_final_expr <- unique(dt_final_expr[,.(target,wildtype,position,mutant,mutation,expr_tot,delta_expr_tot,n_bc_expr_tot,n_libs_expr_tot)])

#merge together
dt_final <- merge(dt_final_bind, dt_final_expr)
setkey(dt_final,target,position,mutant)

#add the rep1 and rep2 bind and expr averages
dt_final[,bind_rep1 := dt_mutant_bind[library=="pool1", mean_bind]]
dt_final[,bind_rep2 := dt_mutant_bind[library=="pool2", mean_bind]]
dt_final[,expr_rep1 := dt_mutant_expr[library=="pool1", mean_expr]]
dt_final[,expr_rep2 := dt_mutant_expr[library=="pool2", mean_expr]]


#rename some of the columns
setnames(dt_final,"bind_tot","bind")
setnames(dt_final,"delta_bind_tot","delta_bind")
setnames(dt_final,"n_bc_bind_tot","n_bc_bind")
setnames(dt_final,"n_libs_bind_tot","n_libs_bind")
setnames(dt_final,"expr_tot","expr")
setnames(dt_final,"delta_expr_tot","delta_expr")
setnames(dt_final,"n_bc_expr_tot","n_bc_expr")
setnames(dt_final,"n_libs_expr_tot","n_libs_expr")
```

Censor any measurements that are from \<3 bc or only sampled in a single
replicate? Don’t do this for now.

``` r
# min_bc <- 2
# min_lib <- 2
# 
# dt_final[n_bc_bind < min_bc & n_libs_bind < min_lib, c("bind","delta_bind","n_bc_bind","n_libs_bind") := list(NA,NA,NA,NA)]
# dt_final[n_bc_expr < min_bc & n_libs_expr < min_lib, c("expr","delta_expr","n_bc_expr","n_libs_expr") := list(NA,NA,NA,NA)]
```

Coverage stats on n_barcodes for different measurements in the final
pooled measurements

``` r
par(mfrow=c(3,2))
#BA286
hist(dt_final[wildtype!=mutant & target %in% c("BA286"), n_bc_bind],col="gray50",main=paste("mutant bind score,\nmedian ",median(dt_final[wildtype!=mutant & target %in% c("BA286"), n_bc_bind],na.rm=T),sep=""),right=F,breaks=max(dt_final[wildtype!=mutant & target %in% c("BA286","BA286"), n_bc_bind])/2,xlab="number barcodes", xlim=c(0,100))
hist(dt_final[wildtype!=mutant & target %in% c("BA286"), n_bc_expr],col="gray50",main=paste("mutant expr score,\nmedian ",median(dt_final[wildtype!=mutant & target %in% c("BA286"), n_bc_expr],na.rm=T),sep=""),right=F,breaks=max(dt_final[wildtype!=mutant & target %in% c("BA286","BA286"), n_bc_expr])/2,xlab="", xlim=c(0,100))

#EG5
hist(dt_final[wildtype!=mutant & target %in% c("EG5"), n_bc_bind],col="gray50",main=paste("mutant bind score,\nmedian ",median(dt_final[wildtype!=mutant & target %in% c("EG5"), n_bc_bind],na.rm=T),sep=""),right=F,breaks=max(dt_final[wildtype!=mutant & target %in% c("EG5","EG5"), n_bc_bind])/2,xlab="number barcodes", xlim=c(0,100))
hist(dt_final[wildtype!=mutant & target %in% c("EG5"), n_bc_expr],col="gray50",main=paste("mutant expr score,\nmedian ",median(dt_final[wildtype!=mutant & target %in% c("EG5"), n_bc_expr],na.rm=T),sep=""),right=F,breaks=max(dt_final[wildtype!=mutant & target %in% c("EG5","EG5"), n_bc_expr])/2,xlab="", xlim=c(0,100))

#FLip
hist(dt_final[wildtype!=mutant & target %in% c("FLip"), n_bc_bind],col="gray50",main=paste("mutant bind score,\nmedian ",median(dt_final[wildtype!=mutant & target %in% c("FLip"), n_bc_bind],na.rm=T),sep=""),right=F,breaks=max(dt_final[wildtype!=mutant & target %in% c("FLip","FLip"), n_bc_bind])/2,xlab="number barcodes", xlim=c(0,100))
hist(dt_final[wildtype!=mutant & target %in% c("FLip"), n_bc_expr],col="gray50",main=paste("mutant expr score,\nmedian ",median(dt_final[wildtype!=mutant & target %in% c("FLip"), n_bc_expr],na.rm=T),sep=""),right=F,breaks=max(dt_final[wildtype!=mutant & target %in% c("FLip","FLip"), n_bc_expr])/2,xlab="", xlim=c(0,100))
```

<img src="collapse_scores_files/figure-gfm/n_barcode_plots-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/histogram_n_bc_per_geno_pooled-libs.pdf",sep="")))
```

## Correlations among backgrounds and to prior Wuhan-Hu-1 DMS data

Look at correlations in mutation effects between each background, for
bind phenotype

``` r
par(mfrow=c(2,2))

x <- dt_final[target=="FLip",bind]; y <- dt_final[target=="EG5",bind]; plot(x,y,pch=16,col="#00000020",xlab="FLip",ylab="EG5",main="", xlim=c(5,14.5),ylim=c(5,14.5));model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

x <- dt_final[target=="BA286",bind]; y <- dt_final[target=="EG5",bind]; plot(x,y,pch=16,col="#00000020",xlab="BA286",ylab="EG5",main="", xlim=c(5,14.5),ylim=c(5,14.5));model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

plot(0,type='n',axes=FALSE,ann=F)

x <- dt_final[target=="BA286",bind]; y <- dt_final[target=="FLip",bind]; plot(x,y,pch=16,col="#00000020",xlab="BA286",ylab="FLip",main="", xlim=c(5,14.5),ylim=c(5,14.5));model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")
```

<img src="collapse_scores_files/figure-gfm/plot_correlations_by_bg_bind-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/background_correlations_bind.pdf",sep=""),useDingbats=F))
```

Look at correlations in mutation effects between each background, for
expr phenotype

``` r
par(mfrow=c(2,2))

x <- dt_final[target=="FLip",expr]; y <- dt_final[target=="EG5",expr]; plot(x,y,pch=16,col="#00000020",xlab="FLip",ylab="EG5",main="", xlim=c(5,11.5),ylim=c(5,11.5));model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

x <- dt_final[target=="BA286",expr]; y <- dt_final[target=="EG5",expr]; plot(x,y,pch=16,col="#00000020",xlab="BA286",ylab="EG5",main="", xlim=c(5,11.5),ylim=c(5,11.5));model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

plot(0,type='n',axes=FALSE,ann=F)

x <- dt_final[target=="BA286",expr]; y <- dt_final[target=="FLip",expr]; plot(x,y,pch=16,col="#00000020",xlab="BA286",ylab="FLip",main="", xlim=c(5,11.5),ylim=c(5,11.5));model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")
```

<img src="collapse_scores_files/figure-gfm/plot_correlations_by_bg_expr-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/background_correlations_expr.pdf",sep=""),useDingbats=F))
```

Merge with the prior VOC DMS experiments

``` r
dt_voc <- data.table(read.csv(file=config$mut_bind_expr,stringsAsFactors=F))

setkey(dt_voc,target,position,mutant)

#rename targets
dt_final[target=="BA286",target:="Omicron_BA286"]
dt_final[target=="EG5",target:="Omicron_EG5"]
dt_final[target=="FLip",target:="Omicron_FLip"]

dt_final$bind_rep3 <- as.numeric(NA)

dt_final <- rbind(dt_final[,.(target,wildtype,position,mutant,mutation,bind,delta_bind,n_bc_bind,n_libs_bind,bind_rep1,bind_rep2,bind_rep3,expr,delta_expr,n_bc_expr,n_libs_expr,expr_rep1,expr_rep2)],dt_voc[,.(target,wildtype,position,mutant,mutation,bind,delta_bind,n_bc_bind,n_libs_bind,bind_rep1,bind_rep2,bind_rep3,expr,delta_expr,n_bc_expr,n_libs_expr,expr_rep1,expr_rep2)])

setkey(dt_final,target,position,mutant)

rm(dt_voc)
```

## Heatmaps!

Order factor variables for plotting

``` r
#order targets in plotting order
dt_final$target <- factor(dt_final$target,levels=c("Wuhan-Hu-1","Eta","Alpha","Beta","Delta","Omicron_BA1","Omicron_BA2","Omicron_BA286","Omicron_BQ11","Omicron_XBB15","Omicron_EG5","Omicron_FLip"))
#order mutant as a factor for grouping by rough biochemical grouping
dt_final$mutant <- factor(dt_final$mutant, levels=c("-","C","P","G","V","M","L","I","A","F","W","Y","T","S","N","Q","E","D","H","K","R"))
#add character vector indicating wildtype to use as plotting symbols for wt
dt_final[,wildtype_indicator := ""]
dt_final[as.character(mutant)==as.character(wildtype),wildtype_indicator := "x"]


#make temp long-form data frame
temp <- data.table::melt(dt_final[target %in% c("Wuhan-Hu-1","Omicron_BA1","Omicron_BA2","Omicron_BA286","Omicron_XBB15","Omicron_EG5","Omicron_FLip"), .(target,position,mutant,bind,delta_bind,expr,delta_expr,wildtype_indicator)],id.vars=c("target","position","mutant","wildtype_indicator"),measure.vars=c("bind","delta_bind","expr","delta_expr"),variable.name="measurement",value.name="value")

#for method to duplicate aa labels on right side of plot https://github.com/tidyverse/ggplot2/issues/3171
guide_axis_label_trans <- function(label_trans = identity, ...) {
  axis_guide <- guide_axis(...)
  axis_guide$label_trans <- rlang::as_function(label_trans)
  class(axis_guide) <- c("guide_axis_trans", class(axis_guide))
  axis_guide
}

guide_train.guide_axis_trans <- function(x, ...) {
  trained <- NextMethod()
  trained$key$.label <- x$label_trans(trained$key$.label)
  trained
}
```

Make heatmaps faceted by target, showing raw affinity and delta-affinity
of muts relative to respective

``` r
p1 <- ggplot(temp[measurement=="bind",],aes(position,mutant))+geom_tile(aes(fill=value),color="black",lwd=0.1)+
  scale_fill_gradientn(colours=c("#FFFFFF","#003366","#003366"),limits=c(5,14.5),values=c(0,6.5/9.5,9.5/9.5),na.value="yellow")+
  #scale_fill_gradientn(colours=c("#FFFFFF","#FFFFFF","#003366"),limits=c(5,12),values=c(0,1/7,7/7),na.value="yellow")+ #three notches in case I want to 'censor' closer to the 5 boundary condition
  scale_x_continuous(expand=c(0,0),breaks=c(331,seq(335,530,by=5)))+
  labs(x="",y="")+theme_classic(base_size=9)+
  coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10),axis.text.y=element_text(face="bold",size=10))+
  facet_wrap(~target,ncol=1)+
  guides(y.sec=guide_axis_label_trans())+
  geom_text(aes(label=wildtype_indicator),size=2,color="gray10")+
  theme(strip.text.x = element_text(size = 18))

p1
```

<img src="collapse_scores_files/figure-gfm/heatmap_DMS_log10Ka-by-target-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/heatmap_SSM_log10Ka-by-target.pdf",sep="")))
```

Second, illustrating delta_log10Ka grouped by SSM position.

``` r
p1 <- ggplot(temp[measurement=="delta_bind",],aes(position,mutant))+geom_tile(aes(fill=value),color="black",lwd=0.1)+
  scale_fill_gradientn(colours=c("#A94E35","#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C","#383C6C"),limits=c(-5.5,4.5),values=c(0/10,1.5/10,3.5/10,5.5/10,6.5/10,7.5/10,10/10),na.value="yellow")+
  scale_x_continuous(expand=c(0,0),breaks=c(331,seq(335,530,by=5)))+
  labs(x="",y="")+theme_classic(base_size=9)+
  coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10),axis.text.y=element_text(face="bold",size=10))+
  facet_wrap(~target,ncol=1)+
  guides(y.sec=guide_axis_label_trans())+
  geom_text(aes(label=wildtype_indicator),size=2,color="gray10")+
  theme(strip.text.x = element_text(size = 18))

p1
```

<img src="collapse_scores_files/figure-gfm/heatmap_DMS_delta-log10Ka-by-target-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/heatmap_SSM_delta-log10Ka-by-target.pdf",sep="")))
```

Make heatmaps faceted by target, showing raw expression and
delta-expression of muts relative to respective wildtype

``` r
p1 <- ggplot(temp[measurement=="expr",],aes(position,mutant))+geom_tile(aes(fill=value),color="black",lwd=0.1)+
  scale_fill_gradientn(colours=c("#FFFFFF","#003366"),limits=c(5,11),na.value="yellow")+
  #scale_fill_gradientn(colours=c("#FFFFFF","#FFFFFF","#003366"),limits=c(5,11.2),values=c(0,1/7,7/7),na.value="yellow")+ #three notches in case I want to 'censor' closer to the 5 boundary condition
  scale_x_continuous(expand=c(0,0),breaks=c(331,seq(335,530,by=5)))+
  labs(x="",y="")+theme_classic(base_size=9)+
  coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10),axis.text.y=element_text(face="bold",size=10))+
  facet_wrap(~target,ncol=1)+
  guides(y.sec=guide_axis_label_trans())+
  geom_text(aes(label=wildtype_indicator),size=2,color="gray10")+
  theme(strip.text.x = element_text(size = 18))

p1
```

<img src="collapse_scores_files/figure-gfm/heatmap_DMS_expression-by-target-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/heatmap_SSM_expression-by-target.pdf",sep="")))
```

Second, illustrating delta_expression grouped by SSM position.

``` r
p1 <- ggplot(temp[measurement=="delta_expr",],aes(position,mutant))+geom_tile(aes(fill=value),color="black",lwd=0.1)+
  scale_fill_gradientn(colours=c("#A94E35","#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-3.5,1),values=c(0/4.5,1.5/4.5,2.5/4.5,3.5/4.5,4/4.5,4.5/4.5),na.value="yellow")+
  scale_x_continuous(expand=c(0,0),breaks=c(331,seq(335,530,by=5)))+
  labs(x="",y="")+theme_classic(base_size=9)+
  coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10),axis.text.y=element_text(face="bold",size=10))+
  facet_wrap(~target,ncol=1)+
  guides(y.sec=guide_axis_label_trans())+
  geom_text(aes(label=wildtype_indicator),size=2,color="gray10")+
  theme(strip.text.x = element_text(size = 18))

p1
```

<img src="collapse_scores_files/figure-gfm/heatmap_DMS_delta-expression-by-target-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/heatmap_SSM_delta-expression-by-target.pdf",sep="")))
```

Save output files.

``` r
dt_final[,.(target,wildtype,position,mutant,mutation,bind,delta_bind,n_bc_bind,n_libs_bind,bind_rep1,bind_rep2,bind_rep3,expr,delta_expr,n_bc_expr,n_libs_expr,expr_rep1,expr_rep2)] %>%
  mutate_if(is.numeric, round, digits=2) %>%
  write.csv(file=config$final_variant_scores_mut_file, row.names=F,quote=F)
```
