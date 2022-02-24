#!/usr/bin/env Rscript
###
# File: merge_RTBD_txt.R
# Created: Tuesday, 6th July 2021 2:56:00 pm
# Author: Zhang Tong (zhangtong516@gmail.com)
# -----
# Last Modified: Tuesday, 6th July 2021 2:56:00 pm
# Modified By: Zhang Tong (zhangtong516@gmail.com)
# -----
# Copyright (c) 2021 GIS
# 
###
library(data.table)
library(stringr)
args = commandArgs(trailingOnly=T)

files_to_merge = args[1] ### files separated by comma (,)
sep = ":"


winSize=10

input_files = unlist(tstrsplit(files_to_merge, sep))

merged_cov_dt = rbindlist(lapply(seq(1,length(input_files)), function(i)  {
    input_file = input_files[i]
    if(input_file!=""){ 
        filename = basename(input_file)
        sampleName = unlist(strsplit(filename, "\\."))[1]
        d = fread(input_file)
        names(d) = c("gene", "length", "pos", "nain3_rt", "nain3_cov", "dmso_rt", "dmso_cov")
        d_win = d[,lapply(.SD, sum, na.rm=T), by=.(gene, length, floor(pos/winSize)*winSize),  .SDcols = c("nain3_rt", "nain3_cov", "dmso_rt", "dmso_cov")]
        # readcounts_win = readcounts[, lapply(.SD, sum, na.rm=T), by = .(gene, floor(pos/winsize)*winsize), .SDcols = colinfo[,uniqueID]]
        d_win[,sample:=sampleName]
        return(d_win)
    }
}))

names(merged_cov_dt)[3] = "pos"
rt_dt =dcast(merged_cov_dt, gene + length + pos ~ sample, value.var=c("nain3_rt", "dmso_rt"))
newColName =data.table("old"=names(rt_dt))
newColName[,c("a","b","c","x","d"):=tstrsplit(old, "_")]
newColName[,new:=ifelse(is.na(d), old, str_c(c, "__", toupper(a), "__", d))]
names(rt_dt) = newColName$new 
fwrite(rt_dt, file = "all.rt.merged.wide.csv.gz")

cov_dt =dcast(merged_cov_dt, gene + length + pos ~ sample, value.var=c("nain3_cov", "dmso_cov"))
newColName2 =data.table("old"=names(cov_dt))
newColName2[,c("a","b","c","x","d"):=tstrsplit(old, "_")]
newColName2[,new:=ifelse(is.na(d), old, str_c(c, "__", toupper(a), "__", d))]
names(cov_dt) = newColName2$new 
fwrite(cov_dt, file = "all.cov.merged.wide.csv.gz")

