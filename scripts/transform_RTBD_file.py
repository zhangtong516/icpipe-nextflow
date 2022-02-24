#!/usr/bin/env python
# -*- coding:utf-8 -*-
###
# File: transform_RTBD_file.py
# Created: Tuesday, 6th July 2021 2:04:51 pm
# Author: Zhang Tong (zhangtong516@gmail.com)
# -----
# Last Modified: Tuesday, 6th July 2021 2:04:51 pm
# Modified By: Zhang Tong (zhangtong516@gmail.com)
# -----
# Copyright (c) 2021 GIS
# 
###

import os,sys
import statistics
## to filter away the gene with 1) mean cov <= minCov; 2) mean RT <= minRT
minCov=100
minRT=2

def read_RTBD(RTBD_file, reformat_rtbd_file, log_file, minCov=minCov, minRT=minRT):
    p1 = os.path.basename(RTBD_file)
    p2 = p1.strip().split(".")
    sampleName = p2[0]

    with open(RTBD_file, "r") as fin, open(reformat_rtbd_file, "w") as fo, open(log_file, "w") as flog:
        for line in fin.readlines():
            p = line.strip().split("\t")
            txn = p[0]
            txn_length = int(p[1])
            nain3_rt = dict()
            nain3_cov = dict()
            dmso_rt = dict()
            dmso_cov = dict()
            nonDetectedNt = 0

            for i in range(2, len(p)):
                cov_string = p[i]
                nt_count = i-2
                if cov_string == "NULL":
                    nonDetectedNt += 1
                    nain3_rt[nt_count] = 'NA'
                    nain3_cov[nt_count] = 'NA'
                    dmso_rt[nt_count] = 'NA'
                    dmso_cov[nt_count] = 'NA'
                else:
                    pp = cov_string.strip().split(",")
                    nain3_rt[nt_count] = int(pp[0])
                    nain3_cov[nt_count] = int(pp[1])
                    dmso_rt[nt_count] = int(pp[2])
                    dmso_cov[nt_count] = int(pp[3])

            valid_cov = [i for i in nain3_cov.values() if i!="NA"]
            valid_rt = [i for i in nain3_rt.values() if i!="NA"]

            gene_meanCoverage = statistics.mean(valid_cov)
            gene_meanRT = statistics.mean(valid_rt)

            if gene_meanCoverage >=minCov and gene_meanRT >= minRT:
                for pos in range(txn_length):
                    fo.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(txn, txn_length, pos+1, nain3_rt[pos], nain3_cov[pos], dmso_rt[pos], dmso_cov[pos]))
            else:
                flog.write("Transcript {} is filtered due to low coverage/RT.\n".format(txn))

RTBD_file = sys.argv[1]
reformat_rtbd_file = sys.argv[2]
log_file = "{}.log".format(reformat_rtbd_file)

read_RTBD(RTBD_file, reformat_rtbd_file, log_file)
