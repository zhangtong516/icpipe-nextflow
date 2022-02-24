#!/usr/bin/env python
# -*- coding:utf-8 -*-
###
# File: shape_to_bedGraph.py
# Created: Thursday, 10th December 2020 6:11:31 pm
# Author: Zhang Tong (zhangtong516@gmail.com)
# -----
# Last Modified: Thursday, 10th December 2020 6:11:32 pm
# Modified By: Zhang Tong (zhangtong516@gmail.com)
# -----
# Copyright (c) 2020 GIS
#
###
#!/usr/bin/env python

import sys,os
import argparse
import gzip

def parse_transcript(line):
    transcript = {}
    p = line.strip().split("\t")
    transcript.update({"transcript_name": p[0]})
    transcript_length = int(p[1])
    transcript.update({"transcript_length": transcript_length})
    transcript['reactivity'] = {}

    for i in range(transcript_length):
        pos_id = i+1
        transcript['reactivity'].update({pos_id: p[i+(len(p)-transcript_length)]})
    return transcript

def main(rt_file, bedgraph_file, name):
    with open(rt_file,"r") as fin, gzip.open(bedgraph_file, "wt") as fo:
        ## bed graph header:
        fo.write('''track type=bedGraph name=" {}" description="BedGraph format" visibility=2 color=200,100,0\n'''.format(rt_file))
        lines = fin.readlines()
        for line in lines:
            transcript = parse_transcript(line)
            transcript_name = transcript['transcript_name']
            transcript_length = transcript['transcript_length']
            for i in range(transcript_length):
                if not transcript['reactivity'][i+1] == "NULL":
                    if float(transcript['reactivity'][i+1]) >= 0:
                        fo.write("{}\t{}\t{}\t{}\n".format(
                            transcript_name, i+1, i+2, transcript['reactivity'][i+1]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Read final RT file from icShape pipeline and write BEDGRAPH file for USCSC genome browser.')
    parser.add_argument('-i', dest='rt_file', type=str, help='Input RT file')
    parser.add_argument('-o', dest='bedgraph_file', type=str, help='output bedgraph file')
    parser.add_argument('-n', dest='name',  type=str, default="name", help='track name')

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)
    main(args.rt_file, args.bedgraph_file, args.name)



