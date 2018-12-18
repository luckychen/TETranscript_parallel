#!/usr/local/bin/python2.7

'''
Created on Jan 27, 2014

@author: Ying Jin
@contact: yjin@cshl.edu
@author: Oliver Tam
@contact tam@cshl.edu
@status: 
@version: 1.5.0
'''
# python module
import sys
import os.path
import math
import operator
import argparse
import traceback

import subprocess
import multiprocessing
from time import time
from pipeproxy import proxy
from multiprocessing import Process
import datetime
import pysam
import cPickle as pickle
import gc

from TEToolkit.IO.ReadInputs import read_opts2
from TEToolkit.TEindex import *
from TEToolkit.EMAlgorithm import *
from TEToolkit.IntervalTree import *
from TEToolkit.GeneFeatures import *

import sys
import pdb


class ForkedPdb(pdb.Pdb):
    """A Pdb subclass that may be used
    from a forked multiprocessing child

    """
    def interaction(self, *args, **kwargs):
        _stdin = sys.stdin
        try:
            sys.stdin = open('/dev/stdin')
            pdb.Pdb.interaction(self, *args, **kwargs)
        finally:
            sys.stdin = _stdin
### Define parameters for program ###
def prepare_parser():
    desc = "Identifying differential transcription of gene and transposable elements."

    exmp = "Example: TEtranscripts -t RNAseq1.bam RNAseq2.bam -c CtlRNAseq1.bam CtlRNAseq.bam --GTF gene_annotation.gtf --TE TE_annotation.gtf --sortByPos --mode multi "

    parser = argparse.ArgumentParser(prog='TEtranscripts',description=desc, epilog=exmp) #'Identifying differential transcription binding/histone modification sites.')

    parser.add_argument('-t','--treatment', metavar='treatment sample', dest='tfiles',nargs='+', required=True,
                        help='Sample files in group 1 (e.g. treatment/mutant)')
    parser.add_argument('-c','--control', metavar='control sample', dest='cfiles',nargs='+', required=False,
                        help='Sample files in group 2 (e.g. control/wildtype)')
    parser.add_argument('--GTF', metavar='genic-GTF-file', dest='gtffile', type=str, required=True,
                        help='GTF file for gene annotations')
    parser.add_argument('--TE', metavar='TE-GTF-file', dest='tefile', type=str, required=True,
                        help='GTF file for transposable element annotations')
    parser.add_argument('--intron', metavar='intron-GTF-file', dest='intronfile', type=str, required=False,
                        help='GTF file for introns flanking transposable elements')
    parser.add_argument('--exonTE', metavar='exon-TE-GTF-file', dest='exontefile', type=str, required=False,
                        help='GTF file for transposable elements overlapping known exons')
    parser.add_argument('--intergenicTE', metavar='intergenic-TE-GTF-file', dest='intergenictefile', type=str, required=False,
                        help='GTF file for transposable elements between genes')
    parser.add_argument('--format', metavar='input file format', dest='format', type=str, nargs='?', default='BAM', choices=['BAM','SAM'],
                        help='Input file format: BAM or SAM. DEFAULT: BAM')
    parser.add_argument('--stranded', metavar='option', dest='stranded', nargs='?', type=str, default="yes", choices=['yes','no','reverse'],
                        help='Is this a stranded library? (yes, no, or reverse). DEFAULT: yes.')
    parser.add_argument('--mode', metavar='TE counting mode', dest='te_mode', nargs='?', type=str, const="multi", default='multi', choices=['uniq','multi'],
                        help='How to count TE: uniq (unique mappers only), or multi (distribute among all alignments).\
                        DEFAULT: multi')
    parser.add_argument('--project', metavar='name', dest='prj_name', nargs='?', default='TEtranscripts_out',
                        help='Name of this project. DEFAULT: TEtranscripts_out')
    parser.add_argument('-p', '--padj', metavar='pvalue', dest='pval', nargs='?', type=float, const=0.1, default=0.05,
                        help='FDR cutoff for significance. DEFAULT: 0.05')
    parser.add_argument('-f', '--foldchange', metavar='foldchange', dest='fc', nargs='?', type=float, const=2.0, default=1.0,
                        help='Fold-change ratio (absolute) cutoff for differential expression. DEFAULT: 1')
    parser.add_argument('--minread',metavar='min_read',dest='min_read',nargs='?',type=int,default=1,
                        help='read count cutoff. genes/TEs with reads less than the cutoff will not be considered. DEFAULT: 1')
    parser.add_argument('-n', '--norm', metavar='normalization', dest='norm', nargs='?', default='DESeq_default', choices=['DESeq_default','TC','quant'],
                        help='Normalization method : DESeq_default (DEseq default normalization method), TC (total annotated counts), quant (quantile normalization). DEFAULT: DESeq_default')
    parser.add_argument('--sortByPos', dest = 'sortByPos', action="store_true",
                        help='Alignment files are sorted by chromosome position.')
    parser.add_argument('-i', '--iteration', metavar='iteration', dest='numItr', nargs='?', type=int,  default=10,
                        help='number of iteration to run the optimization. DEFAULT: 10')
    parser.add_argument('--maxL',metavar='maxL',dest='maxL',nargs='?',type=int,default=500,help='maximum fragment length. DEFAULT:500')
    parser.add_argument('--minL',metavar='minL',dest='minL',nargs='?',type=int,default=0,help='minimum fragment length. DEFAULT:0')
    parser.add_argument('-L', '--fragmentLength', metavar='fragLength', dest='fragLength', nargs='?', type=int,  default=0,
                        help='average fragment length for single end reads. For paired-end, estimated from the input alignment file. DEFAULT: for paired-end, estimate from the input alignment file; for single-end, ignored by default.')
    parser.add_argument('--pickle', dest = 'b_pickle', action="store_true",
                        help='te Idx is pickled to be loaded in the next run')
    parser.add_argument('--origcnt', dest = 'b_origcnt', action="store_true",
                        help='count reads using original method from hammell lab')
    parser.add_argument('--pairedreadname', dest = 'b_pairedreadname', action="store_true",
                        help='read1 and read2 are encoded in readname as /1 and /2 instead of in FLAG: for TCGA mapsplice alignment files')
    parser.add_argument('--verbose', metavar='verbose', dest='verbose', type=int, nargs='?', default=2,
                        help='Set verbose level. 0: only show critical message, 1: show additional warning message, 2: show process information, 3: show debug messages. DEFAULT:2')
    parser.add_argument('--version', action='version', version='%(prog)s 1.5.0')

    return parser


class UnknownChrom(Exception):
    pass

# Reading files containing alignments (BAM of SAM)
def count_reads(samples, format, geneIdx, teIdx, intronIdx, stranded, te_mode, sortByPos, prj_name,numItr,fragLength,maxL, b_origcnt):
    cnt_ele_tbl = {}
    cnt_instance_tbl = {}
    discounted_ele_tbl = {}
    discounted_instance_tbl = {}
    libsize = []
    # check input files exist or not
    for filename in samples :
        if not os.path.isfile(filename):
            sys.stderr.write("File %s does not exist or is not a file.\n" % (filename))
            sys.exit(1)
    try:
        for filename in samples :
            (gene_counts,te_instance_counts,intron_counts) = count_transcript_abundance(filename,format,geneIdx,teIdx,intronIdx, stranded,te_mode, sortByPos, prj_name,numItr,fragLength,maxL, b_origcnt)
           
            #create dict from te count arrays 
            te_i_cnts_dict = teIdx.countByName(te_instance_counts)
            #summarize on elements
            te_ele_counts = teIdx.groupByEle(te_instance_counts)
            # save gene counts and TE counts into count table
            if intronIdx:
                intron_cnts_dict = intronIdx.countByName(intron_counts)
                intron_ele_counts = intronIdx.groupByEle(intron_counts)
                cnt_ele_tbl[filename] = dict(gene_counts.items() + te_ele_counts.items() + intron_ele_counts.items())
                cnt_instance_tbl[filename] = dict(gene_counts.items() + te_i_cnts_dict.items() + intron_cnts_dict.items())
               
                discounted_te_instance_counts = intronIdx.discountByFlankingIntrons(te_instance_counts, intron_counts)
                discounted_te_ele_counts = teIdx.groupByEle(discounted_te_instance_counts)
                discounted_instance_tbl[filename] = teIdx.countByName(discounted_te_instance_counts)
                discounted_ele_tbl[filename] = dict(discounted_te_ele_counts.items())

                total_ele_reads = sum(gene_counts.values())+sum(te_ele_counts.values())+sum(intron_ele_counts.values())
                total_instance_reads = sum(gene_counts.values())+sum(te_i_cnts_dict.values())+sum(intron_cnts_dict.values())
            else:
                cnt_ele_tbl[filename] = dict(gene_counts.items() + te_ele_counts.items())
                cnt_instance_tbl[filename] = dict(gene_counts.items() + te_i_cnts_dict.items())
                total_ele_reads = sum(gene_counts.values())+sum(te_ele_counts.values())
                total_instance_reads = sum(gene_counts.values())+sum(te_i_cnts_dict.values())
                
            #print "total ele reads: ",total_ele_reads
            #print "total instance reads: ", total_instance_reads
            libsize.append(total_ele_reads)
    except:
        sys.stderr.write("Error: %s\n" % str(sys.exc_info()[1]))
        sys.stderr.write( "[Exception type: %s, raised in %s:%d]\n" %
                          ( sys.exc_info()[1].__class__.__name__,
                            os.path.basename(traceback.extract_tb( sys.exc_info()[2] )[-1][0]),
                            traceback.extract_tb( sys.exc_info()[2] )[-1][1] ) )
        sys.exit(1)

    return (cnt_ele_tbl,cnt_instance_tbl,discounted_ele_tbl,discounted_instance_tbl,libsize)


#read assignment
def fetch_exon(chrom, st, cigar,direction,format):
    ''' fetch exon regions defined by cigar. st must be zero based
    return list of tuple of (chrom,st, end)
    '''
    #match = re.compile(r'(\d+)(\D)')
    chrom_st = st
    if format == "BAM" :
        chrom_st += 1
    exon_bound =[]

    for c,s in cigar:    #code and size
        if c==0:        #match
            if direction == 0 :
                exon_bound.append([chrom, chrom_st,chrom_st + s-1,"."])
            if direction == 1 :
                exon_bound.append([chrom, chrom_st,chrom_st + s-1,"+"])
            if direction == -1 :
                exon_bound.append([chrom, chrom_st,chrom_st + s-1,"-"])
            chrom_st += s
        elif c==1:        #insertion to ref
            continue
        elif c==2:        #deletion to ref
            chrom_st += s
        elif c==3:        #gap or intron
            chrom_st += s
        elif c==4:        #soft clipping. We do NOT include soft clip as part of exon
            chrom_st += s
        else:
            continue
    return exon_bound



def ovp_annotation(references,alignments_per_read,geneIdx,teIdx,intronIdx,stranded,format) :
    annot_gene = []
    annot_TE = []
    annot_intron = []

    #loop over every alignment
    for r1,r2 in alignments_per_read :
        TEs = []
        introns = []
        genes = []
        itv_list = []
        direction = 1
        if r1 is not None and r1.is_reverse :
            direction = -1
        if r2 is not None and not r2.is_reverse :
            direction = -1
        chrom1 = ""
        if r1 is not None :
            chrom1 = references[r1.tid]
            #sys.stderr.write("readname: %s \n" % (r1.query_name))
        chrom2 = ""
        if r2 is not None :
            chrom2 = references[r2.tid]
            #sys.stderr.write("readname: %s \n" % (r2.query_name))

        if stranded == "no":
            direction = 0
        if stranded == "reverse":
            direction = direction *(-1)

        #fetch all mapping intervals
        # intervals: genomic region corresponding to mapped portion of the read
        itv_list = []
        if r1 is not None :
            itv_list = fetch_exon(chrom1,r1.pos,r1.cigartuples,direction,format)
   
        if chrom2 != "" : #paired-end, both mates are mapped
            itv_list2 = fetch_exon(chrom2,r2.pos,r2.cigar,direction,format)
            itv_list.extend(itv_list2)
        try:
            #print "itv_list ",itv_list
            TEs = teIdx.TE_annotation(itv_list)
            #print "TEs: ",TEs
            if intronIdx:
                introns = intronIdx.intron_annotation(itv_list)
                #print "introns: ",introns
            genes = geneIdx.Gene_annotation(itv_list)
            #print "genes: ",genes

            if len(TEs) > 0   :
                annot_TE.append(TEs)  # annot_TE: 2-dim array; 1-dim=alignment 2-dim=overlapping_te_idx for each interval
                # changed into a 3-dim array with additional dim of [te_idx, ovp_len]

            if len(introns) > 0   :
                annot_intron.append(introns)

            if len(genes) > 0 :     
                #annot_gene.append(list(set(genes)))
                annot_gene.append(genes)
        except:
            sys.stderr.write("Error occurred during read assignments\n")
            raise

    return (annot_gene,annot_TE,annot_intron)


def readInAlignment(filename, format, sortByPos,prj_name):
    try:
        if format == "BAM" :
            if not sortByPos  :
                samfile = pysam.AlignmentFile(filename,'rb')
                if len(samfile.header) ==0:
                    print >>sys.stderr, "BAM/SAM file has no header section. Exit!"
                    sys.exit(1)

            else :
                #psort = subprocess.Popen(["samtools", "sort", "-n", "-o", filename, "%s_%s_tmp" % (filename, prj_name)], stdout=subprocess.PIPE)
                #samtools_in = subprocess.Popen(["samtools", "view", "-"], stdin=psort.stdout, stdout=subprocess.PIPE)
                #psort.stdout.close()  # Allow p to receive a SIGPIPE if psort exits.

                cur_time = time.time()
                bam_tmpfile = '.'+str(cur_time)
                sys.stderr.write("INFO  @ %s: \n" % (datetime.datetime.now()))
                sys.stderr.write("Sorting bam file ......\n")
                #pysam.sort("-n",filename,bam_tmpfile)
                sortcommand = "sambamba sort -n %s -o %s.bam 2>> sambamba.err" % (filename, bam_tmpfile)
                sortproc = subprocess.Popen(sortcommand, stdout=subprocess.PIPE, shell=True)
                sortproc.wait()
                sys.stderr.write("INFO  @ %s: \n" % (datetime.datetime.now()))
                sys.stderr.write("Done sorting bam file ......\n")
                #samfile = pysam.AlignmentFile(samtools_in.stdout,'r')
                samfile = pysam.AlignmentFile(bam_tmpfile+".bam",'r')
                if len(samfile.header) ==0:
                    print >>sys.stderr, "BAM/SAM file has no header section. Exit!"
                    sys.exit(1)
                #subprocess.call(['rm -f ' + bam_tmpfile+'.bam' ],shell=True)

        else :
                samfile = pysam.AlignmentFile(filename,'r')
                if len(samfile.header) ==0:
                    print >>sys.stderr, "BAM/SAM file has no header section. Exit!"
                    sys.exit(1)

    except:
        sys.stderr.write("Error occured when reading first line of sample file %s.\n" % filename)
        raise
    return samfile


def count_transcript_abundance(filename, format, geneIdx, teIdx, intronIdx, stranded, te_mode, sortByPos, prj_name,numItr,fragLength,maxL,b_origcnt):
    #read in BAM/SAM file
    samfile = readInAlignment(filename, format, sortByPos,prj_name)

    references = samfile.references

    empty = 0
    nonunique = 0
    uniq_reads = 0

    i = 0
    prev_read_name =''

    alignments_per_read = []
    leftOver_gene = []
    leftOver_te = []
    leftOver_intron = []
    avgReadLength =  0
    tmp_cnt = 0
    multi_read1 = []
    multi_read2 = []
    paired = False
    gene_counts = dict(zip(geneIdx.getFeatures(),[0.0]*len(geneIdx.getFeatures())))
    te_counts = [0.0] * teIdx.numInstances()
    te_multi_counts = [0.0]*len(te_counts)
    intron_counts = None
    intron_multi_counts = None
    multi_reads_te = []
    multi_reads_intron = []
    cc = 0
    if intronIdx: 
        intron_counts = [0.0] * intronIdx.numInstances()
        intron_multi_counts = [0.0]*len(intron_counts)

    try:
        while(1):
            i += 1
            cc += 1
            aligned_read = samfile.next()

            if aligned_read.is_unmapped or aligned_read.is_duplicate or aligned_read.is_qcfail :
                continue

            cur_read_name = aligned_read.query_name

            if aligned_read.is_paired :  # this is not reliable if read mate is unmapped
                added_flag_pos = aligned_read.query_name.find('/')
                if added_flag_pos == -1 :
                    added_flag_pos = len(aligned_read.query_name)
                if aligned_read.query_name.endswith(".1") or aligned_read.query_name.endswith(".2") :
                    added_flag_pos = len(aligned_read.query_name) - 2 
    
                cur_read_name = aligned_read.query_name[:added_flag_pos]
 
                paired = True
                #sys.stderr.write("cur_read_name %s == " % (cur_read_name))
                #sys.stderr.write("prev_read_name %s \n" % (prev_read_name))
                if cur_read_name == prev_read_name or prev_read_name == "":
                    prev_read_name = cur_read_name
                    #if aligned_read.is_read1 & (aligned_read.query_name[added_flag_pos+1]=="2"):
                    #    print "WRONG is_read1", aligned_read
                    #if aligned_read.is_read2 & (aligned_read.query_name[added_flag_pos+1]=="1"):
                    #    print "WRONG is_read2", aligned_read
                    if ('b_pairedreadname' in vars()): 
                        if aligned_read.is_read1 :
                                multi_read1.append(aligned_read)
                        if aligned_read.is_read2 :
                                multi_read2.append(aligned_read)
                    else:
######## pair info in read name #####################
                        if (aligned_read.query_name[added_flag_pos+1]=="1"):
                                multi_read1.append(aligned_read)
                        if (aligned_read.query_name[added_flag_pos+1]=="2"):
                                multi_read2.append(aligned_read)
################################################
                    continue

                else : # sort out all reads collected with same name so far
                    if len(multi_read1) <= 1 and len(multi_read2) <= 1 :
                        #sys.stderr.write("unique read\n")
                        uniq_reads += 1
                        read1 = None
                        read2 = None
                        if len(multi_read1) == 1:
                            read1 = multi_read1[0]
                        if len(multi_read2) == 1:
                            read2 = multi_read2[0]
                        if read1 is not None and read1.is_proper_pair :
                            #sys.stderr.write("read1 %s " % (read1.query_name))
                            if read2 is None :
                                sys.stderr.write("******NOT COMPLETE*******\n")
                                sys.stderr.write("If the BAM file is sorted by coordinates, please specify --sortByPos and re-run!\n")
                                #sys.exit(0)
                                continue
                            if tmp_cnt < 10000 :
                                #sys.stderr.write("read2 %s\n" % (read2.query_name))
                                pos1 = read1.reference_start
                                pos2 = read2.reference_start
                                if abs(pos1-pos2) <= maxL :
                                    avgReadLength += abs(pos1-pos2)+read2.query_length
                                    tmp_cnt += 1
                        alignments_per_read.append((read1,read2)) # append paired read to alignments_per_read
                    else :
                        #sys.stderr.write("non-unique\n") # multiple reads next to each other that are more than one pair 
                        nonunique += 1
                        if te_mode == 'uniq' :
                            empty += 1
                            alignments_per_read = []
                            multi_read1 = []
                            multi_read2 = []
                            prev_read_name = cur_read_name
                            if ('b_pairedreadname' in vars()): 
                                if aligned_read.is_read1 :
                                    multi_read1.append(aligned_read)
                                if aligned_read.is_read2 :
                                    multi_read2.append(aligned_read)
                            else:
######## pair info in read name #####################
                                if (aligned_read.query_name[added_flag_pos+1]=="1"):
                                    multi_read1.append(aligned_read)
                                if (aligned_read.query_name[added_flag_pos+1]=="2"):
                                    multi_read2.append(aligned_read)
#####################################################
                            continue

                        else:
                            if (b_origcnt): 
                                #singlton
                                if len(multi_read2) == 0 : 
                                    for r in multi_read1 :
                                        alignments_per_read.append((r,None))
                                if len(multi_read1) == 0:
                                    for r in multi_read2 :
                                        alignments_per_read.append((None,r))
                                if len(multi_read2) == len(multi_read1) :
                                  for i in range(len(multi_read1)) :
                                    read1 = multi_read1[i]
                                    read2 = multi_read2[i]
                                    alignments_per_read.append((read1,read2))
                            else: 
######## new counting ###############################
                                #sys.stderr.write("save alignments per read\n")
                                #print "multi_read1: ", len(multi_read1)
                                #print "multi_read2: ", len(multi_read2)
                              
                                larger_cnt =  max(len(multi_read1), len(multi_read2))
                                #print larger_cnt
                                for c in range(larger_cnt): 
                                    #singlton
                                    if len(multi_read2) == 0 :
                                        r1 = multi_read1.pop(0)
                                        alignments_per_read.append((r1,None))
                                        #print "r1 appended"
                                    elif len(multi_read1) == 0:
                                        r2 = multi_read2.pop(0)
                                        alignments_per_read.append((None, r2))
                                        #print "r2 appended"
                                    else: 
                                        r1 = multi_read1.pop(0)
                                        r2 = multi_read2.pop(0)
                                        alignments_per_read.append((r1,r2))
                                        #print "r1 r2 appended"

                                #print "SAVED ALL COLLECTED READS"
                                #for r1,r2 in alignments_per_read :
                                    #print "r1: ", r1
                                    #print "r2: ", r2
###################################################

            else : #single end read
                    if cur_read_name == prev_read_name or prev_read_name == "":
                        alignments_per_read.append((aligned_read,None))
                        prev_read_name = cur_read_name
                        continue
                    else : # a new read
                        if tmp_cnt < 10000 :
                            avgReadLength += aligned_read.query_length
                            tmp_cnt += 1

                        if len(alignments_per_read) == 1 :
                            uniq_reads += 1
                        else :
                            #sys.stderr.write("non-unique\n")
                            nonunique += 1
                            if te_mode == 'uniq' : #ignore multi-reads
                                empty += 1
                                alignments_per_read = []
                                prev_read_name = cur_read_name
                                alignments_per_read.append((aligned_read,None))
                                continue
            try :

                (annot_gene,annot_TE,annot_intron) = ovp_annotation(references,alignments_per_read, geneIdx, teIdx, intronIdx, stranded,format)

                if len(alignments_per_read) > 1 : #multi read, priority to TE
                    #print "multiple alignments per read" 
                    #for r1,r2 in alignments_per_read :
                        #print "r1: ", r1
                        #print "r2: ", r2
                    no_annot_te = parse_annotations_TE(multi_reads_te,annot_TE, te_counts, te_multi_counts, leftOver_te)
                    if intronIdx: 
                        no_annot_intron = parse_annotations_intron(multi_reads_intron, annot_intron, intron_counts, intron_multi_counts, leftOver_intron)
                        if no_annot_te and no_annot_intron:
                            no_annot_gene = parse_annotations_gene(annot_gene,gene_counts,leftOver_gene)
                            if no_annot_gene :
                                #sys.stderr.write("unannotated\n")
                                empty += 1
                    else:
                        if no_annot_te:
                            no_annot_gene = parse_annotations_gene(annot_gene,gene_counts,leftOver_gene)
                            if no_annot_gene :
                                #sys.stderr.write("unannotated\n")
                                empty += 1
 

                else : #uniq read , priority to gene
                    no_annot_gene = parse_annotations_gene(annot_gene,gene_counts,leftOver_gene)
                    if no_annot_gene :
                        no_annot_te = parse_annotations_TE(multi_reads_te,annot_TE, te_counts, te_multi_counts, leftOver_te)
                        if intronIdx:
                            no_annot_intron = parse_annotations_intron(multi_reads_intron, annot_intron, intron_counts, intron_multi_counts, leftOver_intron)
                            if no_annot_te and no_annot_intron:
                                    #sys.stderr.write("unannotated\n")
                                    empty += 1
                        else:
                            #sys.stderr.write("unannotated\n")
                            empty += 1



            except:
                sys.stderr.write("Error occurred when processing annotations of %s in file %s.\n" % (prev_read_name, filename))
                raise

            if i % 1000000 == 0 :
                sys.stderr.write("%d %s processed.\n" % (i, "alignments " ))

            alignments_per_read = []
            multi_read1 = []
            multi_read2 = []
            prev_read_name = cur_read_name
            if not aligned_read.is_paired :
                alignments_per_read.append((aligned_read,None))
            else :
                if ('b_pairedreadname' in vars()): 
                    if aligned_read.is_read1 :
                        multi_read1.append(aligned_read)
                    if aligned_read.is_read2 :
                        multi_read2.append(aligned_read)
                else:
######## pair info in read name #####################
                    if (aligned_read.query_name[added_flag_pos+1]=="1"):
                        multi_read1.append(aligned_read)
                    if (aligned_read.query_name[added_flag_pos+1]=="2"):
                        multi_read2.append(aligned_read)
#####################################################

    except StopIteration:
            pass
                #the last read
    try:
            #resolve ambiguity
            if len(leftOver_gene) > 0 :
                resolve_annotation_ambiguity(gene_counts,leftOver_gene)
            if intronIdx and len(leftOver_intron) > 0 :
                resolve_annotation_ambiguity(intron_counts,leftOver_intron)
            if len(leftOver_te) > 0 :
                resolve_annotation_ambiguity(te_counts,leftOver_te)

            ss = sum(te_counts)
            sys.stderr.write("uniq te counts = %s \n" % (str(ss)))

            # handling multi-mapped te counts
            te_tmp_counts = [0]*len(te_counts)

            if numItr > 0 :
              try :
                ''' iterative optimization on TE reads '''
                sys.stderr.write(".......start iterative optimization ..........\n")
                if not paired and fragLength > 0:
                    avgReadLength = fragLength

                elif avgReadLength > 0 :
                    avgReadLength = int(avgReadLength/tmp_cnt)
                else :
                    sys.stderr.write("There are not enough reads to estimate fragment length. \n" )
                    raise 
                new_te_multi_counts = [0] *len(te_counts)  
                if len(multi_reads_te) > 0 :
                    new_te_multi_counts = EMestimate(teIdx,multi_reads_te,te_tmp_counts,te_multi_counts,numItr,avgReadLength)

              except :
                sys.stderr.write("Error in optimization\n")
                raise
              te_counts = map(operator.add,te_counts,new_te_multi_counts)
            else :
              te_counts = map(operator.add,te_counts,te_multi_counts)
        
            if intronIdx:
                ss = sum(intron_counts)
                sys.stderr.write("uniq intron counts = %s \n" % (str(ss)))

                # handling multi-mapped intron counts
                intron_tmp_counts = [0]*len(intron_counts)

                if numItr > 0 :
                  try :
                    ''' iterative optimization on TE reads '''
                    sys.stderr.write(".......start iterative optimization ..........\n")
                    if not paired and fragLength > 0:
                        avgReadLength = fragLength

                    elif avgReadLength > 0 :
                        avgReadLength = int(avgReadLength/tmp_cnt)
                    else :
                        sys.stderr.write("There are not enough reads to estimate fragment length. \n" )
                        raise 
                    new_intron_multi_counts = [0] *len(intron_counts)  
                    if len(multi_reads_intron) > 0 :
                        new_intron_multi_counts = EMestimate(intronIdx,multi_reads_intron,intron_tmp_counts,intron_multi_counts,numItr,avgReadLength)

                  except :
                    sys.stderr.write("Error in optimization\n")
                    raise
                  intron_counts = map(operator.add,intron_counts,new_intron_multi_counts)
                else :
                  intron_counts = map(operator.add,intron_counts,intron_multi_counts)


    except:
            sys.stderr.write("Error occurred when assigning read gene/TE annotations in file %s.\n" % (filename))
            raise
    st = sum(te_counts)
    si = 0
    if intronIdx: si = sum(intron_counts)
    sg = sum(gene_counts.values())
    num_reads = st + si + sg

    sys.stderr.write("TE counts total %s\n" % (st))
    if intronIdx: sys.stderr.write("intron counts total %s\n" % (si))
    sys.stderr.write("Gene counts total %s\n" % (sg))
    sys.stderr.write("\nIn library %s:\n" % (filename))
    sys.stderr.write("Total annotated reads = %s\n" % ( str(num_reads) ) )
    sys.stderr.write("Total non-uniquely mapped reads = %s\n" % (str(int(nonunique))))
    sys.stderr.write("Total unannotated reads = %s\n\n" %(str(int(empty))))

    return (gene_counts,te_counts,intron_counts)



def parse_annotations_gene(annot_gene, gene_counts,leftOver_gene):
    #print "parsing gene annotations "
    no_annot = True

    #print "len(annot_gene) ", len(annot_gene)
    if len(annot_gene) > 1 :
        leftOver_gene.append((annot_gene,1.0))
            # if te_mode == "uniq" :
            #return no_annot
    elif len(annot_gene) == 1 :
        genes = annot_gene[0]
        #print "genes: ",genes
        if len(genes) == 1 :
            gene_counts[genes[0][0]] += 1
        else :
            if genes[0][0] == genes[1][0] :
                gene_counts[genes[0][0]] += 1
            else :
                gene_counts[genes[0][0]] += 0.5
                gene_counts[genes[1][0]] += 0.5
    else :
        return True

    return False

# Assign ambiguous genic reads mapped to a location with multiple annotations
def resolve_annotation_ambiguity(counts, leftOvers ) : # w is weight?
    for (annlist,w) in leftOvers :
        readslist = {}
        total = 0.0
        size = len(annlist)
        ww = 1.0 * w
        if size > 1:
            ww = ww/size

        #print annlist
        for ann in annlist :
            #print ann
            for (a,ovp_len) in ann:
                if a not in readslist :
                    readslist[a] = counts[a]
                    total += counts[a]

        if total > 0.0 :
            for a in readslist :
                v = ww * readslist[a] / total
                counts[a] += v
        else :
            for a in readslist :
                counts[a] = ww/len(readslist)


def parse_annotations_TE(multi_reads, annot_TE, uniq_counts, multi_counts, leftOver_list):
    #print "parsing TE annotations "
    if len(annot_TE) == 0 :
        return True
    
    #print "len(annot_TE) ", len(annot_TE)
    #print "len(annot_TE[0]) ", len(annot_TE[0])


    if len(annot_TE) == 1 and len(annot_TE[0]) == 1 :
        #print annot_TE[0][0]
        uniq_counts[annot_TE[0][0][0]] += 1

    if len(annot_TE) == 1 and len(annot_TE[0]) > 1 :
        #for i in range(len(annot_TE[0])):
            #print annot_TE[0][i]
        leftOver_list.append((annot_TE,1.0))

    if len(annot_TE) > 1:
            multi_algn = []
            for i in range(len(annot_TE)):
                for (te_idx, ovp_len) in annot_TE[i] :
                    #print te_idx, " ", ovp_len
                    multi_counts[te_idx] += 1.0 / (len(annot_TE) * len(annot_TE[i]))
                    multi_algn.append(te_idx)

            multi_reads.append(multi_algn)


    return False


# save intron counts for uniq_counts, multi_count and leftOvers
# later discount the TE counts by the read depth ratios of intron counts vs TE counts.
def parse_annotations_intron(multi_reads,annot_intron, uniq_counts,multi_counts, leftOver_list):
    #print "parsing intron annotations "
    if len(annot_intron) == 0 :
        return True

    #print "len(annot_intron) ", len(annot_intron)
    #print "len(annot_intron[0]) ", len(annot_intron[0])

    if len(annot_intron) == 1 and len(annot_intron[0]) == 1 :
        #print annot_intron[0][0]
        uniq_counts[annot_intron[0][0][0]] += 1

    if len(annot_intron) == 1 and len(annot_intron[0]) > 1 :
        #for i in range(len(annot_intron[0])):
            #print annot_intron[0][i]
        leftOver_list.append((annot_intron,1.0))

    if len(annot_intron) > 1:
            multi_algn = []
            for i in range(len(annot_intron)):
                for (intron_idx, ovp_len) in annot_intron[i] :
                    #print intron_idx, " ", ovp_len
                    multi_counts[intron_idx] += 1.0 / (len(annot_intron) * len(annot_intron[i]))
                    multi_algn.append(intron_idx)

            multi_reads.append(multi_algn)

    return False




def output_res(res, ann, smps, prj):
    fname = prj+".png"

    plotHeatmap(res, ann, smps, fname)
    return

def output_count_tbl(t_tbl, c_tbl, fname):
    try:
        f = open(fname, 'w')
    except IOError:
        error("Cannot create report file %s !\n" % (fname))
        sys.exit(1)
    else:
        cnt_tbl = {}
        header = "gene/TE"
        keys = set([])                              # keys: set of all elements or instances 
        for tsmp in t_tbl.keys():                   # tsmp: treatment sample filenames
            keys = keys.union(t_tbl[tsmp].keys())
            header +="\t"+tsmp+".T"
        if c_tbl:
            for csmp in c_tbl.keys():                   # csmp: treatment sample filenames
                keys = keys.union(c_tbl[csmp].keys())
                header +="\t"+csmp+".C"

        for tsmp in t_tbl.keys():
            cnts = t_tbl[tsmp]
            for k in keys:
                val = 0
                if k in cnts :
                    val = cnts[k]
                    if val is None or val is "":
                        sys.exit(1)
                if cnt_tbl.has_key(k) :
                    cnt_tbl[k].append(int(round(val)))
                else :
                    vallist = []
                    vallist.append(int(round(val)))
                    cnt_tbl[k] = vallist

        if c_tbl:
            for csmp in c_tbl.keys():
                cnts = c_tbl[csmp]
                for k in keys:
                    val = 0
                    if k in cnts :
                        val = cnts[k]
                        if val is None or val is "":
                            sys.exit(1)
                    if cnt_tbl.has_key(k) :
                        cnt_tbl[k].append(int(round(val)))
                    else :
                        print "we should never come here, right?"
                        sys.exit(1)
                        vallist = []
                        vallist.append(int(round(val)))
                        cnt_tbl[k] = vallist
        #output
        f.write(header + "\n")
        for gene in sorted(cnt_tbl.keys()) :
            vals = cnt_tbl[gene]
            vals_str = gene if gene else "None"
            for i in range(len(vals)) :
                vals_str +="\t"+str(vals[i])
            f.write(vals_str + "\n")

        f.close()

    return

def output_norm(sf, name, error):
    fname = name + ".norm"
    try:
        f = open(fname, 'w')
    except IOError:
        error("Cannot create report file %s !\n" % (fname))
        sys.exit(1)
    else:
        cnt = 1
        for b in sf:
            desc = "treat" + str(cnt)
            for i in range(len(b)):
                desc += "\t"+str(round(b[i], 2))
            f.write(desc + "\n")
            cnt +=1
        f.close()


def write_R_code(f_cnt_tbl, tfiles, cfiles, prj_name, norm, pval, fc, rpm_val,min_read):
    # Assembling R-code for analysis
    rscript = ''
    rscript += '\n'
    rscript += 'data <- read.table("%s",header=T,row.names=1)\n' % (f_cnt_tbl) # load counts table
    rscript += 'groups <- factor(c(rep("TGroup",%s),rep("CGroup",%s)))\n' % (len(tfiles),len(cfiles)) # generate groups for pairwise comparison
    rscript += 'min_read <- %s\n' % (min_read)
    # Counts filtering (hard coded to 20)
    rscript += 'data <- data[apply(data,1,function(x){max(x)}) > min_read,]\n'
    # Quantile normalization to calculate fold change
    if norm == 'quant':
        rscript += 'colnum <- length(data)\n'
        rscript += 'rownum <- length(data[,1])\n'
        rscript += 'ordMatrix <- matrix(nrow=rownum,ncol=colnum)\n'
        rscript += 'ordIdx <- matrix(nrow=rownum,ncol=colnum)\n'
        rscript += 'for (i in 1:colnum){\n'
        rscript += '  a.sort <- sort(data[,i],index.return=T)\n'
        rscript += '  ordMatrix[,i] <- a.sort$x\n'
        rscript += '  ordIdx[,i] <- a.sort$ix\n'
        rscript += '}\n'
        rscript += 'rowAvg <- rowMeans(ordMatrix)\n'
        rscript += 'data.q.norm <- matrix(nrow=rownum,ncol=colnum)\n'
        rscript += 'for (i in 1:colnum){\n'
        rscript += '  data.q.norm[,i] <- rowAvg[order(ordIdx[,i])]\n'
        rscript += '}\n'
        rscript += 'colnames(data.q.norm) <- colnames(data)\n'
        rscript += 'rownames(data.q.norm) <- rownames(data)\n'
        if len(tfiles) > 1:
            rscript += 'sample1Mean <- rowMeans(data.q.norm[,1:%s],na.rm=T)\n' % (len(tfiles))
        else:
            rscript += 'sample1Mean <- data.q.norm[,1]\n'
        group2_start = len(tfiles) + 1
        group2_stop = group2_start + len(cfiles)
        if len(cfiles) > 1:
            rscript += 'sample2Mean <- rowMeans(data.q.norm[,%s:%s,na.rm=T)\n' % (group2_start, group2_stop)
        else:
            rscript += 'sample2Mean <- data.q.norm[,%s]\n' % (group2_start)
        rscript += 'FoldChange <- (sample2Mean/sample1Mean)\n'
        rscript += 'log2FoldChange <- log2(FoldChange)\n'

    # Normalize by RPM (reads per million mapped)
    if norm == 'TC'  :
        min_libSize = min(rpm_val)
        rpm_vec = ','.join(str(x/min_libSize) for x in rpm_val)
        rscript += 'tc <- c(%s)\n' % (rpm_vec)


    # Performing differential analysis using DESeq
    rscript += 'library(DESeq, quietly=T)\n'
    rscript += 'cds <- newCountDataSet(data,groups)\n'
    if norm == 'TC':
        rscript += 'cds$sizeFactor = tc\n'
    else:
        rscript += 'cds <- estimateSizeFactors(cds)\n'
    if(len(tfiles)==1 and len(cfiles) ==1):
        rscript += 'cds <- estimateDispersions(cds,method="blind",sharingMode="fit-only",fitType="local")\n'
    elif(len(tfiles) > 1 and len(cfiles) > 1):
        rscript += 'cds <- estimateDispersions(cds,method="per-condition")\n'
    else :
        rscript += 'cds <- estimateDispersions(cds,method="pooled")\n'

    rscript += 'res <- nbinomTest(cds,"CGroup","TGroup")\n'

    # Generating output table
    if norm == 'quant':
        rscript += 'res_fc <- cbind(res$id,sample1Mean,sample2Mean,FoldChange,log2FoldChange,res$pval,res$padj)\n'
        rscript += 'colnames(res_fc) = c("id","sample1Mean","sample2Mean","FoldChange","log2FoldChange","pval", "padj")\n'
    else:
        rscript += 'res_fc <- res\n'
    rscript += 'write.table(res_fc, file="%s_gene_TE_analysis.txt", sep="\\t",quote=F,row.names=F)\n' % (prj_name)

    # Generating table of "significant" results

    l2fc = math.log(fc,2)
    if norm == 'quant':
        rscript += 'resSig <- res_fc[(!is.na(res_fc[,7]) & (res_fc[,7] < %f) & (abs(as.numeric(res_fc[,5])) > %f)), ]\n' % (pval, l2fc)
    else:
        rscript += 'resSig <- res_fc[(!is.na(res_fc$padj) & (res_fc$padj < %f) & (abs(res_fc$log2FoldChange)> %f)), ]\n' % (pval, l2fc)
    rscript += 'write.table(resSig, file="%s_sigdiff_gene_TE.txt",sep="\\t", quote=F, row.names=F)\n' % (prj_name)

    return rscript

class GetSetter(object):
    def __init__(self):
        self.var = None

    def set(self, value):
        self.var = value

    def get(self):
        return self.var


def subprocessWorker():
    # Read sample files make count table
    global args
    global geneidx
    global teidx
    global introidx

    info = args.info
    info("\nReading sample files ...\n")
    
    csamples_ele_tbl = None
    csamples_instance_tbl = None
    discounted_c_ele_tbl = None
    discounted_c_i_tbl = None
    csamples_rpm = 0
    procIdx = multiprocessing.current_process()._identify[0]
    tfiles  = tfileArry[procIdx] 

    (tsamples_ele_tbl, tsamples_instance_tbl, discounted_t_ele_tbl, discounted_t_i_tbl, tsamples_rpm) = count_reads(tfiles, args.parser, geneIdx, teIdx, intronIdx, args.stranded, args.te_mode, args.sortByPos, args.prj_name,args.numItr,args.fragLength,args.maxL, args.b_origcnt)
    if args.cfiles:
        (csamples_ele_tbl, csamples_instance_tbl, discounted_c_ele_tbl, discounted_c_i_tbl, csamples_rpm) = count_reads(cfiles, args.parser, geneIdx, teIdx, intronIdx, args.stranded, args.te_mode, args.sortByPos, args.prj_name,args.numItr,args.fragLength,args.maxL, args.b_origcnt)

    info("Finished processing sample files")

    info("Generating counts table")

    ForkedPdb().set_trace()
    f_ele_tbl = args.prj_name + multiprocessing.current_process()._name + ".ele.cntTable"
    f_instance_tbl = args.prj_name + multiprocessing.current_process()._name + ".instance.cntTable"
    f_discounted_ele_tbl = args.prj_name + multiprocessing.current_process()._name + ".discount.ele.cntTable"
    f_discounted_instance_tbl = args.prj_name + multiprocessing.current_process()._name + ".discount.instance.cntTable"
    output_count_tbl(tsamples_ele_tbl, csamples_ele_tbl, f_ele_tbl)
    output_count_tbl(tsamples_instance_tbl, csamples_instance_tbl, f_instance_tbl)
    output_count_tbl(discounted_t_ele_tbl, discounted_c_ele_tbl, f_discounted_ele_tbl)
    output_count_tbl(discounted_t_i_tbl, discounted_c_i_tbl, f_discounted_instance_tbl)
    
def readFileArray(fileHandler, np)
    fileArray = []
    for x in range(0, np):
        line = tfileHandler.readline()  
    if "" == line:
        break
    fileArry.append(line)
    return fileArray

args = None
teIdx = None
geneIdx = None
intronIdx = None
tfileArray = None


# Main function of script
def main():
    """Start TEtranscripts......parse options......"""
    global args
    global teIdx
    global geneIdx
    global intronIdx
    global tfileArray
    args=read_opts2(prepare_parser())
    
    ForkedPdb().set_trace()
    info = args.info
    #warn = args.warn
    #debug = args.debug
    error = args.error

    # Output arguments used for program
    info("\n" + args.argtxt + "\n")

    info("Processing GTF files ...\n")

    #(features, genes) = read_features(args.gtffile, args.stranded, "exon", "gene_id")
    info("Building gene index .......\n")
    geneIdx = GeneFeatures(args.gtffile,args.stranded,"exon","gene_id")
    info("Done building gene index ......\n")
    #TE index
#try :
    teIdx = TEfeatures()
    cur_time = time.time()
    te_tmpfile = '.'+str(cur_time)+'.te.gtf'
    if (args.b_pickle): 
        te_picklefile = '.te.HI.pickle.2.7'
        if os.path.isfile(te_picklefile) and os.path.getsize(te_picklefile) > 0:
            info("\nLoading TE pickle .......\n")
            with open(te_picklefile, 'rb') as fp:
                gc.disable()    # disable garbage collector
                teIdx = pickle.load(fp)
                gc.enable()
            info("\nDone loading TE pickle .......\n")
        else:
            info("\nBuilding TE index .......\n")
            #subprocess.call(['sort -k 1,1 -k 4,4g '+ args.tefile+ ' >'+ te_tmpfile],shell=True )
            #teIdx.build(te_tmpfile,args.te_mode)
            teIdx.build(args.tefile,args.te_mode)
            #subprocess.call(['rm -f ' + te_tmpfile ],shell=True)
            info("Done building TE index ......\n")
            if args.intronfile:
                if args.exontefile:
                    info("\nReading exonic TEs .......\n")
                    teIdx.dictExonicTEs(args.exontefile) 
                    info("\nDone reading exonic TEs .......\n")
                else: 
                    sys.stderr.write("if intron gtf is provided exonicTEs and intergeneTEs have to be provided as well. \n")
                if args.intergenictefile:
                    info("\nReading intergenic TEs .......\n")
                    teIdx.dictIntergenicTEs(args.intergenictefile) 
                    info("\nDone reading intergenic TEs .......\n")
                else:
                    sys.stderr.write("if intron gtf is provided exonicTEs and intergeneTEs have to be provided as well. \n")
            info("\nDumping TE pickle .......\n")
            with open(te_picklefile, 'wb') as fp:
                pickle.dump(teIdx, fp, protocol=pickle.HIGHEST_PROTOCOL)
            info("\nDone dumping TE pickle .......\n")
    else: 
        info("\nBuilding TE index .......\n")
        #subprocess.call(['sort -k 1,1 -k 4,4g '+ args.tefile+ ' >'+ te_tmpfile],shell=True )
        #teIdx.build(te_tmpfile,args.te_mode)
        teIdx.build(args.tefile,args.te_mode)
        #subprocess.call(['rm -f ' + te_tmpfile ],shell=True)
        info("Done building TE index ......\n")
        if args.intronfile:
            if args.exontefile:
                info("\nReading exonic TEs .......\n")
                teIdx.dictExonicTEs(args.exontefile) 
                info("\nDone reading exonic TEs .......\n")
            else: 
                sys.stderr.write("if intron gtf is provided exonicTEs and intergeneTEs have to be provided as well. \n")
            if args.intergenictefile:
                info("\nReading intergenic TEs .......\n")
                teIdx.dictIntergenicTEs(args.intergenictefile) 
                info("\nDone reading intergenic TEs .......\n")
            else:
                sys.stderr.write("if intron gtf is provided exonicTEs and intergeneTEs have to be provided as well. \n")


#except :
#    sys.stderr.write("Error in building TE index \n")
#    sys.exit(1)

    #TEflanking introns index
    if args.intronfile:
    #try :
        intronIdx = IntronFeatures(teIdx)
        cur_time = time.time()
        if (args.b_pickle): 
            intron_picklefile = '.intron.HI.pickle.2.7'
            if os.path.isfile(intron_picklefile) and os.path.getsize(intron_picklefile) > 0:
                info("\nLoading intron pickle .......\n")
                with open(intron_picklefile, 'rb') as fp:
                    gc.disable()    # disable garbage collector
                    intronIdx = pickle.load(fp)
                    gc.enable()
                info("\nDone loading intron pickle .......\n")
            else:
                info("\nBuilding intron index .......\n")
                intronIdx.build(args.intronfile,args.te_mode)
                info("Done building intron index ......\n")
                info("\nDumping intron pickle .......\n")
                with open(intron_picklefile, 'wb') as fp:
                    pickle.dump(intronIdx, fp, protocol=pickle.HIGHEST_PROTOCOL)
                info("\nDone dumping intron pickle .......\n")
        else: 
            info("\nBuilding intron index .......\n")
            intronIdx.build(args.intronfile,args.te_mode)
            info("Done building intron index ......\n")

    #except :
    #    sys.stderr.write("Error in building intron flanking intron index \n")
    #    sys.exit(1)




    # Read sample files make count table
    info("\nReading sample files ...\n")

    #teIdxProxy, teIdxProxyListener = proxy.createProxy(teIdx)
    #argsProxy, argsProxyListener = proxy.createProxy(args)
    #geneIdxProxy, geneIdxProxyListener = proxy.createProxy(geneIdx)
    #intronIdxProxy, intronIdxProxyListener = proxy.createProxy(intronIdx)
    #finished = GetSetter()
    #finishedProxy, finishedProxyListener = proxy.createProxy(finished)
    end = False
    while end == False:
        fp = open(args.tfile, 'r')
        tfileArray = readTfile(args.tfile, np) 
        if args.cfile
            fp = open(args.cfile, 'r')
            cfileArray = readCfile(args.cfile. np)
        if tfileArray.size() == 0:
            end = True 
            break
        p = Pool(tfileArray.size())
        Pool.map(subprocessWorker)

    p2 = Process(target=subprocessWorker, args=())

    p2.start()
    #p3.start()
    #p4.start()
    #while finished.get() < 99:
    #    argsProxyListener.listen()
    #    geneIdxProxyListener.listen()
    #    teIdxProxyListener.listen()
    #    intronIdxProxyListener.listen()
    #    finishedProxyListener.listen()

    p2.join()
    #p3.join()
    #p4.join()
    fp.close()

        
    #csamples_ele_tbl = None
    #csamples_instance_tbl = None
    #discounted_c_ele_tbl = None
    #discounted_c_i_tbl = None
    #csamples_rpm = 0
    #(tsamples_ele_tbl, tsamples_instance_tbl, discounted_t_ele_tbl, discounted_t_i_tbl, tsamples_rpm) = count_reads(args.tfiles, args.parser, geneIdx, teIdx, intronIdx, args.stranded, args.te_mode, args.sortByPos, args.prj_name,args.numItr,args.fragLength,args.maxL, args.b_origcnt)
    #if args.cfiles:
    #    (csamples_ele_tbl, csamples_instance_tbl, discounted_c_ele_tbl, discounted_c_i_tbl, csamples_rpm) = count_reads(args.cfiles, args.parser, geneIdx, teIdx, intronIdx, args.stranded, args.te_mode, args.sortByPos, args.prj_name,args.numItr,args.fragLength,args.maxL, args.b_origcnt)
    #
    #info("Finished processing sample files")

    #info("Generating counts table")

    #f_ele_tbl = args.prj_name + ".ele.cntTable"
    #f_instance_tbl = args.prj_name + ".instance.cntTable"
    #f_discounted_ele_tbl = args.prj_name + ".discount.ele.cntTable"
    #f_discounted_instance_tbl = args.prj_name + ".discount.instance.cntTable"
    #output_count_tbl(tsamples_ele_tbl, csamples_ele_tbl, f_ele_tbl)
    #output_count_tbl(tsamples_instance_tbl, csamples_instance_tbl, f_instance_tbl)
    #output_count_tbl(discounted_t_ele_tbl, discounted_c_ele_tbl, f_discounted_ele_tbl)
    #output_count_tbl(discounted_t_i_tbl, discounted_c_i_tbl, f_discounted_instance_tbl)
    #if args.cfiles:
    #    rpm_val = tsamples_rpm + csamples_rpm
    #else:
    #    rpm_val = tsamples_rpm 

    sys.exit(0)
    #info("Calculating differential expression ...\n")
    # Obtaining R-code for differential analysis

    rscript = write_R_code(f_ele_tbl, args.tfiles, args.cfiles, args.prj_name, args.norm, args.pval, args.fc, rpm_val,args.min_read)

    f_rscript = args.prj_name + '_DESeq.R'
    rcode = open('%s' % (f_rscript) , 'w')
    rcode.write(rscript)
    rcode.close()

    # Running R-code for differential analysi
    try:
        sts = subprocess.call(['Rscript', f_rscript])
    except:
        error("Error in running differential analysis!\n")
        error("Error: %s\n" % str(sys.exc_info()[1]))
        error( "[Exception type: %s, raised in %s:%d]\n" % ( sys.exc_info()[1].__class__.__name__,
                os.path.basename(traceback.extract_tb( sys.exc_info()[2] )[-1][0]),
                traceback.extract_tb( sys.exc_info()[2] )[-1][1] ) )
        sys.exit(1)

    info("Done \n")


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt !\n")
        sys.exit(0)
