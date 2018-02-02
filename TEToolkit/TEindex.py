

"""Module Description

Copyright (c) 2014, Ying Jin <yjin@cshl.edu >


This code is free software; you can redistribute it and/or modify it
under the terms of the Artistic License (see the file COPYING included
with the distribution).

@author:  Ying Jin
@contact: yjin@cshl.edu
"""
import sys, time
import logging
from math import ceil,floor
from ItvTree import ItvTree
from Constants import TEindex_BINSIZE
from Constants import ReadLength

#TEindex_BINSIZE = 200

class TEfeatures:
    """index of TE annotations.
    """
    def __init__ (self):

        self.indexlist = {}
        self._length = []
        self._nameIDmap = []
        self._name2IDdict = {}
        self._elements = []
        self._leftIntronIDmap = [] 
        self._rightIntronIDmap = [] 
        self._intergenicTEs = {}
        self._exonicTEs = {}

    def getNames(self) :
        names = []
        return self._nameIDmap

    def numInstances(self) :
        return len(self._nameIDmap)

    def getElements(self) :
        return self._elements

    def getIdxbyName(self, name):
        return self._name2IDdict[name]

    def initLeftIntronIDmap(self,maxSize) :
        self._leftIntronIDmap = [None]*maxSize

    def initRightIntronIDmap(self,maxSize) :
        self._rightIntronIDmap = [None]*maxSize

    def getLeftIntrons(self) :
        return self._leftIntronIDmap

    def getRightIntrons(self) :
        return self._rightIntronIDmap

    def setLeftIntronIdx(self,idx, val) :
        self._leftIntronIDmap[idx] = val

    def setRightIntronIdx(self,idx, val) :
        self._rightIntronIDmap[idx] = val

    def getStrand(self,idx) :
        f_name = self._nameIDmap[idx]
        return f_name[len(f_name)-1]

    def getInstanceName(self,idx) :
        if idx >= len(self._nameIDmap) or idx < 0 :
            return None
        else :
            full_name =  self._nameIDmap[idx]
            pos = full_name.find(':')
            val = full_name[:pos]
            return val


    def getEleName(self,idx) :
        full_name = None
        if idx >= len(self._nameIDmap) or idx < 0 :
            return None
        else :
            full_name =  self._nameIDmap[idx]
        if full_name is not None:
            pos = full_name.find(':')
            val = full_name[pos+1:(len(full_name)-2)]
            return val
        else :
            return None

    def getFullName(self,idx) :
        if idx >= len(self._nameIDmap) or idx < 0 :
            return None
        else :
            return self._nameIDmap[idx]

    def getLength(self,TE_name_idx) :
        if TE_name_idx < len(self._length) :
            return self._length[TE_name_idx]
        else :
            return -1

    def getFamilyID(self,chr,start,end):
        binID = start/TEindex_BINSIZE
        endbinID = end/TEindex_BINSIZE + 1

        if self.indexlist.has_key(chr) :
            index = self.indexlist[chr]
            (node,RBnode) = index.lookup(binID,index._root,None,None)

            if node is not None and node.overlaps(binID,endbinID) :
                full_name = (node.getName()).split(':')
                famid = full_name[2]
                return famid
            else :
                return None
        else :

            return None

    # looks up the TE index interval tree with start end, 
    # find the nodes corresponding to the intervals 
    # returns all overlapping TE name indices within the node

    def findOvpTE(self,chrom,start,end):
        startbinID = start/TEindex_BINSIZE
        endbinID = end/TEindex_BINSIZE
        if start == startbinID * TEindex_BINSIZE :
           startbinID -= 1
        name_idx_list = []

        if  self.indexlist.has_key(chrom) :
               index = self.indexlist[chrom]
        else :
            return None

        (LBnode,RBnode) = index.lookup_r(startbinID,endbinID,index._root)

        if LBnode is not None :
            telist = LBnode.overlaps(start,end)
            name_idx_list.extend(telist)

        if RBnode is not None :
            telist = RBnode.overlaps(start,end)
            name_idx_list.extend(telist)

        return name_idx_list

    def TE_annotation(self,iv_seq):
        TEs = []
        for iv in iv_seq :
            chromo = iv[0]
            start = iv[1]
            end = iv[2]
            strand = iv[3]
            name_idx_list  = self.findOvpTE(chromo,start,end)
            #print "name_idx_list: ", name_idx_list

            if name_idx_list is not None :
                for (t,ovp_len) in name_idx_list :
                    if strand != "." : #stranded
                        if (strand == self.getStrand(t)) and (ovp_len >= 5) :   # strand matches and overlap is at least 5 bases
                            if t not in TEs :
                                #print self.getFullName(t), " appended"
                                TEs.append([t,ovp_len])
                            #else:
                                #print self.getFullName(t), " already in the list" 
                        #else: 
                            #print "strand ", strand, " != ", self.getStrand(t)
                    else :#not stranded
                        if ovp_len >= 5 :   # overlap is at least 5 bases
                            if t not in TEs :
                                TEs.append([t,ovp_len])

        return TEs

    #group by element
    def groupByEle(self,te_inst_counts) :

        TEs = self.getElements()
        te_ele_counts = dict(zip(TEs,[0]*len(TEs)))

        for i in range(len(te_inst_counts)) :
            ele_name = self.getEleName(i)

            if ele_name is None:
                sys.stderr.write("TE out of index boundary!\n")
                sys.exit(1)

            if ele_name in te_ele_counts :
                te_ele_counts[ele_name] += te_inst_counts[i]
            else :
                sys.stderr.write("TE inconsistency! "+ele_name+"\n")
                sys.exit(1)

        return te_ele_counts

    # count by Instance name
    def countByName(self,te_inst_counts) :
        #print "countByName"
        TEnames = self._nameIDmap
        #print TEnames
        #print te_inst_counts
        assert len(TEnames) == len(te_inst_counts)
        te_name_counts = dict(zip(TEnames,te_inst_counts))
        return te_name_counts




    def build (self,filename,te_mode):
        self.__srcfile = filename

        try:
            f = open(self.__srcfile,'r')
        except:
            logging.error("cannot open such file %s !\n" %(self.__srcfile))
            sys.exit(1)

        name_idx = 0
        linenum = 0
        for line in f :
            line = line.strip()
            items = line.split('\t')
            chrom = items[0]
            start = int(items[3])
            end = int(items[4])
            strand = items[6]
            items[8] = items[8].replace("; ",";")
            desc = items[8].split(';')
            name = ""
            family_id = ""
            ele_id = ""
            class_id = ""
            tlen = end - start + 1
            linenum += 1

            for i in range(len(desc)) :
                desc[i] = desc[i].replace("\"","")
                pos = desc[i].find(" ")
                tid = desc[i][:pos]
                val = desc[i][pos+1:len(desc[i])]

                if tid == "gene_id" :
                    ele_id = val
                if tid == "transcript_id" :
                    name = val
                if tid == "family_id" :
                    family_id = val
                if tid == "class_id" :
                    class_id = val

            if ele_id == "" or name == "" or family_id == "" or class_id == "" :
                sys.stderr.write(line+"\n")
                sys.stderr.write("TE GTF format error! There is no annotation at line %s.\n" % (linenum))
                raise

            # name of TE is e.g.
            # AluSp:AluSp:Alu:SINE:+
            full_name = name+':'+ele_id+':'+family_id+':'+class_id+':'+strand
            ele_name = ele_id+':'+family_id+':'+class_id
            if ele_name not in self._elements :
                self._elements.append(ele_name)

            self._length.append(tlen)
            assert len(self._nameIDmap) == name_idx
            self._nameIDmap.append(full_name)
            self._name2IDdict[name] = name_idx

            if self.indexlist.has_key(chrom) :
                    index = self.indexlist[chrom]

                    bin_startID = start/TEindex_BINSIZE
                    bin_endID = end/TEindex_BINSIZE
                    if start == bin_startID * TEindex_BINSIZE :
                        bin_startID -= 1
                    while bin_startID <= bin_endID :
                        end_pos = min(end,(bin_startID+1) * TEindex_BINSIZE )
                        start_pos = max(start,bin_startID * TEindex_BINSIZE+1)

                        index.insert(start_pos,end_pos,name_idx)
                        bin_startID += 1

            else :
                    index = ItvTree()
                    bin_startID = start/TEindex_BINSIZE
                    bin_endID = end/TEindex_BINSIZE
                    if start == bin_startID * TEindex_BINSIZE :
                        bin_startID -= 1
                    while bin_startID <= bin_endID :
                        end_pos = min(end,(bin_startID+1) * TEindex_BINSIZE )
                        start_pos = max(start,bin_startID * TEindex_BINSIZE+1)
                        index.insert(start_pos,end_pos,name_idx)
                        bin_startID += 1

                    self.indexlist[chrom] = index

            name_idx += 1

        f.close()


    def dictIntergenicTEs(self,filename):
        newdict = {}
        try:
            f = open(filename,'r')
        except:
            logging.error("cannot open such file %s !\n" %(filename))
            sys.exit(1)

        linenum = 0
        for line in f :
            line = line.strip()
            items = line.split('\t')
            strand = items[6]
            items[8] = items[8].replace("; ",";")
            desc = items[8].split(';')
            name = ""
            family_id = ""
            ele_id = ""
            class_id = ""
            linenum += 1

            for i in range(len(desc)) :
                desc[i] = desc[i].replace("\"","")
                pos = desc[i].find(" ")
                tid = desc[i][:pos]
                val = desc[i][pos+1:len(desc[i])]

                if tid == "gene_id" :
                    ele_id = val
                if tid == "transcript_id" :
                    name = val
                if tid == "family_id" :
                    family_id = val
                if tid == "class_id" :
                    class_id = val

            if ele_id == "" or name == "" or family_id == "" or class_id == "" :
                sys.stderr.write(line+"\n")
                sys.stderr.write("TE GTF format error! There is no annotation at line %s.\n" % (linenum))
                raise

            full_name = name+':'+ele_id+':'+family_id+':'+class_id+':'+strand
            newdict[full_name] = 1

        #print "intergenicTEs"
        #print newdict
        self._intergenicTEs = newdict


    def dictExonicTEs(self,filename):
        newdict = {}
        try:
            f = open(filename,'r')
        except:
            logging.error("cannot open such file %s !\n" %(filename))
            sys.exit(1)

        linenum = 0
        for line in f :
            line = line.strip()
            items = line.split('\t')
            strand = items[6]
            items[8] = items[8].replace("; ",";")
            desc = items[8].split(';')
            name = ""
            family_id = ""
            ele_id = ""
            class_id = ""
            linenum += 1

            for i in range(len(desc)) :
                desc[i] = desc[i].replace("\"","")
                pos = desc[i].find(" ")
                tid = desc[i][:pos]
                val = desc[i][pos+1:len(desc[i])]

                if tid == "gene_id" :
                    ele_id = val
                if tid == "transcript_id" :
                    name = val
                if tid == "family_id" :
                    family_id = val
                if tid == "class_id" :
                    class_id = val

            if ele_id == "" or name == "" or family_id == "" or class_id == "" :
                sys.stderr.write(line+"\n")
                sys.stderr.write("TE GTF format error! There is no annotation at line %s.\n" % (linenum))
                raise

            full_name = name+':'+ele_id+':'+family_id+':'+class_id+':'+strand
            newdict[full_name] = 1

        #print "exonicTEs"
        #print newdict
        self._exonicTEs = newdict
        







####################################
#intron annotation

class IntronFeatures(TEfeatures, object):
    """index of intron annotations flanking TEs.
    """
    def __init__(self, TEidx):
        self._start = []
        self._end = []
        self._TEidx = TEidx # TEfeatures idx object needs to be created already

        super(IntronFeatures, self).__init__()
        self._TEidx.initLeftIntronIDmap(self._TEidx.numInstances())
        self._TEidx.initRightIntronIDmap(self._TEidx.numInstances())
        

    def getTEidx(self) :
        return self._TEidx
        
    def getStart(self,TE_name_idx) :
        if TE_name_idx < len(self._start) :
            return self._start[TE_name_idx]
        else :
            return -1

    def getEnd(self,TE_name_idx) :
        if TE_name_idx < len(self._end) :
            return self._end[TE_name_idx]
        else :
            return -1


    # name of intron is e.g.
    # AluSp_dup2_right:AluSp:Alu:SINE:+;L2b_dup12455_left:L2b:L2:LINE:+"
    def getEleName(self,idx) :
        full_name = None
        if idx >= len(self._nameIDmap) or idx < 0 :
            return None
        else :
            full_name =  self._nameIDmap[idx]
        if full_name is not None:
            vals = []
            names = full_name.split(";")
            for i in range(len(names)):
                pos = names[i].find(':')
                vals.append( names[i][pos+1:(len(names[i])-2)])
            return vals
        else :
            return None

    #group by element
    def groupByEle(self,te_inst_counts) :

        TEs = self.getElements()
        #print TEs
        #print self.getNames() 
        #print te_inst_counts
        te_ele_counts = dict(zip(TEs,[0]*len(TEs)))

        for i in range(len(te_inst_counts)) :
            ele_names = self.getEleName(i)
            for j in range(len(ele_names)): 
                ele_name = ele_names[j]
                #print ele_name
                if ele_name is None:
                    sys.stderr.write("TE out of index boundary!\n")
                    sys.exit(1)

                if ele_name in te_ele_counts :
                    #print te_ele_counts[ele_name]
                    te_ele_counts[ele_name] += te_inst_counts[i]
                else :
                    sys.stderr.write("TE inconsistency! "+ele_name+"\n")
                    sys.exit(1)

        return te_ele_counts




    # discount TEcounts by flanking intron depth
    def discountByFlankingIntrons(self,te_inst_counts, intron_counts) :
        new_te_inst_counts = list(te_inst_counts) 
        TEidx = self._TEidx
        TEnames = TEidx.getNames()
        leftIntrons = TEidx.getLeftIntrons()
        rightIntrons = TEidx.getRightIntrons()
        assert len(TEnames) == len(te_inst_counts)
        assert len(TEnames) == len(leftIntrons)
        assert len(TEnames) == len(rightIntrons)
        #print "exonicTEs"
        #print TEidx._exonicTEs
        #print "intergenicTEs"
        #print TEidx._intergenicTEs
        for i in range(len(TEnames)):
            if TEidx._exonicTEs.has_key(TEnames[i]) :        # cannot reliably count TEs in exons (including lncRNA). 
                new_te_inst_counts[i] = 0
                continue
            if TEidx._intergenicTEs.has_key(TEnames[i]) :    # no flanking intron or pre-mRNA to worry about. 
                new_te_inst_counts[i] = te_inst_counts[i]
                continue

            # remaining TEs are within introns
            TEcnt = te_inst_counts[i]
            TElen = TEidx.getLength(i)
            if (TElen == 0):
                print "length zero!"
            TEefflen = TElen - ReadLength + 1
            if TEefflen <= 0:
                TEefflen = 1
            TEdepth = TEcnt/float(TEefflen)
            leftdepth = 0
            rightdepth = 0
            leftintronidx = leftIntrons[i]
            if leftintronidx:
                leftintronlen = self.getLength(leftintronidx)
                leftefflen = leftintronlen - ReadLength + 1
                if (leftefflen <= 0): 
                    leftefflen = 1
                leftintroncnt = intron_counts[leftintronidx]
                #print "leftcnt", leftintroncnt, " leftlen", leftintronlen    
                leftdepth = leftintroncnt/float(leftefflen)
            rightintronidx = rightIntrons[i]
            if rightintronidx:
                rightintronlen = self.getLength(rightintronidx)
                rightefflen = rightintronlen - ReadLength + 1
                if (rightefflen <= 0): 
                    rightefflen = 1
                rightintroncnt = intron_counts[rightintronidx]
                #print "rightcnt", rightintroncnt, " rightlen", rightintronlen    
                rightdepth = rightintroncnt/float(rightefflen)

            # assign new counts based on intron information
            if ((leftintronidx is None) and (rightintronidx is None)):      # maybe merged TEs
                print "missing introns: ", TEnames[i]
                new_te_inst_counts[i] = 0
            else:
                meanintrondepth = float(leftdepth + rightdepth)/((leftintronidx is not None) + (rightintronidx is not None))
                if (TEdepth < 1/float(ReadLength)):
                    if (meanintrondepth == 0):
                        new_te_inst_counts[i] = te_inst_counts[i]
                    else:
                        new_te_inst_counts[i] = 0
                else:
                    discount = meanintrondepth/float(TEdepth)
                    print "discounted: ", TEnames[i], " ",discount
                    new_te_inst_counts[i] = TEcnt - TEcnt*discount
                    if new_te_inst_counts[i] < 0:
                        new_te_inst_counts[i] = 0
                #print leftdepth, " ", rightdepth, " ", TEdepth," ", TEcnt, " ", new_te_inst_counts[i] 
        return new_te_inst_counts




    def build(self,filename,te_mode):
        self.__srcfile = filename

        try:
            f = open(self.__srcfile,'r')
        except:
            logging.error("cannot open such file %s !\n" %(self.__srcfile))
            sys.exit(1)

        name_idx = 0
        linenum = 0
        for line in f :
            line = line.strip()
            items = line.split('\t')
            chrom = items[0]
            start = int(items[3])
            end = int(items[4])
            strand = items[6]
            items[8] = items[8].replace("; ",";")
            desc = items[8].split(';')
            name = []
            family_id = [] 
            ele_id = []
            class_id = []
            tlen = end - start + 1
            linenum += 1
            full_name = ""
            left_TE_names = [] 
            right_TE_names = [] 

            for i in range(len(desc)) :
                desc[i] = desc[i].replace("\"","")
                #print desc[i]
                pos = desc[i].find(" ")
                tid = desc[i][:pos]
                val = desc[i][pos+1:len(desc[i])]

                if tid == "gene_id" :
                    ele_id.append (val)
                if tid == "transcript_id" :
                    name.append (val)
                if tid == "family_id" :
                    family_id.append (val)
                if tid == "class_id" :
                    class_id.append (val)

            if ele_id == [] or name == [] or family_id == [] or class_id == [] :
                sys.stderr.write(line+"\n")
                sys.stderr.write("TE GTF format error! There is no annotation at line %s.\n" % (linenum))
                raise

            # name of intron is e.g.
            # AluSp_dup2_right:AluSp:Alu:SINE:+;L2b_dup12455_left:L2b:L2:LINE:+"
            for i in range(len(name)) :
                #print name[i]
                pos = name[i].rfind("_")
                TEinstancename = name[i][:pos]
                TEnextto = name[i][pos+1:]
                if full_name != "":
                    full_name += ";"
                full_name += name[i]+':'+ele_id[i]+':'+family_id[i]+':'+class_id[i]+':'+strand
                ele_name = ele_id[i]+':'+family_id[i]+':'+class_id[i]
                if ele_name not in self._elements :
                    self._elements.append(ele_name)

                # find out adjacent TE and save the idx
                if (TEnextto == "left"): 
                    #TEidx_list  = self._TEidx.findOvpTE(chrom,end+1,end+1)
                    ##print "intron_left2TE: ", TEidx_list
                    #if TEidx_list: 
                    #    #print self._TEidx.getInstanceName(TEidx_list[0][0]), "==", TEinstancename
                    #    #if (self._TEidx.getInstanceName(TEidx_list[0][0]) == TEinstancename):
                    #        self._TEidx.setLeftIntronIdx(TEidx_list[0][0], name_idx)
                    self._TEidx.setLeftIntronIdx(self._TEidx.getIdxbyName(TEinstancename), name_idx)
                if (TEnextto == "right"): 
                    #TEidx_list  = self._TEidx.findOvpTE(chrom,start-1,start-1)
                    ##print "intron_right2TE: ", TEidx_list
                    #if TEidx_list: 
                    #    #print self._TEidx.getInstanceName(TEidx_list[0][0]), "==", TEinstancename
                    #    if (self._TEidx.getInstanceName(TEidx_list[0][0]) == TEinstancename):
                    #        self._TEidx.setRightIntronIdx(TEidx_list[0][0], name_idx)
                    self._TEidx.setRightIntronIdx(self._TEidx.getIdxbyName(TEinstancename), name_idx)

            assert len(self._nameIDmap) == name_idx
            self._start.append(start)
            self._end.append(end)
            self._length.append(tlen)
            self._nameIDmap.append(full_name)

            # find interval tree for chr in self.indexlist
            if self.indexlist.has_key(chrom) :
                index = self.indexlist[chrom]
            else :
                index = ItvTree()
                self.indexlist[chrom] = index

            # save interval in the interval tree
            bin_startID = start/TEindex_BINSIZE
            bin_endID = end/TEindex_BINSIZE
            if start == bin_startID * TEindex_BINSIZE :
                bin_startID -= 1
            while bin_startID <= bin_endID :
                end_pos = min(end,(bin_startID+1) * TEindex_BINSIZE )
                start_pos = max(start,bin_startID * TEindex_BINSIZE+1)

                index.insert(start_pos,end_pos,name_idx)
                bin_startID += 1

            name_idx += 1

        f.close()

    
    def findOvpIntron(self,chrom,start,end):
        startbinID = start/TEindex_BINSIZE
        endbinID = end/TEindex_BINSIZE
        if start == startbinID * TEindex_BINSIZE :
           startbinID -= 1
        name_idx_list = []

        if  self.indexlist.has_key(chrom) :
               index = self.indexlist[chrom]
        else :
            return None

        (LBnode,RBnode) = index.lookup_r(startbinID,endbinID,index._root)

        if LBnode is not None :
            intronlist = LBnode.overlaps(start,end)
            name_idx_list.extend(intronlist)

        if RBnode is not None :
            intronlist = RBnode.overlaps(start,end)
            name_idx_list.extend(intronlist)

        return name_idx_list

    def get_flankingTE(self,intron_idx):
        intron_name = self.getFullName(intron_idx)
        pos = intron_name.find(':')
        val = intron_name[:pos]
        pos = val.rfind('_')        
        te_name = val[:pos]
        direction = val[pos+1:]
        #print te_name, " ", direction


    def intron_annotation(self,iv_seq):
        introns = []
        for iv in iv_seq :
            chromo = iv[0]
            start = iv[1]
            end = iv[2]
            strand = iv[3]
            name_idx_list  = self.findOvpIntron(chromo,start,end)
            #print "name_idx_list: ", name_idx_list

            if name_idx_list is not None :
                for (t,ovp_len) in name_idx_list :
                    if strand != "." : #stranded
                        if strand == self.getStrand(t)  :
                            if t not in introns :
                                #print self.getFullName(t), " appended"
                                self.get_flankingTE(t)
                                introns.append([t,ovp_len])
                            #else:
                                #print self.getFullName(t), " already in introns"
                        #else:
                            #print "strand ", strand, " != ", self.getStrand(t)
                    else :#not stranded
                        if t not in introns :
                            introns.append([t,ovp_len])

        return introns





if __name__ == '__main__':
    try:

        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt !\n")
        sys.exit(0)
