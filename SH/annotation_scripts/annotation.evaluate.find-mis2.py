#!/usr/bin/env python
# encoding: utf-8

'''
Created on Sep 16, 2017
@author: jiao@mpip.mpg.de

## splitting ERROR  (two genes blastp with one Ref-0 gene, both has high identity)
## merging ERROR (one gene blastp with two Ref-0 gene, both has high identity)
## missing genes (good blastn hit, and protein exonerate hit or RNA-seq support)
## wrong gene structure (if scipo gene model has ortholog??)

'''

import sys
import os
import getopt
import re
from Bio import SeqIO

def main(argv):
    groupFile = ""
    outdir = ""
    ref = ""
    alt = ""
    blastnRes = ""
    blastpRes = ""
    refBed = ""
    altBed = ""
    refProtFas = ""
    altProtFas = ""
    #ara11bed = ""
    blastnResRef = ""
    LofRefFile = ""
    LofAltFile = ""
    rnaGff = ""
    try:
        opts, args = getopt.getopt(argv,"g:o:n:p:s:q:x:y:m:r:",["group=","outdir=","blnAlt=","blastp=","refBed=","altBed=","refProt=","altProt=","blnRef=","rna=",])
    except getopt.GetoptError:
        print 'annotation.evaluate.find-mis.py -g <group> -o <outdir> -n <blastn> -p <blastp> -s <refBed> -q <altBed> -x <refProt> -y <altProt>  -m <blastnRef> -r <rna>'
        sys.exit(2)
    if len(opts) == 0 :
        print 'annotation.evaluate.find-mis.py -g <group> -o <outdir> -n <blastn> -p <blastp>  -s <refBed> -q <altBed> -x <refProt> -y <altProt> -m <blastnRef>'
        sys.exit()
    for opt, arg in opts:
        if opt == '-h':
            print 'annotation.evaluate.find-mis.py -g <group> -o <outdir>  -n <blastn> -p <blastp> -s <refBed> -q <altBed> -x <refProt> -y <altProt>  -m <blastnRef>'
            sys.exit()
        elif opt in ("-g", "--group"):
            groupFile = arg
        elif opt in ("-o", "--out"):
            outdir = arg
        elif opt in ("-n", "--blnAlt"):
            refBlnAltRes = arg
        elif opt in ("-p", "--blastp"):
            blastpRes = arg
        elif opt in ("-s", "--refBed"):
            refBed = arg
        elif opt in ("-q", "--altBed"):
            altBed = arg
        elif opt in ("-x", "--refProt"):
            refProtFas = arg
        elif opt in ("-y", "--altProt"):
            altProtFas = arg
        elif opt in ("-m"," --blnRef") :
            altBlnRefRes = arg
        elif opt in ("-r","--rna"):
            rnaGff = arg
    ##some cutoffs
    minId = 70
    maxIdf = 10
    maxSplitCov = 0.85

    if not os.path.isdir(outdir) :
        os.system("mkdir -p " + outdir)

    refLen = getProtLeng(refProtFas)
    altLen = getProtLeng(altProtFas)

    [refIdx, refPos]= getGeneIndex(refBed)
    [altIdx, altPos]= getGeneIndex(altBed)

    refGeneBed = getBed(refBed)
    altGeneBed = getBed(altBed)

    refBest = getBesthit(refBlnAltRes)
    altBest = getBesthit(altBlnRefRes)

    ## get specific/un-grouped genes
    ref=os.path.basename(refProtFas).split(".")[0].upper()
    alt=os.path.basename(altProtFas).split(".")[0].upper()
    print(ref,"ref")
    [ungrRef,ungrAlt,ort, cnv] = getSpecGene(groupFile,refLen,altLen,ref,alt)

    '''
      Method 1: evaluation based on blastp [Ref-Acc] result and Ref-gene-acc_assembly blastn result
    '''
    #check mis-merged genes
    blastpFlt1 = outdir + "/blastp.flt.out1"
    fltBlastp(blastpRes, blastpFlt1, refIdx, altIdx, altLen, refLen)
    [altMisMerBlp,misMerBlpRef,misMerGeneAlt]= findMisMerBlp(blastpFlt1,outdir,minId,maxSplitCov,maxIdf,refIdx,altIdx,altLen,altPos, refLen)

    #check mis-splitting genes
    [altMisSplBlp,misSplBlpRef,misSplBlpAlt]= findMisSplitBlp(blastpFlt1,outdir,minId,maxSplitCov,maxIdf,refIdx,altIdx,refLen,altPos)

    '''
      Method 2: evaluation based on  Ref-gene-acc_assembly blastn result
    1. genes: mis-merged
    2. genes: mis-split
    3. genes: mis-exon-intron stucture
    4. genes: mis-annotated (not annotated in Ref)
    5. genes: not annotated in accession assembly, to be added
    '''

    blastpFlt2 = outdir + "/blastp.flt.out2"
    [blastp, bestBlastp] = getBlastp(blastpRes, blastpFlt2, refLen, altLen)

    [unAssRef,unAnnRef,misExRef,misExAlt] = findAltMM(refBlnAltRes,altBed,outdir,ort,ungrRef, ungrAlt,refIdx,altIdx, blastp)
    
    [unAssAlt, unAnnAlt, rnaSup] = findRefMM(altBlnRefRes, refBed, altBed, outdir, ungrAlt, ungrRef, rnaGff, blastp)

    inFile = outdir + "/alt-genome.mis-merged.ref-gene.txt"
    [misMerBln, misMerBlnRef] = findMerBlastn(inFile, ungrRef, ungrAlt)

    inFile = outdir + "/alt-genome.mis-split.ref-gene.txt"
    [misSplBln, misSplBlnRef] = findSplitBlastn(inFile, ungrRef, ungrAlt)

    inFile = outdir + "/alt-genome.mis-exon.ref-gene.txt"
    [misExon, misExRef,misExAlt, missingGrAlt]  = findExBlastn(inFile, ungrRef, ungrAlt)

    accBlastnRef = getAccBlastnRef(altBlnRefRes)
    grNum = getGroupNum(groupFile,ref)




    ''' step 2a: ungrouped Ref gene analysis '''

    fi = open(refBed,"r")
    refGenes = {}
    while True :
        line = fi.readline()
        if not line : break
        if not re.search("protein", line) : continue
        t = line.strip().split("\t")
        refGenes[t[3]] = 1
    fi.close()

    fo = open(outdir + "/ref-genome.ungrouped.gene.analysis.txt","w")
    fos = open(outdir + "/ref-genome.ungrouped.gene.analysis.stats","w")
    refUnknown={}
    ss = ["misMerBlp", "misMerBln", "misSplBlp", "misMerBln", "misExon","misUnGr","unAnn", "unAss", "paAss9090", "paAss9060", "paAss9060","cnv-NN","cnv-NM","unknown"]
    refGeneStat = {
            "misMerBlp":0, "misMerBln":0, "misSplBlp":0, "misMerBln":0, "misExon":0,"misUnGr":0,"unAnn":0, "unAss":0, "paAss9090":0, "paAss9060":0, "paAss9060":0,
            "cnv-NN":0,"cnv-NM":0,"unknown":0
            }
    for gene in sorted(refGenes.keys()) :
        flag = "xxx"
        if gene in misMerBlpRef.keys() and gene in misMerBlnRef.keys():
            flag = "misMer-BlpBln"
            refGeneStat["misMerBlp"] += 1
            refGeneStat["misMerBln"] += 1
        elif gene in misMerBlpRef.keys ()  :
            flag = "misMer-Blp"
            refGeneStat["misMerBlp"] += 1
        elif gene in misMerBlnRef.keys ():
            flag = "misMer-Bln"
            refGeneStat["misMerBln"] += 1

        elif gene in misSplBlpRef.keys() and gene in misSplBlpRef.keys ():
            flag = "misSpl-BlpBln"
            refGeneStat["misSplBlp"] += 1
            refGeneStat["misSplBln"] += 1
        elif gene in misSplBlpRef.keys() :
            flag = "misSpl-Blp"
            refGeneStat["misSplBlp"] += 1
        elif gene in misSplBlnRef.keys():
            flag = "misSpl-Bln"
            refGeneStat["misSplBln"] += 1

        elif gene in misExRef.keys():
            flag = "misExon"
            refGeneStat["misExon"] += 1

        elif gene in misGrEef.keys():
            flag = "misUnGr"
            refGeneStat["misUnGr"] += 1

        elif gene in unAssRef.keys() : ## deleted and partially deleted
            flag = unAssRef[gene]
            refGeneStat[flag] += 1

        elif gene in unAnnRef.keys(): # unannotated
            flag = "unAnn"
            refGeneStat[flag] += 1

        elif bestBlastp.has_key(gene) and bestBlastp[gene][0] > 80 and bestBlastp[gene][1] > 0.80  and bestBlastp[gene][2] > 0.80 :
            flag = "misUnGr"
            refGeneStat[flag] += 1
        elif gene in cnv.keys()  :
            flag = cnv[gene]
            n1 = int(flag.split("-")[0])
            n2 = int(flag.split("-")[1])
            if n1 == n2 :
                refGeneStat["cnv-NN"] += 1
            else :
                refGeneStat["cnv-NM"] += 1
        else:
            flag = "unknown"
            refUnknown[gene]=1
        bestp = "xx\txx\txx\txx"
        if bestBlastp.has_key(gene) :
            bestp = "\t".join(str(x) for x in bestBlastp[gene])
        fo.write("{}\t{}\t{}\t{}\t{}\n".format(refGeneBed[gene], gene, flag, refBest[gene], bestp))
    fo.close()

    for k in ss :
        fos.write("{}\t{}\n".format(k, refGeneStat[k]))
    fos.close()


    '''step 2b: analysis of ungrouped genes of query/alt genome '''
    altGeneStat = {
            "misMerBlp":0, "misMerBln":0, "misSplBlp":0, "misMerBln":0, "misExon":0,"misUnGr":0,"unAnn":0, "unAss":0, "paAss9090":0, "paAss9060":0, "paAss9060":0,
            "cnv-NN":0,"cnv-NM":0,"unknown":0
            }
    fo = open(outdir + "/alt-genome.ungrouped.gene.analysis.txt","w")
    fos = open(outdir + "/alt-genome.ungrouped.gene.analysis.stats","w")


    for gene in sorted(altGenes.keys()) :
        flag = "xxx"
        if gene in misMerBlpAlt.keys() and gene in misMerBlnAlt.keys():
            flag = "misMer-BlpBln"
            altGeneStat["misMerBlp"] += 1
            altGeneStat["misMerBln"] += 1
        elif gene in misMerBlpAlt.keys ()  :
            flag = "misMer-Blp"
            altGeneStat["misMerBlp"] += 1
        elif gene in misMerBlnAlt.keys ():
            flag = "misMer-Bln"
            altGeneStat["misMerBln"] += 1

        elif gene in misSplBlpAlt.keys() and gene in misSplBlpAlt.keys ():
            flag = "misSpl-BlpBln"
            altGeneStat["misSplBlp"] += 1
            altGeneStat["misSplBln"] += 1
        elif gene in misSplBlpAlt.keys() :
            flag = "misSpl-Blp"
            altGeneStat["misSplBlp"] += 1
        elif gene in misSplBlnAlt.keys():
            flag = "misSpl-Bln"
            altGeneStat["misSplBln"] += 1

        elif gene in misExAlt.keys():
            flag = "misExon"
            altGeneStat["misExon"] += 1

        elif gene in misGrAlt.keys():
            flag = "misUnGr"
            altGeneStat["misUnGr"] += 1

        elif gene in unAssAlt.keys() :
            flag = unAssAlt[gene]
            altGeneStat[flag] += 1
        elif gene in unAnnAlt.keys():
            flag = "unAnn"
            altGeneStat["unAnn"] += 1

        elif bestBlastp.has_key(gene) and bestBlastp[gene][0] > 80 and bestBlastp[gene][1] > 0.80  and bestBlastp[gene][2] > 0.80 :
            flag = "misUnGr"
            altGeneStat[flag] += 1
        elif gene in cnv.keys()  :
            flag = cnv[gene]
            n1 = int(flag.split("-")[0])
            n2 = int(flag.split("-")[1])
            if n1 == n2 :
                altGeneStat["cnv-NN"] += 1
            else :
                altGeneStat["cnv-NM"] += 1
        else:
            flag = "unknown"
            refUnknown[gene]=1

        bestp = "xx\txx\txx\txx"
        if bestBlastp.has_key(gene) :
            bestp = "\t".join(str(x) for x in bestBlastp[gene])
        fo.write("{}\t{}\t{}\t{}\t{}\n".format(altGeneBed[gene], gene, flag, altBest[gene], bestp))
    fo.close()

    for k in ss :
        fos.write("{}\t{}\n".format(k,altGeneStat[k]))
    fos.close()

    ''' step 3: final stats: '''
    #acc gene Ref-sp Acc-sp
    st1 = [0,0,0]
    st1[0] = len(altLen.keys())
    st1[1] = len(ungrRef.keys())
    st1[2] = len(ungrAlt.keys())
    #s-blastp s-blastn split m-blastp m-blastn merge ex-in acc.unann acc.unass acc.pa.ass col.unass ref.par.ass ref.unann # Ref-LoF Acc-LoF
    st2 = [0,0,0, 0,0,0, 0, 0,0,0, 0,0,0]
    st2[0] = len(altMisSplit.keys())
    st2[1] = len(misSplit.keys())
    tmp = misSplit
    tmp.update(altMisSplit)
    st2[2] = len(tmp.keys())

    st2[3] = len(altMisMer.keys())
    st2[4] = len(misMer.keys())
    tmp = misMer
    tmp.update(altMisMer)
    st2[5] = len(tmp.keys())

    st2[6] = len(misExon.keys())
    st2[7] = len(unAnnRef.keys())
    for x in unAssRef.keys() :
        if re.search(r"unass",unAssRef[x]) :
            st2[8] +=1
        else :
            st2[9] +=1
    for x in unAssAlt.keys() :
        if re.search(r"unass", unAssAlt[x]) :
            st2[10] +=1
        else :
            st2[11] +=1
    st2[12] = len(unAnnAlt.keys())

    #Ref-sp split merge ex-in acc.unann acc.unass acc.pa.ass
    st3 = [0, 0, 0, 0, 0, 0, 0 ]
    st3[0] = st1[1]
    for x in altMisSplitRef.keys() :
        if altMisSplitRef[x] == "ungrouped" :
            st3[1] += 1
    for x in misSplitRef.keys() :
        if not altMisSplitRef.has_key(x) and misSplitRef[x] == "ungrouped":
            st3[1] += 1


    for x in misMerGeneRef.keys() :
        if misMerGeneRef[x] == "ungrouped" :
            st3[2] += 1
    for x in misMerRef.keys() :
        if not misMerGeneRef.has_key(x) and misMerRef[x] == "ungrouped":
            st3[2] += 1
    st3[3] = misExRef
    for x  in unAnnRef.keys() :
        if ungrRef.has_key(x) :
            st3[4] += 1
    for x in unAssRef.keys() :
        if unAssRef[x] == "ungrouped-unass" :
            st3[5] += 1
        elif unAssRef[x] == "ungrouped-parass" :
            st3[6] += 1
    #Acc-sp split merge ex-in ref.unann ref.unass ref.pa.ass  ungr-rnaSup
    st4 = [0, 0, 0, 0, 0, 0, 0 ,0 ]
    st4[0] = st1[2]
    for x in altMisSplitAlt.keys() :
        if altMisSplitAlt[x] == "ungrouped" :
            st4[1] += 1
    for x in misMerGeneAlt.keys() :
        if misMerGeneAlt[x] == "ungrouped" :
            st4[2] += 1
    st4[3] = misExAlt
    for x in unAnnAlt.keys() :
        if ungrAlt.has_key(x) :
            st4[4] += 1
    for x in unAssAlt.keys() :
        if unAssAlt[x] == "ungrouped-unass" :
            st4[5] += 1
        elif unAssAlt[x] == "ungrouped-parass" :
            st4[6] += 1
    #st4[7] = ana2[9]
    for x in rnaSup :
        if ungrAlt.has_key(x) :
            st4[7]+=1
    fo = open(outdir + "/annotation.evaluation.stats","w")
    fo.write("all\t{}\t{}\n".format( "\t".join(str(x) for x in(st1)), "\t".join(str(x) for x in(st2)) ))
    fo.write("\n")
    fo.write("Ref-sp\t{}\tAcc-sp\t{}\n".format( "\t".join(str(x) for x in(st3)), "\t".join(str(x) for x in(st4)) ))






    '''step 4: combine the above information to get query gene to be updated or added '''
    fi = open(altBed,"r")
    outfile1 = outdir + "/alt-genes.to.be.updated.added.txt"
    fo = open(outfile1,"w")
    k0 = 0
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        id = t[3].split(";")[0]
        id = id.split("=")[1]
        blast = ";".join(accBlastnRef[id])

        flag = ""
        gr = "grouped"
        if ungrAlt.has_key(id) :
            gr = "ungrouped"
        else :
            gr = grNum[id]
        if misMer.has_key(id) :
            refGene = misMer[id][0]
            grRef = misMer[id][1]
            for k in range(len(misMer[id])) :
                if k == 0 : continue
                if k % 2 == 0 :
                    refGene = refGene + "," + misMer[id][k]
                    grRef = grRef + "," + misMer[id][k+1]
            flag = "mis-merge" + "\t" + refGene + "\t" + grRef
        elif altMisMer.has_key(id) :
            refGene = altMisMer[id][0]
            grRef = altMisMer[id][1]
            for k in range(len(altMisMer[id])) :
                if k == 0 : continue
                if k % 2 == 0 :
                    refGene = refGene + "," + altMisMer[id][k]
                    grRef = grRef + "," + altMisMer[id][k+1]
            flag = "mis-merge" + "\t" + refGene + "\t" + grRef

        elif misSplit.has_key(id) :
            refGene = misSplit[id][0]
            flag = "mis-split" + "\t" + refGene + "\t" + misSplit[id][1]

        elif altMisSplit.has_key(id):
            refGene = altMisSplit[id][0]
            flag = "mis-split" + "\t" + refGene + "\t" + altMisSplit[id][1]

        elif missingGrAlt.has_key(id) :
            refId = missingGrAlt[id]
            flag = "mis-ungr\t" + refId + "\t-"

        elif misExon.has_key(id) :
            #flag = "mis-ex-int"
            refGene = misExon[id][0]
            flag = "mis-exon" + "\t" + refGene + "\t" + misExon[id][2]

        elif unAssAlt.has_key(id) :
            if unAssAlt[id] == "ungrouped-unass" :
                flag = "Ref-Not-ass" + "\t-\t-"
            else :
                flag = "Ref-Par-ass" + "\t-\t-"

        elif unAnnAlt.has_key(id) :
            if rnaSup.has_key(id) :
                flag = "RefNo-RNAsup" + "\t" + ",".join(unAnnAlt[id]) + "\t-"
            else :
                flag = "Ref-No-ann" + "\t" + ",".join(unAnnAlt[id]) + "\t-"

        elif gtype.has_key(id) :
            flag = "Ref-diff-ann" + "\t" + gtype[id] + "\t-"

        #elif ungrAlt.has_key(id) :
        #    if accLoF.has_key(id) :
        #        flag = "Ref-LoF" + "\t"+ accLoF[id] + "\t-"
        #    else :
        #        flag = "Ref-No-ann" + "\t-\t-"
        #        k0 += 1
        else :
            flag = "unchange\t-\t-"

        fo.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(t[0],t[1],t[2],id,gr,blast,flag))
    fi.close()
    #print "ungrouped acc , ref no LoF and no ann, ", k0
    #not annotated in acc but in Ref
    missRefUngr = 0
    for k in unAnnRef :
        [chrom,start,end,prot,cov,iden] = unAnnRef[k]
        flag = "add"
        if ungrRef.has_key(k) :
            flag = flag + "\t" + k + "\tungrouped"
            missRefUngr += 1
        else  :
            flag = flag + "\t" + k + "\tgrouped"

        fo.write("{}\t{}\t{}\t{}\t{}\t{};{}\t{}\n".format(chrom,start,end,k,"**",cov,iden,flag))


    fi = open(outdir + "/ref-gene.blastn.intersect.alt-gene.txt","r")
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        refId = t[3].split(".")[0]
        if colUnknown.has_key(refId) :
            flag = "add"
            if ungrRef.has_key(refId) :
                flag = flag + "\t" + refId + "\tungrouped"
            else  :
                flag = flag + "\t" + refId + "\tgrouped"
            fo.write("{}\t{}\t{}\t{}\t{}\t{};{}\t{}\n".format(t[0],t[1],t[2],refId,"**",t[4],t[5],flag))
    fi.close()
    fo.close()
    outfile2 = outdir + "/alt-genes.to.be.updated.added.srt.txt"
    os.system("sort -k1,1 -k2,2n " + outfile1  + " > " + outfile2)

def getBesthit(inFile):
    besthit = {}
    fi = open(inFile,"r")
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        besthit[t[3]] = line.strip()
    fi.close()
    return besthit


def getRefLof(inFile):
    colLof = {}
    fi = open(inFile,"r")
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        if t[5] != "No" :
            colLof[t[0]] = t[5]
    fi.close()
    return colLof

def getAccLof(inFile):
    accLof = {}
    fi = open(inFile,"r")
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        if t[5] != "No" :
            accLof[t[1]] = t[5]
    fi.close()
    return accLof

def getBlastp(inFile,outFile, refLen, altLen):
    fi = open(inFile,"r")
    fo = open(outFile,"w")
    blastp = {}
    bestBlastp = {}
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        if t[2] < 60 : continue
        #t[1]=t[1].split(".")[0]
        cov1 = abs(float(t[7])-float(t[6]))/altLen[t[0]]
        cov2 = abs(float(t[9])-float(t[8]))/refLen[t[1]]
        if cov1 < 0.60 or cov2 < 0.60 : continue
        fo.write("{}\t{}\t{}\n".format(line.strip(), cov1, cov2))
        #refId = t[1].split(".")[0]
        refId = t[1]
        blastp[t[0] + refId] = [t[2], cov1, cov2]
        if not bestBlastp.has_key(t[0]) :
            bestBlastp[t[0]] = [refId, float(t[2]), cov1, cov2]
        else :
            if bestBlastp[t[0]][1]*bestBlastp[t[0]][2]*bestBlastp[t[0]][3] < float(t[2])*cov1*cov2 :
                bestBlastp[t[0]] = [refId, float(t[2]), cov1, cov2]

        if not bestBlastp.has_key(refId) :
            bestBlastp[refId] = [t[0], float(t[2]), cov1, cov2]
        else :
            if bestBlastp[refId][1]*bestBlastp[refId][2]*bestBlastp[refId][3] < float(t[2])*cov1*cov2 :
                bestBlastp[refId] = [t[0], float(t[2]), cov1, cov2]
    fi.close()
    return [blastp, bestBlastp]


def getAccBlastnRef (blastnResRef):
    fi = open(blastnResRef,"r")
    accBlastnRef = {}
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        #id = t[3].split(".")[0]
        id = t[3]
        #if re.search("evm",t[3]) :
        #    id = re.sub("model","TU",t[3])
        accBlastnRef[id] = [t[0],t[1],t[2],t[4],t[5]]
    fi.close()
    return accBlastnRef

def getGroupNum (groupFile,ref):
    fi = open(groupFile,"r")
    grNum = {}
    #ref = "AT"
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split()
        refIds = []
        altIds = []
        for k in range(len(t)) :
            if k == 0 : continue
            if re.search(ref,t[k]) :
                id = t[k].split(".")[0]
                refIds.append(id)
            else :
                id = t[k]
                altIds.append(id)
        for k in altIds :
            grNum[k] = str(len(refIds)) + "-" + str(len(altIds))
    fi.close()
    return grNum


def findMerBlastn(inFile, ungrRef, ungrAlt):
    fi = open(inFile,"r")
    misMer={}
    misMerRef = {}
    k1 = 0
    k2 = 0
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        #misMer[t[0]]=[]
        altId = t[9].split(";")[0]
        altId = altId.split("=")[1]
        #refId = t[3].split(".")[0]
        refId = t[3]
        gr = "grouped"
        if ungrRef.has_key(refId) :
            gr = "ungrouped"
            k1 += 1
        misMerRef[refId] = gr

        if ungrAlt.has_key(altId) :
            k2 += 1

        if not misMer.has_key(altId) :
            misMer[altId]=[refId,gr]
        else :
            misMer[altId].append(refId)
            misMer[altId].append(gr)
    fi.close()
    print "\n find mis-merged gene by blastn"
    print "mis-merged genes alt number: ", len(misMer)
    print "mis-merged genes ref number: ", len(misMerRef)
    print "mis-merged genes resulted in ungroupd Ref ", k1
    print "mis-merged genes resulted in ungroupd Alt ", k2
    return misMer, misMerRef

def findSplitBlastn(inFile, ungrRef, ungrAlt):
    fi = open(inFile,"r")
    misSplit = {}
    misSplitRef = {}
    k1 = 0
    k2 = 0
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        altId = t[9].split(";")[0]
        altId = altId.split("=")[1]
        #refId = t[3].split(".")[0]
        refId = t[3]
        misSplitRef[refId] = 1
        gr = "grouped"
        if ungrRef.has_key(refId) :
            gr = "ungrouped"
            k1 += 1
        if ungrAlt.has_key(altId) :
            k2 += 1
        misSplit[altId]=[refId,gr]
        #misSplit[x]=[t[4],xx2[1],t[6]]
    fi.close()
    print "\n # find mis-split gene by blastn "
    print "mis-split gene acc number :" , len(misSplit)
    print "mis-split gene result in ungrouped Ref gene :" , k1
    print "mis-split gene result in ungrouped Acc gene :" , k2
    return [misSplit,misSplitRef]


def findExBlastn(inFile, ungrRef, ungrAlt):
    fi = open(inFile,"r")
    misExon={}
    missingGrAlt = {}
    k = 0
    k1 = 0
    k2 = 0
    k31 = 0
    k32 = 0
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        #misMer[t[0]]=[]
        id = t[9].split(";")[0]
        id = id.split("=")[1]
        t[11] = re.sub("ref\-|alt\-","",t[11])
        t[12] = re.sub("ref\-|alt\-","",t[12])
        refId = t[3]
        if t[-1] == "not-ortholog" :
            misExon[id] = [refId,t[12],t[11]]
        if t[-1] == "missing-ortholog" :
            missingGrAlt[id] = refId
            if ungrRef.has_key(refId) :
                k31 += 1
            if ungrAlt.has_key(id) :
                k32 += 1
        else :
            if ungrAlt.has_key(id) :
                k2 += 1
            if ungrRef.has_key(refId) :
                k1 += 1
        #print id,t[12],t[11]
        #sys.exit()
    fi.close()
    print "\n # find gene with wrong exon-intron structure based on blastn"
    print "mis-exon gene acc number", len(misExon)
    print "mis-exon genes resulted in ungrouped Ref gene", k1
    print "mis-exon genes resulted in ungrouped Acc gene", k2

    print "mis-ungrouped Ref gene", k31
    print "mis-ungrouped Acc gene", k32
    print "len(missingGrAlt):", len(missingGrAlt)
    misExRef = k1
    misExAlt = k2
    return [misExon, misExRef,misExAlt, missingGrAlt]

def getBed(bedFile):
    fi = open(bedFile,"r")
    geneBed = {}
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        id = t[3].split(";")[0]
        id = id.split("=")[1]
        geneBed[id]="\t".join(t[0:3])
    fi.close()
    return geneBed



def findRefMM(blastnRes,refBed,altBed,outdir, ungrAlt, ungrRef, rnaGff, blastp):
    unAssAlt = {}
    unAnnAlt = {}
    geneBed = getBed(accBed)
    '''
      awk '{ if (/None|chl|mito/) print $2"\t"$7"\t"$8"\t"$1"\t"$9"\t"$10 ; else print "Chr"$2"\t"$7"\t"$8"\t"$1"\t"$9"\t"$10}'
      ../blastnRef/gene.blastn.besthit.out >query.prot.besthit.out2
    '''
    fi = open(blastnRes,"r")
    fo = open(outdir + "/ref-genome.un-assembled.alt-gene.txt", "w")
    [k100,k9090,k9060,k9030] = [0,0,0,0]
    [w,x,y,z] = [0,0,0,0]
    #k9030 iden>90, cov<30
    while True :
        line = fi.readline()
        if not line :break
        t = line.strip().split("\t")
        id = t[3]
        if t[0] == "None" :
        	unAssAlt[id] = "unAss"
        	k100 += 1
        	if ungrAlt.has_key(id) :
        		fo.write("{}\t{}\t{}\n".format(geneBed[id],id,"ungrouped-unass"))
        		w += 1
        	else :
        		fo.write("{}\t{}\t{}\n".format(geneBed[id],id,"grouped-unass"))

        elif float(t[5]) > 90 and float(t[4]) < 90 :
            unAssAlt[id] = 2
            k9090 += 1
            if float(t[4]) < 60 :
            	k9060 += 1
            if float(t[4]) < 30:
            	k9039 += 1
            	
            if ungrAlt.has_key(id) :
                fo.write("{}\t{}\t{}\n".format(geneBed[id],id,"ungrouped-parass"))
                x += 1
                if float(t[4]) < 60 :
                	y+=1
                if float(t[4]) < 30 :
                	z+=1	
            else :
                unAssAlt[id] = "grouped-parass"
                fo.write("{}\t{}\t{}\n".format(geneBed[id],id,"grouped-parass"))
    fi.close()
    fo.close()
    
    print("\n#findRefMM unassembled or partially assembled alt gene \n #del100\tI90C90\tI9060\tI9030\n")
    print("  {}\t{}\t{}\t{}\n".format(k100,k9090,k9060,k9030))
    print("  {}\t{}\t{}\t{}\n".format(w,x,y,z))

    blastnRes2 = outdir + "/alt-gene.blastn.Ref.besthit.out2"
    os.system("grep -v None " + blastnRes + " |sort -k1,1 -k2,2n -k3,3n > " + blastnRes2)
    os.system("intersectBed -a " + blastnRes2 + " -b " + refBed + " -wao |sort -k1,1 -k2,2n -k3,3n > " + outdir + "/alt-gene.blastn.intersect.ref-gene.txt")
    fi = open(outdir + "/alt-gene.blastn.intersect.ref-gene.txt", "r" )
    fo = open(outdir + "/ref-genome.missing.alt-gene.txt", "w" )

    unAnnAltSt1 = [0,0,0,0] # both group and ungroup
    unAnnAltSt2 = [0,0,0,0] # only ungrouped
    while True :
    	line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        id = t[3]
        fo.write("{}\t{}".format(geneBed[id],line.strip()))
        if float(t[5]) < 85: continue
        	
        if ungrAlt.has_key(id) :
            fo.write("\tungrouped\n")
            k += 1
            unAnnAlt[id] = t
            unAnnAltSt2[0] += 1
            if float(t[4]) > 30: 
            	unAnnAltSt2[1] += 1
            if float(t[4]) > 60: 
            	unAnnAltSt2[2] += 1
            if float(t[4]) > 90:
            	unAnnAltSt2[3] += 1	
        else :
            fo.write("\tgrouped\n")
            unAnnAlt[id] = t
            
        unAnnAltSt1[0] += 1
        if float(t[4]) > 30: 
        	unAnnAltSt1[1] += 1
        if float(t[4]) > 60: 
           	unAnnAltSt1[2] += 1
        if float(t[4]) > 90:
            unAnnAltSt1[3] += 1	
            	
    fi.close()
    fo.close()
    print "\n###"
    print "Gene annotated in accession and blasted to Ref assembly but not annotated in Ref :" , len(unAnnAlt)
    print "		unAnnAltSt1 " , unAnnAltSt1 
    print "		unAnnAltSt1 " , unAnnAltSt2
    
    
    
    os.system("intersectBed -a " + outdir + "/ref-genome.missing.alt-gene.txt -b " + rnaGff + " -wo > " + outdir + "/ref-genome.missing.alt-gene.intersectBed.RNA-seq.txt")
    fi = open(outdir + "/ref-genome.missing.alt-gene.RNA-seq.supported.txt", "r")
    rnaGene = {}
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        if not rnaGene.has_key(t[6]) :
            rnaGene[t[6]] = int(t[-1])
        else :
            rnaGene[t[6]] += int(t[-1])
    fi.close()

    fi = open(outdir + "/ref-genome.missing.alt-gene.intersectBed.RNA-seq.txt", "r" )
    fo = open(outdir + "/ref-genome.missing.alt-gene.RNA-seq.supported.txt", "w" )
    rnaSup = {}
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        flag = "xxx"
        if rnaGene.has_key(t[6]) :
            flag = float(rnaGene[t[6]])/(float(t[2]) - float(t[1]))
            if flag > 0.9 :
                rnaSup[t[6]] = 1
        fo.write("{}\t{}\n".format(line.strip(), flag))
    fi.close()
    fo.close()
    print "potential Ref missing gene with RNA-sup", len(rnaSup)

    os.system("intersectBed -a " + blastnRes2 + " -b " + refBed + " -wo > " + outdir + "/alt-gene.blastn.intersect.ref-gene.txt")
    ###if the paired gene is not ortholog, then the gene may have wrong exon-intron structure. Also check the protein exonerate alignment
    fi = open(outdir + "/alt-gene.blastn.intersect.ref-gene.txt", "r")
    fo = open(outdir + "/alt-gene.blastn.ref-ass.mis-exon.alt-gene.txt", "w")
    prevAltId = ""
    gtype = {}
    k1 = 0
    k2 = 0
    k3 = 0
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        altId = t[3]
        if not ungrAlt.has_key(altId) :
            continue
        if  float(t[-1])/(float(t[2]) - float(t[1])) < 0.8 or  float(t[-1])/(float(t[8]) - float(t[7])) < 0.8:
            continue

        refId = t[9]
        k1 += 1
        if not ungrRef.has_key(refId) :
            k2 += 1
        if blastp.has_key(t[3] + refId) :
            k3+=1
            continue

        if prevAltId == "" :
            prevAltId = altId
            fo.write(line)
            
        elif prevAltId != altId :
            prevAltId = altId
            fo.write(line)
        else :
            prevAltId = altId
    fi.close()
    fo.close()
    print "ungroupd Acc Genes annotated in accession and blast hit a gene in Ref assembly but they are not orthologs: " , k1
    print "ungrouped Acc Genes annotated in accession and blast hit a gene in Ref assembly, they are not orthologs, but the Ref has ortholog: " , k2
    print "ungroupded Genes annotated in accession ad blast hit a gene in Ref assembly, they are not orthologs, but they have high similarity: " , k3
    return [unAssAlt,unAnnAlt, rnaSup]

def getGeneIndex (geneBed):

    fi = open(geneBed,"r")
    geneIdx = {}
    n = 0
    genePos = {}
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        id = t[3].split(";")[0]
        id = id.split("=")[1]
        geneIdx[id]=n
        genePos[id]=[t[0],t[1],t[2]]
        n += 1
    fi.close()
    return [geneIdx,genePos]

    ### step 1: find mis-merging genes
def fltBlastp(blastpRes,outFile, refIdx,altIdx,altLen, refLen):
    fi = open(blastpRes,"r")
    fo = open(outFile,"w")
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        if float(t[2]) < 60 : continue
        aIdx =  altIdx[t[0]]
        rId = t[1]
        rIdx = refIdx[rId]
        len1 = altLen[t[0]]
        len2 = refLen[rId]
        cov1 = abs((float(t[7]) - float(t[6])))/float(len1)
        cov2 = abs((float(t[9]) - float(t[8])))/float(len2)
        fo.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(line.strip(), aIdx, rIdx, len1, len2, cov1, cov2))
    fi.close()
    fo.close()

def findMisMerBlp (blastpFlt,outdir,minId,maxSplitCov,maxIdf,refIdx,altIdx,altLen,altPos, refLen):

    blastpFltSrt = outdir + "/blastp.flt.srt1.out"
    os.system("sort -k13,13n -k14,14n -k7,7n " + blastpFlt + "  > " + outdir + "/blastp.flt.srt1.out") ## sorted by ref gene-idx

    fi = open(blastpFltSrt,"r")
    fo = open(outdir + "/potential.Qry.mis-merged.gene.by.blastp.txt","w")
    aIdx = 0
    rIdx = 0
    [prevIden,preCov] = [0,0]
    prevRefGene = ""
    prevAltEnd = 0
    minDist = 10
    misMerGene = {}
    misMerGeneRef = {}
    misMerGeneAlt = {}
    preLine = ""
    k1 =0
    k2 =0
    ref = "G"
    alt = "G"
    print "\n### find mis-merging genes in accession based on the blastp results "
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        if re.search(alt,t[0]) and re.search(ref,t[1]) :
            if float(t[2]) <= minId : continue

            aId = t[0]
            rId = t[1]

            if (aIdx == altIdx[aId] and refIdx[rId] == rIdx + 1 and abs(float(t[2])- prevIden) < maxIdf and preCov < maxSplitCov
                    and float(t[3])/altLen[aId] < maxSplitCov and abs(int(t[6]) - prevAltEnd) < minDist ) :
                [chrom,start,end]=altPos[aId]

                if aId not in misMerGene.keys():
                    misMerGene[aId]  = {}
                else :
                    misMerGene[aId][prevRefGene] = 1
                    misMerGene[aId][rId] = 1
                misMerGeneAlt[aId] = 1
                misMerGeneRef[prevRefGene] = aId
                misMerGeneRef[rId] = aId

                fo.write("{}\t{}\t{}\t{}\t{},{}\n".format(aId,chrom,start,end,prevRefGene,rId))
                fo.write("\t#{}\n\t#{}\n".format(preLine,line.strip()))
                aIdx = altIdx[aId]
                rIdx = refIdx[rId]
                prevId = float(t[2])
                preCov =  float(t[3])/altLen[aId]
                prevRefGene = rId
                prevAltEnd = int(t[7])
                preLine = line.strip()
            else :
                aIdx = altIdx[aId]
                rIdx = refIdx[rId]
                prevId = float(t[2])
                preCov =  float(t[3])/altLen[aId]
                prevRefGene = rId
                prevAltEnd = int(t[7])
                preLine = line.strip()

    fi.close()
    fo.close()
    print "mismerged gene acc : ",len(misMerGene)
    print "mismerged gene col : ",len(misMerGeneRef)
    return [misMerGene,misMerGeneRef,misMerGeneAlt]

    ## step 2: find mis-spliting genes
def findMisSplBlp (blastpFlt,outdir,minId,maxSplitCov,maxIdf,refIdx,altIdx,refLen,altPos):
    print "### find mis-split genes in accession based on the blastp results "
    os.system("sort  -k14,14n -k13,13n -k9,9n " + blastpFlt + "  > " + outdir + "/blastp.flt.srt2.out") ## sorted by alt gene-index
    fi = open(outdir + "/blastp.flt.srt2.out","r")
    fo = open(outdir + "/potential.Qry.mis-split.gene.by.blastp.txt","w")
    aIdx = -100
    rIdx = -100
    ref = "G"
    alt = "G"
    [prevIden,preCov] = [0,0]
    preAltGene = ""
    preRefEnd = 0
    altMisSplit = {}
    altMisSplitRef = {}
    altMisSplitAlt = {}
    k1 = 0
    k2 = 0
    #k1k2 = 0
    preLine = ""
    preAltCov = 0
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        if re.search(alt,t[0]) and re.search(ref,t[1]) :
            if float(t[2]) < minId : continue
            aId = t[0]
            rId = t[1]
            if (aIdx == altIdx[aId] - 1 and refIdx[rId] == rIdx and abs(float(t[2])- prevIden) < maxIdf and preCov < maxSplitCov and
                float(t[3])/refLen[rId] < maxSplitCov and  float(t[-2]) > maxSplitCov and preAltCov > maxSplitCov) :

                [chrom,start,end]=[altPos[aId][0],altPos[preAltGene][1],altPos[aId][2]]

                altMisSplitAlt[preAltGene] = rId
                altMisSplitAlt[aId] = rId
                if rId not in altMisSplitRef.keys() :
                    altMisSplitRef[rId] = 1
                else :
                    altMisSplitRef[rId] += 1
                altMisSplit[aId] = rId
                fo.write("{},{}\t{}\t{}\t{}\t{}\n".format(preAltGene,aId,chrom,start,end,rId))
                fo.write("\t#{}\n\t#{}\n".format(preLine,line.strip()))
                aIdx = altIdx[aId]
                rIdx = refIdx[rId]
                prevIden = float(t[2])
                preCov =  float(t[3])/refLen[rId]
                preAltGene = aId
                preAltCov = float(t[-2])
                preLine = line.strip()
            else :
                aIdx = altIdx[aId]
                rIdx = refIdx[rId]
                prevIden = float(t[2])
                preCov =  float(t[3])/refLen[rId]
                preAltGene = aId
                preAltCov = float(t[-2])
                preLine = line.strip()

    fi.close()
    fo.close()

    print "Blastp mis-split gene alt-genome :", len(altMisSplitAlt)
    print "Blastp mis-split gene ref-genome :", len(altMisSplitRef)
    print "\n"
    return [altMisSplit , altMisSplitRef, altMisSplitAlt]

def findAltMM (refBlnAltRes,altBed,outdir,ort,ungrRef,ungrAlt,refIdx,altIdx, blastp):
    ##blastnRes ref gene has a blastn besthit region in accession genome, but no gene was annotated
    unAnnRef = {} # unannotated in acc but annotated in ref
    misExRef = {}  # ref-gene has a blastn best-hit to a alt-gene (cov>90, iden>90), but they are not ortholog due to diff exon-intron structure
    misExAlt = {}
    unAssRef = {}  # genes assembled in ref but deleted or partially deleted in alt genome

    print "\n #findAltMM: find Ref genes which are not assembled (could be deleted) in other accession based on result of Ref-gene blastn against alt genome "

    fi = open(refBlnAltRes,"r")
    fo = open(outdir + "/alt-genome.un-assembled.ref-gene.txt","w")
    [k100,k9090,k9060,k9030] = [0,0,0,0]
    [w,x,y,z] = [0,0,0,0]
    #k9030 iden>90, cov<30
    while True :
        line = fi.readline()
        if not line :break
        t = line.strip().split("\t")
        id = t[3]
        if t[0] == "None" :
            flag = "ungrouped"
            if not ungrRef.has_key(id) :
                flag = "grouped"
            else :
                k100 += 1
            fo.write("{}\t{}\t{}\t{}\n".format(t[3],"un-assembled",flag, line.strip() ))
            unAssRef[id] = "unAss"
            w += 1
        elif float(t[5]) > 90 and float(t[4]) < 90 :
            x += 1
            flag = "ungrouped"
            unAssRef[id] = "paAss9090"
            if not ungrRef.has_key(id) :
                flag = "grouped"
                if float(t[4]) < 60 :
                    y += 1
                    k9060 += 1
                    unAssRef[id] = "paAss9060"
                if float(t[4]) < 30 :
                    z += 1
                    k9030 += 1
                    unAssRef[id] = "paAss9030"
            else :
                if float(t[4]) < 60 :
                    k9060 += 1
                    unAssRef[id] = "paAss9060"
                if float(t[4]) < 30 :
                    k9030 += 1
                    unAssRef[id] = "paAss9030"
            fo.write("{}\t{}\t{}\t{}\n".format(t[3],"partially-assembled",flag, line.strip() ))
            
    fi.close()
    fo.close()
    print("unassembled or partially assembled ref gene \n #del100\tI90C90\tI9060\tI9030\n")
    print("  {}\t{}\t{}\t{}\n".format(k100,k9090,k9060,k9030))
    print("  {}\t{}\t{}\t{}\n".format(w,x,y,z))

    blastnRes2 = outdir + "/ref.gene.blastn.besthit.bed"
    os.system("grep -v None " + refBlnAltRes + " > " + blastnRes2)

    os.system("intersectBed -a " + blastnRes2 + " -b " + altBed + " -wao > " + outdir + "/ref-gene.blastn.intersect.alt-gene.txt")
    os.system("intersectBed -a " + blastnRes2 + " -b " + altBed + " -wo |awk '{if ($6>90) print}' |cut -f 10 |sort -k1,1 |uniq -d -c |sed 's/ \+/\t/g' > " + outdir + "/alt-gene.multiple.blastn-hit.txt")
    os.system("intersectBed -a " + blastnRes2 + " -b " + altBed + " -wo |awk '{if ($6>90) print}' |cut -f 4 |sort -k1,1 |uniq -d -c  |sed 's/ \+/\t/g' > " + outdir + "/ref-gene.multiple.blastn-hit.txt")

    fi = open (outdir + "/ref-gene.multiple.blastn-hit.txt","r")
    repHit1 = {}
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        id = t[1]
        repHit1[id] = int(t[0])
    fi.close()
    print "\n###\nRef-genes have multiple (intersectBed >0) blastn-hit[alt-genes are mis-splitted] : ", len(repHit1)

    fi = open (outdir + "/alt-gene.multiple.blastn-hit.txt","r")
    repHit2 = {}
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        id = t[1].split(";")[0]
        id = id.split("=")[1]
        repHit2[id] = int(t[0])
    fi.close()
    print "Alt-genes have multiple (intersectBed >0) blastn-hit [alt-genes are mis-merged] : ", len(repHit2)

    fi = open( outdir + "/ref-gene.blastn.intersect.alt-gene.txt", "r")
    fo1 = open(outdir + "/alt-genome.un-annotated.ref-gene.txt", "w" )
    fo2 = open(outdir + "/alt-genome.mis-exon.ref-gene.txt", "w")
    fo3 = open(outdir + "/alt-genome.mis-split.ref-gene.txt", "w")
    fo4 = open(outdir + "/alt-genome.mis-merged.ref-gene.txt", "w")
    fo5 = open(outdir + "/alt-genome.m-vs-m.toBeChecked.txt", "w")
    fo6 = open(outdir + "/further.check.list","w")

    #repHit
    repHitM1 = {} ##ref-repeated  mis-split
    repHitM2 = {} ##query-repeated  mis-merge
    unAnnRefSt1 = [0,0,0,0] # both group and ungroup
    unAnnRefSt2 = [0,0,0,0] # only ungrouped
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        id = t[3]
        if t[6] == "." :
        	if float(t[5]) < 85: 
        		continue
        	fo1.write(line.strip())
        	if ungrRef.has_key(id) :
            	unAnnRef[id] = t[0:6]
            	fo1.write("\tungrouped\n")
            	unAnnRefSt2[0] += 1
            	
            	if float(t[4]) > 30: 
            		unAnnRefSt2[1] += 1
            	if float(t[4]) > 60: 
            		unAnnRefSt2[2] += 1
            	if float(t[4]) > 90:
            		unAnnRefSt2[3] += 1	
            else :
            	fo1.write("\tgrouped\n")
            unAnnRefSt1[0] += 1
            
            if float(t[4]) > 30: 
            	unAnnRefSt1[1] += 1
            if float(t[4]) > 60: 
            	unAnnRefSt1[2] += 1
            if float(t[4]) > 90:
            	unAnnRefSt1[3] += 1	
        else :
        	refId = t[3]
        	altId = t[9]
            if float(t[4]) < 80 or float(t[5]) < 90:
                continue
                
            if repHit1.has_key(refId) and repHit2.has_key(altId):
            	if repHitM1.has_key(refId) :
                	repHitM1[refId].append(t)
                else :
                	repHitM1[refId] = []
                	repHitM1[refId].append(t)
                gr1 = "ref-grouped"
                gr2 = "qry-grouped"
                if ungrRef.has_key(refId) :
                    gr1= "ref-ungrouped"
                if ungrAlt.has_key(altId) :
                    gr2= "alt-ungrouped"

                fo5.write("{}\t{}\t{}\n".format(line.strip(), gr1, gr2))

            elif repHit1.has_key(refId) :
                if repHitM1.has_key(refId) :
                    repHitM1[refId].append(t)
                else :
                    repHitM1[refId] = []
                    repHitM1[refId].append(t)

            elif repHit2.has_key(altId) :
                if repHitM2.has_key(altId) :
                    repHitM2[altId].append(t)
                else :
                    repHitM2[altId] = []
                    repHitM2[altId].append(t)
            else :
                if ort.has_key(refId) and (ort[refId].has_key(altId)) :
                    k1 += 1
                    continue
                if ort.has_key(refId) and (not ort[refId].has_key(altId)) :

                    flag1 = "ref-grouped"
                    flag2 = "alt-grouped"
                    if ungrRef.has_key(refId) :
                        flag1 = "ref-ungrouped"

                    if ungrAlt.has_key(altId) :
                        flag2 = "alt-ungrouped"

                    fo2.write("\t".join(t[0:11]))

                    if blastp.has_key(altId + refId) and blastp[altId + refId][0] > 80 and blastp[altId + refId][1] > 0.80 and blastp[altId + refId][2] > 0.80:
                        fo2.write("\t{}\t{}\t{}\n".format(flag1,flag2,"missing-ortholog"))
                    else :
                        fo2.write("\t{}\t{}\t{}\n".format(flag1,flag2,"not-ortholog"))
                        misExRef[refId] = 1
                        misExAlt[altId] = 1
                elif not ort.has_key(refId) :
                    #fo2.write(line.strip())
                    flag1 = "ref-grouped"
                    flag2 = "alt-grouped"
                    if ungrRef.has_key(refId) :
                        flag1 = "ref-ungrouped"

                    if ungrAlt.has_key(altId) :
                        flag2 = "alt-ungrouped"

                    fo2.write("\t".join(t[0:11]))

                    if blastp.has_key(altId + refId) and blastp[altId + refId][0] > 80 and blastp[altId + refId][1] > 0.80 and blastp[altId + refId][2] > 0.80:
                        fo2.write("\t{}\t{}\t{}\n".format(flag1,flag2,"missing-ortholog"))
                    else :
                        fo2.write("\t{}\t{}\t{}\n".format(flag1,flag2,"not-ortholog"))
                        misExRef[refId] = 1
                        misExAlt[altId] = 1

                elif float(t[4]) == 100.0 and float(t[5]) == 100.0 and (t[1] != t[7] or t[2] != t[8] ):
                    flag1 = "ref-grouped"
                    flag2 = "alt-grouped"
                    if ungrRef.has_key(refId) :
                        flag1 = "ref-ungrouped"

                    if ungrAlt.has_key(altId) :
                        flag2 = "alt-ungrouped"

                    fo2.write("\t".join(t[0:11]))

                    if blastp.has_key(altId + refId) and blastp[altId + refId][0] > 80 and blastp[altId + refId][1] > 0.80 and blastp[altId + refId][2] > 0.80:
                        fo2.write("\t{}\t{}\t{}\n".format(flag1,flag2,"missing-ortholog"))
                    else :
                        fo2.write("\t{}\t{}\t{}\n".format(flag1,flag2,"not-ortholog"))
                        misExRef[refId] = 1
                        misExAlt[altId] = 1
    fi.close()
    print "\n### Genes annotated in ref but not in alt, blastn besthit (iden>85) intersect with Al-genome"
    print " Both [cov0, cov30, cov60, cov90]  ", unAnnRefSt1
    print " ungr [cov0, cov30, cov60, cov90], ", unAnnRefSt2


    for k in sorted(repHitM1.keys()) :
        if len(repHitM1[k]) == 2 :
            altId1 = repHitM1[k][0][9].split(";")[0]
            altId1 = altId1.split("=")[1]
            altId2 = repHitM1[k][1][9].split(";")[0]
            altId2 = altId2.split("=")[1]
            if ( float(repHitM1[k][0][8]) - float(repHitM1[k][0][7]) ) * 0.9 < float(repHitM1[k][0][-1]) \
                and ( float(repHitM1[k][1][8]) - float(repHitM1[k][1][7]) ) * 0.9 < float(repHitM1[k][1][-1]) and abs(altIdx[altId1] - altIdx[altId2]) == 1 :
                t1 = repHitM1[k][0]
                t2 = repHitM1[k][1]
                #fo3.write("\t".join(t1[0:11]))
                fo3.write("\t".join(t1[0:10]))
                fo3.write("\n")
                fo3.write("\t".join(t2[0:10]))
                fo3.write("\n")
            else :
                fo6.write("\t".join(repHitM1[k][0]))
                fo6.write("\t".join(repHitM1[k][1]))
                fo6.write("\n")
        elif len(repHitM1[k]) > 2:
            flag = 0

            for i in range(len(repHitM1[k])) :
                if ( float(repHitM1[k][i][8]) - float(repHitM1[k][i][7]) ) * 0.9 < float(repHitM1[k][0][-1]) :
                    flag = 1
                else :
                    flag = 0
                    break
            if flag == 1 :
                for i in range(len(repHitM1[k])) :
                    t = repHitM1[k][i]
                    fo3.write("\t".join(t[0:10]))
                    fo3.write("\n")
            else :
                for i in range(len(repHitM1[k])) :
                    t = repHitM1[k][i]
                    fo6.write("\t".join(t[0:10]))
                    fo6.write("\n")
        else :
            fo6.write("\t".join(repHitM1[k][0]))
            fo6.write("\n")




    for k in sorted(repHitM2.keys()) :
        if len(repHitM2[k]) == 2 :
            if ( float(repHitM2[k][0][2]) - float(repHitM2[k][0][1]) ) * 0.9 < float(repHitM2[k][0][-1]) \
                and ( float(repHitM2[k][0][8]) - float(repHitM2[k][0][7]) ) * 0.9 < float(repHitM2[k][0][-1]) \
                and float(repHitM2[k][1][8]) - float(repHitM2[k][1][7]) *0.9 > float(repHitM2[k][1][-1]) :
                ## two genes annotated in one region, may be +/- or same strand (should not be wrong annotation??)
                t = repHitM2[k][1]
                unAnnRef[id] = t[0:6]
                fo1.write("\t".join(t[0:10]))
                fo1.write("\n")
            elif ( float(repHitM2[k][1][2]) - float(repHitM2[k][1][1]) ) * 0.9 < float(repHitM2[k][1][-1]) \
                and ( float(repHitM2[k][1][8]) - float(repHitM2[k][1][7]) ) * 0.9 < float(repHitM2[k][1][-1]) \
                and float(repHitM2[k][0][8]) - float(repHitM2[k][0][7]) *0.9 > float(repHitM2[k][0][-1]) :
                ## two genes annotated in one region, may be +/- or same strand (should not be wrong annotation??)
                t = repHitM2[k][0]
                unAnnRef[id] = t[0:6]
                fo1.write("\t".join(t[0:10]))
                fo1.write("\n")
            else  :
                refId1 = repHitM2[k][0][3]
                refId2 = repHitM2[k][1][3]
                if abs(refIdx[refId1] - refIdx[refId2]) == 1 and  repHitM2[k][0][1] != repHitM2[k][1][1]:
                    t1 = repHitM2[k][0]
                    t2 = repHitM2[k][1]
                    #if k == "ATAN1-1G16270" :
                    #    print t1,"\n"
                    #    print t2,"\n"
                    fo4.write("\t".join(t1[0:10]))
                    fo4.write("\n")
                    fo4.write("\t".join(t2[0:10]))
                    fo4.write("\n")
                else :
                    if  float(repHitM2[k][0][4]) == 100.0 and float(repHitM2[k][0][5]) == 100.0 and ( repHitM2[k][0][1] != repHitM2[k][0][7] or repHitM2[k][0][2] != repHitM2[k][0][8] ) :
                        flag1 = "ref-grouped"
                        flag2 = "alt-grouped"
                        flag3 = "not-orthlog"
                        refId = repHitM2[k][0][3]
                        altId = k
                        t = repHitM2[k][0]
                        fo2.write("\t".join(t[0:10]))
                        if ort.has_key(refId) and ort[refId].has_key(k) :
                            flag3 = "orthlog"
                        if ungrRef.has_key(refId) :
                            flag1 = "ref-ungrouped"

                        if ungrAlt.has_key(altId) :
                            flag2 = "alt-ungrouped"

                        if blastp.has_key(altId + refId) and blastp[altId + refId][0] > 80 and blastp[altId + refId][1] > 0.80 and blastp[altId + refId][2] > 0.80:
                            fo2.write("\t{}\t{}\t{}\n".format(flag1,flag2,"missing-orthloog"))
                        else :
                            fo2.write("\t{}\t{}\t{}\n".format(flag1,flag2,"not-ortholog"))
                            misExRef[refId] = 1
                            misExAlt[altId] = 1
                        #fo2.write("\t{}\t{}\t{}\n".format(flag1,flag2,flag3))
                    elif  float(repHitM2[k][1][4]) == 100.0 and float(repHitM2[k][1][5]) == 100.0 and ( repHitM2[k][1][1] != repHitM2[k][1][7] or repHitM2[k][1][2] != repHitM2[k][1][8] ) :
                        flag1 = "ref-grouped"
                        flag2 = "alt-grouped"
                        flag3 = "not-orthlog"
                        #refId = repHitM2[k][0][3].split(".")[0]
                        refId = repHitM2[k][0][3]
                        altId = k
                        t = repHitM2[k][1]
                        fo2.write("\t".join(t[0:10]))
                        if ort.has_key(refId) and ort[refId].has_key(k) :
                            flag3 = "orthlog"
                        if ungrRef.has_key(refId) :
                            flag1 = "ref-ungrouped"

                        if ungrAlt.has_key(altId) :
                            flag2 = "alt-ungrouped"

                        if blastp.has_key(altId + refId) and blastp[altId + refId][0] > 80 and blastp[altId + refId][1] > 0.80 and blastp[altId + refId][2] > 0.80:
                            fo2.write("\t{}\t{}\t{}\n".format(flag1,flag2,"missing-ortholog"))
                        else :
                            fo2.write("\t{}\t{}\t{}\n".format(flag1,flag2,"not-ortholog"))
                            misExRef[refId] = 1
                            misExAlt[altId] = 1
                        #fo2.write("\t{}\t{}\t{}\n".format(flag1,flag2,flag3))


        elif len(repHitM2[k]) > 2:
            flag = 0
            for i in range(len(repHitM2[k])) :
                refId = repHitM2[k][i][3]
                #print refId
                if ort.has_key(refId) and ort[refId].has_key(k) and (float(repHitM2[k][i][8]) - float(repHitM2[k][i][7]) ) * 0.8 < float(repHitM2[k][i][-1]):
                    flag = 1  ### multiple vs one, some deleleted in query genome
            if int(repHitM2[k][1][1]) < int(repHitM2[k][0][2]) and  int(repHitM2[k][1][1]) - int(repHitM2[k][0][1]) <  int(repHitM2[k][0][2]) - int(repHitM2[k][1][1]) :
                flag = 1

            if flag == 0 :
                for i in range(len(repHitM2[k])) :
                    t = repHitM2[k][i]
                    fo4.write("\t".join(t[0:10]))
                    fo4.write("\n")
            else :
                for i in range(len(repHitM2[k])) :
                    t = repHitM2[k][i]
                    fo6.write("\t".join(t[0:10]))
                    fo6.write("\n")
        else :
            fo6.write("\t".join(repHitM2[k][0]))
            fo6.write("\n")
    fo1.close()
    fo2.close()
    fo3.close()
    fo4.close()
    fo5.close()
    fo6.close()
    print "\n#false mis-ex as they are orthlog: ",  k1
    print "\n"
    return [unAnnRef,misExRef,misExAlt,unAssRef]

def getAltSpeGene(ungrAlt,altPos,outfile):
    fo = open(outfile,"w")
    #print len(ungrAlt.keys())
    for k in sorted(ungrAlt.keys()):
        [chrom,start,end] = altPos[k]
        fo.write("{}\t{}\t{}\t{}\n".format(chrom,start,end,k))
        #sys.exit()
    fo.close()


def getProtLeng2(protFas):
    fx = open(protFas,"r")
    leng = 0
    l = {}
    id = ""
    while True :
        line = fx.readline()
        if not line : break
        if line[0]==">" :
            if id != "" :
                id = id.split("|")[1]
                if re.search("evm",id) :
                    id = re.sub("model","TU",id)
                else :
                    id = id.split(".")[0]
                l[id] = leng
                leng = 0
            t = line.strip().split()
            id = t[0]
        else :
            leng += len(line.strip())
    id = id.split("|")[1]
    if re.search("evm",id) :
        id = re.sub("model","TU",id)
    else :
        id = id.split(".")[0]
    l[id] = leng
    fx.close()
    return l

def getProtLeng(protFas):
    l = {}
    print(protFas)
    if not os.path.isfile(protFas) :
        print "no such file ", protFas
        sys.exit()
    for seq_record in SeqIO.parse(protFas, "fasta"):
        #print(seq_record.id)
        #print(repr(seq_record.seq))
        #print(len(seq_record))
        l[seq_record.id] = len(seq_record)

    return l


def getSpecGene(groupFile,protRef,protAlt,ref,alt):
    fx = open(groupFile,"r")
    gr = {}
    ort = {}
    cnv = {}
    while True :
        line = fx.readline()
        if not line : break
        if not re.search(ref, line) or not re.search(alt,line) : continue
        t = line.strip().split()
        refIds = []
        altIds = []
        for k in range(len(t)) :
            if k == 0 : continue
            if re.search(ref,t[k]) :
                #id = t[k].split(".")[0]
                id = t[k]
                refIds.append(id)
                gr[id]=1
            else :
                #id = t[k].split(".")[0]
                id = t[k]
                altIds.append(id)
                gr[id]=1
        for k in refIds :
            ort[k]={}
            for j in altIds :
                ort[k][j]=1
            cnv[k] = str(len(refIds)) + "-" + str(len(altIds))
        for k in altIds
            cnv[k] = str(len(refIds)) +"-" + str(len(altIds))
    fx.close()
    ungrRef = {}
    ungrAlt = {}

    for k in sorted(protRef.keys()) :
        if not gr.has_key(k) :
            ungrRef[k]=1
    for k in sorted(protAlt.keys()) :
        if not gr.has_key(k) :
            ungrAlt[k]=1
    print "ungrouped ref gene ", len(ungrRef)
    print "ungrouped alt gene ", len(ungrAlt)
    return [ungrRef,ungrAlt,ort,cnv]


if __name__ == "__main__":
   main(sys.argv[1:])



