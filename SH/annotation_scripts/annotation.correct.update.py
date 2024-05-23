#!/usr/bin/env python
# encoding: utf-8


import sys
import os
import getopt
import re
from Bio import SeqIO
import multiprocessing
import datetime

def main(argv):
    updateFile = ""
    outdir = ""
    gffFile = ""
    sciFile = ""
    augFile = ""
    snapFile = ""
    glimFile = ""
    refBlnFile = ""
    genomeFile = ""
    refProtFile = ""
    refGeneBed = ""
    try:
        opts, args = getopt.getopt(argv,"u:g:o:s:a:n:l:b:f:p:i:",["update=","gff=","outd=","scipio=","augEvi=","glim=","snap=","refBln=","genome=","prot=","refIdx="])
    except getopt.GetoptError:
        print 'update.misann.genes.py -u <update> -g <gff> -o <outd> -s <scipio> -a <augEvi> -n <snap> -l <glim> -b <refBln> -f <genome> -p <prot> -i <refIdx>'
        sys.exit(2)
    if len(opts) == 0 :
        print 'update.misann.genes.py -u <update> -g <gff> -o <outd> -s <scipio>  -a <augEvi> -n <snap> -l <glim> -b <refBln> -f <genome> -p <prot> -i <refIdx> '
        sys.exit()
    for opt, arg in opts:
        if opt == '-h':
            print 'update.misann.genes.py -u <update> -g <gff> -o <outd> -s <scipio>  -a <augEvi> -n <snap> -l <glim> -b <refBln> -f <genome> -p <prot> -i <refIdx>'
            sys.exit()
        elif opt in ("-u", "--update"):
            updateFile = arg
        elif opt in ("-g", "--gff"):
            gffFile = arg
        elif opt in ("-o", "--outd"):
            outdir = arg
        elif opt in ("-s", "--scipio"):
            sciFile = arg
        elif opt in ("-a", "--augEvi"):
            augFile = arg
        elif opt in ("-n", "--snap") :
            snapFile = arg
        elif opt in ("-l", "--glim") :
            glimFile = arg
        elif opt in ("-b", "--refBln") :
            refBlnFile = arg
        elif opt in ("-f", "--genome") :
            genomeFile = arg
        elif opt in ("-p", "--prot") :
            refProtFile = arg
        elif opt in ("-i","--refIdx") :
            refGeneBed = arg

    startTime = datetime.datetime.now()
    print "Starting ======  ",startTime.strftime('%Y-%m-%d %H:%M:%S'), "  ======"

    ## gene model from augustus protein-based, scipio protein-based or stringTie-gff based
    if not os.path.isdir(outdir) :
        os.system("mkdir -p " + outdir)

    ## step 1: add update-action information
    scigene = getScipio(sciFile)
    outfile1 = outdir + "/genes.to.be.updated.txt2"
    checkUpdate(scigene,updateFile,outfile1)

    ## step 2: update gff file
    [augGenes,augLines] = getAugEviGene(augFile)
    [snapGenes,snapLines] = getAbiGene(snapFile)
    [glimGenes,glimLines] = getAbiGene(glimFile)

    blastn = getBlastn(refBlnFile)
    assSeq = getFasta(genomeFile)
    refProt = getFasta(refProtFile)
    print "start updating..."
    augCfg="/netscratch/dep_mercier/grp_schneeberger/bin/augustus/augustus-3.2.3/config/extrinsic/extrinsic.MPE.cfg"
    [geneIdx,genePos] = getGeneIndex(refGeneBed)
    genesAnn =[scigene,augGenes,augLines,snapGenes,snapLines,glimGenes,glimLines]

    upFiles = splitGff(outfile1, outdir)
    print len(upFiles)
    # pool = multiprocessing.Pool(processes = len(upFiles))
    # pool = multiprocessing.Pool(processes = 1) # Change MG
    # Change by MG= removed parallelisation as that was resulting in errors. Can be properly added later if required.
    k = 0
    res = []
    print upFiles
    for upFile in upFiles:
    # for upFile in upFiles[:1] :
        print upFile,"\n"
        outfile2 = outdir + "/updated.chunk." +str(k) + ".gff"
        # upSt = pool.apply_async(updateGff,(gffFile,upFile, outfile2, genesAnn, outdir,blastn,assSeq,augCfg,refProt,geneIdx,genePos,))
        upSt = updateGff(gffFile,upFile, outfile2,genesAnn,outdir,blastn,assSeq,augCfg,refProt,geneIdx,genePos)
        res.append(upSt)
        k += 1
    # pool.close()
    # pool.join()
    print 'Finished UpdateGFF'
    print res
    st1 = [0,0,0,0,0]
    fo = open(outdir + "/update.stats","w")
    fo.write("\t".join(["unchange","mis-merge","mis-split","mis-exon","add",]))
    for i in range(len(res)) :
        fo.write("{}\n".format("\t".join(str(x) for x in res[i][0]) ) )
        for k in range(len(st1)) :
            t = res[i][0]
            st1[k] += t[k]
    fo.write("{}\n".format("\t".join(str(x) for x in st1 ) ) )
    fo.write("\n")
    st2 = [0,0,0,0,0]
    for i in range(len(res)) :
        fo.write("{}\n".format("\t".join(str(x) for x in res[i][1]) ) )
        for k in range(len(st2)) :
            t = res[i][1]
            st2[k] += t[k]
    fo.write("{}\n".format("\t".join(str(x) for x in st2 ) ) )
    fo.close()

    outfile2 = outdir + "/updated.gff"
    os.system("cat " + outdir + "/updated.chunk.*.gff > " + outfile2)
    outfile3 = outdir + "/updated.rmdup.gff"
    outfile4 = outdir + "/updated.highConf.gff"
    outfile5 = outdir + "/updated.highConf.prot.fasta"
    rmDup(outfile2,outfile3)
    os.system("grep -v Low " + outfile3 + " > " + outfile4 )
    print("/srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl " + outfile4 + "  " + genomeFile + " prot > " + outfile5)

    endTime = datetime.datetime.now()
    print "Finished extracting reads. Time elapsed: ",(endTime - startTime).seconds, " s"

def splitGff (upFile, outdir):
    fi = open(upFile,"r")
    chrom = ""
    k = 0
    upFiles = []
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        if chrom == "" :
            fo = open(outdir + "/chunk." + str(k) + ".gff", "w")
            upFiles.append(outdir + "/chunk." + str(k) + ".gff")
            fo.write(line)
            chrom = t[0]
        elif chrom != t[0] :
            fo.close()
            k += 1
            chrom = t[0]
            fo = open(outdir + "/chunk." + str(k) + ".gff", "w")
            upFiles.append(outdir + "/chunk." + str(k) + ".gff")
            fo.write(line)
        else :
            fo.write(line)
    fi.close()
    fo.close()
    return upFiles


def getColLof(inFile):
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

def rmDup(gff,out):
    fi = open(gff,"r")
    fo = open(out,"w")
    flag = 1
    genes = {}
    dup = 0
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        #id = t[8].split(";")[0]
        #id = id.split("=")[1]
        if t[2]=="gene" :
            id =  "_".join([t[0],t[3],t[4],t[6]])
            if not genes.has_key(id) :
               fo.write(line)
               flag = 1
               genes[id]=1
            else :
                flag = 0
                dup += 1

        else :
            if flag == 1 :
                fo.write(line)

    fi.close()
    fo.close()
    print "remove duplcated ones : ", dup



def getFasta(fasta):
    fi = open(fasta,"r")
    id = ""
    seq = {}
    k = 0
    for record in SeqIO.parse(fasta, "fasta"):
        id = record.id
        #if re.search(".", id) : id = id.split(".")[0]
        seq[id] = record.seq
    print "refer seq:",len(seq)

    return seq

def getBlastn (file):
    fi = open(file,"r")
    blastn = {}
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        id = t[3]
        blastn[id] = [t[0],t[1],t[2]]
    fi.close()
    print "blastn:",len(blastn),"\n"
    return blastn



def updateGff(gffFile,updateFile,outFile,genesAnn,outdir,blastn,refSeq,augCfg,refProt,geneIdx,genePos):

    [scigene,augGenes,augLines,snapGenes,snapLines,glimGenes,glimLines] = genesAnn
    fi = open(gffFile,"r")
    gff = {}
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        if t[2] == "gene" :
            id = t[8].split(";")[0]
            id = id.split("=")[1]
            id = re.sub(r'evm.TU.',"",id)
            gff[id] = []
            gff[id].append(line.strip())
        else :
            gff[id].append(line.strip())
    fi.close()

    tempdir = os.path.join(outdir,"temp")
    if not os.path.isdir(tempdir) :
        os.mkdir(tempdir)

    fi = open(updateFile,"r")
    fo = open(outFile,"w")
    genesUp = ""
    checkGenes = {}
    prevStart = ""
    prevId = ""
    prevRefId = ""
    prevGr = ""
    misSplit = {}
    upStats = [[0,0,0,0,0], [0,0,0,0,0]] #unchange mis-split mis-merge mis-ex add
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        #print "\n## update : ", line.strip()
        if t[6] == "mis-merge" :  ##mis-merge
            upStats[0][1] += 1
            print "## mis-merge ",line.strip()
            genes = t[7].split(",")
            cc = t[-1].split(",")
            #if len(cc)!=len(genes) :
            #    print "** wrong ", line
            #    sys.exit()
            n = 0
            isFind = 0
            for i in range(len(genes)) :
                refId = genes[i]
                if checkGenes.has_key(refId) : continue
                checkGenes[refId] = 1
                [blastnCh,blastnStart,blastnEnd] = ["","",""]
                if blastn.has_key(refId) :
                    [blastnCh,blastnStart,blastnEnd] = blastn[refId]
                    if blastnCh != t[0] :
                        print " #diff ch :",blastnCh,refId,t[0], "\n"
                else :
                    print " #no blastn :",refId,"\n"

                if cc[i] == "unchange" :
                    gffline = "\n".join(gff[t[3]])
                    fo.write(gffline)
                    fo.write("\n")
                    isFind += 1
                    print line.strip()," #unchange: " , t[3]
                elif cc[i] == "changeSci" or cc[i] == "changeSciAug":
                    flag = checkScipio(scigene[refId],refSeq)
                    if flag == "good" :
                        lines = getSciGene(scigene[refId])
                        if lines != "" :
                            fo.write("{}\n".format(lines))
                            isFind += 1
                            print  " #updated based on Scipio: ", cc[i]
                        else :
                            print  " #updated based on Scipio [blank lines]: ", cc[i]
                    else :
                        prevEnd = int(blastn[genes[0]][2]) - 100
                        if i >0 :
                            prevEnd =  int(blastn[genes[i-1]][1])
                        nextStart = int(blastnEnd) + 200

                        if i <  len(genes) - 1 and blastn[genes[i+1]][0] == t[0] :
                            nextStart =  int(blastn[genes[i+1]][1])
                        elif i == len(genes) -1 :
                            idx = genePos[refId]
                            if idx+1 >= len(geneIdx) :
                                print "  ", refId, "no blastn of next gene:"
                            elif blastn.has_key(geneIdx[idx+1]) :
                                nextG = geneIdx[idx+1]
                                nextStart = blastn[nextG]
                            else :
                            	nextG = geneIdx[idx+1]
                                print "  ", refId, "no blastn of next gene:",nextG
                        if (not augGenes.has_key(t[0])) or (not snapGenes.has_key(t[0])) or (not glimGenes.has_key(t[0])) :
                            lines = ""
                            print "  no gene was predicted in the contig :", t[0]
                        else :
                            print "prevEnd ", prevEnd, " nextstart ", nextStart
                            #lines = checkAugAbiGene(int(t[1]),int(t[2]),augGenes[t[0]],augLines[t[0]],snapGenes[t[0]],snapLines[t[0]],glimGenes[t[0]],glimLines[t[0]] )
                            lines = checkAugAbiGene2(int(blastnStart),int(blastnEnd),augGenes[t[0]],augLines[t[0]],snapGenes[t[0]],snapLines[t[0]],glimGenes[t[0]],glimLines[t[0]], prevEnd,nextStart)
                        if lines != "" :
                            fo.write("{}\n".format(lines))
                            isFind += 1
                            print " updated based on ab initio prediction"

                        else :
                            n += 1
                            protFile = tempdir + "/" + refId + ".prot.fa"
                            foP = open(protFile,"w")
                            protFa = refProt[refId]
                            foP.write(">" + refId + "\n")
                            foP.write(str(protFa) + "\n")
                            foP.close()

                            faFile = tempdir + "/" + refId + ".genome.fa"
                            [tch,start,end,strand]=flag.split("\t")
                            if blastn[genes[i]][0] == t[0]:
                                start = blastnStart
                                end = blastnEnd
                                strand = "both"

                                foF = open(faFile,"w")
                                subseq = refSeq[t[0]][int(start)- 200 -1 : int(end) + 200]
                                foF.write(">" + t[0] + "_" + str(start) + "_" + str(end) + "\n")
                                foF.write(str(subseq) + "\n")
                                foF.close()

                                lines = getGenewise(faFile, protFile, tempdir, augCfg, strand, start)
                                if lines != "" :
                                    isFind += 1
                                    fo.write("{}\n".format(lines))
                                    print " updated based on Genewise", t[3]
                                else :
                                    print "  no gene\t", t[3]
                            else :
                                print "  NoGene\t",tch==t[0],t[0]==blastn[genes[i]][0],"\t",line
                        #if refId == "AT1G10920" or refId == "AT1G10930":
                        #    print "lines:", lines
                        #    sys.exit()
                elif cc[i] == "changeAug" :
                    if blastnEnd == "" :
                        gffline = "\n".join(gff[t[3]])
                        fo.write(gffline)
                        fo.write("\n")
                        isFind += 1
                        print "  unchange: " , t[3]
                        continue
                    nextStart = int(blastnEnd) + 200
                    prevEnd = int(blastn[genes[0]][2]) - 100
                    if i > 0 : prevStart = int(blastn[genes[i-1]][1])
                    if i <  len(genes) - 1 and blastn[genes[i+1]][0] == t[0] :
                        nextStart =  int(blastn[genes[i+1]][1])
                    if (not augGenes.has_key(t[0])) or (not snapGenes.has_key(t[0])) or (not glimGenes.has_key(t[0])) :
                        lines = ""
                        print "  no gene was predicted in the contig :", t[0]
                    else :
                        #lines = checkAugAbiGene(int(t[1]),int(t[2]),augGenes[t[0]],augLines[t[0]],snapGenes[t[0]],snapLines[t[0]],glimGenes[t[0]],glimLines[t[0]] )
                        lines = checkAugAbiGene2(int(blastnStart),int(blastnEnd),augGenes[t[0]],augLines[t[0]],snapGenes[t[0]],snapLines[t[0]],glimGenes[t[0]],glimLines[t[0]] ,prevEnd, nextStart)
                    if lines != "" :
                        isFind += 1
                        fo.write("{}\n".format(lines))
                        print " updated based on  ab initio prediction"
                    else :
                        n += 1
                        #if t[4]!="ungrouped" :
                        #    gffline = "\n".join(gff[t[3]])
                        #    fo.write(gffline)
                        #    fo.write("\n")
                        print "  no Gene\t",t[3]
                elif cc[i] == "LowConf" :
                    for k in gff[t[3]] :
                        k = k + ";Note=LowConfidence"
                        fo.write(k)
                        fo.write("\n")
                        print " updated, LowConf"
            upStats[1][1] += isFind
            if isFind == len(genes) :
                print "Succeed...Mis-merge was solved, gene model was updated to LowConf as the gene was ungrouped\t",t[3]

            elif isFind >0 :
                if t[4]!="ungrouped" :
                    gffline = "\n".join(gff[t[3]])
                    fo.write(gffline)
                    fo.write("\n")
                    print "Fail...Mis-merge gene was partly solved, gene model was not changed as the gene was grouped\t",t[3]
                else :
                    for k in gff[t[3]] :
                        k = k + ";Note=LowConfidence"
                        fo.write(k)
                        fo.write("\n")
                    print "Fail...Mis-merge gene was partly solved, gene model was updated to LowConf as the gene was ungrouped\t",t[3]
            else :
                if t[4]!="ungrouped" :
                    gffline = "\n".join(gff[t[3]])
                    fo.write(gffline)
                    fo.write("\n")
                    print "Fail...Mis-merge not solved, gene model was not changed as the gene was grouped\t",t[3]
                else :
                    for k in gff[t[3]] :
                        k = k + ";Note=LowConfidence"
                        fo.write(k)
                        fo.write("\n")
                    print "Fail...Mis-merge not solved, gene model was updated to LowConf as the gene was ungrouped\t",t[3]

        elif t[-1] == "unchange" :
            gffline = "\n".join(gff[t[3]])
            fo.write(gffline)
            fo.write("\n")
            #print " unchange" , t[3]
            upStats[0][0] += 1
            upStats[1][0] += 1
        elif t[-1] == "changeSci" or t[-1] == "changeSciAug":
            refId = t[7]
            #if re.search("\.",t[7]) : refId = t[7].split(".")[0]
            if checkGenes.has_key(refId) : continue
            checkGenes[refId] = 1
            isFind = 0
            flag = checkScipio(scigene[refId],refSeq)
            print line.strip()
            if t[6] == "mis-exon" :
                upStats[0][3] += 1
            elif t[6] == "add" :
            	upStats[0][4] += 1
            if flag == "good" :
                lines = getSciGene(scigene[refId])
                if lines != "" : fo.write("{}\n".format(lines))
                print " updated based on Scipio"
                isFind = 1
            else :
                if (not augGenes.has_key(t[0])) or (not snapGenes.has_key(t[0])) or ( not glimGenes.has_key(t[0])) :
                    lines = ""
                    print "  no gene was predicted in the contig :", t[0]
                else :
                    isAdd=0
                    if t[6] == "add" : isAdd = 1
                    strand = scigene[refId][0][6]
                    lines = checkAugAbiGene(int(t[1]),int(t[2]),augGenes[t[0]],augLines[t[0]],snapGenes[t[0]],snapLines[t[0]],glimGenes[t[0]],glimLines[t[0]],isAdd,strand)

                if lines != "" :
                    fo.write("{}\n".format(lines))
                    print " updated based on ab initio prediction"
                    isFind = 1
                else :
                    protFile = tempdir + "/" + refId + ".prot.fa"

                    foP = open(protFile,"w")
                    protFa = refProt[refId]

                    foP.write(">" + refId + "\n")
                    foP.write(str(protFa) + "\n")
                    #sys.exit()
                    foP.close()

                    faFile = tempdir + "/" + refId + ".genome.fa"
                    #if refId == "AT1G04105" :
                    #    print flag,"\n"
                    #    sys.exit()
                    [tch,start,end,strand]=flag.split("\t")
                    if abs(int(end)-int(start)) > 100000 :
                        if blastn.has_key(refId) :
                            [blastnCh,blastnStart,blastnEnd] = blastn[refId]
                            end = int(blastnEnd)
                            start = int(blastnStart)
                    foF = open(faFile,"w")
                    subseq = refSeq[tch][int(start)- 200 -1 : int(end) + 200]
                    foF.write(">" + tch + "_" + str(start) + "_" + str(end) + "\n")
                    foF.write(str(subseq) + "\n")
                    foF.close()

                    if tch != t[0] :
                        print "tch != t[0]t", line,"|n"
                        lines = ""
                    else :
                        lines = getGenewise(faFile, protFile, tempdir, augCfg, strand, start)

                    if lines != "" :
                        fo.write("{}\n".format(lines))
                        print " updated based on Genewise"
                        isFind = 1
                    else :
                        if t[6]== "mis-exon" or (t[4] != "ungrouped" and t[4]!= "**"):
                            gffline = "\n".join(gff[t[3]])
                            fo.write(gffline)
                            fo.write("\n")
                            print "Fail...  correct mis-exon " , t[3], "gene model was kept"
                        else :
                            print "Fail update gene",line
            if isFind == 1 :
                if t[6]  == "mis-exon":
                    print "Succeed update mis-exon gene", t[3]
                    upStats[1][3]+=1
                elif t[6] == "add" :
                    print "Succeed add gene", t[3]
                    upStats[1][4]+=1

        elif t[-1] == "changeAug" :
            isFind = 0
            if t[6] == "mis-exon" :
                upStats[0][3] += 1
            elif t[6] == "add" :
            	upStats[0][4] += 1
            refId = t[7]
            print line.strip()
            #if re.search("\.",t[7]) : refId = t[7].split(".")[0]
            if checkGenes.has_key(refId) : continue
            checkGenes[refId] = 1
            if (not augGenes.has_key(t[0])) or (not snapGenes.has_key(t[0])) or ( not glimGenes.has_key(t[0])) :
                lines = ""
                print "  no gene was predicted in the contig :", t[0]
                isFind = 1
            else :
                #lines = checkAugAbiGene(int(t[1]),int(t[2]),augGenes[t[0]],augLines[t[0]],snapGenes[t[0]],snapLines[t[0]],glimGenes[t[0]],glimLines[t[0]] )
                isAdd=0
                if t[6] == "add" : isAdd = 1
                strand = "both"
                lines = checkAugAbiGene(int(t[1]),int(t[2]),augGenes[t[0]],augLines[t[0]],snapGenes[t[0]],snapLines[t[0]],glimGenes[t[0]],glimLines[t[0]],isAdd,strand )

            if lines != "" :
                fo.write("{}\n".format(lines))
                print " updated based on  ab initio prediction"
                isFind = 1
            else :
                protFile = tempdir + "/" + refId + ".prot.fa"
                foP = open(protFile,"w")
                protFa = refProt[refId]
                foP.write(">" + refId + "\n")
                foP.write(str(protFa) + "\n")
                foP.close()

                faFile = tempdir + "/" + refId + ".genome.fa"
                #[start,end,strand]=flag.split("\t")
                start = int(t[1])
                end = int(t[2])
                strand = "both"
                foF = open(faFile,"w")
                subseq = refSeq[t[0]][int(start)- 200 -1 : int(end) + 200]
                foF.write(">" + t[0] + "_" + str(start) + "_" + str(end) + "\n")
                foF.write(str(subseq) + "\n")
                foF.close()
                lines = getGenewise(faFile, protFile, tempdir, augCfg, strand, start)
                if lines != "" :
                    fo.write("{}\n".format(lines))
                    print " updated based on Genewise"
                    isFind = 1
                else :
                    if t[6]== "mis-exon" or (t[4] != "ungrouped" and t[4]!= "**"):
                        gffline = "\n".join(gff[t[3]])
                        fo.write(gffline)
                        fo.write("\n")
                        print "Fail correct mis-exon\t",t[3], "gene model was kept"
                    else :
                        print "Fail update gene",line
                #if t[6]== "mis-exon":
                #    gffline = "\n".join(gff[t[3]])
                #    fo.write(gffline)
                #    fo.write("\n")
                #print "NoGene\t",line
            if isFind == 1 :
                if t[6]  == "mis-exon":
                    print "Succeed update mis-exon gene", t[3]
                    upStats[1][3] += 1
                elif t[6] == "add" :
                    print "Succeed add gene", t[3]
                    upStats[1][4] += 1
        elif t[-1] == "LowConf" :
            for k in gff[t[3]] :
                k = k + ";Note=LowConfidence"
                fo.write(k)
                fo.write("\n")
            print " updated as LowConf"

        elif t[-1] == "not-add" :
            print " updated : not add"
            #upStats[6] += 1
            continue

    fi.close()
    fo.close()
    print "upstats: " ,  upStats
    return upStats

def checkUpdate(sciGene,updateFile,outfile):
    fi = open(updateFile,"r")
    fo = open(outfile,"w")
    genesUp = ""
    misSplitDup = {}
    toUp = [0,0,0,0,0,0,0,0,0,0,0,0]
    toUp2 = ["unchange", "mis-split","mis-exon", "mis-merQry", "mis-merRef","add","mis-ungr","Ref-diff-ann","Ref-Not-ass","Ref-Par-ass","Ref-No-ann","RefNoAnn-RNAsup"]
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        sci = "sci"
        flag = "**"
        if  t[6] == "unchange" or t[6] == "mis-split" :
            sci = "**"
            flag = "unchange"
            genesUp = "unchange"
            if (t[6] == "unchange") :
                toUp[0]+=1
            else :
                toUp[1]+=1
        elif t[6] == "mis-merge" :
            refGenes = t[7].split(",")
            sci = ""
            flag = ""
            genesUp = ""
            toUp[3]+=1
            toUp[4]+=len(refGenes)
            for g in refGenes :
                fsci = ""
                if sciGene.has_key(g) :
                    if sci == "" :
                        sci = "sci"
                    else :
                        sci = sci + "," + "sci"
                    fsci = "sci"
                else :
                    if sci == "" :
                        sci = "None"
                    else :
                        sci = sci + "," + "None"
                    fsci = "None"
                if flag == "" :
                    if  fsci == "sci" :
                        flag = "misMe-Change-sci"
                        genesUp ="changeSci"
                        if sciGene[g][0][0] != t[0] :
                            genesUp = "changeAug"
                    else :
                        flag = "misMe-NoChange-NOsci"
                        genesUp = "changeAug"
                else :
                    if  fsci == "sci" :
                        flag = flag + "," + "misMe-Change-sci"
                        if sciGene[g][0][0] != t[0] :
                            genesUp = genesUp + "," + "changeAug"
                        else :
                            genesUp = genesUp + "," + "changeSci"
                    else :
                        flag = flag + "," +"misMe-NoChange-NOsci"
                        genesUp = genesUp + "," + "changeAug"

        elif t[6] == "mis-exon" :
            #refGene = t[7].split(".")[0]
            refGene = t[7]
            toUp[2]+=1
            if not sciGene.has_key(refGene) :
                sci = "None"
            if t[4] != "ungrouped"  :
                genesUp = "unchange"
            elif sci == "sci" :
                flag = "misEx-Change-sci"
                genesUp = "changeSci"
                if sciGene[refGene][0][0] != t[0] :
                    genesUp = "changeAug"
            else :
                flag = "misEx-NoChange-NOsci"
                genesUp = "changeAug"

        elif t[6] == "mis-ungr" :
            sci = "**"
            flag = "mis-ungr"
            genesUp = "unchange"
            toUp[6] += 1
        elif t[6] == "RefNo-RNAsup" :
            sci = "**"
            flag = "RNA-sup"
            genesUp = "unchange"
            toUp[11] +=1
        elif t[6] == "Ref-No-ann" :
            flag = "Ref-No-ann"
            genesUp = "unchange"
            toUp[10] += 1
        elif t[6] == "Ref-Not-ass" :
            LoF = "**"
            sci = "**"
            flag = "RefNoAss)"
            genesUp = "unchange"
            toUp[8]+=1
        elif t[6] == "Ref-Par-ass" :
            LoF = "**"
            sci = "**"
            flag = "RefParAss"
            genesUp = "unchange"
            if t[4] != "ungrouped" :
                genesUp = "unchange"
            else :
                genesUp = "unchange"
            toUp[9] += 1
        elif t[6] == "add" :
            toUp[5] += 1
            if not sciGene.has_key(t[7]) :
                sci = "None"
            else :
                sci = "sci"
            if sci == "sci" :
                flag = "add-Sci"
                genesUp = "changeSci" ## further check the features of gene(start,stop codon; splicing sites)
                if sciGene[t[7]][0][0] != t[0] :
                    genesUp = "changeAug"
            else:
                genesUp = "changeAug"
        else :
            genesUp = "unchange"
            toUp[7] += 1
        LoF="NaN-LoF"
        fo.write("{}\t{}\t{}\t{}\t{}\n".format(line.strip(),LoF,sci,flag,genesUp))
    fi.close()
    fo.close()
    for i in range(len(toUp)) :
        print toUp2[i], "\t" , toUp[i]
    return

def getScipio(result):
    fi = open(result,"r")
    genes={}
    gene = ""
    flag = 1
    while True :
        line = fi.readline()
        if not line : break
        if re.search("#query",line) :
            t = line.strip().split("\t")
            gene = t[1].split(" ")[0]
            if genes.has_key(gene) :
                flag = 0
            else :
                flag = 1
                genes[gene] = []
        else :
            if flag == 0 : continue
            if line[0]=="#" : continue
            t = line.strip().split("\t")
            genes[gene].append(t)
    fi.close()
    print("loaded the scipio infor...")
    # sys.exit()        # Change by MG
    return genes

def getSciGene (scig):
    scilines = scig
    strand = scilines[0][6]
    outlines = ""
    start = ""
    end = ""
    id  = scilines[0][8].split(";")[1]
    id = id.split()[0]
    id = id.split("=")[1]
    id = "Scipio-" + id
    infor = "ID=" + id
    if len(scilines) > 1 :
        for k in range(len(scilines)) :
            t = scilines[k]
            ll = t[:]
            if k == 0 :
                if strand == "+" :
                    start = t[3]
                else :
                    end = t[4]

            elif k == len(scilines) - 1:
                if strand == "+" :
                    t[4] = int(t[4]) + 3
                    t[4] = str(t[4])
                    end = t[4]
                    ll[4] = t[4]
                else :
                    t[3] = int(t[3]) - 3
                    t[3] = str(t[3])
                    start = t[3]
                    ll[3] = t[3]

            ll[5] = "."
            ll[2] = "CDS"
            ll[8] = "ID=" + id + ".t1.cds" + str(k + 1) + ";Parent=" + id + ".t1"
            cdsLine = "\t".join(ll)
            ll[2] = "exon"
            ll[7] = "."
            ll[8] = "ID=" + id + ".t1.exon" + str(k + 1) + ";Parent=" + id + ".t1"
            exonLine = "\t".join(ll)
            if outlines == "" :
                outlines = exonLine + "\n" + cdsLine
            else :
                outlines = outlines + "\n" + exonLine + "\n" + cdsLine
    else :
        t = scilines[0]
        ll = t[:]
        if strand == "+" :
            ll[4] = str(int(t[4]) + 3)
            start = scilines[0][3]
            end = str(int(scilines[0][4]) + 3)
        else :
            ll[3] = str(int(t[3]) - 3)
            start = str(int(scilines[0][3]) - 3)
            end = scilines[0][4]

        ll[5] = "."
        ll[2] = "CDS"
        ll[8] = "ID=" + id + ".t1.cds" + str(1) + ";Parent=" + id + ".t1"
        cdsLine = "\t".join(ll)
        ll[2] = "exon"
        ll[7] = "."
        ll[8] = "ID=" + id + ".t1.exon" + str(1) + ";Parent=" + id + ".t1"
        exonLine = "\t".join(ll)
        if outlines == "" :
            outlines = exonLine + "\n" + cdsLine
        else :
            outlines = outlines + "\n" + exonLine + "\n" + cdsLine

    geneLine = t[0] + "\t" + "Scipio\t" + "gene" + "\t" + start + "\t" + end + "\t" + ".\t" + strand + "\t.\t" + infor
    mRNALine = t[0] + "\t" + "Scipio\t" + "mRNA" + "\t" + start + "\t" + end + "\t" + ".\t" + strand + "\t.\t" + infor + ".t1;Parent=" + id
    lines = geneLine + "\n" + mRNALine + "\n" + outlines
    return lines

def getGeneWiseGene (gwlines,name, shiftVal):
    scilines = gwlines
    strand = scilines[0][6]
    outlines = ""
    start = ""
    end = ""

    id = "GW-" + name
    infor = "ID=" + id
    if len(scilines) > 1 :
        for k in range(len(scilines)) :
            t = scilines[k]
            ll = t[:]

            if k == 0 :
                if strand == "+" :
                    start = str(int(shiftVal) + int(t[3]) - 200 -1 )
                else :
                    end = str(int(shiftVal) + int(t[4]) - 200 -1 )

            elif k == len(scilines) - 1:
                if strand == "+" :
                    t[4] = int(t[4]) + 3
                    t[4] = str(t[4])
                    end = str(int(shiftVal) + int(t[4]) - 200 -1 )
                else :
                    t[3] = int(t[3]) - 3
                    t[3] = str(t[3])
                    #start = t[3]
                    start = str(int(shiftVal) + int(t[3]) - 200 -1 )
            ll[3] = str(int(shiftVal) + int(t[3]) - 200 -1 )
            ll[4] = str(int(shiftVal) + int(t[4]) - 200 -1 )

            # ll[0] = ll[0].split("_")[0] # Change MG
            ll[0] = '_'.join(ll[0].split("_")[:-2])
            ll[5] = "."
            ll[2] = "CDS"
            ll[8] = "ID=GW-" + name + ".t1.cds" + str(k + 1) + ";Parent=" + id + ".t1"
            cdsLine = "\t".join(ll)

            ll[2] = "exon"
            ll[7] = "."
            ll[8] = "ID=GW-" + name + ".t1.exon" + str(k + 1) + ";Parent=" + id + ".t1"
            exonLine = "\t".join(ll)
            if outlines == "" :
                outlines = exonLine + "\n" + cdsLine
            else :
                outlines = outlines + "\n" + exonLine + "\n" + cdsLine
    else :
        t = scilines[0]
        ll = t[:]
        if strand == "+" :
            start = str(int(shiftVal) + int(t[3]) - 200 -1 )
            end = str(int(shiftVal) + int(t[4]) - 200 -1 + 3)
        else :
            start = str(int(shiftVal) + int(t[3]) - 200 -1 -3)
            end = str(int(shiftVal) + int(t[4]) - 200 -1 )

        ll[3] = start
        ll[4] = end
        # ll[0] = ll[0].split("_")[0] # Change by MG
        ll[0] = '_'.join(ll[0].split("_")[:-2])
        ll[5] = "."
        ll[2] = "CDS"
        ll[8] = "ID=GW-" + name + ".t1.cds" + str(1) + ";Parent=" + id + ".t1"
        cdsLine = "\t".join(ll)

        ll[2] = "exon"
        ll[7] = "."
        ll[8] = "ID=GW-" + name + ".t1.exon" + str(1) + ";Parent=" + id + ".t1"
        exonLine = "\t".join(ll)

        outlines = exonLine + "\n" + cdsLine

    # geneLine = scilines[0][0].split("_")[0] + "\t" + "GeneWise\t" + "gene" + "\t" + start + "\t" + end + "\t" + ".\t" + strand + "\t.\t" + infor # Change by MG
    geneLine = '_'.join(scilines[0][0].split("_")[:-2]) + "\t" + "GeneWise\t" + "gene" + "\t" + start + "\t" + end + "\t" + ".\t" + strand + "\t.\t" + infor
    # mRNALine = scilines[0][0].split("_")[0] + "\t" + "GeneWise\t" + "mRNA" + "\t" + start + "\t" + end + "\t" + ".\t" + strand + "\t.\t" + infor + ".t1;Parent=GW-" + name # Change by MG
    mRNALine = '_'.join(scilines[0][0].split("_")[:-2]) + "\t" + "GeneWise\t" + "mRNA" + "\t" + start + "\t" + end + "\t" + ".\t" + strand + "\t.\t" + infor + ".t1;Parent=GW-" + name
    lines = geneLine + "\n" + mRNALine + "\n" + outlines
    return lines

def checkScipio(genes,seq):
    ch = genes[0][0]
    strand = genes[0][6]
    flag = "xx;"
    startP = 0
    endP = 0
    if strand == "+" :
        start = int(genes[0][3])
        end = int(genes[len(genes)-1][4])
        startP = start
        endP = end
        startCodon = seq[ch][start -1 :start+2]
        stopCodon = seq[ch][end:end+3]
        if startCodon != "ATG" :
            flag = flag + ";start-lost"
        if not (stopCodon == "TAG" or stopCodon == "TGA" or stopCodon == "TAA" ) :
            flag = ";stop-lost"
        cds = seq[ch][int(genes[0][3])-1 : int(genes[0][4])]
        if len(genes) > 1 :
            cds = ""
            for i in range(len(genes)) :
                if cds == "" :
                    cds = seq[ch][int(genes[i][3])-1 : int(genes[i][4])]
                else :
                    cds = cds + seq[ch][int(genes[i][3])-1 : int(genes[i][4])]
                if i == len(genes) - 1 :
                    break
                else :
                    splice5 = seq[ch][int(genes[i][4]) : int(genes[i][4]) + 2]
                    if splice5 != "GT" :
                        flag = ";splice5-lost"
                    splice3 = seq[ch][int(genes[i+1][3]) - 3 : int(genes[i+1][3]) - 1]
                    if splice3 != "AG" :
                        flag = ";splice3-lost"
        for i in range(len(cds)/3) :
            codon = cds[i*3:i*3+3]
            if codon == "TAG" or codon == "TGA" or codon == "TAA" :
                flag = "stop_gained"
                pos = str(i*3)
                print "stop_gained_scipio at", pos, genes[0][8]

    else :
        start = int(genes[0][4])
        end = int(genes[len(genes)-1][3])
        startP = end
        endP = start
        startCodon = seq[ch][start - 3 : start]
        stopCodon = seq[ch][end-4 : end -1]
        if startCodon != "CAT" :
            flag = flag + ";start-lost"
        if not (stopCodon == "CTA" or stopCodon == "TCA" or stopCodon == "TTA" ) :
            flag = ";stop-lost"
        cds = seq[ch][int(genes[0][3])-1 : int(genes[0][4])]

        if len(genes) > 1 :
            cds = ""
            for i in range(len(genes)) :
                if cds == "" :
                    cds = seq[ch][int(genes[i][3])-1 : int(genes[i][4])]
                else :
                    cds = seq[ch][int(genes[i][3])-1 : int(genes[i][4])] + cds

                if i == len(genes) - 1 :
                    break
                else :
                    splice5 = seq[ch][int(genes[i][3]) -3 : int(genes[i][3]) - 1]
                    if splice5 != "AC" :
                        flag = ";splice5-lost"
                    splice3 = seq[ch][int(genes[i+1][4]) : int(genes[i+1][4]) + 2]
                    if splice3 != "CT" :
                        flag = ";splice3-lost"
        for i in range(len(cds)/3) :
            codon = cds[i*3:i*3+3]
            if codon == "CTA" or codon == "TCA" or codon == "TTA" :
                flag = "stop_gained"
                pos = str(i*3)
                print "stop_gained_scipio at", pos, genes[0][8]
    cdsLen = 0
    for k in range(len(genes)) :
        cdsLen += int(genes[k][4]) - int(genes[k][3]) + 1
    if cdsLen % 3 != 0 :
        flag = flag + ";frame-shift"

    if flag == "xx;" :
        flag = "good"
        return flag
    else :
        print genes[0][8], flag
        flag = ch + "\t" + str(startP) + "\t" + str(endP) + "\t" + strand
        return  flag

def checkGeneWise(genes,seq,isATG):
    ch = genes[0][0]
    strand = genes[0][6]
    flag = "xx;"
    startP = 0
    endP = 0
    stopPos = []
    stopCds = []
    if strand == "+" :
        start = int(genes[0][3])
        end = int(genes[len(genes)-1][4])
        startP = start
        endP = end
        startCodon = seq[ch][start -1 :start+2]
        stopCodon = seq[ch][end:end+3]
        if startCodon != "ATG" :
            if isATG == 1 :
                flag = flag + ";start-lost"
        if not (stopCodon == "TAG" or stopCodon == "TGA" or stopCodon == "TAA" ) :
            flag = flag + ";stop-lost"
        cds = seq[ch][int(genes[0][3])-1 : int(genes[0][4])]

        if len(genes) > 1 :
            cds = ""
            for i in range(len(genes)) :
                if cds == "" :
                    cds = seq[ch][int(genes[i][3])-1 : int(genes[i][4])]
                else :
                    cds = cds + seq[ch][int(genes[i][3])-1 : int(genes[i][4])]

                if i == len(genes) - 1 :
                    break
                else :
                    splice5 = seq[ch][int(genes[i][4]) : int(genes[i][4]) + 2]
                    if splice5 != "GT" :
                        flag = ";splice5-lost"
                    splice3 = seq[ch][int(genes[i+1][3]) - 3 : int(genes[i+1][3]) - 1]
                    if splice3 != "AG" :
                        flag = ";splice3-lost"

        for i in range(len(cds)/3) :
            codon = cds[i*3:i*3+3]
            if codon == "TAG" or codon == "TGA" or codon == "TAA" :
                flag = flag + ";stop_gained"
                pos = str(i*3)
                print "stop_gained_genewise at forward strand", pos, genes[0][8]
                stopPos.append(pos)
                #stopCds.append(i)

    else :
        start = int(genes[0][4])
        end = int(genes[len(genes)-1][3])
        startP = end
        endP = start
        startCodon = seq[ch][start - 3 : start]
        stopCodon = seq[ch][end-4 : end -1]
        if startCodon != "CAT" :
            if isATG == 1 :
                flag = flag + ";start-lost"
        if not (stopCodon == "CTA" or stopCodon == "TCA" or stopCodon == "TTA" ) :
            flag = flag + ";stop-lost"

        cds = seq[ch][int(genes[0][3])-1 : int(genes[0][4])]

        if len(genes) > 1 :
            cds = ""
            for i in range(len(genes)) :
                if cds == "" :
                    cds = seq[ch][int(genes[i][3])-1 : int(genes[i][4])]
                else :
                    cds = seq[ch][int(genes[i][3])-1 : int(genes[i][4])] + cds

                if i == len(genes) - 1 :
                    break
                else :
                    splice5 = seq[ch][int(genes[i][3]) -3 : int(genes[i][3]) - 1]
                    if splice5 != "AC" :
                        flag = ";splice5-lost"
                    splice3 = seq[ch][int(genes[i+1][4]) : int(genes[i+1][4]) + 2]
                    if splice3 != "CT" :
                        flag = ";splice3-lost"

        for i in range(len(cds)/3) :
            codon = cds[i*3:i*3+3]
            if codon == "CTA" or codon == "TCA" or codon == "TTA" :
                flag = flag + ";stop_gained"
                pos = str(i*3)
                print "stop_gained_genewise at reversed strand", pos, genes[0][8]
                stopPos.append(pos)

    if flag == "xx;" :
        flag = "good"
        return flag
    elif flag == "xx;;stop_gained" :
        stop_1 = stopPos[0]
        per = float(stop_1)/float(len(cds))
        if strand == "-" : per = 1.00 - float(stop_1)/float(len(cds))
        print "stop_gained position:",genes[0][8].split(";")[1], stop_1,float(stop_1)/float(len(cds))
        if  per > 0.70:
            print "although stop_gained, stop_1 pos still generate a protein with 70% of orginal protein",genes[0][8].split(";")[1]
            flag = ch + "\t" + str(startP) + "\t" + str(endP) + "\t" + strand + "\t" + stopPos[0] + "\tstopEnd"
            #sys.exit()
            return flag
        else :
            flag = ch + "\t" + str(startP) + "\t" + str(endP) + "\t" + strand + "\t" + flag
            return flag
    else :
        print genes[0][8], flag
        flag = ch + "\t" + str(startP) + "\t" + str(endP) + "\t" + strand + "\t" + flag
        return  flag

def getAugEviGene(gff):

    fi = open(gff,"r")
    genes2 = {}
    genes1 = {}
    strand = ""
    start = ""
    end = ""
    id = ""
    genes = []
    ch = ""
    while True :
        line = fi.readline()
        if not line : break
        if line[0]=="#" : continue
        t = line.strip().split("\t")
        if len(t) < 9 : continue
        if not genes2.has_key(t[0]) : genes2[t[0]] = {}
        if not genes1.has_key(t[0]) : genes1[t[0]] = []
        if t[2] == "gene" :
            if ch =="" :
                ch = t[0]
                id = t[0] + "-" + t[8].split("=")[1]
                genes1[t[0]].append([t[3],t[4],id,t[6]])

                genes2[t[0]][id] = []
                t[8] = "ID=" + id
                genes2[t[0]][id].append(t)
                mRNA = t[:]
                mRNA[2] = "mRNA"
                mRNA[8] = "ID=" + id + ".t1;Parent=" + id
                genes2[t[0]][id].append(mRNA)
                strand = t[6]

            elif ch != t[0] :
                genes1[ch][-1][0] = start
                genes1[ch][-1][1] = end

                genes2[ch][id][0][3] = start
                genes2[ch][id][0][4] = end

                genes2[ch][id][1][3] = start
                genes2[ch][id][1][4] = end

                genes2[ch][id][0] = "\t".join(genes2[ch][id][0])
                genes2[ch][id][1] = "\t".join(genes2[ch][id][1])

                ch = t[0]
                id = t[0] + "-" + t[8].split("=")[1]
                genes1[t[0]].append([t[3],t[4],id,t[6]])

                genes2[t[0]][id] = []
                t[8] = "ID=" + id
                genes2[t[0]][id].append(t)
                mRNA = t[:]
                mRNA[2] = "mRNA"
                mRNA[8] = "ID=" + id + ".t1;Parent=" + id
                genes2[t[0]][id].append(mRNA)
                strand = t[6]

            else :
                #print genes1[t[0]][-1]
                #sys.exit()
                genes1[t[0]][-1][0] = start
                genes1[t[0]][-1][1] = end


                genes2[t[0]][id][0][3] = start
                genes2[t[0]][id][0][4] = end

                genes2[t[0]][id][1][3] = start
                genes2[t[0]][id][1][4] = end

                genes2[t[0]][id][0] = "\t".join(genes2[t[0]][id][0])
                genes2[t[0]][id][1] = "\t".join(genes2[t[0]][id][1])

                id = t[0] + "-" + t[8].split("=")[1]
                genes1[t[0]].append([t[3],t[4],id,t[6]])

                genes2[t[0]][id] = []
                t[8] = "ID=" + id

                #id = t[0] + "."+ t[8].split("=")[1]
                #genes1[t[0]].append([t[3],t[4],id])
                #genes2[t[0]][id] = []
                genes2[t[0]][id].append(t)
                mRNA = t[:]
                mRNA[2] = "mRNA"
                mRNA[8] = "ID=" + id + ".t1;Parent=" + id
                genes2[t[0]][id].append(mRNA)
                strand = t[6]
            start = ""
        elif t[2] == "CDS" :
            #if strand == "+" :
            if not re.search(".t1.",t[8]) : continue
            if start == "" : start = t[3]
            end = t[4]
            exon = t[:]
            exon[2]="exon"
            exon[5]="."
            exon[7]="."
            exon[8] = "ID=" + id + ".t1.exon;Parent=" + id + ".t1"
            exonLine = "\t".join(exon)

            t[8] = "ID=" + id + ".t1.cds;Parent=" + id + ".t1"
            cdsLine = "\t".join(t)
            genes2[t[0]][id].append(exonLine)
            genes2[t[0]][id].append(cdsLine)
    if ch != "" :
        genes2[ch][id][0][3] = start
        genes2[ch][id][0][4] = end
        genes2[ch][id][1][3] = start
        genes2[ch][id][1][4] = end
        genes2[ch][id][0] = "\t".join(genes2[ch][id][0])
        genes2[ch][id][1] = "\t".join(genes2[ch][id][1])

    fi.close()
    #print genes1[ch][0],"\n"

    #for i in genes2["chr1"]["g2"] :
    #    print i
    #print "\n"
    #sys.exit()
    return [genes1,genes2]

def getAbiGene(file):
    fi = open(file,"r")
    genes2 = {}
    genes1 = {}
    strand = ""
    start = ""
    end = ""
    id = ""
    ch = ""
    genes = []
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        if len(t) < 9 : continue
        if not genes2.has_key(t[0]) : genes2[t[0]] = {}
        if not genes1.has_key(t[0]) : genes1[t[0]] = []
        if t[2] == "gene" :
            if ch=="" :
                ch = t[0]
                id = t[8].split(";")[0]
                id = id.split("=")[1]
                genes1[t[0]].append([t[3],t[4],id,t[6]])

                genes2[t[0]][id] = []
                genes2[t[0]][id].append(t)

                strand = t[6]
            elif ch != t[0] :
                genes1[ch][-1][0] = start
                genes1[ch][-1][1] = end

                genes2[ch][id][0][3] = start
                genes2[ch][id][0][4] = end

                genes2[ch][id][1][3] = start
                genes2[ch][id][1][4] = end

                genes2[ch][id][0] = "\t".join(genes2[ch][id][0])
                genes2[ch][id][1] = "\t".join(genes2[ch][id][1])

                ch = t[0]
                id = t[8].split(";")[0]
                id = id.split("=")[1]
                genes1[t[0]].append([t[3],t[4],id,t[6]])

                genes2[t[0]][id] = []
                genes2[t[0]][id].append(t)

            else :
                genes1[t[0]][-1][0] = start
                genes1[t[0]][-1][1] = end

                genes2[t[0]][id][0][3] = start
                genes2[t[0]][id][0][4] = end

                genes2[t[0]][id][1][3] = start
                genes2[t[0]][id][1][4] = end

                genes2[t[0]][id][0] = "\t".join(genes2[t[0]][id][0])
                genes2[t[0]][id][1] = "\t".join(genes2[t[0]][id][1])

                id = t[8].split(";")[0]
                id = id.split("=")[1]
                genes1[t[0]].append([t[3],t[4],id,t[6]])
                genes2[t[0]][id] = []
                genes2[t[0]][id].append(t)

                strand = t[6]
            #if strand == "+" :
            start = ""
            end = ""
        elif t[2] == "mRNA" :
            genes2[t[0]][id].append(t)
        elif t[2] == "exon" :
            genes2[t[0]][id].append("\t".join(t))
        elif t[2] == "CDS" :
            if strand == "+" :
                if start == "" : start = t[3]
                end = t[4]
            else :
                if end == "" : end = t[4]
                start = t[3]
            #genes2[t[0]][id].append(exon)
            genes2[t[0]][id].append("\t".join(t))

    fi.close()
    genes2[ch][id][0][3] = start
    genes2[ch][id][0][4] = end
    genes2[ch][id][1][3] = start
    genes2[ch][id][1][4] = end
    genes2[ch][id][0] = "\t".join(genes2[ch][id][0])
    genes2[ch][id][1] = "\t".join(genes2[ch][id][1])

    if genes2.has_key("chr1") and genes2["chr1"].has_key("gene.chr1-snap.2") :
        for i in genes2["chr1"]["gene.chr1-snap.2"] :
            print i
        print "\n"
    ##ignore the last one gene

    return [genes1,genes2]

def checkAugAbiGene(start,end,augGeneCh,augGeneLines,snapGeneCh,snapGeneLines,glimGeneCh,glimGeneLines,isAdd,strand):
    genesCh = augGeneCh
    geneLines = ""
    for i in range(len(genesCh)) :
        if( start >= int(genesCh[i][0]) and start <= int(genesCh[i][1])) or ( end >= int(genesCh[i][0]) and end <= int(genesCh[i][1]) ) or ( start <= int(genesCh[i][0]) and end >= int(genesCh[i][1]) ) :
            srtArr = sorted([start,end,int(genesCh[i][0]), int(genesCh[i][1])],key=int)
            olap = srtArr[2] - srtArr[1]
            if olap >= 0.8*(end-start) and olap >= 0.8 * (int(genesCh[i][1]) - int(genesCh[i][0]) ) :
                if strand != "both" and strand != genesCh[i][3] : continue
                geneId = genesCh[i][2]
                #if len(augGeneLines[geneId])==0 :
                #    print geneId,"\n"
                #    sys.exit()
                #print augGeneLines[geneId]
                geneLines ="\n".join(augGeneLines[geneId])
                break
            elif isAdd == 1 and olap >= 0.6*(end-start) and olap >= 0.8 * (int(genesCh[i][1]) - int(genesCh[i][0]) ):
                if strand != "both" and strand != genesCh[i][3] : continue
                geneId = genesCh[i][2]
                #if len(augGeneLines[geneId])==0 :
                #    print geneId,"\n"
                #    sys.exit()
                #print augGeneLines[geneId]
                geneLines ="\n".join(augGeneLines[geneId])
                print "isadd , 0.6 overlap" , start, end, geneId
                break

    if geneLines == "" :
        genesCh = snapGeneCh
    #genesChr = abiGenes2[t[0]]
        for i in range(len(genesCh)) :
            if( start >= int(genesCh[i][0]) and start <= int(genesCh[i][1])) or ( end >= int(genesCh[i][0]) and end <= int(genesCh[i][1]) ) or ( start <= int(genesCh[i][0]) and end >= int(genesCh[i][1]) ) :
                srtArr = sorted([start,end,int(genesCh[i][0]), int(genesCh[i][1])],key=int)
                olap = srtArr[2] - srtArr[1]
                if olap >= 0.8*(end-start) and olap >= 0.8 * (int(genesCh[i][1]) - int(genesCh[i][0]) ):
                    if strand != "both" and strand != genesCh[i][3] : continue
                    geneId = genesCh[i][2]
                    geneLines ="\n".join(snapGeneLines[geneId])
                    break
                elif isAdd== 1 and olap >= 0.6*(end-start) and olap >= 0.8 * (int(genesCh[i][1]) - int(genesCh[i][0]) ):
                    if strand != "both" and strand != genesCh[i][3] : continue
                    geneId = genesCh[i][2]
                    geneLines ="\n".join(snapGeneLines[geneId])
                    print "isadd , 0.6 overlap" , start, end, geneId
                    break
        if geneLines == "" :
            genesCh = glimGeneCh
            for i in range(len(genesCh)) :
                if( start >= int(genesCh[i][0]) and start <= int(genesCh[i][1])) or ( end >= int(genesCh[i][0]) and end <= int(genesCh[i][1]) ) or ( start <= int(genesCh[i][0]) and end >= int(genesCh[i][1]) ) :
                    srtArr = sorted([start,end,int(genesCh[i][0]), int(genesCh[i][1])],key=int)
                    olap = srtArr[2] - srtArr[1]
                    if olap >= 0.8*(end-start) and olap >= 0.8 * (int(genesCh[i][1]) - int(genesCh[i][0]) ):
                        if strand != "both" and strand != genesCh[i][3] : continue
                        geneId = genesCh[i][2]
                        geneLines ="\n".join(glimGeneLines[geneId])
                        break
                    elif isAdd==1 and  olap >= 0.8*(end-start) and olap >= 0.8 * (int(genesCh[i][1]) - int(genesCh[i][0]) ):
                        if strand != "both" and strand != genesCh[i][3] : continue
                        geneId = genesCh[i][2]
                        geneLines ="\n".join(glimGeneLines[geneId])
                        print "isadd , 0.6 overlap" , start, end, geneId
                        break
    return geneLines

def checkAugAbiGene2(start,end,augGeneCh,augGeneLines,snapGeneCh,snapGeneLines,glimGeneCh,glimGeneLines,prevStart,nextStart):
    genesCh = augGeneCh
    geneLines = ""
    #lines = checkAugAbiGene2(int(blastnStart),int(blastnEnd),augGenes[t[0]],augLines[t[0]],snapGenes[t[0]],snapLines[t[0]],glimGenes[t[0]],glimLines[t[0]], prevStart,nextStart)
    for i in range(len(genesCh)) :
        if( start >= int(genesCh[i][0]) and start <= int(genesCh[i][1])) or ( end >= int(genesCh[i][0]) and end <= int(genesCh[i][1]) ) or ( start <= int(genesCh[i][0]) and end >= int(genesCh[i][1]) ) :
            srtArr = sorted([start,end,int(genesCh[i][0]), int(genesCh[i][1])],key=int)
            olap = srtArr[2] - srtArr[1]
            if olap >= 0.8*(end-start) and int(genesCh[i][1]) <= nextStart and int(genesCh[i][0]) >= prevStart:
                geneId = genesCh[i][2]
                #if len(augGeneLines[geneId])==0 :
                #    print geneId,"\n"
                #    sys.exit()
                #print augGeneLines[geneId]
                geneLines ="\n".join(augGeneLines[geneId])
                break
    if geneLines == "" :
        genesCh = snapGeneCh
    #genesChr = abiGenes2[t[0]]
        for i in range(len(genesCh)) :
            if( start >= int(genesCh[i][0]) and start <= int(genesCh[i][1])) or ( end >= int(genesCh[i][0]) and end <= int(genesCh[i][1]) ) or ( start <= int(genesCh[i][0]) and end >= int(genesCh[i][1]) ) :
                srtArr = sorted([start,end,int(genesCh[i][0]), int(genesCh[i][1])],key=int)
                olap = srtArr[2] - srtArr[1]
                if olap >= 0.8*(end-start) and  int(genesCh[i][1]) <= nextStart and int(genesCh[i][0]) >= prevStart:
                    geneId = genesCh[i][2]
                    geneLines ="\n".join(snapGeneLines[geneId])
                    break
        if geneLines == "" :
            genesCh = glimGeneCh
            for i in range(len(genesCh)) :
                if( start >= int(genesCh[i][0]) and start <= int(genesCh[i][1])) or ( end >= int(genesCh[i][0]) and end <= int(genesCh[i][1]) ) or ( start <= int(genesCh[i][0]) and end >= int(genesCh[i][1]) ) :
                     srtArr = sorted([start,end,int(genesCh[i][0]), int(genesCh[i][1])],key=int)
                     olap = srtArr[2] - srtArr[1]
                     if olap >= 0.8*(end-start) and int(genesCh[i][1]) <= nextStart and int(genesCh[i][0]) >= prevStart:
                        geneId = genesCh[i][2]
                        geneLines ="\n".join(glimGeneLines[geneId])
                        break

    return geneLines

def getGenewise(faFile,proFile,outdir,augCfg,strand,shiftVal):
    name = os.path.basename(proFile)
    outfile1 = outdir + "/" + name + ".genewise.gff"
    #if os.path.isfile(outfile1) :
    #    print outfile1, "exists, skip genewise alignment"
    if strand == "both" :
        os.system("genewise " + proFile + " " + faFile + " -both -gff -silent > " + outfile1)
    elif strand == "+" :
        os.system("genewise " + proFile + " " + faFile + " -tfor -gff -silent > " + outfile1)
    else :
        os.system("genewise " + proFile + " " + faFile + " -trev -gff -silent > " + outfile1)

    outfile2 = outdir + "/" + name + ".genewise.hints"
    fo = open(outfile2,"w")
    fi = open(outfile1,"r")
    geneLines = ""
    genes = []
    ch = ""
    sco = []
    allGenes = []
    score = 0
    while True :
        line = fi.readline()
        if not line : break
        if re.search("cds",line) :
            t = line.strip().split("\t")
            ch = t[0]
            l = t[:]
            if l[6]=="-" :
                ol[3] = t[4]
                l[4] = t[3]
            genes.append(l)
            l[2] = "CDSpart"
            l[8] = "src=P;grp=" + name + ";pri=4"
            fo.write("\t".join(l))
            fo.write("\n")
        elif re.search("intron",line) :
            t = line.strip().split("\t")
            l = t[:]
            if l[6]=="-" :
                l[3] = t[4]
                l[4] = t[3]
            l[8] = "src=P;grp=" + name + ";pri=4"
            fo.write("\t".join(l))
            fo.write("\n")
        elif re.search("match",line) :
            t = line.strip().split("\t")
            score = float(t[5])
        elif line[0]=="/" :
            sco.append(score)
            allGenes.append(genes)
            genes=[]

    fi.close()
    fo.close()
    if len(sco)== 1 :
        genes = allGenes[0]
    elif len(sco) > 1 :
        maxI = 0
        print "mutiple genewise alignments\t", name
        for i in range(len(sco)-1) :
            if sco[i+1] > sco[i] :
                maxI = i+1
        genes = allGenes[maxI]
    else :
        genes = []
    #print genes
    fi = open(faFile,"r")
    id = ""
    seq = ""
    seqFa = {}
    while True :
        line = fi.readline()
        if not line : break
        if line[0]==">" :
            if id== "" :
                id = line.strip()[1:len(line.strip())]
            else :
                id = line.strip()[1:len(line.strip())]
                seqFa[id]=seq
        else :
            seq += line.strip()
    seqFa[id]=seq
    fi.close()
    flag = "None"
    isATG = 1
    protLeng = 0
    if len(genes) != 0 :
        fi = open(proFile,"r")
        while True :
            line = fi.readline()
            if not line : break
            if line[0]!=">" :
                if line.strip()[0] != "M" :
                    isATG = 0
                    print "No-ATG start...", name                    
                protLeng += len(line.strip())
        fi.close()
        gwLeng = 0
        for k in range(len(genes)) :
            gwLeng += int(genes[k][4]) - int(genes[k][3]) + 1
        if float(gwLeng)/float(protLeng*3) > 0.5 :
            flag = checkGeneWise(genes, seqFa, isATG)
        else :
            print "genewise alignment cov < 50% of the protein", name
            return geneLines
    else :
        print "no GeneWise alignment\t",name
        return geneLines
    if flag == "good" :
        geneLines = getGeneWiseGene(genes, name,shiftVal)
        print "Good genewise:\t" + name + "\tisATG=" + str(isATG) + "\tshiftVal:"+str(shiftVal)+"\n"
        if name == "AT2G12905" :
            print geneLines
            #sys.exit()
        return geneLines
    elif flag.split("\t")[-1] == "stopEnd" :
        stop_1 = flag.split("\t")[-2]
        strand = genes[0][6]
        #if strand == "+"
        cdsIdx = 0
        genesNew = []
        pile = 0
        if strand == "+" :
            for i in range(len(genes)) :
                if pile < int(stop_1) and int(stop_1) < pile + int(genes[i][4]) - int(genes[i][3]) + 1 :
                    genesNew.append(genes[i])
                    genesNew[-1][4]= str( int(genesNew[-1][3]) + ( int(stop_1) - int(pile) ) - 1  )
                    break
                else :
                    pile += int(genes[i][4]) - int(genes[i][3]) + 1
                    genesNew.append(genes[i])
        else :
            # Edited by MG
            first = True
            for i in reversed(range(len(genes))):
                if pile + int(genes[i][4]) - int(genes[i][3]) + 1 <= int(stop_1):
                    pile += int(genes[i][4]) - int(genes[i][3]) + 1
                    continue
                elif first:
                    genesNew.append([j for j in genes[i]])
                    genesNew[-1][3]= str( int(genesNew[-1][3]) + ( int(stop_1) - int(pile) ) + 3  )
                    pile += int(genes[i][4]) - int(genes[i][3]) + 1
                    first = False
                else:
                    genesNew.append(genes[i])
            genesNew.reverse()

            #         # genesNew = genes[:]
            # for i in reversed(range(len(genes))) :
            #     if pile < int(stop_1) and int(stop_1) < pile + int(genes[i][4]) - int(genes[i][3]) + 1 :
            #         genesNew.append(genes[i])
            #         # genesNew[-1][3]= str( int(genesNew[-1][3]) + ( int(stop_1) - int(pile) ) + 3  )   # Changed by MG
            #         genesNew[-1][3]= str( int(genesNew[-1][3]) + int(stop_1) + 3  )
            #         #stopId = i
            #         break
            #     else :
            #         pile += int(genes[i][4]) - int(genes[i][3]) + 1
            #         # genesNew.remove(genesNew[i])  # Changed by MG
            #         genesNew.append(genes[i])

        geneLines = getGeneWiseGene(genesNew, name,shiftVal)
        print geneLines
        #sys.exit()
        return geneLines

    else :
        augCmd = "augustus --UTR=on --print_utr=on --exonnames=on --codingseq=on  --genemodel=complete --alternatives-from-evidence=true --gff3=on --species=arabidopsis "
        outfile3 = outdir + "/" + name + ".augustus.genewise.gff"
        strand = flag.split("\t")[-2]
        if strand == "+" :
            augCmd = augCmd + " --strand=forward "
        elif strand == "-" :
            augCmd = augCmd + " --strand=backward "
        else :
            augCmd = augCmd + " --strand=both "


        augCmd = augCmd + " --extrinsicCfgFile=" + augCfg + " --hintsfile=" + outfile2 + " --outfile=" + outfile3 + " " + faFile
        os.system(augCmd)
        [genes1,genes2]= getAugEviGene(outfile3)
        if genes2.has_key(ch) and len(genes2[ch]) == 1:
            k = genes2[ch].keys()
            #geneLines = "\n".join(genes2[ch][k[0]])
            lines = genes2[ch][k[0]]
            linesNew = []
            if int(lines[0].split("\t")[4]) < 250 :
                print "bad prediction in region:\t" + name
                geneLines = ""
                return geneLines
            for i in range(len(lines)) :
                y = lines[i].split("\t")
                z = y[:]
                # z[0]= z[0].split("_")[0] # Change by MG
                z[0]= "_".join(z[0].split("_")[:-2])
                z[3]= str(int(shiftVal) + int(y[3]) - 200 -1 )
                z[4]= str(int(shiftVal) + int(y[4]) - 200 -1 )
                if i == 0 :
                    z[8] = "ID=" + "AUG-" + name
                elif i== 1 :
                    z[8] = "ID=" + "AUG-" + name + ".t1;Parent=AUG-" + name
                else :
                    z[8] = "Parent=AUG-" + name + ".t1"
                linesNew.append("\t".join(z))
            geneLines = "\n".join(linesNew)
            print "Good Augustus based on genewise hints\n"
            return geneLines
            #print geneLines
            #print "\n"
        elif not genes2.has_key(ch) :
            #if name == "AT1G27535" :
            #    print flag,"\n"
            if flag.split("\t")[-1] == "xx;;stop-lost" :
                strand = genes[0][6]

                isFind = 0
                if strand == "+" :
                    print "Gene stop-lost forward : " , name
                    end = int(genes[-1][4])

                    for i in range(30) : ### extend to find stop codon
                        #print seq[ end+i*3:end+i*3 + 3]
                        if end + i*3 + 3 >= len(seq) :
                            geneLines = ""
                            return geneLines
                        if seq[ end+i*3:end+i*3 + 3]== "TAG" or seq[ end+i*3:end+i*3 + 3]== "TGA"  or seq[ end+i*3:end+i*3 + 3]== "TAA" :
                            print "Find new stop_codon at forword ", str(end+i*3), " although stop-lost", name
                            genes[-1][4] = str(end + i*3)
                            isFind = 1
                            break
                else :
                    end = int(genes[-1][3])
                    for i in range(30) : ### extend to find stop codon
                        if end - i*3 - 3 <=0 :
                            geneLines = ""
                            return geneLines
                        if seq[ end - i*3 - 4 : end - i*3 - 1]== "CTA" or seq[end - i*3 - 4 : end - i*3 - 1]== "TCA"  or seq[ end - i*3 - 4 : end - i*3 - 1]== "TTA" :
                            print "Find new stop_codon at reverse ", str(end-i*3), " although stop-lost", name
                            genes[-1][3] = str(end - i*3)
                            #gwLeng = 0
                            #for k in len(genes) :
                            #    gwLeng += int(genes[k][4]) - int(genes[i][3]) + 1
                            #if gwLeng/len(seq)
                            isFind = 1
                            break
                if isFind == 0 :
                    geneLines = ""
                    return geneLines
                else :
                    geneLines = getGeneWiseGene(genes, name,shiftVal)
                    print geneLines

            elif flag.split("\t")[-1] == "xx;;start-lost" :
                isFind = 0
                if strand == "+" :
                    start = int(genes[0][3])
                    for i in range(30) : ### extend to find stop codon
                        if start - i*3 - 3 <= 0 :
                            geneLines = ""
                            # return geneLines
                        if seq[ start - i*3 -4 : start - i*3 - 1 ]== "TAG" or seq[ start - i*3 -4 : start - i*3 - 1 ]== "TGA" or seq[ start - i*3 -4 : start - i*3 - 1 ]== "TAA" :
                            geneLines = ""
                            # return geneLines
                        if seq[ start - i*3 -4 : start - i*3 - 1 ]== "ATG":
                            print "Find new start_codon at forword ", str(start - i*3 -3), " although start-lost", name
                            genes[0][3] = str(start - i*3 -3)
                            isFind = 1
                            break
                else :
                    start = int(genes[0][4])
                    for i in range(30) : ### extend to find stop codon
                        if start + i*3 + 3 >=len(seq) :
                            geneLines = ""
                            return geneLines
                            break
                        if seq[ start + i*3  : start + i*3 + 3]== "CTA"  or seq[ start + i*3  : start + i*3 + 3]== "TCA"  or seq[ start + i*3  : start + i*3 + 3]== "TTA"  :
                            geneLines = ""
                            return geneLines
                            break
                        if seq[ start + i*3  : start + i*3 + 3]== "CAT" :
                            print "Find new stop_codon at reverse ", str(start+i*3+1), " although start-lost", name
                            genes[-1][4] = str(start + i*3)
                            isFind = 0
                            break

                if isFind == 0 :
                    geneLines = ""
                    return geneLines
                else :
                    geneLines = getGeneWiseGene(genes, name,shiftVal)
                    print geneLines
                    return geneLines

            else :
                print "No gene in augustus-genewise", name, "\n"
                geneLines = ""
            return geneLines

        else :
            print "more then one gene in augustus-genewise", name, "\n"
            geneLines = ""
            return geneLines




def getGeneIndex (geneBed):

    fi = open(geneBed,"r")
    geneIdx = []
    n = 0
    genePos = {}
    while True :
        line = fi.readline()
        if not line : break
        #if (not re.search("Note=protein_coding_gene",line) ) and (not re.search("locus_type=protein",line)) : continue
        t = line.strip().split("\t")
        id = t[3]
        geneIdx.append(id)
        genePos[id]=n
        n += 1
    fi.close()

    return [geneIdx,genePos]


if __name__ == "__main__":
   main(sys.argv[1:])






