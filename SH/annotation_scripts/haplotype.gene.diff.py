#!/usr/bin/env python
# encoding: utf-8


import sys
import os
import getopt
import re
import datetime

def main(argv):
    [grFile, svFile, refGffFile, altGffFile, outdir, refP, altP, refLenFile, altLenFile, refBlnFile, altBlnFile] = ["","","","","","","","","","",""]
    
    try:
        opts, args = getopt.getopt(argv,"g:s:r:a:o:x:y:m:n:j:k:",["gr=","sv=","refG=","altG=","outd=","refP=", "altP=", "refL=", "altL=","refB=","altB="])
    except getopt.GetoptError:
        print('haplotype.gene.diff.py -g <gr> -s <sv> -r <refG -a <altG> -o <outd> -x <refP> -y <altP -m <refL> -n <altL>')
        sys.exit(2)
    if len(opts) == 0 :
        print('haplotype.gene.diff.py -g <gr> -s <sv> -r <refG> -a <altG> -o <outd>  -x <refP> -y <altP -m <refL> -n <altL>')
        sys.exit()
    for opt, arg in opts:
        if opt == '-h':
            print('haplotype.gene.diff.py -g <gr> -s <sv> -r <refG> -a <altG> -o <outd>  -x <refP> -y <altP -m <refL> -n <altL>')
            sys.exit()
        elif opt in ("-g", "--gr"):
            grFile = arg
        elif opt in ("-s", "--sv"):
            svFile = arg
        #elif opt in ("-b", "--blp"):
        #    blpFile = arg
        elif opt in ("-r", "--refG"):
            refGffFile = arg
        elif opt in ("-a", "--altG"):
            altGffFile = arg
        elif opt in ("-m", "--refL"):
            refLenFile = arg
        elif opt in ("-n", "--altL"):
            altLenFile = arg
        elif opt in ("-x", "--refP"):
            refP = arg
        elif opt in ("-y", "--altP"):
            altP = arg
        elif opt in ("-j", "--refP"):
            refBlnFile = arg
        elif opt in ("-k", "--altP"):
            altBlnFile = arg
        elif opt in ("-o", "--outd") :
            outdir = arg
        
    startTime = datetime.datetime.now()
    print("Starting ======  ",startTime.strftime('%Y-%m-%d %H:%M:%S'), "  ======")
    print(refGffFile, altGffFile)
    ## gene model from augustus protein-based, scipio protein-based or stringTie-gff based
    if not os.path.isdir(outdir) :
        os.system("mkdir -p " + outdir)
    
    refGeneBed = outdir + "/ref.genes.bed"
    refGenes = getGeneCoord(refGffFile, refGeneBed)
    altGeneBed = outdir + "/alt.genes.bed"
    altGenes = getGeneCoord(altGffFile, altGeneBed)
 
    ## step 1: gene CNV stats
    getGrCNV (grFile, outdir, refGenes, altGenes, refP, altP)   
    
    ## step 2: get syri output, annotate functional effect of snps and small indels
    syriParser (svFile,refLenFile, altLenFile, outdir)
    
    inFile = outdir + "/ref-coord.snp.indel.out"
    outfile1 = outdir + "/ref-coord.snp.indel.vcf"
    svToVCF(inFile, outfile1)
    refEffOutFile = outdir + "/ref-coord.snp.indel.snpEff.vcf"
    runSNPeff(outfile1, refP.lower(), refEffOutFile)
    refGeneEff = getVarEff(refEffOutFile,refP)
    #refGeneEff = {}
    
    inFile = outdir + "/alt-coord.snp.indel.out"
    outfile2= outdir + "/alt-coord.snp.indel.vcf"
    svToVCF(inFile, outfile2)
    outfile2 = outdir + "/alt-coord.snp.indel.vcf"
    altEffOutFile = outdir + "/alt-coord.snp.indel.snpEff.vcf"
    runSNPeff(outfile2,altP.lower(), altEffOutFile)
    altGeneEff = getVarEff(altEffOutFile,altP)
    #altGeneEff = {}
    
    ## step 3: check haplotype specific genes
    refSpFile = outdir + "/ref.specific.gene.bed"
    delRef = outdir + "/ref-coord.alt.long.del.out"
    outp = outdir + "/ref"
    refBln = getBln(refBlnFile) #[chr,start,end,iden,cov]
    analyzeHapSpGene (refSpFile, delRef, refGeneEff, refBln, outp) 
    infile = outp + ".spGene.noDel.noLoF.blastn.bed"
    outfile = outp + ".spGene.noDel.noLoF.blastn.intersect.alt-genes.bed"
    statFile1 = outp + ".spGene.noDel.noLoF.blastn.stats"
    checkAnn (infile, altGeneBed, outfile, statFile1) 
    
    altSpFile = outdir + "/alt.specific.gene.bed"
    delAlt = outdir + "/alt-coord.ref.long.del.out"
    outp = outdir + "/alt"
    altBln = getBln(altBlnFile)
    analyzeHapSpGene (altSpFile, delAlt, altGeneEff, altBln, outp) 
    infile = outp + ".spGene.noDel.noLoF.blastn.bed"
    outfile = outp + ".spGene.noDel.noLoF.blastn.intersect.ref-genes.bed"
    statFile2 = outp + ".spGene.noDel.noLoF.blastn.stats"
    checkAnn (infile, refGeneBed, outfile, statFile2) 
    
    ## haplotype Copy diff genes [SYN, INV, ITX, CTX, DUP-A, DUP-B]
    
    ## haplotype Copy same genes [SYN, INV, ITX, CTX, DUP-A, DUP-B]
    
def getBln(inFile) :
    fi = open(inFile,"r")
    bln = {}
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        t[0] = t[0].split(".")[0]
        bln[t[0]] = [t[1], t[6], t[7], t[8], t[9]]
    fi.close()
    return bln
    
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
    
    
    '''find gene gain/loss and copy gain/loss between query species and subject species  '''
def getGrCNV (grFile, outdir, refGenes, altGenes,ref, alt):
    #[grNvN,gr1v1,gr2v2,gr3v3,grMvM] = [0,0,0,0,0]
    #[grNv0,gr1v0,gr2v0,gr3v0,grMv0] = [0,0,0,0,0]
    #[gr0vN,gr0v1,gr0v2,gr0v3,gr0vM] = [0,0,0,0,0]
    #[grMvN,gr2v1,gr3v1,gr3v2,grMvX] = [0,0,0,0,0]
    #[grNvM,gr1v2,gr1v3,gr2v3,grXvM] = [0,0,0,0,0]
    
    stats = [[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]]
    fi = open(grFile,"r")
    fo1 = open(outdir + "/gene.cnv.stats", "w")
    fo2 = open(outdir + "/gene.copy.diff.groups.txt", "w")
    fo5 = open(outdir + "/gene.copy.same.groups.txt", "w")
    fo3 = open(outdir + "/ref.specific.gene.bed", "w")
    fo4 = open(outdir + "/alt.specific.gene.bed", "w")
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split()
        [refN, altN] = [0, 0]
        [refId, altId] = [[],[]]
        for k in range(len(t)) :
            if k == 0 : continue
            if re.search(ref,t[k]) :
                refN += 1
                refId.append(t[k])
            else :
                altN += 1
                altId.append(t[k])
        if refN == 1 and altN == 1 :
            stats[0][0] += 1
            stats[0][1] += 1
        elif refN == 2 and altN == 2:
            stats[0][0] += 1
            stats[0][2] += 1
        elif refN == 3 and altN == 3:
            stats[0][0] += 1
            stats[0][3] += 1
        elif refN > 3 and refN == altN :
            stats[0][0] += 1
            stats[0][4] += 1
        elif refN > altN :
            if refN == 1 and altN == 0 :
                stats[1][0] += 1
                stats[1][1] += 1
            elif refN == 2 and altN == 0:
                stats[1][0] += 1
                stats[1][2] += 1
            elif refN == 3 and altN == 0:
                stats[1][0] += 1
                stats[1][3] += 1
            elif refN > 3 and altN == 0:
                stats[1][0] += 1
                stats[1][4] += 1
            elif refN == 2 and altN == 1:
                stats[3][0] += 1
                stats[3][1] += 1
            elif refN == 3 and altN == 1:
                stats[3][0] += 1
                stats[3][2] += 1 
            elif refN == 3 and altN == 2:
                stats[3][0] += 1
                stats[3][3] += 1
            elif refN > 3 :
                stats[3][0] += 1
                stats[3][4] += 1
        elif refN < altN:
            if altN == 1 and refN == 0 :
                stats[2][0] += 1
                stats[2][1] += 1
            elif altN == 2 and refN == 0:
                stats[2][0] += 1
                stats[2][2] += 1
            elif altN == 3 and refN == 0:
                stats[2][0] += 1
                stats[2][3] += 1
            elif altN > 3 and refN == 0:
                stats[2][0] += 1
                stats[2][4] += 1
            elif altN == 2 and refN == 1:
                stats[4][0] += 1
                stats[4][1] += 1
            elif altN == 3 and refN == 1:
                stats[4][0] += 1
                stats[4][2] += 1 
            elif altN == 3 and refN == 2:
                stats[4][0] += 1
                stats[4][3] += 1
            elif altN > 3 :
                stats[4][0] += 1
                stats[4][4] += 1
        t[0]=t[0][0:len(t[0])-1]
        if refN > 0 and altN == 0:
            for g in refId :
                fo3.write("{}\t{}\t{}\t{}\t{}\n".format("\t".join(refGenes[g]), g, t[0], refN, altN ) ) 
        elif refN == 0 and altN > 0:
            for g in altId :
                fo4.write("{}\t{}\t{}\t{}\t{}\n".format("\t".join(altGenes[g]), g, t[0], altN, refN ) ) 
        elif refN != altN :
            fo2.write(line)	 
        elif refN == altN :
            fo5.write(line)
    
    fo2.close()
    fo3.close()
    fo4.close()
    fo5.close()
    fi.close()
    
    os.system("sort -k1,1 -k2,2n -k3,3n " + outdir + "/ref.specific.gene.bed > " + outdir + "/ref.specific.gene.srt.bed")
    os.system("sort -k1,1 -k2,2n -k3,3n " + outdir + "/alt.specific.gene.bed > " + outdir + "/alt.specific.gene.srt.bed")
    for i in range(len(stats)) :
        fo1.write("{}\n".format("\t".join( str(x) for x in stats[i] ) ))	
    fo1.close()
    #return [geneLossS,geneLossQ,geneGainS,geneGainQ,copyLossS,copyLossQ,copyGainS,copyGainQ]
  
def analyzeHapSpGene (spGeneFile, delFile, geneEff, blastn, outp) :
    ##check del
    outfile1 = outp + ".spGene.del.intersect.out"
    cmd = "intersectBed -a " + spGeneFile + " -b " + delFile + " -wao > " + outfile1 
    print(cmd)
    os.system(cmd)
    fi = open(outfile1, 'r')
    outfile2 = outp + ".spGene.del.lof.out"
    fo = open(outfile2,"w")
    delGenes = {}
    minOlap = 5
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        if int(t[-1]) <= minOlap : continue
        olap = float(t[-1])/(float(t[2]) - float(t[1]))
        olap = '{:.2f}'.format(olap)
        #t.append(str(olap))
        #fo.write("{}\n".format("\t".join(t)))
        f = ""
        if int(t[1]) >= int(t[9]) and int(t[2]) <= int(t[10]) : ##100% del
            f = "d100"
        elif float(t[-1]) >= 0.9*(float(t[2]) - float(t[1])) : ## 90% del
            f = "d90"
        elif float(t[-1]) >= 0.6*(float(t[2]) - float(t[1])) : ## 60% del
            f = "d60"
        elif float(t[-1]) >= 0.3*(float(t[2]) - float(t[1])) : ## 30% del
            f = "d30"
        else :
            f = "d00"
        delGenes[t[4]] = str(f) + "-"+ str(olap)
    fi.close()
    
    ##check LoF
    fi = open(spGeneFile, "r")
    spGenes = {}
    lofGenes = {}
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        spGenes[t[4]] = t
        if t[4] in geneEff.keys() :
            lofGenes[t[4]] = geneEff[t[4]]
            #t.append("LoF-" + geneEff[t4])
            #fo.write("{}\n".format("\t".join(t)))
    fi.close()
    
    
    ##check Ann
    outfile5 = outp + ".spGene.noDel.blastn.bed"
    outfile3 = outp + ".spGene.noDel.noLoF.blastn.bed"
    outfile4 = outp + ".spGene.noDel.noLoF.noblastn.bed"
    outfile6 = outp + ".spGene.blastn.bed"
    fo1 = open(outfile3,"w")
    fo2 = open(outfile4,"w")
    fo5 = open(outfile5,"w")
    fo6 = open(outfile6,"w")
    for gene in sorted(spGenes.keys()) :
        ll = "**"
        if gene in lofGenes.keys():
            lofs = []
            for lof in lofGenes[gene]:
                #print(lof)
                lofs.append(":".join(lof))
            ll = ";".join(lofs)
        fo6.write("{}\t{}\t{}\n".format("\t".join(blastn[gene]),"\t".join(spGenes[gene]), ll ))
        if gene in delGenes and gene in lofGenes:
            dd = delGenes[gene]
            ll = lofGenes[gene][-1][-1]
            #print(ll)
            fo.write("{}\t{}\t{}\n".format("\t".join(spGenes[gene]), dd, ll ))
        elif gene in delGenes :
            dd = delGenes[gene]
            ll = "--"
            fo.write("{}\t{}\t{}\n".format("\t".join(spGenes[gene]), dd, ll ))
        elif gene in lofGenes:
            dd = blastn[gene][3] + "-" + blastn[gene][4]
            ll = lofGenes[gene][-1][-1]
            #print(ll)
            fo.write("{}\t{}\t{}\n".format("\t".join(spGenes[gene]), dd, ll ))
            if gene in blastn.keys() :
                if float(blastn[gene][3]) > 80 and float(blastn[gene][4]) > 80:
                    fo5.write("{}\t{}\t{}\t{}\n".format("\t".join(blastn[gene]), "\t".join(spGenes[gene]), dd, ll ))
        elif gene in blastn.keys() and blastn[gene][0] == "None" :
            dd = "bln-None"
            ll = "--"
            fo.write("{}\t{}\t{}\n".format("\t".join(spGenes[gene]), dd, ll ))
        elif gene in blastn.keys() :
            #[ch, start, end, cov, iden] = blastn[gene]
            if float(blastn[gene][3]) < 80 or float(blastn[gene][4]) < 80 :
                dd = blastn[gene][3] + "-" + blastn[gene][4]
                ll = "--"
                fo.write("{}\t{}\t{}\n".format("\t".join(spGenes[gene]), dd, ll ))
            else :
                fo1.write("{}\t{}\n".format("\t".join(blastn[gene]), "\t".join(spGenes[gene]) ))
                fo5.write("{}\t{}\t{}\t{}\n".format("\t".join(blastn[gene]),"\t".join(spGenes[gene]), dd, ll ))
        else :
            fo2.write("{}\n".format("\t".join(spGenes[gene])))
            
    fo1.close()
    fo2.close()
    fo.close()
    fo5.close()
    fo6.close()
    os.system("sort -k1,1 -k2,2n -k3,3n " + outfile3 + " > " + outp + ".spGene.noDel.noLoF.blastn.srt.bed")
    os.system("sort -k1,1 -k2,2n -k3,3n " + outfile5 + " > " + outp + ".spGene.noDel.blastn.srt.bed")
    return
    
    
def checkAnn (inFile, annBed, outfile, statFile) :
    #inFile : chr start end  cov iden chr start end +/- gene gr ref-N alt-N  
    cmd = "intersectBed -a " + inFile + " -b " + annBed + " -wao > " + outfile
    print(cmd)
    os.system(cmd)
    fi = open(outfile,"r")
    print(outfile)
    annStat = {}
    st = [0,0,0,0]
    while True:
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        if t[-1] == "0" :
            if t[9] in annStat.keys() : continue
            annStat[t[9]] = "Ref-No-Ann"
            st[0] += 1
        elif float(t[3]) < 95: 
            if t[9] in annStat.keys() : continue
            annStat[t[9]] = "par-del"
            st[1] += 1
        elif float(t[4]) < 95 :
            if t[9] in annStat.keys() : continue
            #print(t[9])
            annStat[t[9]] = "HDR"
            st[2] += 2
        else :
            if t[9] in annStat.keys() : continue
            annStat[t[9]] = "misAnn"
            st[3] += 3
    fi.close()
    fo = open(statFile, "w")
    fo.write("Ref-NoAnn\t{}\nPar-Del\t{}\nHDR\t{}\nmisAnn\t{}\n".format(st[0],st[1],st[2],st[3]))
    fo.close()
    return annStat
    
def runSNPeff (svVcf, genome, outfile) :
    snpEff = "java -Xmx4G -jar /srv/netscratch/dep_mercier/grp_schneeberger/bin/snpEff/snpEff4.3p/snpEff/snpEff.jar"
    cmd = snpEff + " -noStats -formatEff -q -lof -no-downstream -no-upstream -no-intergenic " + genome + " " + svVcf + " > " + outfile
    print(cmd)
    #os.system(cmd)
    return

def svToVCF(inFile, outFile) :
    fi = open(inFile,"r")
    fo = open(outFile,"w")
    fo.write("##fileformat=VCFv4.2\n")
    fo.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t/xx.bam\n")
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        fo.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(t[0], t[1], ".", t[3], t[4], 40, ".", "DP=10", "GT", "0/0"))
    fi.close()
    fo.close()
    
def getVarEff (inFile,p) :  #output from snpEff
    fi = open(inFile, "r")
    effStats = {} #gene pos ref alt effect
    #return effStats 
    while True :
        line = fi.readline()
        if not line : break
        if line[0] == "#" : continue
        
        t = line.strip().split("\t")
        if not re.search(r"EFF", t[7]) : continue
        
        effStr = t[7].split(";")[1]
        effs = effStr.split("=")[1].split(",")
        for eff in effs:
            if not re.search(r'HIGH', eff) : continue
            mm =  re.search(r'(\dG\d+)\|',eff)
            gene = mm.groups()[0]
            #print(gene)
            ty = eff.split("(")[0]
            if gene in effStats.keys() :
                gene = p + "-" + gene 
                effStats[gene].append([gene, t[1], t[3], t[4], ty])
            else :
                gene = p + "-" + gene
                effStats[gene] = []
                effStats[gene].append([gene, t[1], t[3], t[4], ty])
    fi.close()
    print("total number of genes with LoF variants{}\n".format(len(effStats.keys())))
    return effStats 
    
def syriParser (inFile,refLenFile, altLenFile, outdir):
    ##
    outfile1 = outdir + "/ref-coord.snp.indel.out"
    outfile2 = outdir + "/alt-coord.snp.indel.out"
    fo1 = open(outfile1, "w")
    fo2 = open(outfile2, "w")
    
    refBlk = outdir + "/ref-coord.align.block.out"
    fo3 = open(refBlk, "w")
    altBlk = outdir + "/alt-coord.align.block.out"
    fo4 = open(altBlk, "w")
    fi = open(inFile,"r")
    delRef = outdir + "/ref-coord.alt.long.del.out"
    delAlt = outdir + "/alt-coord.ref.long.del.out"
    fo5 = open(delRef, "w")
    fo6 = open(delAlt, "w")
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        if t[10] != "NOTAL" and int(t[6]) > int(t[7]) :
            [t[7], t[6]] = [t[6],t[7]]
        ttAlt = [t[5], t[6], t[7], t[0], t[1], t[2], t[8], t[9], t[10]]
        ttRef = [t[0], t[1], t[2], t[5], t[6], t[7], t[8], t[9], t[10]]
        
        if t[10] == "SNP" :  
            ttRef = [t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8], t[9], t[10]]
            ttAlt = [t[5], t[6], t[7], t[4], t[3], t[0], t[1], t[2], t[8], t[9], t[10]]
            fo1.write("{}\n".format("\t".join(ttRef)))
            fo2.write("{}\n".format("\t".join(ttAlt)))
        elif t[10] == "INS" :
            ttRef = [t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8], t[9], t[10]]
            ttAlt = [t[5], t[6], t[7], t[4], t[3], t[0], t[1], t[2], t[8], t[9], t[10]]
            ttAlt[-1] == "DEL"
            fo1.write("{}\n".format("\t".join(ttRef)))
            fo2.write("{}\n".format("\t".join(ttAlt)))
        elif t[10] == "DEL" :
            ttRef = [t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8], t[9], t[10]]
            ttAlt = [t[5], t[6], t[7], t[4], t[3], t[0], t[1], t[2], t[8], t[9], t[10]]
            ttAlt[-1] == "INS"
            fo1.write("{}\n".format("\t".join(ttRef)))
            fo2.write("{}\n".format("\t".join(ttAlt)))
        elif re.search(r'AL', t[10]) and t[10] != "NOTAL":
            fo3.write("{}\n".format("\t".join(ttRef)))
            fo4.write("{}\n".format("\t".join(ttAlt)))
        elif t[10] == "NOTAL" :
            if t[0] != "-" :
                nal = "\t".join(t[0:3])
                fo5.write(nal)
                fo5.write("\n")
            else :
                nal = "\t".join(t[5:8])
                fo6.write(nal)
                fo6.write("\n")
            
    fi.close()
    fo1.close()
    fo2.close()
    fo3.close()
    fo4.close()
    fo5.close()
    fo6.close()
    
    
    outfile11 = outdir + "/ref-coord.snp.indel.srt.out"
    outfile21 = outdir + "/alt-coord.snp.indel.srt.out"
    refBlkSrt = outdir + "/ref-coord.align.block.srt.out"
    altBlkSrt = outdir + "/alt-coord.align.block.srt.out"
    
    os.system("sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n " + outfile1 + " > " + outfile11)
    os.system("sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n " + outfile2 + " > " + outfile21)
    
    os.system("sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n " + refBlk + " > " + refBlkSrt)
    os.system("sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n " + altBlk + " > " + altBlkSrt)
    
    cmd = "bedtools complement -i " + refBlkSrt + " -g " + refLenFile + " > " + delRef
    print(cmd)
    #os.system(cmd)
    cmd = "bedtools complement -i " + altBlkSrt + " -g " + altLenFile + " > " + delAlt
    print(cmd)
    return
    
def getGeneCoord (inFile,outFile) : #input: annotation gff file
    geneCoo = {}
    print(inFile)
    fi = open(inFile,"r")
    fo = open(outFile,"w")
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        if t[2] == "gene" :
            id = t[8].split(";")[0].split("=")[1]
            #print(id)
            #sys.exit()
            geneCoo[id] = [t[0], t[3], t[4], t[6]]
            fo.write("{}\n".format("\t".join([t[0], t[3], t[4], t[6], id])))
    fi.close()
    fo.close()
    return geneCoo
    
def checkCopySame(geneBedFile, blkBed, outfile) :
    #geneBedFile : chr start end +/- id group ref-N alt-N
    #refBlkBed : chr start end chr start end SV-type
    cmd = "intersectBed -a " + geneBedFile + " -b " + blkBed + " -wao >" + outfile
    fi = open(outfile, "r")
    geneBlk = {}
    geneOlap = {}
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        ## inside a block
        if t[2] >= t[9] and t[3] <= t[10] :
            if t[4] not in geneBlk.keys() :
                geneBlk[t[4]] = t[14]
            elif (geneBlk[t[4]] == "DUPAL" and (t[14] == "SYNAL" or t[14] == "INVAL" 
                  or t[14] == "TRANSAL" or t[14] == "INVTRAL") ):
                geneBlk[t[4]] = t[14]
            elif (geneBlk[t[4]] == "INVDPAL" and (t[14] == "SYNAL" or t[14] == "INVAL" 
                  or t[14] == "TRANSAL" or t[14] == "INVTRAL") ):
                geneBlk[t[4]] = t[14]
            else :
                continue
        else :
            olap = float(t[-1])/(float(t[3]) - float(t[2]))
            if t[4] in geneOlap.keys() :
                #
                continue
            else :
                geneOlap[t[4]] = olap
        ## 
        
    fi.close()
    
    return

def checkCopyDiff(geneBedFile, blkBed, outfile):

    return
    
if __name__ == "__main__":
    main(sys.argv[1:])






