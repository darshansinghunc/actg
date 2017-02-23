
do_index='zf_graphs'
do_index='pa_u937A'
do_index='pa_u937A2'
do_index='pa_u937A'
do_index='multitx'
do_index='zf_graphs'
do_index='zf_html'
do_index='zf_graphs'
do_index='test1'
do_index='testfdm2'
do_index='mergebam'
do_index='combp'
do_index='mergebam4'
do_index='confighelper'
do_index='testwgsdict'
do_index='testactcorrect'
do_index='testdom'
do_index='combp2'
do_index='testpval'
do_index='paul_dom'
do_index='davidplot'
do_index='davidmaqdata'
#do_index='davidACTdata'
do_index='davidplot'
do_index='davidmaqdata'
do_index='davidplot'
do_index='davidACTdata'
do_index='testactcorrect'
do_index='davidACTdata'
do_index='testactcorrect'
do_index='davidACTdata'
do_index='davidplot'
do_index='thesis_pdfs' 
do_index='testpval'
do_index='T1JS'
do_index='KATE1'
do_index='KATEbox1'
do_index='KATEbox2'
do_index='expchange'
do_index='ATHZ1'
do_index='ATHZ2'
do_index='ATHZ3'
do_index='ATHZ2'
do_index='ATHZ4'
do_index='ATHZ2'
do_index='ATHZ5'
do_index='ATHZ6'
do_index='ATHZ7'
do_index='ATHZ8'
do_index='ATHZ7'
do_index='ATHZ8'
do_index='ATHZ5'
do_index='ATHZ9'
do_index='ATHZ5'
do_index='ATHZ9'
do_index='ATHZ5'
do_index='ATHZ10'
do_index='ATHZ5'
do_index='ATHREP1'
do_index='ATHZ2'
do_index='ATHZ5'
do_index='ATHZ10'
do_index='ATHZ2'
do_index='ATHZ5'
do_index='ATHZ10'
do_index='ATHZ2'
do_index='ATHZ5'
do_index='ATHZ2'
do_index='ATHZ5'
do_index='ATHZ2'
do_index='ATHZ5'
do_index='ATHZ10'
do_index='ATHZ2'
do_index='gtf1'
do_index='ATHimage'
do_index='a'
do_index='ATHZ5_9'
do_index='ATHREP2'
do_index='ATHZ2'
do_index='ATHZ5'
do_index='ATHZ11'
do_index='ATHZ2'
do_index='ATHZ5' 
do_index='ATHZ2'
do_index='ATHZ5' 
do_index='ATHZ11a'
do_index='ATHZ9'
do_index='ATHZ5'
do_index='ATHZ11a'
do_index='ATHZ2'
do_index='RMOUSE'

if do_index=='zf_graphs':
    import actgraph
    import ast
    import os 
    
    infile='/playpen/rootFDM/project/ZF_GZ0/ffull/TOP_SIG__data2.txt'
    worklist=[]
    for lntxt in open(infile):
        ln=lntxt.strip('\n').split('\t')
        if ln[0]=='gene':
            continue
        l2=ast.literal_eval(ln[15])
        l2.append(ast.literal_eval(ln[14]))
        minx=min([min(x) for x in l2])
        maxx=max([max(x) for x in l2])
        worklist.append([max(float(ln[2]),float(ln[5]),float(ln[8])),ln[0],int(ln[1]),[minx,maxx]])
    
    worklist.sort()
    worklist.reverse()
    worklist=worklist[1900:2000]
    tslist=['ZF_CTRL1_A','ZF_RAP_A','ZF_TG3_A']
    project_name='ZF_GZ0'
    rootdir='/playpen/rootFDM'
    actfilelist=['%s/dataout/%s/%s_%s.act'%(rootdir,ts,project_name,ts) for ts in tslist]
    imagedir='/playpen/rootFDM/project/%s/report/image'%project_name
    for actfilename in actfilelist:
        act0=actgraph.actFile(actfilename)
        for maxfdm,gene,pos,genomerange in worklist:
            if gene in ['CR623802','OGG1 type 1e']:
                continue
            print maxfdm,gene,pos
            act0.Toimage(gene, imagedir,genomerange,[pos],'pdf')
    
    for maxfdm,gene,pos,genomerange in worklist:
        os.system('montage -geometry 760x1600 -tile 3x1 "%s/*%s*%d-%d.pdf" "%s/all_%s_%s-%s_joined.pdf"'%
                  (imagedir,gene.replace('/','_'),genomerange[0],genomerange[1],imagedir,gene.replace('/','_').replace(' ','_'),genomerange[0],genomerange[1]))

if do_index=='zf_html':
    import actgraph
    import ast
    import os 
    import common
    
    infile='/playpen/rootFDM/project/ZF_GZ0/ffull/TOP_SIG__data2.txt'
    csvfile='/playpen/rootFDM/project/ZF_GZ0/report/TOP_SIG_REPORT.txt'
    fout=open(csvfile,'w')
    htmlfile=open('/playpen/rootFDM/project/ZF_GZ0/report/TOP_SIG_REPORT.html','w')
    outtxtlist=[]
    for lntxt in open(infile):
        ln=lntxt.strip('\n').split('\t')
        if ln[0]=='gene':
            fout.write('%s\t%s\n'%(lntxt.strip(),'act graph'))
            continue
        l2=ast.literal_eval(ln[15])
        l2.append(ast.literal_eval(ln[14]))
        minx=min([min(x) for x in l2])
        maxx=max([max(x) for x in l2])
        imgname='image/all_%s_%d-%d_joined.pdf'%(ln[0],minx,maxx)
        urlname='<a href="%s">act graph</a>'%imgname
        outtxtlist.append([2.0-max(float(ln[2]),float(ln[5]),float(ln[8])),'%s\t%s\n'%(lntxt.strip(),urlname)])
    outtxtlist.sort()
    for maxfdm,txt in outtxtlist:
        fout.write(txt)
    fout.close()
    common.csv2html(open(csvfile),htmlfile)
    htmlfile.close()

    
if do_index=='pa_u937A':
    genelist=['NPM1','RUNX2','HOXA9','GAS6','SF1','SF3B2','SIGLEC6','SMYD3','RAVER1','HNRNPM','BAFF','CEBPA',
              'ELANE','GFI1', 'HOXA9', 'NOTCH2NL','ASXL1','U2AF1']
    genelist.sort()
    tslist=['STEM-U937A']
    project_name='STEM2'
    rootdir='/playpen/rootFDM'
    actfilelist=['%s/dataout/%s/%s_%s.act'%(rootdir,ts,project_name,ts) for ts in tslist]
    imagedir='/playpen/rootFDM/project/%s/image'%project_name
    for actfilename in actfilelist:
        act0=actgraph.actFile(actfilename)
        for gene in genelist:
            print gene
            act0.Toimage(gene, imagedir)
    
    
if do_index=='pa_u937A2':
    import actgraph    
    tslist=['STEM-U937A']
    project_name='STEM2'
    rootdir='/playpen/rootFDM'
    actfilelist=['%s/dataout/%s/%s_%s.act'%(rootdir,ts,project_name,ts) for ts in tslist]
    imagedir='/playpen/rootFDM/project/%s/image'%project_name
    for actfilename in actfilelist:
        act0=actgraph.actFile(actfilename)
        act0.Toimage('CHD4',imagedir,[6679249,6716551])    
        
if do_index=='ATHimage':
    """
            #temphack
        exonlist=[[82031576,82032336],[82032336,82032691],[82032691,82032718],[82032718,82032720],[82032720,82033176],
                  [82033176,82033215],[82033215,82033256],[82033256,82033639],[82034276,82040548],[82041452,82043672],
                  [82043672,82043794],[82045268,82045340],[82045340,82045345],[82049089,82049434]]
        exonwtlist=[1779.31,2076.53,1695.21,1699.76,1285.04,449.62,492.84,700.25,800.6,2.42,726.17,1028.76,1305.32,953.39]
        
        intronlist=[[82033639,82034276]]
        intronwtlist=[5.37]

        splicelist=[[82032336,82034276],[82032718,82033176],[82032720,82033215],[82033639,82034276],[82040548,82041452],
                    [82040548,82043672],[82040548,82045268],[82043794,82045268],[82045340,82049089],[82045345,82049089]]
        
        splicewtlist=[13.0,148.0,10.0,726.0,6.0,805.0,57.0,934.0,8.0,1218.0]
        
        startnodelist=[82031576]
        
        endnodelist=[82049434]
        
        novelnodelist=[82032336,82032718,82033176,82045340]
    """
    import actgraph    
    tslist=['ATHFCD15_T01']
    project_name='ATHFCD15'
    tslist=['STEM-U937A']
    project_name='STEM2'    
    tslist=['ATHFCD052_T01']
    project_name='ATHFCD052'       
    rootdir='/playpen/rootFDM'
    actfilelist=['%s/dataout/%s/%s_%s.act'%(rootdir,ts,project_name,ts) for ts in tslist]
    imagedir='/playpen/rootFDM/project/%s/report/image'%project_name
    print imagedir
    for actfilename in actfilelist:
        act0=actgraph.actFile(actfilename)
        act0.Toimage('MAT1A',imagedir,genomerange=[82031575,82049435])     

if do_index=='multitx':
    import cPickle        
    x=cPickle.load(open('/playpen/rootFDM/annotation/fdm_human_hg19_knowngene.gtf/genetranscriptdict.pck'))
    print len(x)
    y=[a for a in x if len(x[a])>1]
    print len(y)
    for z in y:
        print z

if do_index=='test1':
    import common
    ln=['a','b','c','a:1','a:2','b:1','b:2','c:1','c:2','ab:3','ac:3','bc:3']
    hdict=common.findheaderdict(ln,['3'])
    print hdict
    
    
if do_index=='testfdm2':
    import common
    group1=[1,1,1,1,1,1]
    group2=[0,0,0,1,0,0,0,1,1]
    group1=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
    group2=[0,1,0,0,0,1,0,0,1,1,0,1]    
    print common.getgroupshufflepvalue(group1,group2,numiter=10000)



    
    
if do_index=='mergebam':
    import os
    
    def mergebamfiles(pathsamtools,mergedict):
        for out in mergedict:
            cmd1='mkdir -p %s'%('/'.join(out.split('/')[:-1]))
            print cmd1
            cmd2='%s/samtools merge %s %s'%(pathsamtools,out,' '.join(mergedict[out]))
            print cmd2
            #os.system(cmd1) 
            #os.system(cmd2)
        
    def twoexptojsdfile(in1,in2,out):
        return
    
    def mergesregenmetadata(mergedict,jsd_pair_dict):
        #new expression.txt
        for out in jsd_pair_dict:
            twoexptojsdfile(jsd_pair_dict[out][0],jsd_pair_dict[out][1],out)
    
    pathsamtools='/home/darshan/nextgen/tools/samtools-0.1.8'
    rootdir='/playpen/sregen'
    outdir='%s/TYPE3'%rootdir
    indirs=['%s/TYPE3A'%rootdir,'%s/TYPE3B'%rootdir,'%s/TYPE3C'%rootdir]
    #the first one dictates the set
    dirs=dict([(x,os.listdir('%s/data'%x)) for x in indirs])
    for x in dirs:
        dirs[x]=['%s/data/%s/alignments.bam'%(x,y) for y in dirs[x]]
        dirs[x].sort()
    #print dirs
    allfiles=[ dirs[x] for x in indirs]
    mergefilelist=zip(*allfiles)
    #print mergefilelist
    out=['%s/data/%s/alignments.bam'%(outdir,x.split('/')[-2]) for x in dirs[indirs[0]]]
    #print out
    mergedict=dict(zip(out,mergefilelist))
    print mergedict
    mergebamfiles(pathsamtools,mergedict)
    
if do_index=='mergebam4':
    import os
    
    def mergebamfiles(pathsamtools,mergedict):
        for out in mergedict:
            cmd1='mkdir -p %s'%('/'.join(out.split('/')[:-1]))
            #print cmd1
            cmd2='%s/samtools merge %s %s'%(pathsamtools,out,' '.join(mergedict[out]))
            print cmd2
            os.system(cmd1) 
            os.system(cmd2)            
        
    def twoexptojsdfile(in1,in2,out):
        return
    
    def mergesregenmetadata(mergedict,jsd_pair_dict):
        #new expression.txt
        for out in jsd_pair_dict:
            twoexptojsdfile(jsd_pair_dict[out][0],jsd_pair_dict[out][1],out)
    
    pathsamtools='/home/darshan/nextgen/tools/samtools-0.1.8'
    rootdir='/playpen/sregen'
    outdir='%s/TYPE4'%rootdir
    indirs=['%s/TYPE4%s'%(rootdir,x) for x in ['A','B','C','D','E','F']]
    #the first one dictates the set
    dirs=dict([(x,os.listdir('%s/data'%x)) for x in indirs])
    #print dirs
    for x in dirs:
        dirs[x]=['%s/data/%s/alignments.bam'%(x,y) for y in dirs[x]]
        dirs[x].sort()
    #print dirs
    allfiles=[ dirs[x] for x in indirs]
    set1=allfiles[0]
    set2=allfiles[1][0:5]+allfiles[2][0:20]+allfiles[1][5:]+allfiles[2][20:]
    set3=allfiles[3][0:5]+allfiles[4][0:10]+allfiles[5][0:10]+allfiles[3][5:]+allfiles[4][10:]+allfiles[5][10:]
    #print len(set1),len(set2),len(set3)
    
    sets=[set1,set2,set3]
    #print sets
    mergefilelist=zip(*sets)
    #print mergefilelist
    out=['%s/data/%s/alignments.bam'%(outdir,x.split('/')[-2]) for x in dirs[indirs[0]]]
    #print out
    mergedict=dict(zip(out,mergefilelist))
    #print mergedict
    mergebamfiles(pathsamtools,mergedict)    

if do_index=='combp':
    import common
    pvaluelist=[0.0001,0.85,0.8,0.85,0.05,0.0002,0.9,0.9]
    print common.combinepvalues(pvaluelist,method='fisher')
    print common.combinepvalues(pvaluelist,method='liptak')
    
if do_index=='confighelper':
    rootdir='/playpen/rootFDM'
    sources=['T4T%02dS01'%x for x in range(1,51)]
    sourcepath=['%s/dataout/%s/alignments.bam'%(rootdir,x) for x in sources]
    origsourcepath=['/playpen/sregen/TYPE4/data/T%02dS01/alignments.bam'%x for x in range(1,51)]
    
    mkdircmds=['mkdir -p %s'%'/'.join(x.split('/')[:-1]) for x in sourcepath]
    print '\n'.join(mkdircmds)
    
    lncmds=['ln -s %s %s'%(x,y) for x,y in zip(origsourcepath,sourcepath)]
    print '\n'.join(lncmds)
    
    project_groups_config_line=','.join(sources)
    data_config_lines='\n'.join(['%s=%s'%(x,y) for x,y in zip(sources,sourcepath) ])
    act_config_lines='\n'.join(['%s=0'%x for x in sources])
    
    print 'project_groups=',project_groups_config_line
    print data_config_lines
    print act_config_lines
    
if do_index=='testwgsdict':
    import common
    import cPickle
    wgsdict=cPickle.load(open('/playpen/rootFDM/dataout/DMx/DM_DM1.wgs'))
    genelist=wgsdict.keys()
    genelist=['AF085995']
    for gene in genelist[0:5]:
        wgs=wgsdict[gene]
        print gene
        print wgs
        sparsedict,nodedict=common.wgs2sparsegraphdict(wgs)
        print sparsedict
        print nodedict
        print

if do_index=='testactcorrect':
    import cPickle
    import actgraph
    import os
    import time
    
    ts='DM2'
    
    rpkfile='/playpen/rootFDM/dataout/%s/DMD_%s.rpk'%(ts,ts)
    actfile='/playpen/rootFDM/dataout/%s/DMD_%s.act'%(ts,ts)
    outactfile='/playpen/rootFDM/dataout/%s/DMDCORR_%s_subset3.act'%(ts,ts)
    fileptr=open(outactfile,'w')
    
    rpkdict=dict([lntxt.split('\t')[0],float(lntxt.split('\t')[4])] for lntxt in open(rpkfile))
    
    act0=actgraph.actFile(actfile)
    imagedir1='/playpen/rootFDM/dataout/%s/image'%ts
    imagedir2='/playpen/rootFDM/dataout/%s/imagecorr'%ts
    os.system('mkdir -p %s'%imagedir1)
    os.system('mkdir -p %s'%imagedir2)
    
    printall=0
    
    genelist=['FLJ00075','HES4']
    genelist=['CYP1A1','EGF','CIRBP','DAPK3','DHX15','DNMT1','EIF4G2','EREG','FAT2','CYR61','CXorf2','COX7A1','CPA1',
              'CPB1','CTPS','DCT','ECGF1','EYA3','EXOSC10','CLCN5','COL1A1','CYP1B1','CYP2D6','ERCC3','ERCC5','FAH',
              'F12','DAF','CYP1A2','CYP2B6','CYP2F1','CYP4A11','EIF2AK2','ENTPD1','CTNNB1','DSG3','FGB','CX3CL1',
              'EPO','ELAVL1','ERBB4','POLR2A','EZH2','F13A1','CRH','CHGB','EDN1','DEFB4','DGKG','DLG4','ERN1','CSNK2A2',
              'CSK','EPHA7','EIF2AK3','EPHA3','FGR','DCC','CHRNA4','FGF9','CHRNG','FANCG','CHST1','DUSP3','FKBP5','COX15',
              'CSPG3','DFFA','CRLF1','EGFR','E2F2','DAP','CSF1R','CP','CYP3A5','DYSF','CRYGB','FER','FGF2','DNAJC7','CYP2J2',
              'FADD','CRP','DRD5','CYP3A7','CYP2A13','CYP2C9','FBP2','ELK1','DPYD','FABP4','CUL4A','UBC','CKS2','18S','HPRT1',
              'ABCC8','POLR2A','CNIH','ABCB6','ARMET','CELSR1','ABCB1','ABCC6','ABCG2','CENPF','ABCC5','ABCC4','ALKBH',
              'CACNA2D2','CARD4','ARFGEF1','CACNG3','AP4B1','CD86','ADAT1','BRDG1','CLCA4','C11orf9','APH1A','C14orf122','CGI-119',
              'ASB4','ACTL6B','CLST11240','BET1L','BET1L','PHF21A','BAIAP2','C10orf92','BHMT2','C14orf103','C10orf3','CKAP2','AGPAT5',
              'ChGn','ABHD10','CACNA2D3','BTNL2','ABCC1','AGPAT4','C21orf7','CA10','CABC1','ANKRD2','C11orf15','ANGPTL7','BACH2',
              'ABCG5','CENTD3','APG3L','CDCP1','C20orf116','ASPSCR1','APG9L1','CARD14','ADIPOR2','ALG9','C14orf161','CLUAP1',
              'ACTR5','C9orf76','ASRGL1','C13orf23','ATF6','CDH19','C6orf47','CAP1','ABCC11','CGI-115','C2GNT3','BEX1','CCNL2',
              'ADAMTSL1','C15orf16','ABC1','C16orf33','BM88','AD023','C10orf93','C6orf210','SCCPDH','LMBRD1','C9orf121','C6orf78',
              'C3orf18','C12orf8','AICDA','UBC','CHCHD2','18S','HPRT1','DDX41','POLR2A','GNA13','DELGEF','CRISP3','C1orf61','GNAI3',
              'CYP46A1','GALNT6','EML2','FBXO22','DNAJB9','GLS2','EPN1','GGA1','DKFZP586A0522','GABRP','GAPDHS','DHDH','DLG7','DRD1IP',
              'DCXR','GNG13','FKBP11','DKFZp434H2215','FLJ20130','FAM83E','CASZ1','FLJ20449','EPN3','IARS2','FLJ10374','FLJ10560','DPPA4',
              'FLJ10922','GIMAP5','FOXJ2','FLJ11127','GK001','FSCN3','FGF22','DAB1','FLJ13149','FLJ13220','ELAVL4','GAN','EPB41L4A','ELSPBP1',
              'MOSC1','FLJ13855','GDAP1L1','CUEDC2','DCC1','GLB1L','FYCO1','FLJ12650','EFTUD1','UGT2A3','FLJ13848','CSPP','FLJ12443','FLJ13236',
              'FLJ21019','CCDC15','WDR71','ELL3','FOXB1','DKFZp762E1312','DUSP21','FLJ21901','FLJ14627','FLJ14816','DNAJC4','DCD','CRLF3','DIRAS2',
              'C1orf113','CHMP2B','DNAJB1','CYP3A4','FLJ20014','RNF186','FLJ23584','FBXO39','DEFB125','ESX1L','FLJ37549','FLJ35894','FBXL12',
              'DLEU1','FLJ25168','FXYD7','FLJ20244','UBC','18S','HPRT1','TP53','SLC25A4','TNFRSF9','S100A5','SCN9A','SCP2','SLC5A2','SPIB',
              'CLEC3B','UBE2B','RHAG','SLC12A3','TGFBI','TGM1','TYR','SLC10A2','SOD1','SLC2A4','SLC12A2','TFF1','THBS1','RECQL4','POLR2A',
              'RAD52','TRIM21','RORC','TNFRSF8','TNFSF8','TFRC','SCAP1','SRPK1','TK1','TJP2','TIE1','SHH','TIAM1','TPD52','SRF','TCOF1',
              'TAP1','SERPINB5','USP5','SERPINB7','SNX4','TNFSF14','SNAP23','SGCE','SDHA','SNRPD3','TRIP13','TRFP','TRAF4','SEPT2','SFRS10',
              'SPARCL1','SYNGR1','SFRS11','STOML1','SERPINI1','RPL3L','SFPQ','SLC1A6','SLC17A1','SCO2','THRB','STAT4','TFAP2C','RUNX1T1','TP73',
              'SERPINB2','STAT5A','UGCG','TPST2','SLC9A5','SCGB2A2','SLC22A4','THRSP','TJP1','TNNC2','TNNI2','ROM1','TNP1','SLC16A3','S100A8','SUCLG1',
              'TRA1','TGFBR2','STAT6','SRD5A1','TRPC1','TDG','SFRS3','TMSB4Y','UBC','18S','HPRT1','MYC','RAD51','MDM4','MMP7','NEUROD1','NVL',
              'PEX13','PGM1','PLS1','PMP2','POLB','POLE2','PPP1CC','PRPSAP1','PSD','PSMA5','PSMB10','PTPRB','PZP','MYO5A','PEX1','NQO2','P4HB',
              'PNLIP','PRLR','HTRA1','PPFIA2','POLR2A','MCM3','MGMT','MCM3AP','MLNR','PPP3CA','NFKB2','POMC','PGLYRP1','MVK','NTRK3','PCTK2','PIK3C3','PRKCG',
              'MAPKAPK3','MPP2','PLA2G2A','PTCH','PPIC','MX1','NGFR','MTRF1','PEX14','MGAM','PIGB','MRPL33','NDUFAB1','NDUFS1','PDE1C','MYB','NEF3','POU1F1',
              'NR2C2','MDS1','PPARA','NR1I3','MMP1','MMP13','MMP2','MMP9','MIF','PAX3','MDM2','PER1','PPP6C','NME1','MAGEB3','MYH8','OMG','P2RY11','PITX1','PKP4',
              'NUMB','NOL4','PABPC4','NAP1L3','POLR2G','PSAP','RAB11A','OASL','PCNA','PSMD7','RAC2','PFDN1','NAT2','P2RY1','UBC','18S','HPRT1','CCNA2','CFLAR','ACTN2',
              'ADPRH','ALOX12B','ANXA5','AOX1','CD68','AHSG','AMBP','APOF','ARF5','ATP1B2','BTN1A1','CD58','CDA','CENPA','ACO1','ABCD3','ACADM','ACADS','ADA','APOH',
              'ASPA','BTD','ABCC2','AMT','ANXA1','CACNG1','ANGPT2','CDH1','CCR7','CCL20','ARHGDIB','BLM','POLR2A','CD3D','CD81','ADCYAP1','CD80','ATM','CDK4','CDK8',
              'CDK9','BRS3','APC','APOB','ADCY8','ADM','BARD1','ATP6AP1','ABCB11','BAG1','CDC45L','ALDH4A1','BECN1','B2M','BCL2A1','ANXA13','BPHL','CAD','CALR','BRPF1','CALB1','CCS','CDX4','CASP5','CEACAM5',
              'CDK5R1','ATOH1','ADRB1','ATP1A2','ATP6V1B1','CCR3','CCNB2','CDKN1A','BIRC2','C16orf3','ALPI','C11orf13','CDK5','C4BPB','BAG4','BAT1','ADORA2B',
              'ACAT1','BYSL','BRCA2','ALOX15','BMP3','AP1S1','AMD1','ACADVL','UBC','18S','HPRT1','PLAU','PAX5','POLR2A','RBM5','NRL','RRH','NTSR2','NTS','PTK2','PKIA','NRGN',
              'PDGFRA','OPTN','S100A12','PPIF','PLXNC1','S100P','NELL1','ORC2L','PPP2R5B','RPP30','NPC2','PSMC4','RPIP8','PDIA5','RAB32','NPM3','PWP1','PICALM','POLI','MMD',
              'PADI4','PSD4','MRPL13','ORC6L','RAB30','RB1CC1','RNF10','POLR1D','PADI3','PRRX2','DTL','MASK','NUDT11','NEIL3','RSAD1','C1orf103','PBK','POLR1B','NS3TP1','PELI1',
              'NDUFV3','NGB','MRPL46','RAB38','OSBPL11','MOSPD3','MLPH','PSTPIP2','PRKRIP1','NIP30','RUFY1','PUS3','POU6F2','PPARG','RRAGB','NCK1',
              'REV1L','MS4A4A','OR5V1','RNASE6','P2RY5','PAICS','OR1F1','PURG','PNMA3','RAPGEFL1','RLN3R1','SALL4',
              'NEDD8','PES1','MTCBP-1','NKG7','RAD21','RNUT1','RNASE6','PDIA2','NURIT','PPARD','MPP4','RPP38','NONO','UBC','PREI3','18S','HPRT1','HIF1A','GBAS','INPPL1','GFAP',
              'GFPT1','GLP1R','GYS1','HSD17B2','ITIH2','ITIH3','KCNJ3','KCNS3','KIF3C','LIFR','LOXL2','LUM','FURIN','HGD','ICAM1','G6PD','IL4R','IGFBP2','GABRA1','GRM1',
              'GRM3','GSTP1','IL1R1','ITGB7','JAK3','HOXD13','IGF2','INSL4','POLR2A','GRLF1','LPL','GP1BA','KIT','IL6','IL13','IL18R1','KDR','IGFBP5','GABRA6','GABRD','GRIA2',
              'IGF2R','HBEGF','GSTT1','LGI1','KCNK5','KYNU','HABP2','HADH2','KCNC3','GOLGA5','GIF','GRPR','GZMM','HIP2','FLI1','HCFC1','ICAM3','JAK1',
              'LRP1','JAK2','LTA','HS3ST1','GAS1','KCNA6','GPR65','HOXB7','KCNQ2','GPR8','GPR22','HIST1H1D','HSPA6','JUN','KNG1','JUNB','ITIH4','HRG','IGFBP1','KCNC1',
              'HIST1H3A','FSCN1','HIST3H3','HDAC1','GABRA4','GGH','G6PC','GLI3','KALRN','HBZ','UBC','18S',
              'HPRT1','LMNA','GULP1','HUNK','POLR2A','LOXL1','GPR64','IRS1','MAP3K5','LRRC17','MLLT4','IQGAP2',
              'ITGB1BP2','INPP5A','CD180','KLF10','KIF20A','INADL','KNTC2','GZMA','IRF6','MASP2','HEAB','KIF3A','GUCA2B','MYLPF','HSU79274','LR8','IL1F6','JMJD2A','KIAA0773',
              'KIAA0152','KIAA0101','KIAA0652','MELK','KIAA0196','VASH1','CLCC1','LEPROTL1','HS747E2A','LOC51066','LOC51337','LOC55831','MDS032','HCA112','IL17RB','KRT24','IL1F9',
              'LXN','MEPE','LZTFL1','HB-1','LGR7','ITM2B','TINAGL1','LHX5','MFSD1','March7','HTR5A','MGC4504','MAPKAP1','LENG4','KLHL18','IQCG','HEY2','LHX6','C17orf40','IARS','HOXB8',
              'LOC81691','MEGF11','KCNE1L','BTBD15','TAAR2','HYDIN','C16orf24','LRRC7','MGC2963','MIA2','MAB21L1','KCTD5','HYDIN',
              'HOXA5','KIAA1900','HTRA4','LOC132321','HOXA7','HDHD1A','KIF4A','LIN28','MGC17337','KCNE4','UBC','LDHA','LYPLA2','18S','HPRT1','SPDEF','POLR2A','SEMA3A','WIF1','TRAF2',
              'SLC16A2','TRAF1','SLC17A2','TSPAN31','SIX1','TGOLN2','SLC2A1','SLC22A7','SLC38A3','SOX15','STMN2','VIL1','TU3A','SLC6A14','SLCO2B1','TUSC2','SCRG1','TAGLN3','ST8SIA5','SLC16A8','TMEM5','TTLL3','TTC15','VPS24',
              'TRIM17','TRPC4','VRK3','SPTBN5','UBN1','SPOCK3','SURF2','TSC','SLC30A6','TMEM19','SEC5L1','UEV3','TRIM36','XAB2','SLC15A2','SEZ6L','TNFAIP1','TRIB2','UCP1','TEX27','SRR','SLC13A1','TPK1','SLC30A5','SLC25A23',
              'SHCBP1','TREML2','TMC7','WDR19','SPTBN4','WNT4','STMN4','TCF7L1','TFE3','TBR1','TRPS1','TCF8','TNFAIP3','STAT1','SEPT9','SLC4A1AP','SPATA20','WNT2B','TOE1','SYNC1','ST3GAL4','TLE4',
              'SUGT1','USP28','TMSB10','TEX13B','TBN','SHC3','CD300LG','STAG3','THUMPD1','TAIP-2','TFG','ST18','SMAP','SEMA4G','TXNDC10','TMSL8','STARD8','UBC','18S','HPRT1','RB1','BIRC5','CHAD','CHGA','ABCD1','CD3E','VWF','ZNF262',
              'ZNF258','WRN','POLR2A','ZNF146','AFP','IL1B','VCAM1','FLT3','CGA','CHEK1','MAPK14','TGFA','MAP3K14','MAP3K8','WNT1','RARRES3','XK','ZP2','FPGT','ABCD2','TROAP','ABCA1','DMXL1','MLANA','HMGCS2','COG2','COMMD3',
              'GTF3C5','ZNF225','SEC24D','CRELD1','ZBTB20','GTF2A1','ZNF331','ZNF167','MOSPD1','ZNF350','APOBEC3G','NR0B2','ZNF557','SE57-1','COASY','RARA','ZNF140','ZNF192','ABCA2','P2RXL1',
              'FOXA1','ZFP36L2','ZNF358','FMO5','ZNF101','ZIM3','MMP21','ZSCAN4','ZNF434','ZNF548','WNT2','VAMP3','ZMYM1','UBC','18S','HPRT1','SELE','LDLR','IL8','PTGS2','TYMS','DAD1','FABP1','DPP4','MX2','TNFRSF11B','DNCL1','MAP2K6','AES']
    genelist=['DDRGK1','C12orf24','FAM105A','CCDC79','GUF1','ZEB1','NCAN','HEATR6','HSP90B1','FNDC8','C22orf31','CEP55','CCDC94','UBE2Z','PAAF1','ATG2B',
              'NOD1','CCDC47','CSRNP3','SKAP1','TBCCD1','NAT11','TMEM9B','C8orf55','DLGAP5','CALY','ASNSD1','MBOAT7','ZMYM4','STAP1','USE1','ZFAND3','HJURP',
              'ARAP3','FAM192A','ALKBH1','TAF8','MLEC','NDC80','RUNDC3A','DNAJC22','TMEM176B','MARCH7','CSGALNACT1','TMBIM4','SYNC','MIF4GD','HIGD1B','CCDC68',
              'FASTKD5','OTUD7A','FAM107A','HSD17B10','NUP62CL','TRMT1','SPERT','DSCC1','HMHB1','SNRNP25','NEFM','C4orf33','NKAIN1','RASSF7','RMI1','RRP15',
              'CLP1','CSPP1','REV1','FAM131B','MANF','CHAC1','TMEM176A','FAM173A','CATSPERB','TESC','MST4','ESX1','SNUPN','TMEM93','SOHLH2','LPCAT1','FAM158A',
              'ATG3','PDSS2','EXOC2','ZBTB44']
    genelist=['P2RX6', 'HMFT1638', 'CICE', 'fls485', 'METTL7A', 'HCA66', 'DelGEF', 'LIG', 'TMX3', 'DKFZp434G1310', 'FASTKD1', 'KDM4A', 'C9orf30', 'PTCH1']
    #genelist=[]
    allislandlist=cPickle.load(open('/playpen/rootFDM/annotation/fdm_human_hg19_knowngene.gtf/islandlist.pck'))
    allislandlist.sort()
    if len(genelist)>0:
        islandlist=[island for island in allislandlist if island[3] in genelist]
        print len(genelist), len(islandlist)
    else:
        islandlist=allislandlist
    islandlist=islandlist
    lncnt=0
    for island in islandlist:
        lncnt+=1
        gene=island[3]
        time1=time.time()
        wgrstuple,error=act0.ACT2Corrected(gene,5)
        timetaken=time.time()-time1
        meanerr=error/(rpkdict[gene]+0.0001)
        print '###\t%d/%d\t%s\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%d'%(lncnt,len(islandlist),gene,error,rpkdict[gene],meanerr,timetaken, 
                                                        len(wgrstuple[0])+len(wgrstuple[1])+len(wgrstuple[2])+len(wgrstuple[3])+len(wgrstuple[4]))
        if lncnt%200==1 or (rpkdict[gene]>10 and meanerr<0.1) or printall==1:  
            #print 'printing...',gene 
            act0.Toimage(gene, imagedir1) 
            act0.Toimage(gene, imagedir2,inwgs=wgrstuple)  
        try:
            actgraph.wgrstuple2file(wgrstuple,island,fileptr)
        except:
            print 'Problem in writing act',[len(x) for x in wgrstuple]
    fileptr.close()

if do_index=='testdom':
    import project
    
    flow_file='/playpen/rootFDM/project/TYPE301/flows/TYPE301_ALL_flows.txt'
    adict=project.extractNdominantsplices(flow_file,N=2)
    k=adict.keys()
    print adict[k[0]]

if do_index=='combp2':
    import common
    pvaluelist=[0.00,0.031200,0.000000,0.0,0.393500,0.000000,0.002200,0.00]
    print common.combinepvalues(pvaluelist,method='fisher')
    print common.combinepvalues(pvaluelist,method='liptak')

if do_index=='testpval':    
    import project
    distrlist=[0,0,0,0,0,0,0,0,0,0,0]
    x=0
    print project.distr2pvalue(distrlist,x)
    import common
    csvfile='/playpen/nextgen/exam/OralExamList.txt'
    htmlfile=open('/playpen/nextgen/exam/OralExamList.html','w')
    common.csv2html(open(csvfile),htmlfile,',')
    htmlfile.close()

if do_index=='paul_dom':
    import project
    minwt=5.0
    bsiggenefile='/playpen/rootFDM/project/STEM2/report/BSIG_gene.txt.csv'
    bsiggeneposlist=[(lntxt.split(',')[0],int(lntxt.split(',')[2])) for lntxt in open(bsiggenefile).readlines()[1:]]    
    print len(bsiggeneposlist)
    dominanttslist=['STEM-DISC-L','STEM-DOWC-L','STEM-DURM-L','STEM-GRAA-L','STEM-HENE-L','STEM-KERR-L','STEM-SHOL-L','STEM-VARV-L']
    flowfile='/playpen/rootFDM/project/STEM2/flows/T3LSC_ALL_flows.txt'
    domdict=project.extractNdominantsplices(flowfile,dominanttslist,N=2)
    domkeys=domdict.keys()
    edges,datadict=domdict[domkeys[0]]
    tslist=datadict.keys()
    tslist.sort()

    outfile1='/playpen/rootFDM/project/STEM2/report/BSIG_main_flows.txt'
    outfile2='/playpen/rootFDM/project/STEM2/report/BSIG_main_flows_wt_annotated.txt'
    fout1=open(outfile1,'w')
    fout2=open(outfile2,'w')
    header1=['gene','pos','edge1','edge2']
    header2=['gene','pos','edge1','edge2']
    h1cols=['edge1_frac','edge2_frac','edge_tot_wt']
    h2cols=['edge1_frac','edge2_frac']
    for ts in tslist:
        for h1col in h1cols:
            header1.append('%s__%s'%(ts,h1col))
    for ts in tslist:
        for h2col in h2cols:
            header2.append('%s__%s'%(ts,h2col))
    fout1.write('%s\n'%'\t'.join(header1))
    fout2.write('%s\n'%'\t'.join(header2))
    for key in bsiggeneposlist:
        if key not in domdict:
            print key
        else:
            edges,datadict=domdict[key]
        if len(edges)==2:
            line1=[key[0],'%d'%key[1],'%s'%str(edges[0]),'%s'%str(edges[1])]
            line2=[key[0],'%d'%key[1],'%s'%str(edges[0]),'%s'%str(edges[1])]
        else:
            line1=[key[0],'%d'%key[1],'%s'%str(edges[0]),'']
            line2=[key[0],'%d'%key[1],'%s'%str(edges[0]),'']            
        for ts in tslist:
            line1+=['%6.4f'%datadict[ts][0][0],'%6.4f'%datadict[ts][0][1],'%7.2f'%datadict[ts][1]]
            if datadict[ts][1]>minwt:
                line2+=['%6.4f'%datadict[ts][0][0],'%6.4f'%datadict[ts][0][1]]
            else:
                line2+=['_%6.4f_'%datadict[ts][0][0],'_%6.4f_'%datadict[ts][0][1]]
        fout1.write('%s\n'%'\t'.join(line1))
        fout2.write('%s\n'%'\t'.join(line2))
    fout1.close()
    fout2.close()

if do_index=='paul_novelcntplot':
    def splicecount(filename,type,mincnt):
        distrdict={}
        genedict={}
        for lntxt in open(filename):
            ln=lntxt.rstrip('\n').split('\t')
            if ln[2]!='splice':
                continue
            if ln[1]==type and float(ln[5])>=mincnt:
                genedict[ln[8]]=genedict.get(ln[8],0)+1
        for gene in genedict:
            distrdict[genedict[gene]]=distrdict.get(genedict[gene],0)+1
        return distrdict
    
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pylab as plt
    mincnt=3.0
    totgene=28092
    tsprefixlist=['STEM-DISC','STEM-DOWC','STEM-DURM','STEM-GRAA','STEM-HENE','STEM-KERR','STEM-SHOL','STEM-VARV']
    for tsprefix in tsprefixlist:
        for suffix in ['L','N']:
            filename='/nas02/home/d/a/darshan/prinsrootFDM/dataout/%s-%s/T3LSC_%s-%s.act'%(tsprefix,suffix,tsprefix,suffix)
            imagefilename='/nas02/home/d/a/darshan/prinsrootFDM/project/T3LSC/scripts/%s_%s.pdf'%(tsprefix,suffix)
            fig = plt.figure(figsize=(8, 8))
            plt.title('%s_%s'%(tsprefix,suffix))
            plt.ylabel('Number of Genes')
            plt.xlabel('Count of Splices')
            distrdict1=splicecount(filename,'annot',mincnt)
            distrdict2=splicecount(filename,'novel',mincnt)
            xlist=range(1,max(distrdict1.values()+distrdict2.values())+1)
            #y1list=[distrdict1.get(x,0) for x in xlist[1:]]
            y1list.insert(0,totgene-sum(y1list))
            #y2list=[distrdict2.get(x,0) for x in xlist[1:]]
            y2list.insert(0,totgene-sum(y2list))
            plt.plot(xlist,y1list,'k.')
            plt.plot(xlist,y2list,'c.')
            plt.savefig(imagefilename)
            print '%s_%s'%(tsprefix,suffix),'annot',sum([x*y for x,y in zip(xlist,y1list)])
            print '%s_%s'%(tsprefix,suffix),'novel',sum([x*y for x,y in zip(xlist,y2list)])



def davidarrowplot(title,xdata,y1y2data,slope1,slope2,axes,imagefilename):
    fig = plt.figure(figsize=(16, 16))
    ax = fig.add_subplot(111)
    for i in range(len(y1y2data)):
        plt.plot(xdata[i],y1y2data[i][0],'ro')
        plt.plot(xdata[i],y1y2data[i][1],'bs')
        expy=slope1*xdata[i]
        if abs(expy-y1y2data[i][0])<abs(expy-y1y2data[i][1]):
            improveflg=0
        else:
            improveflg=1
        if improveflg==0:
            plt.annotate("",
                        xy=(xdata[i],y1y2data[i][0]), xycoords='data',
                        xytext=(xdata[i],y1y2data[i][1]), textcoords='data',
                        arrowprops=dict(arrowstyle="<|-",linewidth=0.5,connectionstyle="arc3",ec='Orange'),
                        ) 
        if improveflg==1:
            plt.annotate("",
                        xy=(xdata[i],y1y2data[i][0]), xycoords='data',
                        xytext=(xdata[i],y1y2data[i][1]), textcoords='data',
                        arrowprops=dict(arrowstyle="<|-",linewidth=0.5,connectionstyle="arc3",ec='Lime'),
                        )             
    plt.title(title)
    plt.plot([0,100],[0,slope1*100],'r-')
    plt.plot([0,100],[0,slope2*100],'b-')
    plt.xlabel('MAQC')
    plt.ylabel('Splice Coverage')
    plt.axis(axes)
    try:
        plt.savefig(imagefilename)
    except:
        message='Problem in writing %s'%imagefilename
        common.printstatus(message,'W',common.func_name())            

if do_index=='davidmaqdata':
    import common
    filelist=['/playpen/rootFDM/project/DMD/scripts/junction_maqc_values_%d.txt'%i for i in range(1,5)]
    outfileptr=open('/playpen/rootFDM/project/DMD/scripts/maqc_values_all.txt','w')
    common.mergefilelist(filelist,[5],outfileptr,keycollist=[0,1,2,3,4],sep='\t')
    
if do_index=='davidACTdata':
    act1='/playpen/rootFDM/dataout/DM1/DMD_DM1.act'
    actc1='/playpen/rootFDM/dataout/DM1/DMDCORR_DM1_subset.act'
    act2='/playpen/rootFDM/dataout/DM2/DMD_DM2.act'
    actc2='/playpen/rootFDM/dataout/DM2/DMDCORR_DM2_subset.act'
    
    actfilelist=[act1,actc1,act2,actc2]
    
    maqcfile=open('/playpen/rootFDM/project/DMD/scripts/maqc_values_all.txt')
    probelist=[]
    probedict={}
    for lntxt in maqcfile:
        ln=lntxt.rstrip('\n').rstrip('\r').split('\t')
        probelist.append((ln[2],int(ln[3]),int(ln[4])+1))
        probedict[(ln[2],int(ln[3]),int(ln[4])+1)]=ln[5:9]
    probelist.sort()
    
    
    
    
    allactdict={}
    for act in actfilelist:
        actdict={}
        for lntxt in open(act):
            ln=lntxt.rstrip('\n').split('\t')
            if ln[2]!='splice':
                continue
            if ln[8] in ['CR592260', 'BC014606', 'gis5', 'DKFZp313O1925','TGN','DKFZp686J07132','ScRG-1','FUS1','NR1C1'
                         , 'DL491360', 'N-Tes','BCYRN1', 'PC-1', 'KIAA0755', 'Rh50', 'HOYS6', 'Elk1', 'MTABC', 'TMEM143', 'MRP6',
                          'hTPK1','PSF','MDMX','C5orf32', 'AX748200','GM2A','DKFZp564P0662','AK092143','LGMD2B','BC036909','BC062349',
                          'BC019904','TIE','DKFZp762K197','CDD','NF1','UTP6', 'NME1-NME2', 'KIAA1515', 'CD55', 'patched a isoform', 'AX747756', 
                          'RIE2 sid2705', 'TMEFF1', 'SERGEF', 'TS', 'ALS2CR10', 'C19orf23', 'ARC20', 'RB3', 'PIF', 'KIAA1638', 'NR1A2', 'TMEM111', 
                          'KIAA1830', 'Pex14', 'ENOSF1', 'ADI1', 'SUGT1P', 'INS-IGF2', 'PIG33', 'C3orf32', 'RXFP1', 'MRP5', 'PLEKHH3', 
                          'DKFZp586A0522', 'DKFZp761D0712', 'GFPP', 'ARHGDIG', 'CERKL', 'hPer', 'MED20', 'SAS', 'KIAA0795', 'HP2XM', 
                          'DKFZp564O123', 'TRA2B', 'MAP3K7IP1', 'DKFZp761J1810', 'KLHL32', 'DPE2', 'MPDU1', 'PPFIA2 variant protein', 
                          'FLJ00383', 'KIAA0535', 'BMI1', 'MELKv3', 'MELKv2', 'BAF53b', 'pp9974', 'ATB0+', 'AK127087', 'PKCG', 'KIAA0677', 
                          'hINADL', 'smap-4', 'DKFZp666O0110', 'TNNI3K', 'SYNGR1C', 'UBE2K', 'IDH3GL', 'NME2', 'AK096625', 'hAO', 'LPAAT-delta', 
                          'DM119463', 'CR612522', 'EFA6B', 'ZNF547', 'DKFZp434B103', 'PleA 1', 'KIAA1800', 'UNC24', 'DM004282','CICE']:
                continue
            splice=(ln[0],int(ln[3]),int(ln[4]))
            if splice in probelist:
                if splice not in allactdict:
                    allactdict[splice]=[[],[]]
                if allactdict[splice][0]!=[ln[8]]:
                    allactdict[splice][0]+=[ln[8]]
                allactdict[splice][1]+=[ln[5]]
    
    moregenelist=[];iggenelist=[]
    for splice in probedict:
        if splice in allactdict:
            if len(allactdict[splice][0])!=1:
                if len(allactdict[splice][0])==4:
                    print '####', allactdict[splice]
                    moregenelist.append(allactdict[splice][0][0])
                    iggenelist.append(allactdict[splice][0][1])
                else:
                    xlist=allactdict[splice][0]
                    for x in xlist:
                        if xlist.count(x)!=4:
                            iggenelist.append(x)
                    print '#####',splice,len(allactdict[splice][0]),allactdict[splice] 
#            
#                print splice, probedict[splice],allactdict[splice][0],allactdict[splice][1]
#            else:
                
        
    print moregenelist
    print list(set(iggenelist))
    
    print len(allactdict),len(probedict)
    
    actmqqcdata=open('/playpen/rootFDM/project/DMD/scripts/maqc_act_data.txt','w')
    dataline=['chr','splice1','splice2','gene','maqc1','maqc2','maqc3','maqc4','DM1','DMCORR1','DM2','DMCORR2']
    actmqqcdata.write('%s\n'%'\t'.join(dataline))
    for splice in allactdict:
        outline=[splice[0],'%d'%splice[1],'%d'%splice[2]]+allactdict[splice][0]+probedict[splice]+allactdict[splice][1]
        actmqqcdata.write('%s\n'%'\t'.join(outline))
    actmqqcdata.close()
            
if do_index=='davidplot':
    import matplotlib.pylab as plt
    import numpy as np
    
    infile='/playpen/rootFDM/project/DMD/scripts/maqc_act_data.txt'
    mincov=5
    
    
    
    

    lnlist=[]
    lncnt=0
    for lntxt in open(infile):
        lncnt+=1
        if lncnt==1:
            continue
        ln=lntxt.rstrip('\n').split('\t')
        lnlist.append(ln)
    
    datacollist=[[4,8,9],[4,10,11],[5,8,9],[5,10,11],[6,8,9],[6,10,11],[7,8,9],[7,10,11],[0,0,0]]
    #datacollist=[[0,0,0]]
    
    namedict={4:'MAQC1',5:'MAQC2',6:'MAQC3',7:'MAQC4',8:'DM1',10:'DM2',0:'AVG'}
        
    for datacols in datacollist:
        xdata=[];y1data=[];y2data=[]
        filxdata=[];fily1data=[];fily2data=[]
        print datacols
        for ln in lnlist:
            if len(ln)!=12:
                #print ln
                continue
            if datacols==[0,0,0]:
                xdata.append(0.25*(float(ln[4])+float(ln[5])+float(ln[6])+float(ln[7])))
                
                y1data.append(0.5*(float(ln[8])+float(ln[10])))
                y2data.append(0.5*(float(ln[9])+float(ln[11])))
            else:
                xdata.append(float(ln[datacols[0]]))
                
                y1data.append(float(ln[datacols[1]]))
                y2data.append(float(ln[datacols[2]]))
    
            if y1data[-1]>mincov:
                filxdata.append(xdata[-1])
                fily1data.append(y1data[-1])
                fily2data.append(y2data[-1])
        print len(xdata),len(filxdata)
            
        imagefilename1='/playpen/rootFDM/project/DMD/scripts/image/MAQCplot_%s_%s_60.pdf'%(namedict[datacols[0]],namedict[datacols[1]])
        imagefilename2='/playpen/rootFDM/project/DMD/scripts/image/MAQCplot_%s_%s_06.pdf'%(namedict[datacols[0]],namedict[datacols[1]])
        imagefilename3='/playpen/rootFDM/project/DMD/scripts/image/MAQCplot_%s_%s_01.pdf'%(namedict[datacols[0]],namedict[datacols[1]])
        
        x=np.array(xdata)
        y=np.array(y1data)
        print 'uncorrected_all corr',np.corrcoef(x,y)
        x = x[:,np.newaxis]
        a, resid= np.linalg.lstsq(x, y)[:2]
        rsq1 = 1 - resid / (y.size * y.var())
        print a,rsq1
        slope1=a[0]
        
        x=np.array(xdata)
        y=np.array(y2data)
        print 'corrected_all corr',np.corrcoef(x,y)
        x = x[:,np.newaxis]
        a, resid= np.linalg.lstsq(x, y)[:2]
        rsq2 = 1 - resid / (y.size * y.var())
        print a,rsq2
        slope2=a[0]
    
        x=np.array(filxdata)
        y=np.array(fily1data)
        print 'uncorrected_fil corr',np.corrcoef(x,y)
        x = x[:,np.newaxis]
        a, resid= np.linalg.lstsq(x, y)[:2]
        rsq = 1 - resid / (y.size * y.var())
        print a,rsq
        
        x=np.array(filxdata)
        y=np.array(fily2data)
        print 'corrected_fil corr',np.corrcoef(x,y)
        x = x[:,np.newaxis]
        a, resid= np.linalg.lstsq(x, y)[:2]
        rsq = 1 - resid / (y.size * y.var())
        print a,rsq
        
        title='MAQC_comparison\nR2 raw=%5.3f\nR2 corr=%5.3f'%(rsq1,rsq2)
        
        y1y2data=zip(y1data,y2data)
#        davidarrowplot(title,xdata,y1y2data,slope1,slope2,[0,60,0,1000],imagefilename1)
#        davidarrowplot(title,xdata,y1y2data,slope1,slope2,[0,6,0,100],imagefilename2)
#        davidarrowplot(title,xdata,y1y2data,slope1,slope2,[0,1,0,40],imagefilename3)

if do_index=='thesis_pdfs':
    import os
    infile='/tmp/images.txt'
    outdir='/home/darshan/workspace/athesis/mainmatter/figures'
    for f in open(infile):
        if f=='\n':
            continue
        f=f.rstrip('\n')
        cmd='echo %s > %s/%s.txt'%(f,outdir,f)
        print cmd
        os.system(cmd)
        cmd='enscript -f Times-Roman80 -B  -r %s/%s.txt -o - | ps2pdf - %s/%s.pdf'%(outdir,f,outdir,f)
        print cmd
        os.system(cmd)
        cmd='rm %s/%s.txt'%(outdir,f)
        print cmd
        os.system(cmd)

if do_index=='T1JS':     
    import common   
    infile='/playpen/nextgen/project/T1JS/report/unfiltered/top_report.txt'
    csvfile='/playpen/nextgen/project/T1JS/report/top_report.txt'
    fout=open(csvfile,'w')
    bdict={'SIG':1,'NOT':0}
    wdict={'SIG':-1,'NOT':0}
    lncnt=0
    for lntxt in open(infile):
        ln=lntxt.rstrip('\n').split('\t')
        if ln[0]=='#':
            fout.write('%s\n'%'\t'.join(ln))
            continue
        diffcnt=bdict[ln[21]]+bdict[ln[24]]+bdict[ln[27]]+bdict[ln[30]]+wdict[ln[18]]+wdict[ln[33]]
        if diffcnt<4:
            print ln[18], ln[33], ln[21],ln[24],ln[27],ln[30]
            continue
        lncnt+=1
        ln[0]='%d'%lncnt
        fout.write('%s\n'%'\t'.join(ln))
    fout.close()
    htmlfile=open('/playpen/nextgen/project/T1JS/report/top_report.html','w')
    common.csv2html(open(csvfile),htmlfile)
    htmlfile.close()

if do_index=='nothing':    
    x=[i for i in range(20,94) if i%3!=1]


if do_index=='KATE1':
    import actgraph
    import ast
    import os 
    
    genelist=['NCBP2','TNRC6A','HIF1A']
    tslist=['JSCTL1','JSCTL2','JSKO01','JSKO02']
    project_name='T1JS'
    rootdir='/playpen/rootFDM'
    actfilelist=['%s/dataout/%s/%s_%s.act'%(rootdir,ts,project_name,ts) for ts in tslist]
    imagedir='/playpen/rootFDM/project/%s/report/image'%project_name
    for actfilename in actfilelist:
        act0=actgraph.actFile(actfilename)
        for gene in genelist:
            act0.Toimage(gene, imagedir,ext='pdf')

if do_index=='KATEbox1':
    import common
    import ast
    fast_min_cov=5.0
    root_dir='/playpen/rootFDM'
    project_name='KATE'
    flowfile='%s/project/%s/flows/%s_ALL_flows.txt'%(root_dir,project_name,project_name)
    FEFfile1='%s/project/%s/flows/%s_ALL_FEFs.txt'%(root_dir,project_name,project_name)
    FEFfile2='%s/project/%s/flows/%s_ALL_FEF_processed.txt'%(root_dir,project_name,project_name)
    fout1=open(FEFfile1,'w')
    fout2=open(FEFfile2,'w')
    lncnt=0
    for lntxt in open(flowfile):
        lncnt+=1
        ln=lntxt.rstrip('\n').split('\t')
        if lncnt==1:
            hdict=common.findheaderdict(ln)
            tslist=hdict.keys()
            tslist.sort()
            fout1.write('\t'.join(ln[0:7]))
            fout1.write('\tedge used')
            for ts in tslist:
                fout1.write('\t%s:FEF\t%s:wt'%(ts,ts))
            fout1.write('\n')
            fout2.write('\t'.join(ln[0:7]))
            fout2.write('\tedge used')
            for ts in tslist:
                fout2.write('\t%s:FEF'%(ts))
            fout2.write('\n')                
        else:
            edges=ast.literal_eval(ln[6])
            numrows=max(1,len(edges)-1)
            for i in range(numrows):
                fout1.write('\t'.join(ln[0:7]))
                fout1.write('\t%s'%str(edges[i]))
                fout2.write('\t'.join(ln[0:7]))
                fout2.write('\t%s'%str(edges[i]))                    
                for ts in tslist:
                    wt1=float(ln[hdict[ts][0]])
                    wt2=float(ln[hdict[ts][1]])
                    flow=ast.literal_eval(ln[hdict[ts][2]])
                    FEF=flow[i]
                    fout1.write('\t%6.4f\t%6.4f'%(FEF,wt2))
                    if wt2<fast_min_cov:
                        FEF=-1
                    fout2.write('\t%6.4f'%(FEF))
                fout1.write('\n')
                fout2.write('\n')         
    fout1.close()
    fout2.close()

if do_index=='KATEbox2':
    import common
    import ast
    import scipy.stats
    import numpy
    import random
    import plot
    project_groups=[['JSREN%03d'%i for i in range(1,39)],['JSREN%03d'%i for i in range(39,417)]]
    root_dir='/playpen/rootFDM'
    project_name='KATE'    
    outlinelist=[]
    mingroupsz=10
    numiter=101
    minstd=-0.01
    FEFfile2='%s/project/%s/flows/%s_ALL_FEF_processed.txt'%(root_dir,project_name,project_name)
    twogroupdiff='%s/project/%s/report/%s_Two_Group_diff.txt'%(root_dir,project_name,project_name)    
    fout=open(twogroupdiff,'w')
    lncnt=0
    for lntxt in open(FEFfile2):
        lncnt+=1
        ln=lntxt.rstrip('\n').split('\t')
        if lncnt==1:
            hdict=common.findheaderdict(ln,coltype=['FEF'])
            tslist=hdict.keys()
            tslist.sort()
            fout.write('\t'.join(ln[0:8]))
            fout.write('Size Group1\tSize Group2\tGroup1 Median\tGroup2 Median\tRandom Median1\tRandom Media2\tSample Variance')
            fout.write('\tSample Range\tKruskal pvalue\tMedian Random pvalue\tGroup1 Median ACT\tGroup2 Median ACT\n')
        else:
            allgroup1=[float(ln[hdict[ts][0]]) for ts in tslist if ts in project_groups[0]]
            allgroup2=[float(ln[hdict[ts][0]]) for ts in tslist if ts in project_groups[1]]
            group1=[x for x in allgroup1 if x>=0.0]
            group2=[x for x in allgroup2 if x>=0.0]
            group1.sort()
            group2.sort()
            if len(group1)<mingroupsz or len(group2)<mingroupsz:
                continue 
            medgroup1=common.median(group1)
            medgroup2=common.median(group2)
            medts1=[ts for ts in tslist if ts in project_groups[0] and float(ln[hdict[ts][0]])-medgroup1<0.001 ][0]
            medts2=[ts for ts in tslist if ts in project_groups[1] and float(ln[hdict[ts][0]])-medgroup2<0.001 ][0]
            bothgroup=group1+group2
            samplestd=scipy.stats.nanstd(bothgroup)
            #print lncnt, samplestd
            if samplestd<0.0001:
                adddummy=max(group1)-0.0001
                group1.append(adddummy)
                group2.append(adddummy)
                bothgroup=group1+group2
            if samplestd>minstd:
                
                #print group1, group2
                _,grouppval=scipy.stats.kruskal(numpy.array(group1),numpy.array(group2))
                rndpvaluelist=[]
                for j in range(numiter):
                    #print j
                    random.shuffle(bothgroup)
                    rndgrp1=bothgroup[0:len(group1)]
                    rndgrp2=bothgroup[len(group1):]
                    _,rndpval=scipy.stats.kruskal(numpy.array(rndgrp1),numpy.array(rndgrp2))
                    rndpvaluelist.append(rndpval)
                
                outline=ln[0:8]+['%d'%len(group1),'%d'%len(group2),'%5.4f'%medgroup1,'%5.4f'%medgroup2,'%5.4f'%common.median(rndgrp1),'%5.4f'%common.median(rndgrp2),
                         '%5.4f'%samplestd,'%5.4f'%(max(bothgroup)-min(bothgroup)),'%5.4f'%grouppval,'%5.4f'%common.median(rndpvaluelist),medts1,medts2,group1,group2]
                outlinelist.append(outline)
                fout.write('%s\n'%'\t'.join(outline[0:-2]))        
    fout.close()     
    
    reportdir='%s/project/%s/report'%(root_dir,project_name)
    csvhlist=['chr','pos','gene','outflag','exonstartflag','divexon','flowedges','edge used']
    csvhlist+=['Size Group1','Size Group2','Group1 Median','Group2 Median','Random Median1','Random Median2','Sample STD']
    csvhlist+=['Sample Range','Kruskal pvalue','Median Random pvalue','Group1 Median ACT','Group2 Median ACT','Boxplot']
    csvhline='\t'.join(csvhlist)
    worklist=[]
    for line in outlinelist:
        worklist.append([[float(line[16]),1-float(line[15])],line])
    worklist.sort()
    toplist=worklist[:]
    csvfile='%s/top_twolargegroup_report.txt'%reportdir
    csvptr=open(csvfile,'w')
    csvptr.write('#\t%s\n'%csvhline)
    lncnt=0
    for line in toplist:
        lncnt+=1
        l2=ast.literal_eval(line[1][6])
        l2.append(ast.literal_eval(line[1][5]))
        minx=min([min(x) for x in l2])
        maxx=max([max(x) for x in l2])
        #print line
        #print line[1][18],line[1][0],line[1][2],minx,maxx
        urlact1='urlact1'
        urlact2='urlact2'
        imagefilename='%s/image/%s_jitter__%s__%d-%d.pdf'%(reportdir,project_name,line[1][2],minx,maxx)
        ylabel='Fractional Edge Weight'
        xlabel='Group1 vs Group2'
        title='%s__%d-%d\np-value=%s'%(line[1][2],minx,maxx,line[1][16])
        plot.makejitter(line[1][20],line[1][21],[],ylabel,xlabel,title,imagefilename)
        urlimg='%s_jitter__%s__%d-%d.pdf'%(project_name,line[1][2],minx,maxx)
        csvline=line[1][0:18] 
        csvline.append(urlact1)
        csvline.append(urlact2)
        csvline.append(urlimg)
        csvptr.write('%d\t%s\n'%(lncnt,'\t'.join(csvline)))
        genepos=[line[1][2],[minx,maxx],int(line[1][1])]
        mediantslist=[line[1][18],line[1][19]]
        #print mediantslist

    csvptr.close()
    htmlfile=open('%s/top_twolargegroup_report.html'%reportdir,'w')    
    common.csv2html(open(csvfile),htmlfile)
    htmlfile.close()    


if do_index=='expchange':
    import random
    import matplotlib.pylab as plt
    
    exp1='/playpen/sregen/TH_1/metadata/expression_T02S01.txt'
    outdir='/playpen/sregen/TH_2_2_expression'
    fliptype=3
    # flip any two
    # flip largst-any
    # flip largest smallest
    # flip largst-second
    # make fold change histogram
    outfile='%s/exp2_%d.txt'%(outdir,fliptype)
    foldchangelist=[]
    fout=open(outfile,'w')
    outlines=[]
    prevgene=''
    for lntxt in open(exp1):
        ln=lntxt.split('\t')
        if ln[0]==prevgene:
            outlines.append(ln)
        else:
            if prevgene!='':
                if fliptype==1:
                    flippers=random.sample(range(len(outlines)),2)
                if fliptype in [2,3,4]:
                    explist=[float(outlines[i][4]) for i in range(len(outlines))]
                    largest=explist.index(max(explist))
                    seclarge=explist.index(max([explist[x] for x in range(len(outlines)) if x!=largest]))
                    lany=random.choice([x for x in range(len(outlines)) if x!=largest])
                    smallest=explist.index(min(explist))  
                if fliptype==2:
                    flippers=[largest,lany]
                if fliptype==3:
                    flippers=[largest,smallest]    
                if fliptype==4:
                    flippers=[largest,seclarge]    
                outlines[flippers[0]][4],outlines[flippers[1]][4]=outlines[flippers[1]][4],outlines[flippers[0]][4]
                foldchange=float(outlines[flippers[0]][4])/float(outlines[flippers[1]][4])
                if foldchange<1.0:
                    foldchange=1/foldchange
                foldchangelist.append(foldchange)
                for outline in outlines:
                    fout.write('%s\n'%'\t'.join(outline))                
            prevgene=ln[0]
            outlines=[ln]
                
    plt.subplot(111)
    num_bins=1000
    foldchangelist.sort()
    print foldchangelist

    n, bins, patches = plt.hist(foldchangelist, num_bins, normed=1, facecolor='blue', alpha=0.3)
    plt.xlim([1,3])      
    plt.savefig('%s/exp2_%d.pdf'%(outdir,fliptype))

def getffullsig(ffullfile,psig=0.01,fdmth=0.05,covth=1.1):
    genedict={}
    prevgene=''; genedata=[]
    for lntxt in open(ffullfile):
        ln=lntxt.rstrip('\n').split('\t')
        if ln[2]=='gene':
            continue
        gene=ln[2];fdmvalue=float(ln[11]);pvalue=float(ln[12]); wt1=min(float(ln[7]),float(ln[8]));wt2=min(float(ln[9]),float(ln[10]))
        if gene==prevgene or prevgene=='':
            genedata.append((fdmvalue,pvalue,wt1,wt2))
        else:
            genedata.sort()
            genedata.reverse()
            maxfdm,pvalue,wt1,wt2=genedata[0]
            while (min(wt1,wt2)<covth) and len(genedata)>0:
                del genedata[0]
                maxfdm,pvalue,wt1,wt2=genedata[0]
            if min(wt1,wt2)<covth:
                pvalue=0.5
                maxfdm=0.0
            if maxfdm<fdmth:
                pvalue=0.5
            fdmcnt=len(genedata)
            if pvalue*fdmcnt<psig:
                genedict[prevgene]=[wt1,wt2,maxfdm,fdmcnt,pvalue,'DIFF']
            else:
                genedict[prevgene]=[wt1,wt2,maxfdm,fdmcnt,pvalue,'NOT']
            genedata=[(fdmvalue,pvalue,wt1,wt2)]
        prevgene=gene
    return genedict
       
def getgeneexpjsd(jsdfile):
    genejsddict={}
    for lntxt in open(jsdfile):
        ln=lntxt.rstrip('\n').split('\t')
        if ln[0]=='gene':
            continue
        gene=ln[0]; cov1=ln[1] ; cov2=ln[2] ; jsd=ln[5]
        genejsddict[gene]=[cov1,cov2,jsd]
    return genejsddict
 
                
if do_index=='ATHZ1':
    import matplotlib.pylab as plt
    import math
    sroot='/playpen/sregen'
    proot='/playpen/rootFDM/project'
    outdir='/playpen/rootFDM/analysis/ATH'
    jsdfiles=['%s/TH_1/metadata/expression_jsd.txt'%sroot]*9
    fclist=['ATHFC%s'%suf for suf in ['01','05','10','20','20H']]
    fclist=['ATHFC%s'%suf for suf in ['01','05','10','20']]
    fclist=['ATHFC%s'%suf for suf in ['01']]
    jsdfiles+=['%s/%s/metadata/expression_jsd.txt'%(sroot,fc) for fc in fclist]
    suflist=['ATH1_%03d_%05d'%(part,perm) for part in [15,30,60] for perm in [2000,5000,10000]]
    ffullfiles=['%s/ATHpermutation/%s_ALL_ffull.txt'%(proot,suf) for suf in suflist]
    ffullfiles+=['%s/%s/ffull/%s_ALL_ffull.txt'%(proot,fc,fc) for fc in fclist]
    exptup=zip(suflist+fclist,jsdfiles,ffullfiles)
    print '\n'.join(['%s,%s,%s'%x for x in exptup])
    for tup in exptup:
        outfile='%s/%s_combo.txt'%(outdir,tup[0])
        imfile='%s/%s_combo.pdf'%(outdir,tup[0])
        fout=open(outfile,'w')
        genesigdict=getffullsig(tup[2])
        geneexpdict=getgeneexpjsd(tup[1])
        genelist=sorted(geneexpdict.keys())
        x1list=[];y1list=[];x2list=[];y2list=[]
        for gene in genelist:
            outline=geneexpdict[gene][:]
            outline.insert(0,gene)
            if gene in genesigdict:
               x=genesigdict[gene]
               outline+=['%6.4f'%x[0],'%6.4f'%x[1],'%6.4f'%x[2],'%d'%x[3],'%6.4f'%x[4],x[5]]
               xval=min(x[0],x[1]);yval=float(outline[3])
            else:
               outline+=['','','','','','NDAT']
               xval=min(float(outline[1])/25,float(outline[2])/25)  #hack
               yval=float(outline[3])
            if outline[-1]=='DIFF':
                x1list.append(math.log10(max(xval,0.01))); y1list.append(math.log10(yval))
            elif outline[-1]=='NOT':
                x2list.append(math.log10(max(xval,0.01))); y2list.append(math.log10(yval))
            fout.write('%s\n'%'\t'.join(outline))
        plt.clf()
        plt.xlim([-0.1,4])
        msize=3
        plt.plot(x2list,y2list,linestyle='.',marker='o',markerfacecolor='white',markersize=msize)
        plt.plot(x1list,y1list,linestyle='.',marker='o',markerfacecolor='k',fillstyle='full',markersize=msize)
        plt.xlabel('log min coverage')
        plt.ylabel('log jsd/foldchange')        
        plt.savefig(imfile)
        fout.close()


if do_index=='ATHZ2':
    import matplotlib.pylab as plt
    import math
    import os
    sroot='/playpen/sregen'
    proot='/playpen/rootFDM/project'
    outdir='/playpen/rootFDM/analysis/ATH/EXP2'

    suflist=['ATHJSD_ALL_ffull_%02d_%05d'%(part,perm) for part in [15,30,60] for perm in [2000,5000,10000]]
    
    jsdfiles=[]
    ffullfiles=[]
    fclist=[]
    
#    for suf in suflist:
#        ffullfile='%s/ATHJSD/ffull/comparative/%s.txt'%(proot,suf)
#        #print ffullfile
#        if os.path.exists(ffullfile):
#            jsdfiles.append('%s/ATHJSD01/metadata/expression_jsd.txt'%sroot)
#            ffullfiles.append(ffullfile)
#            fclist.append(suf)
    

#    fc='ATHFCD10'
#    fclist.append(fc)
#    jsdfiles.append('%s/%s/metadata/expression_jsd.txt'%(sroot,fc))
#    ffullfiles.append('%s/%s/ffull/%s_ALL_ffull.txt'%(proot,fc,fc))
#    
#    fc='ATHFCD15'
#    fclist.append(fc)
#    jsdfiles.append('%s/%s/metadata/expression_jsd.txt'%(sroot,fc))
#    ffullfiles.append('%s/%s/ffull/%s_ALL_ffull.txt'%(proot,fc,fc))
#    
#    fc='ATHJSD15'
#    fclist.append(fc)
#    jsdfiles.append('%s/%s/metadata/expression_jsd.txt'%(sroot,fc))
#    ffullfiles.append('%s/%s/ffull/%s_ALL_ffull.txt'%(proot,fc,fc))
#            
#    fc='ATHJSD03'
#    fclist.append(fc)
#    jsdfiles.append('%s/%s/metadata/expression_jsd.txt'%(sroot,fc))
#    ffullfiles.append('%s/%s/ffull/%s_ALL_ffull.txt'%(proot,fc,fc))
#
#    fc='ATHJSD18'
#    fclist.append(fc)
#    jsdfiles.append('%s/%s/metadata/expression_jsd.txt'%(sroot,fc))
#    ffullfiles.append('%s/%s/ffull/%s_ALL_ffull.txt'%(proot,fc,fc))
    
    fc='ATHREP'
    fc='ATHREP2'
    fclist.append(fc)
    jsdfiles.append('%s/%s/metadata/expression_jsd.txt'%(sroot,fc))
    ffullfiles.append('%s/%s/ffull/%s_ALL_ffull.txt'%(proot,fc,fc))    
    
    fclist=['ATHJSD18_ALL_ffull_%02d_%05d'%(part,perm) for part in [15,30,60] for perm in [2000,5000,10000]]
    jsdfiles=['%s/%s/metadata/expression_jsd.txt'%(sroot,fc) for fc in ['ATHJSD18']*len(fclist)]
    ffullfiles=['/playpen/rootFDM/project/ATHJSD18/ffull/Many/%s.txt'%fc for fc in fclist]
    

    exptup=zip(fclist,jsdfiles,ffullfiles)
    print '\n'.join(['%s,%s,%s'%x for x in exptup])
    for tup in exptup:
        outfile='%s/%s_combo.txt'%(outdir,tup[0])
        imfile='%s/%s_combo.pdf'%(outdir,tup[0])
        fout=open(outfile,'w')
        genesigdict=getffullsig(tup[2])
        geneexpdict=getgeneexpjsd(tup[1])
        genelist=sorted(geneexpdict.keys())
        x1list=[];y1list=[];x2list=[];y2list=[]
        for gene in genelist:
            outline=geneexpdict[gene][:]
            outline.insert(0,gene)
            if gene in genesigdict:
                x=genesigdict[gene]
                outline+=['%6.4f'%x[0],'%6.4f'%x[1],'%6.4f'%x[2],'%d'%x[3],'%6.4f'%x[4],x[5]]
                xval=min(x[0],x[1])
                yval=float(outline[3])
            else:
                outline+=['','','','','','NDAT']
                xval=min(float(outline[1])/25,float(outline[2])/25)  #hack
                yval=float(outline[3])
            if 'JSD' in tup[0]:
                yval=math.sqrt(yval)
                ylabel='SQRT JSD'
                ylimit=[0.0,1.0]
            if 'FC' in tup[0]:
                ylabel='Transcript usage change'
                ylimit=[0.0,1.0]
            if 'REP' in tup[0]:
                ylabel='Jitter' 
                ylimit=[0.0,0.1]               
            if outline[-1]=='DIFF':
                x1list.append(math.log10(max(xval,0.01))); y1list.append(yval)
            elif outline[-1] != 'DIFF':
                x2list.append(math.log10(max(xval,0.01))); y2list.append(yval)
            fout.write('%s\n'%'\t'.join(outline))
        plt.clf()
        plt.xlim([0.0,4])
        
        msize=3
        plt.plot(x1list,y1list,linestyle='.',marker='o',markerfacecolor='k',fillstyle='full',markersize=msize)
        plt.plot(x2list,y2list,linestyle='.',marker='o',markerfacecolor='white',markersize=msize)
        plt.grid()
        plt.yticks([i/10.0 for i in range(0,11)])   
        if 'REP' in tup[0]:
            plt.yticks([i/10.0 for i in range(0,1)])    
        plt.xticks([i/2.0 for i in range(0,9)])  
        plt.ylim(ylimit)
        plt.xlabel('log min coverage')
        plt.ylabel(ylabel)        
        plt.savefig(imfile)
        print tup,len(x1list),len(x2list)

        fout.close()
        
 
#        xthreshold=0.80
#        ythreshold=0.27
#        xylist=[]
#        numn=0; nump=0
#        for x,y in zip(x1list,y1list):
#            if x<xthreshold:
#                continue
#            if y<ythreshold:
#                xylist.append([x,y,1,0])
#                numn+=1
#            else:
#                xylist.append([x,y,1,1])
#                nump+=1
#        for x,y in zip(x2list,y2list):
#            if x<xthreshold:
#                continue
#            if y<ythreshold:
#                xylist.append([x,y,0,0])
#                numn+=1
#            else:
#                xylist.append([x,y,0,1])        
#                nump+=1
#        xylist.sort()
#        xylist.reverse()
#        fprlist=[]
#        tprlist=[]
#        numfp=0;numtp=0
#        auc=0
#        fpr=0;tpr=0
#        fprlist.append(fpr)
#        tprlist.append(tpr)
#        for x,y,result,real in xylist:
#            if result==1:
#                if real==1:
#                    numtp+=1
#                else:
#                    numfp+=1
#            prevfpr=fpr
#            fpr=numfp*1.0/numn
#            tpr=numtp*1.0/nump   
#            auc+=(fpr-prevfpr)*tpr
#            fprlist.append(fpr)
#            tprlist.append(tpr)
#            print result,real,nump,numn,numtp,numfp,fpr,tpr
#        rocfile='%s/%s_roc.pdf'%(outdir,tup[0])
#        plt.clf()
#        plt.xlim([0.0,1.0])
#        plt.ylim([0.0,1.0])            
#        plt.plot(fprlist,tprlist,'b-')
#        plt.xlabel('False Positive Rate')
#        plt.ylabel('True Positive Rate')   
#        plt.title('%s: AUC=%5.3f'%(tup[0],auc))  
#        plt.savefig(rocfile)
        
        
if do_index=='ATHZ3': 
    import os
    sroot='/playpen/sregen'   
    proot='/playpen/rootFDM/project'   
    
    outpref='ATHFCD10'
    outdir1='%s/%s/metadata'%(sroot,outpref)
    outdir2='%s/%s/ffull'%(proot,outpref)
    
    os.system('mkdir -p %s'%outdir1)
    os.system('mkdir -p %s'%outdir2)
     
    file1='%s/ATHFCD05/metadata/expression_jsd.txt'%sroot
    file2='%s/ATHFCD052/metadata/expression_jsd.txt'%sroot

    genelist1=[x.split('\t')[0] for x in open(file1).readlines()]
    genelist2=[x.split('\t')[0] for x in open(file2).readlines()]
    
    common=[x for x in genelist1 if x in genelist2]
    print common, len(common)
    
    
    outfile='%s/expression_jsd.txt'%outdir1
    fout=open(outfile,'w')
    for lntxt in open(file1):
        fout.write(lntxt)
    for lntxt in open(file2):
        ln=lntxt.split('\t')
        if ln[0] not in genelist1:
            fout.write(lntxt)
    fout.close()    
    
         
    file1='%s/ATHFCD05/ffull/ATHFCD05_ALL_ffull.txt'%proot
    file2='%s/ATHFCD052/ffull/ATHFCD052_ALL_ffull.txt'%proot
    outfile='%s/%s_ALL_ffull.txt'%(outdir2,outpref)
    fout=open(outfile,'w')
    for lntxt in open(file1):
        fout.write(lntxt)
    for lntxt in open(file2):
        ln=lntxt.split('\t')
        if ln[2] not in genelist1:
            fout.write(lntxt)
    fout.close()    
    
if do_index=='ATHZ4': 
    import os
    sroot='/playpen/sregen'   
    proot='/playpen/rootFDM/project'   
    
    outpref='ATHJSD18'
    outdir1='%s/%s/metadata'%(sroot,outpref)
    outdir2='%s/%s/ffull'%(proot,outpref)
    
    os.system('mkdir -p %s'%outdir1)
    os.system('mkdir -p %s'%outdir2)
     
    file1='%s/ATHJSD15/metadata/expression_jsd.txt'%sroot
    file2='%s/ATHJSD03/metadata/expression_jsd.txt'%sroot

    genelist1=[x.split('\t')[0] for x in open(file1).readlines()]
    genelist2=[x.split('\t')[0] for x in open(file2).readlines()]
    
    common=[x for x in genelist1 if x in genelist2]
    print common, len(common)
    
    
    outfile='%s/expression_jsd.txt'%outdir1
    fout=open(outfile,'w')
    for lntxt in open(file1):
        fout.write(lntxt)
    for lntxt in open(file2):
        ln=lntxt.split('\t')
        if ln[0] not in genelist1:
            fout.write(lntxt)
    fout.close()    
    
         
    file1='%s/ATHJSD15/ffull/ATHJSD15_ALL_ffull.txt'%proot
    file2='%s/ATHJSD03/ffull/ATHJSD03_ALL_ffull.txt'%proot
    outfile='%s/%s_ALL_ffull.txt'%(outdir2,outpref)
    fout=open(outfile,'w')
    for lntxt in open(file1):
        fout.write(lntxt)
    for lntxt in open(file2):
        ln=lntxt.split('\t')
        if ln[2] not in genelist1:
            fout.write(lntxt)
    fout.close()    
    
               
def gettprfpr(aroot,sroot,ts,rangemin,rangemax,jsdthreshold,rangetype,covrange):
    import ast

    comfile='%s/%s_combo.txt'%(aroot,ts)
    expfile='%s/%s/metadata/expression_jsd.txt'%(sroot,ts.split('_')[0])
    generesultdict={}
    for lntxt in open(comfile):
        ln=lntxt.rstrip('\n').split('\t')
        gene=ln[0]
        jsdfold=float(ln[3])
        if ln[9]=='NDAT':
            deno=25
            cov1=float(ln[1])/deno
            cov2=float(ln[2])/deno
        else:
            cov1=float(ln[4])
            cov2=float(ln[5])
        if ln[9]=='DIFF':
            diff_flg=1
        else:
            diff_flg=0
        generesultdict[gene]=[jsdfold,cov1,cov2,diff_flg]
    if rangetype=='mincov':
        realdata=[]
        rsltdata=[]
        for gene in generesultdict:
            jsdfold,cov1,cov2,diff_flg=generesultdict[gene]
            mincov=min(cov1,cov2)
            if mincov>=10**rangemin and mincov<10**rangemax:
                if jsdfold>jsdthreshold:
                    realdata.append(1)
                else:
                    realdata.append(0)
                rsltdata.append(diff_flg)
#                if diff_flg==1 and realdata[-1]==0:
#                    print '~~~~~~~~#$#$#$$#$#$#$$#$#$#$',gene,generesultdict[gene]
#                if diff_flg==0 and realdata[-1]==1:
#                    print '********#$#$#$$#$#$#$$#$#$#$',gene,generesultdict[gene]
        combdata=zip(realdata,rsltdata)
        tpr=1.0*combdata.count((1,1))/(realdata.count(1)+1e-20)
        fpr=1.0*combdata.count((0,1))/(realdata.count(0)+1e-20)

    if rangetype=='numtx':
        genetxdict={}
        for lntxt in open(expfile):
            ln=lntxt.split('\t')
            gene=ln[0];txexp=ast.literal_eval(ln[3])
            genetxdict[gene]=len(txexp)
        realdata=[]
        rsltdata=[]
        for gene in generesultdict:
            jsdfold,cov1,cov2,diff_flg=generesultdict[gene]
            mincov=min(cov1,cov2)
            if mincov<covrange[0] or mincov>covrange[1]:
                continue
            numtx=genetxdict[gene]
            if numtx>=rangemin and numtx<rangemax:
                if jsdfold>jsdthreshold:
                    realdata.append(1)
                else:
                    realdata.append(0)
                rsltdata.append(diff_flg)
        combdata=zip(realdata,rsltdata)
        tpr=1.0*combdata.count((1,1))/(realdata.count(1)+1e-20)
        fpr=1.0*combdata.count((0,1))/(realdata.count(0)+1e-20)      
        
        
    if rangetype=='coverageratio':
        realdata=[]
        rsltdata=[]
        for gene in generesultdict:
            jsdfold,cov1,cov2,diff_flg=generesultdict[gene]
            mincov=min(cov1,cov2)
            if mincov<covrange[0] or mincov>covrange[1]:
                continue
            covratio=max(cov1,cov2)/min(cov1,cov2)
            if covratio>=rangemin and covratio<rangemax:
                if jsdfold>jsdthreshold:
                    realdata.append(1)
                else:
                    realdata.append(0)
                rsltdata.append(diff_flg)
        combdata=zip(realdata,rsltdata)
        tpr=1.0*combdata.count((1,1))/(realdata.count(1)+1e-20)
        fpr=1.0*combdata.count((0,1))/(realdata.count(0)+1e-20)               
    
    return [len(combdata),tpr,fpr]
    
if do_index=='ATHZ5':  
    # 'ATHFCD10','ATHFCD15','ATHJSD18' 
    """
    """
    aroot='/playpen/rootFDM/analysis/ATH/EXP2'
    sroot='/playpen/sregen'  
    jsdthreshold=0.02
    foldthreshold=0.1
    
    tstuplist1=[
               ('ATHJSD18','mincov',[x/2.0 for x in range(9)],jsdthreshold),
               ('ATHFCD10','mincov',[x/2.0 for x in range(9)],foldthreshold)]
#    ,
#               ('ATHFCD15','mincov',[x/2.0 for x in range(9)],foldthreshold)]
    tstuplist2=[         
               ('ATHFCD10','numtx',[2,3,4,7,11,1000],foldthreshold),
               ('ATHFCD10','coverageratio',[1.0,1.5,3.0,6.0,10,10000],foldthreshold)           
               ]
    tstuplist3=[         
               ('ATHFCD15','numtx',[2,3,4,7,11,1000],foldthreshold),
               ('ATHFCD15','coverageratio',[1.0,1.5,3.0,6.0,10,10000],foldthreshold)    
               ]

    tstuplist1=[
               ('ATHREP2','mincov',[x/2.0 for x in range(9)],foldthreshold)]

    tstuplist1=[
               ('ATHJSD18_ALL_ffull_60_10000','mincov',[x/2.0 for x in range(9)],jsdthreshold),
               ('ATHJSD18_ALL_ffull_30_10000','mincov',[x/2.0 for x in range(9)],jsdthreshold),
               ('ATHJSD18_ALL_ffull_15_10000','mincov',[x/2.0 for x in range(9)],jsdthreshold),
               ('ATHJSD18_ALL_ffull_60_05000','mincov',[x/2.0 for x in range(9)],jsdthreshold),
               ('ATHJSD18_ALL_ffull_30_05000','mincov',[x/2.0 for x in range(9)],jsdthreshold),
               ('ATHJSD18_ALL_ffull_15_05000','mincov',[x/2.0 for x in range(9)],jsdthreshold),
               ('ATHJSD18_ALL_ffull_60_02000','mincov',[x/2.0 for x in range(9)],jsdthreshold),
               ('ATHJSD18_ALL_ffull_30_02000','mincov',[x/2.0 for x in range(9)],jsdthreshold),
               ('ATHJSD18_ALL_ffull_15_02000','mincov',[x/2.0 for x in range(9)],jsdthreshold)
               ]

    tstuplist1=[
               ('ATHJSD18_ALL_ffull_60_02000','mincov',[x/2.0 for x in range(9)],jsdthreshold),
               ]
    tstuplist1=[
               ('ATHJSD18','mincov',[x/2.0 for x in range(9)],jsdthreshold),
               ('ATHFCD10','mincov',[x/2.0 for x in range(9)],foldthreshold)]

    for i in range(len(tstuplist1)):
        X=[]
        for covrange in [[5,10],[10,30],[30,10000]]:
            ts,rangetype,rangelist,jsdthreshold=tstuplist1[i]
            print '###\n\n',ts,rangetype,rangelist,covrange,'\n'
            rangetup=zip(rangelist[:-1],rangelist[1:])
            t=[];f=[]
            for rangemin,rangemax in rangetup:
                numsampe,tpr,fpr=gettprfpr(aroot,sroot,ts,rangemin,rangemax,jsdthreshold,rangetype,covrange)     
                #print rangemin,rangemax, [numsampe,tpr,fpr] 
                t.append(tpr); f.append(fpr)
            X.append([t,f])
        print X[0]

if do_index=='ATHZ5_9':  
    # 'ATHFCD10','ATHFCD15','ATHJSD18' 
    """
    """
    aroot='/playpen/rootFDM/project/ATHJSD18/ffull/Many'
    sroot='/playpen/sregen'  
    
    jsdthreshold=0.02
    foldthreshold=0.1
    
    tstuplist1=[
               ('ATH1_030_10000','mincov',[x/2.0 for x in range(9)],jsdthreshold),
               ('ATH1_060_10000','mincov',[x/2.0 for x in range(9)],jsdthreshold),
    ]
#    ,
    for i in range(len(tstuplist1)):
        X=[]
        for covrange in [[5,10]]:
            ts,rangetype,rangelist,jsdthreshold=tstuplist1[i]
            print '###\n\n',ts,rangetype,rangelist,covrange,'\n'
            rangetup=zip(rangelist[:-1],rangelist[1:])
            t=[];f=[]
            for rangemin,rangemax in rangetup:
                numsampe,tpr,fpr=gettprfpr(aroot,sroot,ts,rangemin,rangemax,jsdthreshold,rangetype,covrange)     
                print rangemin,rangemax, [numsampe,tpr,fpr] 
                t.append(tpr); f.append(fpr)
            X.append([t,f])
        print X[0]

                
            
def plotnormaldist(outfile,datalist,fdm,title,pvalue):
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib.mlab as mlab
    import math
    import random
    
    plt.clf()
    
    mean=np.mean(datalist)    
    variance=np.var(datalist)
    sigma = math.sqrt(variance)
    x = np.linspace(0,1,100)
    plt.plot(x,mlab.normpdf(x,mean,sigma)*100)
    
    if '10000' in outfile:
        yinc=1
    if '5000' in outfile:
        yinc=2
    if '2000' in outfile:
        yinc=5
            
    
    xlist=datalist[:]
    ylist=[];y=0;xprevbin=0
    for x in xlist:
        xbin=int(x*100)
        if xbin==xprevbin:
            y+=yinc
        else:
            y=yinc
        ylist.append(y)
        xprevbin=xbin
    xlist,ylist=zip(*random.sample(zip(xlist,ylist),500))
    plt.xlim([0,1])
    plt.ylim([0,1000])
    plt.plot(xlist,ylist,'k.',markersize=1)
    plt.plot([fdm,fdm],[0,250],'r-')
    plt.text(0.6,700,'mean=%5.3f\nsigma=%5.3f\npval=%6.4f'%(mean,sigma,pvalue),fontsize=24)
    plt.title(title)
    plt.savefig(outfile)        

if do_index=='ATHZ6': 
    import random
    
    datalist=[random.gauss(0.3,0.1) for _ in range(10000)]
    datalist.sort()
    outfile='/playpen/rootFDM/analysis/ATH/EXP2/plot/test.pdf'
    fdm=0.5
    title='test'
    plotnormaldist(outfile,datalist,fdm,title)

if do_index=='ATHZ7':
    import os
    import cPickle
    
    dir='/playpen/rootFDM/analysis/ATH/EXP2'
    ts='ATHFCD15'
    for partition in [15,30,60]:
        for permutation in [2000,5000,10000]:
            nullpck='%s/nulldist_%s_%2d_%05d.pck'%(dir,ts,partition,permutation)
            if os.path.isfile(nullpck):
                d=cPickle.load(open(nullpck))
                for genepos in d:
                    fdmnulllist,part1,part2,fdmvalue,pvalue=d[genepos]
                    outfile='%s/plot/%s_%s_%2d_%05d.pdf'%(dir,ts,genepos.replace(':','_'),partition,permutation)
                    title='%s_%2d_%05d'%(genepos.replace(':','_'),partition,permutation)
                    print title
                    plotnormaldist(outfile,fdmnulllist,fdmvalue,title,pvalue)

if do_index=='ATHZ8':
    import os
    import cPickle
    
    dir='/playpen/rootFDM/analysis/ATH/EXP2'
    ts='ATHFCD15'
    filelist=[]; preflist=[]
    for partition in [15,30,60]:
        for permutation in [2000,5000,10000]:
            nullpck='%s/nulldist_%s_%2d_%05d.pck'%(dir,ts,partition,permutation)
            if os.path.isfile(nullpck):
                d=cPickle.load(open(nullpck))
                for genepos in d:
                    filelist.append('%s/plot/%s_%s_%2d_%05d.pdf'%(dir,ts,genepos.replace(':','_'),partition,permutation)) 
                    preflist.append('%s_%s'%(ts,genepos.replace(':','_')))
    z=zip(preflist,filelist)
    z.sort()
    preflist,filelist=zip(*z)
    os.system('cd %s/plot/final'%dir)
    for i in range(len(filelist)/9):
        outfile='%s/plot/montage/%s.pdf'%(dir,preflist[i*9])
        infilestr=' '.join(filelist[i*9:(i+1)*9])
        cmdstr1='gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=%s %s'%(outfile,infilestr)
        print cmdstr1
        cmdstr2="pdfnup --nup 3x3 --suffix 'join' --batch %s"%(outfile)
        print cmdstr2        
        os.system(cmdstr1)
        os.system(cmdstr2)
  
if do_index=='ATHZ9':
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib.mlab as mlab
    import math
    import random
    
    plt.figure(figsize=(12,8))
    plt.clf()
    ax=plt.subplot(311)
    
    #numtx
    X=[[[0.5, 0.4444444444444444, 0.45161290322580644, 0.5, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0]], 
       [[0.5865384615384616, 0.5161290322580645, 0.625, 0.6666666666666666, 0.6666666666666666],[0.1, 0.0, 0.07407407407407407, 0.0, 0.0]], 
       [[0.9008620689655172, 0.8866666666666667, 0.9205607476635514, 0.9090909090909091, 0.8571428571428571], [0.18, 0.23333333333333334, 0.14893617021276595, 0.2, 0.0]]]
#    X=[[[0.47619047619047616, 0.425531914893617, 0.5, 0.6470588235294118, 0.0], [0.0, 0.0, 0.14285714285714285, 0.0, 0.0]], 
#       [[0.7207207207207207, 0.6794871794871795, 0.5934959349593496, 0.6086956521739131, 0.75], [0.4, 0.25, 0.5, 0.16666666666666666, 0.0]], 
#       [[0.8787878787878788, 0.8785046728971962, 0.89171974522293, 0.8928571428571429, 0.75], [0.19047619047619047, 0.2857142857142857, 0.3076923076923077,0.25, 0.6666666666666666]]]
    xtickslabels=('2','3','4-6','7-10','11+')
    titles=['coverage between 5-10','coverage between 10-30','coverage >30']
    xlabel='Number of Transcripts'
    ylabel='percentage'
    barlabels=('True +ve','False +ve')
    imgf='/playpen/rootFDM/analysis/ATH/EXP2/plot/realfinal/numtx_fc.pdf'
    
#    #covratio
#    X=[[[0.3, 0.30952380952380953, 0.32, 0.625, 0.7804878048780488], [0.0, 0.0, 0.0, 0.0, 0.0]], 
#       [[0.4166666666666667, 0.6067415730337079, 0.5875, 0.6829268292682927, 0.6333333333333333], [0.0, 0.0, 0.1, 0.2, 0.0]], 
#       [[0.900709219858156, 0.9069767441860465, 0.8895348837209303, 0.9107142857142857, 1.0], [0.21875, 0.32432432432432434, 0.08333333333333333, 0.1, 0.0]]]
##    X=[[[0.0, 0.40425531914893614, 0.5, 0.4666666666666667, 0.7], [0.0, 0.0, 0.0, 0.0, 0.5]], 
##       [[0.46774193548387094, 0.6571428571428571, 0.6804123711340206, 0.75, 0.8285714285714286], [0.5, 0.0, 0.2, 0.6, 0.5]], 
##       [[0.9, 0.8739495798319328, 0.8842105263157894, 0.84375, 0.9107142857142857], [0.1111111111111111, 0.34146341463414637, 0.2857142857142857, 0.2, 0.0]]]
#    xtickslabels=('1.0-1.5','1.5-3.0','3.0-6.0','6.0-10.0','10.0+')
#    titles=['coverage between 5-10','coverage between 10-30','coverage >30']
#    xlabel='Coverage Ratio between samples'
#    ylabel='percentage'
#    barlabels=('True +ve','False +ve')
#    imgf='/playpen/rootFDM/analysis/ATH/EXP2/plot/realfinal/covratio_fc.pdf'



#    X=[[[0.39, 0.42, 0.69, 0.84, 0.91, 0.93,0.98,1.0], [0.11, 0.09, 0.13, 0.15, 0.08,0.54,0,0]],
#       [[0.37, 0.44, 0.72, 0.84, 0.92, 0.94,0.98,1.0], [0.09, 0.08, 0.12, 0.14, 0.08,0.54,0,0]],
#       [[0.40, 0.41, 0.70, 0.81, 0.88, 0.90,0.95,0.92], [0.09, 0.08, 0.12, 0.12, 0.07,0.47,0,0]]]
#    imgf='/playpen/rootFDM/analysis/ATH/EXP2/plot/realfinal/covjsd_partitions.pdf'
#    xtickslabels=('1.0-3.2','3.2-10.0','10.0-31.6','31.6-100','100-316','316-1000','1000-3162','3162+')
#    titles=['Num partitions=15','Num partitions=30','Num partitions=60']
#    xlabel='minimum coverage'
#    ylabel='percentage'
#    barlabels=('True +ve','False +ve')
#
#[[0.33760683760683763, 0.4269230769230769, 0.596551724137931, 0.8668941979522184, 0.9653465346534653, 0.9310344827586207, 0.8695652173913043, 1.0], 
#[0.08163265306122448, 0.0, 0.05063291139240506, 0.11764705882352941, 0.20588235294117646, 0.2413793103448276, 0.3333333333333333, 0.0]]
    
    x=X[0]
    ind = np.arange(len(x[0]))  # the x locations for the groups    
    width = 0.35       # the width of the bars
    bar1=ax.bar(ind, x[0], width, color='g')
    bar2=ax.bar(ind+width, x[1], width, color='r')
    ax.set_ylabel(ylabel)
    ax.set_title(titles[0])
    ax.set_xticks(ind+width)
    ax.set_xticklabels( () )
    ax.set_ylim([0,1])
    ax.yaxis.grid(True)
    ax.set_xlim([0,5])
    ax.legend( (bar1,bar2),  barlabels)

    ax=plt.subplot(312)
    x=X[1]
    ind = np.arange(len(x[0]))  # the x locations for the groups    
    width = 0.35       # the width of the bars
    bar1=ax.bar(ind, x[0], width, color='g')
    bar2=ax.bar(ind+width, x[1], width, color='r')
    ax.set_ylabel(ylabel)
    ax.set_title(titles[1])
    ax.set_xticks(ind+width)
    ax.set_ylim([0,1])
    ax.yaxis.grid(True)
    ax.set_xlim([0,5])
    ax.set_xticklabels( () )

    
    ax=plt.subplot(313)
    x=X[2]
    ind = np.arange(len(x[0]))  # the x locations for the groups    
    width = 0.35       # the width of the bars
    bar1=ax.bar(ind, x[0], width, color='g')
    bar2=ax.bar(ind+width, x[1], width, color='r')
    ax.set_ylabel(ylabel)
    ax.set_title(titles[2])
    ax.set_xticks(ind+width)
    ax.set_xticklabels( xtickslabels )
    ax.set_xlabel(xlabel,fontsize=20)
    ax.set_ylim([0,1])
    ax.yaxis.grid(True)
    ax.set_xlim([0,5])
    

    plt.savefig(imgf)
    
if do_index=='ATHZ10':
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib.mlab as mlab
    import math
    import random
    
    plt.figure(figsize=(12,12))
    plt.clf()
    
    X1=[[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.03977272727272727, 0.03543307086614173, 0.011494252873563218, 0.0, 0.012396694214876033, 0.03968253968253968, 0.09375, 0.0]]
    X2=[[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.01320408163265306, 0.03543307086614173, 0.011494252873563218, 0.0, 0.008178672985781991, 0.03968253968253968, 0.06666666666666667, 0.0]]
    imgf='/playpen/rootFDM/analysis/ATH/EXP2/plot/realfinal/covrepl.pdf'
    xtickslabels=('1.0-3.2','3.2-10.0','10.0-31.6','31.6-100','100-316','316-1000','1000-3162','3162+')
    titles=['']
    xlabel='minimum coverage'
    ylabel='percentage'
    barlabels=('True +ve','False +ve')
    
    ax=plt.subplot(211)
    x=X1[:]
    ind = np.arange(len(x[0]))  # the x locations for the groups    
    width = 0.35       # the width of the bars
    bar1=ax.bar(ind, x[0], width, color='g')
    bar2=ax.bar(ind+width, x[1], width, color='r')
    ax.set_ylabel(ylabel)
    ax.set_title(titles[0])
    ax.set_xticks(ind+width)
    ax.set_xticklabels( () )
    ax.set_ylim([0,0.15])
    ax.set_xlim([0,8])
    ax.yaxis.grid(True)
    
    ax.legend( (bar1,bar2),  barlabels)

    ax=plt.subplot(212)
    x=X2[:]
    ind = np.arange(len(x[0]))  # the x locations for the groups    
    width = 0.35       # the width of the bars
    bar1=ax.bar(ind, x[0], width, color='g')
    bar2=ax.bar(ind+width, x[1], width, color='r')
    ax.set_ylabel(ylabel)
    ax.set_title(titles[0])
    ax.set_xticks(ind+width)
    ax.set_xticklabels( xtickslabels )
    ax.set_ylim([0,0.15])
    ax.set_xlim([0,8])
    ax.yaxis.grid(True)
    
    ax.set_xlabel(xlabel,fontsize=20)

    

    plt.savefig(imgf)   


if do_index=='ATHZ11':
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib.mlab as mlab
    import math
    import random
    
    plt.figure(figsize=(12,6))
    plt.clf()
    
    #X1=[[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.03977272727272727, 0.03543307086614173, 0.011494252873563218, 0.0, 0.012396694214876033, 0.03968253968253968, 0.09375, 0.0]]
    X2=[[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.026785714285714284, 0.007142857142857143, 0.01824817518248175, 0.02, 0.05904059040590406, 0.075, 0.0423728813559322]]
    imgf='/playpen/rootFDM/analysis/ATH/EXP2/plot/realfinal/covreplFC40.pdf'
    xtickslabels=('1.0-3.2','3.2-10.0','10.0-31.6','31.6-100','100-316','316-1000','1000-3162','3162+')
    titles=['']
    xlabel='minimum coverage'
    ylabel='percentage'
    barlabels=('True +ve','False +ve')
    
#    ax=plt.subplot(211)
#    x=X1[:]
#    ind = np.arange(len(x[0]))  # the x locations for the groups    
#    width = 0.35       # the width of the bars
#    bar1=ax.bar(ind, x[0], width, color='g')
#    bar2=ax.bar(ind+width, x[1], width, color='r')
#    ax.set_ylabel(ylabel)
#    ax.set_title(titles[0])
#    ax.set_xticks(ind+width)
#    ax.set_xticklabels( () )
#    ax.set_ylim([0,0.15])
#    ax.set_xlim([0,8])
#    ax.yaxis.grid(True)
#    
#    ax.legend( (bar1,bar2),  barlabels)

    ax=plt.subplot(111)
    x=X2[:]
    ind = np.arange(len(x[0]))  # the x locations for the groups    
    width = 0.35       # the width of the bars
    bar1=ax.bar(ind, x[0], width, color='g')
    bar2=ax.bar(ind+width, x[1], width, color='r')
    ax.set_ylabel(ylabel)
    ax.set_title(titles[0])
    ax.set_xticks(ind+width)
    ax.set_xticklabels( xtickslabels )
    ax.set_ylim([0,0.10])
    ax.set_xlim([0,8])
    ax.yaxis.grid(True)
    ax.legend( (bar1,bar2),  barlabels)
    ax.set_xlabel(xlabel,fontsize=20)

    

    plt.savefig(imgf)   


if do_index=='ATHZ11a':
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib.mlab as mlab
    import math
    import random
    
    plt.figure(figsize=(12,6))
    plt.clf()
    
    #X1=[[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.03977272727272727, 0.03543307086614173, 0.011494252873563218, 0.0, 0.012396694214876033, 0.03968253968253968, 0.09375, 0.0]]
    X2=[[0.3038461538461538, 0.4157303370786517, 0.570957095709571, 0.8410596026490066, 0.9375, 0.84375, 0.851063829787234, 1.0], 
        [0.07692307692307693, 0.0, 0.04878048780487805, 0.11594202898550725, 0.20588235294117646, 0.22580645161290322, 0.3333333333333333, 0.0]]

    imgf='/playpen/rootFDM/analysis/ATH/EXP2/plot/realfinal/covfc10_tpfp.pdf'
    xtickslabels=('1.0-3.2','3.2-10.0','10.0-31.6','31.6-100','100-316','316-1000','1000-3162','3162+')
    titles=['']
    xlabel='minimum coverage'
    ylabel='percentage'
    barlabels=('True +ve','False +ve')
    
#    ax=plt.subplot(211)
#    x=X1[:]
#    ind = np.arange(len(x[0]))  # the x locations for the groups    
#    width = 0.35       # the width of the bars
#    bar1=ax.bar(ind, x[0], width, color='g')
#    bar2=ax.bar(ind+width, x[1], width, color='r')
#    ax.set_ylabel(ylabel)
#    ax.set_title(titles[0])
#    ax.set_xticks(ind+width)
#    ax.set_xticklabels( () )
#    ax.set_ylim([0,0.15])
#    ax.set_xlim([0,8])
#    ax.yaxis.grid(True)
#    
#    ax.legend( (bar1,bar2),  barlabels)

    ax=plt.subplot(111)
    x=X2[:]
    ind = np.arange(len(x[0]))  # the x locations for the groups    
    width = 0.35       # the width of the bars
    bar1=ax.bar(ind, x[0], width, color='g')
    bar2=ax.bar(ind+width, x[1], width, color='r')
    ax.set_ylabel(ylabel)
    ax.set_title(titles[0])
    ax.set_xticks(ind+width)
    ax.set_xticklabels( xtickslabels )
    ax.set_ylim([0,1.0])
    ax.set_xlim([0,8])
    ax.yaxis.grid(True)
    ax.legend( (bar1,bar2),  barlabels)
    ax.set_xlabel(xlabel,fontsize=20)

    

    plt.savefig(imgf)  



if do_index=='ATHREP1':
    import common
    import random
    file1='/playpen/sregen/ATHFCD11/metadata/expression_T01S01.txt'
    file2='/playpen/sregen/ATHFCD11/metadata/expression_T01S02.txt'
    fout=open('/playpen/sregen/ATHFCD11/metadata/expression_jsd.txt','w')
    genedict1={}
    for lntxt in open(file1):
        ln=lntxt.split('\t')
        gene=ln[0].split(':')[0]
        if gene not in genedict1:
            genedict1[gene]=[[],[]]
        genedict1[gene][0].append(float(ln[4]))
        genedict1[gene][1].append(int(ln[5]))
    genedict2={}
    for lntxt in open(file2):
        ln=lntxt.split('\t')
        gene=ln[0].split(':')[0]
        if gene not in genedict2:
            genedict2[gene]=[[],[]]
        genedict2[gene][0].append(float(ln[4]))
        genedict2[gene][1].append(int(ln[5]))
    
    for gene in sorted(genedict1.keys()):
        fout.write('%s\t%d\t%d\t%s\t%s\t%6.4f\n'%
                   (gene,sum(genedict1[gene][1]),sum(genedict2[gene][1]),
                    common.fl2str(genedict1[gene][0]),
                    common.fl2str(genedict1[gene][0]),random.random()/10.0))

if do_index=='ATHREP2':
    import common
    import random
    file1='/playpen/sregen/ATHFCD40/metadata/expression_T01S01.txt'
    file2='/playpen/sregen/ATHFCD40/metadata/expression_T01S02.txt'
    fout=open('/playpen/sregen/ATHFCD40/metadata/expression_jsd.txt','w')
    genedict1={}
    for lntxt in open(file1):
        ln=lntxt.split('\t')
        gene=ln[0].split(':')[0]
        if gene not in genedict1:
            genedict1[gene]=[[],[]]
        genedict1[gene][0].append(float(ln[4]))
        genedict1[gene][1].append(int(ln[5]))
    genedict2={}
    for lntxt in open(file2):
        ln=lntxt.split('\t')
        gene=ln[0].split(':')[0]
        if gene not in genedict2:
            genedict2[gene]=[[],[]]
        genedict2[gene][0].append(float(ln[4]))
        genedict2[gene][1].append(int(ln[5]))
    
    for gene in sorted(genedict1.keys()):
        fout.write('%s\t%d\t%d\t%s\t%s\t%6.4f\n'%
                   (gene,sum(genedict1[gene][1]),sum(genedict2[gene][1]),
                    common.fl2str(genedict1[gene][0]),
                    common.fl2str(genedict1[gene][0]),random.random()/10.0))

   
if do_index=='gtf1':
    import common
    def attrtxt2attrdict(attrtxt):
        attrlist=attrtxt.split(';')
        attrtuplist=[attr.strip(' ').split(' ',1) for attr in attrlist]
        attrtuplist=[(attrtup[0],attrtup[1][1:-1]) for attrtup in attrtuplist if len(attrtup)==2]
        return dict(attrtuplist)
    
    
    def gtf2genetranscriptdict(gtffile,type):
        if type=='E':
            geneidx='gene_name'
        elif type=='F':
            geneidx='gene_id'
        genetranscriptdict={}
        lnum=0
        for linetxt in open(gtffile):
            lnum+=1
            if lnum%200000==1:
                message='Processing GTF file for transcripts at line %d'%lnum
                common.printstatus(message,'S',common.func_name())
            line=linetxt.rstrip('\n').split('\t')
            if line[2]!='exon':
                continue
            attrdict=attrtxt2attrdict(line[8])
            if geneidx in attrdict:
                gene=attrdict[geneidx]
            else:
                continue
            transcript=attrdict['transcript_id']
            chrnm=line[0]
            if line[6]=='+':
                strand=1
            else:
                strand=0
            #print gene, transcript, chrnm, strand
            if gene not in genetranscriptdict.keys():
                genetranscriptdict[gene]={}
            if transcript not in  genetranscriptdict[gene].keys():
                genetranscriptdict[gene][transcript]=(chrnm,strand,[])
            genetranscriptdict[gene][transcript][2].append((int(line[3]),int(line[4])))
        for gene in genetranscriptdict:
            for transcript in genetranscriptdict[gene]:
                genetranscriptdict[gene][transcript][2].sort()
        message='Completed processing GTF file for transcripts'
        common.printstatus(message,'S',common.func_name())
        
        genetranscriptdict_file='%s/genetranscriptdict.pck'%self.dir
        cPickle.dump(genetranscriptdict,open(genetranscriptdict_file,'w'))
        genetranscriptdict=self._fixgenetranscriptdict()
        return genetranscriptdict       
            
    gtffile='/playpen/rootFDM/annotation/Arabidopsis_thaliana.TAIR10.21.gtf'
    type='E'
    gtx=gtf2genetranscriptdict(gtffile,type)
    genelist=gtx.keys()
    print genelist[0:10]
    print gtx[genelist[0]]
    
if do_index=='a':
    import math
    flow1=[0.5,0.5,0]
    flow2=[0,0,1]

    fdm=0
    if len(flow1)!=len(flow2):
        message='Flow sizes are different %s:%s'%(common.fl2str(flow1),common.fl2str(flow2))
        common.printstatus(message,'E',common.func_name())   
    else:
        for i in range(len(flow1)):
            fdm+=1.0/2*math.fabs(flow1[i]-flow2[i]) 
    print fdm
    
if do_index=='RMOUSE1':
    """
    Consistent differences
    """                 
    #50,53,56;59,62,65;68,71,74,77
    projroot='/home/darshan/rootFDM/project/MOUSE01/report'
    top_file='%s/top_report.txt'%projroot
    
    lncnt=0
    for lntxt in open(top_file):
        if lncnt==0:
            continue
        lncnt+=1
        ln=lntxt.rstrip('\n').split('\t')
        if ln[50]=='SIG' and ln[53]=='SIG' and ln[56]=='SIG':
            print ' 3: ',lncnt, lntxt
        if ln[59]=='SIG' and ln[62]=='SIG' and ln[65]=='SIG':
            print ' 5: ',lncnt, lntxt         
        if ln[68]=='SIG' and ln[71]=='SIG' and ln[74]=='SIG' and ln[74]=='SIG':
            print '10: ',lncnt, lntxt    
        
if do_index=='RMOUSE1':
    """
    Consistent differences
    """                 
    #50,53,56;59,62,65;68,71,74,77
    projroot='/home/darshan/rootFDM/project/MOUSE01/report'
    top_file='%s/top_report.txt'%projroot
    
    lncnt=0
    for lntxt in open(top_file):
        if lncnt==0:
            continue
        lncnt+=1
        ln=lntxt.rstrip('\n').split('\t')
        if ln[50]=='SIG' and ln[53]=='SIG' and ln[56]=='SIG':
            print ' 3: ',lncnt, lntxt
        if ln[59]=='SIG' and ln[62]=='SIG' and ln[65]=='SIG':
            print ' 5: ',lncnt, lntxt         
        if ln[68]=='SIG' and ln[71]=='SIG' and ln[74]=='SIG' and ln[74]=='SIG':
            print '10: ',lncnt, lntxt        



if do_index=='RMOUSE2':
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pylab as plt
    
    def plotfile(xlist,ylist, ptitle, xlabel, ylabel, img_file):
        corr = statistics.correlation(xlist, ylist, method = "Pearson") 
        plt.ioff() 
        plt.plot(xlist,ylist,'.')
        plt.xlabel(xlabel,{'fontsize': 'x-large'})
        plt.ylabel(ylabel,{'fontsize': 'x-large'})
        plt.title('%s\nCorrelation: %4.3f'%(ptitle,corr),{'fontsize': 'x-large'})
        plt.axis([0.0,1.0,0.0,1.0])
        plt.savefig(img_file)
        plt.clf()     
        return corr
    """
    Cufflinks-FDM
    """    
    cuffdir='/data/Porter_RNAseq/Brian/Analysis_temp'
    #8,12
    generpkdict={}
    generpkdict['FDM']={}
    generpkdict['CUFF']={}
    lncnt=0
    for lntxt in open('%s/CuffDiff-WIII_vs_MIII/genes.fpkm.tracking'%cuffdir):
        if lncnt==0:
            continue
        ln=lntxt.rstrip('\n').split('\t')
        generpkdict['CUFF'][ln[0]]=[float(ln[8]),float(ln[12])]
    lncnt=0
    for lntxt in open('%s/CuffDiff-WV_vs_MV/genes.fpkm.tracking'%cuffdir):
        if lncnt==0:
            continue
        ln=lntxt.rstrip('\n').split('\t')
        generpkdict['CUFF'][ln[0]]+=[float(ln[8]),float(ln[12])]        
    lncnt=0
    for lntxt in open('%s/CuffDiff-WX_vs_MX/genes.fpkm.tracking'%cuffdir):
        if lncnt==0:
            continue
        ln=lntxt.rstrip('\n').split('\t')
        generpkdict['CUFF'][ln[0]]+=[float(ln[8]),float(ln[12])]        
    fdmdir='/home/darshan/rootFDM/dataout'
    tslist=['W03_1','W03_2','W03_3','W05_1','W05_2','W05_3','W10_1','W10_2','W10_3','W10_4',
            'M03_1','M03_2','M03_3','M05_1','M05_2','M05_3','M10_1','M10_2','M10_3','M10_4']
    for ts in tslist:
        infile='%s/%s/MOUSE01_%s.rpk'%(fdmdir,ts,ts)
        lncnt=0
        for lntxt in open('%s/CuffDiff-WV_vs_MV/genes.fpkm.tracking'%cuffdir):
            if lncnt==0:
                continue
            ln=lntxt.rstrip('\n').split('\t')
            if ln[0] not in generpkdict['FDM']:
                generpkdict['FDM'][ln[0]]=[]
            generpkdict['FDM'][ln[0]]+=[float(ln[5])]                
    genelist=list(set(generpkdict['FDM'].keys()).intersection(set(generpkdict['CUFF'].keys())))
    print len(genelist)
    
    imagedir='/home/darshan/rootFDM/project/MOUSE01/analysis/image'
    for num in range(20):
        ts=tslist[num]
        type=num/10
        rep=num%10
        lob=min(rep%3,2)
        cuffi=lob*2+type
        xlist=[]
        ylist=[]
        for gene in genelist:
            xlist.append(math.log10(generpkdict['CUFF'][gene][cuffi]))
            ylist.append(math.log10(generpkdict['FDM'][gene][num]))
        xlabel='CUFFLINKS'; ylabel='FDM'
        ptitle=ts
        img_file='%s/%s.pdf'%(imagedir,ts)
        corr=plotfile(xlist,ylist, ptitle, xlabel, ylabel, img_file)
        print ts, corr
        
if do_index=='RMOUSE3':
    import numpy as np
    threshold=10
    ffastdir='/home/darshan/rootFDM/project/MOUSE01/ffast'
    tslist=['W03_1','W03_2','W03_3','W05_1','W05_2','W05_3','W10_1','W10_2','W10_3','W10_4',
            'M03_1','M03_2','M03_3','M05_1','M05_2','M05_3','M10_1','M10_2','M10_3','M10_4']
    for i in range(19):
        for j in range(i+1,20):
            ffastfile='%s/MOUSE01__%s__%s_ALL_fdm.txt'%(ffastdir,tslist[i],tslist[j])
            fdmlist=[]
            lncnt=0
            for lntxt in open(ffastfile):
                if lncnt==0:
                    continue
                lncnt+=1
                ln=lntxt.rstrip('\n').split('\t')
                wt=float(ln[-2]); fdm=float(ln[-1])
                if wt>threshold:
                    fdmlist.append(fdm)
            mean,var=np.mean(fdmlist),np.var(fdmlist)
            print tslist[i],tslist[j], mean, var
            
            
            
    
    
    
    
    
    
    
    
    
    
    
   
print "done"        

