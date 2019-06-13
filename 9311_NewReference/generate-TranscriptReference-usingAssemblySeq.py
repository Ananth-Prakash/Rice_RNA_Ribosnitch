# Z:\YD\ANALYSIS\OSA\HONG-RICE\New-Ref-9311-12022018
"""
 NOTE: Genes must have isoforms !!
"""
import sys
import itertools
import re
from Bio.Seq import Seq
from collections import defaultdict
import pprint 
from itertools import groupby
import gffutils

def getEasy(f):
    for header, group in itertools.groupby(f, key = lambda x: x.startswith(">") ):
      if header:
         line = next(group)
         tag = line[1:].strip()
      else:
         sequence = ''.join(line.strip() for line in group)
         yield tag, sequence.strip()

def loadWholeSeq( fastaFiler = "93-11_AltReferenceGenome.fasta" ):
   seqD = {}
   with open(fastaFiler, 'r' ) as inpFast:
    for h, seq in getEasy(inpFast):
        seqD[h]= seq
   return seqD

def orderExons(isoforms_exon):
    isoforms_exon = list(isoforms_exon)
    isok = [] # [ s.split(':exon_')[0]: s for s in isoforms_exon } #
    ## keep the order 
    sign = set()
    pipeD = defaultdict(list)
    for pack in isoforms_exon:
        s, info = pack 
        if ':' not in s:
           print 'Error:', pack
           exit()    
        gene_id, exid  = s.split(':')
        sign.add(info[0])  # this to ensure we have same sign exons
        pipeD[gene_id].append(pack)

        if gene_id  not in isok: # these is to keep order in Fasta file 
           isok.append(gene_id)               

    assert len(pipeD)  == len(isok), ' failed number of exons check gff file again!'
    assert len(sign) == 1,   'exon direction is bad for this gene: ' +  str(gene_id)
    return  sign, isok, pipeD      
#-------------------------------------------

def fetchRegion(pipe,  seqD):
    #child.strand, child.seqid, child.start, child.stop
    st=''
    for head, reg in pipe:
        strand, seqid, istart, istop =  reg
        st += seqD[seqid][istart -1: istop]
    return st  


#------ driver --------------------------------
rice_gff3_massive = "all_gff3_MSU7.gff3"
G = gffutils.FeatureDB('rice-sanit.db', keep_order = True)
#============== Loading seq ==============================
seqD = loadWholeSeq( fastaFiler = "93-11_AltReferenceGenome.fasta" ) # 12 chromosomes
minus_stranded_genes = 0
plus_stranded_genes  = 0

megaSign = defaultdict(list)
iso_form_order = {}  # for deferred printing 

with open('rice-9311-chr1-12-refJitu-isformOrdered.fasta','w') as  outf:
 for gene in G.features_of_type('gene', order_by='start'):

    if gene.seqid not in seqD: # we dont need gene which are not present in Fasta refence
       continue
  
    # get all grandchildren, only keeping the exons adding
    isoforms_exon = []
    for child in G.children(gene.id, 2,  order_by='start'):
       if child.featuretype == 'exon':
          isoforms_exon.append([child.id, [child.strand, child.seqid, child.start, child.stop] ] ) 
    #=============================================================================================     
    sign, gen_ord, pipeD = orderExons(isoforms_exon) #:=
    megaSign[gene.id].append(list(sign)[0])
    
    for g in gen_ord:
        st = fetchRegion(pipeD[g], seqD)  

        if '-' in sign:# - gene we are taking rc
           st = str( Seq(st).reverse_complement() )   
           iso_form_order[g] = st
           minus_stranded_genes += 1 

        elif '+' in sign:# + gene
           iso_form_order[g] = st
           plus_stranded_genes  += 1

        else:
           print "strange sign: ", sign 

#--------------------------------
 iso_names = iso_form_order.keys()
 iso_names.sort() 
 for k in iso_names: # lexiographic order
        seqa = iso_form_order[k]
        outf.write('>' + k + '\n')
        outf.write(   seqa + '\n')
       
print "minus_stranded_genes:", minus_stranded_genes 
print "plus_stranded_genes:",  plus_stranded_genes 
print "Total isoformed Gene: ", plus_stranded_genes + minus_stranded_genes 

#---------isoform sign intergrity check --------------
print "Pure Genes:", len(megaSign)
for puregene, signal in megaSign.items():
    if len(set(signal)) > 1:
       print 'Bad gene:', puregene
print ("DONE")
#--------------------------------------------------
