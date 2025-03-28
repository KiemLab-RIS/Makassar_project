import sys
import os
from collections import defaultdict
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
#
# assay dna file from liftover of off-target amplicon
# Assay_MacacaMakas.bed
# the assay data is always + strand and matches fastq
#
# 0	chr10	10633821	10634030	GGAAGTTTACAAGAACCATCGTTGTGCCAAAATTATCCTATATAAGGATATTAGAAAACCAAAAAACTAGTTCAA...
# 1	chr10	12691171	12691361	GGATTACCGATGTGAGCCGCTGTGCCTGGCCCAGCCATGTTTTTTGGGTTGAATGGAACTCACTGCTTTTCTTCT...
# 2	chr10	13689929	13690143	AGGCTAACGCTGGAGGATCTCTTGATCCTAGGATTTTGAGACCATCCTGGACAACAGAGTGAGATTCCATCTCTT...
#
# off target sequences also from liftover hglft_genome_2b26_eea760.bed
# but this sequence also has strand info added based on whether the
# off-target sequence aligns to + or - strand
# 
#
# GTGCATCTGACTCCTGAGGAGAA	-
# TCTCATCTGACTCCTGAGGACAA	-
# TTATCCTCAGGAGTCAGATGCTG	+
# TTCTTCTCTGGAGTCAGATGGGG	+
# TTCTCCACAGAAGTCAGATGGGA	+
#---------------------------------------------------------
#
# reverse complement
#
#---------------------------------------------------------
def dnacomp(ch):
  if ch == 'A': return 'T'
  if ch == 'T': return 'A'
  if ch == 'C': return 'G'
  return 'C'
#---------------------------------------------------------
#
# reverse complement
#
#---------------------------------------------------------
def revcomp(s):
  l = [dnacomp(s[x-1]) for x in range(len(s),0,-1)]
  l = ''.join(l)
  return(l)
#
# protospace or guideseq
# the last three bases are PAM  'CAC'
#
# if this aligns to + strand then protospacer location is [guideSeqLoc:guideSeqLoc+20]
# if this aligns to - strand then protospacer location is [guideSeqLoc+3:guideSeqLoc+23]
#
#  for + strand example
#                              TTCTCCTCAGGAGTCAGATGCAC
#   GGAAGTTTACAAGAACCATCGTTGTGCTTCTCCTCAGGAGTCAGATGCACAAAATTATCCTATATTTCTCCTCAGGAGTCAGATGCAC
#                              ********************
#
# for - strand example
#                              GTGCATCTGACTCCTGAGGAGAA
#   GGAAGTTTACAAGAACCATCGTTGTGCGTGCATCTGACTCCTGAGGAGAAAAATTATCCTATATTTCTCCTCAGGAGTCAGATGCAC
#                                *********************
#
#
guideSeq = 'TTCTCCTCAGGAGTCAGATGCAC'
guideSeqR = revcomp(guideSeq)
#
# params
#
if len(sys.argv) != 5:
  print(f'Useage: rhAmp.py primerList offTargetList fastqFile outDir')
  quit()

fPrimerList = sys.argv[1]
fOffTarget  = sys.argv[2]
r1          = sys.argv[3]
outDir      = sys.argv[4]
log         = 'log/' + outDir + '_prog_log.txt'
#
# read off-target seq 
#
tf = [x.strip().split('\t') for x in open(fOffTarget,'r')] 
if len(tf) == 0:
  print(f'bad offTarget file {fOffTarget}')
  quit()
#
#
#
offTarget = []
#
#
#
for index,offLine in enumerate(tf):
  offTarget.append((offLine[0],offLine[1],index))

print(f'read {len(offTarget)} off target seqs')
#
# read primers
#
primersT = [x.strip().split('\t') for x in open(fPrimerList,'r')]
if len(primersT) == 0:
  print(f'bad primer file {fPrimerList}')
  quit()
#
#
# why doesn't base 1 align?
#
ampliconsT = [x[4][1:] for x in primersT]
primersT = [(x[4][1:16],x[4][-15:]) for x in primersT]
#
# add revcomp,index,f/r
#
validPrimersF = {}
validPrimersR = {}
amplicons = {}
print('process amplicon and primers')
for index,(priF,priR) in enumerate(primersT):
  validPrimersF[priF] = (index,'+',priR)
  validPrimersR[priR] = (index,'+',priF)
  amplicons[ampliconsT[index]] = (index,'+')

print(f'loaded {len(validPrimersF)} primers from {fPrimerList}')

#
# primary storage is ampSeqD
#
# key = [forwardPrimer:reversePrimer]
#
# value = list of all sequences in a defaultdict(int)
#
ampSeqD = {}
#
# ampSeqRef is the amplicon sequence for mathcing priF : priR
#
ampSeqRef = {}
#
# init every ampSeq with reference sequence
#
for keySeqA in amplicons.keys():
  key = keySeqA[:15] + keySeqA[-15:]
  d = defaultdict(int)
  d[keySeqA] = 1
  ampSeqRef[key] = keySeqA
  ampSeqD[key] = d
  
badSeq = defaultdict(int)

#--------------------------------------------------------------------------------
#
# compare two sequences
#
#--------------------------------------------------------------------------------
def compareSeq(sa,sb):
  if len(sa) != len(sb):
    return(100)
  error = 0
  for i in range(0,len(sa)):
    if sa[i] != sb[i]:
      error += 1
  return error
#--------------------------------------------------------------------------------
#
# checkQuality: how many bases are below Q
#
#--------------------------------------------------------------------------------
def checkQuality(qualSeq):
  #
  # count how many bases are below low quality
  #
  fail = 0
  i = 0
  for char in qualSeq:
    q = ord(char) - 33   #ord('!')
    if (q < 14):         # change from 15
      fail += 1
  return fail
#
#
#
lines = 0
goodLines = 0
#
# read fastq files
#
with open(r1,'r') as fh1:
  while True:
    header = fh1.readline()
    if not header:
      break

    seq   = fh1.readline().strip()
    plus  = fh1.readline()
    q     = fh1.readline()

    #e1 = checkQuality(q)
    
    primerF = seq[:15]
    primerR = seq[-15:]
    #
    # look  for 'primer' seq in white(valid) list
    #
    sequenceMatch = False
    primerWhiteF = ''
    primerWhiteR = ''
    if primerF in validPrimersF.keys() and primerR == validPrimersF[primerF][2]:
      #print(f'exact match')
      #
      # exact match
      #
      primerWhiteF = primerF
      primerWhiteR = primerR
      priIndex,priStrand,t_priR = validPrimersF[primerF]
      sequenceMatch = True
    elif primerF in validPrimersF.keys():
      #
      # found forward primer but reverese may have error
      #
      priIndex,priStrand,priReverse = validPrimersF[primerF]
      error = compareSeq(primerR,priReverse)
      if error <= 1:
        primerWhiteF = primerF
        primerWhiteR = priReverse
        sequenceMatch = True
        #print(f' repair reverse, error = {error}')
        #print(f'{primerF}')
        #print(f'{primerR}')
        #print(f'{priReverse}')
        #print(f'{seq}')
    elif primerR in validPrimersR.keys():
      #
      # found reverse but forward may have error
      #
      priIndex,priStrand,priForward = validPrimersR[primerR]
      error = compareSeq(primerF,priForward)
      if error <= 1:
        primerWhiteF = priForward
        primerWhiteR = primerR
        sequenceMatch = True
        #print(f'repair forward')
    #
    # is there is a primer match, store sequence in dict where
    # key is fPrimer + rPrimer
    #
    if sequenceMatch:
      goodLines += 1
      saveSeq = f'{seq}_{priIndex}_{priStrand}'
      primerKey = primerWhiteF + primerWhiteR
      if primerKey in ampSeqD.keys():
        ampSeqD[primerKey][saveSeq] += 1
      else:
        print(f'primer key {primerKey} not found in ampSeq')
        d = defaultdict(int)
        d[saveSeq] += 1
        ampSeqD[primerKey] = d
    else:
      bad = primerF + " " + primerR
      #print(f' no match for {primerF}    {primerR}')
      #print(f'{primerF in validPrimersF.keys()}  {primerR in validPrimersR.keys()}')
      badSeq[bad] += 1
    lines += 1
    if lines % 1000 == 0:
      print(f'{lines} number of amps = {len(ampSeqD)}  good lines = {goodLines}')

    if lines > 1000000000:
      print("End on line count")
      break
#
# record top bad
#
sd = sorted(badSeq.items(),key=lambda x:x[1],reverse=True)
for i in range(0,40):
  print(f'bad seq[{i}]   {sd[i][0]} {sd[i][1]}')
#
#  generate output
#
offTargetSeqFound = [False for x in offTarget]

with open(log,'w') as fLog:
  fLog.write(f'Log file for {outDir}\n')
  fLog.write(f'PrimerKey : SeqMatchRef\n')
  #
  # for each primer set...a group of sequences is stored
  #
  for pSeq,ampD in ampSeqD.items():
    #
    # if only one seq, it is jsut reference
    # 
    if len(ampD) == 1:
      continue
    #
    # make dir tostore
    #
    if not os.path.exists(outDir):
      os.mkdir(outDir)

    outFile = outDir + '/' + pSeq + '.tsv'
 
    totalSeq = 0
    #
    #
    #
    sd = sorted(ampD.items(),key=lambda x:x[1],reverse=True)
    # sum total reads
    for k,v in sd:
      totalSeq += v
    #
    # get ref sequence
    #
    ref = ampSeqRef[pSeq]
    #
    # are any off-targets in ref
    #
    offLoc = -1
    offStrand = '+'
    offIndex = -1
    for tIndex,offLine in enumerate(offTarget):
      offSeq = offLine[0]
      strand = offLine[1]
      oIndex = offLine[2]
      if offSeq in ref:
        offLoc = ref.find(offSeq)
        offTargetSeqFound[tIndex] = True
        offStrand = offStrand = strand
        offIndex = oIndex
        break
    #
    # write file
    #
    with open(outFile,'w') as fh:
      rsp = len(ref)-len(pSeq)
      # primers
      fh.write(f'pri\t{pSeq[:15]}{" "*rsp}{pSeq[-15:]}\n')
      fh.write(f'SeqCount\t{totalSeq}\n')
      # off-target protospacer
      if offStrand == '+':
        fh.write(f'pro.f\t{" " * offLoc}{guideSeq}\n')
      else:
        fh.write(f'pro.r\t{" " * offLoc}{guideSeqR}\n')
      if offLoc != -1:
        if offStrand == '+':
          fh.write(f'OL\t{offLoc}\n')
          fh.write(f'off\t{" " * offLoc}{offSeq}   index {offIndex}  strand  {offStrand}\n')
        else:
          fh.write(f'OL\t{offLoc+3}\n')
          fh.write(f'off\t{" " * offLoc}{offSeq}   index {offIndex}  strand  {offStrand}\n')
      # reference (largest seq)
      fh.write(f'ref\t{ref}\n')
      #
      # write all seq
      #
      for k,v in sd:
        # get rid of _index_strand
        k = k[:k.find('_')]
        # string to list
        # for putting (.) in seq
        # 
        #
        # seq = [dna for dna in k]
        # comp to ref
        #  iMax = min(len(ref),len(seq))
        #  for i in range(0,iMax):
        #    if seq[i] == ref[i]:
        #      seq[i] = '.'
        #  seq = ''.join(seq)
        seq = k
        fh.write(f'seq\t{seq}\t{v}\t{v/totalSeq:2.4f}\n')
    rMatch = sd[0][1] / totalSeq
  
    fLog.write(f'{pSeq}\t{rMatch}\n')
  #
  # offTarget looks like
  #------------------------
  #  seq           index0
  #  revComp(seq)  index0
  #  seq           index1
  #  revComp(seq)  index1
  #
  fLog.write(f'OffTargetSequences found\n')
  found = 0
  index = 0
  for index in range(0,len(offTargetSeqFound),2):
    if offTargetSeqFound[index]:
      fLog.write(f'found forward   {offTarget[index]}\n')
      found += 1
    elif offTargetSeqFound[index+1]:
      fLog.write(f'found reverse   {offTarget[index+1]}\n')
      found += 1

  fLog.write(f'found total of {found} offTargets in this run\n\n\n')
  fLog.write(f'OffTargetSequences not found\n')

  notFound = 0
  for index in range(0,len(offTargetSeqFound),2):
    if offTargetSeqFound[index] == False and offTargetSeqFound[index+1] == False:
      fLog.write(f'Not Found    {offTarget[index]}\n')
      notFound += 1

  fLog.write(f'did not find total of {notFound} offTargets in this run\n')

