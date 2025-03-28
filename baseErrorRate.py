from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import Align
import sys
import os
from collections import defaultdict
inputDir = sys.argv[1]
#
#
#
aligner = Align.PairwiseAligner()
aligner.mode = 'global'
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -2.5
aligner.extend_gap_score = -0.1
aligner.target_end_gap_score = 0.0
aligner.query_end_gap_score = 0.0
#
#  example aligner output
#
#TATCTCCTCCTGCTTGTGGATCAGCTGGATAGGTACCTAAGCTATCTGTGGATTGGCCCA
#TATCTCCTCATCTTGTGGATCAGCTGGATAGGTACCTAAGCTATCTGTGGATTGGCCCAA
#TATCTCCTCCTGCTTGTGGATCAGCTGGATAGGTACCTAAGCTATCTGTGGATTGGCCCA
#|||||||||.|-||||||||||||||||||||||||||||||||||||||||||||||||
#TATCTCCTCAT-CTTGTGGATCAGCTGGATAGGTACCTAAGCTATCTGTGGATTGGCCCA
#
#(((0, 11), (12, 192)), ((0, 11), (11, 191))
#
# in this case, indel at loc 11
#


#
# classify indel
#    inputRanges = (((0, 139), (140, 209)), ((0, 139), (139, 208)))
#
def classifyIndel(inputRanges):
  src = inputRanges[0]
  dst = inputRanges[1]
  lSrc = 0
  lDst = 0
  for x1,x2 in src:
     lSrc = x2

  for x1,x2 in dst:
    lDst = x2

  #
  # dst > src = (+) deletion
  # src > dst = (-) insertion
  #
  return(lDst - lSrc)
#
#
#
#
def compSeq(s1,s2,maxError):
  # may have been an indel?
  if len(s1) != len(s2):
    return None
  # how many base mis-match
  baseMismatch = []
  if s1 == s2:
    return []
  errors = 0
  for i in range(0,len(s1)):
    if s1[i] != s2[i]:
      baseMismatch.append((i,s1[i],s2[i]))
      if len(baseMismatch) >= maxError:
        return baseMismatch

  return baseMismatch


def processFile(fIn,fOut,fOut2,fOutIndel):

  fh = open(fOut,'w')
  fh2 = open(fOut2,'w')
  fhIndel = open(fOutIndel,'w')

  seqL = [f.strip().split('\t') for f in open(fIn,'r')]
  
  ref = ''
  totalSeqCount = 0
  refCount = 0
  # store mismatch errors and indel by location
  perBaseError = []
  # store all base counts
  baseCount = []
  # found proto?
  protoSpacerLoc = -1
  #
  # record indels
  #
  indelDict = defaultdict(list)
  #
  # for each unique sequence
  #
  for index in range(0,len(seqL)):
    line = seqL[index]

    if line[0] == "OL":
      protoSpacerLoc = int(line[1])
    elif line[0] == 'ref':
      ref = line[1]
      for x in range(0,len(ref)):
        perBaseError.append(defaultdict(int))
        baseCount.append(defaultdict(int))
      #fh.write(f'{ref}\n')
      #fh.write(f'sequence count = {totalSeqCount}\n')
    elif line[0] == 'seq':
      seq = seqL[index][1]
      seqCount = int(seqL[index][2])
      totalSeqCount += seqCount
      #fh.write(f'{ref}\n')
      #fh.write(f'{seq}\n')
      #fh.write(f'sequence count = {seqCount}\n')
      #
      # get a list of base errors [loc,error] from a direct comparison
      # 
      # a large number of errors = prob indel
      #
      # empty list means lengths different = prob indel
      #
      # could also run alignment first? just slower
      #
      Max_Error = 6
      baseErrors = compSeq(ref,seq,Max_Error)
      if baseErrors == None or len(baseErrors) >= Max_Error:
        #
        # possible indel...run alignet
        #
        #print(f'possible indel, count = {seqCount}')
        # is an alignment possible?
        if abs(len(seq) - len(ref)) > 10:
          # very bad seq
          if seqCount > 100:
            print(f'Unalignable sequence\n')
            print(f'{ref}')
            print(f'{seq}')
            print(f'{seqCount}')
            pass
        else:
          alignments = aligner.align(ref, seq)
          if len(alignments) > 1000:
            #print(f'WARNING:  # alignments = {len(alignments)}')
            pass
          sd = sorted(alignments)
          # a is highest ranked alignment
          a = sd[0]
          #fh.write(f'{a}\n')
          #fh.write(f'{a.aligned}\n')
          # find index of indexl
          indel = a.aligned[0]
          if len(indel) > 1:
            indelStart = indel[0][1]
            indelType = classifyIndel(a.aligned)
            indelDict[indelStart].append((indelType,seqCount))
            #print(f'a = ')
            #print(f'{a}')
            #print(f'a.aligned = {a.aligned}')
            #print(f'indel = {indel}')
            #print(f'indexStart = {indelStart}')
            #
            #  this will add to the base error rate of this position so
            #  an indel will be counted twice  [once for indel and once for '-' in aligned seq]
            #
            perBaseError[indelStart]['I'] += seqCount
            #
            # all base stats
            # 
            # example
            # a =
            # GCTGGTGAGTTTTAATTGATGA
            # ||||||||||||||||||||||
            # GCTGGTGAGTTTTAATTGATGA
            #a.aligned = (((0, 139), (140, 209)), ((0, 139), (139, 208)))
            #indel source = ((0, 139), (138, 208))
            #indexStart = 139
            dstAlign = a.aligned[0]
            srcAlign = a.aligned[1]
            
            alignSeq = ['-' for x in range(0,len(ref))]
            
            for alIndex in range(0,len(srcAlign)):
              srcRange = srcAlign[alIndex]
              dstRange = dstAlign[alIndex]
              #print(f' src range = {srcRange}')
              #print(f' dst range = {dstRange}')
              srcIndex = srcRange[0]
              dstIndex = dstRange[0]
              while srcIndex < srcRange[1]:
                base = seq[srcIndex]
                #print(srcIndex,dstIndex,base,len(seq),len(ref))
                # alignment should never be longer than ref
                if dstIndex > len(alignSeq):
                  print(f'error, dstIndex {dstIndex} longer than ref {len(ref)}')
                else:
                  alignSeq[dstIndex] = base;
                  srcIndex += 1
                  dstIndex += 1
            #print(ref)
            #t = ''.join(alignSeq)
            #print(t)
            for dstIndex,base in enumerate(alignSeq):
                baseCount[dstIndex][base] += seqCount
            #
            # re-calc baseErrors...no max error
            #
            baseErrors = compSeq(ref,alignSeq,len(ref))
            #print(f'ref {ref[:160]}')
            #t = ''.join(alignSeq)
            #print(f'als {t[:160]}')
            #for i,br,bs in baseErrors:
            #  print(i,br,bs)
            #quit()


            if baseErrors == None:
              print(f'error in re-calculatre base errors after indel alignment')
            
            for bi,bRef,bSeq in baseErrors:
              perBaseError[bi][bSeq] += seqCount
      else:
        #
        # non-indel
        #
        # base errors
        #
        for bi,bRef,bSeq in baseErrors:
          perBaseError[bi][bSeq] += seqCount
        #
        # all base stats
        #
        for seqIndex in range(0,len(seq)):
          base = seq[seqIndex]
          baseCount[seqIndex][base] += seqCount
          
          

  print(fIn,totalSeqCount)

  baseReturn = [(-1,0.0) for x in range(0,len(ref))]
  fh.write(f'#total sequence count = {totalSeqCount}  ref count = {refCount}  match = {refCount/totalSeqCount:2.8f}\n')
  fh2.write(f'#totalSequenceCount,{totalSeqCount}\n')
  fh2.write(f'index,protoLoc,reference,A,C,G,T,-\n')
  fh.write(f'#protospacer location = {protoSpacerLoc}\n')
  #
  # summarie per-base errors for out.csv
  #
  for index in range(0,len(perBaseError)):
    errorD = perBaseError[index]
    if len(errorD) > 0:
      fh.write(f'#--------------------------\n')
      fh.write(f'#[{index}]\n')
      fh.write(f'#ref = {ref[index]}\n')
      baseErrors = 0
      for bSeq,bCount in errorD.items():
        fh.write(f'#    {bSeq}    {bCount}\n')
        baseErrors += bCount
      rError = baseErrors/totalSeqCount
      fh.write(f'#{index} total base error rate = {rError:2.5f}\n')
      fh.write(f'{index},{rError:2.5f}\n')
      baseReturn[index] = (index,rError)
    else:
      fh.write(f'{index},0.0\n')
      baseReturn[index] = (index,0.0)
  #
  #  r_out.csv just base composition
  #
  for index in range(0,len(baseCount)):
    baseDict = baseCount[index]

    b_a = baseDict['A']
    b_c = baseDict['C']
    b_g = baseDict['G']
    b_t = baseDict['T']
    b__ = baseDict['-']

    fh2.write(f'{index},{protoSpacerLoc},{ref[index]},{b_a},{b_c},{b_g},{b_t},{b__}\n')
  fh2.close()
  fh.close()



  #
  # summarize indels
  #
  fhIndel.write('location,type,length,count,fraction\n')
  # sort by key(location)
  sd = sorted(indelDict.items(),key=lambda x:x[0],reverse=False)
  for loc,indelList in sd:
    tmpD = defaultdict(int)
    for iType,count in indelList:
      tmpD[iType] += count
    sortTmp = sorted(tmpD.items(),key=lambda x:x[0],reverse=False)
    for iType,iCount in sortTmp: 
      iText = 'insertion'
      if iType > 0:
        iText = 'deletion'
      fhIndel.write(f'{loc},{iText},{abs(iType)},{iCount},{iCount/totalSeqCount:2.5f}\n')

  fhIndel.close()
  #
  # return value is summary list 
  #
  return baseReturn
#------------------------------------------------------------------------------------------
#
#    Main code point
#
#
#------------------------------------------------------------------------------------------
#
# gather all files
#
tsv_files = [f for f in os.listdir(inputDir) if f.endswith('.tsv')]
print(f'process {len(tsv_files)} files from {inputDir}')
#
# process each file
#
mainList = []
for f in tsv_files:
  fin = inputDir + '/' + f
  fout = inputDir + '/' + f[:-4] + ".csv"
  fout2 = inputDir + '/r_' + f[:-4] + ".csv"
  fIndel = inputDir + '/i_' + f[:-4] + ".csv"
  #print(f'Process file {fin} to {fout}')
  errorList = processFile(fin,fout,fout2,fIndel)
  #
  # errorList is base error rate for each position
  #
  #
  # make a dictionary for sorting while keeping location
  #
  # also make a list of just the error rates
  #
  #
  errorLocDict = defaultdict(float)
  numericList = []
  for i_l,e_l in errorList:
    errorLocDict[i_l] = e_l
    numericList.append(e_l)
  #
  #
  #
  meanV = sum(numericList) / len(numericList)
  #
  # find top error locations
  #
  print(f'Error summary for {f}')
  sd = sorted(errorLocDict.items(),key=lambda x:x[1],reverse=True)
  #
  # print up to first 10
  #
  for index,(loc,e) in enumerate(sd):
    print(f'[{loc}]   e =  {e:2.5f}')
    if index >= 10:
      break
  #
  #if len(sd) >= 4:
  #  mainList.append((fin,meanV,[sd[0],sd[1],sd[2],sd[3]]))
  #else:
  #  print(f'{f} errorList = {len(errorList)}   {len(errorLocDict)}    {len(sd)}')
  #  for a in errorList:
  #    print(a)
  #for f,m,n in mainList:
  #s = ''
  #for entry in n:
  #  s = f'{s} ({entry[0]} {entry[1]:2.5f})'
  #print(f'{f}  {m:2.5f}   {s}')


