import random
import collections
import dnatools.getools as gu
import numpy as np
from datetime import datetime
import time
import math
from dnatools.dnautils import BWT



dnaAlphabet = 'ACGT'
rnaAlphabet = 'ACGU'

dnaSeqType = 'dna'
rnaSeqType = 'rna'


class Aligner:
    def __init__(self):
        self.pattern = ''
        self.text    = ''
        self.matchPositions = []
        self.matchComputations = 0
        self.minComputations = 0
        self.maxComputations = 0
        self.preProcessExecutionTime = 0
        self.T = ''
  
    def processT(self,P,T,revCompl,matchAr,debug=False):
        
        return -1
    
    def splitP(self,P,m):
        
        splits = []
        p = P.seq
        
        splitSpan = math.floor(P.seqLen() / (m+1))
        for i in range(0,m):
            split = p[i*splitSpan:(i+1)*splitSpan]
            splits.append((split,i*splitSpan,len(p)- i*splitSpan - len(split)))
            
        split =  p[m*splitSpan:]
        splits.append((split,m*splitSpan,len(p)-m*splitSpan - len(split)))
        return splits    
            
     
    def hamDist(self,p,t):
        if (len(p) != len(t)):
           return len(p) + len(t)
    
        hamDist = 0
        
        for i,pCh in enumerate(p):
            if (pCh == t[i]):
               pass
            else:
               hamDist +=1 
        return hamDist    
 
    def showMatches(self,l=100):
        # can be called after a find to show the actual matching strings
        matches = []
        for entry in self.matchPositions:
            matches.append(self.T.seq[entry:entry+l])
        return matches                   
            
            

    def numMismatches(self,st,P,T,extraDict):
        return self.hamDist(P.seq,T.seq[st:st+P.seqLen()])

    def findMatches(self,P,T,revCompl=False,maxMismatches = 0,debug=False):
       # print ('in base findmatch')
        matchAr = []
        self.matchComputations = 0
        
        startFind = time.time()
        
       # print ('finding matches in T: ',T.seq[:200])
        if (maxMismatches == 0):
           self.matchComputations = self.processT(P,T,revCompl,matchAr,debug)
        else:
           splits = self.splitP(P,maxMismatches)
           #print ('splits:',splits)
           for split,splitPos,revSplitPos in splits:
               potMatchAr = [] 
               #print ('split: ',split)
               partP = P.copy(split)
               #print ('part p: ',partP.seq) 
               #print ('P copy split: ',P.copy(split).seq)
               self.matchComputations+= self.processT(P.copy(split),T,revCompl,potMatchAr,debug)
               #print ('pot match ar: ',potMatchAr)
               for entry in potMatchAr:
                   
                   st = entry - splitPos
                    #print ('entry: ',entry, ' splitpos: ',splitPos, ' st: ',st)
                   #print ('potential: ',entry, 'P: ',P.seq, ' T: ',T.seq[st:st+P.seqLen()])
                   #if (self.hamDist(P.seq,T.seq[st:st+P.seqLen()]) <= maxMismatches):
                   if (self.numMismatches(st,P,T)  <= maxMismatches):
                           matchAr.append(st)
                   if (revCompl):
                      st = entry - revSplitPos
                      pR = P.reverseComplement().seq
                      #print ('pR: ',pR, 't: ',T.seq[st:st+P.seqLen()], ' st: ',st ) 
                      #if (self.hamDist(P.reverseComplement().seq,T.seq[st:st+P.seqLen()]) <= maxMismatches):
                      if self.numMismatches(st,P.reverseComplement(),T) <= maxMismatches:
                            matchAr.append(st)
                
               
           
        self.T = T        
 
        matchAr = list(set(matchAr))
        self.matchPositions = matchAr
        
        endFind = time.time()
        self.findExecutionTime = endFind - startFind
        
        if (self.matchComputations == -1):
           return -1 #Not implemented
        else:
           return len(matchAr)  
    
    
class BMAligner(Aligner):

      def __init__(self):
         super(BMAligner, self).__init__()
   
      def preProcess(self,P,debug=False):
          #print ('pre st: ',str(datetime.now()))
          startPreProcess = time.time()  
          badCharDists = P.allDistsToNearestBaseOnLeft()
          goodSuffixDists = P.allSufficesDistsToNearestOnLeft() 
          if (debug):
             print (badCharDists)
             print (goodSuffixDists)
          #print ('pre end: ',str(datetime.now()))  
          self.badCharDists = badCharDists
          self.goodSuffixDists = goodSuffixDists 
          endPreProcess = time.time()
          self.preProcessExecutionTime = endPreProcess - startPreProcess
        
          return badCharDists,goodSuffixDists

      def onePass(self,P,T,matchAr,debug):
        
          comps = 0 
        
          badCharDists,goodSuffixDists = self.preProcess(P,debug)
            
          currStart = 0
          currEnd = currStart + P.seqLen() - 1
          done = False
          p = P.seq
          t = T.seq  
        
          #print ('main st: ',str(datetime.now()))
          while not done:
              if (debug):
                 print ('curr pos: ',currStart, 'txt: ',t[currStart:currEnd+1])  
              currPrefix = ''
              match = True
              for pos in range(currEnd,currStart-1,-1):
                  comps  +=1
                  pPos = len(p) - 1 - (currEnd - pos)  
                  currPrefix = t[pos] + currPrefix
                  if currPrefix == p[pPos:]:
                     pass
                  else:
                     match = False
                     break
              if match:
                 if (debug):
                    print ('match at pos: ',currStart)
                 matchAr.append(currStart)   
                 currStart +=1
                 currEnd   +=1 
              else:
                 badChar = t[pos]
                 goodSuffix = p[pPos+1:] 
  
                 if (badChar in badCharDists[pPos][1]):   
                    skipBad = badCharDists[pPos][1][badChar]
                 else:
                    skipBad = pPos + 1
                 skipGood =  goodSuffixDists[goodSuffix]
                 if (debug):
                     print ('bad char: ',badChar, 'move: ',skipBad,'  good suff:',goodSuffix,'move: ',skipGood)
                 #print ('skip bad: ',skipBad,'skip good: ',skipGood)   
                 currStart +=max(skipGood,skipBad)
                 currEnd +=max(skipGood,skipBad)   
                
              if (currEnd >= len(t)):
                  done = True   
          #print ('main end: ',str(datetime.now()))            
          return comps
        
      def processT(self,P,T,revCompl,matchAr,debug=False):
            
          comps = self.onePass(P,T,matchAr,debug)
        
          if (revCompl):
             PRev = P.reverseComplement()
             comps += self.onePass(PRev,T,matchAr,debug) 
                

          return comps
    

                
class IndexAligner(Aligner):

      def __init__(self,k):
         super(IndexAligner, self).__init__()
         self.k = k
         self.T = None
   
      def preProcess(self,T,debug=False):
          print ('Preprocessing')
          startPreProcess = time.time()  
          t = T.seq
          self.iDict = {}
          for i in range(len(t) - self.k + 1):
              kmer = t[i:i+self.k]
              if (kmer in self.iDict):
                 self.iDict[kmer].append(i)
              else:
                 self.iDict[kmer] = [i]
                    
          self.T = T
        
          endPreProcess = time.time()
          self.preProcessExecutionTime = endPreProcess - startPreProcess
          print ('End Preprocessing')
 
        
      def verify(self,P,T,potPos,debug=False):
         p = P.seq
         t = T.seq
        
         if p == t[potPos:potPos+len(p)]:
            return True
        
         return False


      def lookUp(self,P,T,matchAr,debug=False):
        
          p = P.seq
          t = T.seq
        
          comps = 0
            
          kmer = p[:self.k]
          
          if (kmer in self.iDict):
            
             for potEntry in self.iDict[kmer]:
                 comps += 1   
                 if self.verify(P,T,potEntry,debug):
                    matchAr.append(potEntry)
          return comps
            
      def  processT(self,P,T,revCompl,matchAr,debug=False):
         
          if (self.k > P.seqLen()):
              print ('Index too large for pattern. k = ',self.k,' Please decrease k to at most: ',P.seqLen())
              return -1
        
          if self.T and (self.T.seq ==  T.seq):
              pass #already preprocessed
          else:
              self.preProcess(T,debug)
        
        
          comps = self.lookUp(P,T,matchAr,debug)
          
          if (revCompl):
             PRev = P.reverseComplement()
             comps +=self.lookUp(PRev,T,matchAr,debug)
                
          return comps      
    
class BWTAligner(Aligner):
     def __init__(self):
         super(BWTAligner, self).__init__()
         self.T = None
 
     def numMismatches(self,st,P,T):

        #inverse=False
        #if ('inverse' in extraDict):
        #    inverse=extraDict['inverse']
            
        #if inverse:
        #   return self.hamDist(P.seq,T.seq[st:st+P.seqLen()]) 
        #else:
        #   return self.hamDist(P.seq,T.seq[st:st+P.seqLen()])
        if ((st + P.seqLen()) >= len(self.bwtObj.t)):
            return  99999
        else:    
            return self.hamDist(P.seq,self.bwtObj.t[st:st+P.seqLen()]) 
           
     def preProcess(self,T,debug=False):
          print ('Preprocessing')

          startPreProcess = time.time()  
       
          self.bwtObj = BWT(debug=debug)
          #self.bwtObj = BWT(T.seq + '$',debug=debug)
         
          if (self.inverse):
              self.bwtObj.loadBWT(T.seq,psd=self.psd,psa=self.psa,reconstruct=self.reconstruct)
                  
          else:
              self.bwtObj.loadT(T.seq + '$')
            
          self.T = T
            
          endPreProcess = time.time()
          self.preProcessExecutionTime = endPreProcess - startPreProcess
          print ('End Preprocessing')
   

     def findMatches(self, *args, **kwargs):
        inverseParams = {}
        argopt2 = kwargs.pop('inverseParams', None)
        # remove the extra arg so the base class doesn't complain. 
        #del kwargs['argopt2']
        if (argopt2 is None):
            pass
        else:
            inverseParams = argopt2
        #print (argopt2)
        self.inverse=False
        self.reconstruct = False
        self.psa = None
        self.psd = None
        if (inverseParams is None):
           pass
        else: 
            #print ('inv params: ',inverseParams)
            if ('inverse' in inverseParams):
               self.inverse = inverseParams['inverse']

            if ('reconstruct' in inverseParams):
               self.reconstruct = inverseParams['reconstruct']


            if ('psa' in inverseParams):
               self.psa = inverseParams['psa']


            if ('psd' in inverseParams):
               self.psd = inverseParams['psd']  
         
        res = super().findMatches(*args, **kwargs) 
        return res   
            
     def  processT(self,P,T,revCompl,matchAr,debug=False):
         
          if self.T and (self.T.seq ==  T.seq):
              pass #already preprocessed
          else:
              self.preProcess(T,debug)
   
          _,_,_,_,_,matchPositions = self.bwtObj.find(P.seq)
    
          #print ('matchposiitions: ',matchPositions)
    
          matchPositionsRev = []
          if (revCompl):
            _,_,_,_,_,matchPositionsRev = self.bwtObj.find(P.reverseComplement().seq)
            #print ('matchpositionsrev: ',matchPositionsRev)
    
          matchAr.extend(matchPositions)
          matchAr.extend(matchPositionsRev)
          return 0
          
   

class SuffixArrayAligner(Aligner):

      def __init__(self,maxSuffArrSortLen=101):
         super(SuffixArrayAligner, self).__init__()
         self.T = None
         self.maxSuffArrSortLen = 51 #maximum length of suffix of T used in sorting
         self.matchSuffArrayIndices = []
            
   
      def preProcess(self,T,debug=False):
          print ('Preprocessing')
          startPreProcess = time.time()  
          t = T.seq
          self.suffArray = [i for i in range(len(t))]
          
          self.suffArray = sorted(self.suffArray, key=lambda x: t[x:x+self.maxSuffArrSortLen] + '$')
          
          self.T = T
            
          endPreProcess = time.time()
          self.preProcessExecutionTime = endPreProcess - startPreProcess
          print ('End Preprocessing')
   
        
 
 
      def intToSuff(self,i,T):
          t = T.seq
          return t[i:] + '$'
          
      def posToSuff(self,pos,T):
          if (pos >= len(self.suffArray)):
             print('Error: Pos: ',pos, ' Suf Arr len: ',len(self.suffArray))   
          return self.intToSuff(self.suffArray[pos],T)

      def  findInRegion(self,P,T,stReg,endReg):
           
           #print ('stReg: ',stReg,'endReg: ',endReg) 
           p = P.seq
           t = T.seq
            
           if (stReg == endReg):
               if (p == self.posToSuff(stReg,T)[:len(p)]):
                   return self.suffArray[stReg],stReg
               else:
                   return -1,-1
            
           else:
               midPoint = stReg + int((endReg - stReg)/2)
               #print ('midpoint:',midPoint) 
               midStr =  self.posToSuff(midPoint,T)
               if (p == midStr[:len(p)]):
                  return self.suffArray[midPoint],midPoint
               elif (p < midStr):
                    return self.findInRegion(P,T,stReg,midPoint)
               else:
                    return self.findInRegion(P,T,midPoint+1,endReg)
                 
      def  processT(self,P,T,revCompl,matchAr,debug=False):
         
          if self.T and (self.T.seq ==  T.seq):
              pass #already preprocessed
          else:
              self.preProcess(T,debug)
        
        
          tPos,suffArrPos = self.findInRegion(P,T,0,T.seqLen()-1)
          #print ('result: ',tPos)  
        
          if (tPos == -1):
             pass
          else:
             self.spread(P,T,tPos,suffArrPos,matchAr)
            
            
          if (revCompl):
             PRev = P.reverseComplement()
             tPos,suffArrPos = self.findInRegion(PRev,T,0,T.seqLen()-1)
             if (tPos == -1):
                pass
             else:
                self.spread(PRev,T,tPos,suffArrPos,matchAr)  
          #comps = self.lookUp(P,T,matchAr,debug)
          
          #if (revCompl):
          #   PRev = P.reverseComplement()
          #   comps +=self.lookUp(PRev,T,matchAr,debug)
                
          return 0      
 
      def spread(self,P,T,hitPos,hitSuffArrPos,matchAr):
          p = P.seq
          currPos = hitSuffArrPos
            
          while (currPos >= 0) and (p == self.posToSuff(currPos,T)[:len(p)]):
               currPos-=1
                
          currPos +=1
          while (currPos < T.seqLen()) and (p == self.posToSuff(currPos,T)[:len(p)]):
                matchAr.append(self.suffArray[currPos])
                self.matchSuffArrayIndices.append(currPos)
                currPos+=1
                
 
        
                
                
            
        
       
    
    
class NaiveExactAligner(Aligner):
    
      def __init__(self):
         super(NaiveExactAligner, self).__init__()
        

      def checkExactMatch(self,P,T):
          matchFound = True
          comps = 0
         
          for i in range(0,len(P)):
             comps +=1
             if P[i] != T[i]:
                matchFound = False
                break
          return matchFound,comps
    
      def processT(self,P,T,revCompl,matchAr,debug=False):
           
         comps = 0
            
         if (revCompl):
              rCompP = P.reverseComplement()
     
         for i in range(0,T.seqLen()- P.seqLen()+1):
              
              compF = 0
              compR = 0
                
              matchFound,compF = self.checkExactMatch(P.seq,T.seq[i:i+P.seqLen()])
 
              if matchFound:
                 #print ('pos match found: ',i)
                 pass
            
 
              if (not matchFound) and (revCompl):
                 
                 matchFound,compR =   self.checkExactMatch(rCompP.seq,T.seq[i:i+P.seqLen()])
                 if matchFound:
                    pass    
                    #print ('neg match found: ',i)
                    #print ('checked neg p: ',rCompP.seq, ' T: ',T.seq[i:i+P.seqLen()])
                    
 
              if matchFound:
                 #print ('match found: ',i)
                 matchAr.append(i)
                    
              comps += compF
              comps += compR
         
         self.maxComputations = P.seqLen() * (T.seqLen() - P.seqLen() + 1)
         self.minComputations = T.seqLen() - P.seqLen() + 1
        
         return comps      
        
            
class Assembler():
    def __init__(self,reads,refGenome,aligner,maxMismatches=0):
        self.reads = reads
        self.refGenome = refGenome
        self.aligner = aligner
        self.maxMismatches = maxMismatches
        self.assembledSeqArr = ['-' for ch in self.refGenome.seq]
       
        

    def addToAssembled(self,matchPos,read):
        for i in range(0, read.seqLen()):
            self.assembledSeqArr[matchPos + i] = read.seq[i]
        
        
    def assemble(self):
        
        for read in self.reads:
            numMatches = self.aligner.findMatches(read,self.refGenome,maxMismatches=self.maxMismatches)
            if (numMatches == 1):
                #print ('single match: ',numMatches)
                self.addToAssembled(self.aligner.matchPositions[0],read)
            elif (numMatches == 0):
                print ('No match: ',numMatches)
            else:
                print('Multiple matches!', numMatches)
        
        
        
        return ''.join(self.assembledSeqArr)
                
        
if __name__ == '__main__':
    p = 'ACG'
    t = 'CGATGCTGCTGAGCAGCTACGGATGCTGGCTAGCTGATGC'
    al = NaiveExactAligner()
    print(al.findMatches(p,t))
