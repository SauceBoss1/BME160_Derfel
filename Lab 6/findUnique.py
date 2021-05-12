#!/usr/bin/env python3
# Name: Derfel Terciano (dtercian)
# Group Members: List full names (CATS usernames) or “None”
'''NOTES: code runs in 10-11 seconds'''
# Design for findUnique

# 1) input all sequences into class
# 2) find power set for all of the tRNA sequences and store it into the class
# 3) have a function that takes in a header, and do all the set comparison based on the header and sequence
# 4) have the final output be created based on calling the header of each sequence           
import sys
import time
start_time = time.time()
from sequenceAnalysis import FastAreader #type: ignore
class FindUnique:
    def __init__(self, header='', seq=''):
        self.tRnaInfo ={} #the form of the dictionary is in the following => header: [seq, powerset(seq)] (this is good for comparing current sets)
        self.headers = [] #store headers
        self.allTrnaPset = [] #store all pSets from all 22 sequences (useful for unioning without using for loops)
        self.addSeq(seq, header)    

    def addSeq(self, seq, header):
        header = header.replace(' ','')
        self.seq= seq.upper().replace('-','').replace('_','').replace('.','')
        if header != '':
            self.headers.append(header)
            seqPset = self.powerset(self.seq)
            self.tRnaInfo[header] = [self.seq, seqPset]
            self.allTrnaPset.append(seqPset)

    def powerset(self,inSeq):
        pSet= set()
        for i in range(len(inSeq)):
            for j in range(i+1,len(inSeq)+1):
                pSet.add(inSeq[i:j])
        return pSet

    def uniqueFinder(self, curHeader):
        allSets = self.allTrnaPset.copy()
        allSets.remove(self.tRnaInfo[curHeader][1])
        return set(self.tRnaInfo[curHeader][1]) - set().union(*allSets)

    def essential(self, curHeader):
        nonEssentials = set()
        uniqueSet = self.uniqueFinder(curHeader)
        compareSet = uniqueSet.copy()
        for subset in uniqueSet:
            compareSet.remove(subset)
            for compareSub in compareSet:
                if compareSub in subset:
                    nonEssentials.add(subset)
            compareSet.add(subset)
        return set(uniqueSet) - nonEssentials

    def findAllOccurances(self, subSeq, mainSeq):
        posFound = []
        index = mainSeq.find(subSeq)
        while index != -1:
            posFound.append(index)
            index = mainSeq.find(subSeq, index+1)
        return posFound
    
    def findUnique(self): #problem: some essentials have multiple positions on sequence
        sys.stdout.reconfigure(encoding='utf-8')
        with sys.stdout as outFile:
            for currentHeader in sorted(self.headers, key = lambda a:a[5:8]): #adjust order of sequences here
                print(f'{currentHeader}\n{self.tRnaInfo[currentHeader][0]}')
                essentials = self.essential(currentHeader)
                essentialsIndex = {}
                for essentialSubs in essentials:
                    for indexes in self.findAllOccurances(essentialSubs, self.tRnaInfo[currentHeader][0]):
                        essentialsIndex[indexes] = essentialSubs
                for index, seq in sorted(essentialsIndex.items(), key = lambda a:a[0]):
                    print(f'{index * "."}{seq}')
            #print('run time (in seconds) : ', time.time()-start_time) #uncomment this if you want to see the printed run time of the program at the end of output file

########################################################################
# Main
# Here is the main program
# 
########################################################################

def main(inCL=None):
    inFile = FastAreader()
    genes = FindUnique()
    for header,seq in inFile.readFasta():
        genes.addSeq(seq, header)
    genes.findUnique()
if __name__ == "__main__":
    main()  