#!/usr/bin/env python3
# Name: Your full name (CATS account username)
# Group Members: List full names (CATS usernames) or “None”

import sys
from sequenceAnalysis import FastAreader #type: ignore

class FindUnique:
    def __init__(self, seq=''):
        self.tRnaInfo ={}
        self.addSeq(seq)

    def addSeq(self, seq, header):
        seq= seq.upper().replace('-','').replace('_','').replace('.','')
        self.tRnaInfo[header] = [seq, self.powerset(self.seq)]
        
    def powerset(self,inSeq):
        pSet= set()
        for i in range(len(inSeq)):
            for j in range(i+1,len(inSeq)+1):
                pSet.add(inSeq[i:j])
        return pSet
'''
Design for findUnique

1) input all sequences into class
2) find power set for all of the tRNA sequences and store it into the class
3) have a function that takes in a header, and do all the set comparison based on the header and sequence
4) have the final output be created based on calling the header of each sequence
'''
            
########################################################################
# Main
# Here is the main program
# 
########################################################################

def main(inCL=None):
    ''' '''
    pass
if __name__ == "__main__":
    main()  
