#!/usr/bin/env python3
# Name: Your full name (CATS account username)
# Group Members: List full names (CATS usernames) or “None”

import sys
from sequenceAnalysis import FastAreader #type: ignore

class FindUnique:
    def __init__(self, seq=''):
        self.seq = seq
        self.tRnaPowersets =[]
        self.addSeq(seq)

    def addSeq(self,seq):
        self.seq= seq.upper().replace('-','').replace('_','').replace('.','')
        

    def powerset(self):
        pSet= set()
        for i in range(len(self.seq)):
            for j in range(i+1,len(self.seq)+1):
                pSet.add(self.seq[i:j])
        return pSet
    
            
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
