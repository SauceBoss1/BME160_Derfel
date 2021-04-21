#!/usr/bin/env python3
# Name: Derfel Terciano (dtercian)
# Group Members: None

import sequenceAnalysis
import sys
#from importlib import reload
#reload(sequenceAnalysis)

class GenomeAnalyzer:
    def __init__(self,fileName=''):
        self.myReader=sequenceAnalysis.FastAreader(fileName)
        self.myNuc=sequenceAnalysis.NucParams()
        for head,seq in self.myReader.readFasta():
            head = head
            self.myNuc.addSequence(seq)

    def analyzeGenome(self):
        numBases = self.myNuc.nucCount()/1000000
        contentGC = ((self.myNuc.nucComposition()['C'] + self.myNuc.nucComposition()['G']) / self.myNuc.nucCount()) *100
        print(f'sequence length = {numBases:0.2f} Mb\n')
        print(f'GC content = {contentGC:0.1f}%\n')


        for codon, aa in sorted(self.myNuc.rnaCodonTable.items(), key=lambda a:(a[1],a[0])):
            codonFreq = (self.myNuc.codonComposition()[codon] / self.myNuc.aaComposition()[aa]) * 100
            print (f'{codon} : {aa} {codonFreq:5.1f} ({self.myNuc.codonComposition()[codon]:6d})')



def main ():
    genome = GenomeAnalyzer()
    genome.analyzeGenome()
main()
