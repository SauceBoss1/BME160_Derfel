#!/usr/bin/env python3
# Name: Derfel Terciano (dtercian)
# Group Members: None

import sequenceAnalysis


class GenomeComparator:
    def __init__(self,fileName1='',fileName2=''):
        self.fileName1 = fileName1
        self.fileName2 = fileName2
        self.genome1 = self.fileInitiator(fileName1)
        self.genome2 = self.fileInitiator(fileName2)

    def fileInitiator(self,fileName):
        inputFile = sequenceAnalysis.FastAreader(fileName)
        genome = sequenceAnalysis.NucParams()
        for header, seq in inputFile.readFasta():
            genome.addSequence(seq)

        return genome
    def genomeComparator(self):
        gcContent1 = ((self.genome1.nucComposition()['C'] + self.genome1.nucComposition()['G']) / self.genome1.nucCount()) *100
        gcContent2 = ((self.genome2.nucComposition()['C'] + self.genome2.nucComposition()['G']) / self.genome2.nucCount()) *100
        
        print (f'Genomes we are comparing: => {self.fileName1} | {self.fileName2}\n\n')

        print(f'GC content comparison of: => {self.fileName1} | {self.fileName2}')
        if gcContent1 > gcContent2:
            print(f'{gcContent1:2.1f}% | > | {gcContent2:2.1f}%\n')
        elif gcContent1 < gcContent2:
            print(f'{gcContent1:2.1f}% | < | {gcContent2:2.1f}%\n')
        else:
            print(f'{gcContent1:2.1f}% | = | {gcContent2:2.1f}%\n')

        print(f'\nNucleotide comparison of: => {self.fileName1} | {self.fileName2}')
        for aa in sorted(self.genome1.nucComposition()):
            aa1 = self.genome1.nucComposition()[aa]
            aa2 = self.genome2.nucComposition()[aa]

            if aa1 > aa2:
                print(f'{aa} : {aa1:8d} | > | {aa2:8d}')
            elif aa1 < aa2:
                print(f'{aa} : {aa1:8d} | < | {aa2:8d}')
            else:
                print(f'{aa} : {aa1:8d} | = | {aa2:8d}')

        print(f'\nCodon usage comparison of: => {self.fileName1} | {self.fileName2}')
        for codon, aa in sorted(self.genome1.rnaCodonTable.items(), key=lambda a:(a[1],a[0])):
            genome1CodonCount = self.genome1.codonComposition()[codon]
            genome2CodonCount = self.genome2.codonComposition()[codon]
            codonFreq1 = (genome1CodonCount/self.genome1.aaComposition()[aa])*100
            codonFreq2 = (genome2CodonCount/self.genome2.aaComposition()[aa])*100

            if codonFreq1 > codonFreq2:
                print(f'{codon} : {aa} => {codonFreq1:5.1f} ({genome1CodonCount:6d}) | > | {codonFreq2:5.1f} ({genome2CodonCount:6d})')
            elif codonFreq1 < codonFreq2:
                print(f'{codon} : {aa} => {codonFreq1:5.1f} ({genome1CodonCount:6d}) | < | {codonFreq2:5.1f} ({genome2CodonCount:6d})')
            else:
                print(f'{codon} : {aa} => {codonFreq1:5.1f} ({genome1CodonCount:6d}) | = | {codonFreq2:5.1f} ({genome2CodonCount:6d})')
def main():
    genomes = GenomeComparator('testGenome.fa','haloVolc1_1-genes.fa')
    #genomes = GenomeComparator() #i dont know how to input 2 files
    genomes.genomeComparator()
main()