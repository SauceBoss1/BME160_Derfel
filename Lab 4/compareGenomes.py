#!/usr/bin/env python3
# Name: Derfel Terciano (dtercian)
# Group Members: None
'''
compareGenomes.py

Written by: Derfel Terciano

This program compares two genomes from two given FastA files and returns a formatted output of its results.

The results are:
    -GC content comparison of both files
    -Nucleotide content comparison of both genomes
    -A comparison of codon usage of both genomes 

Class(es) used:
    -GenomeComparator: => handles the formatted output of the comparison of both files

Assumptions:
    -ONLY two files are given no more and no less.
'''
import sequenceAnalysis


class GenomeComparator:
    '''
    Formats an output of the comparison between 2 genomes (in FastA format)

    Compares:
        -GC Conent
        -Nucleotide count
        -Codon usage
    
    public method(s):
        -genomeComparator: => prints a formatted output of the above comparisons

    private method(s):
        -__fileInitiator(fileName): => helps convert a FastA file to a readable genome by returning a NucParams object of the file
        -__init__(fileName1,fileName2): => takes in two files
    '''
    def __init__(self,fileName1='',fileName2=''):
        '''Takes in two FastA file and extracts a readable genome'''
        self.fileName1 = fileName1
        self.fileName2 = fileName2
        self.genome1 = self.__fileInitiator(fileName1) #extracts the NucParams object of the first genome
        self.genome2 = self.__fileInitiator(fileName2) #extracts the NucParams object of the second genome

    def __fileInitiator(self,fileName):
        '''
        Using the FastA object of sequenceAnalysis, this function will extract all genomic data
        and stores it into the NucParams class with all necessary calculations
        '''
        inputFile = sequenceAnalysis.FastAreader(fileName) #read file input with FastAreader
        genome = sequenceAnalysis.NucParams()
        for header, seq in inputFile.readFasta(): #get only the genomic data of the FastA file
            genome.addSequence(seq)

        return genome #return the NucParams object of the inputted file

    def genomeComparator(self):
        '''Prints a formatted output of the genomic comparison of the 2 genomes'''

        #the following if elif statment makes sure that there are valid codons to compare if not, terminate program.
        if self.genome1.nucCount()==0 or (sum(self.genome1.nucComposition().values())==0):
            print(f'WARNING NO VALID CODONS OR NUCLEOTIDES FOUND IN {self.fileName1}. NOTHING TO COMPARE TO {self.fileName2}')
            exit() #terminate program
        elif self.genome2.nucCount()==0 or (sum(self.genome2.nucComposition().values())==0):
            print(f'WARNING NO VALID CODONS OR NUCLEOTIDES FOUND IN {self.fileName2}. NOTHING TO COMPARE TO {self.fileName1}')
            exit() #terminate program

    
        #the 2 lines below finds the GC content of both genomes
        gcContent1 = ((self.genome1.nucComposition()['C'] + self.genome1.nucComposition()['G']) / self.genome1.nucCount()) *100
        gcContent2 = ((self.genome2.nucComposition()['C'] + self.genome2.nucComposition()['G']) / self.genome2.nucCount()) *100
        
        print (f'Genomes we are comparing: => {self.fileName1} | {self.fileName2}\n\n') #announce to the users the files they have given the program

        print(f'GC content comparison of: => {self.fileName1} | {self.fileName2}') #announce the GC content comparison is taking place
        if gcContent1 > gcContent2: #shows that genome1 has a greater GC content
            print(f'{gcContent1:2.1f}% | > | {gcContent2:2.1f}%\n')
        elif gcContent1 < gcContent2: #shows that genome2 has a greater GC content
            print(f'{gcContent1:2.1f}% | < | {gcContent2:2.1f}%\n')
        else: #shows that both genomes have the same GC content
            print(f'{gcContent1:2.1f}% | = | {gcContent2:2.1f}%\n')

        print(f'\nNucleotide comparison of: => {self.fileName1} | {self.fileName2}') #announce the nucleotide comparison is taking place
        for aa in sorted(self.genome1.nucComposition()): #iterate through the nucComp dictionary of both genomes (Note: it doesn't matter which genome to iterate, both are NucParams)
            aa1 = self.genome1.nucComposition()[aa]
            aa2 = self.genome2.nucComposition()[aa]

            if aa1 > aa2: #shows that genome1 has a greater count of the current AA
                print(f'{aa} : {aa1:8d} | > | {aa2:8d}')
            elif aa1 < aa2: #shows that genome2 has a greater count of the current AA
                print(f'{aa} : {aa1:8d} | < | {aa2:8d}')
            else: #shows that both have the same AA count of the current AA
                print(f'{aa} : {aa1:8d} | = | {aa2:8d}')

        print(f'\nCodon usage comparison of: => {self.fileName1} | {self.fileName2}') #announce that codon usage comparison is taking place
        for codon, aa in sorted(self.genome1.rnaCodonTable.items(), key=lambda a:(a[1],a[0])): #iterate through the sorted codonComp dictionary of both genomes (Note: it doesn't matter which genome to iterate, both are NucParams)
            currentAAcomp1 = self.genome1.aaComposition()[aa]
            currentAAcomp2 = self.genome2.aaComposition()[aa]

            #the following code prevents a division by zero error if the current AA has no counts
            if self.genome1.aaComposition()[aa] == 0:
                currentAAcomp1 = 1
            elif self.genome2.aaComposition()[aa] == 0:
                currentAAcomp2 = 1

            genome1CodonCount = self.genome1.codonComposition()[codon] 
            genome2CodonCount = self.genome2.codonComposition()[codon]
            codonFreq1 = (genome1CodonCount/currentAAcomp1)*100 #calculates codon freq of genome1
            codonFreq2 = (genome2CodonCount/currentAAcomp2)*100 #calculates codon freq of genome2

            if codonFreq1 > codonFreq2: #shows that codon freq of current freq is greater in genome1
                print(f'{codon} : {aa} => {codonFreq1:5.1f} ({genome1CodonCount:6d}) | > | {codonFreq2:5.1f} ({genome2CodonCount:6d})')
            elif codonFreq1 < codonFreq2: #shows that codon freq of current freq is greater in genome2
                print(f'{codon} : {aa} => {codonFreq1:5.1f} ({genome1CodonCount:6d}) | < | {codonFreq2:5.1f} ({genome2CodonCount:6d})')
            else: #shows that codon freq of current freq is the same in both genomes
                print(f'{codon} : {aa} => {codonFreq1:5.1f} ({genome1CodonCount:6d}) | = | {codonFreq2:5.1f} ({genome2CodonCount:6d})')
def main():
    '''executes main output sequence of program'''
    genomes = GenomeComparator('testGenome.fa','haloVolc1_1-genes.fa') #give the 2 files here
    #genomes = GenomeComparator() #i dont know how to input 2 files through command line
    genomes.genomeComparator() #print outputs

if __name__ == '__main__':
    main()