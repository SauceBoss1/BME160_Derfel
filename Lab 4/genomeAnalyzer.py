#!/usr/bin/env python3
# Name: Derfel Terciano (dtercian)
# Group Members: None
'''
genomeAnalyzer.py

Written by: Derfel Terciano

This program prints the sequence length (in Mb), GC content (in %), and the relative codon usage from a given FastA file's genomic info.

Class(es)/object(s):
    -GenomeAnalyzer: => handles the output formatting of a given genome (from the FastA file)

Assumptions:
    - FastA files are the only files that is to be used
'''
import sequenceAnalysis #type: ignore
import sys
#from importlib import reload
#reload(sequenceAnalysis)

class GenomeAnalyzer:
    '''
    Formats and prints relvant genomic info about the given FastA file.

    Written by: Derfel Terciano

    Relevant info includes:
        -Sequence length (in Mb)
        -GC content (in %)
        -codon usage
    
    Notes:
        -Files can be specified directly when initializing the class or through the command line's STDIN

    public methods:
        -analyzeGenome(): => prints a formatted output of the relevant FastA info
    '''
    def __init__(self,fileName=''):
        '''
        Initalizes the object by reading FastA file and updating internal info about the data in the file

        Parameter(s) (optional): [fileName] (files can be specified through initialization or by STDIN)
        '''
        self.myReader=sequenceAnalysis.FastAreader(fileName) #creates a FastA object in order to read the inputted file
        self.myNuc=sequenceAnalysis.NucParams() #creates a NucParams object in order to store genomic data from the FastA file
        for head,seq in self.myReader.readFasta(): #after reading the FastA file, update the NucParams object
            self.myNuc.addSequence(seq)
        

    def analyzeGenome(self):
        '''
        Prints a formatted output of the relevant info from the genomic content of the FastA file
        '''
        if (self.myNuc.nucCount() == 0) or (sum(self.myNuc.codonComposition().values()) == 0): #if there are no nucleotides or valid codons then end program
            print('WARNING NO VALID CODON OR NUCLEOTIDES FOUND. TERMINATING PROGRAM')
            exit() #from sys

        numBases = self.myNuc.nucCount()/1000000 #Mega Bases is the nucleotide length divded by 1 million
        contentGC = ((self.myNuc.nucComposition()['C'] + self.myNuc.nucComposition()['G']) / self.myNuc.nucCount()) *100 #GC content is the sum of Gs and Cs in nucleotide div. by nucCount

        print(f'sequence length = {numBases:0.2f} Mb\n') #prints seq length to 2 decimal places
        print(f'GC content = {contentGC:0.1f}%\n') #returns the percentage to 1 decimal place

        #the loop below prints the codon usage of each codon
        for codon, aa in sorted(self.myNuc.rnaCodonTable.items(), key=lambda a:(a[1],a[0])): #ordered by AA and then by the codon itself

            #the following if/else statment prevents a divide by zero error if there was no codon count for the current codon
            if self.myNuc.aaComposition()[aa] == 0: #if there are no counts of the current aa, set to 1 to preven divide by zero error
                currentAAcomp = 1
            else:
                currentAAcomp = self.myNuc.aaComposition()[aa]

            
            codonFreq = (self.myNuc.codonComposition()[codon] / currentAAcomp) * 100 #finds codon freq from the # of codons divided by # of amino acid of that codon
            print (f'{codon} : {aa} {codonFreq:5.1f} ({self.myNuc.codonComposition()[codon]:6d})')



def main ():
    '''Executes the main sequence of the program'''
    genome = GenomeAnalyzer()
    genome.analyzeGenome() #calls the formatted input

if __name__ == '__main__':
    main()
