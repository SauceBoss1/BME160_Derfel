#!/usr/bin/env python3
# Name: Derfel Terciano (dtercian)
# Group Members: None
'''
sequenceAnalysis.py v1

written by: Derfel Terciano (with some source code given by David Bernik)

This module provides tools for biological data. These tools are mainly helpful for the bioinformatics field but can be used by anyone else. \n
So far, the main functions of the module include: Nucleotide parameters, Protein parameters, and a FastA reader. This module will become heavily developed further down the line.
Each of the objects in this module will also have its own design specs, overview, and assumptions in their respective docstrings.

Objects in this class:
    - NucParams([inString]) -> a collection of methods that can determine amino acid composition, nucleotide composition, codon composition, and nucletide count.
        -Initialized with a string of nucleotides
    - ProteinParam(protein(string type)) -> a collection of methods that can determine aacount, pI, aa composition, molar & mass extinc. coeff., and molec. weight.
        -Initialized with a string of amino acids (also known as a protein)
    - FastAreader([fname (a file name or path of type string)]) -> a collection of methods that can read and return FastA file info through either a file name or STDIN
        -Initialized with a filename, path (both of these are string types) or through STDIN from the console.
    -More features coming soon?
'''
import numpy
import sys


#TODO: DOCUMENTATION AND EXCEPTIONAL CASES (CONSIDER ONLY Ns AS AN INPUT)
class NucParams:
    '''
    Calculates compositional properties of a given string of nucleotides.

    Written by Derfel Terciano

    With a given string of nucleotides, the public methods will do the following:
        -addSequence(sequence): updates all relevant tables and data structures within the object with the inputted string of nucs
                -returns nothing
        -aaComposition(): returns a dictionary of the composition of amino acids that are translated in the string of nucs
        -nucComposition(): returns a dictionary of the composition of nucleotides with the string. This is based off the {ATCGUN} alphabet
        -codonComposition(): returns the composition of valid codons in a given string of nucleotides. 
            -(This means bases with N will be ignored and not be counted)
        -nucCount(): returns the count (sum) of all nucleotides that are in the {ACTGUN} alphabet

    Assumptions:
        - {ACTGN} is the only valid alphabet we will be using
        - Codons (nucleotides of 3) will be ignored if they contain at least one 'N'
    '''
    rnaCodonTable = {
    # RNA codon table
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}

    def __init__ (self, inString=''):
        '''
        Initializes the NucParams object
            parameter (optional): inString (string type)
        '''
        self.nucComp = { 'A':0, 'T':0, 'C':0, 'G':0, 'U':0, 'N':0 } # nucleotide composition dictionary
        self.aaComp = {aa:0 for aa in self.rnaCodonTable.values()} # amino acid composition dictionary
        self.codonComp = {aa:0 for aa in self.rnaCodonTable} # codon composition dictionary
        self.addSequence(inString)#updates the dictionaries if the object is initialized with a string (which is optional)


    def addSequence (self, inSeq):
        '''
        updates the object's corresponding dictionaries with relevant info from an entered string
            -parameters: inSeq (a string type)
        '''
        inSeq = inSeq.replace(' ','').upper() #this makes sure that theres no spaces in sequence

        for codonIndex in range(0,len(inSeq),3): #deals with the codon comp and increments in multiples of 3
            codonString = inSeq[codonIndex:codonIndex+3].replace('T','U') #extracts the string from the codonIndex and converts DNA to RNA
            if codonString in self.codonComp: #if the codon is valid that means there's a matching AA as well
                self.codonComp[codonString] += 1 #adds codonString to the valid key in codonComp
                self.aaComp[self.rnaCodonTable[codonString]] += 1 #adds the valid AA to the valid key in aaComp
            

        for character in inSeq: #deals with the nucleotide composition
            character = character #makes the current char uppercase
            if character in self.nucComp: #checks if the character in inSeq is a valid nucleotide in the {ATCGUN} alphabet
                self.nucComp[character]  += 1



    def aaComposition(self):
        '''returns a dictionary of the composition of amino acids that are translated in the string of nucs'''
        return self.aaComp

    def nucComposition(self):
        '''returns a dictionary of the composition of nucleotides with the string. This is based off the {ATCGUN} alphabet'''
        return self.nucComp

    def codonComposition(self):
        '''returns the composition of valid codons in a given string of nucleotides.'''
        return self.codonComp

    def nucCount(self):
        '''returns the count (sum) of all nucleotides that are in the {ACTGUN} alphabet'''
        return sum(self.nucComp.values())

class ProteinParam :
    '''
    Calculate the physical-chemical properties of a given protein string (a string of amino acids)

    Written by: Derfel Terciano\n
    
    With a given string from the cosntruction of the class, the following public methods will return the following:
            -aaCount(): the number (integer or float) of VALID amino acids in the given string
            -pI(): the theortical isoelectric point (a float) of the given protein
            -aaComposition(): a dictionary of all of the number of amino acids present (0 means that the specific amino acid is not present)
            -molarExtinction(,[cystine]): a float of the molar extinction coefficent of the protein (the default input is true; if False, Cystine calculations will be ignored)
            -massExtinction(,[cystine]): a float of the the mass extinction coefficient of the protein (the default input is true; if False, Cystine calculations will be ignored)
            -molecularWeight(): the TOTAL molecular weight of the given protein (a float) 
    Notes: 
    The pI is returned under the Lehninger method and not with the IPC protein method.
    Therefore, use isoelectric.org to see if the pI returned is correct

    Input: String\n
    Output of all methods: Floats

    ProteinParam also uses certain private methods to complete operations in the public methods:
            -_charge_(pH): calculates the net charge on the protein at a specific pH (a float)
            -__sumOfPosCharges(pH): sums up the all the positive charges in the given protein with the number of K, R, and, H present with N-Terminus
            -__sumOfNegCharges(pH): sums up all the negative charges given in the protein with number of D, E, C, and Y present with the C-Terminus
            -__leftMost(listX,target, precision): Used in duplicate binary searches and finds the first occurance of a duplicate value (returns an index)
            -__rightMost(listX,target,precision): Used in duplicate binary searches and finds the last occurance of a duplicate value (returns an index)
            -__dupBinSearch(listX,target,precision): returns a list of duplicate items that was found in a given list based on its precision value (default precision search is 2 decimal places) 


    '''
        # These tables are for calculating:
    #     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
    #     absorbance at 280 nm (aa2abs280)
    #     pKa of positively charged Amino Acids (aa2chargePos)
    #     pKa of negatively charged Amino acids (aa2chargeNeg)
    #     and the constants aaNterm and aaCterm for pKa of the respective termini
    #  Feel free to move these to appropriate methods as you like

    # As written, these are accessed as class attributes, for example:
    # ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }
    
    
    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__ (self, protein):
        '''
        Takes in a string for the class's methods to use and stores the valid amino acids in a dictionary
        -Example input: 'VLSPADKTNVKAAW'
        '''
        self.protein = protein.upper() #makes sure that the passed in string is all capitals (makes identifying valid amino acids easier)
        self.aaComp = {key:0 for key in self.aa2mw} #creates a dictionary of all AA keys but set to 0 (uses aa2mw as a reference of keys)
        for char in self.protein: #iterates through each character in given input
            if char in self.aaComp: #checks to see if the character from the given input is a valid
                self.aaComp[char] += 1
        # if self.aaCount() == 0: #If nothing valid is entered, then warn the user that their input is invalid
        #     print('WARNING: No valid Amino Acid(s) were entered. All methods in the class will now return zero (0) :)')
    #private helper methods 
    #These methods are not meant to be accessed by anyone else therefore, I will be using name mangling to make them private.
    def __sumOfPosCharges(self,pH):
        '''
        Takes in a pH can calculates the sum of positive charges in the given protein sequence
        -Note: K,R,H, and the N-Terminus are the only Amino Acids (AAs) with a positive charge
        '''
        nTerm = (10**self.aaNterm)/(10**self.aaNterm + 10**pH) #since the N-terminus is not in the aaComp dictionary, it's easier to calculate this independently
        sumOfAA = sum((self.aaComp[aa])*(10**self.aa2chargePos[aa]/(10**self.aa2chargePos[aa]+10**pH)) for aa in self.aa2chargePos) #sums the charged Amino Acids (with a specific formula) listed in the docstring with the given pH 
        return nTerm + sumOfAA

    def __sumOfNegCharges(self,pH):
        '''
        Takes in a pH can calculates the sum of negative charges in the given protein sequence
        -Note: D, E, C, Y and the C-Terminus are the only Amino Acids (AAs) with a negative charge
        '''
        cTerm = (10**pH)/(10**self.aaCterm + 10**pH) #since the C-terminus is not in the aaComp dictionary, it's easier to calculate this independently
        sumOfAA = sum((self.aaComp[aa])*(10**pH/(10**self.aa2chargeNeg[aa]+10**pH)) for aa in self.aa2chargeNeg) #sums the charged Amino Acids (with a specific formula) listed in the docstring with the given pH 
        return cTerm + sumOfAA


    #private method(s) for binary search
    #These methods are not meant to be accessed by anyone else therefore, I will be using name mangling to make them private

    #binary search
    def __binarySearchOfCharges(self,leftBound, rightBound ,targetValue ,precision=2):
        '''
        Searches for the charge that is closes to thee target value then it returns the pH closest to the target value

                -inputs: a leftbound(float), a rightbound(float), a target value(float), and an optional precision parameter (2 by default)
                -Output: a float of the pH value

            -Example inputs: 0,14,0
            -Output of example: some pH where the net charge is closest to the target

        '''
        left=leftBound
        right=rightBound
        while left < right: #this makes the loop run until both the right and left bound are equal to each other, any more and the loop will break
            midPoint = (left+right)/2 #take the 2 bounds and find the midpoint of them
            chargeATpH = self._charge_(midPoint) # identifies the charge at the pH
            if chargeATpH > targetValue: 
                left = midPoint + (.1**precision) #based on a standard pH and charge graph, if the charge at the midpoint is bigger than target, we need to move the leftbound to the right
            else:
                right = midPoint - (.1**precision) #if the charge at the midpoint is smaller, move the rightbound to the left
        return right #returns the final result (returning either the left or right bound should be ok)


    #public methods
    def aaCount (self):
        '''returns the sum of valid AAs in the aaComposition() dictionary'''
        return sum(self.aaComp[aa] for aa in self.aaComp) #since invalid characters have already been accounted for, it is ok to just sum the value of each key in aaComposition()

    def pI (self):
        '''
        returns the theoretical isoelectric point (Lehninger method) of the given protein string \n
        Uses binary search in order to find the smallest net charges of the protein across the pH range
            -Why binary search? In short its faster and more efficient!
        '''
        if self.aaCount() == 0: #if there are no valid AAs, then return basically nothing
            return 0
        return self.__binarySearchOfCharges(0,14,0) #Searches for the charge that is closes to thee target value then it returns the pH closest to the target value

    def aaComposition (self) :
        '''returns a dictionary with Amino Acids as keys and the number of their respective amino acids found'''
        return self.aaComp #this was done in __init__

    def _charge_ (self,pH):
        '''
        Calculates the net charge of all the AAs with a given pH
        -Note: Only K, R, H, D, E, C, Y have associated charges
        '''
        return self.__sumOfPosCharges(pH) - self.__sumOfNegCharges(pH) #the net charge is found by the sum of all the positive charges - the sum of all negative charges 

    def molarExtinction (self, cystine = True):
        '''
        Calculates the molar extinction coefficent of the protein \n
        Optional parameter: cystine (True by default) If cystine is False, then Cystine and Cysteine is ignored.
        '''
        cysResidue = 0
        if cystine: #if cystine is False, then the cystine part of the calculation will be ignored
            cysResidue = (self.aaComp['C']*self.aa2abs280['C'])
        return (self.aaComp['Y']*self.aa2abs280['Y'])+(self.aaComp['W']*self.aa2abs280['W'])+cysResidue #calculates the molar extinction coeff. of the protein with Y,W,and,C


    def massExtinction (self, cystine = True):
        '''Calculates the mass extinction coefficient of the protein'''
        myMW =  self.molecularWeight() #gets molecular weight of protein
        return self.molarExtinction(cystine) / myMW if myMW else 0.0 #uses molar extinction to get mass extinction

    def molecularWeight (self):
        '''Calculates the molecular weight of the protein with H2O taken into account'''
        if self.aaCount() > 0: #if there are no AAs, there should be no calculations taking place
            return self.mwH2O + sum((self.aaComp[aa]*(self.aa2mw[aa]-self.mwH2O)) for aa in self.aaComp) #uses the sum of the number of the AA times the (AA weight - H2O weight) then adds mwH2O at the end
        else:
            return 0


class FastAreader :
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence
