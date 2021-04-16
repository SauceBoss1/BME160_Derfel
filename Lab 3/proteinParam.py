#!/usr/bin/env python3
# Name: Derfel Terciano (dtercian)
# Group Members: None
'''
proteinParam.py
This program takes in a string of Amino Acids (also known as a protein) and the program returns the caluclated 
physical-chemical properties of a protein sequence similar to what ProtParam outputs. 

Notes: 
    The pI is returned under the Lehninger method and not with the IPC protein method.
    Therefore, use isoelectric.org to see if the pI returned is correct

classes: ProteinParam (methods will be revealed under ProteinParam docstring)

Assumptions:
    -The input sequence is not guaranteed to be uppercase and might contain unexpected characters.
    -Only count the following (A, C, D, E, F, G, H, I, L, K, M, N, P, Q, R, S, T, V, Y, W) or the lower-case equivalents, and ignore anything else
'''

import numpy
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


#private methods for binary search
#These methods are not meant to be accessed by anyone else therefore, I will be using name mangling to make them private

#more documentation on binary search will be in method __dupBinSearch
    def __leftMost(self,listX,target, precision=2):
        '''
        Takes in a list and searches and returns the first occurance of the duplicate target (in the list) based of the precision

        optional parameter: precision (by default the algorithm withh search up to 2 decimal places)

            -input: a list, target value, [precision (set to 2 by default)]
            -output: index (an integer) of the first occurance fo the duplicate item
        
        Example Input: [-1.00, 0.0005, 0.00007, 5.00]
        Output (of the example): 1
        
        '''
        target = target*(.1**(-1*precision)) #since we're searching with a precision value, we need to change the decimal place of the target value
        leftBound=0 #set beginning index at first index of list
        rightBound=len(listX)-2 #set last index at the last element of the list
        while leftBound < rightBound: 
            mid = (leftBound+rightBound)//2 #in order for binary search to work properly, we need to find the middle element which is found by floor dividing the lest and right bounds by 2
            midVal = int(listX[mid]*(0.1**(-1*precision))) #this makes the value of the middle element the correct precision (by default this will only compare up to 2 decimal places)

            #the following lines checks to see if the value of the element in the middle of the 2 bounds is less than or greater than the targeted value
            if midVal < target:
                leftBound = mid +1 #if the midpoint value is less than the targeted value, move the left bound of the search to the right of the mid point
            else:
                rightBound = mid #if the midpoint value is greater than the targeted value, move the right bound to the position of the mid point
        return leftBound
    
    def __rightMost(self,listX,target,precision=2):
        '''
        Takes in a list and searches and returns the last occurance of the duplicated target (in the list) based of the precision

        optional parameter: precision (by default the algorithm withh search up to 2 decimal places)

            -input: a list, target value, [precision (set to 2 by default)]
            -output: index (an integer) of the lasr occurance fo the duplicate item
        
        Example Input: [-1.00, 0.0005, 0.00007, 5.00]
        Output (of the example): 1
        
        '''
        #for the 3 lines below refer to the first 3 comments in __leftMost for the setup
        target = target*(.1**(-1*precision)) 
        leftBound=0
        rightBound=len(listX)-1
        while leftBound < rightBound:
            mid = (leftBound+rightBound)//2 #in order for binary search to work properly, we need to find the middle element which is found by floor dividing the lest and right bounds by 2
            midVal = int(listX[mid]*(0.1**(-1*precision))) #this makes the value of the middle element the correct precision (by default this will only compare up to 2 decimal places)
            if midVal > target: 
                rightBound = mid #if the value of the midpoint is bigger than the targeted value, then make the right bound the current position of the mid point
            else:
                leftBound = mid + 1 #if the value of the midpoint is smaller than the targeted value, then move the left bound to the right of the midpoint
        return rightBound-1

    def __dupBinSearch(self,listX,target,precision=2):
        '''
        Does a duplicate(in the context of this program it's an approximate) binary search of a given list and a target value.
        The method then returns a list of values that it found were duplicate values based on the precision value (2 decimal places is the default).\n
        Why would we use a binary search? Binary searches are typically faster than a linear search. This creates low run times and high efficency!

            -inputs: a list, a target value,[precison (by default the precision is 2 decimal places)]
            -outputs: a list of item(s) that were found to be duplicates based on the precision

        -Example input: [-1, 0.00004, 0.0003, 10], target = 0,
        -Output of Example: [0.00004, 0.0003] (this is because since our precision is 2 decimal places (by default), 0.00004 and 0.0003 is actually 0.00 in the algorithm)

        Credits: my readings on the binary searches come mainly from Wikipedia (https://en.wikipedia.org/wiki/Binary_search_algorithm)
        '''
        listX=sorted(listX) #in order for binary searches to work, arrays or lists must be sorted from most negative to most positive items
        leftIndex=self.__leftMost(listX,target,precision) #this uses binary search to find the first index of the duplicate element based on precision
        rightIndex=self.__rightMost(listX,target,precision) #this uses binary search to find the last index of the duplicate element based on precision
        return [listX[items] for items in range(leftIndex,rightIndex+1)] #returns a list of what the algorithm found were duplicates based on the the precision value 
    
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
        if self.aaCount() == 0: #if there are no valid AAs in the given string, then just return 0
            return 0

        chargeOfAA = [self._charge_(pH) for pH in numpy.arange(0,14.01,0.01)] #cosntructs a list net charges from pH range 0-14 in increments of 0.01 (there will be 1400 elements in list!!)
        #note: each index in the list is correlated to a pH

        smallestCharge = self.__dupBinSearch(chargeOfAA,0) #returns list of all the charges closest to zero (this is because -0.007 and 0.005 is the same as 0.00 in the algorithm)
        #I could've written the enitre pI method in 4 lines of code using the min() function but as I state before:
        #there are 1400 elements in the list chargeOfAA, its much more efficient if I find the min() of a list of 10 elements than a list of 1400 elements.

        smallestAbsCharge = [abs(netCharge) for netCharge in smallestCharge] #in order to find the smallest charge closest to zero, I have to make all elements in smallestCharge positive/

        if min(smallestAbsCharge) not in chargeOfAA: #sometimes, the charge closest to zero is a negative number and keep in mind I made all values positive in smallestAbsCharge
            return chargeOfAA.index(min(smallestAbsCharge)*-1) * 0.01 #this coverts the smallest charge back to negative number
        return chargeOfAA.index(min(smallestAbsCharge)) * 0.01 #returns the index of the charge closest to zero then I multiply the index by 0.01 in order to get the correct decimal place

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
        myMW =  self.molecularWeight() #gets moleculat weight of protein
        return self.molarExtinction(cystine) / myMW if myMW else 0.0 #uses molar extinction to get mass extinction

    def molecularWeight (self):
        '''Calculates the molecular weight of the protein with H2O taken into account'''
        if self.aaCount() > 0: #if there are no AAs, there should be no calculations taking place
            return self.mwH2O + sum((self.aaComp[aa]*(self.aa2mw[aa]-self.mwH2O)) for aa in self.aaComp) #uses the sum of the number of the AA times the (AA weight - H2O weight) then adds mwH2O at the end
        else:
            return 0

# Please do not modify any of the following.  This will produce a standard output that can be parsed
    
import sys
def main():
    '''Main execution sequence that produces a formatted input and output of the protein and its properties'''
    inString = input('protein sequence?')
    while inString :
        myParamMaker = ProteinParam(inString)
        myAAnumber = myParamMaker.aaCount()
        print ("Number of Amino Acids: {aaNum}".format(aaNum = myAAnumber))
        print ("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
        print ("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
        print ("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
        print ("Theoretical pI: {:.2f}".format(myParamMaker.pI()))
        print ("Amino acid composition:")
        
        if myAAnumber == 0 : myAAnumber = 1  # handles the case where no AA are present 
        
        for aa,n in sorted(myParamMaker.aaComposition().items(), 
                           key= lambda item:item[0]):
            print ("\t{} = {:.2%}".format(aa, n/myAAnumber))
    
        inString = input('protein sequence?')

if __name__ == "__main__":
    main()