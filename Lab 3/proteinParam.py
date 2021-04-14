#!/usr/bin/env python3
# Name: Derfel Terciano (dtercian)
# Group Members: None


import numpy
class ProteinParam :
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
        self.protein = protein.upper()
        self.aaComp = {key:0 for key in self.aa2mw} #creates a dictionary of all AA keys but set to 0
        for char in self.protein:
            if char in self.aaComp:
                self.aaComp[char] += 1

#private helper methods 
#These methods are not meant to be accessed by anyone else therefore, I will be using name mangling to make them private
    def __sumOfPosCharges(self,pH):
        nTerm = (10**self.aaNterm)/(10**self.aaNterm + 10**pH)
        sumOfAA = sum((self.aaComp[aa])*(10**self.aa2chargePos[aa]/(10**self.aa2chargePos[aa]+10**pH)) for aa in self.aa2chargePos)
        return nTerm + sumOfAA
    def __sumOfNegCharges(self,pH):
        cTerm = (10**pH)/(10**self.aaCterm + 10**pH)
        sumOfAA = sum((self.aaComp[aa])*(10**pH/(10**self.aa2chargeNeg[aa]+10**pH)) for aa in self.aa2chargeNeg)
        return cTerm + sumOfAA


#private methods for binary search
#These methods are not meant to be accessed by anyone else therefore, I will be using name mangling to make them private
    def __leftMost(self,listX,target, precision=2):
        target = target*(.1**(-1*precision))
        leftBound=0
        rightBound=len(listX)-2
        while leftBound < rightBound:
            mid = (leftBound+rightBound)//2
            midVal = int(listX[mid]*(0.1**(-1*precision)))
            if midVal < target:
                leftBound = mid +1
            else:
                rightBound = mid
        return leftBound
    
    def __rightMost(self,listX,target,precision=2):
        target = target*(.1**(-1*precision))
        leftBound=0
        rightBound=len(listX)-1
        while leftBound < rightBound:
            mid = (leftBound+rightBound)//2
            midVal = int(listX[mid]*(0.1**(-1*precision)))
            if midVal > target:
                rightBound = mid
            else:
                leftBound = mid + 1
        return rightBound-1

    def __dupBinSearch(self,listX,target,precision=2):
        listX=sorted(listX)
        left=self.__leftMost(listX,target,precision)
        right=self.__rightMost(listX,target,precision)
        duplicates = [abs(listX[items]) for items in range(left,right+1)]
        smallestValue=2**256
        for item in duplicates:
            if item < smallestValue:
                smallestValue = item
        return smallestValue
    
#public methods
    def aaCount (self):
        return sum(self.aaComp[aa] for aa in self.aaComp)

    def pI (self,useBinary=False):
        if self.aaCount() == 0:
            return 0
        if useBinary:
            chargeOfAA = [self._charge_(pH) for pH in numpy.arange(0,14.01,0.01)]
            smallestCharge = self.__dupBinSearch(chargeOfAA,0)
            return chargeOfAA.index(smallestCharge) * 0.01
        else:
            chargeAA = [abs(self._charge_(pH)) for pH in numpy.arange(0,14.01,0.01)]
            return chargeAA.index(min(chargeAA))*.01

    def aaComposition (self) :
        return self.aaComp

    def _charge_ (self,pH):
        return self.__sumOfPosCharges(pH) - self.__sumOfNegCharges(pH)

    def molarExtinction (self):
        return (self.aaComp['Y']*self.aa2abs280['Y'])+(self.aaComp['W']*self.aa2abs280['W'])+(self.aaComp['C']*self.aa2abs280['C'])

    def massExtinction (self):
        myMW =  self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0

    def molecularWeight (self):
        if self.aaCount() > 0:
            return self.mwH2O + sum((self.aaComp[aa]*(self.aa2mw[aa]-self.mwH2O)) for aa in self.aaComp)
        else:
            return 0

# Please do not modify any of the following.  This will produce a standard output that can be parsed
    
import sys
def main():
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