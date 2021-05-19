#!/usr/bin/env python3
# Name: Derfel Terciano (dtercian)
# Group Members: List full names (CATS usernames) or “None”
'''
findUnique.py

Find all unique and essential subsequence of each tRNA from a population from other tRNA sequnces

written by Derfel Terciano

Summary:
    -Stream in a FastA file of all tRNA sequnces and return a formatted output of all the unique and essential subsequence of each sequence

NOTES:  
    -code runs in 0.2 seconds
    -Unicode character are processed differently on different machines so unicode characters may look different
    -Without comments, the code is around ~70 lines
'''
# Design for findUnique

# 1) input all sequences into class
# 2) find power set for all of the tRNA sequences and store it into the class
# 3) have a function that takes in a header, and do all the set comparison based on the header and sequence
# 4) have the final output be created based on calling the header of each sequence           
import sys
import codecs
import time
start_time = time.time() #initialize a time object (for determing run-time)
from sequenceAnalysis import FastAreader #type: ignore

class FindUnique:
    '''
    Find all unique and essential subseqences of all tRNAs
    
    written by: Derfel Terciano

    methods:
        -addSeq(seq, header) => add a tRNA sequence to the tRNA dictionary and list for the class
        -powerSet(inSeq) => return a powersey
        -uniqueFinder(curHeader) => Based on the specified header, return the unique subsequences of that tRNA
        -essential(curHeader) => Based on the specified header, return the essential subsequences of that tRNA
        -findAllOccurances(subSeq, mainSeq) => find all instances of a specified subsequence in a sequence
        -findUnique() => find the unique and essential subsequences of ALL tRNAs
    '''
    def __init__(self, header='', seq=''):
        '''Instantiate the object with an empty list of headers, tRNA pSets, and an empty dictionary of tRNA information.'''
        self.tRnaInfo ={} #the form of the dictionary is in the following => header: [seq, powerset(seq)] (this is good for comparing current sets)
        self.headers = [] #store headers
        self.allTrnaPset = [] #store all pSets from all 22 sequences (useful for unioning without using for loops)
        self.addSeq(seq, header) #this was an idea that was inspired from lab 4

    def addSeq(self, seq, header):
        '''add a header and a sequence to their respective data structures in the object'''
        header = header.replace(' ','') #remove all spaces from header
        self.seq= seq.upper().replace('-','').replace('_','').replace('.','') #remove all alignment characters
        if header != '': #since the default instance is an empty header and sequence, we must ignore that
            self.headers.append(header) #store header
            seqPset = self.powerset(self.seq) #find power set of seq
            self.tRnaInfo[header] = [self.seq, seqPset] #add the header, and all its info into the dictionary
            self.allTrnaPset.append(seqPset) #add the pSet of this header to the allTrnaPset list

    def powerset(self,inSeq):
        '''find the power set of a given sequence'''
        pSet= set()
        for i in range(len(inSeq)): #iterate through each character in the sequence
            for j in range(i+1,len(inSeq)+1): #find all logical subsequences for each character in the sequence
                pSet.add(inSeq[i:j]) 
        return pSet #remeber that a powerset is a set

    def uniqueFinder(self, curHeader):
        '''Find all unique subsequences of a given header'''
        allSets = self.allTrnaPset.copy() #retrieve a copy of all the tRNApSets we have found (we do this so we don't accidentally tamper with the list)
        allSets.remove(self.tRnaInfo[curHeader][1]) #remove the pSet of the given header (this will avoid unioning all 22 pSets since we only want to union 21 pSets)
        return set(self.tRnaInfo[curHeader][1]) - set().union(*allSets) #(current pSet) - (union of all other pSets)

    def essential(self, curHeader):
        '''Find all essential subsequences of each tRNA sequence'''
        nonEssentials = set() #we need to make a set of all non-essential and unique subsequences
        uniqueSet = self.uniqueFinder(curHeader) 
        for rawUnique in uniqueSet:
            if (rawUnique[:-1] in uniqueSet) or (rawUnique[1:] in uniqueSet):
                nonEssentials.add(rawUnique)
        return set(uniqueSet) - nonEssentials #remove all nonessentials from the unique set

    def findAllOccurances(self, subSeq, mainSeq):
        '''Find all occurances of a given subsequence in a sequence'''
        posFound = []
        index = mainSeq.find(subSeq) #finds the first instance of the occuring subsequence
        while index != -1: #.find() returns negative 1 if nothing can be found
            posFound.append(index) #append the current instance of the subsequence
            index = mainSeq.find(subSeq, index+1) # find another instance of the subsequence by move the start position to the right by 1 posisiton
        return posFound
    
    def findUnique(self): 
        sys.stdout.reconfigure(encoding='utf-8') #we need to encode for utf8 since we will get an error if we don't encode it
        with sys.stdout as outFile:
            for currentHeader in sorted(self.headers): #adjust order of sequences here (iterate through the list of headers)
                print(f'{currentHeader}\n{self.tRnaInfo[currentHeader][0]}') #print general tRNA info
                essentials = self.essential(currentHeader) #find the essentials of the current header
                essentialsIndex = {} #in the form of index:essentialSubSequence (it's more convient to do it this way so it's easier to print the dictionary)
                for essentialSubs in essentials: #iterate through all essential subsequences
                    for indexes in self.findAllOccurances(essentialSubs, self.tRnaInfo[currentHeader][0]): #find all occurances of the essential subsequence
                        essentialsIndex[indexes] = essentialSubs #add the occurance and the essential sub to the dinctionary
                for index, seq in sorted(essentialsIndex.items(), key = lambda a:a[0]): #print out the essentialsIndex dictionary
                    print(f'{index * "."}{seq}')
            #print('run time (in seconds) : ', time.time()-start_time) #uncomment this line if you want to see the printed run time of the program at the end of output file

def main(inCL=None):
    '''Execute main program, stream in file, stream out an output file'''
    inFile = FastAreader()
    genes = FindUnique()
    for header,seq in inFile.readFasta(): #concept of using a for loop to add sequences into the object taken from Lab 4
        genes.addSeq(seq, header)
    genes.findUnique()
if __name__ == "__main__":
    main()  