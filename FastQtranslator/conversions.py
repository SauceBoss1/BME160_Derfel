#!/usr/bin/env python3
# Collaborators: Sagarika Kannoly (skannoly) & Derfel Terciano (dtercian)
'''Functions for all phred conversions'''

# ########################################
# HOW TO APPROACH THE PROJECT

# 1) Make functions that corresponds to each CL input and extract it's qScores in the (0:40) range
# 2)Take those (0:40) scores and convert it to either p33 or p64 (using seperate functions)

# EXCEPTIONS:
#    B offset requires change of actual bases
#    1) if a base has a score of B, then change the base to N

# Edits we may need to make:
#     1) everything (except for Solexa) can be done mathmatically but we may need to make dictionaries if program becomes too slow
# ########################################

from solexaDictionary import solexaToPhred #remember to have this file in the same folder

class PhredConversions:
    def __init__(self, qScores, seq, inCl):
        self.inOptions = ['p33','p64','p64b','pSol'] #list allowed options
        self.orderedQscores = []
        self.newSeq = seq #this is only needed if we use B offset
    
    def p33toPhred (self, qScores):
        '''convert p33 to qScores that range from (0:40)'''
        pass

    def p64toPhred (self, qScores):
        '''convert p64 to qScores that range from (0:40)'''
        pass
    
    def p64bToPhred (self, qScore, seq):
        '''convert p64b offset to qScores that range from (0:40)'''
        newSeq = seq
        phredQscore = []
        charPos = 0
        for char in qScore:
            rawQscore = ord(char)-66
            phredQscore.append(rawQscore)
            if rawQscore == 0:
                newSeq = newSeq[:charPos] + 'N' + newSeq[charPos+1:]
            charPos += 1
        self.newSeq = newSeq
        self.orderedQscores = phredQscore

    def pSolToPhred (self, qScore):
        '''convert p64Solexa offset to qScores that range from (0:40)'''
        phredQscore = []
        for char in qScore:
            phredQscore.append(solexaToPhred[ord(char)])
        self.orderedQscores

    def toP33 (self):
        '''convert raw qScores that range from (0:40) to p33'''
        newQScore = ''
        for qScore in self.orderedQscores:
            newQScore += newQScore.join(chr(qScore+33))
        return newQScore

    def toP64 (self):
        '''convert raw qScores that range from (0:40) to p64'''
        newQScore = ''
        for qScore in self.orderedQscores:
            newQScore += newQScore.join(chr(qScore+64))
        return newQScore