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

from solexaDictionary import solexaToPhred, asciiToSolexa  # remember to have this file in the same folder


class PhredConversions:
    def __init__(self, qScores, seq, clOut=False):
        self.clOut = clOut
        self.qScore = qScores
        self.newSeq = seq  # this is only needed if we use B offset

    def getNewSeq(self):
        return self.newSeq

    def outputConverter(self, clOut):
        '''Return the p33 or p64 ascii offset'''
        if clOut:
            return 64
        else:
            return 33

    def p33toPhred(self):
        '''convert p33 to qScores that range from (0:40)'''
        phredQscore = ''
        for char in self.qScore:
            rawQscore = ord(char) - 33
            phredQscore += phredQscore.join(chr(rawQscore + self.outputConverter(self.clOut)))
        return phredQscore

    def p64toPhred(self):
        '''convert p64 to qScores that range from (0:40)'''
        phredQscore = ''
        for char in self.qScore:
            rawQscore = ord(char) - 64
            phredQscore += phredQscore.join(chr(rawQscore + self.outputConverter(self.clOut)))
        return phredQscore

    def p64bToPhred(self):
        '''convert p64b offset to qScores that range from (0:40)'''
        newSeq = self.newSeq
        phredQscore = ''
        charPos = 0
        for char in self.qScore:
            if ord(char) == 66:
                rawQscore = 0
            else:
                rawQscore = ord(char) - 64

            if rawQscore == 0:
                newSeq = newSeq[:charPos] + 'N' + newSeq[charPos + 1:]
            charPos += 1
            phredQscore += phredQscore.join(chr(rawQscore + self.outputConverter(self.clOut)))
        return phredQscore, newSeq

    def pSolToPhred(self):
        '''convert p64Solexa offset to qScores that range from (0:40)'''
        phredQscore = ''
        for char in self.qScore:
            solQ = ord(char)
            rawQ = solexaToPhred[asciiToSolexa[solQ]]
            phredQscore += phredQscore.join(chr(rawQ + self.outputConverter(self.clOut)))
        return phredQscore