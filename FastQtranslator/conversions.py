#!/usr/bin/env python3
# Collaborators: Sagarika Kannoly (skannoly) & Derfel Terciano (dtercian)
'''
Convert PHRED scores to different PHRED formats

Written By: Derfel Terciano & Sagarika Kannoly

Summary: 
    -take in a quality line, sequence, and return the newly converted PHRED quality line from the specified output

Objects:
    -PhredConversions (qScore, seq, clOut = False) -> converts all quality scores to a spcifed PHRED format

Notes:
    -Solexa conversions run in around 3 seconds
    -All other conversions run in 0.3 seconds
'''

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
    '''
    Convert all quality scores to a specifed PHRED format

    Written by: Derfel Terciano & Sagarika Kannoly

    Methods:
        -getNewSeq() -> returns the converted sequence (is converted at all)
        -outputConverter() -> returns the PHRED33 or PHRED64 offset
        -p33toPhred() -> converts PHRED33 to either phred 33 or phred 64
        -p64toPhred() -> converts PHRED64 to either phred 33 or phred 64
        -p64bToPhred() -> converts PHRED64 with B offset to either phred 33 or phred 64
        -pSolToPhred() -> converts PHRED64 with Solexa scoring to either phred 33 or phred 64

    '''
    def __init__(self, qScores, seq, clOut=False):
        '''Initialize object with a required quality line, sequence, and an optional output of either p33[default] or p64 format'''
        self.clOut = clOut #True -> use p64 format // False -> use p33 format
        self.qScore = qScores
        self.newSeq = seq  # this is only needed if we use B offset

    def getNewSeq(self):
        '''return the sequence'''
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
        for char in self.qScore: #iterate through each quality character
            rawQscore = ord(char) - 33 #find its quality score in the (0:40) range
            phredQscore += phredQscore.join(chr(rawQscore + self.outputConverter(self.clOut))) #covert each rawQscore to either p33 or p64 and then find its ascii char
        return phredQscore

    def p64toPhred(self):
        '''convert p64 to qScores that range from (0:40)'''
        phredQscore = ''
        for char in self.qScore: #iterate through each quality character
            rawQscore = ord(char) - 64 #find its quality score in the (0:40) range
            phredQscore += phredQscore.join(chr(rawQscore + self.outputConverter(self.clOut))) #covert each rawQscore to either p33 or p64 and then find its ascii char
        return phredQscore

    def p64bToPhred(self):
        '''convert p64b offset to qScores that range from (0:40)'''
        newSeq = self.newSeq
        phredQscore = ''
        charPos = 0
        for char in self.qScore: #iterate through each quality character
            if ord(char) == 66: #if we see 'B' then it's considered the lowest score so assign it to 0
                rawQscore = 0
            else:
                rawQscore = ord(char) - 64 #else find its qScore in the (0:40) range

            if rawQscore == 0: #when we find 'B' or the lowest qScore, change it's corresponding base position in the seq to N
                newSeq = newSeq[:charPos] + 'N' + newSeq[charPos + 1:]
            charPos += 1 #increase the character position
            phredQscore += phredQscore.join(chr(rawQscore + self.outputConverter(self.clOut))) #covert each rawQscore to either p33 or p64 and then find its ascii char
        return phredQscore, newSeq

    def pSolToPhred(self):
        '''convert p64Solexa offset to qScores that range from (0:40)'''
        phredQscore = ''
        for char in self.qScore:
            solQ = ord(char) #find the ascii integer of the current character
            rawQ = solexaToPhred[asciiToSolexa[solQ]] #using dictionary mappings, find the qScore in the (0:40) range
            phredQscore += phredQscore.join(chr(rawQ + self.outputConverter(self.clOut))) #covert each rawQscore to either p33 or p64 and then find its ascii char
        return phredQscore