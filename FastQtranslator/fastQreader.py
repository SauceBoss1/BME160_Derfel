#!/usr/bin/env python3
# Collaborators: Sagarika Kannoly (skannoly) & Derfel Terciano (dtercian)


import sys
class FastQReader:
    def __init__(self, fname = ''):
        self.fname = fname
    
    def doOpen(self):
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFastQ(self):

        header = ''
        sequence = ''
        qDesc = ''
        qScore = ''

        with self.doOpen() as f:
            header = ''
            sequence = ''
            qDesc = ''
            qScore = ''

            lineNum = 0
            for line in f:
                if line.rstrip():
                    lineNum += 1
                    lineNumMod = lineNum % 4 
                    if lineNumMod == 1:
                        header = line[1:].rstrip()
                    elif lineNumMod == 2:
                        sequence = line.rstrip().replace('.','N').replace('*', 'N').upper()
                    elif lineNumMod == 3:
                        qDesc = line[1:].rstrip()
                    elif lineNumMod ==0:
                        qScore = line.rstrip()
                        yield header, sequence, qDesc, qScore
                        header = ''; sequence = ''; qDesc = ''; qScore = ''
                    