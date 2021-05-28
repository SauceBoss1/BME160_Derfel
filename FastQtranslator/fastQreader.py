#!/usr/bin/env python3
# Name: Derfel Terciano (dtercian)
# Group Members: Sagarika Kannoly (skannoly)


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

        # header = ''
        # sequence = ''
        # qDesc = ''
        # qScore = ''

        with self.doOpen() as f:
            #reset headings
            header = ''
            sequence = ''
            qDesc = ''
            qScore = ''

            #skip to first fastQ header
            line = f.readline()
            while not line.startswith('@'):
                line = f.readline()
            header = line[1:].rstrip()

            for line in f:
                if line.rstrip():
                    if line.startswith('@'):
                        yield header, sequence, qDesc, qScore
                        header = line[1:].rstrip()
                        sequence = '' ; qDesc = '' ; qScore = ''
                    elif line.startswith('+'):
                        qDesc = line[1:].rstrip()
                    elif sequence != '':
                        qScore = line.rstrip()
                    else:
                        sequence = line.rstrip()
        yield header, sequence, qDesc, qScore