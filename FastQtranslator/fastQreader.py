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

        header = ''
        sequence = ''
        qDesc = ''
        qScore = ''

        with self.doOpen() as f:
            #reset headings
            header = ''
            sequence = ''
            qDesc = ''
            qScore = ''
            
            line = f.readline()
            while not line.startswith('@'):
                line = f.readline()
            header = line[1:].rstrip()

            for line in f:
                if line.startswith('@') and line[1:len(line)].isupper() and qScore != '':
                    yield header, sequence, qDesc, qScore
                    header = line[1:].rstrip()
                    sequence =''; qDesc = ''; qScore = ''
                elif line.startswith('+'):
                    qDesc = line[1:].rstrip()
                elif (qDesc == '') and (qScore == ''):
                    sequence = line.rstrip()
                elif (sequence != '') and (qDesc != ''):
                    qScore = line.rstrip()
        yield header, sequence, qDesc, qScore