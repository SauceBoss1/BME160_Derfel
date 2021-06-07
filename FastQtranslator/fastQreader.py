#!/usr/bin/env python3
# Collaborators: Sagarika Kannoly (skannoly) & Derfel Terciano (dtercian)

'''
Read a FastQ file and yield all its fields

Written by: Derfel Terciano & Sagarika Kannoly

Summary: 
    -Stream in (or specify) a fastQ file and yield each read
    -Each read contains a header, sequence, quality description, and the quality scores
'''

import sys
class FastQReader:
    '''
    Read a FastQ file and yield each read

    Written by: Derfel Terciano & Sagarika Kannoly

    instantiation:
        - fastQfile = FastQReader('someFile.fastq')

    usage:
        for header, seq, qDesc, qScore for readFastQ():
            do something
    '''
    def __init__(self, fname = ''):
        '''Save the file name for doOpen()'''
        self.fname = fname
    
    def doOpen(self):
        '''Determine if a file is specified or input is from stdin'''
        if self.fname == '': #if nothing has been specified, then assume input is from sys.stdin
            return sys.stdin
        else:
            return open(self.fname)

    def readFastQ(self):
        '''Read a fastQ file and yield the header, sequence, quality description, and the quality score'''
        header = ''
        sequence = ''
        qDesc = ''
        qScore = ''

        with self.doOpen() as f: #with either a specified file or from stdin
            header = ''
            sequence = ''
            qDesc = ''
            qScore = ''

            lineNum = 0 
            for line in f:
                if line.rstrip(): #skip all empty lines
                    lineNum += 1
                    lineNumMod = lineNum % 4 #since each read is in groups of 4, you can just calculate what each line corresponds to by using modulo 4
                    
                    #check if the line is either the header (line 1), seq (line 2), qDesc (line 3), qScore (line 4)
                    if lineNumMod == 1 and line.startswith('@'): #headers are line 1 and start with an '@'
                        header = line[1:].rstrip() #get rid of the '@' and eny white space
                    elif lineNumMod == 2: #the second line is always the sequence
                        sequence = line.rstrip().replace('.','N').replace('*', 'N').upper() #get rid of any invalid characters
                    elif lineNumMod == 3 and line.startswith('+'): #line 3 is a quality description if it starts with a '+'
                        qDesc = line[1:].rstrip() #get rid of any empty whitespace
                    elif lineNumMod == 0: #line 4 (or 0 since 4%4 = 0) is always the quality scores
                        qScore = line.rstrip()

                        #below is the error handling
                        if (u'\x00' in qScore) or (u'\x00' in sequence): #if we see null characters, then this is a bad read
                            print ('warning bad record', '@'+header, file=sys.stderr)
                        elif header == '': #if we don't see a header since the numbering may have become out of sync, then tell the user its a bad read
                            print('warning bad read at line:', lineNum , file=sys.stderr)
                        else:
                            yield header, sequence, qDesc, qScore #yield all fastQ fields of current read
                        header = ''; sequence = ''; qDesc = ''; qScore = '' #reset fields