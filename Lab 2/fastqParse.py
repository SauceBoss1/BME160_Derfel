#!/usr/bin/env python3 
# Name: Derfel Terciano (dtercian)
# Group Members: None
'''Sequence information parsing
In this exercise, you will create a program to “parse” sequence name information from a single line of a FASTQ formatted file. Your task is to create a Python script called fastqParse that:

asks for and collects the seqname line of a FASTQ file using input()
parses out each field of the run information from the string and displays each of them on a new line For example, if I enter the FASTQ seqname line:
@EAS139:136:FC706VJ:2:2104:15343:197393 then the program will output:
Instrument = EAS139
Run ID = 136
Flow Cell ID = FC706VJ
Flow Cell Lane = 2
Tile Number = 2104
X-coord = 15343
Y-coord = 197393

Hints:
The input string is guaranteed to have 7 fields.
The first character of the FASTQ seqname line is “@” and each field of the run information is separated by a colon“:”.
A reasonable solution would be around 16 lines of code excluding comments.
To get full credit on this assignment, your code needs to:

Run properly (execute and produce the correct output)
Contain documentation/comments
Include any assumptions or design decisions you made in writing your code
Include an overview describing what your program does with expected inputs and outputs
'''

'''
Program docstring goes here.
'''

class FastqString (str):
    ''' Class docstring goes here.'''
    def parse(self):
        ''' Method docstring goes here.'''
        qString = self
        qString = qString.replace('@','')
        qString = qString.split(':')


        instrument = 'Intrument = '+ f'{qString[0]}' +'\n'
        runID = 'Run ID = ' + f'{qString[1]}' +'\n'
        flowCellID = 'Flow Cell ID = ' + f'{qString[2]}' +'\n'
        flowCellLane = 'Flow Cell Lane = ' + f'{qString[3]}' +'\n'
        tileNumber = 'Tile Number = ' + f'{qString[4]}' +'\n'
        xCoord = 'X-coord = ' + f'{qString[5]}' +'\n'
        yCoord = 'Y-coord = ' + f'{qString[6]}'
        return instrument + runID + flowCellID + flowCellLane + tileNumber + xCoord + yCoord
def main():
    ''' Function docstring goes here.'''
    fastQ = FastqString(input('FASTQ Seq name: '))
    print(fastQ.parse())

main()