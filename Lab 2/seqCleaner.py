#!/usr/bin/env python3 
# Name: Derfel Terciano (dtercian) 
# Group Members: None

'''
Logistics: make sure it contain docstrings and specific line or block comments that explain the semantics of your implementation.
Read a DNA string from user input and return a collapsed substring of embedded Ns to: {count}.

Example: 
 input: AaNNNNNNGTC
output: AA{6}GTC

Any lower case letters are converted to uppercase
'''

class DNAstring (str):
    def length (self):
        return (self.length()) #edit this
    
    def purify(self):
        ''' Return an upcased version of the string, collapsing a single run of Ns.'''
        upperString = self.upper()
        firstInstanceOfN = upperString.index('N')
        counter = upperString.count('N')
        upperString = upperString.replace('N','')
        upperString = upperString[:firstInstanceOfN] + '{' + f'{counter}' + '}' +upperString[firstInstanceOfN:]
        return upperString
    
def main():
    ''' Get user DNA data and clean it up.'''
    data = input('DNA data?')
    thisDNA = DNAstring (data)
    pureData = thisDNA.purify()
    print (pureData)
    
main()