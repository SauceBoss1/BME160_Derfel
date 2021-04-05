#!/usr/bin/env python3 
# Name: Derfel Terciano (dtercian) 
# Group Members: None

'''
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
        return "this is empty" #edit this as well
    
def main():
    ''' Get user DNA data and clean it up.'''
    data = input('DNA data?')
    thisDNA = DNAstring (data)
    pureData = thisDNA.purify()
    print (pureData)
    
main()