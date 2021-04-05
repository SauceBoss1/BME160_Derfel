#!/usr/bin/env python3 
# Name: Your full name  
# Group Members: List full names or “None”

'''
Read a DNA string from user input and return a collapsed substring of embedded Ns to: {count}.

Example: 
 input: AaNNNNNNGTC
output: AA{6}GTC

Any lower case letters are converted to uppercase
'''

class DNAstring (str):
    def length (self):
        return (length(self))
    
    def purify(self):
        ''' Return an upcased version of the string, collapsing a single run of Ns.'''
        pass # this is just a placeholder, so dont leave "pass" here
    
def main():
    ''' Get user DNA data and clean it up.'''
    data = input('DNA data?')
    thisDNA = DNAstring (data)
    pureData = thisDNA.purify()
    print (pureData)
    
main()