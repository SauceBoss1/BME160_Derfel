#!/usr/bin/env python3 
# Name: Derfel Terciano (dtercian) 
# Group Members: None


# Hints:

# The input sequence is not guaranteed to be uppercase, but should be interpreted as though it is all uppercase.
# Only the letters (A,C,G,T,N) will be included in the input.
# Only the DNA characters (A,C,G,T) should remain after the cleanup.
# The input may include, at most, one block of 'N' characters.

# To get full credit on this assignment, your code needs to:

# Run properly (execute and produce the correct output)
# Contain docstrings and specific line or block comments that explain the semantics of your implementation.
# Include any assumptions or design decisions you made in writing your code
# Include an overview describing what your program does with expected inputs and outputs. This should be in the form of a program level dosctring



'''
Read a DNA string from user input and return a collapsed substring of embedded Ns to: {count}.

Example: 
 input: AaNNNNNNGTC
output: AA{6}GTC

Any lower case letters are converted to uppercase

Assumptions:
    -Only (A,C,G,T,N) are allowed characters
    -An input without N (i.e. a nucleotide sequence that does not contain the character N) will not be tested
'''

class DNAstring (str):
    '''
    determines length of string and 'cleans up' given nucleotide sequence

    input: string
    output/return values:
        -method length(): length of given string
        -method purify(): returns an upcased version of the string, collapsing a single run of Ns.
    '''
    def length (self):
        '''This returns the length of the inputted string'''
        return (self.length()) 
    
    def purify(self):
        ''' Return an upcased version of the string, collapsing a single run of Ns.'''
        '''
        -Note: There are many ways to insert the formatted { (number if Ns) } in our string. 
            I found that inserting the formatted output at the first instance of N to be the easiest since we are assuming there 
            is one block (or chain) of Ns AT MOST. (i.e. ATNNNNNNGCCCCTTT is a valid test case)
            If there were multiple blocks of Ns (i.e. ATNNNGCGCGNNNNTTAATATANNNN) then I would have used a different approach.
        '''
        upperString = self.upper() #makes all characters of string input uppercase (this is the equivalence of a caps lock on a keyboard)
        firstInstanceOfN = upperString.index('N') #this is important for the formatted output later
        counter = upperString.count('N') #counts the number of in the inputted string
        upperString = upperString.replace('N','') #removes all instanced of N in the string

         
        # The below line inserts { (number of Ns) } at the index of the first instance of N
        upperString = upperString[:firstInstanceOfN] + '{' + f'{counter}' + '}' +upperString[firstInstanceOfN:] 

        return upperString
    
def main():
    ''' Get user DNA data and clean it up.'''
    
    '''
    Asks for nucleotide data from user and returns a cleaned string with {Number of Ns} replacing the block of Ns from orginal input

    Input: a string of  characters
    Output: a collapsed substring of embedded Ns to: {count}.

    '''
    data = input('DNA data?')
    thisDNA = DNAstring (data)
    pureData = thisDNA.purify()
    print (pureData)
    
main()