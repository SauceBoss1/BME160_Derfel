#!/usr/bin/env python3 
# Name: Derfel Terciano (dtercian)
# Group Members: None

# For extra-credit, provide a direct replacement for the Triad class. 
# The external methods that calculate angles, distances, and points (tuples) p,q and r must be maintained such that either version of the Triad class can be used.

# You could use the cosine law to calculate angles instead of the dot product. 
# You might make use of the numpy module. You might recode each of the methods to avoid using zip. You might consider using list iterations.
# Your Triad replacement must reimplement all of Triad public function, without using zip and without being a trivial rewrite.
# Your implementation need not be as compact as the current implementation, and it needs to be correct and fully documented to receive full credit. 
import math

class Triad:
    '''
    TODO: ADD DOCSTRING DOCUMENTATION
    '''
    def __init__(self,p,q,r) :
        """ Construct a Triad.  
        
        Example object construction:
            p1 = Triad( p=(1.,0.,0.), q=(0.,0.,0.), r=(0.,0.,0.) ). 
        """
        self.p = p
        self.q = q
        self.r = r
#private helper methods
        
    def __lawOfCosineAngle(self,a,b,c):
        '''returns an angle by taking in three points and returning a radian'''
        x = math.acos(((a**2)+(b**2)-(c**2))/(2*a*b))
        return x

    def __bondLength(self,a,b):
        summationOfPoints = 0
        for i in range(len(a)):
            summationOfPoints += (a[i]-b[i])**2
        return math.sqrt(summationOfPoints)
    
    def __angleCalculation(self,a,b,c):
        lengthA = self.__bondLength(a,b)
        lengthB = self.__bondLength(b,c)
        lengthC = self.__bondLength(a,c)
        theta = self.__lawOfCosineAngle(lengthA,lengthB,lengthC)
        return theta 

#Calculates lengths(distances) of PQ, PR, QR
    def dPQ (self):
        """ Provides the distance between point p and point q """
        return self.__bondLength(self.p,self.q)
    
    def dPR (self):
        """ Provides the distance between point p and point r """
        return self.__bondLength(self.p,self.r)
    
    def dQR (self):
        """ Provides the distance between point q and point r """
        return self.__bondLength(self.q,self.r)
    
    def angleP (self) :
        """ Provides the angle made at point p by segments pq and pr (radians). """
        return self.__angleCalculation(self.q,self.p,self.r)
    
    def angleQ (self) :
        """ Provides the angle made at point q by segments qp and qr (radians). """
        return self.__angleCalculation(self.p,self.q,self.r)
        
    def angleR (self) :
        """ Provides the angle made at point r by segments rp and rq (radians). """
        return self.__angleCalculation(self.p,self.r,self.q)


def convertToCoordinates (coordInput):
    '''This function converts the string of inputs into 3 lists of floats that the Triad class can then use

        input: a string
        output: 3 float lists
    '''
    coordInput = coordInput.replace(' ','').replace('(','').replace(')','') # removes any unnecessary data like spaces and parenthesis

    # the following line removes all idenitfiers and replaces 'N=' and 'Ca=' with ':' (Any unique character that isn't related to the string input will work)
    # I used the approach of inserting unique character because I need to eventually split the string into a 3-element list that
    #   correspond to the triad class inputs
    coordInput = coordInput.replace('C=', '').replace('N=', ':').replace('Ca=', ':')
    # using the example string above, the current result would look like the following string (keep in mind, the following is still one string):
    #   '39.447, 94.657, 11.824:39.292, 95.716, 11.027:39.462, 97.101, 11.465'

    # I create the 3-element list in the below line where each element will correspond to the a certain input of the Triad class
    coordInput = coordInput.split(':')
    # now the example string would look like the following 3-element list:
    #   ['39.447, 94.657, 11.824','39.292, 95.716, 11.027','39.462, 97.101, 11.465']


    # the 3 lines below turns each element into its own individual list and is then formatted into a float
    # in terms of the triad class inputs: coordC = p, coordN = q, coordCa = r
    coordC = list(map(float,coordInput[0].split(',')))
    coordN = list(map(float,coordInput[1].split(',')))
    coordCa = list(map(float,coordInput[2].split(',')))

    return coordC, coordN, coordCa


def main():
    '''asks for 3 atomic coordinates and outputs bond lengths and angles

        -input: string
        -outputs: bond lengths and bond angles
        Example: C = (39.447, 94.657, 11.824) N = (39.292, 95.716, 11.027) Ca = (39.462, 97.101, 11.465)
    '''
    coordInput = input('Enter coordinates: ')
    #coordInput = 'C = (39.447, 94.657, 11.824) N = (39.292, 95.716, 11.027) Ca = (39.462, 97.101, 11.465)'
    p,q,r = convertToCoordinates(coordInput) #since the convertToCoordinates function returns 3 lists, 3 names are assigned simultaneously to the function 
    bondInfo = Triad(p,q,r) #the Triad class is constructed

    # the 3 lines below formats the calculations to the correst number of sig figs 
    ncBondLength = 'N-C bond length = ' +f'{bondInfo.dPQ():.3}\n'
    ncaBondLength = 'N-Ca bond length = ' + f'{bondInfo.dQR():.3}\n'
    nBondAngle = 'C-N-Ca bond angle = ' +f'{float(math.degrees(bondInfo.angleQ())):.4}'
    print (f'{ncBondLength}{ncaBondLength}{nBondAngle}') #prints results of calculations

main()