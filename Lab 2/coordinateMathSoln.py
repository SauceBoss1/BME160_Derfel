#!/usr/bin/env python3 
# Name: Derfel Terciano (dtercian)
# Group Members: None

#from triad import *
'''
coordinateMathSoln.py
This program takes in three sets of atomic coordinates, all provided on a single line. The program then calculates and returns the bond lengths and angles.

Example:
    -input: C = (39.447, 94.657, 11.824) N = (39.292, 95.716, 11.027) Ca = (39.462, 97.101, 11.465)
    -output: 
        N-C bond length = 1.33
        N-Ca bond length = 1.46
        C-N-Ca bond angle = 124.0 (in degrees)

Assumptions:
    -the input for the three coordinates are in the following order and uses only the following characters : C N Ca
    -all numbers are floats and all outputs are floats only
    -Each coordinate will contain only 3 numeric values.
    -all inputs will be formatted exactly like how the example input is 
'''

import math
class Triad :
    """
    Calculate angles and distances among a triad of points.
 
    Author: David Bernick
    Date: March 21, 2013
    Points can be supplied in any dimensional space as long as they are consistent.
    Points are supplied as tupels in n-dimensions, and there should be three
    of those to make the triad. Each point is positionally named as p,q,r
    and the corresponding angles are then angleP, angleQ and angleR.
    Distances are given by dPQ(), dPR() and dQR()
 
    Required Modules: math
    initialized: 3 positional tuples representing Points in n-space
             p1 = Triad( p=(1,0,0), q=(0,0,0), r=(0,1,0) )
    attributes: p,q,r the 3 tuples representing points in N-space
    methods:  angleP(), angleR(), angleQ() angles measured in radians
          dPQ(), dPR(), dQR() distances in the same units of p,q,r
 
    """
 
    def __init__(self,p,q,r) :
        """ Construct a Triad. 
        
        Example construction:
            p1 = Triad( p=(1.,0.,0.), q=(0.,0.,0.), r=(0.,0.,0.) ). 
        """
        self.p = p
        self.q = q
        self.r = r
# private helper methods
    def d2 (self,a,b) : # calculate squared distance of point a to b
        return float(sum((ia-ib)*(ia-ib)  for  ia,ib in zip (a,b)))
    
    def dot (self,a,b) : # dotProd of standard vectors a,b
        return float(sum(ia*ib for ia,ib in zip(a,b)))
    
    def ndot (self,a,b,c) : # dotProd of vec. a,c standardized to b
        return float(sum((ia-ib)*(ic-ib) for ia,ib,ic in zip (a,b,c)))
    
# calculate lengths(distances) of segments PQ, PR and QR
    def dPQ (self):
        """ Provides the distance between point p and point q """
        return math.sqrt(self.d2(self.p,self.q))
    
    def dPR (self):
        """ Provides the distance between point p and point r """
        return math.sqrt(self.d2(self.p,self.r))
    
    def dQR (self):
        """ Provides the distance between point q and point r """
        return math.sqrt(self.d2(self.q,self.r))
    
    def angleP (self) :
        """ Provides the angle made at point p by segments pq and pr (radians). """
        return math.acos(self.ndot(self.q,self.p,self.r) /   math.sqrt(self.d2(self.q,self.p)*self.d2(self.r,self.p)))
    
    def angleQ (self) :
        """ Provides the angle made at point q by segments qp and qr (radians). """
        return math.acos(self.ndot(self.p,self.q,self.r) /  math.sqrt(self.d2(self.p,self.q)*self.d2(self.r,self.q)))
 
    def angleR (self) :
        """ Provides the angle made at point r by segments rp and rq (radians). """
        return math.acos(self.ndot(self.p,self.r,self.q) /  math.sqrt(self.d2(self.p,self.r)*self.d2(self.q,self.r)))


def convertToCoordinates (coordInput):
    '''This function converts the string of inputs into 3 lists of floats that the Triad class can then use

        -input: a string
        -output: 3 float lists
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