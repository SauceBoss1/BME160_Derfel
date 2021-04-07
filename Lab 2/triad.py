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
    Original Author: David Bernick \n
    Improved by: Derfel Terciano \n
    Date: April 6, 2021\n
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
        '''returns an angle by taking in three side lengths and returning a radian'''
        theta = math.acos(((a**2)+(b**2)-(c**2))/(2*a*b))
        return theta

    def __bondLength(self,a,b):
        '''Calculates the bond length from two 3-elemented tuples'''
        #this function uses a for loop in order to calculate the summation of the distance between two 3D points
        summationOfPoints = 0
        for i in range(len(a)):
            summationOfPoints += (a[i]-b[i])**2
        return math.sqrt(summationOfPoints)
    
    def __angleCalculation(self,a,b,c):
        '''
        calculates the angle of a tuple by calculating 3 bond lengths between 3 given points then returning the angle (in radians) using the law of cosines
        '''

        # the 3 lines below calculate the bond lengths of the 3 given points
        lengthA = self.__bondLength(a,b)
        lengthB = self.__bondLength(b,c)
        lengthC = self.__bondLength(a,c)

        #the line below calculates the angle (in radians) based on the calculated bond lengths above
        theta = self.__lawOfCosineAngle(lengthA,lengthB,lengthC )
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