#!/usr/bin/env python3 
# Name: Derfel Terciano (dtercian)
# Group Members: None

# For extra-credit, provide a direct replacement for the Triad class. 
# The external methods that calculate angles, distances, and points (tuples) p,q and r must be maintained such that either version of the Triad class can be used.

# You could use the cosine law to calculate angles instead of the dot product. 
# You might make use of the numpy module. You might recode each of the methods to avoid using zip. You might consider using list iterations.
# Your Triad replacement must reimplement all of Triad public function, without using zip and without being a trivial rewrite.
# Your implementation need not be as compact as the current implementation, and it needs to be correct and fully documented to receive full credit. 

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
    def dPQ (self):
        """ Provides the distance between point p and point q """
        pass
    
    def dPR (self):
        """ Provides the distance between point p and point r """
        pass
    
    def dQR (self):
        """ Provides the distance between point q and point r """
        pass
    
    def angleP (self) :
        """ Provides the angle made at point p by segments pq and pr (radians). """
        pass
    
    def angleQ (self) :
        """ Provides the angle made at point q by segments qp and qr (radians). """
        pass
        
    def angleR (self) :
        """ Provides the angle made at point r by segments rp and rq (radians). """
        pass