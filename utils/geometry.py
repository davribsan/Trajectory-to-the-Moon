'''
(a) Geometric functions for an interplanetary trajectory to the Moon
(b) Version number: 1
(c) Autors: davribsan
(d) Date of initializaition: 19/10/25
(e) Inputs:

   r: Radius of the orbit at a specific position                        [m]
   theta: True anomaly of the orbit at a specific position              [º]
   positions: Array containing positions in an orbit                    [m,m,m]
   a: Semimajor axis of the orbit                                       [m]
   e: Eccentricity of the orbit                                         [ ]
   R: Radius of a spherere                                              [m]

(f) Outputs:

    p: Coordenates with the position in an orbit                       [m,m,m]                                                        [m,m,m]
    r: Array with the radii in the orbit associated to each position   [m,m,...,m]                                   [m]
    unique: List of touples with information related to the 
            intersection of an ellipsis and a sphere                   [º,m,m,m] 

(g) Software version: Python 3.12.4
'''

# ------------------------------------------------------------------------------------
# Libraries 

import math
import numpy as np

# ------------------------------------------------------------------------------------
# Functions

def compute_position(r,theta):
    "Computes the position in an orbit given the radius and true anomaly"
    x = r*math.cos(math.radians(theta))
    y = r*math.sin(math.radians(theta))
    z = 0
    p = np.array([x, y, z])
    return p.copy()

def orbit_radii_from_coordinates(positions):
    "Computes all the radii corresponding to each coordenate in the orbit"
    r = np.linalg.norm(positions, axis=1) 
    return r

def intersections(a,e,R):
    "Computes the intersection point between the Hohmann transfer orbit and the Moon's sphere of influence"
    p = a*(1-e*e)
    r_a = a*(1+e)
    A = 2*p*r_a*e + (r_a*r_a - R*R)*e*e
    B = 2*p*r_a + 2*(r_a*r_a - R*R)*e
    C0= p*p + (r_a*r_a - R*R)
    D = B*B - 4*A*C0
    sols = []
    if D < 0:
        return sols
    for sign in (+1,-1):
        if A == 0:
            if abs(B) > 1e-12:
                C = -C0/B
            else:
                continue
        else:
            C = (-B + sign*math.sqrt(abs(D)))/(2*A)
        if abs(C) <= 1+1e-12:
            C = max(min(C,1.0),-1.0)
            theta1 = math.acos(C)
            theta2 = 2*math.pi - theta1
            for theta in (theta1, theta2):
                r = p/(1+e*math.cos(theta))
                x = r*math.cos(theta)
                y = r*math.sin(theta)
                sols.append((theta, r, x, y))
    unique = []
    for s in sols:
        if not any(abs(s[0]-u[0])<1e-8 for u in unique):
            unique.append(s)
    return unique