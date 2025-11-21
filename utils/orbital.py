'''
(a) Geometric functions for an interplanetary trajectory to the Moon
(b) Version number: 1
(c) Autors: davribsan
(d) Date of initializaition: 19/10/25
(e) Inputs:

   t: Instant of time                                                   [s]
   w: State vector of a body in the space                               [m,m,m,m/s,m/s,m/s]
   mu: Standard gravitational parameter                                 [m^3,s^2]
   a: Semimajor axis of the orbit                                       [m]
   rp: Radius at periapsis                                              [m]
   e: Eccentricity of the orbit                                         [ ]
   r: Radius of the orbit at a specific position                        [m]
   ra: Radius at apoapsis                                               [m]
   Mi: Mean anomaly at position in the orbit                            [rad]
   n: Mean motion of an orbit                                           [rad/s]
   theta: True anomaly of the orbit at a specific position              [º]
   E: Eccentric anomaly of the orbit at a specific position             [rad]
   v: Linear velocity of the body at time t                             [m/s]

(f) Outputs:

    w_dot: State derivative vector of a body in the space               [m/s,m/s,m/s,m/s^2,m/s^2,m/s^2]
    v: Velocity at periapsis                                            [m/s] 
    theta: True anomaly of the orbit at a specific position             [º]
    e: Eccentricity of the orbit                                        [ ]
    a: Semimajor axis of the orbit                                      [m]
    b: Semiminor axis of the orbit                                      [m]
    T: Period of an orbit                                               [s]
    tof: Time of flight between true anomalies                          [s]
    n: Mean motion of an orbit                                          [rad/s]
    E: Eccentric anomaly of the orbit at a specific position            [rad]
    M: Mean anomaly at position in the orbit                            [rad]
    vc: orbital velocity of a circular orbit                            [m/s]

(g) Software version: Python 3.12.4
'''

# ------------------------------------------------------------------------------------
# Libraries 

import math
import numpy as np

# ------------------------------------------------------------------------------------
# Functions

def two_body(t, w, mu):
    "Derivative of the state vector w for a two-body problem"
    r = w[:3]
    v = w[3:]
    r_norm = np.linalg.norm(r)
    a = -mu * r / r_norm**3
    w_dot = np.concatenate([v, a])
    return w_dot

def velocity_periapsis(mu,a,rp):
    "Velocity at periapsis according to vis-viva equation"
    v = math.sqrt(mu*((2/rp) - (1/a)))
    return v

def compute_true_anomaly_from_r(a,e,r):
    "Computes the true anomaly in orbit given the radius"
    theta = math.acos(((a*(1-e**2)/r)-1)/e)
    return math.degrees(theta)

def compute_eccentricty_ellipsis(a,rp):
    "Computes the orbit eccentricuty from the semimajor axis and the periapsis"
    e = 1-(rp/a)
    return e

def compute_semimajor_axis_from_positions(rp,ra):
    "Computes the semimajor axis of the orbit from periapsis and apoapsis"
    a = (rp+ra)/2
    return a

def semiminor_axis(a,e):
    "Computes the semiminor axis of the orbit"
    b = a*math.sqrt(1-e**2)
    return b

def orbital_period(a,mu):
    "Computes the orbital period"
    T = 2*math.pi*math.sqrt(a**3/mu)
    return T

def time_of_flight(M1,M2,n):
    "Time of flight between true anomalies"
    tof = (M2-M1)/n
    return tof

def mean_motion(mu,a):
    "Computes the mean motion in a Keplerian orbit"
    n = math.sqrt(mu/(a**3))
    return n

def eccentric_anomaly(e,theta):
    "Computes the eccentric anomaly in a Keplerian orbit for a given true anomaly"
    theta = math.radians(theta)
    E = 2.0 * math.atan2(math.sqrt(1 - e) * math.sin(theta/2),
                         math.sqrt(1 + e) * math.cos(theta/2)) 
    return E

def mean_anomaly(E,e):
    "Computes the mean anomaly in a Keplerain orbit corresponding to certain true anomaly"
    M = E-e*math.sin(E)
    return M

def compute_semimajor_axis(r,v,mu):
    "Computes the semimajor axis of the orbit"
    a = (r*mu)/(2*mu-r*(v**2))
    return a

def circular_veocity(mu,r):
    "Computes the orbital velocity of a circular orbit"
    vc = math.sqrt(mu/r)
    return vc
