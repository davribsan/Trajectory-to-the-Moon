'''
(a) Script containing the inputs for an interplanetary trajectory to the Moon

(b) Inputs:

    mass_launcher: Mass of the launcher                                 [kg] 
    h_orbit: Altitude of the initial LEO orbit                          [m]
    h_c: Altitude of the final circular orbit at Moon                   [m]
    peri_m: Periapsis of the final elliptical orbit at Moon             [m]
    apo_m: Apoapsis of the final elliptical orbit at Moon               [m]
    r_sphere_influence: Radius of the Moon's sphere of influence        [m]
    maxiter: Number of iterations the solver will do                    [ ]
    atol: Absolute tolerance of the Lambert solver                      [ ]
    rtol: Relative tolerance of the Lambert solver                      [ ]

'''

# ------------------------------------------------------------------------------------

mass_launcher = 30e3 
h_orbit = 200e3
h_c = 490e3
peri_m = 300e3 
apo_m = 1658e3
maxiter = 50              
atol = 1e-7                
rtol = 1e-7  
r_sphere_influence = 60000e3 