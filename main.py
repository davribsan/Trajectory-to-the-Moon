'''
(a) Interplanetary trajectory to the Moon
(b) Version number: 1
(c) Autors: davribsan
(d) Date of initializaition: 19/10/25
(e) Description of the program: 
    Design the interplanetary trajectory for a spacecraft from a LEO orbit to an specific orbit 
    in the Moon minimizing Δv budget. The solution must provide the total Δv. Adittionally, in 
    order to minimize Δv, the plane inclination in LEO is the same as the final one in the Moon 
    leading to a 2D case.   
(f) Lambert sover: Izzo 
    References:
    [1] Izzo, D. (2015). Revisiting Lambert's problem. Celestial Mechanics
           and Dynamical Astronomy, 121(1), 1-15.

    [2] Lancaster, E. R., & Blanchard, R. C. (1969). A unified form of
           Lambert's theorem (Vol. 5368). National Aeronautics and Space
           Administration.
       
g) Range of validity expected of the parameters and range tested:
   - Expected: An average error in the velocity vector of 10^-13 is expected, with a maximum error of 10-8.
   For the single revolution case, the expected average of iterations is 2.1 iterations while, in the 
   multiple revolution case 3.3. No convergence when the number of iterations is too small. Aditionally, 
   it should work properly with elliptical orbits.
   - Tested: It has been tested for an interplanetary Earth-Moon trajectory in 2D with a vehicle of 30 tones, 
   assuming the mean distance between the bodies, considering only one gravitational field at a time, a 
   departure point on the side of the Earth opposite the Moon along the line connecting the centres of the two 
   bodies, launcher constant mass and instant impulse. The final orbits were circular with an altitude of 490 
   km and inclination of 49º and a elliptical with 300 km of periapsis, 1658 km of apoapsis and 157º of inclination. 

(h) Inputs:

    G: Universal gravitational constant (constant)                                          [m3⋅kg-1⋅s-2]
    M_EARTH: Mass of the Earth (constant)                                                   [kg]
    M_MOON: Mass of the Moon (constant)                                                     [kg]
    R_EARTH: Earth radius (constant)                                                        [m]
    R_MOON: Moon radius (constant)                                                          [m]
    r_MOON: Distance Earth's centre - Moon's centre (constant)                              [m]
    V_EQ: Velocity on the Earth's surface at the equator (constant)                         [m/s]
    mass_launcher: Mass of the launcher                                                     [kg] 
    h_orbit: Altitude of the initial LEO orbit                                              [m]
    h_c: Altitude of the final circular orbit at Moon                                       [m]
    peri_m: Periapsis of the final elliptical orbit at Moon                                 [m]
    apo_m: Apoapsis of the final elliptical orbit at Moon                                   [m]
    r_sphere_influence: Radius of the Moon's sphere of influence                            [m]
    maxiter: Number of iterations the solver will do                                        [ ]
    atol: Absolute tolerance of the Lambert solver                                          [ ]
    rtol: Relative tolerance of the Lmabert solver                                          [ ]

(i) Outputs:

    Time of flight between LEO and the Moon's sphere of incluence                           [h]
    Time of flight between the boundary and the insertion point for final circular orbit    [h]
    Time of flight between the boundary and the insertion point for final elliptical orbit  [h]
    Total time of fligth for final circular orbit                                           [h]
    Total time of fligth for final elliptical orbit                                         [h]
    Δv for leaving LEO                                                                      [m/s]
    Δv passig trhought the sphere of influence for final circular orbit                     [m/s]
    Δv passig trhought the sphere of influence for final elliptical orbit                   [m/s]
    Δv to circularize the circular orbit                                                    [m/s]
    Δv to circularize the elliptical orbit                                                  [m/s]
    Δv total to get to the Moon's circular orbit                                            [m/s]
    Δv total to get to the Moon's elliptical orbit                                          [m/s]
    Trajectory to get to the Moon's circular orbit
    Trajectory to ger to the Moon's elliptical orbit
    
(j) List of dependencies:
    - This program requires the izzo Lambert solver from the solver folder
    - This program requires the orbital mechanics constants from the data folder
    - This program requires the inputs from the data folder
    - This program requires the orbital equations from the utils folder
    - This program requires the geometric functions from the utils folder

(k) Software version: Python 3.12.4
'''

# ------------------------------------------------------------------------------------
# Libraries 

import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from matplotlib.patches import Circle
from solvers import izzo
from utils.geometry import * 
from utils.orbital import *
from data.orbital_mechanics_constants import *
from data.inputs import *

# ------------------------------------------------------------------------------------
# Program

# 0. Save all results printed on the screen
sys.stdout = open("Trajectory_delta_v.txt", "w")

# 1. Initial orbital parameters (poiting to Moon)
r_sat = h_orbit + R_EARTH
a_initial = compute_semimajor_axis_from_positions(r_sat,r_MOON)
e_initial = compute_eccentricty_ellipsis(a_initial,r_sat)
b_initial = semiminor_axis(a_initial,e_initial)

# 2. Find intersection point to the Moon's sphere of influence  
intersection = intersections(a_initial,e_initial,r_sphere_influence)
p_intersection = intersection[0][2:]
p_intersection = np.concatenate((p_intersection, [0]))
p_intersection = p_intersection.reshape(1, 3)
r_intersection_point = orbit_radii_from_coordinates(p_intersection)
r_intersection_point = r_intersection_point[0]
t_anomaly_intersection_point = compute_true_anomaly_from_r(a_initial,e_initial,r_intersection_point)

# 3. Time of flight between LEO and the intersection point
# 3.1. Eccenctric anomaly
E1 = eccentric_anomaly(e_initial,0)
E2 = eccentric_anomaly(e_initial,t_anomaly_intersection_point)

# 3.2. Mean anomaly
M1 = mean_anomaly(E1,e_initial)
M2 = mean_anomaly(E2,e_initial)

# 3.3. Mean motion
mu_earth = G*M_EARTH 
n = mean_motion(mu_earth,a_initial)

# 3.4. Time of flight
tof_initial = time_of_flight(M1,M2,n)
print(f"\nTof 1: {tof_initial/3600} h")

# 4. Coordenates of initial and final true anomalies
p1 = compute_position(r_sat,0)
p_intersection = p_intersection[0]

# 5. Solve Lambert's problem
v1_, v2_ = izzo.izzo2015(mu_earth, p1, p_intersection, tof_initial, maxiter, atol, rtol)
v1 = np.linalg.norm(v1_)
v2 = np.linalg.norm(v2_)

# 6. Δv1 
v_LEO_circ = circular_veocity(mu_earth,r_sat)
v_LEO = ([0,v_LEO_circ,0])
delta_v1 = np.linalg.norm(v1_-v_LEO)
print(f"Δv1: {delta_v1} m/s\n")

# 7. Change the reference frame to the Moon
moon_position_from_Earth = ([-r_MOON,0,0])
relative_distance_moon_intersection = p_intersection-moon_position_from_Earth

# 8. Flight between sphere of influence boundary and Moon circular orbit
mu_moon = G*M_MOON
p_moon_c_orbit = np.array([-R_MOON-h_c, 0, 0], dtype=np.float64)
p_influnce_sphere = relative_distance_moon_intersection

delta_v2 = delta_v1
tof_optimized = 0
v3_optimized = 0
v4_optimized = 0

# 8.1. Test different time of flights and select the one that minimizes Δv2
for tof in range(10*3600,50*3600):
    
    v3_, v4_ = izzo.izzo2015(mu_moon, p_influnce_sphere, p_moon_c_orbit, tof, maxiter, atol, rtol)

    delta_v_candidate = np.linalg.norm(v3_-v2_)

    if delta_v_candidate < delta_v2:
        tof_optimized = tof
        delta_v2 = delta_v_candidate
        v3_optimized = v3_
        v4_optimized = v4_

print('Circular orbit')
print(f'∆v2: {delta_v2} m/s')

# 8.2. Orbit circularization
r_moon_c_orbit = R_MOON + h_c
v_moon_circ = circular_veocity(mu_moon,r_moon_c_orbit) 
v_moon_circ = np.array([0, -v_moon_circ, 0], dtype=np.float64)

# 8.2.1. Δv3 
delta_v3 = np.linalg.norm(v_moon_circ-v4_optimized)
print(f'∆v3: {delta_v3} m/s')
total_detlta_v = delta_v1 + delta_v2 + delta_v3
print(f'∆v total: {total_detlta_v} m/s')
print(f'Tof 2: {tof_optimized/3600} h')
print(f'Tof total: {(tof_initial+tof_optimized)/3600} h\n')

# 9. Flight between sphere of influence boundary and Moon elliptical orbit
delta_v2 = delta_v1
tof_optimized = 0
v3_optimized = 0
v4_optimized = 0
p_moon_e_orbit = np.array([-R_MOON-peri_m, 0, 0], dtype=np.float64)

# 9.1. Test different time of flights and select the one that minimizes Δv2
for tof in range(10*3600,50*3600): 
    
    v3_, v4_ = izzo.izzo2015(mu_moon, p_influnce_sphere, p_moon_e_orbit, tof, maxiter, atol, rtol)
    delta_v_candidate = np.linalg.norm(v3_-v2_)
    if delta_v_candidate < delta_v2:
        tof_optimized = tof
        delta_v2 = delta_v_candidate
        v3_optimized = v3_
        v4_optimized = v4_

print('Elliptical orbit:')
print(f'∆v2: {delta_v2} m/s')

# 9.2. Orbit circularization (insertion in the periapsis to save Δv)
r_peri_moon = peri_m + R_MOON
r_apo_moon = apo_m + R_MOON
a_moon = compute_semimajor_axis_from_positions(r_peri_moon,r_apo_moon)
v_moon_peri = velocity_periapsis(mu_moon,a_moon,r_peri_moon)
v_moon_peri = np.array([0, -v_moon_peri, 0], dtype=np.float64)


# 9.2.1. Δv3 
delta_v3 = np.linalg.norm(v_moon_peri-v4_optimized)
print(f'∆v3: {delta_v3} m/s')
total_detlta_v = delta_v1 + delta_v2 + delta_v3
print(f'∆v total: {total_detlta_v} m/s')
print(f'Tof 2: {tof_optimized/3600} h')
print(f'Tof total {(tof_initial+tof_optimized)/3600} h\n')

# 10. Save results in a txt file
sys.stdout.close()
sys.stdout = sys.__stdout__

# 11. Construct trajectory
# 11.1. Segment 1
w0 = np.concatenate([p1, v1_]) 
t_span = [0, tof_initial]
sol = solve_ivp(lambda t, w: two_body(t, w, mu_earth), t_span, w0, rtol=1e-9, atol=1e-9, max_step=10)
x_segment_1 = sol.y[0]  
y_segment_1 = sol.y[1]

# 11.2. Segment 2
w0 = np.concatenate([p_influnce_sphere, v3_optimized]) 
t_span = [0, tof_optimized]
sol = solve_ivp(lambda t, w: two_body(t, w, mu_moon), t_span, w0, rtol=1e-9, atol=1e-9, max_step=10)
x_segment_2 = sol.y[0] - r_MOON # Return to the Earth's reference frame 
y_segment_2 = sol.y[1]

# 11.3. Segement 3, final circular orbit
T_orbit = orbital_period(r_moon_c_orbit,mu_moon)
w0 = np .concatenate([p_moon_c_orbit, v_moon_circ]) 
t_span = [0, T_orbit]
sol = solve_ivp(lambda t, w: two_body(t, w, mu_moon), t_span, w0, rtol=1e-9, atol=1e-9, max_step=10)
x_segment_3 = sol.y[0] - r_MOON  
y_segment_3 = sol.y[1]

# 11.3.1. Earth
angle = np.linspace(0, 2*np.pi, 200)
x_circle = R_EARTH * np.cos(angle)
y_circle = R_EARTH * np.sin(angle)
plt.fill(x_circle, y_circle, color='b', alpha=1, label='Earth')

# 11.3.2. LEO initial orbit
LEO_circle = Circle((0,0), radius=r_sat, fill=False, edgecolor='magenta', linewidth=1, label='LEO')
plt.gca().add_patch(LEO_circle)

plt.plot(x_segment_1,y_segment_1, color='green', label='Segment 1')
plt.plot(x_segment_2,y_segment_2, color='darkorange', label='Segment 2')
plt.plot(x_segment_3,y_segment_3,color='black', label='Moon orbit')

# Markers for visualization
x_soi = -3.34234085e+08
y_soi = 3.35164688e+07
x_soi_2 = -3.85e8
y_soi_2 = 6e7
plt.plot(x_soi, y_soi, 'o', markerfacecolor='none', markeredgecolor='red', markersize=20)
plt.text(x_soi_2, y_soi_2, 'Beginning of the Sphere of Influence', color='red', fontsize=15, ha='left', va='bottom')

plt.legend(loc = 'best')
plt.title('Preliminary interplanetary trajectory to circular orbit', fontweight='bold', fontsize=20)
plt.xlabel('X[m]', fontsize=18)
plt.ylabel('Y[m]', fontsize=18)
plt.axis("equal")
plt.show()

# 11.4. Segment 3, final elliptical orbit
T_orbit = orbital_period(a_moon,mu_moon)
w0 = np .concatenate([p_moon_e_orbit, v_moon_peri])
t_span = [0, T_orbit]
sol = solve_ivp(lambda t, w: two_body(t, w, mu_moon), t_span, w0, rtol=1e-9, atol=1e-9, max_step=10)
x_segment_3 = sol.y[0] - r_MOON  
y_segment_3 = sol.y[1]

# 11.4.1. Earth
angle = np.linspace(0, 2*np.pi, 200)
x_circle = R_EARTH * np.cos(angle)
y_circle = R_EARTH * np.sin(angle)
plt.fill(x_circle, y_circle, color='b', alpha=1, label='Earth')

# 11.4.2. LEO initial orbit
LEO_circle = Circle((0,0), radius=r_sat, fill=False, edgecolor='magenta', linewidth=1, label='LEO')
plt.gca().add_patch(LEO_circle)

plt.plot(x_segment_1,y_segment_1, color='green', label='Segment 1')
plt.plot(x_segment_2,y_segment_2, color='darkorange', label='Segment 2')
plt.plot(x_segment_3,y_segment_3,color='black', label='Moon orbit')

# Markers for visualization
x_soi = -3.34234085e+08
y_soi = 3.35164688e+07
x_soi_2 = -3.85e8
y_soi_2 = 6e7
plt.plot(x_soi, y_soi, 'o', markerfacecolor='none', markeredgecolor='red', markersize=20)
plt.text(x_soi_2, y_soi_2, 'Beginning of the Sphere of Influence', color='red', fontsize=15, ha='left', va='bottom')

plt.legend(loc = 'best')
plt.title('Preliminary interplanetary trajectory to elliptical orbit', fontweight='bold', fontsize=20)
plt.xlabel('X[m]', fontsize=18)
plt.ylabel('Y[m]', fontsize=18)
plt.axis("equal")
plt.show()




