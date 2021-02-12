# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 21:38:03 2020

@author: manoj
"""


import numpy as np
import PyDEM
import MDEM_3d
import matplotlib.pyplot as plt
import matplotlib
import os
import logging
import time
import math
time_initial = time.time()
matplotlib.use('Qt5Agg')

test_name = input('Name of the Simulation:')
output_forlder_name = test_name + '_output'
os.mkdir(output_forlder_name)
original_dir = os.getcwd()
os.chdir(output_forlder_name)
current_dir = os.getcwd()
os.mkdir('particle_data')
os.chdir('particle_data')
save_par_dir = current_dir+'\particle_data'
log_filename = 'log_filename.log'
for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)
logging.basicConfig(filename=log_filename,level=logging.INFO)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Input parameters

Total_particles_for_consolidation = 100
Time_for_consolidation = 5
Timestep_for_consolidation = 1e-5
rmax = 0.02
rmin = 0.005
Radius_of_particles = [rmin,rmax]
Density_of_particles = 2203
Youngs_modulus_of_particle = 7.2e9
Shear_modulus_of_particle = 31.3e9
Poissons_ratio_of_paticle = 0.01
Coeficient_of_restitution_of_particle = 0.3#
friction_coef_particle = 0.2
Poissons_ratio_of_wall = 0.01
Youngs_modulus_of_wall = 7.2e9
Shear_modulus_of_wall = 31.3e9
Bottom_wall_position = 0.25
Top_wall_position = 50+0.25
Cylinder_radius = .2
mass_scale = 1
accleration_due_to_gravity = 9.81
artificial_rolling_damping = 0
mu_r = 0.02

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

"""
Computation part
"""
P = np.array(0)
N_tot = Total_particles_for_consolidation
N = 0

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_overlap_hc = np.zeros([N_tot,N_tot,2])

T_overlap_z = np.zeros(N_tot)
T_overlap_x = np.zeros(N_tot)

T_overlap_horizontal = np.zeros(N_tot)
T_overlap_vertical = np.zeros(N_tot)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iteration = 0
times = 0
timestep = Timestep_for_consolidation
time_end = Time_for_consolidation

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

KE = []
Ti = []
t = []
u = []
v = []
w = []
x = []
y = []
z = []
ou = []
ow = []
ov = []
FZ = []
FY = []
vsv = []
Tq = []
sld = []
kec = []
KE_conv = 0
KE_conv_tol = 1e-3
disp_max = 0
dis = []
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P,E = PyDEM.ini_sim(Total_particles_for_consolidation,1)
P,E = PyDEM.Prop_gen(P,E,0,'all',Youngs_modulus_of_particle,Shear_modulus_of_particle,Poissons_ratio_of_paticle,
                     Coeficient_of_restitution_of_particle,friction_coef_particle,mu_r,Density_of_particles)
P = PyDEM.particle_gen_cyl(P,E,'all',[rmin,rmax],.5,Bottom_wall_position,Top_wall_position,Cylinder_radius)
tot_hf = []
Top_wall_position = PyDEM.place_top_plate(P)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""
Binning the domain
"""
resolution = 1.1
ngrid,x_array,y_array,z_array,del_d = PyDEM.gridit_cuboid(rmax,-0.35,0.35,-0.35,0.35,0,Top_wall_position,resolution)
Bin_array,Index_array = PyDEM.binit(P,ngrid,x_array,y_array,z_array,del_d)
Neiboughers = PyDEM.neibs_array(P,Bin_array,Index_array,ngrid,x_array,y_array,z_array)

t_overlap_hc = PyDEM.t_overlap_ip_init(P, Neiboughers)


t_overlap_cyl = PyDEM.t_overlap_wall(P)

t_overlap_bottom_wall = PyDEM.t_overlap_wall(P)


#wall1
nv1 = (-1,1,0)
p1 = (0,0,0)
hc_paras_wall1 = PyDEM.hertz_contact_parameters_wall(E,Youngs_modulus_of_wall,Shear_modulus_of_wall,Poissons_ratio_of_wall)
t_overlap_wall1 = PyDEM.t_overlap_wall(P)

#wall2
nv2 = (1,1,0)
p2 = (0,0,0)
hc_paras_wall2 = PyDEM.hertz_contact_parameters_wall(E,Youngs_modulus_of_wall,Shear_modulus_of_wall,Poissons_ratio_of_wall)
t_overlap_wall2 = PyDEM.t_overlap_wall(P)

#wall1
nv7 = (0,1,-1)
p7 = (0,0,0)
hc_paras_wall7 = PyDEM.hertz_contact_parameters_wall(E,Youngs_modulus_of_wall,Shear_modulus_of_wall,Poissons_ratio_of_wall)
t_overlap_wall7 = PyDEM.t_overlap_wall(P)

#wall2
nv8 = (0,1,1)
p8 = (0,0,0)
hc_paras_wall8 = PyDEM.hertz_contact_parameters_wall(E,Youngs_modulus_of_wall,Shear_modulus_of_wall,Poissons_ratio_of_wall)
t_overlap_wall8 = PyDEM.t_overlap_wall(P)

#wall2
nv9 = (0,1,0)
p9 = (0,0,0)
hc_paras_wall9 = PyDEM.hertz_contact_parameters_wall(E,Youngs_modulus_of_wall,Shear_modulus_of_wall,Poissons_ratio_of_wall)
t_overlap_wall9 = PyDEM.t_overlap_wall(P)


#wall1
nv3 = (0,0,1)
p3 = (0,0,-0.25)
hc_paras_wall3 = PyDEM.hertz_contact_parameters_wall(E,Youngs_modulus_of_wall,Shear_modulus_of_wall,Poissons_ratio_of_wall)
t_overlap_wall3 = PyDEM.t_overlap_wall(P)

#wall1
nv4 = (0,0,-1)
p4 = (0,0,0.25)
hc_paras_wall4 = PyDEM.hertz_contact_parameters_wall(E,Youngs_modulus_of_wall,Shear_modulus_of_wall,Poissons_ratio_of_wall)
t_overlap_wall4 = PyDEM.t_overlap_wall(P)

#wall1
nv5 = (-1,0,0)
p5 = (0.25,0,0)
hc_paras_wall5 = PyDEM.hertz_contact_parameters_wall(E,Youngs_modulus_of_wall,Shear_modulus_of_wall,Poissons_ratio_of_wall)
t_overlap_wall5 = PyDEM.t_overlap_wall(P)

#wall1
nv6 = (1,0,0)
p6 = (-0.25,0,0)
hc_paras_wall6 = PyDEM.hertz_contact_parameters_wall(E,Youngs_modulus_of_wall,Shear_modulus_of_wall,Poissons_ratio_of_wall)
t_overlap_wall6 = PyDEM.t_overlap_wall(P)



wall = [t_overlap_wall9]


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ite_delay_hc = 10
P_old = P

hc_paras_ip = PyDEM.hertz_contact_parameters_interparticle(E)
hc_paras_wall = PyDEM.hertz_contact_parameters_wall(E,Youngs_modulus_of_wall,Shear_modulus_of_wall,Poissons_ratio_of_wall)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while times <= 7:


   P = PyDEM.ini_ite(P)




   if iteration%ite_delay_hc == 0:
      Bin_array,Index_array = PyDEM.binit(P,ngrid,x_array,y_array,z_array,del_d)
      Neiboughers = PyDEM.neibs_array(P,Bin_array,Index_array,ngrid,x_array,y_array,z_array)
      t_overlap_hc = PyDEM.t_overlap_update(P, Neiboughers, t_overlap_hc)


   P,t_overlap_hc = PyDEM.Hertz_contact_interparticle(P,E,hc_paras_ip,timestep,t_overlap_hc,Neiboughers)



   # P,t_overlap_wall1,normal_force1 = PyDEM.Hertz_contact_arb_wall(P, E, nv1, p1 , hc_paras_wall1, timestep,t_overlap_wall1)
   # P,t_overlap_wall2,normal_force2 = PyDEM.Hertz_contact_arb_wall(P, E, nv2, p2 , hc_paras_wall2, timestep,t_overlap_wall2)
   P,t_overlap_wall3,normal_force3 = PyDEM.Hertz_contact_arb_wall(P, E, nv3, p3 , hc_paras_wall3, timestep,t_overlap_wall3)
   P,t_overlap_wall4,normal_force4 = PyDEM.Hertz_contact_arb_wall(P, E, nv4, p4 , hc_paras_wall4, timestep,t_overlap_wall4)
   P,t_overlap_wall5,normal_force5 = PyDEM.Hertz_contact_arb_wall(P, E, nv5, p5 , hc_paras_wall5, timestep,t_overlap_wall5)
   P,t_overlap_wall6,normal_force6 = PyDEM.Hertz_contact_arb_wall(P, E, nv6, p6 , hc_paras_wall6, timestep,t_overlap_wall6)
   # P,t_overlap_wall7,normal_force7 = PyDEM.Hertz_contact_arb_wall(P, E, nv7, p7 , hc_paras_wall7, timestep,t_overlap_wall7)
   # P,t_overlap_wall8,normal_force8 = PyDEM.Hertz_contact_arb_wall(P, E, nv8, p8 , hc_paras_wall8, timestep,t_overlap_wall8)
   P,t_overlap_wall9,normal_force9 = PyDEM.Hertz_contact_arb_wall(P, E, nv9, p9 , hc_paras_wall9, timestep,t_overlap_wall9)


   if times <= 5:
         p9 = (0,0.03*math.sin(70*times),0)

   P = PyDEM.timestep_update_constant_accleration(P,Total_particles_for_consolidation,timestep,[0,-10,0])

   KEt = PyDEM.total_KE(P)
   t.append(times)
   KE.append(KEt)
   if iteration%500 == 0:
       print(times)
       PC = PyDEM.comprehensive_array(P,E)
       PyDEM.save_to_csv(PC,nnew,iteration)
       logging.info("________________________________")
       logging.info("time elapsed = %f"%times)
       logging.info("Kinetic Energy = %f"%KEt)
       KE_conv = max(KE[-1000:len(KE)])
       disp_max = MDEM_3d.disp_array(P,P_old)
   dis.append(disp_max)
   kec.append(KE_conv)
   iteration = iteration + 1
   times = times + timestep


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
os.chdir(current_dir)
np.save('final_particle_data',P)
np.savetxt("final_particle_data.csv", P, delimiter=",")

plt.figure(1)
plt.plot(t,KE)
plt.savefig('KE_covergance.png',dpi = 2000)

P = PyDEM.comprehensive_array(P,E)
if Top_wall_position-Bottom_wall_position >= Cylinder_radius:
 	rlim = Top_wall_position-Bottom_wall_position
 	xlim = [-(rlim/2),(rlim/2)]
 	ylim = [0,rlim]
 	zlim = [-(rlim/2),(rlim/2)]
else:
 	rlim = Cylinder_radius
 	xlim = [-(rlim/2),(rlim/2)]
 	ylim = [0,rlim ]
 	zlim = [-(rlim/2),(rlim/2)]
MDEM_3d.plot_particle(P,xlim,ylim,zlim,rlim,Cylinder_radius,Bottom_wall_position,Top_wall_position)

Top_wall_position = PyDEM.place_top_plate(P)
Boundary_conditions = np.array([Cylinder_radius, Bottom_wall_position, Top_wall_position, Youngs_modulus_of_wall, Shear_modulus_of_wall, Poissons_ratio_of_wall])
np.save('boundary_condition_geometry',Boundary_conditions)

plt.savefig('final_particle_matrix.png',dpi = 600)
logging.shutdown()
os.chdir(original_dir)

print(MDEM_3d.packing_fraction(P,Cylinder_radius))
time_final = time.time()
print(time_final-time_initial)
