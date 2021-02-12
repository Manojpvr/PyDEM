#PyDEM v2.0(previously MDEM)
#author: P. V. R. Manoj
#Date: 15-12-2020
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
import numpy as np
import random
import os
import math
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Function to generate Particle array and Property array to initiate iteration

def ini_sim(no_par, no_type):
    #np is the number of particles and nt is the number of types of particles(particle materials)
    P = np.zeros([no_par, 18])#[0.Property_number, 1.Radius, 2.X, 3.Y, 4.Z, 5.Fx,
    #6.Fy, 7.Fz, 8.Tx, 9.Ty, 10.Tz, 11.Ux, 12.Uy, 13.Uz, 14.Wx, 15.Wy, 16.Wz, 17.Mass]
    E = np.zeros([no_type, 8])#[0.Youngs_modulus, 1.Shear_modulus, 2.Poissons ratio,
    #3.Coefficient_of_restitution, 4.Friction_coefficient, 5.Coeficient_of_rolling_resistance, 6.Density, 7.Beta_coef]
    return(P, E)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Property Generator
def Prop_gen(P, E, property_number, Par_range, Youngs_modulus, Shear_modulus, Poissons_ratio, Coeficient_of_restitution, Coeficient_of_friction, Coeficient_of_rolling_resistance, Density):
    no_par = np.shape(P)[0]
    if Par_range == 'all':
        Par_range = [par for par in range(no_par)]
    E[property_number, 0] = Youngs_modulus
    E[property_number, 1] = Shear_modulus
    E[property_number, 2] = Poissons_ratio
    E[property_number, 3] = Coeficient_of_restitution
    E[property_number, 4] = Coeficient_of_friction
    E[property_number, 5] = Coeficient_of_rolling_resistance
    E[property_number, 6] = Density
    E[property_number, 7] = 1.825742*math.log(Coeficient_of_restitution)/(math.sqrt((math.log(Coeficient_of_restitution))**2+9.8696))
    for idx in Par_range:
        P[idx, 0] = property_number
    return(P, E)

#Particle Generator
def particle_gen_cyl(P, E, Par_range, Radius_range, max_seed_vel, Cyl_base, Cyl_height, Cyl_radius):
    no_par = np.shape(P)[0]
    if Par_range == 'all':
        Par_range = [par for par in range(no_par)]
    else:
        no_par = len(Par_range)
    rmax = max(Radius_range)
    rmin = min(Radius_range)
    if rmin <= 0:
        raise ValueError("Radius of particle cannot be less than or equal to zero")
    if Cyl_radius <= 0:
        raise ValueError("Radius of cylinder cannot be less than or equal to zero")
    if Cyl_radius <= rmax:
        raise ValueError("Radius smaller than maximum allowed particle radius")
    if Cyl_height-Cyl_base <= rmax:
        raise ValueError("Cylinder height smaller than maximum allowed particle radius")
    r_sep = rmax*1.3
    nlp = int((Cyl_radius-r_sep)/(2*r_sep))#no. of particles in one radial extension
    nap = 8*nlp+1#no. of particles in one planar extension
    nhl = int(no_par/nap)+1#no. of vertical layers
    ang = 0
    h = Cyl_base+r_sep
    particle = 0
    par_idx = Par_range[particle]
    for idx in range(nhl):
        u_seed = -max_seed_vel + 2*max_seed_vel*random.random()
        w_seed = -max_seed_vel + 2*max_seed_vel*random.random()
        radius_seed = rmin + (rmax-rmin)*random.random()
        vol_par = (4/3)*3.14159*radius_seed**3
        density = E[int(P[par_idx, 0]), 6]
        mass_par = vol_par*density
        if particle < no_par-1:
            P[par_idx, 1:] = [radius_seed, 0, h, 0, 0, 0, 0, 0, 0, 0, u_seed, 0, w_seed, 0, 0, 0, mass_par]
            particle += 1
            par_idx = Par_range[particle]
        elif particle == no_par-1:
            P[par_idx, 1:] = [radius_seed, 0, h, 0, 0, 0, 0, 0, 0, 0, u_seed, 0, w_seed, 0, 0, 0, mass_par]
            particle += 1
        ang2 = ang
        for jdx in range(8):
            dia_sep = 2*r_sep
            for kdx in range(nlp):
                u_seed = -max_seed_vel + 2*max_seed_vel*random.random()
                w_seed = -max_seed_vel + 2*max_seed_vel*random.random()
                radius_seed = rmin + (rmax-rmin)*random.random()
                vol_par = (4/3)*3.14159*radius_seed**3
                density = E[int(P[par_idx, 0]), 6]
                mass_par = vol_par*density
                if particle < no_par-1:
                    P[par_idx, 1:] = [radius_seed, dia_sep*np.cos(ang2), h, dia_sep*np.sin(ang2), 0, 0, 0, 0, 0, 0, u_seed, 0, w_seed, 0, 0, 0, mass_par]
                    particle += 1
                    par_idx = Par_range[particle]
                    dia_sep = dia_sep+2*r_sep
                elif particle == no_par-1:
                    P[par_idx, 1:] = [radius_seed, dia_sep*np.cos(ang2), h, dia_sep*np.sin(ang2), 0, 0, 0, 0, 0, 0, u_seed, 0, w_seed, 0, 0, 0, mass_par]
                    particle += 1
            ang2 += np.pi/4
        ang += np.pi/6 + random.random()
        h += (2*r_sep)
    return(P)
def add_par(P,particle):
    P = np.row_stack((P,particle))
    return(P)
def remove_par(P,idx):
    P = np.delete(P,idx,axis = 0)
    return(P)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#From Particle array and Property array to comprehensive array"""

#List of attributes in comprehensive array"""

#0. Radius
#1. X Position
#2. Y Position
#3. Z Position
#4. X Orientation(not used)
#5. Y Orientation(not used)
#6. Z Orientation(not used)
#7. Youngs Modulus
#8. Shear Modulus
#9. Poissons Ration
#10. Coefficient of Restitution
#11. Traction zone tensile stress(not updated as of now)
#12. Mass
#13. Friction coeficient
#14. Initial Radius(not used as of now)
#15. U
#16. V
#17. W
#18. Ox
#19. Oy
#20. Oz
#21. Fx
#22. Fy
#23. Fz
#24. Tx
#25. Ty
#26. Tz


def comprehensive_array(P, E):
    no_par = np.shape(P)[0]
    P_comp = np.zeros([no_par, 27])
    for idx in range(no_par):
        prop_no = int(P[idx, 0])
        P_comp[idx] = [P[idx, 1], P[idx, 2], P[idx, 3], P[idx, 4], 0, 0, 0, E[prop_no, 0], E[prop_no, 1], E[prop_no, 2],
                       E[prop_no, 3], 0, P[idx, 17], E[prop_no, 4], P[idx, 1], P[idx, 11], P[idx, 12], P[idx, 13], P[idx, 14],
                       P[idx, 15], P[idx, 16], P[idx, 5], P[idx, 6], P[idx, 7], P[idx, 8], P[idx, 9], P[idx, 10]]
    return(P_comp)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Save comprehensive array
def prepend_line(file_name, line):
    #Insert given string as a new line at the beginning of a file
    # define name of temporary dummy file
    dummy_file = file_name + '.bak'
    # open original file in read mode and dummy file in write mode
    with open(file_name, 'r') as read_obj, open(dummy_file, 'w') as write_obj:
        # Write given line to the dummy file
        write_obj.write(line + '\n')
        # Read lines from original file one by one and append them to the dummy file
        for line in read_obj:
            write_obj.write(line)
    # remove original file
    os.remove(file_name)
    # Rename dummy file as the original file
    os.rename(dummy_file, file_name)

def save_to_csv(P, N, i):
    np.savetxt('at_ite_%d.xyz'%(i), P, delimiter=" ")
    prepend_line('at_ite_%d.xyz'%(i), '"particle data"')
    prepend_line('at_ite_%d.xyz'%(i), str(N))
    return()

def save_to_csv2(P, i):
    np.savetxt('at_ite_%d.txt'%(i), P, delimiter=", ")
    prepend_line('at_ite_%d.txt'%(i), 'r, x, y, z, orx, ory, orz, ym, sm, pr, cr, sxxo, mass, fc, ir, U, V, W, Ox, Oy, Oz, Fx, Fy, Fz, Tx, Ty, Tz')
    return()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Binning and neighbour finding

def gridit_cyl(par_radius_max, Cyl_radius, Cyl_height, Cyl_base, bin_size):#bin size as a fraction of maximum particle radius
    del_d = 2*par_radius_max*bin_size

    n_lin_array = int(2*Cyl_radius/del_d) + 1   #number of bins along horizontal axis
    n_lin_array_h = int((Cyl_height-Cyl_base)/del_d) + 1   #number of bins along vertical axis

    x_array = np.arange(-Cyl_radius, (n_lin_array)*del_d/2, del_d).tolist()
    z_array = np.arange(-Cyl_radius, (n_lin_array)*del_d/2, del_d).tolist()
    y_array = np.arange(Cyl_base, (n_lin_array_h)*del_d, del_d).tolist()

    n_array_x = len(x_array)
    n_array_y = len(y_array)
    n_array_z = len(z_array)

    n_grid = n_array_x*n_array_y*n_array_z

    return(n_grid, x_array, y_array, z_array,del_d)


def gridit_cuboid(par_radius_max, xmin, xmax, ymin, ymax, zmin, zmax, bin_size):#bin size as a fraction of maximum particle radius
    del_d = 2*par_radius_max*bin_size

    x_array = np.arange(xmin, xmax + del_d, del_d).tolist()
    z_array = np.arange(zmin, zmax + del_d, del_d).tolist()
    y_array = np.arange(ymin, ymax + del_d, del_d).tolist()

    n_array_x = len(x_array)
    n_array_y = len(y_array)
    n_array_z = len(z_array)

    n_grid = n_array_x*n_array_y*n_array_z

    return(n_grid, x_array, y_array, z_array,del_d)


def binit(P, n_grid, x_array, y_array, z_array, del_d):

    npar = np.shape(P)[0]
    # N_x = len(x_array)
    N_y = len(y_array)
    N_z = len(z_array)
    N_yz = N_y*N_z

    min_x = x_array[0]
    min_y = y_array[0]
    min_z = z_array[0]

    index_array = [0 for particle in P]

    bin_array = [[] for bin in range(n_grid)]

    for idx in range(npar):
        x = P[idx, 2]
        y = P[idx, 3]
        z = P[idx, 4]

        xn = int((x-min_x)/del_d)
        yn = int((y-min_y)/del_d)
        zn = int((z-min_z)/del_d)

        n_ar = xn*N_yz + yn*N_z + zn

        bin_array[n_ar].append(idx)
        index_array[idx] = n_ar

    return(bin_array, index_array)

def find_neibs(par_no, bin_array, index_array, n_grid, x_array, y_array, z_array):

    bin_number = index_array[par_no]

    N_y = len(y_array)
    N_z = len(z_array)
    Nyz = N_y*N_z

    nb_0 = bin_number
    nb_1 = bin_number - 1
    nb_2 = bin_number + 1
    nb_3 = bin_number - N_z
    nb_4 = bin_number - N_z + 1
    nb_5 = bin_number - N_z - 1
    nb_6 = bin_number + N_z
    nb_7 = bin_number + N_z + 1
    nb_8 = bin_number + N_z - 1
    nb_9 = bin_number + (Nyz)
    nb_10 = bin_number + (Nyz) + 1
    nb_11 = bin_number + (Nyz) - 1
    nb_12 = bin_number - (Nyz)
    nb_13 = bin_number - (Nyz) + 1
    nb_14 = bin_number - (Nyz) - 1
    nb_15 = bin_number + (Nyz) - N_z
    nb_16 = bin_number + (Nyz)- N_z - 1
    nb_17 = bin_number + (Nyz) - N_z + 1
    nb_18 = bin_number + (Nyz) + N_z
    nb_19 = bin_number + (Nyz) + N_z - 1
    nb_20 = bin_number + (Nyz) + N_z + 1
    nb_21 = bin_number - (Nyz) - N_z
    nb_22 = bin_number - (Nyz)- N_z - 1
    nb_23 = bin_number - (Nyz) - N_z + 1
    nb_24 = bin_number - (Nyz) + N_z
    nb_25 = bin_number - (Nyz) + N_z - 1
    nb_26 = bin_number - (Nyz) + N_z + 1

    nb_array = [nb_0, nb_1, nb_2, nb_3, nb_4, nb_5, nb_6, nb_7, nb_8, nb_9, nb_10, nb_11, nb_12, nb_13, nb_14, nb_15, nb_16, nb_17, nb_18, nb_19, nb_20, nb_21, nb_22, nb_23, nb_24, nb_25, nb_26]

    for bin in nb_array:
        if bin < 0 or bin >= n_grid:
            nb_array.remove(bin)
    neibs_list = []
    for bin in nb_array:
        particles_in_bin = bin_array[bin]
        for particle in particles_in_bin:
            if particle > par_no:
                neibs_list.append(particle)
    return(neibs_list)

def neibs_array(P, bin_array, index_array, n_grid, x_array, y_array, z_array):

    N = np.shape(P)[0]
    neibs = [find_neibs(idx, bin_array, index_array, n_grid, x_array, y_array, z_array) for idx in range(N)]
    return(neibs)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def ini_ite(P):
    N = np.shape(P)[0]
    for idx in range(N):
        P[idx,5:11] = [0,0,0,0,0,0]
    return(P)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def timestep_update_constant_accleration(P,N_par,t,a_g):#Z corresponds to the particle array and t corresponds to time step duration

    tt = 0.5*(t**2)
    for i in range(N_par):
        par = P[i].tolist()

        mass_par = P[i][17]

        Ax = par[5]/(mass_par) + a_g[0]
        Ay = par[6]/(mass_par) + a_g[1]
        Az = par[7]/(mass_par) + a_g[2]


        par[2] += t*par[11] + (tt*Ax) #-\
        par[3] += t*par[12] + (tt*Ay) #- |-->update position
        par[4] += t*par[13] + (tt*Az) #-/

        par[11] += Ax*t
        par[12] += Ay*t
        par[13] += Az*t
        MoI = (2/5)*mass_par*(par[1]**2) #compute moment of inertia of the sphere

        Ang_ax = par[8]/MoI #-\
        Ang_ay = par[9]/MoI #- |-->compute angular accleration
        Ang_az = par[10]/MoI #-/

        par[14] += Ang_ax*t
        par[15] += Ang_ay*t
        par[16] += Ang_az*t
        P[i] = par
    return(P)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#total kinetic energy
def total_KE(P):
    
    N = np.shape(P)[0]
    total_KE = 0
    for idx in range(N):
        mass_par = P[idx][17]
        coef = 0.5*mass_par
        total_KE += coef*(P[idx][11]**2+P[idx][12]**2+P[idx][13]**2)
    return(total_KE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def place_top_plate(P):

   N = np.shape(P)[0]

   top_plt = 0

   i = 0

   while i < N:

      if P[i][1]+P[i][3] > top_plt:

         top_plt = P[i][1]+P[i][3]

      i = i + 1

   return(top_plt)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Compute parameters used for hertzian contact

def yesn(y1, y2, e1, e2):
    y = y1*y2
    e1s = e1**2
    e2s = e2**2
    return(2*y/(y2-y2*e1s+y1-y1*e2s))

def gest(y1, y2, e1, e2):
    y = y1*y2
    e1s = e1**2
    e2s = e2**2
    return(4*y/(2*y2+y2*e1-y2*e1s+2*y1+y1*e2-y1*e2s))

def knc(y1, y2, e1, e2):
    return(0.66667*yesn(y1, y2, e1, e2))

def ktc(y1, y2, e1, e2):
    return(gest(y1, y2, e1, e2))


def hertz_contact_parameters_interparticle(E):
    npar = np.shape(E)[0]
    hc_paras = []
    yesnmat = np.zeros([npar, npar])
    for idx in range(npar):
        for jdx in range(npar):
            yesnmat[idx, jdx] = yesn(E[idx, 0], E[jdx, 0], E[idx, 2], E[jdx, 2])
    hc_paras.append(yesnmat)

    gestmat = np.zeros([npar, npar])
    for idx in range(npar):
        for jdx in range(npar):
            gestmat[idx, jdx] = gest(E[idx, 0], E[jdx, 0], E[idx, 2], E[jdx, 2])
    hc_paras.append(gestmat)

    kncmat = np.zeros([npar, npar])
    for idx in range(npar):
        for jdx in range(npar):
            kncmat[idx, jdx] = knc(E[idx, 0], E[jdx, 0], E[idx, 2], E[jdx, 2])
    hc_paras.append(kncmat)

    ktcmat = np.zeros([npar, npar])
    for idx in range(npar):
        for jdx in range(npar):
            ktcmat[idx, jdx] = ktc(E[idx, 0], E[jdx, 0], E[idx, 2], E[jdx, 2])
    hc_paras.append(ktcmat)

    return(hc_paras)

def hertz_contact_parameters_wall(E, Youngs_modulus_wall, Shear_modulus_wall, Poissons_ratio_wall):
    npar = np.shape(E)[0]
    hc_paras = []
    hc_paras = [[yesn(E[idx, 0], Youngs_modulus_wall, E[idx, 2], Poissons_ratio_wall),
    gest(E[idx, 0], Youngs_modulus_wall, E[idx, 2], Poissons_ratio_wall),
    knc(E[idx, 0], Youngs_modulus_wall, E[idx, 2], Poissons_ratio_wall),
    ktc(E[idx, 0], Youngs_modulus_wall, E[idx, 2], Poissons_ratio_wall)] for idx in range(npar)]
    return(hc_paras)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def numsign(n):
    if n > 0:
        return(1)
    if n < 0:
        return(-1)
    if n == 0:
        return(0)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#compute forces on partcles due to collision with each other
#hertzian contact

def t_overlap_ip_init(P, neibs):
    N = np.shape(P)[0]
    t_overlap_ip_list = [{j:[0, 0] for j in neibs[i]} for i in range(N)]
    return(t_overlap_ip_list)

def t_overlap_ip_init_pbc(P):
    N = np.shape(P)[0]
    t_overlap_ip_list = [{j:[0, 0] for j in range(N)} for i in range(N)]
    return(t_overlap_ip_list)

def t_overlap_update(P, neibs, t_overlap_hc):
    N = np.shape(P)[0]
    t_overlap_ip_list = [{j:[0, 0] for j in neibs[i]} for i in range(N)]
    for i in range(N):
        common_keys = t_overlap_hc[i].keys() & t_overlap_ip_list[i].keys()
        for j in common_keys:
            t_overlap_ip_list[i][j] = t_overlap_hc[i][j]
    return(t_overlap_ip_list)
        
        
def Hertz_contact_interparticle(P, E, hc_parameters_interparticle, timestep, t_overlap_ip_list, neibs):

    N = np.shape(P)[0]
  
    for idx in range(N):
        t_overlap = t_overlap_ip_list[idx]
        neibs_idx = neibs[idx]
        par_idx = P[idx].tolist()
        [pi, ri, xi, yi, zi, fxi, fyi, fzi, txi, tyi, tzi, ui, vi, wi, oxi, oyi, ozi, mass_i] = par_idx
        pi = int(pi)
        b_coef = E[pi, 7]
        mu = E[pi, 4]
        mu_r = E[pi, 5]
        fxg = 0.0
        fyg = 0.0
        fzg = 0.0
        txg = 0.0
        tyg = 0.0
        tzg = 0.0

        for jdx in neibs_idx:
            par_jdx = P[jdx].tolist()
            [pj, rj, xj, yj, zj, fxj, fyj, fzj, txj, tyj, tzj, uj, vj, wj, oxj, oyj, ozj, mass_j] = par_jdx
            rad_sum = ri + rj
            if xi-xj > rad_sum or yi-yj >rad_sum or zi-zj>rad_sum:
                t_overlap[jdx] = [0, 0]
            else:
                cd = ((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)
                if cd >= rad_sum**2:
                    t_overlap[jdx] = [0, 0]
                else:
                    pj = int(pj)
                    cd = math.sqrt(cd)
                    overlap = ri+rj-cd
                    contact_angle_zx = math.atan2((xj-xi), (zj-zi))
                    ang_ver = math.asin((yj-yi)/cd)
                    RE = ri*rj/(ri+rj)
                    ME = mass_i*mass_j/(mass_i+mass_j)
                    yesn = hc_parameters_interparticle[0][pi, pj]
                    gest = hc_parameters_interparticle[1][pi, pj]
                    knc = hc_parameters_interparticle[2][pi, pj]
                    ktc = hc_parameters_interparticle[3][pi, pj]
                    [tox, toy] = t_overlap[jdx]
                    
                    olph = overlap/2
                    mali = ri - olph
                    malj = rj - olph

                    vel_g_x = ui-uj
                    vel_g_y = vi-vj
                    vel_g_z = wi-wj

                    reo = math.sqrt(RE*overlap)

                    SN, ST, kn, kt = yesn*(reo), gest*(reo), knc*(reo), ktc*(reo)
                    bn, bt = b_coef*math.sqrt(SN*ME), b_coef*math.sqrt(ST*ME)

                    S1, C1 = math.sin(contact_angle_zx), math.cos(contact_angle_zx)
                    S2, C2 = math.sin(ang_ver), math.cos(ang_ver)

                    coef1 = S1*S2
                    coef2 = S2*C1
                    coef3 = S1*C2
                    coef4 = C1*C2

                    U_rel_loc = C1*vel_g_x - S1*vel_g_z
                    V_rel_loc = - coef1*vel_g_x + C2*vel_g_y - coef2*vel_g_z
                    W_rel_loc = coef3*vel_g_x + S2*vel_g_y + coef4*vel_g_z

                    oxil = C1*oxi - S1*ozi
                    oyil = - coef1*oxi + C2*oyi - coef2*ozi
                    ozil = coef3*oxi + S2*oyi + coef4*ozi

                    oxjl = C1*oxj - S1*ozj
                    oyjl = - coef1*oxj + C2*oyj - coef2*ozj
                    ozjl = coef3*oxj + S2*oyj + coef4*ozj

                    oxl_rel = oxil-oxjl
                    oyl_rel = oyil-oyjl
                    ozl_rel = ozil-ozjl

                    Fz_local = (-kn*overlap + bn*W_rel_loc)

                    relative_sliding_velocity_x = U_rel_loc + oyil*ri + oyjl*rj
                    relative_sliding_velocity_y = V_rel_loc - (oxil*ri + oxjl*rj)

                    muf = np.abs(Fz_local)*mu

                    Fx_local = -kt*tox + bt*relative_sliding_velocity_x
                    if abs(Fx_local) > abs(muf):
                        Fx_local = -muf*np.sign(relative_sliding_velocity_x)

                    Fy_local = -kt*toy + bt*relative_sliding_velocity_y
                    if abs(Fy_local) > abs(muf):
                        Fy_local = -muf*np.sign(relative_sliding_velocity_y)

                    t_overlap[jdx] = [tox+relative_sliding_velocity_x*timestep, toy+relative_sliding_velocity_y*timestep]

                    T_x_loc = -Fy_local*mali
                    T_y_loc = Fx_local*mali

                    T_x_locj = -Fy_local*malj
                    T_y_locj = Fx_local*malj

                    trcoef = mu_r*RE*Fz_local

                    tmurx = trcoef*numsign(oxl_rel)
                    tmury = trcoef*numsign(oyl_rel)
                    tmurz = trcoef*numsign(ozl_rel)

                    fxgj = - C1*Fx_local + coef1*Fy_local - coef3*Fz_local
                    fygj = - C2*Fy_local - S2*Fz_local
                    fzgj = + S1*Fx_local + coef2*Fy_local - coef4*Fz_local

                    txgi =  C1*T_x_loc - coef1*T_y_loc
                    tygi = C2*T_y_loc
                    tzgi = -S1*T_x_loc - coef2*T_y_loc

                    txgj =  C1*T_x_locj - coef1*T_y_locj
                    tygj = C2*T_y_locj
                    tzgj = -S1*T_x_locj - coef2*T_y_locj

                    fxg -= fxgj
                    fyg -= fygj
                    fzg -= fzgj

                    txg += txgi+tmurx
                    tyg += tygi+tmury
                    tzg += tzgi+tmurz

                    P[jdx, 5:11] += [fxgj, fygj, fzgj, txgj-tmurx, tygj-tmury, tzgj-tmurz]

        t_overlap_ip_list[idx] = t_overlap
        P[idx, 5:11] += [fxg, fyg, fzg, txg, tyg, tzg]

    return(P, t_overlap_ip_list)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def t_overlap_wall(P):
    N = np.shape(P)[0]
    t_overlap_wall_array = [[0, 0] for idx in range(N)]
    return(t_overlap_wall_array)

def t_overlap_wall_ap(tovlpw,tovlp_par):

    tovlpw.append(tovlp_par)
    return(tovlpw)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#hertz contact for cylindrical wall
def Hertz_contact_cylinder_wall(P, E, Cyl_radius, hc_parameter_cylinder, timestep, t_overlap_cylinder):
    normal_force = 0
    N = np.shape(P)[0]

    for idx in range(N):
        par_idx = P[idx]
        [pi, ri, xi, yi, zi, fxi, fyi, fzi, txi, tyi, tzi, ui, vi, wi, oxi, oyi, ozi, mass_i] = par_idx
        yaxis_dist = xi**2 + zi**2
        if yaxis_dist >= (Cyl_radius-ri)**2:
            pi = int(pi)
            [toh, tov] = t_overlap_cylinder[idx]
            yaxis_dist = math.sqrt(yaxis_dist)
            normal_overlap_magnitude = yaxis_dist + ri - Cyl_radius
            [yesn, gest, knc, ktc] = hc_parameter_cylinder[pi]
            mu = E[pi, 4]
            mu_r = E[pi, 5]
            reo = math.sqrt(ri*normal_overlap_magnitude)
            SN, ST, kn, kt = yesn*reo, gest*reo, knc*reo, ktc*reo
            bet_coef = E[pi, 7]
            bn, bt = bet_coef*math.sqrt(SN*mass_i), bet_coef*math.sqrt(ST*mass_i)
            y_axis_angle = math.atan2(xi, zi)
            C = np.cos(y_axis_angle)
            S = np.sin(y_axis_angle)
            Vz_local = wi*C + ui*S
            Vx_local = -wi*S + ui*C
            Wx_loc = - ozi*S + oxi*C
            horizontal_sliding_velocity = Vx_local + oyi*ri #Vx_local + Angular velocity along y axis
            vertical_sliding_velocity = vi - Wx_loc*ri
            F_z_local = -kn*normal_overlap_magnitude+bn*Vz_local
            normal_force += kn*normal_overlap_magnitude
            F_x_local = (-toh*kt + horizontal_sliding_velocity*bt) #F_tangential_horizontal(friction)
            if np.abs(F_x_local) >= np.abs(F_z_local)*mu:
               F_x_local = np.abs(F_z_local)*mu*np.sign(horizontal_sliding_velocity)
            F_vertical = -tov*kt + vertical_sliding_velocity*bt#F_y_global
            if np.abs(F_vertical) >= np.abs(F_z_local)*mu:
               F_vertical = np.abs(F_z_local)*mu*np.sign(vertical_sliding_velocity)
            fy = F_vertical - mu_r*ri*F_z_local*np.sign(oyi)#F_y_global#Torque due to kinetic friction and rolling resistance
            ty = F_x_local*ri #torque along y-axis
            T_x_local = -F_vertical*mu-mu_r*ri*F_z_local*np.sign(Wx_loc)#Torque due to kinetic friction and rolling resistance
            fz = F_z_local*np.cos(-y_axis_angle) - F_x_local*np.sin(-y_axis_angle)#F_Z_global
            fx = - F_z_local*np.sin(- y_axis_angle) - F_x_local*np.cos(- y_axis_angle)#F_X_global
            tz = - np.sin(- y_axis_angle)*T_x_local #torque along z-axis
            tx = - np.cos(- y_axis_angle)*T_x_local
            t_overlap_cylinder[idx] = [toh+horizontal_sliding_velocity*timestep, tov+vertical_sliding_velocity*timestep]
            P[idx, 5:11] += [fx, fy, fz, tx, ty, tz]
        else:
            t_overlap_cylinder[idx] = [0, 0]
    return(P, t_overlap_cylinder, normal_force)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#hertz contact bottom wall
def Hertz_contact_bottom_wall(P, E, base, hc_parameter_bottom_wall, timestep, t_overlap_bottom_wall):
    normal_force = 0
    N = np.shape(P)[0]

    for idx in range(N):
        par_idx = P[idx]
        [pi, ri, xi, yi, zi, fxi, fyi, fzi, txi, tyi, tzi, ui, vi, wi, oxi, oyi, ozi, mass_i] = par_idx
        if yi < base+ri:
            pi = int(pi)
            [tox,toz] = t_overlap_bottom_wall[idx]
            normal_overlap_magnitude = ri-(yi-base)
            [yesn, gest, knc, ktc] = hc_parameter_bottom_wall[pi]
            mu = E[pi, 4]
            mu_r = E[pi, 5]
            reo = math.sqrt(ri*normal_overlap_magnitude)
            SN, ST, kn, kt = yesn*reo, gest*reo, knc*reo, ktc*reo
            bet_coef = E[pi, 7]
            bn, bt = bet_coef*math.sqrt(SN*mass_i), bet_coef*math.sqrt(ST*mass_i)
            F_vertical =  kn*normal_overlap_magnitude + bn*vi
            normal_force += -kn*normal_overlap_magnitude
            sliding_velocity_z = wi - oxi*ri
            sliding_velocity_x = ui + ozi*ri
            Fx = -(kt*tox - bt*sliding_velocity_x)
            if np.abs(Fx) >= np.abs(F_vertical)*mu:
               Fx = -np.abs(F_vertical)*mu*np.sign(sliding_velocity_x)
            Fz = -(kt*toz - bt*sliding_velocity_z)
            if np.abs(Fz) >= np.abs(F_vertical)*mu:
               Fz = -np.abs(F_vertical)*mu*np.sign(sliding_velocity_z)
            tz =  (Fx/ri) - mu_r*ri*np.abs(F_vertical)*np.sign(ozi)
            tx = -(Fz/ri) - mu_r*ri*np.abs(F_vertical)*np.sign(oxi)
            t_overlap_bottom_wall[idx] = [tox+sliding_velocity_x*timestep,toz+sliding_velocity_z*timestep]
            P[idx,5:11] += [Fx,F_vertical,Fz,tx,0,tz]
        else:
            t_overlap_bottom_wall[idx] = [0,0]
    return(P,t_overlap_bottom_wall,normal_force)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#hertz contact arbitrary wall defined by unit normal vector and a point on the
#plane
def Hertz_contact_arb_wall(P, E, norm_vect, ref_point ,hc_parameter_arb_wall, timestep, t_overlap_arb_wall):
    normal_force = 0
    N = np.shape(P)[0]
    
    a = norm_vect[0]
    b = norm_vect[1]
    c = norm_vect[2]
    
    xr = ref_point[0]
    yr = ref_point[1]
    zr = ref_point[2]
    
    d = -(a*xr + b*yr + c*zr)
    
    the1 = -math.atan2(a, b)
    the2 = math.atan2(c, math.sqrt((a**2)+(b**2)))
    
    C1 = math.cos(the1)
    C2 = math.cos(the2)
    S1 = math.sin(the1)
    S2 = math.sin(the2)
    
    coef1 = -S1*C2
    coef2 = S1*S2
    coef3 = C1*C2
    coef4 = -C1*S2
    for idx in range(N):
        par_idx = P[idx]
        [pi, ri, xi, yi, zi, fxgi, fygi, fzgi, txgi, tygi, tzgi, ugi, vgi, wgi, oxgi, oygi, ozgi, mass_i] = par_idx
        cdsqr = np.abs(((a*xi+b*yi+c*zi+d)**2)/((a**2)+(b**2)+(c**2)))
        if cdsqr < ri**2:
            normal_overlap_magnitude = ri - math.sqrt(cdsqr)
            ui = C1*ugi + S1*vgi
            vi = coef1*ugi + coef3*vgi  + S2*wgi
            wi = coef2*ugi + coef4*vgi + C2*wgi
            oxi = C1*oxgi + S1*oygi
            ozi = coef2*oxgi + coef4*oygi + C2*ozgi
            pi = int(pi)
            [tox,toz] = t_overlap_arb_wall[idx]
            [yesn, gest, knc, ktc] = hc_parameter_arb_wall[pi]
            mu = E[pi, 4]
            mu_r = E[pi, 5]
            reo = math.sqrt(ri*normal_overlap_magnitude)
            SN, ST, kn, kt = yesn*reo, gest*reo, knc*reo, ktc*reo
            bet_coef = E[pi, 7]
            bn, bt = bet_coef*math.sqrt(SN*mass_i), bet_coef*math.sqrt(ST*mass_i)
            
            Fly =  kn*normal_overlap_magnitude + bn*vi
            
            normal_force += -kn*normal_overlap_magnitude
            
            sliding_velocity_z = wi - oxi*ri
            sliding_velocity_x = ui + ozi*ri
            
            Flx = -(kt*tox - bt*sliding_velocity_x)
            if np.abs(Flx) >= np.abs(Fly)*mu:
               Flx = -np.abs(Fly)*mu*numsign(sliding_velocity_x)
            Flz = -(kt*toz - bt*sliding_velocity_z)
            if np.abs(Flz) >= np.abs(Fly)*mu:
               Flz = -np.abs(Fly)*mu*numsign(sliding_velocity_z)
               
            tlz =  (Flx/ri) - mu_r*ri*np.abs(Fly)*np.sign(ozi)
            tlx = -(Flz/ri) - mu_r*ri*np.abs(Fly)*np.sign(oxi)
            
            t_overlap_arb_wall[idx] = [tox+sliding_velocity_x*timestep,toz+sliding_velocity_z*timestep]
            
            Fx = C1*Flx + coef1*Fly + coef2*Flz
            Fy = S1*Flx + coef3*Fly + coef4*Flz
            Fz = S2*Fly + C2*Flz
            tx = C1*tlx + coef2*tlz
            ty = S1*tlx + coef4*tlz
            tz = C2*tlz
        
            P[idx,5:11] += [Fx,Fy,Fz,tx,ty,tz]
        else:
            t_overlap_arb_wall[idx] = [0,0]
    return(P,t_overlap_arb_wall,normal_force)
        
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#horizontal 4-wall periodic boundary condition(cubic)

def pbc_find_neibs(bin_number, n_grid, bin_array, x_array, y_array, z_array):
    N_y = len(y_array)
    N_z = len(z_array)
    Nyz = N_y*N_z

    nb_0 = bin_number
    nb_1 = bin_number - 1
    nb_2 = bin_number + 1
    nb_3 = bin_number - N_z
    nb_4 = bin_number - N_z + 1
    nb_5 = bin_number - N_z - 1
    nb_6 = bin_number + N_z
    nb_7 = bin_number + N_z + 1
    nb_8 = bin_number + N_z - 1
    nb_9 = bin_number + (Nyz)
    nb_10 = bin_number + (Nyz) + 1
    nb_11 = bin_number + (Nyz) - 1
    nb_12 = bin_number - (Nyz)
    nb_13 = bin_number - (Nyz) + 1
    nb_14 = bin_number - (Nyz) - 1
    nb_15 = bin_number + (Nyz) - N_z
    nb_16 = bin_number + (Nyz)- N_z - 1
    nb_17 = bin_number + (Nyz) - N_z + 1
    nb_18 = bin_number + (Nyz) + N_z
    nb_19 = bin_number + (Nyz) + N_z - 1
    nb_20 = bin_number + (Nyz) + N_z + 1
    nb_21 = bin_number - (Nyz) - N_z
    nb_22 = bin_number - (Nyz)- N_z - 1
    nb_23 = bin_number - (Nyz) - N_z + 1
    nb_24 = bin_number - (Nyz) + N_z
    nb_25 = bin_number - (Nyz) + N_z - 1
    nb_26 = bin_number - (Nyz) + N_z + 1

    nb_array = [nb_0, nb_1, nb_2, nb_3, nb_4, nb_5, nb_6, nb_7, nb_8, nb_9, nb_10, nb_11, nb_12, nb_13, nb_14, nb_15, nb_16, nb_17, nb_18, nb_19, nb_20, nb_21, nb_22, nb_23, nb_24, nb_25, nb_26]

    for bin in nb_array:
        if bin < 0 or bin >= n_grid:
            nb_array.remove(bin)
    neibs_list = []
    for bin in nb_array:
        particles_in_bin = bin_array[bin]
        for particle in particles_in_bin:
            neibs_list.append(particle)
    return(neibs_list)

def pbc_hc(P, E, hc_parameters_interparticle, timestep, t_overlap_ip_list, imaginary_particles, bin_array, index_array, n_grid, x_array, y_array, z_array):
    N = np.shape(imaginary_particles)[0]
  
    for idx in range(N):
        original_par = imaginary_particles[idx][0]
        
        bin_no = imaginary_particles[idx][4]
        
        neibs_idx = pbc_find_neibs(bin_no, n_grid, bin_array, x_array, y_array, z_array)
        # for par in neibs_idx:
        #     if par not in t_overlap_ip_list[original_par].keys():
        #         t_overlap_ip_list[original_par][par] = [0,0]
        t_overlap = t_overlap_ip_list[original_par]    
        par_idx = P[original_par].tolist()
        [pi, ri, xi, yi, zi, fxi, fyi, fzi, txi, tyi, tzi, ui, vi, wi, oxi, oyi, ozi, mass_i] = par_idx
        xi =  imaginary_particles[idx][1]
        yi =  imaginary_particles[idx][2]
        zi =  imaginary_particles[idx][3]
        pi = int(pi)
        b_coef = E[pi, 7]
        mu = E[pi, 4]
        mu_r = E[pi, 5]
        fxg = 0.0
        fyg = 0.0
        fzg = 0.0
        txg = 0.0
        tyg = 0.0
        tzg = 0.0

        for jdx in neibs_idx:
            par_jdx = P[jdx].tolist()
            [pj, rj, xj, yj, zj, fxj, fyj, fzj, txj, tyj, tzj, uj, vj, wj, oxj, oyj, ozj, mass_j] = par_jdx
            rad_sum = ri + rj
            if xi-xj > rad_sum or yi-yj >rad_sum or zi-zj>rad_sum:
                t_overlap[jdx] = [0, 0]
            else:
                cd = ((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)
                if cd >= rad_sum**2:
                    t_overlap[jdx] = [0, 0]
                else:
                    pj = int(pj)
                    cd = math.sqrt(cd)
                    overlap = ri+rj-cd
                    contact_angle_zx = math.atan2((xj-xi), (zj-zi))
                    ang_ver = math.asin((yj-yi)/cd)
                    RE = ri*rj/(ri+rj)
                    ME = mass_i*mass_j/(mass_i+mass_j)
                    yesn = hc_parameters_interparticle[0][pi, pj]
                    gest = hc_parameters_interparticle[1][pi, pj]
                    knc = hc_parameters_interparticle[2][pi, pj]
                    ktc = hc_parameters_interparticle[3][pi, pj]
                    [tox, toy] = t_overlap[jdx]
                    
                    olph = overlap/2
                    mali = ri - olph
                    malj = rj - olph

                    vel_g_x = ui-uj
                    vel_g_y = vi-vj
                    vel_g_z = wi-wj

                    reo = math.sqrt(RE*overlap)

                    SN, ST, kn, kt = yesn*(reo), gest*(reo), knc*(reo), ktc*(reo)
                    bn, bt = b_coef*math.sqrt(SN*ME), b_coef*math.sqrt(ST*ME)

                    S1, C1 = math.sin(contact_angle_zx), math.cos(contact_angle_zx)
                    S2, C2 = math.sin(ang_ver), math.cos(ang_ver)

                    coef1 = S1*S2
                    coef2 = S2*C1
                    coef3 = S1*C2
                    coef4 = C1*C2

                    U_rel_loc = C1*vel_g_x - S1*vel_g_z
                    V_rel_loc = - coef1*vel_g_x + C2*vel_g_y - coef2*vel_g_z
                    W_rel_loc = coef3*vel_g_x + S2*vel_g_y + coef4*vel_g_z

                    oxil = C1*oxi - S1*ozi
                    oyil = - coef1*oxi + C2*oyi - coef2*ozi
                    ozil = coef3*oxi + S2*oyi + coef4*ozi

                    oxjl = C1*oxj - S1*ozj
                    oyjl = - coef1*oxj + C2*oyj - coef2*ozj
                    ozjl = coef3*oxj + S2*oyj + coef4*ozj

                    oxl_rel = oxil-oxjl
                    oyl_rel = oyil-oyjl
                    ozl_rel = ozil-ozjl

                    Fz_local = (-kn*overlap + bn*W_rel_loc)

                    relative_sliding_velocity_x = U_rel_loc + oyil*ri + oyjl*rj
                    relative_sliding_velocity_y = V_rel_loc - (oxil*ri + oxjl*rj)

                    muf = np.abs(Fz_local)*mu

                    Fx_local = -kt*tox + bt*relative_sliding_velocity_x
                    if abs(Fx_local) > abs(muf):
                        Fx_local = -muf*np.sign(relative_sliding_velocity_x)

                    Fy_local = -kt*toy + bt*relative_sliding_velocity_y
                    if abs(Fy_local) > abs(muf):
                        Fy_local = -muf*np.sign(relative_sliding_velocity_y)

                    t_overlap[jdx] = [tox+relative_sliding_velocity_x*timestep, toy+relative_sliding_velocity_y*timestep]

                    T_x_loc = -Fy_local*mali
                    T_y_loc = Fx_local*mali

                    T_x_locj = -Fy_local*malj
                    T_y_locj = Fx_local*malj

                    trcoef = mu_r*RE*Fz_local

                    tmurx = trcoef*numsign(oxl_rel)
                    tmury = trcoef*numsign(oyl_rel)
                    tmurz = trcoef*numsign(ozl_rel)

                    fxgj = - C1*Fx_local + coef1*Fy_local - coef3*Fz_local
                    fygj = - C2*Fy_local - S2*Fz_local
                    fzgj = + S1*Fx_local + coef2*Fy_local - coef4*Fz_local

                    txgi =  C1*T_x_loc - coef1*T_y_loc
                    tygi = C2*T_y_loc
                    tzgi = -S1*T_x_loc - coef2*T_y_loc

                    txgj =  C1*T_x_locj - coef1*T_y_locj
                    tygj = C2*T_y_locj
                    tzgj = -S1*T_x_locj - coef2*T_y_locj

                    fxg -= fxgj
                    fyg -= fygj
                    fzg -= fzgj

                    txg += txgi+tmurx
                    tyg += tygi+tmury
                    tzg += tzgi+tmurz

                    P[jdx, 5:11] += [fxgj, fygj, fzgj, txgj-tmurx, tygj-tmury, tzgj-tmurz]

        t_overlap_ip_list[original_par] = t_overlap
        P[original_par, 5:11] += [fxg, fyg, fzg, txg, tyg, tzg]

    return(P, t_overlap_ip_list)
    



def pbc_cube_rearange(xmin,xmax,zmin,zmax,boundary_thickness,P,E,t_overlap_hc,timestep,x_array,y_array,z_array,bin_array,index_array,neibs,max_par_size,bin_size,hc_parameters_interparticle):
    
    del_d = 2*max_par_size*bin_size
    xmin_in = xmin + boundary_thickness
    xmin_out = xmin - boundary_thickness
    xmax_in = xmax - boundary_thickness
    xmax_out = xmax + boundary_thickness
    zmin_in = zmin + boundary_thickness
    zmin_out = zmin - boundary_thickness
    zmax_in = zmax - boundary_thickness
    zmax_out = zmax + boundary_thickness
    xwidth = xmax-xmin
    zwidth = zmax-zmin
    bin_width_z = math.floor(zwidth/del_d)
    N_x = len(x_array)
    N_y = len(y_array)
    N_z = len(z_array)
    min_x = x_array[0]
    min_y = y_array[0]
    min_z = z_array[0]
    N_yz = N_y*N_z
    bin_width_x = math.floor(xwidth/del_d)*N_yz
    n_grid = N_x*N_y*N_z
    
    N = np.shape(P)[0]
    imaginary_particles = []
    
    for idx in range(N):
        [pi, ri, xi, yi, zi, fxi, fyi, fzi, txi, tyi, tzi, ui, vi, wi, oxi, oyi, ozi, mass_i] = P[idx]
        
        if zi <= zmin_out and xmin < xi < xmax:
            P[idx][4] += zwidth
        
        elif zi >= zmax_out and xmin < xi < xmax:
            P[idx][4] -= zwidth
            
        elif xi <= xmin_out and zmin < zi < zmax:
            P[idx][2] += xwidth
        
        elif xi >= xmax_out and zmin < zi < zmax:
            P[idx][2] -= xwidth
            
        elif zi <= zmin_out and xi <= xmin_out:
            P[idx][2] += xwidth
            P[idx][4] += zwidth
            
        elif zi <= zmin_out and xi >= xmax_out:
            P[idx][2] -= xwidth
            P[idx][4] += zwidth
            
        elif zi >= zmax_out and xi >= xmax_out:
            P[idx][2] -= xwidth
            P[idx][4] -= zwidth
            
        elif zi >= zmax_out and xi <= xmin_out:
            P[idx][2] += xwidth
            P[idx][4] -= zwidth
            
        
        
        
    return(P)
        
        
def pbc_cube_imaginary(xmin,xmax,zmin,zmax,boundary_thickness,P,E,t_overlap_hc,timestep,x_array,y_array,z_array,bin_array,index_array,neibs,max_par_size,bin_size,hc_parameters_interparticle):
    
    del_d = 2*max_par_size*bin_size
    xmin_in = xmin + boundary_thickness
    xmin_out = xmin - boundary_thickness
    xmax_in = xmax - boundary_thickness
    xmax_out = xmax + boundary_thickness
    zmin_in = zmin + boundary_thickness
    zmin_out = zmin - boundary_thickness
    zmax_in = zmax - boundary_thickness
    zmax_out = zmax + boundary_thickness
    xwidth = xmax-xmin
    zwidth = zmax-zmin
    bin_width_z = math.floor(zwidth/del_d)
    N_x = len(x_array)
    N_y = len(y_array)
    N_z = len(z_array)
    min_x = x_array[0]
    min_y = y_array[0]
    min_z = z_array[0]
    N_yz = N_y*N_z
    bin_width_x = math.floor(xwidth/del_d)*N_yz
    n_grid = N_x*N_y*N_z
    
    N = np.shape(P)[0]
    imaginary_particles = []
        
    for idx in range(N):
        [pi, ri, xi, yi, zi, fxi, fyi, fzi, txi, tyi, tzi, ui, vi, wi, oxi, oyi, ozi, mass_i] = P[idx]
        
        if zi <= zmin_in and xmin_in < xi < xmax_in:
            xim = xi
            zim = zi + zwidth
            yim = yi
            xn = int((xim-min_x)/del_d)
            yn = int((yim-min_y)/del_d)
            zn = int((zim-min_z)/del_d)

            bin_im = xn*N_yz + yn*N_z + zn
            imaginary_particles.append([idx,xim,yim,zim,bin_im])
            
        
        elif zi >= zmax_in and xmin_in < xi < xmax_in:
            xim = xi
            zim = zi - zwidth
            yim = yi
            xn = int((xim-min_x)/del_d)
            yn = int((yim-min_y)/del_d)
            zn = int((zim-min_z)/del_d)

            bin_im = xn*N_yz + yn*N_z + zn
            imaginary_particles.append([idx,xim,yim,zim,bin_im])
        
        elif xi <= xmin_in and zmin_in < zi < zmax_in:
            xim = xi + xwidth
            zim = zi
            yim = yi
            xn = int((xim-min_x)/del_d)
            yn = int((yim-min_y)/del_d)
            zn = int((zim-min_z)/del_d)

            bin_im = xn*N_yz + yn*N_z + zn
            imaginary_particles.append([idx,xim,yim,zim,bin_im])
            
        
        elif xi >= xmax_in and zmin_in < zi < zmax_in:
            xim = xi - xwidth
            zim = zi 
            yim = yi
            xn = int((xim-min_x)/del_d)
            yn = int((yim-min_y)/del_d)
            zn = int((zim-min_z)/del_d)

            bin_im = xn*N_yz + yn*N_z + zn
            imaginary_particles.append([idx,xim,yim,zim,bin_im])
            
            
        elif zi <= zmin_in and xi <= xmin_in:
            xim = xi + xwidth
            zim = zi 
            yim = yi
            xn = int((xim-min_x)/del_d)
            yn = int((yim-min_y)/del_d)
            zn = int((zim-min_z)/del_d)

            bin_im = xn*N_yz + yn*N_z + zn
            imaginary_particles.append([idx,xim,yim,zim,bin_im])
            
            xim = xi 
            zim = zi + zwidth
            yim = yi
            xn = int((xim-min_x)/del_d)
            yn = int((yim-min_y)/del_d)
            zn = int((zim-min_z)/del_d)

            bin_im = xn*N_yz + yn*N_z + zn
            imaginary_particles.append([idx,xim,yim,zim,bin_im])
            
            
            xim = xi + xwidth
            zim = zi + zwidth
            yim = yi
            xn = int((xim-min_x)/del_d)
            yn = int((yim-min_y)/del_d)
            zn = int((zim-min_z)/del_d)

            bin_im = xn*N_yz + yn*N_z + zn
            imaginary_particles.append([idx,xim,yim,zim,bin_im])
            
        elif zi <= zmin_in and xi >= xmax_in:
            xim = xi - xwidth
            zim = zi 
            yim = yi
            xn = int((xim-min_x)/del_d)
            yn = int((yim-min_y)/del_d)
            zn = int((zim-min_z)/del_d)

            bin_im = xn*N_yz + yn*N_z + zn
            imaginary_particles.append([idx,xim,yim,zim,bin_im])
            
            xim = xi 
            zim = zi + zwidth
            yim = yi
            xn = int((xim-min_x)/del_d)
            yn = int((yim-min_y)/del_d)
            zn = int((zim-min_z)/del_d)

            bin_im = xn*N_yz + yn*N_z + zn
            imaginary_particles.append([idx,xim,yim,zim,bin_im])
            
            
            xim = xi - xwidth
            zim = zi + zwidth
            yim = yi
            xn = int((xim-min_x)/del_d)
            yn = int((yim-min_y)/del_d)
            zn = int((zim-min_z)/del_d)

            bin_im = xn*N_yz + yn*N_z + zn
            imaginary_particles.append([idx,xim,yim,zim,bin_im])
            
        elif zi >= zmax_in and xi >= xmax_in:
            xim = xi - xwidth
            zim = zi 
            yim = yi
            xn = int((xim-min_x)/del_d)
            yn = int((yim-min_y)/del_d)
            zn = int((zim-min_z)/del_d)

            bin_im = xn*N_yz + yn*N_z + zn
            imaginary_particles.append([idx,xim,yim,zim,bin_im])
            
            xim = xi 
            zim = zi - zwidth
            yim = yi
            xn = int((xim-min_x)/del_d)
            yn = int((yim-min_y)/del_d)
            zn = int((zim-min_z)/del_d)

            bin_im = xn*N_yz + yn*N_z + zn
            imaginary_particles.append([idx,xim,yim,zim,bin_im])
            
            
            xim = xi - xwidth
            zim = zi - zwidth
            yim = yi
            xn = int((xim-min_x)/del_d)
            yn = int((yim-min_y)/del_d)
            zn = int((zim-min_z)/del_d)

            bin_im = xn*N_yz + yn*N_z + zn
            imaginary_particles.append([idx,xim,yim,zim,bin_im])
            
        elif zi >= zmax_in and xi <= xmin_in:
            xim = xi + xwidth
            zim = zi 
            yim = yi
            xn = int((xim-min_x)/del_d)
            yn = int((yim-min_y)/del_d)
            zn = int((zim-min_z)/del_d)

            bin_im = xn*N_yz + yn*N_z + zn
            imaginary_particles.append([idx,xim,yim,zim,bin_im])
            
            xim = xi 
            zim = zi - zwidth
            yim = yi
            xn = int((xim-min_x)/del_d)
            yn = int((yim-min_y)/del_d)
            zn = int((zim-min_z)/del_d)

            bin_im = xn*N_yz + yn*N_z + zn
            imaginary_particles.append([idx,xim,yim,zim,bin_im])
            
            
            xim = xi + xwidth
            zim = zi - zwidth
            yim = yi
            xn = int((xim-min_x)/del_d)
            yn = int((yim-min_y)/del_d)
            zn = int((zim-min_z)/del_d)

            bin_im = xn*N_yz + yn*N_z + zn
            imaginary_particles.append([idx,xim,yim,zim,bin_im])
            
        
    P,t_overlap_hc = pbc_hc(P, E, hc_parameters_interparticle, timestep, t_overlap_hc, imaginary_particles, bin_array, index_array, n_grid, x_array, y_array, z_array)
    
    return(P,t_overlap_hc)
    
    
    
            
      
        