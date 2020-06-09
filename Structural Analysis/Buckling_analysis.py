import numpy as np
import matplotlib.pyplot as plt

#Input parameters
stringer_type = 2   #value 0 means hat, 1 means I, 2 means z, 3 means t (that is with the horizontal part of t not attached to skin)
E = 70000.0         #N/mm^2
v = 0.3             #- Poisson ratio
k_ss = 4.0          #value for simply supported plate
k_ssf = 0.43        #value for 3 sides simply supported and one edge free
t_skin = 2.5        #mm
t_stiff = 2.5       #mm
h_stiff = 50        #mm height of the stiffener
w_stiff = 30        #mm (please note that this is half the total width for z stiffener and only includes the top (of hat) width for hat stiffener)
d_stiff = 5         #mm (set this 0 for everything except hat, for hat it is the small part of the hat stiffener being in touch with the skin)
a = 1000            #mm distance between two ribs
b = 120             #mm (stringer spacing)
c = 3.0             #value 4 for bolt (HiLok/Jobolts), 3.5 for spot welds, 3 for universal head rivets, 1.5 for countersunk rivets (tension head) and 1 for countersunk rivets (shear head)
s = 90              #mm (spacing between rivets)

Limit_load = 600            #N/mm
Ultimate_load = 900         #N/mm

#first failure -> Local skin buckling
sigma_crit_skin = (np.pi**2 * k_ss * E) / (12 * (1 - v**2)) * (t_skin/b)**2

#Second failure -> Local stiffener buckling
if stringer_type == 0: #for hat stringer
    sigma_cr_baseflange = (np.pi**2 * k_ssf * E) / (12*(1-v**2)) * (t_stiff/d_stiff)**2
    sigma_cr_web = (np.pi**2 * k_ss * E) / (12*(1-v**2)) * (t_stiff/h_stiff)**2
    sigma_cr_topflange = (np.pi**2 * k_ss * E) / (12*(1-v**2)) * (t_stiff/w_stiff)**2
    
    sigma_cr_stiff = min(sigma_cr_baseflange, sigma_cr_web, sigma_cr_topflange)
    
elif stringer_type == 1: #for I stringer
    sigma_cr_flange = (np.pi**2 * k_ssf * E) / (12*(1-v**2)) * (t_stiff/(w_stiff/2))**2
    sigma_cr_web = (np.pi**2 * k_ss * E) / (12*(1-v**2)) * (t_stiff/h_stiff)**2
    
    sigma_cr_stiff = min(sigma_cr_flange, sigma_cr_web)
    
elif stringer_type == 2: #for Z-stringer
    sigma_cr_flange = (np.pi**2 * k_ssf * E) / (12*(1-v**2)) * (t_stiff/(w_stiff))**2
    sigma_cr_web = (np.pi**2 * k_ss * E) / (12*(1-v**2)) * (t_stiff/h_stiff)**2
    
    sigma_cr_stiff = min(sigma_cr_flange, sigma_cr_web)
    
else: #for T_stringer
    sigma_cr_flange = (np.pi**2 * k_ssf * E) / (12*(1-v**2)) * (t_stiff/(w_stiff/2))**2
    sigma_cr_web = (np.pi**2 * k_ss * E) / (12*(1-v**2)) * (t_stiff/h_stiff)**2
    
    sigma_cr_stiff = min(sigma_cr_flange, sigma_cr_web)
    
    
#Third failure -> Inter rivet buckling
sigma_rivet = (np.pi**2 * c * E) / (12*(1-v**2)) * (t_skin/s)**2

#Fourth failure -> Column stiffness buckling
if stringer_type == 0: #for hat stringer
    A_skin = 30*t_skin*t_skin
    y_skin = 0
    A_str_flange_base = t_stiff * 2 * d_stiff
    y_str_flange_base = 0
    A_str_web = 2*t_stiff*h_stiff
    y_str_web = 0.5*h_stiff
    A_str_flange_top = t_stiff * w_stiff
    y_str_flange_top = h_stiff
    
    y_bar = (A_skin * y_skin + A_str_flange_base*y_str_flange_base + A_str_web*y_str_web + A_str_flange_top*y_str_flange_top) / (A_skin + A_str_flange_base + A_str_web + A_str_flange_top)
    I_yy_skin = A_skin * (y_bar)**2
    I_yy_str_flange_base = A_str_flange_base * (y_bar)**2
    I_yy_str_web = 2*(t_stiff * h_stiff**3)/12 + A_str_web * (y_bar - y_str_web)**2
    I_yy_str_flange_top = A_str_flange_top*(y_bar - y_str_flange_top)**2
    
    I_yy = I_yy_skin + I_yy_str_flange_base + I_yy_str_web + I_yy_str_flange_top
    A_stiff = A_str_flange_base + A_str_web + A_str_flange_top
    
elif stringer_type == 1: #for I stringer
    A_skin = 30*t_skin*t_skin
    y_skin = 0
    A_str_flange_base = t_stiff * w_stiff
    y_str_flange_base = 0
    A_str_web = t_stiff*h_stiff
    y_str_web = 0.5*h_stiff
    A_str_flange_top = t_stiff * w_stiff
    y_str_flange_top = h_stiff
    
    y_bar = (A_skin * y_skin + A_str_flange_base*y_str_flange_base + A_str_web*y_str_web + A_str_flange_top*y_str_flange_top) / (A_skin + A_str_flange_base + A_str_web + A_str_flange_top)
    I_yy_skin = A_skin * (y_bar)**2
    I_yy_str_flange_base = A_str_flange_base * (y_bar)**2
    I_yy_str_web = (t_stiff * h_stiff**3)/12 + A_str_web * (y_bar - y_str_web)**2
    I_yy_str_flange_top = A_str_flange_top*(y_bar - y_str_flange_top)**2
    
    I_yy = I_yy_skin + I_yy_str_flange_base + I_yy_str_web + I_yy_str_flange_top
    A_stiff = A_str_flange_base + A_str_web + A_str_flange_top
    
elif stringer_type == 2: #for Z-stringer
    A_skin = 30*t_skin*t_skin
    y_skin = 0
    A_str_flange_base = t_stiff * w_stiff
    y_str_flange_base = 0
    A_str_web = t_stiff*h_stiff
    y_str_web = 0.5*h_stiff
    A_str_flange_top = t_stiff * w_stiff
    y_str_flange_top = h_stiff
    
    y_bar = (A_skin * y_skin + A_str_flange_base*y_str_flange_base + A_str_web*y_str_web + A_str_flange_top*y_str_flange_top) / (A_skin + A_str_flange_base + A_str_web + A_str_flange_top)
    I_yy_skin = A_skin * (y_bar)**2
    I_yy_str_flange_base = A_str_flange_base * (y_bar)**2
    I_yy_str_web = (t_stiff * h_stiff**3)/12 + A_str_web * (y_bar - y_str_web)**2
    I_yy_str_flange_top = A_str_flange_top*(y_bar - y_str_flange_top)**2
    
    I_yy = I_yy_skin + I_yy_str_flange_base + I_yy_str_web + I_yy_str_flange_top
    A_stiff = A_str_flange_base + A_str_web + A_str_flange_top
    
else: #for T_stringer
    A_skin = 30*t_skin*t_skin
    y_skin = 0
    A_str_web = t_stiff*h_stiff
    y_str_web = 0.5*h_stiff
    A_str_flange_top = t_stiff * w_stiff
    y_str_flange_top = h_stiff
    
    y_bar = (A_skin * y_skin + A_str_web*y_str_web + A_str_flange_top*y_str_flange_top) / (A_skin + A_str_web + A_str_flange_top)
    I_yy_skin = A_skin * (y_bar)**2
    I_yy_str_web = (t_stiff * h_stiff**3)/12 + A_str_web * (y_bar - y_str_web)**2
    I_yy_str_flange_top = A_str_flange_top*(y_bar - y_str_flange_top)**2
    
    I_yy = I_yy_skin + I_yy_str_web + I_yy_str_flange_top
    A_stiff = A_str_web + A_str_flange_top

P_critital = np.pi**2 * E * I_yy / (a**2)

#Normalizing everything
t_star = (A_stiff/b) + t_skin
Load_crit_skin = sigma_crit_skin*t_star
Load_cr_stiff = sigma_cr_stiff*t_star
Load_rivet = sigma_rivet*t_star
Load_critical = P_critital/b

load_names = ('Skin critical', 'Stiffener critical', 'Rivet critical', 'Load critical')
y_pos = np.arange(len(load_names))
performance = [Load_crit_skin, Load_cr_stiff, Load_rivet, Load_critical]

plt.bar(y_pos, performance, align='center', alpha=0.5)
plt.xticks(y_pos, load_names)
plt.xlabel('Different scenarios')
plt.ylabel('Load [N/mm]')
plt.hlines(y=Limit_load, xmin=-0.5, xmax=3.5, linestyles='dashed')
plt.text(3.5, Limit_load, 'Limit load', ha='left', va='center')

plt.hlines(y=Ultimate_load, xmin=-0.5, xmax=3.5, colors = 'r')
plt.text(3.5, Ultimate_load, 'Ultimate load', ha='left', va='center')
plt.title('Experienced load per scenario compared to allowable loads')

plt.show()











