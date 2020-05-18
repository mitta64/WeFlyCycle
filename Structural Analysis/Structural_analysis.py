# -*- coding: utf-8 -*-
"""
Created on Mon May 11 17:03:15 2020

@author: Matthew
"""

import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
"""
 
    
Assumptions for wing & fuselage loading:
    - Wing is a cantilever beam
    - Shape of the wing is a tapered wingbox
    - Location of front spar
    - Location of rear spar
    - Centre of gravity is at 50% of the MAC

Height_airfoil is for a block, this height stays constant 
Upward forces are positive -> Lift
Downward forces are negative -> Weight

self.data
Column 1: Start of a section in terms of the span
Column 2: End of a section in terms of the span
Column 3: Average distributed lift force per section, positive since upward
Column 4: Moment arm, middle of each section
Column 5: Total moment per section
Column 6: Force of the engine, negative since downward
Column 7: 
Column 8: Force of the fuel per section, assumed to be located in the middle of a section
Column 9: Volume per section
Column 10: Chord length per section
Column 11: Minimal skin thickness, root is most critical

self.fus_data
Column 1: Start of a section in terms of the fuselage length
Column 2: End of a section in terms of the fuselage length
Column 3: Force of pax + luggage per section
Column 4: Force due to mass tank and hydrogen per section
Column 5: Force due to the engines per section
Column 6: Force due to the forward and aft cargo compartment per section
Column 7: Force due to the fuselage
Column 8: Force due to the wing
Column 9: Total force per section
Column 10: Moment arm from section to centre of gravity
Column 11: Moment per section about critical point
"""
class structural:
    
    def __init__(self, lift, wingspan, root_chord, tip_chord, sections, height_airfoil, material_yield,
                 safety_factor, mass_engine_1, loc_engine_1, mass_engine_2, loc_engine_2,
                 mass_fuel, front_spar, rear_spar, t_spar,
                 density, material_name, metal_wing, composite_wing, fus_sections,
                 fus_len, start_wing, mass_Htank, mass_H, start_Htank,
                 end_Htank, seat_pitch, seat_rows, seat_abreast, mass_pax, total_len_rows,
                 mass_fwd_cargo, loc_fwd_cargo, mass_aft_cargo, loc_aft_cargo, mass_engine, loc_engine,
                 mass_fus, mass_wing, pressure_difference, diameter_normal, diameter_double_bottom,
                 diameter_double_top, width_floor, width_beam, metal, metal_name, metal_yield, composite,
                 composite_name, composite_ultimate, composite_safety, t_beam, t_floor, loc_floor, loc_beam,
                 composite_symmetry_plies, height_double_bubble, metal_density, composite_density):
        # Wing
        self.lift = lift
        self.wingspan = wingspan
        self.root_chord = root_chord
        self.tip_chord = tip_chord
        self.sections = sections
        self.height = height_airfoil
        self.material_yield = material_yield
        self.safety_factor = safety_factor
        self.mass_engine_1 = mass_engine_1
        self.loc_engine_1 = loc_engine_2
        self.mass_engine_2 = mass_engine_2
        self.loc_engine_2 = loc_engine_2
        self.mass_fuel = mass_fuel
        self.front_spar = front_spar
        self.rear_spar = rear_spar
        self.t_spar = t_spar
        self.density = density
        self.material_name = material_name
        self.metal_wing = metal_wing
        self.composite_wing = composite_wing
        
        # Fuselage
        self.fus_sections = fus_sections
        self.fus_len = fus_len
        self.start_wing = start_wing
        self.mass_Htank = mass_Htank
        self.mass_H = mass_H
        self.start_Htank = start_Htank
        self.end_Htank = end_Htank
        self.seat_pitch = seat_pitch
        self.seat_rows = seat_rows
        self.seat_abreast = seat_abreast
        self.mass_pax = mass_pax
        self.total_len_rows = total_len_rows
        self.mass_fwd_cargo = mass_fwd_cargo
        self.loc_fwd_cargo = loc_fwd_cargo
        self.mass_aft_cargo = mass_aft_cargo
        self.loc_aft_cargo = loc_aft_cargo
        self.mass_engine = mass_engine
        self.loc_engine = loc_engine
        self.mass_fus = mass_fus
        self.mass_wing = mass_wing
        self.pressure_difference = pressure_difference
        self.diameter_normal = diameter_normal
        self.diameter_double_bottom = diameter_double_bottom
        self.diameter_double_top = diameter_double_top
        self.width_floor = width_floor
        self.width_beam = width_beam
        self.metal = metal
        self.metal_yield = metal_yield
        self.metal_name = metal_name
        self.composite = composite
        self.composite_name = composite_name
        self.composite_ultimate = composite_ultimate
        self.composite_safety = composite_safety
        self.t_beam = t_beam
        self.t_floor = t_floor
        self.loc_floor = loc_floor
        self.loc_beam = loc_beam
        self.composite_symmetry_plies = composite_symmetry_plies
        self.height_double_bubble = height_double_bubble
        self.metal_density = metal_density
        self.composite_density = composite_density
        
        # Gravity
        self.g = 9.81
        
    def wing_loading(self):
        # Define array size
        self.data = np.zeros((int(self.sections), 11))
                
        # Generate sections
        start_end = np.linspace(0, 0.5 * int(self.wingspan), int(self.sections) + 1)
        sec_length = start_end[1] - start_end[0]
        
        # Define highest value of elliptical lift distribution
        height = (self.lift * 4) / (np.pi * self.wingspan)
                       
        # Force(s) due to engine(s)
        force_engine_1 = - self.mass_engine_1 * self.g
        force_engine_2 = - self.mass_engine_2 * self.g
        
        # Force due to fuel, per section!
        force_fuel = (-(0.5 *self.mass_fuel / self.sections) * self.g)
        
        # Define chord length for each section depending on spar locations
        self.root_chord = (1 - ((1 - self.rear_spar) + self.front_spar)) * self.root_chord
        self.tip_chord = (1 - ((1 - self.rear_spar) + self.front_spar)) * self.tip_chord
        chord = np.linspace(self.root_chord, self.tip_chord, int(self.sections))
        
        # Generate moment arm
        arm_start = np.linspace(int(self.sections), 1, int(self.sections))
        self.arm = np.zeros((int(self.sections), int(self.sections)))
        
        # Max stress location
        y_max = self.height / 2
        
        for i in range(int(self.sections)):
            
            # Add start and end location of each section spanwise
            self.data[i, 0] = start_end[i]
            self.data[i, 1] = start_end[i + 1]
                                   
            # Compute lift force for each section
            sec_lift_start = height * math.sqrt(1 - (self.data[i,0]**2 / (0.5 * self.wingspan)**2))
            sec_lift_end = height * math.sqrt(1 - (self.data[i,1]**2 / (0.5 * self.wingspan)**2))
            sec_lift_avg = round((sec_lift_start + sec_lift_end) / 2, 3)
            self.data[i, 2] = sec_lift_avg
            
            # Compute moment arm from root
            self.data[i, 3] = (self.data[i, 0] + self.data[i, 1]) / 2
            
            # Add force due to the engine(s)
                # Engine 1
            if self.loc_engine_1 >= self.data[i, 0] and self.loc_engine_1 < self.data[i, 1]:
                self.data[i, 5] = force_engine_1
                             
                # Engine 2
            if self.loc_engine_2 >= self.data[i, 0] and self.loc_engine_2 < self.data[i, 1]:
                self.data[i, 5] = force_engine_2
            
            # Compute force due to the fuel
            self.data[i, 7] =  force_fuel
         
        for i in range(int(self.sections)):
            # Generated the correct moment arms per section
            b = int(arm_start[i])
            self.arm[i:int(self.sections), i] = self.data[:, 3][0:b]
            
            
            # Compute moment per section due to lift force, engine weight, fuel weight, SAFETY FACTOR USED
            self.data[i, 4] = sum((self.data[:, 2] * sec_length + self.data[:, 5] + self.data[:, 7]) 
                                  * self.arm[:, i]) * self.safety_factor
                    
            # Change in chord along the span
            self.data[i, 9] = chord[i]
                        
            # Compute minimal thickness of the skin
            # data[0,10] is the most critical value, does not matter if the rest is negative
            self.data[i, 10] = ((((abs(self.data[i, 4]) * y_max) / 
                                   self.material_yield) - ((1/6) * self.t_spar * self.height**3)) 
                                / ((1/2) * self.height**2 * (self.data[i, 9] - 2 * self.t_spar))) 
        # Final minimum thickness required for the skin    
        if self.metal_wing == 'yes':
            self.plate_thickness = max(self.data[:,10])
        
        if self.composite_wing == 'yes':
            self.plate_thickness = max(self.data[:,10]) * self.composite_symmetry_plies
        
        # Compute the volume of the material per section
        for i in range(int(self.sections)):
            self.data[i, 8] = (self.data[i, 9] * self.height - ((self.data[i, 9] - 2 * self.t_spar) *
                              (self.height - 2 * self.plate_thickness))) * sec_length 
        # Total material volume for a half wing
        self.volume_material = sum(self.data[:, 8])
        
        # Compute the weight of the wing
        # Multiplication by 2: wing = 2 * halfwing
        self.wing_mass = 2 * self.volume_material * self.density
        
        # Print results
        # print('The weight of the wing is:', self.wing_mass, '[kg]')
        # print('The thickness of the top and bottom panel is', 1000 * round(self.plate_thickness, 4), '[mm]')
        
        # Plot moment distribution over the half wingspan
        # plt.scatter(np.linspace(1, self.sections, self.sections), self.data[:, 4])
        # plt.show()
    def fuselage_loading(self):
        # Define array size
        self.fus_data = np.zeros((int(self.fus_sections), 11))
    
        # Generate sections
        start_end = np.linspace(0, int(self.fus_len), int(self.fus_sections) + 1)
        sec_length = start_end[1] - start_end[0]
        
        # Create critical point for moment analysis (assume centre of gravity is at the middle of the root chord)
        self.crit_point = self.start_wing + 0.5 * self.root_chord
        
         # Force pax per section
        force_pax_section = self.mass_pax * self.seat_abreast * 9.81
        
        # Force hydrogen tank + hydrogen per section
        force_hydrogen = (((self.mass_H + self.mass_Htank) * self.g) / 
                          ((math.floor(self.end_Htank) - math.floor(self.start_Htank)) / (sec_length)))
        
        # Force due to the engines attached to the fuselage (always 2 engines)
        force_engines = 2 * self.mass_engine * self.g
        
        # Force due to cargo
        force_fwd_cargo = self.mass_fwd_cargo * self.g
        force_aft_cargo = self.mass_aft_cargo * self.g
        
        # Force due to mass of fuselage from class II weight estimation
        force_fuselage = (self.mass_fus * self.g ) / (self.fus_len / sec_length)
        
        # Force due to mass of wing from class II weight estimation
        force_wing = self.mass_wing * self.g
        
        # Force due to pressurisation
         
            # Add pax + luggage force per section
        for i in range(int(seat_rows)):
            self.fus_data[i, 2] = force_pax_section
            
        for i in range(int(self.fus_sections)):
            # Add start and end location of each section, starting at the forward part of the cabin
            self.fus_data[i, 0] = start_end[i]
            self.fus_data[i, 1] = start_end[i + 1]
            
            # Add force hydrogen tank + hydrogen
            if self.fus_data[i, 0] >= math.floor(self.start_Htank):
                self.fus_data[i, 3] = force_hydrogen
            
            # Add force engine
            if self.loc_engine >= self.fus_data[i, 0] and self.loc_engine < self.fus_data[i, 1]:
                self.fus_data[i, 4] = force_engines
       
            # Add force cargo
            if self.loc_fwd_cargo >= self.fus_data[i, 0] and self.loc_fwd_cargo < self.fus_data[i, 1]:
                self.fus_data[i, 5] = force_fwd_cargo
                
            if self.loc_aft_cargo >= self.fus_data[i, 0] and self.loc_aft_cargo < self.fus_data[i, 1]:
                self.fus_data[i, 5] = force_aft_cargo
                
            # Add force fuselage
            self.fus_data[i, 6] = force_fuselage
            
            # Add force wing
            if self.crit_point >= self.fus_data[i, 0] and self.crit_point < self.fus_data[i, 1]:
                self.fus_data[i, 7] = force_wing
            
            # Add total force per section
            self.fus_data[i, 8] = (self.fus_data[i, 2] + self.fus_data[i, 3] + self.fus_data[i, 4] +
                                   self.fus_data[i, 5] + self.fus_data[i, 6] + self.fus_data[i, 7]) * self.safety_factor
            
            # Moment arm from section to critical point
            self.fus_data[i, 9] = abs(self.crit_point - ((self.fus_data[i, 1] + self.fus_data[i, 0]) / 2))
            
            # Moment per section about critical point
            self.fus_data[i, 10] = self.fus_data[i, 8] * self.fus_data[i, 9]
        
        # Total moment about critical point        
        moment_crit_point = sum(self.fus_data[:, 10])
        
        # Required moment of inertia, fuselage thicknesses and weight fuselage
        # Metal conventional fuselage
        if self.metal == 'yes' and self.width_beam == 0:
            self.Ixx = (moment_crit_point * (self.diameter_normal / 2)) / self.metal_yield
            self.t_fus = ((8 * self.Ixx - (self.t_floor * self.width_floor * (self.loc_floor - self.diameter_normal/2)**2))
                          / (np.pi * self.diameter_normal**3))
            self.fuselage_area = (2 * np.pi * (0.5 * self.diameter_normal - self.t_fus) * self.t_fus 
                                    + self.t_floor * self.width_floor)
            self.fuselage_volume = self.fuselage_area * self.fus_len
            self.fuselage_mass = self.fuselage_volume * metal_density
            
        # Metal double bubble fuselage
        if self.metal == 'yes' and self.width_beam != 0:
            self.Ixx = (moment_crit_point * ((self.diameter_double_bottom + self.diameter_double_top) / 2)) / self.metal_yield
            self.t_fus = ((self.Ixx - ((self.t_beam * self.width_beam) * (0.5 * self.height_double_bubble - self.loc_beam)**2))
                          / ((1/8) * np.pi * (self.diameter_double_bottom + self.diameter_double_top)))
            self.fuselage_area = ((2 * np.pi * (0.5 * self.diameter_double_bottom - self.t_fus) * self.t_fus) 
                                  + (2 * np.pi * (0.5 * self.diameter_double_top - self.t_fus) * self.t_fus) 
                                  + self.t_floor * self.width_floor + self.t_beam * self.width_beam)
            self.fuselage_volume = self.fuselage_area * self.fus_len
            self.fuselage_mass = self.fuselage_volume * metal_density
            
        #Composite conventional fuselage
        if self.composite == 'yes' and self.width_beam == 0:
            self.Ixx = (moment_crit_point * (self.diameter_normal / 2)) / (self.composite_ultimate * self.composite_safety)
            self.t_fus = ((8 * self.Ixx - (self.t_floor * self.width_floor * (self.loc_floor - self.diameter_normal/2)**2))
                          / (np.pi * self.diameter_normal**3)) * self.composite_symmetry_plies
            self.fuselage_area = (2 * np.pi * (0.5 * self.diameter_normal - self.t_fus) * self.t_fus 
                                    + self.t_floor * self.width_floor)
            self.fuselage_volume = self.fuselage_area * self.fus_len
            self.fuselage_mass = self.fuselage_volume * composite_density
        
        # Composite double bubble fuselage
        if self.composite == 'yes' and self.width_beam != 0:
            self.Ixx = ((moment_crit_point * ((self.diameter_double_bottom + self.diameter_double_top) / 2)) / 
                        (self.composite_ultimate * self.composite_safety))
            self.t_fus = ((self.Ixx - ((self.t_beam * self.width_beam) * (0.5 * self.height_double_bubble - self.loc_beam)**2))
                          / ((1/8) * np.pi * (self.diameter_double_bottom + self.diameter_double_top))) * self.composite_symmetry_plies
            self.fuselage_area = ((2 * np.pi * (0.5 * self.diameter_double_bottom - self.t_fus) * self.t_fus) 
                                  + (2 * np.pi * (0.5 * self.diameter_double_top - self.t_fus) * self.t_fus) 
                                  + self.t_floor * self.width_floor + self.t_beam * self.width_beam)
            self.fuselage_volume = self.fuselage_area * self.fus_len
            self.fuselage_mass = self.fuselage_volume * composite_density
                   
        # Print results
        print('The weight of the fuselage is:', self.fuselage_mass, '[kg]')
        print('The thickness of skin panels is:', 1000 * round(self.t_fus, 4), '[mm]')
            
# Import excel data
    # Wing
lift = pd.read_excel('parameters structures.xlsx', sheet_name = 'Wing')['value'][14]
wingspan = pd.read_excel('parameters structures.xlsx', sheet_name = 'Wing')['value'][0]
root_chord = pd.read_excel('parameters structures.xlsx', sheet_name = 'Wing')['value'][4]
tip_chord = pd.read_excel('parameters structures.xlsx', sheet_name = 'Wing')['value'][5]
sections = pd.read_excel('parameters structures.xlsx', sheet_name = 'Wing')['value'][15]
height_airfoil = pd.read_excel('parameters structures.xlsx', sheet_name = 'Wing')['value'][11]
material_yield = pd.read_excel('parameters structures.xlsx', sheet_name = 'Wing')['value'][16]
safety_factor = pd.read_excel('parameters structures.xlsx', sheet_name = 'Wing')['value'][18]
mass_engine_1 = pd.read_excel('parameters structures.xlsx', sheet_name = 'Wing')['value'][19]
loc_engine_1 = pd.read_excel('parameters structures.xlsx', sheet_name = 'Wing')['value'][20]
mass_fuel = pd.read_excel('parameters structures.xlsx', sheet_name = 'Wing')['value'][21]
front_spar = pd.read_excel('parameters structures.xlsx', sheet_name = 'Wing')['value'][8]
rear_spar = pd.read_excel('parameters structures.xlsx', sheet_name = 'Wing')['value'][9]
t_spar = pd.read_excel('parameters structures.xlsx', sheet_name = 'Wing')['value'][22]
density = pd.read_excel('parameters structures.xlsx', sheet_name = 'Wing')['value'][23]
material_name = pd.read_excel('parameters structures.xlsx', sheet_name = 'Wing')['value'][24]
mass_engine_2 = pd.read_excel('parameters structures.xlsx', sheet_name = 'Wing')['value'][25]
loc_engine_2 = pd.read_excel('parameters structures.xlsx', sheet_name = 'Wing')['value'][26]
metal_wing = pd.read_excel('parameters structures.xlsx', sheet_name = 'Wing')['value'][27]
composite_wing = pd.read_excel('parameters structures.xlsx', sheet_name = 'Wing')['value'][28]


    # Fuselage
fus_sections =  pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][0]
fus_len = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][17]
start_wing = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][24]
mass_Htank = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][13]
mass_H = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][14]
start_Htank = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][15]
end_Htank = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][16]
seat_pitch = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][18]
seat_rows = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][19]
seat_abreast = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][22]
mass_pax = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][23]
total_len_rows = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][20]
mass_fwd_cargo = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][25]
loc_fwd_cargo = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][26]
mass_aft_cargo = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][27]
loc_aft_cargo = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][28]
mass_engine = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][29]
loc_engine = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][30]
mass_fus = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][31]
mass_wing = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][32]
pressure_difference = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][5]
diameter_normal = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][9]
diameter_double_bottom = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][7]
diameter_double_top = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][6]
width_floor = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][12]
width_beam = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][33]
metal = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][34]
metal_name = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][35]
metal_yield = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][36]
composite = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][37]
composite_name = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][38]
composite_ultimate = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][39]
composite_safety = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][40]
t_beam = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][41]
t_floor = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][42]
loc_floor = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][43]
loc_beam = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][44]
composite_symmetry_plies = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][47]
height_double_bubble = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][10]
metal_density = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][45]
composite_density = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][46]

"""Use code"""


loading = structural(lift, wingspan, root_chord, tip_chord, sections, height_airfoil, material_yield,
                     safety_factor, mass_engine_1, loc_engine_1, mass_engine_2, loc_engine_2,
                     mass_fuel, front_spar, rear_spar, t_spar, density,
                     material_name, metal_wing, composite_wing, fus_sections, fus_len, start_wing, mass_Htank, mass_H, start_Htank,
                     end_Htank, seat_pitch, seat_rows, seat_abreast, mass_pax, total_len_rows,
                     mass_fwd_cargo, loc_fwd_cargo, mass_aft_cargo, loc_aft_cargo, mass_engine, loc_engine,
                     mass_fus, mass_wing, pressure_difference, diameter_normal, diameter_double_bottom,
                     diameter_double_top, width_floor, width_beam, metal, metal_name, metal_yield, composite,
                     composite_name, composite_ultimate, composite_safety, t_beam, t_floor, loc_floor,
                     loc_beam, composite_symmetry_plies, height_double_bubble, metal_density, composite_density)
loading.wing_loading()   
loading.fuselage_loading()
data = loading.data
arm = loading.arm

fus_data = loading.fus_data
    


# Check whether lift distribution is elliptical
# lift = 10000
# wingspan = 60
# height = (lift * 4) / (np.pi * wingspan)    
# y = []
# x = np.arange(0, 30.1, 0.1)
# for i in x:
#     z = height * math.sqrt(1 - (i**2 / (0.5 * wingspan)**2))
#     y.append(z)
# plt.plot(x, y)
# plt.show()