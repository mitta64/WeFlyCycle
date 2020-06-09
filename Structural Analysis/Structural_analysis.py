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
Column 7: Force of the wing mass from class II weight estimation
Column 8: Force of the fuel per section, assumed to be located in the middle of a section
Column 9: Volume per section
Column 10: Chord length per section
Column 11: Minimal skin thickness, root is most critical
Column 12: Ixx, x axis from the centre of the cross section to the right
Column 13: Iyy, y axis from the centre of the cross section to the bottom
Column 14: Izz, z axis from the centre of the fuselage to the tip of the half wing

self.fus_data
Column 1: Start of a section in terms of the fuselage length
Column 2: End of a section in terms of the fuselage length
Column 3: Force of pax + luggage per section
Column 4: Force due to mass tank and hydrogen per section
Column 5: Force due to the engines and nose landing gear per section
Column 6: Force due to the forward and aft cargo compartment per section
Column 7: Force due to the fuselage
Column 8: Force due to the wing
Column 9: Total force per section
Column 10: Moment arm from section to centre of gravity
Column 11: Moment per section about critical point
Column 12: Moment analysis per section, verification
"""
class structural:
    
    def __init__(self, lift, wingspan, root_chord, tip_chord, sections, height_airfoil, material_yield,
                 safety_factor, mass_engine_1, loc_engine_1, mass_engine_2, loc_engine_2,
                 mass_fuel, front_spar, rear_spar, t_spar,
                 density, material_name, metal_wing, composite_wing, wing_class2, 
                 mass_landing_gear, loc_landing_gear, fus_sections,
                 fus_len, start_wing, mass_Htank, mass_H, start_Htank,
                 end_Htank, seat_pitch, seat_rows, seat_abreast, mass_pax, total_len_rows,
                 mass_fwd_cargo, loc_fwd_cargo, mass_aft_cargo, loc_aft_cargo, mass_engine, loc_engine,
                 mass_fus, mass_wing, pressure_difference, diameter_normal, diameter_double_bottom,
                 diameter_double_top, width_floor, width_beam, metal, metal_name, metal_yield, composite,
                 composite_name, composite_ultimate, composite_safety, t_beam, t_floor, loc_floor, loc_beam,
                 composite_symmetry_plies, height_double_bubble, metal_density, composite_density,
                 mass_nose_gear, loc_nose_gear):
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
        self.wing_class2 = wing_class2
        self.mass_landing_gear = mass_landing_gear
        self.loc_landing_gear = loc_landing_gear
        
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
        self.mass_nose_gear = mass_nose_gear
        self.loc_nose_gear = loc_nose_gear
        
        # Gravity
        self.g = 9.81
        
    def wing_loading(self):
        # Define array size
        self.data = np.zeros((int(self.sections), 14))
                
        # Generate sections
        start_end = np.linspace(0, 0.5 * int(self.wingspan), int(self.sections) + 1)
        sec_length = start_end[1] - start_end[0]
        
        # Define highest value of elliptical lift distribution
        height = (self.lift * 4) / (np.pi * self.wingspan)
                       
        # Force(s) due to engine(s)
        force_engine_1 = - self.mass_engine_1 * self.g
        force_engine_2 = - self.mass_engine_2 * self.g
        
        # Force due to the landing gear
        force_landing_gear = -self.mass_landing_gear * self.g
        
        # Force due to fuel, per section!
        force_fuel = (-(0.5 *self.mass_fuel / self.sections) * self.g)
        
        # Force wing mass per section
        force_wing = (- self.wing_class2 * self.g ) / (self.sections)
        
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
            
            # Add force due to the engine(s) and landing gear
                # Engine 1
            if self.loc_engine_1 >= self.data[i, 0] and self.loc_engine_1 < self.data[i, 1]:
                self.data[i, 5] = force_engine_1
                             
                # Engine 2
            if self.loc_engine_2 >= self.data[i, 0] and self.loc_engine_2 < self.data[i, 1]:
                self.data[i, 5] = force_engine_2
           
                # Landing gear
            if self.loc_landing_gear >= self.data[i, 0] and self.loc_landing_gear < self.data[i, 1]:
                self.data[i, 5] = force_landing_gear
            
            #Compute force due to the wing
            self.data[i, 6] = force_wing
            # Compute force due to the fuel
            self.data[i, 7] =  force_fuel
          
        for i in range(int(self.sections)):
            # Generated the correct moment arms per section
            b = int(arm_start[i])
            self.arm[i:int(self.sections), i] = self.data[:, 3][0:b]
            
            
            # Compute moment per section due to lift force, engine weight, fuel weight, SAFETY FACTOR USED
            self.data[i, 4] = sum((self.data[:, 2] * sec_length + self.data[:, 5] + self.data[:,6] + self.data[:, 7]) 
                                  * self.arm[:, i]) * self.safety_factor
                    
            # Change in chord along the span
            self.data[i, 9] = chord[i]
                        
            # Compute minimal thickness of the skin and Ixx, Iyy
            # data[0,10] is the most critical value, does not matter if the rest is negative
            self.data[i, 10] = ((((abs(self.data[i, 4]) * y_max) / 
                                   self.material_yield) - ((1/6) * self.t_spar * self.height**3)) 
                                / ((1/2) * self.height**2 * (self.data[i, 9] - 2 * self.t_spar))) 
            # Compute moment of inertia per section: Ixx, Iyy, Izz
            self.data[i, 11] =  ((1/6) * self.t_spar * self.height**3 + (1/2) * self.height**2 * self.data[0, 10]
                                * (self.data[i, 9] - 2 * self.t_spar))
            self.data[i, 12] = ((1/6) * self.t_spar**3 * self.height + (1/6) * (self.data[i, 9] - self.t_spar)**3 
                                * self.data[0, 10] + (1/2) * self.height**2 * self.data[0, 10] * 
                                (self.data[i, 9] - 2 * self.t_spar))
            #self.data[i, 13] = 
            
        # Final minimum thickness required for the skin    
        if self.metal_wing == 'yes':
            self.plate_thickness = max(self.data[:,10])
        
        if self.composite_wing == 'yes':
            self.plate_thickness = max(self.data[:,10]) * 2 * self.composite_symmetry_plies
        
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
        print('The weight of the wing is:', self.wing_mass, '[kg]')
        print('The thickness of the top and bottom panel is', 1000 * round(self.plate_thickness, 4), '[mm]')
        
        # # Check for fatigue analysis
        # tensile_stress_root = (self.data[0, 4] / self.safety_factor * y_max) / self.data[0, 11]
        # thickness = ((((abs(self.data[0, 4]) * y_max) / 
        #                            290726196.3) - ((1/6) * self.t_spar * self.height**3)) 
        #                         / ((1/2) * self.height**2 * (self.data[0, 9] - 2 * self.t_spar))) 
        # print(1000 * thickness)
        #==============================================================================
        # Verification of wing internal moments
        #Plot moment distribution over the half wingspan
        # plt.plot(np.linspace(1, self.sections, self.sections), self.data[:, 4], marker ='o', linestyle = '-' )
        # plt.plot([1,20], [0,0])
        # plt.title('Moment distribution over the half wing with starting position the root')
        # plt.xlabel('Section number [-]')
        # plt.ylabel('Internal moment [Nm]')
        # plt.show()
        #==============================================================================
    def fuselage_loading(self):
        # Define array size
        self.fus_data = np.zeros((int(self.fus_sections), 12))
    
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
        
        # Force due to nose landing gear
        force_nose_gear = self.mass_nose_gear * self.g
         
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
       
            # Add force nose landing gear
            if self.loc_nose_gear >= self.fus_data[i, 0] and self.loc_nose_gear < self.fus_data[i, 1]:
                self.fus_data[i, 4] = force_nose_gear
        
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
        
        #==============================================================================
        # # Verification of fuselage internal moments
        
        # # Amount of left and right sections of critical point
        # self.left_sections_crit = np.where(self.fus_data[:,0] < self.crit_point)
        # self.right_sections_crit = np.where(self.fus_data[:,0] > self.crit_point)
        
        # right_fus = int(self.right_sections_crit[0][-1] - self.right_sections_crit[0][0] ) + 1
        # fus_arm_right = np.linspace(right_fus, 1, right_fus)
        # self.fus_arm_right = np.zeros((right_fus, right_fus))
        # arm_right_end = int(self.right_sections_crit[0][-1])
        # for i in range(right_fus):
        #     # Generated the correct moment arms per section
        #     b = int(fus_arm_right[i])
        #     self.fus_arm_right[i:arm_right_end, i] = self.fus_data[:, 0][0:b]
        
        # left_fus = int(self.left_sections_crit[0][-1] - self.left_sections_crit[0][0] ) 
        # fus_arm_left = np.linspace(left_fus, 1, left_fus)
        # self.fus_arm_left = np.zeros((left_fus, left_fus))
        # arm_left_end = int(self.left_sections_crit[0][-1])
        # for i in range(left_fus):
        #     # Generated the correct moment arms per section
        #     b = int(fus_arm_left[i])
        #     self.fus_arm_left[i:arm_left_end, i] = self.fus_data[:, 0][0:b]
        
        # for i in range(int(self.fus_sections)):
        #     if i < self.left_sections_crit[0][-1]:
        #         self.fus_data[i, 11] = sum( self.fus_data[:, 8][0:arm_left_end] * np.flip(self.fus_arm_left[:, (-1-i)]))
        #     if i > self.left_sections_crit[0][-1]:
        #         self.fus_data[i, 11] = sum( self.fus_data[:, 8][self.right_sections_crit[0][0] : arm_right_end + 1] * self.fus_arm_right[:, i -( self.left_sections_crit[0][-1] + 1) ])
        # self.fus_data[self.right_sections_crit[0][0] - 1, 11] = self.fus_data[self.right_sections_crit[0][0], 11]
        # plt.plot(np.linspace(1, self.fus_sections, self.fus_sections), self.fus_data[:,11],  marker ='o', linestyle = '-')
        # plt.plot([1,62], [0,0])
        # plt.title('Moment distribution over the fuselage with starting position the nose')
        # plt.xlabel('Section number [-]')
        # plt.ylabel('Internal moment [Nm]')
        # plt.show()
        #==============================================================================
            
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
            
            # Moment of inertia Izz, Iyy, Ixx 
            Izz = self.Ixx
            Iyy = self.Ixx - (self.t_floor * self.width_floor * (self.loc_floor - self.diameter_normal/2)**2)
            self.I = [Izz, Iyy]
            
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
            
            
            # Moment of inertia Izz, Iyy, Ixx 
            Izz = self.Ixx
            Iyy = self.Ixx - (self.t_floor * self.width_floor * (self.loc_floor - self.diameter_normal/2)**2)
            self.I = [Izz, Iyy]
            
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
        # print('The weight of the fuselage is:', self.fuselage_mass, '[kg]')
        # print('The thickness of skin panels is:', 1000 * round(self.t_fus, 4), '[mm]')
            
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
wing_class2 = pd.read_excel('parameters structures.xlsx', sheet_name = 'Wing')['value'][29]
mass_landing_gear = pd.read_excel('parameters structures.xlsx', sheet_name = 'Wing')['value'][30]
loc_landing_gear = pd.read_excel('parameters structures.xlsx', sheet_name = 'Wing')['value'][31]

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
mass_nose_gear = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][48]
loc_nose_gear = pd.read_excel('parameters structures.xlsx', sheet_name = 'Fuselage')['value'][49]

"""Use code"""


loading = structural(lift, wingspan, root_chord, tip_chord, sections, height_airfoil, material_yield,
                     safety_factor, mass_engine_1, loc_engine_1, mass_engine_2, loc_engine_2,
                     mass_fuel, front_spar, rear_spar, t_spar, density,
                     material_name, metal_wing, composite_wing, wing_class2, mass_landing_gear, loc_landing_gear,
                     fus_sections, fus_len, start_wing, mass_Htank, mass_H, start_Htank,
                     end_Htank, seat_pitch, seat_rows, seat_abreast, mass_pax, total_len_rows,
                     mass_fwd_cargo, loc_fwd_cargo, mass_aft_cargo, loc_aft_cargo, mass_engine, loc_engine,
                     mass_fus, mass_wing, pressure_difference, diameter_normal, diameter_double_bottom,
                     diameter_double_top, width_floor, width_beam, metal, metal_name, metal_yield, composite,
                     composite_name, composite_ultimate, composite_safety, t_beam, t_floor, loc_floor,
                     loc_beam, composite_symmetry_plies, height_double_bubble, metal_density, composite_density,
                     mass_nose_gear, loc_nose_gear)
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