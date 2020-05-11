import pandas as pd
import numpy as np
from math import *
import matplotlib.pyplot as plt
import sys

            ##### Create loading diagram #####

#Loading passengers

    # idx = FlightData[FlightData['time']==time].index.values[0]
    # th0.append(FlightData['Pitch_angle'][FlightData[FlightData['time']==time].index.values[0]]*(pi/180))

class loading:
    def __init__(self, Title, OEW, Payload, Wingloc, seat_pitch = 1, n_rows = 120/6):
        self.oew_masses = OEW['element mass']
        self.oew_locs = OEW['element location']
        
        self.payload_masses = Payload['payload mass']
        self.payload_locs = Payload['payload location']
        
        self.oew_tot = sum(self.oew_masses)
        self.cg_oew = np.dot(np.array(self.oew_masses), np.array(self.oew_locs))/self.oew_tot
        
        self.n_rows = int(n_rows)
        self.seat_pitch = seat_pitch
        
        self.pax_mass = Payload['payload mass'][Payload[Payload['payload element']=='Passenger'].index.values[0]]
        self.first_seat_loc = Payload['payload location'][Payload[Payload['payload element']=='Passenger'].index.values[0]]
        self.seat_locs = [self.first_seat_loc + seat_pitch*i for i in range(self.n_rows)]
        
        self.fwd_mass = Payload['payload mass'][Payload[Payload['payload element']=='fwd cargo'].index.values[0]]
        self.fwd_loc = Payload['payload location'][Payload[Payload['payload element']=='fwd cargo'].index.values[0]]
        
        self.aft_mass = Payload['payload mass'][Payload[Payload['payload element']=='aft cargo'].index.values[0]]
        self.aft_loc = Payload['payload location'][Payload[Payload['payload element']=='aft cargo'].index.values[0]]
        
        self.fuel_mass = Payload['payload mass'][Payload[Payload['payload element']=='Fuel'].index.values[0]]
        self.fuel_loc = Payload['payload location'][Payload[Payload['payload element']=='Fuel'].index.values[0]]
        
        self.cargo_mass = [self.fwd_mass, self.aft_mass]
        self.cargo_locs = [self.fwd_loc, self.aft_loc]
        
        self.wingloc = Wingloc['length'][Wingloc[Wingloc['distance']=='x lemac'].index.values[0]]
        self.mac = Wingloc['length'][Wingloc[Wingloc['distance']=='mac'].index.values[0]]
        
        self.seat_locs = [self.payload_locs[0] + seat_pitch*i for i in range(self.n_rows)]
        
        self.title = 'Loading diagram of {} aircraft'.format(Title)
        
    def checkargs(self, *args):
        pars = [i for i in args]
        if len(pars) == 0:
            mass_list = [self.oew_tot]
            cg_list = [self.cg_oew]
            return mass_list, cg_list
        elif len(pars) == 2:
            mass_list = [pars[0]]
            cg_list = [pars[1]]
            return mass_list, cg_list
        else:
            print("please enter two values: mass and cg")
            sys.exit()
            
    def convert_to_mac(self, cg_list):
        if type(cg_list[0]) == list:
            converted_cg_list = []
            for cgs in cg_list:
                converted_cgs = []
                for cg in cgs:
                    converted_cgs.append((cg - self.wingloc)/self.mac*100)
                converted_cg_list.append(converted_cgs)
        else:
            converted_cg_list = []
            for cgs in cg_list:
                converted_cg_list.append((cgs - self.wingloc)/self.mac*100)
        return converted_cg_list
        
            
        

    
    def ftbcargo(self, *args):        
        mass_list, cg_list = self.checkargs(*args)
        
        for i in range(2):
            mass_list.append(mass_list[-1]+self.cargo_mass[i])            
            cg_list.append((cg_list[-1]*mass_list[i] + self.cargo_locs[i]*self.cargo_mass[i])/mass_list[i+1])
        
        new_mass = mass_list[-1]
        new_cg = cg_list[-1]
        cg_list = self.convert_to_mac(cg_list)
        return [mass_list, cg_list], new_mass, new_cg
    
    def btfcargo(self, *args):        
        mass_list, cg_list = self.checkargs(*args)
        
        for i in range(2):
            mass_list.append(mass_list[-1]+self.cargo_mass[i])            
            cg_list.append((cg_list[-1]*mass_list[i] + self.cargo_locs[-i-1]*self.cargo_mass[-i-1])/mass_list[i+1])
            
        new_mass = mass_list[-1]
        new_cg = cg_list[-1]
        cg_list = self.convert_to_mac(cg_list)
        return [mass_list, cg_list], new_mass, new_cg
    
    def ftbpax(self, *args):       
        endmass, endcg = self.checkargs(*args)
        endmass = endmass[0]
        endcg = endcg[0]
        
        mass_list = []
        cg_list = []
        
        for i in range(3):
            mass_list.append([endmass + 2*self.pax_mass*i for i in range(self.n_rows+1)])
            cg_new = [endcg]
            
            for j in range(self.n_rows):
                cg_new.append((endcg*mass_list[i][j] + self.seat_locs[j]*2*self.pax_mass)/mass_list[i][j+1])
                endcg = cg_new[-1]
                
            cg_list.append(cg_new)
            endmass = mass_list[-1][-1]
            
        new_mass = mass_list[-1][-1]
        new_cg = cg_list[-1][-1]
        cg_list = self.convert_to_mac(cg_list)
        return [mass_list, cg_list], new_mass, new_cg
    
    def btfpax(self, *args):
        endmass, endcg = self.checkargs(*args)
        endmass = endmass[0]
        endcg = endcg[0]
        
        mass_list = []
        cg_list = []
        
        for i in range(3):
            mass_list.append([endmass + 2*self.pax_mass*i for i in range(self.n_rows+1)])
            cg_new = [endcg]
            
            for j in range(self.n_rows):
                cg_new.append((endcg*mass_list[i][j] + self.seat_locs[-j-1]*2*self.pax_mass)/mass_list[i][j+1])
                endcg = cg_new[-1]
                
            cg_list.append(cg_new)
            endmass = mass_list[-1][-1]
        
        new_mass = mass_list[-1][-1]
        new_cg = cg_list[-1][-1]
        cg_list = self.convert_to_mac(cg_list)
        return [mass_list, cg_list], new_mass, new_cg
    
    def fuel(self, *args):
        mass_list, cg_list = self.checkargs(*args)
        
        mass_list.append(mass_list[-1] + self.payload_masses[3])
        cg_list.append((mass_list[0]*cg_list[-1] + self.payload_masses[3]*self.payload_locs[3])/mass_list[1])
        
        new_mass = mass_list[-1]
        new_cg = cg_list[-1]
        cg_list = self.convert_to_mac(cg_list)
        return [mass_list, cg_list], new_mass, new_cg
    
    def maxima(self, margin, *args):
        cg_list = []
        for arg in args:
            if type(arg[1][0]) == list:
                for cgs in arg[1]:
                    for cg in cgs:
                        cg_list.append(cg)
            else:
                for cg in arg[1]:
                    cg_list.append(cg)
                    
        minimum = min(cg_list)-margin
        maximum = max(cg_list)+margin
        return [minimum, maximum]
    
    def diagram(self, first, second, third):
        if first == 'Passengers':
            ftbpax, junk, junk = self.ftbpax()
            btfpax, new_mass, new_cg = a=self.btfpax()
        elif first == 'Cargo':
            ftbcargo, junk, junk = self.ftbcargo()
            btfcargo, new_mass, new_cg = self.btfcargo()
        elif first == 'Fuel':
            fuel, new_mass, new_cg = self.fuel()
        else:
            print('please enter: "Passengers", "Cargo" or "Fuel"')
            sys.exit()
        
        if second == 'Passengers':
            ftbpax, junk, junk = self.ftbpax(new_mass, new_cg)
            btfpax, new_mass, new_cg = a=self.btfpax(new_mass, new_cg)
        elif second == 'Cargo':
            ftbcargo, junk, junk = self.ftbcargo(new_mass, new_cg)
            btfcargo, new_mass, new_cg = self.btfcargo(new_mass, new_cg)
        elif second == 'Fuel':
            fuel, new_mass, new_cg = self.fuel(new_mass, new_cg)
        else:
            print('please enter: "Passengers", "Cargo" or "Fuel"')
            sys.exit()
            
        if third == 'Passengers':
            ftbpax, junk, junk = self.ftbpax(new_mass, new_cg)
            btfpax, new_mass, new_cg = a=self.btfpax(new_mass, new_cg)
        elif third == 'Cargo':
            ftbcargo, junk, junk = self.ftbcargo(new_mass, new_cg)
            btfcargo, new_mass, new_cg = self.btfcargo(new_mass, new_cg)
        elif third == 'Fuel':
            fuel, new_mass, new_cg = self.fuel(new_mass, new_cg)
        else:
            print('please enter: "Passengers", "Cargo" or "Fuel"')
            sys.exit()
        
        cg_maxima = self.maxima(0.5, ftbcargo, btfcargo, ftbpax, btfpax, fuel)
        print(cg_maxima)
        
        plt.plot(ftbpax[1][0], ftbpax[0][0])
        plt.plot(ftbpax[1][1], ftbpax[0][1])
        plt.plot(ftbpax[1][2], ftbpax[0][2])
        
        plt.plot(btfpax[1][0], btfpax[0][0])
        plt.plot(btfpax[1][1], btfpax[0][1])
        plt.plot(btfpax[1][2], btfpax[0][2])
        
        plt.plot(ftbcargo[1], ftbcargo[0])
        plt.plot(btfcargo[1], btfcargo[0])
        
        plt.plot(fuel[1], fuel[0])
        
        plt.plot([cg_maxima[0], cg_maxima[0]], [self.oew_tot, new_mass], color = 'black', linewidth = 0.5)
        plt.plot([cg_maxima[1], cg_maxima[1]], [self.oew_tot, new_mass], color = 'black', linewidth = 0.5)
        
        plt.title(self.title)
        plt.show()
        
        
        
hybrid_oew = pd.read_excel('Stability and Control - Input Values.xlsx', sheet_name = 'Conventional Hybrid', usecols = 'A:C').dropna()
hybrid_pay = pd.read_excel('Stability and Control - Input Values.xlsx', sheet_name = 'Conventional Hybrid', usecols = 'E:G').dropna()
hybrid_wingloc = pd.read_excel('Stability and Control - Input Values.xlsx', sheet_name = 'Conventional Hybrid', usecols = 'I:J').dropna()
hydrogen_oew = pd.read_excel('Stability and Control - Input Values.xlsx', sheet_name = 'Conventional Hydrogen', usecols = 'A:C').dropna()
hydrogen_pay = pd.read_excel('Stability and Control - Input Values.xlsx', sheet_name = 'Conventional Hydrogen', usecols = 'E:G').dropna()
hydrogen_wingloc = pd.read_excel('Stability and Control - Input Values.xlsx', sheet_name = 'Conventional Hydrogen', usecols = 'I:J').dropna()
dubble_bubble_oew = pd.read_excel('Stability and Control - Input Values.xlsx', sheet_name = 'Dubble Bubble', usecols = 'A:C').dropna()
dubble_bubble_pay = pd.read_excel('Stability and Control - Input Values.xlsx', sheet_name = 'Dubble Bubble', usecols = 'E:G').dropna()
dubble_bubble_wingloc = pd.read_excel('Stability and Control - Input Values.xlsx', sheet_name = 'Dubble Bubble', usecols = 'I:J').dropna()


hybrid = loading('hybrid', hybrid_oew, hybrid_pay, hybrid_wingloc)
hydrogen = loading('hydrogen', hydrogen_oew, hydrogen_pay, hydrogen_wingloc)
dubble_bubble = loading('dubble bubble', dubble_bubble_oew, dubble_bubble_pay, dubble_bubble_wingloc)

hybrid.diagram('Cargo', 'Passengers', 'Fuel')
hydrogen.diagram('Cargo', 'Passengers', 'Fuel')
dubble_bubble.diagram('Cargo', 'Passengers', 'Fuel')
    