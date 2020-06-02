#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 31 19:39:05 2020

@author: suyiwen
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

inputs=pd.read_excel(r'inputs_outputs.xlsx', sheet_name="Inputs") #Read input excel
print(inputs)

class wing_loading_diagram:
    
    def __init__(self, inputs):
        self.begincruise_MTOW=0.99*0.99*0.995*0.98
        self.T_0=288.15 #K
        self.rho_0=1.225 #kg/m^3
        self.P_0 =101325 #Pa
        self.R= 287   #gas constant in J/kgK
        self.g=9.80665 #m/s^2
        self.gamma=1.4
        self.a=-0.0065
        self.MTOW=inputs["Dimensions"][0]
        self.W_f=inputs["Dimensions"][1]
        self.W_L=self.MTOW-self.W_f
        self.V_stall_landing=inputs["Dimensions"][2]
        self.M_cruise= inputs["Dimensions"][3]
        self.H_cruise=inputs["Dimensions"][4]
        self.N=inputs["Dimensions"][5]     #number of engines
        self.A=inputs["Dimensions"][6]
        self.CL_max_clean=inputs["Dimensions"][7]
        self.CL_max_takeoff=inputs["Dimensions"][8]
        self.CL_max_landing=inputs["Dimensions"][9]
        self.s_FL=inputs["Dimensions"][10]
        self.cj=inputs["Dimensions"][11]/3600/self.g
        self.ROC=inputs["Dimensions"][12]
        
        self.wingloadings=self.wing_loading(inputs)
        
    @staticmethod
    def wingloading(V,rho,CL):
        W_to_S = 0.5*rho*V**2*CL
        return W_to_S
    
    @staticmethod
    def ISA(T_0,P_0, rho_0, gamma,R, g,a,h):
        T=T_0+a*h
        P=P_0*(T/T_0)**(-g/a/T_0)
        rho=P/T/R
        return T,P,rho
    
    @staticmethod
    def cruise_velocity(T,M,gamma,R):
        a=np.sqrt(R*T*gamma)
        V=a*M
        return V
    
    @staticmethod
    def drag_polar(self, e,CL,CD0):
        CD=CD0+CL**2/np.pi/self.A/e    
        return CD
    
    @staticmethod
    def drag(self,A,CL,S):
        kg_to_lbs=2.20462
        c=0.0119
        d=0.7531
        S_wet=10**(c+d*np.log10(self.MTOW*kg_to_lbs))/10.7639   #in m^2
        
        # assume skin friction coefficient c_f=0.003       
        a=-2.5229
        b=1
        f=10**(a+b*np.log10(S_wet*10.7639))/10.7639
        print(S_wet,f)
        
        CD_0_clean=f/S
        CD_0_to=0.015+CD_0_clean
        CD_0_togear=CD_0_to+0.02
        CD_0_L=CD_0_clean+0.065
        CD_0_Lgear=CD_0_L+0.02
        CD_0_A=(CD_0_togear+CD_0_Lgear)/2
        e=0.8
        e_to=0.75
        e_L=0.7
        
        CD_clean=self.drag_polar(self,e,CL,CD_0_clean)
        CD_to=self.drag_polar(self,e_to,CL,CD_0_to)
        CD_togear=self.drag_polar(self,e_to,CL,CD_0_togear)
        CD_L=self.drag_polar(self,e_L,CL,CD_0_L)
        CD_Lgear=self.drag_polar(self,e_L,CL,CD_0_Lgear)
        CD_A=self.drag_polar(self,e_L,CL,CD_0_A)
        return CD_0_togear, CD_clean, CD_to,CD_togear,CD_L,CD_Lgear,CD_A
    
    @staticmethod
    def _rect_inter_inner(x1, x2):
        n1 = x1.shape[0]-1
        n2 = x2.shape[0]-1
        X1 = np.c_[x1[:-1], x1[1:]]
        X2 = np.c_[x2[:-1], x2[1:]]
        S1 = np.tile(X1.min(axis=1), (n2, 1)).T
        S2 = np.tile(X2.max(axis=1), (n1, 1))
        S3 = np.tile(X1.max(axis=1), (n2, 1)).T
        S4 = np.tile(X2.min(axis=1), (n1, 1))
        return S1, S2, S3, S4
    
    @staticmethod
    def _rectangle_intersection_(self, x1, y1, x2, y2):
        S1, S2, S3, S4 = self._rect_inter_inner(x1, x2)
        S5, S6, S7, S8 = self._rect_inter_inner(y1, y2)
    
        C1 = np.less_equal(S1, S2)
        C2 = np.greater_equal(S3, S4)
        C3 = np.less_equal(S5, S6)
        C4 = np.greater_equal(S7, S8)
    
        ii, jj = np.nonzero(C1 & C2 & C3 & C4)
        return ii, jj
    @staticmethod
    def intersection(self, x1, y1, x2, y2):
        x1 = np.asarray(x1)
        x2 = np.asarray(x2)
        y1 = np.asarray(y1)
        y2 = np.asarray(y2)
    
        ii, jj = self._rectangle_intersection_(self, x1, y1, x2, y2)
        n = len(ii)
    
        dxy1 = np.diff(np.c_[x1, y1], axis=0)
        dxy2 = np.diff(np.c_[x2, y2], axis=0)
    
        T = np.zeros((4, n))
        AA = np.zeros((4, 4, n))
        AA[0:2, 2, :] = -1
        AA[2:4, 3, :] = -1
        AA[0::2, 0, :] = dxy1[ii, :].T
        AA[1::2, 1, :] = dxy2[jj, :].T
    
        BB = np.zeros((4, n))
        BB[0, :] = -x1[ii].ravel()
        BB[1, :] = -x2[jj].ravel()
        BB[2, :] = -y1[ii].ravel()
        BB[3, :] = -y2[jj].ravel()

        for i in range(n):
            try:
                T[:, i] = np.linalg.solve(AA[:, :, i], BB[:, i])
            except:
                T[:, i] = np.Inf
    
        in_range = (T[0, :] >= 0) & (T[1, :] >= 0) & (
            T[0, :] <= 1) & (T[1, :] <= 1)
    
        xy0 = T[2:, in_range]
        xy0 = xy0.T
        return xy0[:, 0], xy0[:, 1]
    
    def wing_loading(self,inputs):
        m_to_ft=3.28084
        Pa_to_psf=0.02088543
        kts_to_mps=0.514444
        n=100
        TW_vertical =np.linspace(0.2,0.6, num=n)
        
        #sizing to stall speed req
        Wingloadings1= self.wingloading(self.V_stall_landing, self.rho_0, self.CL_max_landing)/(self.W_L/self.MTOW)
        
        #sizing to take-off
        TOP_imp=self.s_FL* m_to_ft/37.5
        Wingloadings_step=np.linspace(3000,6000,num=n)
        Wingloadings_step_imp=Wingloadings_step*Pa_to_psf
        CL_max_takeoff=np.linspace(1.8,2.4,num=4)
        Thrustloading_TO1=np.array([])
        plt.figure(figsize=(8,6))
        for i in CL_max_takeoff:
            Thrustloading_TO=Wingloadings_step_imp/i/TOP_imp
            plt.plot(Wingloadings_step, Thrustloading_TO, label='CL_max_TO=%s '%i)
            Thrustloading_TO1=np.append(Thrustloading_TO1,Thrustloading_TO)
            
        #Sizing to landing
        #Roskam
        V_A=np.sqrt(self.s_FL*m_to_ft/0.3)
        V_SL=V_A/1.3*kts_to_mps
        CL_max_landing=np.array([2.2,2.4,2.6,2.8])
        Wingloadings_landing=np.array([])
        for i in CL_max_landing:
            WS_L=self.wingloading(V_SL,self.rho_0,i)/(self.W_L/self.MTOW)
            Wingloadings_landing=np.append(Wingloadings_landing,WS_L)
            plt.plot(WS_L*np.ones(len(TW_vertical)),TW_vertical,label='C_L_max_L=%s '%i)
        print(Wingloadings_landing)
        #Sizing to EASA climb
        N=int(self.N)
        S=125 #assume in m^2
        CGR=0
        
        if N==2:
            CGR=0.021
        if N==3:
            CGR=0.024
        if N==4:
            CGR=0.027
              
        CL_A=2.4
        CL=CL_A/1.5**2
        CD_0_togear, CD_clean, CD_to,CD_togear,CD_L,CD_Lgear,CD_A=self.drag(self,self.A,CL,S)
        
        TW_EASA=N/(N-1)*(CD_A/CL+CGR)/0.8*0.84
        
        #Sizing to cruise speed
        T,P,rho=self.ISA(self.T_0,self.P_0, self.rho_0, self.gamma,self.R, self.g,self.a,self.H_cruise)    
        e=0.8
        CD_0=0.018
        V_cruise=self.cruise_velocity(T,self.M_cruise,self.gamma,self.R)
        print(V_cruise)
        A=np.linspace(9.5,11.5,num=5)
        for i in A:
            TW_TO = (CD_0*0.5*rho*V_cruise**2/Wingloadings_step+Wingloadings_step*2/np.pi/i/e/rho/V_cruise**2)/0.23  
            plt.plot(Wingloadings_step, TW_TO, label='Cruise speed A=%s '%i)
        
        #Sizing to climb rate requirement
        #Assum climb takes place at max L/D and gear down and take-off flaps
        e_to=0.75
        Thrustloading_TO2=np.array([])
        for i in A:
            V_ROC=np.sqrt(2*Wingloadings_step/self.rho_0* np.sqrt( CD_0_togear*np.pi*i*e_to))
            TW_TO_climb=self.ROC/V_ROC+2/np.sqrt(np.pi*i*e_to/CD_0_togear)
            plt.plot(Wingloadings_step, TW_TO_climb, label='ROC req A=%s '%i)
            Thrustloading_TO2=np.append(Thrustloading_TO2, TW_TO_climb)
        
        
        plt.plot(Wingloadings1*np.ones(len(TW_vertical)),TW_vertical,label="Stall speed landing")
        plt.plot(Wingloadings_step,TW_EASA*np.ones(len(Wingloadings_step)),label="EASA climb" )
        plt.title('Thrus-Wing Loading Diagram',fontsize=22)
        plt.xlabel("Wing loading [$W/m^2$]",fontsize=18)
        plt.ylabel("Thrust loading [-]",fontsize=18)
        plt.tight_layout()
        plt.legend(loc="upper left",bbox_to_anchor=(1, 1))
        plt.tick_params(labelsize=15)
                
        #for CL_max_L=2.8, CL_max_to=2.2
        print(len(Wingloadings_step),len(Thrustloading_TO1))
        x, y = self.intersection(self, Wingloadings_step, Thrustloading_TO1[2*len(Wingloadings_step):3*len(Wingloadings_step)], Wingloadings_landing[-1]*np.ones(len(TW_vertical)), TW_vertical)
        plt.plot(x, y, '*g',markersize=30)
        plt.show()
        #plt.savefig("Wing_loading_diagram")
        
        T_TO=y*self.MTOW*self.g/1000
        
        return x,y, T_TO
        
if __name__ == "__main__":
    results=wing_loading_diagram(inputs)
    print(results.wingloadings)
    import openpyxl
    
    file = openpyxl.load_workbook('inputs_outputs.xlsx')
    sheet= file['Outputs'] 
    for i in range(2,5):
        cell= sheet.cell(row=i, column=2)    
        cell.value=results.wingloadings[i-2].item()
    file.save('inputs_outputs.xlsx')
    
    
        
        
