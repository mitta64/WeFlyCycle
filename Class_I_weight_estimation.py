# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 13:37:03 2020

@author: Matthew
"""
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

"""
Imports the NO-ENG till the Take-off length columns
from the spreadsheet:'List of 100-120 pax aircraft flying on fossil fuel'

Column 0: Number of engines
Column 1: Number of seats
Column 2: Range [km]
Column 3: Operational empty weight [kg]
Column 4: Maximum take-off weight [kg]
Column 5: Ratio of OEW/MTOW
Column 6: Thrust [kN]
Column 7: Take-off length [m]
Column 8: Wingspan [m]
Column 9: List price [M$]
Column 10: Program cost [M$]

"""
data = np.genfromtxt('Reference_aircraft.csv', delimiter = ",", skip_header=1, usecols = (3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13), invalid_raise = False)


"""
By performing Class I weight estimation the following parameters are the outputs:
    - Maximum take-off weight W_TO
    - Operational empty weight W_OE
    - Fuel used during mission F_used
    - Reserve fuel              F_res
    - Payload weight            

"""
# Calculate linear regression line, y = a * x + b
# Inputs are: operational empty weight and maximum take off weight of reference aircraft
# y = stats_OEW
# x = stats_MTOW
# Outputs are: the slope a, interception point b, correlation coefficient R, coefficient of determination R**2
stats_OEW = data[:, 3]
stats_MTOW = data[:, 4]
slope, interception, R, p_value, ssd = stats.linregress(stats_MTOW, stats_OEW)

# Plots the OEW versus the MTOW plus the linear regression line
plt.plot(stats_MTOW, stats_OEW, 'o')
plt.plot(stats_MTOW, slope * stats_MTOW + interception )
plt.xlabel('MTOW')
plt.ylabel('OEW')
# plt.show()


# Determine payload weight
pax = 120
mass_pax_luggage = 100          #[kg]       Initial value
g = 9.80665                     #[m/s^2]
W_p = pax * mass_pax_luggage * g



