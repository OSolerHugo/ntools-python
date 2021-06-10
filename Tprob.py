import numpy as np #Numeric packge
from scipy.integrate import odeint #Solve DE using Rugge-Kutta
import matplotlib.pyplot as plt #Plotting graphics
######
from mpmath import * #Math functions
import time as time #Function for counting time
from scipy.special import legendre #Legendre polinomials
####
import scipy.constants as sc #Physical constants
from sympy.physics.quantum.cg import CG
from sympy import S
###


def ClebshGordan (j1, m1, j2, m2, j3, m3):
  return (float(CG(S(2*j1)/2, S(2*m1)/2, S(2*j2)/2, S(2*m2)/2, j3, m3).doit()))


file1 = open('Tprob.in')
line = file1.readline()

lread = np.fromstring(line, dtype=float, sep=',')
A = lread[0]
Z = lread[1]
r0= lread[2]

line = file1.readline()
lamb = np.fromstring(line, dtype=float, sep=',')


line = file1.readline()
stat = np.fromstring(line, dtype=float, sep=',')

line = file1.readline()
nBe2uw= np.fromstring(line, dtype=float, sep=',')

Ia=stat[0]

if (int(Ia + 0.5) != int(Ia) ):
  k = 0.5
else:
  k = 1.0

R = r0*A**(1/3)


Ia=stat[0]
i=1
while i<= (len(stat)-1):
    Ib=stat[i]
    Be2uw= nBe2uw[i-1]

    print(Ib)

    cbab = ClebshGordan(Ia,k,2,0,Ib,k)

    Be2 = (1.2)**(2*lamb)/(4*np.pi)*(3/(lamb + 3))**(2)*A**(2/3*lamb)*Be2uw
    Bab = (2*Ib + 1 )/(2*Ia + 1)*Be2
    ME = np.sqrt((2*Ia + 1)*Bab)

    MnE = np.sqrt(Bab)/(cbab)*(-1)**((Ib - Ia - np.abs(Ia - Ib))/2)

    beta = 4*np.pi/3*(MnE)/(Z*R**(lamb))
    Def = R*beta
    Rdef = 4*np.pi/(3*Z)*R**(1-lamb)*ME

    print(u'\u03B2 = %8.4f' %beta)
    print('M(E%i) ='%lamb, '%8.4f' %ME)
    print('Mn(E%i) ='%lamb, '%8.4f' %MnE)
    print('DEF(%i) ='%lamb, '%8.4f' %Def)
    print('RDEF(%i) ='%lamb, '%8.4f' %Rdef)
    i+=1
