import numpy as np 
from matplotlib import pyplot as plt
import random as rand
from scipy.optimize import fsolve
import scipy
import math
import matplotlib
from scipy.optimize import root
from scipy.special import lambertw
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy import integrate
from scipy.integrate import cumulative_trapezoid 
from ParamsECO1 import *
from Functionsdiurnal1 import * 


#for main jupyter notebook, edit this for the day length!
length =  60 #days (length of simulation)
dt     = .01 #time step

ws0   = 0.5   #w at time = 0 

#atmospheric variable vectors

ws_vec           = np.zeros(int(length/dt))
hh           = np.zeros(int(length/dt))
qqq          = np.zeros(int(length/dt))
thetavir     = np.zeros(int(length/dt))
theta        = np.zeros(int(length/dt))
phi          = np.zeros(int(length/dt))
lcl          = np.zeros(int(length/dt))
CAPE         = np.zeros(int(length/dt))
CIN          = np.zeros(int(length/dt))
CT           = []
CAC          = []





#stochastic variables 
lamb2   =  0.0 #frequency of precipitation (1/day)
alpha  = 12.47 #precip (mm)

#main loop - rainfall 

def rain(n,Zr,sfc,alpha,lamb2,dt,length):
    rain_vec = np.zeros(int(length/dt))

    for i in range (0,int(length/dt)-1):

    
        gamma  = (n*Zr)/alpha

        r = rand.uniform(0,1)
                
        if r < (lamb2*dt):
        
            #generate rainfall value
            rain_vec[i] = np.random.exponential(1/gamma)
        
        else:
            rain_vec[i] = 0
            
    return rain_vec

#zlcl function
def zlcl(T, q, Rv=461, g=9.81, cp=1004.0, Pa=100000.0, Tref=273):
    
    R = 287
    lambda_w = 2.501e6  # Latent heat of vaporization for water (J/kg) 
    rel = 0.263 * Pa*q * np.exp(17.67 * (T - Tref) / (T - 29.65)) ** -1
    #could also use 0.622*ARDbuck(Td)/P
    
    #for checking with relative humidity
    
    #q = (rel/((0.263*101000)))*np.exp((17.67*(theta-Tref)/(theta-29.65)))

    s = lambda_w*R
    
    w1 = ((lambda_w * R * (-g / cp)) / (g * Rv * T)) * (((q*Pa) / (0.622*613.75))**(R / cp))
    
    w2 = np.exp((lambda_w * R * (-g / cp)) / (g * Rv * Tref))

    j = s / ((g*Rv* lambertw(w1*w2 ,k=-1)))
    z = T / (g / cp) + j.real
    
    return z



#simple LCL function - also an option 

def simp_LCL(theta, rel, Rv=461, g=9.81, cp=1004.0, Pa=100000.0, Tref=273, a = 125):
    q = (rel/((0.263*101000)))*np.exp((17.67*(theta-Tref)/(theta-29.65)))
    #rel = 0.263 * Pa*q * np.exp(17.67 * (theta - Tref) / (theta - 29.65)) ** -1
    T = np.arange(230,351,1)
    qs = qsat(T)
    Td = interp1d(qs,T)(q)
    zlcl = a*(theta-Td)
    
    return zlcl

#find zlcl, convert to plcl and tlcl
def tp(zlcl,T0, Rd=287):

    #assume free temp is dry adiabatic process below lcl
    tlcl=T0+zlcl*-g/cp
    
    #assume pressure is dry adiabatic process below lcl
    plcl=Pa*(tlcl/T0)**(cp/Rd)
    
    return tlcl,plcl


def adiabat(tlcl,plcl,zlcl,leng, gam_q,qf0, Rv=461, Rd=287,Tref=273):
    
    
    Z   = np.linspace(zlcl, leng, 500)
    Tf  = np.empty(Z.shape)
    Tm = np.empty(Z.shape)
    Tmp = np.empty(Z.shape)
    Tfv = np.empty(Z.shape)
    P = np.empty(Z.shape)
    qs = np.empty(Z.shape)
    
    dZ  = Z[1]-Z[0]

    #free temp at lcl
    Tm[0] = tlcl
    
    #converting temperature to 
    #virtual temp
    Tmp[0] = tlcl*(1+0.61*(qsat(tlcl)*Pa/plcl))
    
    #input pressure at lcl (Pascals)
    P[0]  = plcl
    
    #input parameters
    
    #virtual potential temp and q
    thetav = gam_theta*Z + theta_f0
    q      = gam_q*Z+qf0
    theta  = thetav/(1 + 0.61*(q))
    
    #convert back to temperature
    Tf[0]  = theta[0]/((Pa/plcl)**(Rd/cp))
    
    for i in range(len(Z) - 1):
        #use true temperature to compute pressure
        Tf[i]  = theta[i]/((Pa/P[i])**(Rd/cp))
        Tfv[i] = Tf[i]*(1+0.61*q[i])
        dPdz = -P[i]*g/Rd/Tf[i]
        P[i+1] = P[i] + dPdz*dZ
        
        qs[i] = qsat(Tm[i])*Pa/P[i]
        dTdP1=Rd*Tm[i]+lambda_w*qs[i];
        dTdP2=(cp+lambda_w**2*qs[i]/Rv/Tm[i]/Tm[i]) * P[i]
        dTdP=dTdP1/dTdP2
        Tm[i+1]=Tm[i]+dTdP*(P[i+1]-P[i])
        Tmp[i]=Tm[i]*(1+0.61*qs[i])
        #returns potential for both 
    
    return Tfv,Tmp,P,Z,qs


def RAINC(PE, tlcl, plcl, zlcl, leng, gam_q, qf0, h0, tabl, qabl,
        Rv=461, Rd=287, Tref=273,Pa=Pa):
    
        

    P0 = np.linspace(Pa, plcl, 500, endpoint=False)
    Z1 = np.linspace(zlcl, leng, 500)

    n0 = len(P0)
    n1 = len(Z1)

    dP0 = P0[1] - P0[0] 
    dZ1 = Z1[1] - Z1[0]
    
 
    Tf0 = np.zeros(n0)   # temperature (parcel/env, dry)


    Tm1   = np.zeros(n1)   # parcel moist temperature
    Tmp1  = np.zeros(n1)   # parcel virtual T
    Tenv1 = np.zeros(n1)   # environment T
    Tfv1  = np.zeros(n1)   # environment virtual T
    P1    = np.zeros(n1)   # pressure
    w1    = np.zeros(n1)   # updraft velocity
    qs1   = np.zeros(n1)   # saturation specific humidity

    thetav1 = gam_theta * Z1 + theta_f0
    q1      = gam_q * Z1 + qf0
    theta1  = thetav1 / (1.0 + 0.61 * q1)


    T0=tlcl*(P0[0]/plcl)**(R/cp)
    Tf0[0] = T0
    Tf0[:]  = T0 * (P0[:] / P0[0])**(Rd / cp)


    
    # match pressure at base of moist layer with top of dry layer
    P1[0]  = P0[-1]#plcl

    #this shoud be the start of lcl
    Tm1[0] = tlcl  # parcel temperature at LCL

    for i in range(n1 - 1):

        # environment temperature from theta profile
        Tenv1[i] = theta1[i] / ((Pa / P1[i])**(Rd / cp))
        Tfv1[i]  = Tenv1[i] * (1.0 + 0.61 * q1[i])

        # hydrostatic using environment T
        dPdz      = -P1[i] * g / (Rd * Tenv1[i])
        P1[i+1]   = P1[i] + dPdz * dZ1
        
        # moist parcel thermodynamics
        qs1[i]  = qsat(Tm1[i]) * Pa / P1[i]
        dTdP1   = Rd * Tm1[i] + lambda_w * qs1[i]
        dTdP2   = (cp + (lambda_w**2) * qs1[i] / (Rv * Tm1[i]**2)) * P1[i]
        dTdP    = dTdP1 / dTdP2

        Tmp1[i] = Tm1[i] * (1.0 + 0.61 * qs1[i])

        # update parcel temp along moist adiabat
        Tm1[i+1] = Tm1[i] + dTdP * (P1[i+1] - P1[i])
        
        if i > 0:
            cape_integrand = np.maximum(Tmp1[:i] / Tfv1[:i] - 1.0, 0.0)
            integral       = g * np.trapz(cape_integrand, Z1[:i])
            w1[i]          = np.sqrt(2.0 * integral)
        else:
            w1[i] = 0.0

    # last point for env and parcel virtual T
    Tenv1[-1] = theta1[-1] / ((Pa / P1[-1])**(Rd / cp))
    Tfv1[-1]  = Tenv1[-1] * (1.0 + 0.61 * q1[-1])
    qs1[-1]   = qsat(Tm1[-1]) * Pa / P1[-1]
    Tmp1[-1]  = Tm1[-1] * (1.0 + 0.61 * qs1[-1])

    dqs_dT = np.gradient(qs1, Z1)
    dT_dz  = 1

    integrand_P = (
        w1
        * -dqs_dT
        * dT_dz
        * PE * P1 / (Rd * Tenv1)
    )

    integral_P = np.trapz(np.maximum(integrand_P, 0.0), x=Z1)

    PRECIP = integral_P 

    cape_integrand = np.maximum(g * (Tmp1 / Tfv1 - 1.0), 0.0)
    cin_integrand = np.maximum(-g * (Tmp1 / Tfv1 - 1.0), 0.0)
    CAPE = np.trapz(cape_integrand, x=Z1)
    CIN = np.trapz(cin_integrand, x=Z1)
    PGROUND = h0 * rho_a * qabl * PE

    return PRECIP, PGROUND, Tmp1, Tfv1, CAPE


def PIN(PE, AF, tlcl, plcl, zlcl, leng, gam_q, qf0, h0, tabl, qabl,
        Rv=461, Rd=287, Tref=273,Pa=Pa):
        
    Z   = np.linspace(.1, leng, 1000)
    Tf  = np.empty_like(Z)
    Tm  = np.empty_like(Z)
    Tmp = np.empty_like(Z)
    Tfv = np.empty_like(Z)
    P   = np.empty_like(Z)
    w   = np.empty_like(Z)
    qs  = np.empty_like(Z)
    
    dZ  = Z[1]-Z[0]

    P[0]  = Pa
    Tm[0] = tabl
    Tmp[0]= tabl * (1 + 0.61*qabl)  
    thetav = gam_theta*Z + theta_f0
    q      = gam_q*Z+qf0
    theta  = thetav/(1 + 0.61*(q))
    
    lcli = np.argmin(np.abs(Z - zlcl))

    abli = np.argmin(np.abs(Z - h0))

    #up to abl
    for i in range(abli):
        Tf[i]  = theta[i]/((Pa/P[i])**(Rd/cp))
        Tfv[i] = Tf[i] * (1 + 0.61*qabl)

        # parcel
        Tm[i]  = tabl - g/cp*(Z[i])
        Tmp[i] = Tm[i] * (1 + 0.61*qabl)  

        # hydrostatic pressure
        dPdz   = -P[i]*g/(Rd*Tf[i])
        P[i+1] = P[i] + dPdz*dZ


    for i in range(abli, lcli+1):
        Tf[i]  = theta[i]/((Pa/P[i])**(Rd/cp))
        Tfv[i] = Tf[i] * (1 + 0.61*q[i])

        # parcel
        Tm[i]  = tabl - g/cp*(Z[i])
        Tmp[i] = Tm[i] * (1 + 0.61*qabl)  

        # hydrostatic pressure
        dPdz   = -P[i]*g/(Rd*Tf[i])
        P[i+1] = P[i] + dPdz*dZ



    for i in range(lcli, len(Z)-1):

        # free atmosphere
        Tf[i]  = theta[i]/((Pa/P[i])**(Rd/cp))
        Tfv[i] = Tf[i] * (1 + 0.61*q[i])
    
        # hydrostatic pressure
        dPdz = -P[i]*g/(Rd*Tf[i])
        P[i+1] = P[i] + dPdz*dZ
    
        # moist processes
        qs[i]  = qsat(Tm[i]) * Pa/P[i]
        dTdP1  = Rd*Tm[i] + lambda_w*qs[i]
        dTdP2  = (cp + lambda_w**2*qs[i]/Rv/Tm[i]/Tm[i]) * P[i]
        dTdP   = dTdP1/dTdP2
        Tmp[i] = Tm[i]*(1+0.61*qs[i])
    
        Tm[i+1] = Tm[i] + dTdP*(P[i+1]-P[i])

        
            # CAPE
    integral = g*np.trapz(np.maximum((Tmp[:i+1]/Tfv[:i+1]-1),0), Z[:i+1])
    w[i]     = np.sqrt(2*integral)

    integral_P = w[:-1] * -np.gradient(qs[:-1], Tm[:-1]) * np.gradient(Tm[:-1], Z[:-1]) \
                * PE * AF * P[:-1]/(Rd*Tf[:-1])



    V_profile = integrate.cumulative_trapezoid(g*(Tmp[:-1]/Tfv[:-1] - 1), x=Z[:-1], initial=0)

    

    PRECIP = np.trapz(np.maximum(integral_P,0), x=Z[:-1]) #mm per second from conversion using water density....
    CAPE   = np.trapz(np.maximum(g*(Tmp[:-1]/Tfv[:-1]-1),0), x=Z[:-1])

    return PRECIP, Tfv, Tmp, V_profile
    
