import numpy as np 
from importlib import reload
import ParamsECO1
reload(ParamsECO1)
from ParamsECO1 import *

#NO REVIEW OF FOLLOWING CODE 
def ET_L(ss_in,b,sh,Ks,emax):
    
    if (sh < ss_in < 1):
        
        return emax*((ss_in-sh)/(1-sh)) + Ks*(ss_in)**(2*b+3) 
    
    else:    
        
        return Ks*(ss_in)**(2*b+3)
        
#NO REVIEW OF FOLLOWING CODE    
def Ks_func(ss_in,Ksm,b):

    return Ksm*(ss_in**(2*b+3))

#NO REVIEW OF FOLLOWING CODE 
def gp(psi_l,psi_50): #vulnerability curve
    #changed to 4
    return gp_max/(1+((psi_l/psi_50)**4))


#NO REVIEW OF FOLLOWING CODE 
def gs(ss_in,Zrm,Ksm,b): #soil conductance, conversion from days to seconds, and Pa to MPa
        
        #final units are m/(MPa*s)
        
    return ((Ks_func(ss_in,Ksm,b))*np.sqrt(RAI*ss_in**(-8)))/((np.pi*g*rho_w*Zrm)*1*10**-6)

#NO REVIEW OF FOLLOWING CODE 
def ET_NL(ss_in,sh,emax):
    
    if (sh < ss_in < 1):
        
        return emax*((ss_in-sh)/(1-sh))
    
    else:    
        
        return 0 

#NO REVIEW OF FOLLOWING CODE 
def gsrp(ss_in,psi_l,Zrm,psi_50,Ksm,b): 
    
    return (gs(ss_in,Zrm,Ksm,b)*gp(psi_l,psi_50)/f*LAI)/((gs(ss_in,Zrm,Ksm,b)+LAI*gp(psi_l,psi_50)/f))

#NO REVIEW OF FOLLOWING CODE 
def gw(ws_in,ss_in): #water storage conductance
    
    return 0#gw_max*(ws_in)**4
#NO REVIEW OF FOLLOWING CODE 
def gstom(psi_l,psi_min,gst_max): #stomatal conductance function
 
    if (psi_min <= psi_l <= -0.07):
        
        return gst_max*((psi_l-psi_min)/(-psi_min-0.07))
            
    elif (psi_l > -0.07):
        return gst_max
    
    else:
        
        return 0 
             
 #NO REVIEW OF FOLLOWING CODE       
def fD(Ta,Dx,qa):
            
    return 1/(1+VPD(Ta,qa)/Dx)
    
def fT(Ta):
            
    return 1 - k22*(Ta-Topt)**2

def fPAR(phi):
    
    return 1 - np.exp(-0.005*phi)
    
#NO REVIEW OF FOLLOWING CODE 
#Potentials (MPa)

def psi_w(ws_in,Tl): #storage potential

    
    psi = 0#4.2923*np.log(ws_in) - 0.1634 #from C3 plant, Nobel
            
    return psi #osmotic + turgor
    
#NO REVIEW OF FOLLOWING CODE 
def psi_s(ss_in,b,psi_s0): #soil
    
    if ss_in !=0: 
            
        return psi_s0*(1/(ss_in**b))
        
    else:
        ss_in = 0.000001
        
        return psi_s0*(1/(ss_in**b))
#m/s
def gmsa(psi_l,psi_min,gst_max,Ta,phi,qa):
    
    return (LAI*fD(Ta,Dx,qa)*fT(Ta)*fPAR(phi)*gstom(psi_l,psi_min,gst_max)*ga)/(LAI*fD(Ta,Dx,qa)*fT(Ta)*fPAR(phi)*gstom(psi_l,psi_min,gst_max)+ga)

#surface energy balance (E in mm/d with conversion from m/s)
def Tl_Evap(Tl,EE,phi,Ta):
     
    return -EE+conv*(phi/(lambda_w*rho_w) - (cp*ga*(Tl - Ta)*rho_a)/(lambda_w*rho_w))

def qsat(Tl):

    return (.622/Pa)*(asat*np.exp(((Tl-273)*bsat)/(csat+Tl-273)))

#exponential function within specific humidity
def VPD(Ta,qa):
    
    return ((Pa/.622)*qsat(Ta)-(qa*(Pa/0.622)))

    
def ex_fun(psi_l,Tl):
    
    return np.exp(((vww)*((psi_l*10**6)-g*rho_w*ch))/(R*Tl))

#full Air-leaf potential difference with conductances 

def ALP(psi_l,EE,Tl,qa,psi_min,gst_max,Ta,phi):
    
    return -EE+conv*((rho_a/(rho_w))*((qsat(Tl)*ex_fun(psi_l,Tl))-(qa))*(gmsa(psi_l,psi_min,gst_max,Ta,phi,qa)))

def Evap(ss_in,psi_l,ws_in,EE,Tl,Ksm,b,psi_s0,psi_50,Zrm):

    #mm/d converted from m/s
    if ss_in != 0 and gp(psi_l,psi_50) != 0:
            aa=((1 + ((1 - f)*gsrp(ss_in,psi_l,Zrm,psi_50,Ksm,b))/(LAI*gp(psi_l,psi_50)))) #+ (((1 - f)*(gw(ws_in,ss_in))*LAI)/(LAI*gp(psi_l))))
            bb=((((gsrp(ss_in,psi_l,Zrm,psi_50,Ksm,b))*(psi_s(ss_in,b,psi_s0) - psi_l)))) #+ ((gw(ws_in,ss_in)))*LAI*(psi_w(ws_in,Tl,DwMax) - psi_l)))
            return -EE+conv*(bb/aa)
    else: 
            aa=1
            bb=0
            return -EE+conv*(bb/aa) #added conversion

    #water storage and soil fluxes

def QW(psi_l,ws_in,ss_in,Tl,EE):
    
    return 0#conv*(((gw(ws_in,ss_in))*LAI)*(psi_w(ws_in,Tl)-(((EE/conv)*(1-f))/((gp(psi_l,psi_50))*LAI))-psi_l))

def QS(psi_l,ss_in,Tl,EE,Ksm,b,psi_s0,Zrm,psi_50):
    
    return conv*(gsrp(ss_in,psi_l,Zrm,psi_50,Ksm,b)*(psi_s(ss_in,b,psi_s0)-psi_l)) #(((EE/conv)*(1-f))/((gp(psi_l,psi_50))*LAI))


def psimin(psi_50):
    #from chen et al. 
    psimin = 0.6552*psi_50 - 1.3124

    return psimin




#Do review this code 
def gmax(psi_50):
    stom_50 = ((psimin(psi_50)-.07)/(2)) #50 percent loss of conductance - the gmax portion 

    #power law from supplementary code
    gmax= 0.01143287 * (-psi_50)**-1.585659 
    return gmax
    #calculate initial temp given gam_theta,h_init

def initial_theta(h0,gam_theta,theta_f0):
    
    theta_sqig = gam_theta*(1+Beta)/(1+2*Beta)*h0
    initial_theta = theta_sqig+theta_f0
    return initial_theta


def initial_gamma(h0,gam_theta,gam_q,qf0):
    q_sqig = 0.5 * (gam_q - (cp * gam_theta) / (lambda_w * (1 + 2 * Beta))) * h0
    initial_q = q_sqig+qf0
    
    return initial_q
