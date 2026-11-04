import numpy as np 

def soil_params(x):
    if x == "sand":
        n = 0.35  # soil porosity [-],
        emax = 1  # maximum value of ET [cm/day]
        ET_w = 0.01  #  value of ET [mm/day] at wilting point
        sfc = 0.35  # soil moisture at field capacity
        sw = 0.11  # wilting point
        sh = 0.08  # soil moisture at hygroscopic point
        ss = 0.33  # soil moisture at onset of plant stress
        beta = 12.1  # soil term to calculate K(s) decay as part of leakage term
        Ks = 2000  # saturated hydraulic conductivity of the soil [K(ss=1)] (mm/day)
        psi_s0 = -0.34*10**-3 #MPa 
        b      = 4.05
        

    elif x == "loamy sand":
        n = 0.42  # soil porosity [-],
        emax = 1  # maximum value of ET [cm/day]
        ET_w = 0.01  #  value of ET [mm/day] at wilting point
        sfc = 0.52  # soil moisture at field capacity
        sw = 0.11  # wilting point
        sh = 0.08  # soil moisture at hygroscopic point
        ss = 0.31  # soil moisture at onset of plant stress
        beta = 12.7  # soil term to calculate K(s) decay as part of leakage term
        Ks = 1000 # saturated hydraulic conductivity of the soil [K(ss=1)] (mm/day)
        psi_s0 = -0.17*10**-3 #MPa 
        b      = 4.38
    elif x == "sandy loam":
        n = 0.43  # soil porosity [-],
        emax = 1  # maximum value of ET [cm/day]
        ET_w = 0.01  #  value of ET [mm/day] at wilting point
        sfc = 0.56  # soil moisture at field capacity
        sw = 0.18  # wilting point
        sh = 0.14  # soil moisture at hygroscopic point
        ss = 0.46  # soil moisture at onset of plant stress
        beta = 13.8  # soil term to calculate K(s) decay as part of leakage term
        Ks = 800 # saturated hydraulic conductivity of the soil [K(ss=1)] (mm/day)
        psi_s0 = -0.70*10**-3 #MPa 
        b      = 4.90

    elif x == "loam":
        n = 0.45  # soil porosity [-],
        emax = 1  # maximum value of ET [cm/day]
        ET_w = 0.01  #  value of ET [mm/day] at wilting point
        sfc = 0.65  # soil moisture at field capacity
        sw = 0.24  # wilting point
        sh = 0.19  # soil moisture at hygroscopic point
        ss = 0.57  # soil moisture at onset of plant stress
        beta = 14.8  # soil term to calculate K(s) decay as part of leakage term
        Ks = 200  # saturated hydraulic conductivity of the soil [K(ss=1)] (mm/day)
        psi_s0 = -1.43*10**-3 #MPa 
        b      = 5.39

    elif x == "clay":
        n = 0.5  # soil porosity [-],
        Zr = 300  # Depth of the rooting zone [cm]
        emax = 1  # maximum value of ET [cm/day]
        ET_w = 0.01  #  value of ET [mm/day] at wilting point
        sfc = 1  # soil moisture at field capacity
        sw = 0.52  # wilting point
        sh = 0.47  # soil moisture at hygroscopic point
        ss = 0.78  # soil moisture at onset of plant stress
        beta = 26.8  # soil term to calculate K(s) decay as part of leakage term
        Ks = 100 # saturated hydraulic conductivity of the soil [K(ss=1)] (mm/day)
        psi_s0 = -1.82*10**-3 #MPa 
        b      = 11.4
    else:
        pass
    
    return n, emax, sfc, sw, sh, ss, Ks, psi_s0,b



#Ksm  =  Ks * 1.15741*10**-8   #m/s (saturated hydraulic conductivity)


#stomatal-atmosphere params 

asat  =   613.75 #NU 
bsat  =   17.502 #NU
csat  =   240.97 #NU
vww   =   18*10**-6 #kg/mol (molar density)
R     =   8.314 #J/mol*K (gas constant)
ga    =   0.02  #m/s (conductance to atm)
#gst_max = 0.025 #m/s (conductance through stomata, maximum)
gst_min = 5*10**-8 #m/s (from Qin et al., 2020, Frontiers in Plant Science )
rho_a =   1.225 #kg/m^3 (density of air)
rho_w =   1000  #kg/m^3 (density of water)


#Ta = 290
#qa = 0.01

#rel =(qa*((0.263*101000)))/np.exp((17.67*(Ta-273)/(Ta-29.65)))

###### plant potential params 

gp_max  = 1.17*10**-5 #m/s*MPa (maximum plant conductance)
c       = 2 #NU
d       = 8 #NU
#psi_min = -4.5  #MPa
#psi_50  = -1.7 #Fix

LAI = 8
gw_max = 0#1*10**-4

#soil params

a      = 8 
RAI    = 11.6 #m^2/m^2 #NU (Root area index)
g      = 9.81 #m/s^2 (gravity acc.)

#PAR params
k1     = 0.005 #m^2/W


#storage params
d1     = 0.028   #moles
d2     = 8       #NU
n_w    = 4 #nw2 #time constant variable 


#log layer params

Pa    = 101325   #Pascal (air pressure at z=0)
conv  = 86400000 #m/s to mm/d (conversion)
f     = 1 #NU (height of plant storage node)

#atomospheric variables

cp       = 1005 #specific heat of dry air, J/kg*K
lambda_w = 2.5 * 10**6 #J/kg, latent heat of water 
t_not    = 6   #hr 
delt     = 12  #hr, day length

#atmospheric conductance from Porporato and Yin, and Data from Kaimal and Finnigan
ch         = 10 #m 
k        = 0.40 #Von Karman constant
ua       = 1.2  #m/s, mean wind speed
dis      = 0.75*ch 
e_star   = dis/10

e_star2  = 0.2*e_star

ga          =   0.02 #(ua*(k**2))/((np.log((15-dis)/e_star))*np.log((15-dis)/e_star2))

#other variables regarding humidity 

Dx       = 1250 #Pa, VPD function variable 

#optimum temperature variables 
k22  = 0.0016 #K^-2
Topt = 298 #K

rho_a       = 1.225 #kg/m^3
rho_w       = 1000 #kg/m^3
Rd         = 287.053 # gas constant for dry air J·K–1·kg–1
Rv         = 461.5 #J/(kg·K) gas constant for water vapor 


#CARBON ASSIMILATION VARIABLES (C3 from Daly (MAKE SURE)) 
a1 = 6.4

Kc0 = 302*10**-6 #mol/mol

Hkc = 59430 #J/mol

Hko = 36000 #J/mol

HvV = 116300 #J/mol

HdV = 202900 #J/mol

HdJ = 201000 #J/mol

HvJ = 79500 #J/mol

Vc0   = 50*10**-6 #mol/m^2*s

J0    = 75*10**-6 #mol/m^2*s

Kc0 = 302*10**-6 #mol/mol

Ko0 = 256/1000 #mol/mol

K2  = 0.2  #mol CO2 per mol photons

oi  = 0.209 #mol/mol

T_0 = 293.2 #K

gam_o = 34.6 * 10 ** -6 #mol/mol

gam_1 = 0.0451 #K^-1

gam_2 = 0.000347 #K^-2


Sv    = 649 #J/mol

lambm = 550 * 10 ** -9 #m
h     = 6.6 * 10 ** -34 #J/s
Na    = 6 * 10   **  23 #mol^-1
cc    = 299792458 #m/s
ca    = 420*10**-6 #mol/mol
R     = 8.314 #J/mol*K
#simple defs 
yk = 365  * 10**-6 #mol/mol

nc = 1
rc = 0.7
aw11 = 20.3*10**-6 #mol/m^2*s
alph1 = .25
alph2 = 3
#take intercept and slope from data, due to change
intercept_T = 235.20137558545503 
slope_T     = -9237.297232178247
slope_q         = -2000
intercept_q     = 0.008




#atmospheric boundary params (from Jun)
Beta      = 0.2 #entrainment ratio
gam_T = -0.0033 #K/m
#gam_q     = -2.1*10**-6 #kg/(kg*m)
T_f0  = slope_T*gam_T+intercept_T #K

#qf0       = 0.011#slope_q*gam_q+intercept_q #kg/kg
h0        = 100

#gtt,tff = converT(gam_T,T_f0,gam_q,qf0,50,100000,1000)
gam_theta = 0.0033
theta_f0  = 298 #-3481.8*gam_theta+313.75 
Rv=461
Rd=287
RAT = (Rv/Rd-1)
