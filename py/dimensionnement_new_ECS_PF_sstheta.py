# -*- coding: utf-8 -*-
"""
Created on Fri Nov 21 10:19:52 2014

@author: yleguennec
"""

import numpy as np



Constant_vec = {}
# Simulation parameters values
Constant_vec["Ma"] = 0.0;                                 # Aircraft mach number
Constant_vec["pax"] = 120. ;                              # Number of passengers + aircraft crew
Constant_vec["P_eq"] = 7800.;                             # Equipments + outside flow dissipation
Constant_vec["P_pax"] = (70.+100.)/2.;                    # Thermal dissipation by passenger
Constant_vec["P_HT"] =  Constant_vec.get("pax")* Constant_vec.get("P_pax")+ Constant_vec.get("P_eq") ; # Thermal power to dissipate
Constant_vec["T_consigne"] = 273. + 24.;                    # Cabine temperature
Constant_vec["T5_min"] = 273. + 15.;                        # Min ECS output flow temperature
Constant_vec["T5_max"] = 273. + 30.;                        # Max ECS output flow temperature
Constant_vec["P5_min"] = 101300.;                          # Min ECS output flow pressure
Constant_vec["P5_max"] = 1.05* Constant_vec.get("P5_min");               # Max ECS output flow pressure
Constant_vec["T1"] = 473.;                                # Bleed tempearture in K
Constant_vec["P1"] = 260000.;                             # Bleed pressure in Pa
Constant_vec["T1r"] = 323.;                           # Ram air inlet temperature
Constant_vec["P1r"] = 101300.;                           # Ram air inlet pressure
Constant_vec["nu_is"] = 0.8;                              # Turbomachines isentropic efficiency
Constant_vec["DPhx"] = 40000.;                            # HX pressure losses
#Constant_vec["theta"] = 0.;                                # Valve opening

#constants
Constant_vec["cp"]=1004. # Air cp in J/kg.K
Constant_vec["r"]=287.
Constant_vec["gamma"]=1.4 #air isentropic coefficient
Constant_vec["rho_acier"] = 7800.
Constant_vec["rho"] =  Constant_vec.get("P1")/ Constant_vec.get("r")/ Constant_vec.get("T1")   # Bleed air density
Constant_vec["rho_r"] =  Constant_vec.get("P1r")/ Constant_vec.get("r")/ Constant_vec.get("T1r")   # Ram air density
Constant_vec["etam"] = 1.                                 # Shaft mechanical efficiency
Constant_vec["rp"] = 0.5                                 # Compressor/Turbine power ratio
Constant_vec["Ar"] = 0.20                                # Ram stream cross surface


# Heat exchangers parameters values
Constant_vec["rho_HX"] = 1415.                           # Heat exchanger density
Constant_vec["mu"] = 2.286e-5                            # Viscosity for main stream
Constant_vec["mur"] = 2.286e-5                           # Viscosity for ram stream
Constant_vec["delta"] = 0.102e-3                         # Fin thickness
Constant_vec["tw"] = 6.e-4                               # Wall thickness
Constant_vec["kw"] = 237.                                 # Thermal conductivity
Constant_vec["b"] = 5.21e-3                              # Plate spacing main stream
Constant_vec["br"] = 12.3e-3                             # Plate spacing ram stream
Constant_vec["beta"] = 2231.                              # Ratio between heat transfer area and volume between plates for main stream
Constant_vec["beta_r"] = 1115.                            # Ratio between heat transfer area and volume between plates for ram stream
Constant_vec["Pr"] = 0.70                                # Prandtl number for main stream
Constant_vec["Prr"] = 0.70                               # Prandtl number for ram stream
Constant_vec["Dh"] = 1.537e-3                            # Hydraulic diameter for main stream 
Constant_vec["Dhr"] = 3.41e-3;                            # Hydraulic diameter for ram stream
Constant_vec["lam"] = 0.0350                            # Constant in heat transfer coefficient calculation for main stream
Constant_vec["lamr"] = 0.0350                            # Constant in heat transfer coefficient calculation for main stream
Constant_vec["rA"] = 0.841                               # Ratio between fin and total area for main stream                            
Constant_vec["rAr"] = 0.862                              # Ratio between fin and total area for ram stream


# Compressor parameters values
Constant_vec["sigc"] = 0.9;                               # Compressor slip factor
Constant_vec["etac"] = 0.8;                               # Compressor adiabatic efficiency
Constant_vec["hc"]   = 0.7;                               # Compressor aspect ratio
Constant_vec["ec"]   = 0.01;                              # Compressor blade thickness
Constant_vec["Zc"]   = 21.;                                # Number of blades for compressor

# Turbine parameters values
Constant_vec["sigt"] = 0.8;                               # Turbine slip factor
Constant_vec["etat"] = 0.92;                              # Turbine adiabatic efficiency
Constant_vec["ht"]   = 0.5;                               # Turbine aspect ratio
Constant_vec["et"]   = 0.01;                              # Turbine blade thickness
Constant_vec["Zt"]   = 21.;                                # Number of blades for turbine






#WITH nlopt


def system_solver(X,Constant_vec):
  #  global compteur_infeas_1
 #   global compteur_infeas_2
    
    etac=Constant_vec.get("etac")                           # Compressor adiabatic efficiency
    etat=Constant_vec.get("etat")                            # Turbine adiabatic efficiency
    T1=Constant_vec.get("T1")                              # Bleed tempearture in K
    P1=Constant_vec.get("P1")                        # Bleed pressure in Pa
    T1r=Constant_vec.get("T1r")                            # Ram air inlet temperature
    DPhx=Constant_vec.get("DPhx")                           # HX pressure losses
    cp=Constant_vec.get("cp") # Air cp in J/kg.K
    r=Constant_vec.get("r")
    gamma=Constant_vec.get("gamma")#air isentropic coefficient
    rho=Constant_vec.get("rho")   # Bleed air density
    rho_r=Constant_vec.get("rho_r")    # Ram air density
    etam=Constant_vec.get("etam")                                  # Shaft mechanical efficiency
    Ar=Constant_vec.get("Ar")                                 # Ram stream cross surface
    #theta=0.   ######CETTE VALEUR NE DOIT PAS ETRE CHANGE. LE MODELE AVEC CHANGEMENT DE L ANGLE DOIT ETRE VALIDE / MODIFIE

    qm      = X[0]            # Bleed flowrate
    qmr     = X[1]           # Ram flowrate
    r2p     = X[2]           # Compressor inlet blade foot radius
    r2t     = X[3]           # Compressor inlet blade tip radius
    r3      = X[4]            # Compressor outlet blade radius
    b3      = X[5]            # Compressor outlet blade tip
    beta3   = X[6]         # Compressor outlet blade angle
    r5p     = X[7]           # Turbine outlet blade foot radius
    r5t     = X[8]           # Turbine outlet blade tip radius
    r4      = X[9]            # Turbine inlet blade radius
    b4      = X[10]            # Turbine inlet blade tip
    alpha4  = X[11]        # Turbine stator blade angle
    Lx1     = X[12]           # Heat exchanger 1 x length
    Ly1     = X[13]           # Heat exchanger 1 y length
    Lz1     = X[14]           # Heat exchanger 1 z length
    Lx2     = X[15]           # Heat exchanger 2 x length
    Ly2     = X[16]           # Heat exchanger 2 y length
    Lz2     = X[17]           # Heat exchanger 2 z length
    
    success = 10.    
    
    eps1 = epsilon_calculation (qm, qmr, Lx1, Ly1, Lz1)
    eps2 = epsilon_calculation ((1.)*qm, qmr, Lx2, Ly2, Lz2)
   
    # Shaft rotational speed
    Omega = (1.)*qm/(2.*np.pi*rho*r3**2) * (np.tan(beta3)/b3 + etam*np.tan(alpha4)/b4)
    E = qmr**3/(2.*rho_r**2*Ar**2)

    omega1 = Omega/2. - np.sqrt(Omega**2/4. - qmr**3/(2.*(1.)*qm*r3**2*rho_r**2*Ar**2))
    Wc1 = r3**2 * (1.)*qm * omega1**2 - (1.)**2*qm**2*np.tan(beta3)/(2.*np.pi*rho*b3) * omega1
    Wt1 = (1.)**2*qm**2*np.tan(alpha4)/(2.*np.pi*rho*b4)*omega1
    
    omega2 = Omega/2. + np.sqrt(Omega**2/4. - qmr**3/(2.*(1.)*qm*r3**2*rho_r**2*Ar**2))
    Wc2 = r3**2 * (1.)*qm * omega2**2 - (1.)**2*qm**2*np.tan(beta3)/(2*np.pi*rho*b3) * omega2
    Wt2 = (1.)**2*qm**2*np.tan(alpha4)/(2.*np.pi*rho*b4)*omega2

    # Pick up highest plausible solution
    if (Wc2 > 0) & (np.abs(etam * Wt2 - (Wc2 + E)) < 1e-8):
        omega = omega1
    elif (Wc1 > 0) & (np.abs(etam * Wt1 - (Wc1 + E)) < 1e-8):
        omega = omega2
    else:
        success = 0.
        T2 = -10
        T3 = -10
        T4 = -10
        T5 = -10
        T2r = -10
        T3r = -10
        P2 = -10
        P3 = -10
        P4 = -10
        P5 = -10
        eps1 = -10
        eps2 = -10
        omega = -10
        M2 = -10
        M3 = -10
        M4 = -10
        M5 = -10
        #print X
        #print Constant_vec
     #   compteur_infeas_1 = compteur_infeas_1+1
        return T2, T3, T4, T5, T2r, T3r, P2, P3, P4, P5, eps1, eps2, omega, M2, M3, M4, M5, success

    # Energy received by the fluid at compressor and turbine stages
    wc = r3**2*omega**2 - (1.)*qm*np.tan(beta3)/(2*np.pi*rho*b3)*omega
    wt = (1.)*qm*np.tan(alpha4)/(2.*np.pi*rho*b4)*omega
    
    # Input stagnation conditions (low mach number is assumed)
    T01  = T1
    P01  = P1
    T01r = T1r
    
    # Stagnation temperatures
    L1 = qm/qmr
    L2 = (1.)*qm/qmr
    T02  = ((1.-eps1)*T01 + eps1*(1.-eps2*L2)*T01r + eps1*eps2*L2*(wc/etac/cp)) / (1.-eps1*eps2*L1)         #erreur dans cette equation suspectee
    T03  = T02 + wc/etac/cp
    T04  = T03 + eps2*(T01r-T03)
    T05  = T04 - etat*wt/cp
    T02r = T01r + L2*(T03-T04)
    T03r = T02r + L1*(T01-T02)
    
    # Stagnation pressures
    P02 = P01 - DPhx
    P03 = P02 * (1. + etac*(T03-T02)/T02)**(gamma/(gamma-1.))
    P04 = P03 - DPhx
    P05 = P04 * (1. + 1./etat*(T05-T04)/T04)**(gamma/(gamma-1.))
    
    # Control surfaces
    S2 = np.pi*(r2t**2-r2p**2)
    S3 = 2.*np.pi*r3*b3
    S4 = 2.*np.pi*r4*b4
    S5 = np.pi*(r5t**2-r5p**2)
    
    # Fluid velocity at control surfaces
    C2 = (1.)*qm/rho/S2
    C3 = np.sqrt((r3*omega - (1.)*qm*np.tan(beta3)/(rho*S3))**2 + ((1.)*qm/(rho*S3))**2)
    C4 = (1.)*qm/(rho*S4)/np.cos(alpha4)
    C5 = (1.)*qm/rho/S5

    # Static temperatures
    T2 = T02 - C2**2/2./cp
    T3 = T03 - C3**2/2./cp
    T4 = T04 - C4**2/2./cp
    T5 = T05 - C5**2/2./cp
    T2r = T02r
    T3r = T03r
    
    # Mach at control surfaces
    M2 = C2/np.sqrt(gamma*r*T2)
    M3 = C3/np.sqrt(gamma*r*T3)
    M4 = C4/np.sqrt(gamma*r*T4)
    M5 = C5/np.sqrt(gamma*r*T5)
    
    # Static pressures
    P2 = P02/(1.+(gamma-1.)/2.*M2**2)**(gamma/(gamma-1.))
    P3 = P03/(1.+(gamma-1.)/2.*M3**2)**(gamma/(gamma-1.))
    P4 = P04/(1.+(gamma-1.)/2.*M4**2)**(gamma/(gamma-1.))
    P5 = P05/(1.+(gamma-1.)/2.*M5**2)**(gamma/(gamma-1.))           #on peut calculer P6 avec la loi des gaz parfait
    
    # Sanity check
    if any (~np.isreal([T2, T3, T4, T5, T2r, T3r, P2, P3, P4, P5, eps1, eps2, omega, M2, M3, M4, M5])):
        success = 0. 
     #   compteur_infeas_2 = compteur_infeas_2+1
    return T2, T3, T4, T5, T2r, T3r, P2, P3, P4, P5, eps1, eps2, omega, M2, M3, M4, M5, success




def epsilon_calculation (qm, qmr, Lx, Ly, Lz):

    # Ntu calculation
    Ntu = Ntu_calculation (qm, qmr, Lx, Ly, Lz)
    if Ntu < 0:
        eps = NaN
        return eps
        
    # Epsilon calculation
    L = qm/qmr
    # assert (L < 1); % Modelling assumption in Ntu theory
    eps = 1. - np.exp((1./L)*(Ntu**0.22)*(np.exp(-L*(Ntu**0.78))-1.))

    return eps




def Ntu_calculation (qm, qmr, Lx, Ly, Lz):
    
    mu=Constant_vec.get("mu")                            # Viscosity for main stream
    mur=Constant_vec.get("mur")                      # Viscosity for ram stream
    delta=Constant_vec.get("delta")                         # Fin thickness
    tw=Constant_vec.get("tw")                        # Wall thickness
    kw=Constant_vec.get("kw")                            # Thermal conductivity
    b=Constant_vec.get("b")                        # Plate spacing main stream
    br=Constant_vec.get("br")                          # Plate spacing ram stream
    beta=Constant_vec.get("beta")                      # Ratio between heat transfer area and volume between plates for main stream
    beta_r=Constant_vec.get("beta_r")                       # Ratio between heat transfer area and volume between plates for ram stream
    Pr=Constant_vec.get("Pr")                           # Prandtl number for main stream
    Prr=Constant_vec.get("Prr")                          # Prandtl number for ram stream
    Dh=Constant_vec.get("Dh")                 # Hydraulic diameter for main stream 
    Dhr=Constant_vec.get("Dhr")                            # Hydraulic diameter for ram stream
    lam=Constant_vec.get("lam")                     # Constant in heat transfer coefficient calculation for main stream
    lamr=Constant_vec.get("lamr")                      # Constant in heat transfer coefficient calculation for main stream
    rA=Constant_vec.get("rA")                            # Ratio between fin and total area for main stream                            
    rAr=Constant_vec.get("rAr")                     # Ratio between fin and total area for ram stream
    cp=Constant_vec.get("cp") # Air cp in J/kg.K
    # Exchanger total volume
    V = Lx * Ly * Lz
    
    # Frontal areas for main and ram stream
    Afr = Lx * Lz
    Afrr = Ly * Lz
    
    # Ratios of total heat transfer area to total exchanger volume
    alpha = b*beta / (b+br+2.*tw)
    alphar = br*beta_r / (b+br+2.*tw)
    
    # Total heat transfer areas and wall area
    A = alpha * V
    Ar = alphar * V
    
    # Wall area is taken as the arithmetic mean of the two heat transfer areas
    Aw = (A+Ar) / 2.
    
    # Ratios between free-flow and frontal area
    sigma = alpha * Dh/4.
    sigmar = alphar * Dhr/4.
    
    # Free-flow area
    Ac = sigma * Afr
    Acr = sigmar * Afrr
    
    # Mass velocities
    G = qm / Ac
    Gr = qmr / Acr
    
    # Reynolds numbers
    Re = Dh * G / mu
    Rer = Dhr * Gr / mur
    
    # Heat transfer coefficients (using Gnielinski correlation)
    f = (0.79*np.log(Re)-1.64)**-2
    fr = (0.79*np.log(Rer)-1.64)**-2
    h = f/8.*(Re-1000.)*Pr / (1.+12.7*np.sqrt(f/8.)*(Pr**(2./3.)-1.)) * (lam/Dh)
    hr = fr/8.*(Rer-1000.)*Prr / (1.+12.7*np.sqrt(fr/8.)*(Prr**(2./3.)-1.)) * (lamr/Dhr)
    m = np.sqrt(2.*h/kw/delta)
    mr = np.sqrt(2.*hr/kw/delta)
    
    # Effectiveness calculation
    nuf = np.tanh(m*b/2.) / (m*b/2.)
    nufr = np.tanh(mr*br/2.) / (mr*br/2.)
    nu0 = 1. - rA*(1.-nuf)
    nu0r = 1. - rAr*(1.-nufr)
    
    # Ntu calculation
    Ntu = 1. / (qm*cp * (1./(nu0*h*A) + tw/(kw*Aw) + 1./(nu0r*hr*Ar)))

    return Ntu


#
#def retourval():
#    bb = Xres
#
#    Xopt = np.concatenate(([bb[9]],bb[0:2],[bb[5]*np.sqrt(2/(1+bb[8]**2))],[bb[8]*bb[5]*np.sqrt(2/(1+bb[8]**2))],[bb[2]],[bb[6]],[bb[3]],[bb[6]],[bb[7]],[bb[4]],[bb[7]]))
#    eps1 = HX.eps(HX.NUT(Xres[9],Xres[9]/Xres[0],Xres[6],Xres[3],Xres[6]),Xres[9]/Xres[0])
#    eps2 = HX.eps(HX.NUT(Xres[9],Xres[9]/Xres[0],Xres[7],Xres[4],Xres[7]),Xres[9]/Xres[0])
#    Yopt1 = Oa.Problem5(Xpb,Yini2,Constant_vec)[0] 
#    print((Yopt1[2]-(-Constant_vec[13]/(Constant_vec[0]*bb[9])+Constant_vec[14])-(Yopt1[1]-Yopt1[0]))/(Yopt1[2]))      
#    Yopt = np.concatenate((Yopt1[0:-1],[1./2. * (1-bb[5]**2/bb[1]**2)],[bb[5]],[HX.NUT(Xres[9],Xres[9]/Xres[0],Xres[6],Xres[3],Xres[6])],[HX.NUT(Xres[9],Xres[9]/Xres[0],Xres[7],Xres[4],Xres[7])],[eps1],[eps2]))
#    Yopt=np.insert(Yopt,3,[Constant_vec[14] - Constant_vec[13]/(Constant_vec[0]*bb[9])]) 
#    Da = np.sqrt(1./(ga*C)*(-(1./Yopt1[9]*np.exp(ga/2*(1-(bb[5]/bb[1])**2)*Yopt1[10]**2))**2+1))
#    Da2 = bb[9]/(bb[1]**2*Yopt1[6])*np.sqrt(r*Yopt1[2]/ga)
#    Yopt=np.insert(Yopt,11,Da)    
#    return Xopt,Yopt,Constant_vec[0:10]
#aa = retourval()
#
#
#
#
#fid.write('X for reg: '+repr(aa[0])+'\n\n')
#fid.write('Y for reg: '+repr(aa[1])+'\n\n')
#fid.write('Cons for reg: '+repr(aa[2])+'\n\n')
#fid.close()
#
