function problem = ECS_ENGOPT2016 ()

%-------------------------------------------------------------------------%
% This ECS version implements a scenario with aircraft on ground, full of %
% passengers, and T = 50°C outside. Cabin must be maintained at T5 = 24°C %
%-------------------------------------------------------------------------%

% Load parameters values
param = load_param_A320_sol();

% Variable bounds
qm_min  = param.P_HT / (param.cp*(param.T_consigne-param.T5_min)); 
qm_max = 4*qm_min;
qmr_min = qm_min     ; qmr_max = qm_max;
r2p_min  = 0.03      ; r2p_max  = 0.1;
r2t_min  = 0.04      ; r2t_max  = 0.2;
r3_min  = 0.1        ; r3_max  = 0.3;
b3_min  = 0.01       ; b3_max  = 0.1;
beta3_min  = -pi/3   ; beta3_max  = pi/3;
r5p_min  = 0.03      ; r5p_max  = 0.1;
r5t_min  = 0.04      ; r5t_max  = 0.2;
r4_min  = 0.1        ; r4_max  = 0.3;
b4_min  = 0.01       ; b4_max  = 0.1;
alpha4_min  = 0      ; alpha4_max  = pi/3;
Lx1_min = 0.025      ; Lx1_max = 0.7;
Ly1_min = 0.025      ; Ly1_max = 0.7;
Lz1_min = 0.025      ; Lz1_max = 0.7;
Lx2_min = 0.025      ; Lx2_max = 0.7;
Ly2_min = 0.025      ; Ly2_max = 0.7;
Lz2_min = 0.025      ; Lz2_max = 0.7;

% Problem structure definition
problem.name = 'ECS_v2';
problem.type = 'simu';
problem.simu = @(x)simu(x, param);
problem.F = {@(x)f1(x, param), @(x)f2(x, param)};
problem.C = {@(x)c1(x, param), @(x)c2(x, param), @(x)c3(x, param), @(x)c4(x, param),...
    @(x)c5(x, param), @(x)c6(x, param), @(x)c7(x, param), @(x)c8(x, param),...
    @(x)c9(x, param), @(x)c10(x, param), @(x)c11(x, param), @(x)c12(x, param),...
    @(x)c13(x, param), @(x)c14(x, param), @(x)c15(x, param)};
problem.nonlcon = @(x)nonlcon(x, param);
problem.objfunc = @(x)objfunc(x, param);
problem.nObj = 2;
problem.nConst = 15;

problem.domain.bounds = ...
    [qm_min,   qmr_min,...
    r2p_min,   r2t_min,   r3_min,   b3_min,   beta3_min,...
    r5p_min,   r5t_min,   r4_min,   b4_min,   alpha4_min,...
    Lx1_min,   Ly1_min,   Lz1_min,...
    Lx2_min,   Ly2_min,   Lz2_min;...
    qm_max,    qmr_max,...
    r2p_max,   r2t_max,   r3_max,   b3_max,   beta3_max,...
    r5p_max,   r5t_max,   r4_max,   b4_max,   alpha4_max,...
    Lx1_max,   Ly1_max,   Lz1_max,...
    Lx2_max,   Ly2_max,   Lz2_max];
    
problem.domain.dim = size (problem.domain.bounds, 2);
problem.domain.indicator = @(x)isindomain(x, param);

end % function ECS_base


% -------------------------------------------------------------------------
function is_indom = isindomain(x, param)

    % Design variables
    qm      = x(:,1);            % Bleed flowrate
    qmr     = x(:,2);            % Ram flowrate
    r2p     = x(:,3);            % Compressor inlet blade foot radius
    r2t     = x(:,4);            % Compressor inlet blade tip radius
    r3      = x(:,5);            % Compressor outlet blade radius
    b3      = x(:,6);            % Compressort outlet blade tip    
    beta3   = x(:,7);            % Compressor outlet blade angle
    r5p     = x(:,8);            % Turbine outlet blade foot radius
    r5t     = x(:,9);            % Turbine outlet blade tip radius
    r4      = x(:,10);           % Turbine inlet blade radius
    b4      = x(:,11);           % Turbine inlet blade tip  
    alpha4  = x(:,12);           % Turbine stator blade angle
    Lx1     = x(:,13);           % Heat exchanger 1 x length 
    Ly1     = x(:,14);           % Heat exchanger 1 y length 
    Lz1     = x(:,15);           % Heat exchanger 1 z length 
    Lx2     = x(:,16);           % Heat exchanger 2 x length            
    Ly2     = x(:,17);           % Heat exchanger 2 y length
    Lz2     = x(:,18);           % Heat exchanger 2 z length
    
    % Discrimant of second order polynomial in omega
    a = qm.*r3.^2;
    b = -qm.^2/(2*pi*param.rho).*(tan(beta3)./b3 + param.etam*tan(alpha4)./b4);
    c = qmr.^3/(param.etaf * 2*param.rho^2*param.Ar^2);
    Delta = b.^2 - 4*a.*c;
    
    % Check appartenance to simulation domain
    is_indom = (qm  < qmr) &...
        (r2p + 0.02 < r2t) & (r2t < r3) & (b3 < param.hc*r3) &...
        (r5p + 0.02 < r5t) & (r5t < r4) & (b4 < param.ht*r4) &...
        (Delta > 0) &...
        (tan(beta3)./b3 + param.etam*tan(alpha4)./b4 > 0);
    
end % function isindomain


% -------------------------------------------------------------------------
function [y, c] = simu (x, param)

% Number of experiments to perform
n = size (x, 1);

% Objectives and constraints values are initialized with NaNs
y = nan (n, 2);
c = nan (n, 15);

% Loop over experiments
for i = 1:n

    % Create var structure containing current experiment design variables
    var.qm      = x(i,1);            % Bleed flowrate
    var.qmr     = x(i,2);            % Ram flowrate
    var.r2p     = x(i,3);            % Compressor inlet blade foot radius
    var.r2t     = x(i,4);            % Compressor inlet blade tip radius
    var.r3      = x(i,5);            % Compressor outlet blade radius
    var.b3      = x(i,6);            % Compressort outlet blade tip    
    var.beta3   = x(i,7);            % Compressor outlet blade angle
    var.r5p     = x(i,8);            % Turbine outlet blade foot radius
    var.r5t     = x(i,9);            % Turbine outlet blade tip radius
    var.r4      = x(i,10);           % Turbine inlet blade radius
    var.b4      = x(i,11);           % Turbine inlet blade tip  
    var.alpha4  = x(i,12);           % Turbine stator blade angle
    var.Lx1     = x(i,13);           % Heat exchanger 1 x length 
    var.Ly1     = x(i,14);           % Heat exchanger 1 y length 
    var.Lz1     = x(i,15);           % Heat exchanger 1 z length 
    var.Lx2     = x(i,16);           % Heat exchanger 2 x length            
    var.Ly2     = x(i,17);           % Heat exchanger 2 y length
    var.Lz2     = x(i,18);           % Heat exchanger 2 z length
    
    % Compute the physical values of interest for this experiment
    [T2, T3, T4, T5, T2r, T3r, P2, P3, P4, P5, eps1, eps2, omega, M2, M3, M4, M5] = system_solver (var, param);
    
    % Check for experiment failure
    if any (isnan([T2, T3, T4, T5, T2r, T3r, P2, P3, P4, P5, eps1, eps2, omega, M2, M3, M4, M5]))
        continue
    end
    
    % Mass calculation
    mass_hx1      = param.rho_HX * (var.Lx1*var.Ly1*var.Lz1);
    mass_hx2      = param.rho_HX * (var.Lx2*var.Ly2*var.Lz2);
    mass_blades_c = param.Zc * param.rho_acier * param.ec * ((var.r3-var.r2p)*param.hc*var.r3/2 - (var.r3-var.r2t)*(param.hc*var.r3-var.b3)/2);
    mass_body_c   = param.rho_acier * pi/3 * (param.hc*var.r3^3 + param.hc*var.r2p*var.r3^2 - param.hc*var.r2p^3);
    mass_blades_t = param.Zt * param.rho_acier * param.et * ((var.r4-var.r5p)*param.ht*var.r4/2 - (var.r4-var.r5t)*(param.ht*var.r4-var.b4)/2);
    mass_body_t   = param.rho_acier * pi/3 * (param.ht*var.r4^3 + param.ht*var.r5p*var.r4^2 - param.ht*var.r5p^3);
    mass          = mass_hx1 + mass_hx2 + mass_blades_c + mass_body_c + mass_blades_t + mass_body_t;
    
    % Entropy generation rate calculation
    S = var.qm  * (log(T5/param.Ta)-(param.gamma-1)/param.gamma*log(P5/param.Pa)) + ...
        var.qmr * (log(T3r/param.Ta)-(param.gamma-1)/param.gamma*log(param.rho_r*param.r*T3r/param.Pa));
    
    % Objectives
    y(i,1) = mass;
    y(i,2) = S;
    
    % Constraints
    c(i,1)  = (param.T5_min - T5)/param.T5_min;
    c(i,2)  = (T5 - param.T_consigne)/param.T_consigne;
    c(i,3)  = (param.P5_min-P5)/param.P5_min;
    c(i,4)  = (P5 - param.P5_max)/param.P5_max;
    c(i,5)  = eps1 - 0.9;
    c(i,6)  = 0.5 - eps1;
    c(i,7)  = eps2 - 0.9;
    c(i,8)  = 0.5 - eps2;
    c(i,9)  = M2 - 0.95;
    c(i,10) = M3 - 0.95;
    c(i,11) = M4 - 0.95;
    c(i,12) = M5 - 0.95;
    c(i,13) = var.r3*omega - sqrt(param.gamma*param.r*T3);
    c(i,14) = var.r4*omega - sqrt(param.gamma*param.r*T4);
    c(i,15) = (param.P_HT - var.qm*param.cp*(param.T_consigne - T5))/param.P_HT;
end

end % function simu


%--------------------------------------------------------------------------
function [T2, T3, T4, T5, T2r, T3r, P2, P3, P4, P5, eps1, eps2, omega,...
          M2, M3, M4, M5] = system_solver (var, param)

% Parameters
rho     = param.rho;
rho_r   = param.rho_r;
etam    = param.etam;
etac    = param.etac;
etat    = param.etat;
Ar      = param.Ar;
cp      = param.cp;
T1      = param.T1;
T1r     = param.T1r;
P1      = param.P1;
DPhx    = param.DPhx;
gamma   = param.gamma;
Ar      = param.Ar;
r       = param.r;

% Design variables
qm      = var.qm;            % Bleed flowrate
qmr     = var.qmr;           % Ram flowrate
r2p     = var.r2p;           % Compressor inlet blade foot radius
r2t     = var.r2t;           % Compressor inlet blade tip radius
r3      = var.r3;            % Compressor outlet blade radius
b3      = var.b3;            % Compressor outlet blade tip
beta3   = var.beta3;         % Compressor outlet blade angle
r5p     = var.r5p;           % Turbine outlet blade foot radius
r5t     = var.r5t;           % Turbine outlet blade tip radius
r4      = var.r4;            % Turbine inlet blade radius
b4      = var.b4;            % Turbine inlet blade tip
alpha4  = var.alpha4;        % Turbine stator blade angle
Lx1     = var.Lx1;           % Heat exchanger 1 x length
Ly1     = var.Ly1;           % Heat exchanger 1 y length
Lz1     = var.Lz1;           % Heat exchanger 1 z length
Lx2     = var.Lx2;           % Heat exchanger 2 x length
Ly2     = var.Ly2;           % Heat exchanger 2 y length
Lz2     = var.Lz2;           % Heat exchanger 2 z length

% Epsilon calculation for both heat exchangers
eps1 = epsilon_calculation (qm, qmr, Lx1, Ly1, Lz1, param);
eps2 = epsilon_calculation (qm, qmr, Lx2, Ly2, Lz2, param);

% Shaft rotational speed
Omega = qm/(2*pi*rho*r3^2) * (tan(beta3)/b3 + etam*tan(alpha4)/b4);
E = qmr^3/(2*rho_r^2*Ar^2);

omega1 = Omega/2 - sqrt(Omega^2/4 - qmr^3/(2*qm*r3^2*rho_r^2*Ar^2));
Wc1 = r3.^2 .* qm .* omega1.^2 - qm.^2.*tan(beta3)./(2*pi*rho*b3) .* omega1;
Wt1 = qm^2*tan(alpha4)/(2*pi*rho*b4)*omega1;

omega2 = Omega/2 + sqrt(Omega^2/4 - qmr^3/(2*qm*r3^2*rho_r^2*Ar^2));
Wc2 = r3.^2 .* qm .* omega2.^2 - qm.^2.*tan(beta3)./(2*pi*rho*b3) .* omega2;
Wt2 = qm^2*tan(alpha4)/(2*pi*rho*b4)*omega2;

% Pick up highest plausible solution
if (Wc2 > 0) && (abs(etam * Wt2 - (Wc2 + E)) < 1e-8)
    omega = omega1;
elseif (Wc1 > 0) && (abs(etam * Wt1 - (Wc1 + E)) < 1e-8)
    omega = omega2;
else
    T2 = NaN; T3 = NaN; T4 = NaN; T5 = NaN;
    P2 = NaN; P3 = NaN; P4 = NaN; P5 = NaN;
    M2 = NaN; M3 = NaN; M4 = NaN; M5 = NaN;
    T2r = NaN; T3r = NaN; omega = NaN; 
    return
end
    
% Energy received by the fluid at compressor and turbine stages
wc = r3^2*omega^2 - qm*tan(beta3)/(2*pi*rho*b3)*omega;
wt = qm*tan(alpha4)/(2*pi*rho*b4)*omega;

% Input stagnation conditions (low mach number is assumed)
T01  = T1;
P01  = P1;
T01r = T1r;

% Stagnation temperatures
L = qm/qmr;
T02  = ((1-eps1)*T01 + eps1*(1-eps2*L)*T01r + eps1*eps2*L*(wc/etac/cp)) / (1-eps1*eps2*L);
T03  = T02 + wc/etac/cp;
T04  = T03 + eps2*(T01r-T03);
T05  = T04 - etat*wt/cp;
T02r = T01r + L*(T03-T04);
T03r = T02r + L*(T01-T02);

% Stagnation pressures
P02 = P01 - DPhx;
P03 = P02 * (1 + etac*(T03-T02)/T02)^(gamma/(gamma-1));
P04 = P03 - DPhx;
P05 = P04 * (1 + 1/etat*(T05-T04)/T04)^(gamma/(gamma-1));

% Control surfaces
S2 = pi*(r2t^2-r2p^2);
S3 = 2*pi*r3*b3;
S4 = 2*pi*r4*b4;
S5 = pi*(r5t^2-r5p^2);

% Fluid velocity at control surfaces
C2 = qm/rho/S2;
C3 = sqrt((r3*omega - qm*tan(beta3)/(rho*S3))^2 + (qm/(rho*S3))^2);
C4 = qm/(rho*S4)/cos(alpha4);
C5 = qm/rho/S5;

% Static temperatures
T2 = T02 - C2^2/2/cp;
T3 = T03 - C3^2/2/cp;
T4 = T04 - C4^2/2/cp;
T5 = T05 - C5^2/2/cp;
T2r = T02r;
T3r = T03r;

% Mach at control surfaces
M2 = C2/sqrt(gamma*r*T2);
M3 = C3/sqrt(gamma*r*T3);
M4 = C4/sqrt(gamma*r*T4);
M5 = C5/sqrt(gamma*r*T5);

% Static pressures
P2 = P02/(1+(gamma-1)/2*M2^2)^(gamma/(gamma-1));
P3 = P03/(1+(gamma-1)/2*M3^2)^(gamma/(gamma-1));
P4 = P04/(1+(gamma-1)/2*M4^2)^(gamma/(gamma-1));
P5 = P05/(1+(gamma-1)/2*M5^2)^(gamma/(gamma-1));

% Sanity check
if any (~isreal([T2, T3, T4, T5, T2r, T3r, P2, P3, P4, P5, eps1, eps2, omega, M2, M3, M4, M5])) |...
   any([T2, T3, T4, T5, T2r, T3r, P2, P3, P4, P5, eps1, eps2, omega, M2, M3, M4, M5] <= 0)
    T2 = NaN; T3 = NaN; T4 = NaN; T5 = NaN;
    P2 = NaN; P3 = NaN; P4 = NaN; P5 = NaN;
    M2 = NaN; M3 = NaN; M4 = NaN; M5 = NaN;
    T2r = NaN; T3r = NaN; omega = NaN; 
    return
end

end % function system_solver


%--------------------------------------------------------------------------
function param = load_param_A320_sol()

% Ambient conditions
param.Ta = 323;                                 % Ambient temperature
param.Pa = 101300;                              % Ambient pressure

% Simulation parameters values
param.Ma = 0.0;                                 % Aircraft mach number
param.pax = 120. ;                              % Number of passengers + aircraft crew
param.P_eq = 7800.;                             % Equipments + outside flow dissipation
param.P_pax = (70.+100.)/2.;                    % Thermal dissipation by passenger
param.P_HT = param.pax*param.P_pax+param.P_eq ; % Thermal power to dissipate
param.T_consigne = 273 + 24;                    % Cabine temperature
param.T5_min = 273 + 15;                        % Min ECS output flow temperature
param.T5_max = 273 + 30;                        % Max ECS output flow temperature
param.P5_min = 101300;                          % Min ECS output flow pressure
param.P5_max = 1.05*param.P5_min;               % Max ECS output flow pressure
param.T1 = 473.;                                % Bleed tempearture in K
param.P1 = 260000.;                             % Bleed pressure in Pa
param.T1r = param.Ta;                           % Ram air inlet temperature
param.P1r = param.Pa;                           % Ram air inlet pressure
param.DPhx = 40000.;                            % HX pressure losses
param.theta = 0;                                % Valve opening

% Miscellaneous parameters
param.cp = 1004.;                               % Thermic capacity
param.r = 287.;                                 % Perfect gaz constant
param.gamma = 1.4 ;                             % Air isentropic coefficient
param.rho = param.P1/param.r/param.T1;          % Bleed air density
param.rho_r = param.P1r/param.r/param.T1r;      % Ram air density
param.rho_acier = 7800.;                        % Steel density
param.etam = 1;                                 % Shaft mechanical efficiency
param.Ar = 0.20;                                % Ram stream cross surface
param.etaf = 0.95;                              % Ram fan efficiency

% Heat exchangers parameters values
param.rho_HX = 1415.;                           % Heat exchanger density
param.mu = 2.286e-5;                            % Viscosity for main stream
param.mur = 2.286e-5;                           % Viscosity for ram stream
param.delta = 0.102e-3;                         % Fin thickness
param.tw = 6.e-4;                               % Wall thickness
param.kw = 237;                                 % Thermal conductivity
param.b = 5.21e-3;                              % Plate spacing main stream
param.br = 12.3e-3;                             % Plate spacing ram stream
param.beta = 2231;                              % Ratio between heat transfer area and volume between plates for main stream
param.beta_r = 1115;                            % Ratio between heat transfer area and volume between plates for ram stream
param.Pr = 0.70;                                % Prandtl number for main stream
param.Prr = 0.70;                               % Prandtl number for ram stream
param.Dh = 1.537e-3;                            % Hydraulic diameter for main stream 
param.Dhr = 3.41e-3;                            % Hydraulic diameter for ram stream
param.lam = 0.0350 ;                            % Constant in heat transfer coefficient calculation for main stream
param.lamr = 0.0350 ;                           % Constant in heat transfer coefficient calculation for ram stream
param.rA = 0.841;                               % Ratio between fin and total area for main stream                            
param.rAr = 0.862;                              % Ratio between fin and total area for ram stream

% Compressor parameters values
param.sigc = 0.9;                               % Compressor slip factor
param.etac = 0.8;                               % Compressor adiabatic efficiency
param.hc   = 0.7;                               % Compressor aspect ratio
param.ec   = 0.01;                              % Compressor blade thickness
param.Zc   = 21;                                % Number of blades for compressor

% Turbine parameters values
param.sigt = 0.8;                               % Turbine slip factor
param.etat = 0.92;                              % Turbine adiabatic efficiency
param.ht   = 0.5;                               % Turbine aspect ratio
param.et   = 0.01;                              % Turbine blade thickness
param.Zt   = 21;                                % Number of blades for turbine

end % function load_param_A320_sol


%--------------------------------------------------------------------------
function eps = epsilon_calculation (qm, qmr, Lx, Ly, Lz, param)

% Ntu calculation
Ntu = Ntu_calculation (qm, qmr, Lx, Ly, Lz, param);
if Ntu < 0
    eps = NaN;
    return
end

% Epsilon calculation
L = qm./qmr;
% assert (L < 1); % Modelling assumption in Ntu theory
eps = 1 - exp((1/L).*(Ntu.^0.22).*(exp(-L.*(Ntu.^0.78))-1));

end % function epsilon


% -------------------------------------------------------------------------
function Ntu = Ntu_calculation (qm, qmr, Lx, Ly, Lz, param)

% Parameters
cp = param.cp;                % Thermic capacity
mu = param.mu;                % Viscosity for main stream
mur = param.mur;              % Viscosity for ram stream
delta = param.delta;          % Fin thickness
tw = param.tw;                % Wall thickness
kw = param.kw;                % Thermal conductivity
b = param.b;                  % Plate spacing main stream
br = param.br;                % Plate spacing ram stream
beta = param.beta;            % Ratio between heat transfer area and volume between plates for main stream
betar = param.beta_r;         % Ratio between heat transfer area and volume between plates for ram stream
Pr = param.Pr;                % Prandtl number for main stream
Prr = param.Prr;              % Prandtl number for ram stream
Dh = param.Dh;                % Hydraulic diameter for main stream 
Dhr = param.Dhr;              % Hydraulic diameter for ram stream
lam = param.lam;              % Constant in heat transfer coefficient calculation for main stream
lamr = param.lamr;            % Constant in heat transfer coefficient calculation for ram stream
rA = param.rA;                % Ratio between fin and total area for main stream  (Af/A)                          
rAr = param.rAr;              % Ratio between fin and total area for ram stream (Afr/A)

% Exchanger total volume
V = Lx * Ly * Lz;

% Frontal areas for main and ram stream
Afr = Lx * Lz;
Afrr = Ly * Lz;

% Ratios of total heat transfer area to total exchanger volume
alpha = b*beta / (b+br+2*tw);
alphar = br*betar / (b+br+2*tw);

% Total heat transfer areas and wall area
A = alpha * V;
Ar = alphar * V;

% Wall area is taken as the arithmetic mean of the two heat transfer areas
Aw = (A+Ar) / 2;

% Ratios between free-flow and frontal area
sigma = alpha * Dh/4;
sigmar = alphar * Dhr/4;

% Free-flow area
Ac = sigma * Afr;
Acr = sigmar * Afrr;

% Mass velocities
G = qm / Ac;
Gr = qmr / Acr;

% Reynolds numbers
Re = Dh * G / mu;
Rer = Dhr * Gr / mur;

% Heat transfer coefficients (using Gnielinski correlation)
f = (0.79*log(Re)-1.64).^-2;
fr = (0.79*log(Rer)-1.64).^-2;
h = f/8*(Re-1000)*Pr ./ (1+12.7*sqrt(f/8)*(Pr^(2/3)-1)) * (lam/Dh);
hr = fr/8.*(Rer-1000)*Prr ./ (1+12.7*sqrt(fr/8)*(Prr^(2/3)-1)) * (lamr/Dhr);
m = sqrt (2*h/kw/delta);
mr = sqrt (2*hr/kw/delta);

% Effectiveness calculation
nuf = tanh(m*b/2) ./ (m*b/2);
nufr = tanh(mr*br/2) ./ (mr*br/2);
nu0 = 1 - rA*(1-nuf);
nu0r = 1 - rAr*(1-nufr);

% Ntu calculation
Ntu = 1 / (qm*cp * (1/(nu0.*h.*A) + tw/(kw*Aw) + 1/(nu0r*hr*Ar)));

end % function Ntu_calculation

%--------------------------------------------------------------------------
function y1 = f1(x, param)

y = simu(x, param);
y1 = y (:, 1);

end

%--------------------------------------------------------------------------
function y2 = f2(x, param)

y = simu(x, param);
y2 = y (:, 2);

end

%--------------------------------------------------------------------------
function [c, ceq] = nonlcon(x, param)

[~, c] = simu(x, param);
ceq = [];

end

%--------------------------------------------------------------------------
function y = objfunc(x, param)

[y, ~] = simu(x, param);

end

%--------------------------------------------------------------------------
function y = c1(x, param)

[~, c] = simu(x, param);
y = c(:,1);

end

%--------------------------------------------------------------------------
function y = c2(x, param)

[~, c] = simu(x, param);
y = c(:,2);

end

%--------------------------------------------------------------------------
function y = c3(x, param)

[~, c] = simu(x, param);
y = c(:,3);

end

%--------------------------------------------------------------------------
function y = c4(x, param)

[~, c] = simu(x, param);
y = c(:,4);

end

%--------------------------------------------------------------------------
function y = c5(x, param)

[~, c] = simu(x, param);
y = c(:,5);

end

%--------------------------------------------------------------------------
function y = c6(x, param)

[~, c] = simu(x, param);
y = c(:,6);

end

%--------------------------------------------------------------------------
function y = c7(x, param)

[~, c] = simu(x, param);
y = c(:,7);

end

%--------------------------------------------------------------------------
function y = c8(x, param)

[~, c] = simu(x, param);
y = c(:,8);

end

%--------------------------------------------------------------------------
function y = c9(x, param)

[~, c] = simu(x, param);
y = c(:,9);

end

%--------------------------------------------------------------------------
function y = c10(x, param)

[~, c] = simu(x, param);
y = c(:,10);

end

%--------------------------------------------------------------------------
function y = c11(x, param)

[~, c] = simu(x, param);
y = c(:,11);

end

%--------------------------------------------------------------------------
function y = c12(x, param)

[~, c] = simu(x, param);
y = c(:,12);

end

%--------------------------------------------------------------------------
function y = c13(x, param)

[~, c] = simu(x, param);
y = c(:,13);

end

%--------------------------------------------------------------------------
function y = c14(x, param)

[~, c] = simu(x, param);
y = c(:,14);

end

%--------------------------------------------------------------------------
function y = c15(x, param)

[~, c] = simu(x, param);
y = c(:,15);

end

