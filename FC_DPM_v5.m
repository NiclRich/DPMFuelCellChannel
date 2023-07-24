%% ============ Mechanistic Channel Discrete Particle Model (MCDPM)==================

%------------------------------------------------------------------------------------

% Uses a discrete particle model approach to model two-phase flow in fuel
% cell gas channels. Uses mechanistic assumptions for droplet attachment and merging.

% Published in Niblett, Daniel, Stuart Martin Holmes, and Vahid Niasar. 
% "Discrete-particle model to optimize operational conditions of proton-exchange membrane fuel-cell gas channels." 
% ACS Applied Energy Materials 4.10 (2021): 10514-10533.

% Developed by Dr. Daniel Niblett mainly during 2019 - 2021 (PhD Student 2017 - 2021)

% https://scholar.google.com/citations?user=3fkp6J0AAAAJ&hl=en

%------------------------------------------------------------------------------------
%     Recommended Updates for others to implement:
%
%           - Update algorithms to work for any flow field (instead of straight channel)
%           - Incorporate continuum scale flow model for flow in channels and porous sections combined
%           - Use solved velocity field (taking into account water fraction) to transport O2 + H2O vapor
%           - Estimate droplet/water cluster evaporation rates using lagrangian approach (i.e. Sh and Nu correlations)
%           - Incorporate solving charge transport coupled to local reaction rates (Butler-volmer with mass transfer correction)
%
%------------------------------------------------------------------------------------

% NOTE: I created this model a while ago and so some of the programming efficiency will be low, if I was to write this model again
% I would do it differently.

% If you are interested to develop this model further, please get in touch:

% - daniel.niblett@newcastle.ac.uk
% - daniel.niblett@live.co.uk
% - www.danielniblett.com

%% Model Initialisation

clc
clear

% -------- Saving Data -------- %
  mkdir jul                       % make new folder 
  Savefilename = 'jul/run_';      % make new file name
  mkdir jul VTK                   % make folder to save .vtk outputs
  vtkname = 'jul/VTK/';           % make name of vtk files

% -------- Read external .csv to specify operating conditions ------- %
% OPC = csvread('opc_jul2021.csv');
% repeatmax = length(OPC);                 % number of model iterations
  repeatmax = 1;                          % or specify iterations


%============== Example (varying air and water velocity)===========
% uinr = OPC(:,2);                         %air velocity
% u_wr = OPC(:,1);                         %water velocity
% sens_theta_GDL = 40 - 160                % (sens)itivity analysis
% sens_theta_w = 40 - 160
% sens_uw = 1e-4 - 1
% sens_ua = 0.1 - 30
% sens_rhoa = 0.7 - 2
% sens_npores = [0.1 1 5 10 20];
% sens_sigma = 0.02 - 0.072
% snes_Rpmax = 
% channel_length = [0.01 0.05 0.1 0.2 0.4 0.8];


% Repeat running of model for repeatmax times (set to 1 for single run)
for repeat=1:1:repeatmax 
    
%% Model Visualisation and Saving Settings
sid = repeat;                        
disp(['Starting Run Number: ' num2str(sid)]) % output run number (to see on screen)

clear D Time S RegimeSatSave Regime timetest % must clear to preallocate for next run
vtk = 1;                              % Vtk plotting of results at each time step (1 for on)
vtkplot = 1;                          % how many iterations to skip for vtk saving
serpentine_plot = 0;                  % Vtk plotting in serpentine arrangement (1 for on)
coalescence = 1;                      % Droplet Collision and Merge Algorithm (1 for on)                   
transport = 0;                        % Transport Algorithm on (1 for on) (IN PROGRESS)
evaporation = 0;                      % Evaporation Algorithm on (1 for on) (IN PROGRESS)
condensation = 0;                     % Condensation Algorithm (IN PROGRESS)
film_formation = 1;                   % Film Formation Regime (1 for on)
plug_formation = 1;                   % Plug Formation Regime (1 for on)
ct = 0;                               % iteration counter
saveresults = 0;                      % turn on to save results in .csv and .vtk of final timestep
timeplot = 0;                         % if you want to see transient Saturation,WCR,pressure
vdclimit = 2;                        % After how many iterations to plot the surface
gridRes = 50;                         % maximum resolution for grid (slow plotting)

%% Model Initial & Boundary Conditions

% Channel Dimensions
H = 0.001;                             % channel height (m)
W = 0.001;                             % channel width (m)
L = 0.01;                              % channel length (m)
% L = channel_length(repeat)

xmin = 0;                             % min xaxis
xmax = W;                             % max xaxis
ymin = 0;                             % min yaxis
ymax = H;                             % max yaxis

% GDL and Channel Wetting Properties
theta_GDL = 120*(pi()/180);             % GDL contact angle (radians)
theta_wall = 70*(pi()/180);             % wall contact angle (radians)
theta_hyst = 25*(pi()/180);             % contact angle hysteresis (radians)

% Discretise channel points for transport and regime plotting along channel
nz = 10;                                  % number of nodes
dz = L/nz;                                % cell size
Cz = 0:dz:L;                              % location of nodes

% Electrochemical Operating Parameters
% Sto = OPC(repeat,2);                       % Stochiometry of o2 reduction reaction (read from file)
Sto = 2;                                     % Stochiometry
Active_Area = (L*W)*2;                       % Active Area (L*W) + 2 half channel ribs = 2*L*W
Active_Area = 1.5125e-05;                    % Or specify cell cross-section
%current = OPC(repeat,1)*Active_Area;        % Specify current externally
current = 10000*Active_Area;                 % Operating current (must change numerical value)
current_density = current/Active_Area;       % Current density based on active area

% Operating Parameters needed for transport algorithm
xo2 = 0.21;                                 % fraction of O2
xN2 = 1-xo2;                                % fraction of N2
RHi = 1;                                    % Inlet RH
RHd = 1;                                    % Domain RH
Tin = 70+273;                               % Inlet Air Temperature (K)
Tdin = 70+273;                              % droplet temperature
Tdomain = 70+273;                           % Temperature of Domain
Pin = 120000;                               % Inlet Pressure (Pa)
T = ones(size(Cz)).*Tdomain;                % Nodal temperature of domain
Psat = @(T) exp(34.494 - (4924.99./((T-273) +237.1)))./ (((T-273)+105).^1.57); %Saturated vapour pressure equation (Pa)
Ps = Psat(T);                                 % Initialising Psat in domain
Pv = RHd.*Ps;                                 % Initialising vapour pressure in domain
Pvin = RHi.*Psat(Tin);                        % Calculating water vapor pressure in inlet
Cwin = Pvin./(8.314.*Tin);                    % Calculating water concentration inlet (mol/m3)
Cw = (Pv./(8.314.*T));                        % Determine water concentration in domain (mol/m3)
xw = (Pv./Pin);                               % Water mol fraction
u_sto = (Sto/xo2).*((current_density*W*L)./(4*96487)).*((8.314*Tin)/(Pin-Pvin))*(1/(H*W)); %stochiometric air velocity (m/s)
Oin = ((xo2.*(Pin-Pvin))/(8.314*Tin));        % oxygen concentration at inlet
N2 = (0.79/0.21).*Oin;                        % nitrogen concentration
rhow = 1000;                                  % density of water (kg/m3)
rhoa = 1;                                     % density of air   (kg/m3)
sigma = 0.072;                                % surface tension of water (N/m)
w_c = zeros(size(Cz));                        % water in channel segments
Sw_c = w_c;                                   % Saturation of water in channel segments

%% Two-phase Droplet Injection Parameters

% Number of Injection points (i.e. number of discrete water clusters in GDL)

%npores=L*1000;                                                % Set number of pores as function of length
%npores = round(4.6509e-06.*((current_density/10).^3.108).*L) ;% I in mA cm-2 i.e. 1 A cm2 = 1000 mA cm-2 zang 2012 prediction
%npores = round((2.139.*erf(16.95*((current_density/10000)-0.7078)) + 7.631) * (W*1000*L*1000)); % Published Paper correlation
%npores = 730; %i.e. average cluster density of 7 for 100 mm channel
%npores = sens_npores(repeat).*(W*1000*L*1000);
npores=round(7*(W*L*1000000));                              % set number of injection pores (7 clusters/mm2)

%=============== Specify the pore size and location =============================

Rpmax = (50e-6)/2;                                % maximum pore radius                                
%Rpmax = 100e-6;                                  % maximum pore radius (VoF compare)
xp = rand(npores,1).*(W-(2.*Rpmax)) + Rpmax;      % x coord of pores ensures no pore touches wall due to Rpmax
%xp = rand(npores,1)*(W);                         % x coord of pores (random)
%xp = ones(npores,1)*(W/2);                       % x coord of pores (ordered)
yp = zeros(npores,1);                             % y coord of pores (set to zero)
zp = rand(npores,1)*(L/1.0005);                   % z coord of pores (random along channel)
%zp = ones(npores,1)*(L/2);                       % z coord of pores (ordered)
%set_pore_points                                  % externally set pore locations from .csv
%VoF_Straight_Pores_2                             % Sets pore data for VOF comparison simulation

%======================= Preallocate of droplet data ============================

poresize = ones(npores,1)*Rpmax;      % Assign Rpmax to all the pores
Vd = zeros(size(poresize));           % volume of droplets
Rc = poresize;                        % contact radius
Rd = poresize;                        % radius
h = ones(size(poresize)).*1e-9;       % height
b = zeros(size(poresize));            % medial axis
Lc = zeros(size(poresize));           % characteristic length (wall transition)
Ld = zeros(size(poresize));           % length droplet in z direction
E = zeros(size(poresize));            % Equation number
v_evap = zeros(size(poresize));       % evaporation volume
A_surf = zeros(size(poresize));       % area covering GDL surface
A_sl = zeros(size(poresize));         % solid liquid area coverage
flow_l = zeros(size(poresize));       % flow rate into droplet
u_d = zeros(size(poresize));          % velocity of droplets
u_i = zeros(size(poresize));          % velocity of air above droplets
u_ac = zeros(size(poresize));         % velocity of air in channel section
As = zeros(size(poresize));           % surface area of droplets
Ac = zeros(size(poresize));           % cross section area droplet
Lad = zeros(size(poresize));          % length of adhesion contact line
Fdrag = zeros(size(poresize));        % drag force
Fadhesion = zeros(size(poresize));    % adhesion force
theta_app = zeros(size(poresize));    % apparent contact angle
Rp = ones(size(poresize)).*poresize;  % radius of pore
hlim = ones(size(poresize));          % height limit of slug
Td = ones(size(poresize)).*Tdin;

% Determine initial adhesion force based on advancing receding contact angles

adv = theta_GDL;                      % advancing contact angle
rec = theta_GDL-theta_hyst;           % receding contact angle
Lad = 2.*Rp;                          % contact line normal to flow (adhesion force length)
Fadhesion = (24/(pi()^3))*sigma*2.*Lad.*(cos(rec)-cos(adv));  % adhesion force (N)

% ================= Water Flow rates by current density =================================

%mass_flow_water = 1e-6;% * ;                              %Set water injection rate by mass flow rate (kg/s)                    
mass_flow_water = ((current_density*W*L*0.018)/(2*96487)); %Set water injection rate by current density (m3/s)

volume_flow_water = mass_flow_water./rhow;                 %convert mass to volume flow (m3/s)
%volume_flow_water = ((2400/3600)*10^-6)/1000;             %Set water injection rate by volume flow rate (m3/s)

volume_flow_water_pore = (volume_flow_water./npores);                % volume flow per pore (m3/s/pore)
velocity_flow_water = volume_flow_water_pore./(pi().*(poresize.^2)); % velocity of water in each pore (m/s)
u_w = velocity_flow_water;                                           % velocity of water in each pore (m/s)                                       
u_wav = u_w(1);

% =============== Water Injection by Velocity =================================

u_wav = 0.1; % uncomment to set water velocity per pore explicit

% Used when flow rate is specified
 volume_flow_water_pore = u_wav.*(pi().*(poresize(1).^2));     % m/s/pore
 volume_flow_water = volume_flow_water_pore.*npores;

% Assign a flow rate to each pore
flow_l = ones(size(poresize)).*volume_flow_water_pore; %m3/s/pore

% ================ Preallocate Droplet storage Array ==========================
% with original pore values saved in Dorig.

%     1   2       3  4 5  6 7  8  9 10 11 12 13 14 15 16   17   18    19    20  
Dorig = [Vd flow_l xp yp zp E Rd Rc Rp h b Lc Ld As Ac Lad A_surf A_sl v_evap u_d u_i u_ac Fdrag Fadhesion theta_app hlim Td];
%21 u_i  %22 u_ac %23fdrag %24 fadhesion 25 theta_app 26 hlim 27 Td
D = Dorig; %save pore injection conditions in Dorig

%% Fuel Cell Operation Parameters

% Fuel Cell Dimensions
GDLy = 200e-6;          % Thickness of GDL (unused)
MPLy = 10e-6;           % Thickness of MPL (unused)
CLy = 10e-6;            % Thickness of CL  (unused)

% Physical Properties of cell
Do2 = 2.19e-5;                                    % Diffusion coefficeint of ozygen in nitrogen (m2/s)
Dh2o = 21.6*10^-6 .* (1 + (0.0071 .* (T-273)));   % diffusion coefficient of water vapour (m2/s)
Cpg = (1.00926*10^3)-(4.0403*10^-2 *(T-273)) + (6.1759*10^-4 *(T-273).^2) - (4.097*10^-7 *(T-273).^3);    %heat capacity of air
Cpw = 4181;                                       % heat capacity water (J K-1 g-1)
rhog = Pin./(287.058*T);                          % density of air (kg/m3)
kg = (1.97*10^-4).*(T.^0.858);                    % thermal conductivity of air (W m-1 K-1)
alphag = kg./(Cpg.*rhog);                         % thermal diffusivity of air (W m-1 K-1)
mug = (1.691*10^-5)+(4.984*10^-8.* (T-273)) -(3.187*10^-11 .* (T-273).^2)+(1.319*10^-14 .* (T-273).^3); % kg m-2 s-1
hlv =((-2.5464.*D(:,27))+2503.6).*1000;           % latent heat of vaporisation/condensation (kJ/kg)
muw = 0.0006539;                                  % viscosity of water (kg m-2 s-1)
Mw = 0.018;                                       % molecular mass kg/mol

% =========== Initialise Species Matix (Transport) =====================
%T,O2,H20

%Tdin=60+278;                                 % Temperature of Inlet air
%Tdomain=Tin;                                 % Make domain same as inlet
Co2=ones(size(Cz)).*Oin;                      % Concentration of O2 (mol/m3)
%Co2=zeros(size(Cz)).*Oin;                    % zeros (mol/m3)
%T=ones(size(Cz)).*Tdomain;                   % Preallocate domain T
P=ones(size(Cz)).*Pin;                        % Preallocate domain P

%Droplet Temperature (Evaporation)
%Td=ones(size(D(:,1))).*Tdin;
x_o2=Co2.*(8.314.*T)./(P);
Saturation=zeros(size(Cz));

% =============  Set Air Inlet Velocity (m/s) ==============================

uin = u_sto;              % air velocity based on stochiometry and current density (m/s)
uin = 10;                 % air velocity based on set value (m/s)

% IF air velocity (uin=0), air must flow induced by growth of droplets.
%uin = (volume_flow_water./(H*W)).*dt; %i.e. injection of droplets must cause air flow

u_a=ones(size(Cz)).*uin;   % assign velocity to channel points (unused)
V_air_in = uin.*(H*W);     % volume flow of air in (unused)
V_air = u_a.*(H*W);        % volume flow of air in each segment (unused)


Location_Evap = D(:,5) > Cz ;       % find where the location is less than the channel points
[~,i_evap] = min(Location_Evap,[],2); % find the first 1 value in each row (i_evap)
Cbulk = Cw(i_evap)';                  % converting rows to column array for droplet evaporation assesment
Tbulk = T(i_evap)';
Nu = 2 + 0.6.*(Cpg(1).*mug(1)/kg(1)).^(1/3).*((2.*D(:,7).*(u_a(1)-D(:,20)).*rhoa) ./(mug(1)) ).^(1/2); % nusselt number
h_conv = Nu.*kg(1)./(2.*D(:,7)); % convective heat transfer coefficient
h_conv(isinf(h_conv)|isnan(h_conv)) = 0;

%% Balance Viscous, Inertial and Capillary Shear Forces (Film formation == 1)

Ca = uin.*mug(1)./sigma;
A = (rhoa.*(uin.^2))./(4.*sigma.*(1-cos(theta_wall)));
%coefficients for equation
a = A;
b = 1 + (3.*Ca.*H./(1-cos(theta_wall))) - (2.*A.*H) ;
c = -(3.*Ca.*(H)./(1-cos(theta_wall))) -(2.*H);
d =  (H.^2);

% Newton iteration to solve
f=@(h)  a.*(h.^3) + b.*(h.^2) + (c.*h) + d;
df=@(h)  (a.*h.^2).*3 + (b.*h).*2 + c;

     h=ones(size(uin)).*(0);
     holde = h; % h;
     h1 = h; %h 
     hnew = holde + 10000;
    
     f1_error=1e-9;
iter=0;
w=0.7;
     while sqrt(abs(transpose(holde-hnew)*(holde-hnew)))>1e-12 %solving the height equation equation (Newton iterations)
         hnew=h1 - (w.*df(h1)).\f(h1);
         holde=h1;
         h1=hnew;
         iter=iter+1;
     end
     
     h1(h1>H) = H; %h cannot be above H
     h1(h1<0) = 0; %h cannot be below 0
  
     hlim = ones(size(poresize)).*h1;    %maximum height before deformation in capillary bridge
     D(:,26) = hlim;
     Dorig(:,26) = hlim; %height limit for each droplet based on the external forces
     
%% ______________ Simulation Time Integration Settings_______________________

dt = 1e-7;                      % timestep (may be unstable for too large flow rates)
dtorig_i = 1;                   % original timestep
%dtorig_i = 0.5*H./uin;         % if air velocity is high (0.5 courant number) %remove for limiting time step
dt3 = 0.5*H./u_wav;             % if water velocity is high
dtorig = min([dtorig_i dt3]);   % finds minimum between the two
dt2 = dtorig;
%dtorig = 1e-3;                 % alternatively, set the max time step

Evap = 0;                         % initialise evaporation
V_evap= 0;                        % initialise evaporation
ic = 0;                           % iteration counter
plotlim = 1;                      % plotting limit (how many timesteps to skip)
vtkploti = 0;
plotlimcount = 0;                 % plotting limit counter
write_count = 0;
av_vol_fraction = 0;
maxdroplets = 10000;              % determine exit condition for divergence of model

% Determine simulatuon time defined by 1 channel volume of water injected:
Channel_Volume = H*W*L;                          % Channel Volume (m3)
endtime = 3.*(Channel_Volume/volume_flow_water); % Simulation will stop at this time (now 3 times longer to be able to get better time average)
%endtime = 100;                                  % alternatively, set endtime
t=0;                                             % initial time
time_value=0;                                    % time counter
nwrite=2;                                        % change number to change write interval
write_precision=1*10^(-nwrite);                
timer_tic = 0;
timer_tic_coal = 0;
courant = 0.1;
visualiseDropletCounter = 0;

%% Start Time Integration 
for tt=0:1:100000000000 %maximum time loop

%%============= Initialise simulation with droplets already present ==============
 if t==0
%     Example of introducing artificial droplets
%     Di=Dorig(1,:)
%     Di(:,1) = 5.98e-10
%     Di(:,3) = W/2
%     Di(:,4) = 2*H
%     Di(:,5) = L/2
%     Di(:,7) = 0.001
%     Di(:,6) = 5
%     %Di(:,24) = 
%     Di(:,20) = 0
%     D=[D;Di];
 % D(:,1)=5.98e-11; %m3
 end
 
 %%=========== Adaptive timestep ==============%%
 
 adapt_t=courant.*(0.001)./(D(:,20));
 adapt_t(isinf(adapt_t))=dtorig;
 dt=min(adapt_t);
 
 if max((D(:,23)./D(:,24)))>0.999 
     dt2=1e-3;
 end
 
 dt2(dt2==0)=dtorig;
 dt3 = 0.5*dz./uin; %new addition for transport
 dt3=0;             %remove if transport not needed
 dt3(dt3==0)=dtorig;
 dt = min([dt dt2 dt3]);
 dt2=0;
 
 %dt=1e-3;
 %dt=1e-4;           %activate to set timestep

 t=t+dt;            %advancment of time
 if t > endtime
     break
 end
tic;         
ic = ic + 1;                     % iteration counter
plotlimcount = plotlimcount + 1; % plotting iteration counter



if transport==1 %solves transport equations if transport ==1
%transport of water vapour explicit
Cwn=Cw;
Tn=T;
Co2n=Co2;
V_airn=V_air;
i=2:1:size(Cz,2)-1;

Cw(:,i)=Cwn(:,i) + dt.*((u_a(:,i-1).*Cwn(:,i-1) - u_a(:,i).*Cwn(:,i))./dz)...
                 + dt.*(Dh2o(:,i).*((Cwn(:,i-1)- 2.*Cwn(:,i) + Cwn(:,i+1))./dz.^2));
                 
             
% %transport of temperature             
% T(:,i) = Tn(:,i) + dt.*((u_a(:,i-1).*Tn(:,i-1) - u_a(:,i).*Tn(:,i))./dz) ...
%                  + alphag(:,i).*((Tn(:,i-1) - 2*Tn(:,i) + Tn(:,i+1))./(dz.^2)).*dt;
%                  
                 
% %transport of oxygen
% 
 Co2(:,i)=Co2n(:,i) + dt.*((u_a(:,i-1).*Co2n(:,i-1) - u_a(:,i).*Co2n(:,i))./dz)...
                    + dt.*(Do2.*((Co2n(:,i-1)- 2.*Co2n(:,i) + Co2n(:,i+1))./dz.^2))...
                    - dt.*(((current_density.*W.*dz)/(4*96487))./(H*W*dz));

% droplet temperature update D(:,27)=Td


    
%Boundary Condition Update
%temperature
T(:,1) = Tin;
%T(:,1) = Tn(:,1) + dt.*((uin.*Tin - u_a(:,1).*Tn(:,1))./dz)...
%        - alphag(:,1).*((Tin - 2*Tn(:,1) + Tn(:,2))./(dz.^2)).*dt;

%T(:,end) = Tn(:,end) + dt.*((u_a(:,end-1).*Tn(:,end-1) - u_a(:,end).*Tn(:,end))./dz) ...
                % - alphag(:,end).*((Tn(:,end-1) - 2*Tn(:,end) + Tn(:,end))./(dz.^2)).*dt; ...
T(:,end)= T(:,end-1);
%water vapour          
Cw(:,1)=Cwin;
                 %- dt.*(Dh2o(:,1).*((Cwin - 2.*Cwn(:,1) +
                 %Cwn(:,2))./dz.^2));2
                 
Cw(:,end)=Cwn(:,end-1);
%                               
% Cw(:,end)=Cwn(:,end) + dt.*((u_a(:,end-1).*Cwn(:,end-1) - u_a(:,end).*Cwn(:,end))./dz)...
%                  - dt.*(Dh2o(:,end).*((Cwn(:,end-1) - 2.*Cwn(:,end) + Cwn(:,end))./dz.^2));
%              
%oxygen concentration
%Co2n(:,1) + dt.*((uin.*Oin - u_a(:,1).*Co2n(:,1))./dz)...
%                    - dt.*(Do2.*((Oin - 2.*Co2n(:,1) + Co2n(:,2))./dz.^2))...
%                    - dt.*(((current_density.*W.*dz)/(4*96487))./(H*W*dz))
Co2(:,1) = Oin;              
% Co2(:,1)=Co2n(:,1) + dt.*((uin.*Oin - u_a(:,1).*Co2n(:,1))./dz)...
%                    - dt.*(Do2.*((Oin - 2.*Co2n(:,1) + Co2n(:,2))./dz.^2))...
%                    - dt.*(((current_density.*W.*dz)/(4*96487))./(H*W*dz));
Co2(:,end)=Co2n(:,end-1);             
%  Co2(:,end)=Co2n(:,end) + dt.*((u_a(:,end-1).*Co2n(:,end-1) - u_a(:,end).*Co2n(:,end))./dz)...
%                         - dt.*(Do2.*((Co2n(:,end-1) - 2.*Co2n(:,end) + Co2n(:,end))./dz.^2))...
%                         - dt.*(((current_density.*W.*dz)/(4*96487))./(H*W*dz));
%

%volume balance
 i=2:1:size(Cz,2)-1;
 V_air(:,1) = V_air(:,1) + V_air_in - V_air(:,2)...
             - dt.*(((current_density.*W.*dz)/(4*96487)).*(8.314*T(:,1)./P(:,1)));
         
 V_air(:,i) = V_air(:,i) + V_air(:,i-1) - V_air(:,i+1)...
             - dt.*(((current_density.*W.*dz)/(4*96487)).*(8.314*T(:,i)./P(:,i)));
             
 V_air(:,end) = V_air(:,end) + V_air(:,end-1) - V_air(:,end)...
             - dt.*(((current_density.*W.*dz)/(4*96487)).*(8.314*T(:,end)./P(:,end)));
 
 uaprox=V_air./(H*W);
% 
%Ps=Psat(T);
RH=(Cw.*8.314.*T)./Ps;

end

%% Droplet Transport and update

% Search all injection points to see if droplet is being injected into.
% amend july 2021 (only search those points that are within the search radius, i.e create new search radius function
search_radius = 4.*D(:,7); %i.e. search radius is 4 times the droplet radius

timer_tic = timer_tic+1;
timer_tic = 1;
if timer_tic == 1
    
if size(D(:,1),1)>1
    
    for i=1:1:size(D,1)
        
        %improve speed algorithm, 
    %distance_compare= D(i,7) > sqrt((D(i,5) - zp(:,1)).^2 + (D(i,3) - xp(:,1)).^2 + (D(i,4) - yp(:,1)).^2 );
              %droplet radius + pore radius > distance from droplet center to injection point center (for all injection points)
   % distance_compare= D(i,7)+(poresize(zp<D(i,5)+0.001 & zp>D(i,5)-0.001)) > sqrt((D(i,5) - zp(zp<D(i,5)+0.001 & zp>D(i,5)-0.001 ,1)).^2 + (D(i,3) - xp(zp<D(i,5)+0.001 & zp>D(i,5)-0.001 ,1)).^2 + (D(i,4) - yp(zp<D(i,5)+0.001 & zp>D(i,5)-0.001 ,1)).^2 );
    distance_compare= D(i,7)+(poresize) > sqrt((D(i,5) - zp(:,1)).^2 + (D(i,3) - xp(:,1)).^2 + (D(i,4) - yp(: ,1)).^2 );
    D(i,2)=sum(distance_compare).*volume_flow_water_pore;
    
    if D(i,6)==1
        D(i,2)=volume_flow_water_pore;
    end
    
    if D(i,6)==0
        D(i,2)=volume_flow_water_pore;
    end
    end
   
    
end

 timer_tic=-1;
 
 
 
end
    
%% ============== injection volume balance ============
%D(D(:,1)==0,19)=0;

D(:,1) = D(:,1) + ((D(:,2) - D(:,19)).*dt);
%D(:,1) = D(:,1) + ((D(:,2)).*dt)- D(:,19);

if evaporation==1
%temperature update
D(:,27) = D(:,27) + (h_conv.*D(:,14).*(Tbulk - D(:,27))./(rhow.*D(:,1).*Cpw)).*dt...
                  - (hlv.*(D(:,19).*rhow)./(rhow.*D(:,1).*Cpw)).*dt;
end

% if D(:,1)<0
%     break
% end
%removal of droplets with D(:,1)<0
D=D(D(:,1)>0,:); %new addition 17/06/2020, may help for evaporation

D(D(:,1)<0,1)=0;

%D(:,23)=D(:,23)*0; % remove this section for code to work 

%%===========  momemtum balance (Newton second law) ===================%%

dvdt=subplus((D(:,23)-D(:,24))./(D(:,1).*rhow));            % acceleration update
% move=D(:,20)>0;
% dvdt(move)=(D(move,23)-D(move,24))./(D(move,1).*rhow);    % acceleration update

D(:,20)=D(:,20) + (dvdt.*dt);                       % velocity update
D(D(:,20)<0,20)=0;                                  % ensure no negative velocity
%D((D(:,6)<6) & (D(:,10)>H),20)=uin;                % new addition 08/06/2020
%D((D(:,6)<6) & ((D(:,4) - D(:,7))>H),20)=uin; 
%D((D(:,6)==6) & ((D(:,4) - D(:,7))>H),20)=uin;
D((D(:,6)==7),20)=uin;
D(D(:,20)>uin,20)=uin;                              % new addition 08/06/2020
D(:,5)=D(:,5) + (D(:,20).*dt);                      % z position update
 
%%============ After water injection, update classificiations ==================

        % class 0 - ghost
        %D(D(:,6)==0,7)=Rpmax;
        
        % class 1 - emerging
        D((D(:,1)>0) & (D(:,6)==0),6)=1;                 %finds all droplets with volume and assigns as emerging (1)
        D((D(:,6)==2) & (D(:,25) < theta_GDL),6)=1;
        %D(D(:,6)==1,7)=Rpmax;
        
        % class 2 - isolated
        D((D(:,6)==1) & (D(:,25) > theta_GDL),6)=2; %finds all emerging droplets, and assigns as isolated if constant contact angle regime (2)
        
        % class 3 - side wall (west)
        r3find_minus=find((D(:,6)<3) & ( (D(:,3) - D(:,7) ) < 0));  %finds droplets in class less than 3 and that have a radius that touches left wall
        D(r3find_minus,6)=3;          %assign classification 3
        D(r3find_minus,3)=0;          %x position is left wall
        D(r3find_minus,4)=0;          %y position is zero
        
        % class 3 - side wall (east)
        r3find_plus=find((D(:,6)<3) & ( (D(:,3) + D(:,7) ) > W));  %finds droplets in class less than 3 and that have a radius that touches right wall
        D(r3find_plus,6)=3;          %assign classification 3
        D(r3find_plus,3)=W;          %x position is left wall
        D(r3find_plus,4)=0;          %y position is zero
        
        
        % class 4 - corner droplet (north)
        r4find=find( ( D(:,6)<4 ) & (D(:,10)>H) ); %find droplets class less than 4, and also h greater than the height of channel
        D(r4find,6)=4;
        D(r4find,4)=H;
        D(r4find,2)=0; %new 05/06/2020
        
        
        % class 5 - truncated capillary bridge
        r5find=find( ( D(:,6) <5 ) & ( D(:,12)>W ) ); %finds droplets class less than 5 and where Lc (length across width is greater than the channel width
        D(r5find,6)=5;
        
        
if film_formation==1
        r6find=find((D(:,6) == 5) & D(:,10)>D(:,26));
        D(r6find,6)=6;
        
end

if plug_formation==1
        %class 6 plug flow
        r7find = find((D(:,6)==6) & D(:,10)>(0.99*H));
        D(r7find,6)=7;
%         
        r5find2=find( ( D(:,6)==7 ) & ( D(:,10)<H ) );
        D(r5find2,6)=5;
end 
   
   water_removed=find((D(:,5) - D(:,13))>L);
   water_cond = sum(D(water_removed,1)); %new addition, track water removed
   newmatrix=find((D(:,5) - D(:,13))<L);  % new 16/06/2020
   D = D(newmatrix,:);
   timetest(ic,1)=toc;                      %save time for section 1
   
   
 %% ============= Coalescence Algorithm ================= %%
 
 tic
 if coalescence==1
     
     
timer_tic_coal = timer_tic_coal+1;
timer_tic_coal = 1;
if timer_tic_coal == 1
     
 if size(D(:,1),1)>1
 
 % each droplet tracks all other droplets, build upper triangular matrix, not to self duplicate connection
 G = ones(size(D,1));
 % remove selfconnection from upper tri matrix and make sparse matrix
 K = sparse(triu(G,1));
 % sparse matrix saves memory because 0 are not stored (i.e. half matrix)
 [i,j,~] = find(K);         % access index of sparse matrix
 ij = [i j];                % create connectivity matrix
 k = 1:1:size(ij,1);        % for all potential connections
 % determine distance between centers of all droplets
 s = sqrt((D(ij(k,1),3) - D(ij(k,2),3)).^2 + (D(ij(k,1),4) - D(ij(k,2),4)).^2 + (D(ij(k,1),5) - D(ij(k,2),5)).^2);
 % determine sum of radii
 DR = D(ij(k,1),7) + D(ij(k,2),7);
 % special case for film flow because R is very large
 film_index = find(D(ij(k,1),6)==6 & D(ij(k,1),6)==6);
 DR(film_index) = (D(ij(film_index,1),13)./2) + (D(ij(film_index,2),13)./2);
 
 %mark pair for coalescence
 coal = find(DR > s);  %find row coordination for coalescing pair
 d1 = ij(coal,1);      %coalescing pair d1
 d2 = ij(coal,2);      %coalescing par d2 
 ik=0;
 %for all colliding droplets
 for ik=1:1:size(d1,1)
   
     [imax,jmax] = max([D(d1(ik),1) D(d2(ik),1)]); %find which is largest droplet
     dvolumes=[D(d1(ik),3) D(d2(ik),3);D(d1(ik),4) D(d2(ik),4);D(d1(ik),5) D(d2(ik),5)]; %find largest droplet
     %new droplet is equal to droplet 1
     Dn(ik,:) =  D(d1(ik),:);
     Dn(ik,1) = D(d1(ik),1) + D(d2(ik),1); %volume balance
     %Dn(ik,2) = 0;                 %flow rate into new droplet
     %Dn(ik,3) = dvolumes(1,jmax);  %x position new
     
     Dn(ik,3) = (D(d1(ik),1).*rhow.*D(d1(ik),3) + D(d2(ik),1).*rhow.*D(d2(ik),3))./(Dn(ik,1).*rhow);
     Dn(ik,4) = (D(d1(ik),1).*rhow.*D(d1(ik),4) + D(d2(ik),1).*rhow.*D(d2(ik),4))./(Dn(ik,1).*rhow);
     Dn(ik,5) = (D(d1(ik),1).*rhow.*D(d1(ik),5) + D(d2(ik),1).*rhow.*D(d2(ik),5))./(Dn(ik,1).*rhow);
     
     %Dn(ik,4) = dvolumes(2,jmax);  %y position new
     %Dn(ik,5) = dvolumes(3,jmax);  %z position new
     Dn(ik,6) = max([D(d1(ik),6) D(d2(ik),6)]); %Classifciation is equal to the largest between
     Dn(ik,7) = max([D(d1(ik),7) D(d2(ik),7)]);
     Dn(ik,8) =  max([D(d1(ik),8) D(d2(ik),8)]); %changed 16/06/2020
     Dn(ik,9) = max(poresize);
     Dn(ik,20)=(D(d1(ik),1).*rhow.*D(d1(ik),20) + D(d2(ik),1).*rhow.*D(d2(ik),20))./(Dn(ik,1).*rhow); 
     Dn(ik,11)=(H/2);
     Dn(ik,27)=max([D(d1(ik),27) D(d2(ik),27)]);
     
     
     
     %kill connecting droplets and merge new droplet
     D(d1(ik),:) = zeros(1,size(D,2));      %droplet 1 == 0
     D(d2(ik),:) = zeros(1,size(D,2));      %droplet 2 == 0;
     D(isnan(D))=0;                     %make any nan =0;
     D = [D;Dn(ik,:)];
     %D = D(any(D,2),:); 
     
 
 end
  
  
 if isempty(ik)==0
 %Dnew = [D;Dn];
 %D = Dnew
 clear Dn
 clear dvolumes
 clear coal
 clear d1 d2
 clear s G K ij k
 D = D(any(D,2),:); 
 end
 D((D(:,3)==0 & D(:,4)==0 & D(:,5)==0),:) = D((D(:,3)==0 & D(:,4)==0 & D(:,5)==0),:).*0; %added 26/06/2020 fixes coalescence problem
 D(isnan(D))=0;      
 D = D(any(D,2),:); 
 end
 timer_tic_coal=-1;
end  
 end

  timetest(ic,2)=toc;

  %% emergence algorithm

  tic 
  i_e=0;
  
  for i_p=1:1:size(zp,1)
      
      s2 = sqrt((xp(i_p) - D(:,3)).^2 + (yp(i_p) - D(:,4)).^2 + (zp(i_p) - D(:,5)).^2);
      dr2 = D(:,7) + (poresize(i_p));
      
      
     s2(D(:,6)==1)=0; %new addition 9th june
     s2(D(:,6)==0)=0; %new addition 9th june
      
      pore_cover = s2 < dr2;
      if sum(pore_cover)<1
          i_e=i_e+1;
          D_emerge(i_e,:) = Dorig(i_p,:); 
      end  
  end
 
  if i_e>0
        size(D_emerge,1);
        D = [D;D_emerge];
  end
  
 if size(D,1) > maxdroplets
     disp('Model has diverged, lower time step')
     break 
 end
    D(isnan(D))=0;
    D = D(any(D,2),:); 
    clear D_emerge pore_cover dr2 s2
 timetest(ic,3)=toc;
 
 D(D(:,27) < 273,27)=Dorig(1,27); %added for temperature phase change
 
 
  %% Create Droplet Classes
  tic
  D0 = find(D(:,6)==0); %Find all droplets that use equation 0
  D1 = find(D(:,6)==1); %Find all droplets that use stage 1
  D2 = find(D(:,6)==2); %Find all droplets that use stage 2
  D3 = find(D(:,6)==3); %Find all droplets that use stage 3
  D4 = find(D(:,6)==4); %Final all drops using stage 4
  D5 = find(D(:,6)==5); %Final all drops using stage 5
  D6 = find(D(:,6)==6); %Find all drops using stage 6
  D7 = find(D(:,6)==7); %Find all drops using stage 7
  
  %separate droplets based on class
  F0 = D(D0,:);   
  F1 = D(D1,:);
  F2 = D(D2,:);
  F3 = D(D3,:);
  F4 = D(D4,:);
  F5 = D(D5,:);
  F6 = D(D6,:);
  F7 = D(D7,:);
  
  %% Update Droplet Geometric Parameters
  %F0(:,17)=pi().*(Rpmax).^2;         %surface area GDL coverage
  %% Emerging (Constant Contact Radius) (CCR)
  if isempty(D1)==0
     
  %Equation V=pi()h/6 (3 rc^2 + h^2) ... We know rc and V -> h and theta
  %height (:,10)
  %radius of curvature (:,8)
  %radius of pore (:,9)
  %theta_app (:,25)
  %droplet radius (:,7) used for coalescence
  
     holde = F1(:,10); % h;
     h1 = F1(:,10); %h 
     hnew = holde + 10000;
    
     f1_error=1e-9;
     F1_1=     @(rpore,h,Vd)    h.^3 + 3.*(rpore.^2).*h - 6*Vd/pi();    %Equation for stage 1
     F1_J=     @(rpore,h)       3.*h.^2 + 3.*(rpore.^2);                %jacobian for stage 1
     F1_R=     @(Vd,h)          (((3*Vd)./(pi()*(h.^2)))+h)/3;           %radius as a function of height
     F1_Ac=    @(R,theta_1c)    ((R.^2)/2).*(theta_1c-sin(theta_1c));    %cross section area
     F1_As=    @(R,h)           2.*pi().*R.*h;                          %surface area
     F1_Ac_90= @(R,theta_app)   (pi()*R.^2)-((R.^2)/2).*((2.*(pi()-theta_app))-sin(2.*(pi()-theta_app)));
     
     while sqrt(abs(transpose(holde-hnew)*(holde-hnew)))>1e-9 %solving the height equation equation (Newton iterations)
         hnew=h1 - F1_J(F1(:,9),h1).*eye(length(D1))\F1_1(F1(:,9),h1,F1(:,1));
         holde=h1;
         h1=hnew;
     end
  
     F1(:,10)=hnew;                     %update height
     F1(:,8)=F1_R(F1(:,1),F1(:,10));    %update radius of curvature
     
     b90=find(F1(:,10) - F1(:,8) < 0); %for all droplets where the radius is greater than the height o<90
     a90=find(F1(:,10) - F1(:,8) > 0); %for all droplets where the radius is greater than the height o>90
     
     F1(b90,25)=real(asin(F1(b90,9)./F1(b90,8)));      %apparent contact angle      
     F1(b90,15) = F1_Ac( F1(b90,8) , 2.*F1(b90,25) ); %cross section area
     F1(b90,14) = F1_As( F1(b90,8)   , F1(b90,10)  ); %surface area
     F1(b90,7) = F1(b90,9);                           %radius of curvature equal to rpore for plotting
     
     F1(a90,25)=pi()-real(asin(F1(a90,9)./F1(a90,8))); %apparent contact angle      
     F1(a90,15) = F1_Ac_90( F1(a90,8) , F1(a90,25) ); %cross section area
     F1(a90,14) = F1_As( F1(a90,8)   , F1(a90,10)  ); %surface area
     F1(a90,7) = F1(a90,8);                           %radius of curvature used for plotting
     
     r2find=find(F1(:,25)>theta_GDL);
     F1(r2find,6)=2;
     
     
     %F1(:,4) = F1(:,10) - F1(:,7);        %y position update might need to change radius
     F1(:,4) = F1(:,10) - F1(:,8);        %y position update might need to
     %change radius 11/06/2020
     F1(:,11)= (H - F1(:,10))/2;          %medial axis y location update (b)
     F1(:,13)= 2.*F1(:,7);                %length of droplet in z direction
     F1(:,17)=pi().*(F1(:,7)).^2;         %surface area GDL coverage
     F1(:,16)=F1(:,9);                    %adhesion length
     F1(:,18)=(pi().* F1(:,9).^2);         %liquid solid shear area
     
  end
  
  %% Isolated
  if isempty(D2)==0
  
%Analytical Geometry equation for R as function of V 
%(V = (4/3)*pi()*R^3*(2+cos(theta)(1-cos(theta))^2
%Finding Radius of curvature
          
%new radius of curvature found by analytical geometric equation 2
Rnew=(F2(:,1)./((4/3)*pi() - ((pi()/3)*((2+cos(pi()-theta_GDL))*((1-cos(pi()-theta_GDL)).^2))))).^(1/3); %Finding Radius of curvature      
                        
                                    F2(:,8)=Rnew;                                 %Radius of Curvature
                                    F2(:,9)=F2(:,8)*sin(pi()-theta_GDL);             %Rpore (also known as Rcontact)
                                    F2(:,10)=F2(:,8)*(1-cos(theta_GDL));             %Height of droplet
                                    F2(:,15)=F1_Ac_90(F2(:,8),F2(:,25));             %Cross sectional area of droplet
                                    F2(:,14)=F1_As(F2(:,8),F2(:,10));                %surface area
                                    F2(:,7)=F2(:,8);                                 %R for plotting and coalescence
                                    F2(:,4)=F2(:,10)-F2(:,7);                        %y position
                                    F2(:,11)=(H-F2(:,10))/2;                         %medial axis b
                                    F2(:,13)=2*F2(:,7);                              %length droplet along channel z
                                    F2(:,17)=(F2(:,9)).^2*pi();                      %area coverage of droplet (need to update this)
                                    F2(:,16)=F2(:,9);                                %adhesion length
                                    F2(:,18)=(pi().* F2(:,9).^2);                     %liquid solid shear area
                                    if theta_GDL< (pi()/2)
                                     F2(:,7)=F2(:,9);   
                                    end
                       
      
  end
  clear Rnew
  
  %% Side Wall
  if isempty(D3)==0
      
      alpha_wall = (pi()/2) - theta_wall;
      alpha_GDL = theta_GDL - (pi()/2);
      theta_c = 2*(pi()-theta_GDL);
               
      if theta_GDL<(pi()/2)
          theta_c=theta_GDL;
      end
               
      theta_c_total = theta_wall + alpha_GDL;
      a3 = 2 - (3*sin(alpha_wall)) + (sin(alpha_wall)^3);
      %a3 = 2 - (3*sin(theta_wall)) + (sin(theta_wall)^3);
      %b3 = 1 - (theta_c-sin(theta_c))/2;
      b3 = 1 - (theta_c -  sin(theta_c))/(2*pi());
      %c3 = a3*b3;
      c3 = pi()*a3*b3;
      
      %New Radius of curvature (look into changing this 10/06/2020
      Rnew3=(3.*F3(:,1)./(c3)).^(1/3);
               
      x31 = Rnew3.*sin(alpha_wall);
      x32 = Rnew3-x31;
      x33 = Rnew3.*sin(alpha_GDL);
      x34 = (Rnew3.^2 - x31.^2).^(1/2);
      x35 = (((x33./tan(alpha_GDL))-x31)./(1-cos(theta_wall)));
      
      F3_Ac = @(R,theta_c_total,x5,x3,x4) ((R.^2)/2).*(theta_c_total-sin(theta_c_total)) + (x5.*(x3+x4)/2); %cross section area
      F3_As = @(R,theta_c,x2) 2.*pi().*R.*x2.*(1-(theta_c-sin(theta_c))./2);
      
      F3(:,8) = Rnew3;                                                       %Radius of Curvature
      F3(:,9) = x34;                                                         %Rpore (also known as Rcontact)
      F3(:,16) = x34 + x33 + (F3(:,8)-  x31);                                %adhesion length
      F3(:,10) = x33 + x34;                                                  %Height of droplet
      F3(:,15) = F3_Ac(Rnew3,theta_c_total,x35,x33,x34);                     %Cross sectional area of droplet
      F3(:,14) = F3_As(Rnew3,theta_c,x32);                                   %surface area
      F3(:,7) =  F3(:,8);                                                    %R for coalescence and plotting
      F3(:,11) = (W - x32)/2;                                                %medial axis b
      F3(:,13) = F3(:,8);                                                    %length along channel z % need to update
      %added the ratio of contact angles to work out bottom length
      F3(:,17) = ((((x33./tan(alpha_GDL))-x31)./(1-cos(theta_wall))).^2/2)*(2*theta_wall - sin(2*theta_wall));
      %F3(:,17) = ((x32).^2)/2 * ((2*theta_wall)-sin(2*theta_wall));
      F3(:,18) = pi().*(x34.^2) ;                                             %liquid solid shear area (update)
      
      %Determine x and y location of droplet
      x3min = find((F3(:,3) < (W/2)));                %for droplets at xmin wall
      x3max = find((F3(:,3) > (W/2)));                %for droplets at xmax wall
      F3(x3min,3) = xmin - x31(x3min);                %update x position
      
      F3(x3max,3) = xmax + x31(x3max);                %update x position
      y3min = find((F3(:,4)==ymin));                %update y position
      F3(y3min,4) = ymin + x33(y3min);                %update y position
      F3(:,4) = F3(:,10) - x34;                        %update y position
         
  end
  
  
  %% Corner
  if isempty(D4)==0
      
        alpha_4 = (pi()/2) - theta_wall;
        alpha_42 = theta_wall - (pi()/4); %used to calculate cross section
        a4 = 2 - (3*sin(alpha_4)) + (sin(alpha_4)^3);
        b4 = (2*theta_wall - sin(2*theta_wall));
        c4 = 2.*b4;
         
        Rnew4 = ((6*F4(:,1)./(a4*b4))).^(1/3);
        x41 = Rnew4*sin(alpha_4);
        x42 = Rnew4 - x41;
        F4_As = @(R,alpha,theta_wall) 2*pi().*(R.^2) * (1 - sin(alpha))*(2*theta_wall - sin(2*theta_wall)); %surface area
        %need to fix this cross section equation 19/05/2020
        %F4_Ac = @(R,alpha,theta_wall,x42) ((x42.^2)/2); %+ ((R.^2)/2)*(2*(pi() - alpha) - sin(2*(pi() - alpha))); 
        F4_Ac = @(R,x42) ((x42.^2)/2) + ((R.^2)/2)*((2*alpha_42) - sin(2*alpha_42)); %fixed 11/06/2020

        %update geometry based on curvature                 
        F4(:,8) = Rnew4;                                               %Radius of Curvature
        F4(:,9) = F4(:,8)*(1-sin(alpha_4));                            %Rpore (also known as Rcontact) 
        F4(:,15) = F4_Ac(F4(:,8),x42);                                 %Cross section area
        F4(:,14) = F4_As(F4(:,8),alpha_4,theta_wall);                  %Surface area
        F4(:,7) = Rnew4;                                               %R for coal and plotting
        F4(:,12) = x42;                                                %length in x direction Lc
        F4(:,11) = (W - x42)/2;                                        %medial axis b (updated 7/04/2020)
        F4(:,13) = 2*Rnew4*sin(theta_wall);                            %length in channel z direction
        F4(:,17) = 0;                                                  %no GDL surface coverage for corner droplets
        F4(:,4) = ymax + x41;                                          %new y position
        F4(:,16) = x42.*2;                                             %adhesion length
        F4(:,18) = 2.*pi().*(F4(:,8).^2).*(theta_wall/pi());               %may need to change
        
        x4min = find(F4(:,3) < (W/2));
        x4max = find(F4(:,3) > (W/2));
        F4(x4min,3) = xmin - x41(x4min);
        F4(x4max,3) = xmax + x41(x4max);
   
      
  end
 
  %% Capillary bridge
  if isempty(D5)==0
      
    clear z51
  
    V5 = F5(:,1);
    z51 = ((W/2)./cos(theta_wall)) * (1-cos((pi()/2) - (theta_wall)));
    R51 = (W/2) * sec(theta_wall);
    theta_51 = (pi() - (2*theta_wall));
    theta_52 = theta_wall*2;
    
    phi_51 = (R51.^2).*(theta_51 - sin(theta_51))./(2*z51*W);
    phi_52 = (theta_52 - sin(theta_52))./(2*pi());
    
    a5 = phi_52*W*pi();
    b5 = -2*phi_52*W*pi()*phi_51*z51;
    c5 = phi_52*W*pi()*phi_51*(z51.^2) - F5(:,1);
    
    R52 = (b5 + sqrt(b5.^2 - (4*a5.*c5)))./(2*a5);
    R53 = R52-z51;
    h5 = R52.*(1 - cos(theta_wall));
    h53 = R53.*(1 - cos(theta_wall));
    %new section 11/06/2020
    h5b = (h53 - H);
    h5b(h5b<0)=0;
    c5b = 2.*R52.*sin(acos(1 - (h5b./R52)));
    c5b(c5b<0)=0;
   
      F5(:,8) = R52;                                % Radius of Curvature
      F5(:,7) = R52;                                % R for plotting/coalescence
      F5(:,13) = 2.*R52.*sin(theta_wall);           % length along channel LD
      F5(:,15)= R53.*(1 - cos(theta_wall)).*W;      % Crosssectionarea
      F5(:,15)= h5.*W;                             %may affect 25/06/2020
      F5(F5(:,15) > H*W,15) = H*W;
      F5(:,14) = pi().*(R52).*phi_52*W;             % surface area
      F5(:,11) = (H - h5)./2;                       % medial axis b (need to check)
      F5(:,17) = 0;                                 % GDL surface Coverag
      F5(:,17) = c5b.*W;                            % new estimation 11/06/2020
      F5(:,3) = (W/2);                              %x position
      F5(:,4) = H + (R52 - h5);                     %y positionneed to change this 27/01/2020
      F5(:,9) = W + 2*(R52*(1-cos(theta_wall)));    %Rpore
      %F5(:,16) = F5(:,9);                           %LAD
      F5(:,16) = W;
      F5(:,18) = F5(:,13).*W;                       %liquid solid area
      F5(:,12) = h5;
      %new additions 16/06/2020
      F5(:,10) = h53;                    %new droplet height
     
  end
  clear h5 
  %% plug %empty for now
  if isempty(D6)==0
   %%New Addition 16/07/2020 Inertial,viscous,capillary balance relvelocity
   rel_v6 = uin - F6(:,20);
   
   Ca6 = rel_v6.*mug(1)./sigma;
   A6 = (rhoa.*(rel_v6.^2))./(4.*sigma.*(1-cos(theta_wall)));
   %coefficients for equation
   a6 = A6;
   b6 = 1 + (3.*Ca6.*H./(1-cos(theta_wall))) - (2.*A6.*H) ;
   c6 = -(3.*Ca6.*(H)./(1-cos(theta_wall))) -(2.*H);
   d6 =  (H.^2);

   f6=@(h6)  a6.*(h6.^3) + b6.*(h6.^2) + (c6.*h6) + d6;
   df6=@(h6)  (a6.*h6.^2).*3 + (b6.*h6).*2 + c6;

   h6=ones(size(uin)).*(0);

     holde6 = h6; % h;
     h16 = h6; %h 
     hnew6 = holde6 + 10000;
    
     f1_error6=1e-9;
     iter6=0;
     w6=0.7;

     while sqrt(abs(transpose(holde6-hnew6)*(holde6-hnew6)))>1e-12 %solving the height equation equation (Newton iterations)
         hnew6 = h16 - (w6.*df6(h16)).\f6(h16);
         holde6 = h16;
         h16 = hnew6;
         iter6 = iter6 + 1;
    
     end
    
     h16(h16>H) = H; %h cannot be above H
     h16(h16<0) = 0; %h cannot be below 0
  
     F6(:,26) = h16;    %maximum height before deformation in capillary bridge hlim
      
      %% start of geometry section
      
      
    F6(:,15)=(H*W);
    h6 = F6(:,10);
    h6 = F6(:,26);
    theta_61 = pi()-(2.*(theta_wall));
    theta_62 = 2.*(theta_wall - (pi()/4));
    R61=(W/2).*sec(theta_61);
    x61 = tan(pi()/4).*h6;
    x62 = x61./tan(theta_62);
    %R62 = x62 + h6;\
    R62 = F6(:,8);
    y6 = H + x62;
    
    L6 = (F6(:,1) - (2.*h6.*x61.*W) - ((R62.^2).*(theta_62 - sin(theta_62)).*W))./ ...
         ((h6.*W) -((R61.^2)/2).*(theta_61 - sin(theta_61)));
     
    L6(L6<0)=0;
    
    L6 = L6 + (2.*R62);
    R6plot = (L6.^2)./(8.*h6) + (h6/2);
    y6plot = (H + R6plot - h6);
   

      F6(:,8) = R62;                                   % Radius of Curvature
      F6(:,7) = R6plot;                                % R for plotting/coalescence
      F6(:,13) = L6;                                   % length along channel LD
      %F6(:,15)= R61.*(1 - cos(theta_wall)).*W;        % Crosssectionarea
      F6(:,15)= h6.*W;                                 %may affect
      F6(F6(:,15) > H*W,15) = H*W;
      %F6(:,14) = (W*L6) + (2.*theta_wall.*R62.*W);    % surface area
      F6(:,14) = (W*L6);
      F6(:,11) = (H - h6)./2;                         % medial axis b (need to check)
      %F6(:,17) = 0;                                   % GDL surface Coverag
      % F6(:,17) = c5b.*W;                            % new estimation 11/06/2020
      F6(:,3) = (W/2);                                %x position
      F6(:,4) = y6plot;                              %y positionneed to change this 27/01/2020
      F6(:,9) = W + 2*(R62*(1-cos(theta_wall)));    %Rpore
      F6(:,16) = F6(:,9);                           %LAD
      F6(:,18) = (F6(:,13).*W) + (2.*h6.*L6);                       %liquid solid area
      F6(:,12) = W;
      F6(:,10) = h6;
     

clear h6 R6plot y6plot L6 R62
  end
  % is the same algorithm as the capillary bridge regime
  if isempty(D7)==0
       clear z571
  
    V7 = F7(:,1);
    z71 = ((W/2)./cos(theta_wall)) * (1-cos((pi()/2) - (theta_wall)));
    R71 = (W/2) * sec(theta_wall);
    theta_71 = (pi() - (2*theta_wall));
    theta_72 = theta_wall*2;
    
    phi_71 = (R71.^2).*(theta_71 - sin(theta_71))./(2*z71*W);
    phi_72 = (theta_72 - sin(theta_72))./(2*pi());
    
    a7 = phi_72*W*pi();
    b7 = -2*phi_72*W*pi()*phi_71*z71;
    c7 = phi_72*W*pi()*phi_71*(z71.^2) - F7(:,1);
    
    R72 = (b7 + sqrt(b7.^2 - (4*a7.*c7)))./(2*a7);
    R73 = R72-z71;
    h7 = R72.*(1 - cos(theta_wall));
    h73 = R73.*(1 - cos(theta_wall));
    %new section 11/06/2020
    
    h7b = (h73 - H);
    h7b(h7b<0)=0;
    
    c7b = 2.*R72.*sin(acos(1 - (h7b./R72))); 
    if theta_wall>(pi()/2)
    c7b = 2.*R73.*sin(acos(1 - (h7b./R73))); 
    end
    c7b(c7b<0)=0;
   
    %NEW DERIVATION OCTOBER 2020
    L7 = F7(:,1)./(H*W);
    
    
      F7(:,8) = R72;                                % Radius of Curvature
      F7(:,7) = R72;                                % R for plotting/coalescence
      F7(:,13) = 2.*R72.*sin(theta_wall);           % length along channel LD
      F7(:,15)= R73.*(1 - cos(theta_wall)).*W;      % Crosssectionarea
      F7(:,15)= h73.*W;                             %may affect 25/06/2020
      F7(F7(:,15) > H*W,15) = H*W;
      F7(:,14) = pi().*(R72).*phi_72*W;             % surface area
      F7(:,11) = (H - h7)./2;                       % medial axis b (need to check)
      F7(:,17) = 0;                                 % GDL surface Coverag
      %F7(:,17) = c7b.*W;                            % new estimation 11/06/2020
      %NEW
      F7(:,8) = L7./2;                                % Radius of Curvature
      F7(:,7) = L7./2;                                % R for plotting/coalescence
      F7(:,13) = L7;           % length along channel LD
      %NEW
      F7(:,17) = F7(:,13).*W;                       %amended 20/10/2020
      F7(:,3) = (W/2);                              %x position
      %F7(:,4) = H/2 + (R72 - h7);                     %y positionneed to change this 27/01/2020
      F7(:,4) = H/2;                     %y positionneed to change this 27/01/2020
      F7(:,9) = W + 2*(R72*(1-cos(theta_wall)));    %Rpore
      F7(:,16) = F7(:,9);                           %LAD
      F7(:,18) = F7(:,13).*W;                       %liquid solid area
      F7(:,12) = h7;
      %new additions 16/06/2020
      F7(:,10) = h73;                    %new droplet height
  end
  
  %% re-establish droplet matrix
  if isempty(D1)==1
F1=zeros(1,size(D,2));
  end
if isempty(D2)==1
F2=zeros(1,size(D,2));
end
if isempty(D3)==1
F3=zeros(1,size(D,2));
end
if isempty(D4)==1
F4=zeros(1,size(D,2));
end
if isempty(D5)==1
F5=zeros(1,size(D,2));
end
if isempty(D6)==1
F6=zeros(1,size(D,2));
end
if isempty(D7)==1
F7=zeros(1,size(D,2));
end
if isempty(D0)==1
F0=zeros(1,size(D,2));
end


               D=[F1;F2;F3;F4;F5;F6;F7;F0];
               D = D(any(D,2),:); 
               clear F1 F2 F3 F4 F5 F6 F0 F7
               
%% ===========  force balance ======================

%First find projected cross sectional area defined from classes
D(D(:,15)>(H*W),15)=(H*W);
%D(D(:,10)>H,15)=H*W;
P_A = D(:,15);

%find the average velocity of air in the segment (in progress)
%D(:,22) = u_a( (D(:,5) - (dz/2) < Cz) & (D(:,5) > Cz) )

 D(:,22) = u_a(1)  ;
 D(D(:,20) > D(:,22),20)= D(D(:,20)>D(:,22),22);
 
%Use mass balance (continuity equation) to find velocity of air
U_I = ((D(:,22).*H*W)./((H*W) - P_A)); 
%U_I = (((D(:,22)-D(:,20)).*H*W)./((H*W) - P_A)); 

%Determine pressure force from bernuolli equation (pressure to kinetic)
DP_p = (rhoa/2).*((((U_I - D(:,20)).^2) - (D(:,22)).^2));

%FPressure = 1.22.*DP_p.*P_A; %Pressure Force = Pressure Difference * Cross section area
FPressure = DP_p.*P_A;
FPressure(isnan(FPressure)==0);

%Determine Shear Force gas-liquid on each droplet (parallel plate with b)
DP_v= (3*mug(1).*(U_I - D(:,20))./D(:,11));
FShear = DP_v.*D(:,14); %new shear 27/03/2020
FShear(isnan(FShear)==0);

%Shear Force liquid_solid
FShear_ls=((muw.*D(:,20))./((D(:,7)/2))).*D(:,18);
FShear_ls(D(:,6)==6) = (6.*(muw.*D(D(:,6)==6,20))./(D(D(:,6)==6,10))).*D(D(:,6)==6,18);
FShear_ls(isnan(FShear_ls)==0);

%Total Drag Force
D(:,23) = FShear + FPressure - FShear_ls;
D(D(:,23)<0,23)=0;

% for plug flow
       %index = find(D(:,6)>5 & D(:,1)>0);
       index = find(D(:,6)==7 & D(:,1)>0);
       Fshearplug=((6*muw*uin)/H)*(2*(H*D(index,13)) + 2*(W*D(index,13))); % N
       DPplug=Fshearplug/(H*W);
       DP3(index,1) = DPplug;
       
       %new oct 2020
        bp =(W)/2;
        cp =(H)/2;
        betap =((1/3) - (64/(pi()^5)).*((cp./bp).*tanh((pi().*bp)./ (2.*cp))));
        DPplug =(uin.*muw.*D(index,13))./(cp.^2 .* betap);
        DPplug(isinf(DPplug)|isnan(DPplug)) = 0;
   

%Adhesion Force based on contact line normal to flow and contact angle
ADV=ones(size(D,1),1).*theta_GDL;
ADV(D(:,6)==1) = theta_GDL;     %emerging
ADV(D(:,6)==2) = theta_GDL;     %isolated
ADV(D(:,6)==3) = theta_wall;    %side wall
ADV(D(:,6)==4) = theta_wall;    %corner droplet
ADV(D(:,6)==5) = theta_wall;    %truncated capillary bridge
ADV(D(:,6)==6) = theta_wall;    %plug
ADV(D(:,6)==7) = theta_wall;    %plug
REC = ADV - theta_hyst;         %find receding angle based on constant CAH
REC(D(:,6)==5) = theta_wall-1;  %hysteresis when not attached to gdl
REC(D(:,6)==6) = theta_wall-1;
REC(D(:,6)==7) = theta_wall-1;
REC(D(:,6)==3) = theta_wall-1;

%D(:,24) = (24/(pi()^3))*sigma*D(:,16).*(cos(REC)-cos(ADV));
D(:,24) = sigma*D(:,16).*(cos(REC)-cos(ADV)); %new adhesion 2021 check
D(D(:,6)==6,24) = sigma*D(D(:,6)==6,16).*(cos(REC(D(:,6)==6))-cos(ADV(D(:,6)==6)));

%Record Pressure Drop from Droplets Alone
DP = DP_v + DP_p;
DP(D(:,6)==7)=DPplug;
 
%% Determine Pressure Drop In Channel

%Single Phase Pressure Drop
   bp=(W)/2;
   cp=(H)/2;
   betap=((1/3) - (64/(pi()^5)).*((cp./bp).*tanh((pi().*bp)./ (2.*cp))));
   DPsinglephase=(uin.*mug.*dz)./(cp.^2 .* betap);
   DP(isinf(DP)|isnan(DP)) = 0;
   DP=abs(DP);
   
   TotalP=sum(DPsinglephase) + sum(DP);
   
   Rel_perm = uin*mug(1)/TotalP;
   perm = uin*mug(1)./sum(DPsinglephase);
   
%%  (Dec 2020) Determine pressure drop along the channel

%    pressure_length = [D(:,5) DP];
%    p_ls = sort(pressure_length,1); %sort pressure based on distance
%    maxpressure = max(cumsum(p_ls(:,2)));
%    jp=2:1:size(D,1)';
%    distance_diff(1,1) = 0;
%   distance_diff(jp,1) = abs(p_ls(jp-1,1)-p_ls(jp,1)); %find distance between consec droplets
%   sp_ls = (uin.*mug(1).*distance_diff./(cp.^2 .* betap));
%    %update pressure length 
%    p_ls(:,2) = p_ls(:,2) + sp_ls;
%   maxpressure = max(cumsum(p_ls(:,2)));
%    p_ls(:,2) = p_ls(:,2)*-1;
%    clear PL
%    clear sp_ls
%    clear distance_diff
%    PL(:,1) = p_ls(:,1);
%    PL(:,2) = maxpressure + cumsum(p_ls(:,2));
%    
%    
%   figure(56)
%   plot(PL(:,1),PL(:,2))

%% ==================  Evaporation (IN PROGRESS) ====================

if evaporation==1
% As=A(:,14) surface area of droplets
% Find channel segments they are evaporating to
% D(:,5) is z location
%Td is temperature of droplets


Location_Evap = D(:,5) > Cz ;%find where the location is less than the channel points
[~,i_evap]=min(Location_Evap,[],2); %find the first 1 value in each row (i_evap)

Cbulk=Cw(i_evap)'; %converting rows to column array for droplet evaporation assesment
Tbulk = T(i_evap)';
%surface is at saturated conditions of temperature of droplet
%Psurf=Psat(T(i_evap)); 
Psurf=Psat(D(:,27)); % need to update Td every timestep
%kelvin eqaution
Psurf=exp(2*sigma*(18/1000)./(D(:,7).*8.314.*D(:,27))).*Psurf;

%concentration at surface in mol/m3 by ideal gas law
%Csurf=(Psurf./(8.314.*Tin))';
Csurf=(Psurf./(8.314.*D(:,27)));

%update physical properties for boundary layer properties
kgS=(1.97*10^-4).*((((D(:,27)+Tbulk)/2)).^0.858);
mugS= (1.691*10^-5)+(4.984*10^-8.* ((((D(:,27)+Tbulk)/2))-273)) -(3.187*10^-11 .* ((((D(:,27)+Tbulk)/2))-273).^2)+(1.319*10^-14 .* ((((D(:,27)+Tbulk)/2))-273).^3);
Dh2oS = 21.6*10^-6 .* (1 + (0.0071 .* (((D(:,27)+Tbulk)/2)-273))); %added diffusion coefficient
%Dh2oS=Dh2o(1);
%Cbulk=Cbulk.*0; % remove this 
%Sherwood number estimation
%Sh=2 + 0.6.*(mug(1)./(rhoa*Dh2o(1))).^(1/3).*((2.*D(:,7).*(u_a(1)-D(:,20)).*rhoa) ./(mug(1)) ).^(1/2);
%Sh=2 + 0.6.*(mugS./(rhoa*Dh2oS)).^(1/3).*((2.*D(:,7).*(u_a(1)-D(:,20)).*rhoa) ./(mugS) ).^(1/2);
Sh=2 + 0.6.*(mugS./(rhoa*Dh2oS)).^(1/3).*((2.*D(:,7).*(U_I-D(:,20)).*rhoa) ./(mugS) ).^(1/2);
%km = Sh.*Dh2o(1)./(2.*D(:,7));

Sh(D(:,6)==1)=2;
%Sh=2 + 0.6.*(mugS./(rhoa*Dh2oS).^(1/3).*((2.*D(:,7).*(u_a(1)-D(:,20)).*rhoa) ./(mugS)).^(1/2))
km = Sh.*Dh2oS./(2.*D(:,7));
km(isinf(V_evap)|isnan(V_evap)) = 0;
V_evap = 0.018.*km.*D(:,14).*(Csurf - Cbulk)./(rhow);
%V_evap = 0.018.*km.*D(:,14).*(Csurf - Cbulk)./(rhow);


% Re=(2.*D(:,7).*(u_a(1)-D(:,20)).*rhoa) ./(mug(1));
% C_bl = 2.*D(:,7)./sqrt(Re);
% C_bl(isinf(C_bl)|isnan(C_bl)) = D((isinf(C_bl)|isnan(C_bl)),11);

%volume evaporation rate DV/dt
%V_evap = 0.018.*(Dh2o(i_evap))'.*D(:,14).*(Csurf - Cbulk)./((D(:,11)).*rhow);
%V_evap = 0.018.*(Dh2o(i_evap))'.*D(:,14).*(Csurf - Cbulk)./(C_bl.*rhow);
V_evap(isinf(V_evap)|isnan(V_evap)) = 0;
%V_evap(V_evap<0)=0;
D(:,19)= V_evap;
%Td
hlv =((-2.5464.*D(:,27))+2503.6).*1000;               %latent heat of vaporisation/condensation
BT=(Cpg(1).*(Tbulk-D(:,27))./hlv);

%Nu=(2 + 0.6.*(Cpg(1).*mugS./kgS).^(1/3).*((2.*D(:,7).*(u_a(1)-D(:,20)).*rhoa) ./(mugS) ).^(1/2))./((1+BT).^0.7);
Nu=(2 + 0.6.*(Cpg(1).*mugS./kgS).^(1/3).*((2.*D(:,7).*(U_I-D(:,20)).*rhoa) ./(mugS) ).^(1/2))./((1+BT).^0.7);
Nu(D(:,6)==1)=2;
h_conv = Nu.*kgS./(2.*D(:,7));
h_conv(isinf(h_conv)|isnan(h_conv)) = 0;

%find which segments the droplet is adding mass to
%Evap = (Dh2o(i_evap))'.*D(:,14).*(Csurf - Cbulk)./(D(:,11)) ;
w_c=w_c*0;
for ie=1:1:size(i_evap);
   
Cw(i_evap(ie))=Cw(i_evap(ie)) + ((V_evap(ie).*rhow./0.018).*dt)./(dz*H*W);
%T(i_evap(ie)) = T(i_evap(ie)) - (h_conv(ie).*D(ie,14).*(Tbulk(ie) - D(ie,27)).*dt);
    %+(hlv(ie).*(D(ie,19).*rhow)./(rhog(i_evap(ie)).*Cpg(i_evap(ie))).*(H*W*dz)).*dt;
%
w_c(i_evap(ie)) = w_c(i_evap(ie)) + D(ie,19);
end

%Sw_c = w_c./(dz*W);
%%update properties
%Physical Properties of cell
Dh2o=21.6*10^-6 .* (1 + (0.0071 .* (T-273)));   %diffusion coefficient of water vapour
Cpg=(1.00926*10^3)-(4.0403*10^-2 *(T-273)) + (6.1759*10^-4 *(T-273).^2) - (4.097*10^-7 *(T-273).^3);    %heat capacity of air
rhog=Pin./(287.058*T);                          %density of air
kg=(1.97*10^-4).*(T.^0.858);                    %thermal conductivity of air
alphag=kg./(Cpg.*rhog);                         %thermal diffusivity of air
mug= (1.691*10^-5)+(4.984*10^-8.* (T-273)) -(3.187*10^-11 .* (T-273).^2)+(1.319*10^-14 .* (T-273).^3);
hlv =((-2.5464.*D(:,27))+2503.6).*1000;               %latent heat of vaporisation/condensation                  
AverageRH(ic,:) = RH;

%%Initialise simulation with droplets already present

end

Ps=Psat(T);
RH=(Cw.*8.314.*T)./Ps;

%% =============== Condesnation (IN PROGRESS) =================

 if condensation==1
     
Cc=(RH-1).*Psat(T)./8.314.*T;
Vcond=((RH-1).*Psat(T)*0.018)./(8.314.*T.*rhow);
Vcond(Vcond<0)=0;

     if Vcond>0
     icond=1:1:size(Cz,2);
     Di(icond,:)=Dorig(1,:);
     Di(icond,1) = Vcond(icond);
     Di(icond,3) = round(rand(1).*W,3);
     Di(icond,4) = H;
     Di(icond,5) = Cz(icond);
     Di(icond,7) = ((3.*Di(icond,:))./(4*pi()))^(1/3)
     Di(:,6) = 4
     D=[D;Di];
     end
     
     %Cw=Cw-Cc
 end

timetest(ic,4)=toc;

%% ============== Saving Data =======================
  
   S(ic,1)=t;
   S(ic,2)=TotalP;
   S(ic,3)=sum(V_evap);
   S(ic,4)=mean(D(D(:,1)>0,20));

   Time(ic,1)=toc;
   Time(ic,2)=size(D,1);
   Time(ic,3)=sum(D(:,1))./(L*W*H);
   Time(ic,4)=sum(D(:,17))./(L*W);
   
   %save condensation data
   Water(ic,1)=water_cond;
   % save evaporation data
   E_data(ic,1)=D(1,1);
   E_data(ic,2)=D(1,27);
   
   
%save relative permeability estimation
permsave(ic,1) = Rel_perm./perm;

%visualiseDroplets2
visualiseDropletCounter = visualiseDropletCounter + 1;

if visualiseDropletCounter > vdclimit
visualiseDroplets4
visualiseDropletCounter = 0;
end




   
%% ============== Figure Plotting During Simulation =================

if timeplot==1
%_____SATURATION______
  
  figure(1)
  clf('reset')
  plot(S(:,1),Time(:,3),'k')
  hold on
  plot(S(:,1),Time(:,4),'b')
  ylabel('Saturation/Water coverage ratio')
  yyaxis right
  plot(S(:,1),S(:,2),'r')
  set(gca, 'YScale', 'log')
  ylabel('Channel Pressure (Pa)')
  xlabel('Time (s)')
  set(gca,'fontsize',16)
  legend('Channel Saturation','Channel water coverage ratio','Channel air pressure (Pa)')
  
  
end


%_____Water coverage ratio______
% 
%  if rem(ic,2)==1
% figure(2)
%  plot(S(:,1),Time(:,4))
%  hold on
%  end
 
 %set(gca, 'XScale', 'log')


%_____Number of Droplets______
%  figure(3)
%  plot(S(:,1),Time(:,2))


%_____Channel Pressure drop______
%  figure(4)
%  plot(S(:,1),S(:,2))


%_____relperm______
% figure(4)
%  %plot(S(:,1),permsave(:,1))
%   plot(Time(:,3),permsave(:,1))
%plot(Time(:,3),permsave(:,1),'-s')


%% channel plotting X - Z cross section (uncomment all)
% 
% figure(2)
% plot(zp , xp ,'mo','MarkerFaceColor','m')
% viscircles([D(:,5) D(:,3)],D(:,7),'color','k')
% legend( sprintfc('%g s', t),'FontSize',12 )
% xlim([0-(dz/2) (L+(dz/2))/100]);
% xlim([0 (L+(dz/2))/100]);
% 
% ylim([0 W])
% set(gcf, 'Position',  [0, 500, 2000, 250]) %1400 %100

%% channel plotting X - Y cross section (uncomment all)
% 
% clear figure
% sumdt=0;
% figure(3)
% clf('reset')
% viscircles([D(:,3) D(:,4)],D(:,7),'Color','b','LineWidth',2)
%  line([0,0],[0,W],'LineWidth',0.5,'Color','k')
%  %hold on
%  line([W,0],[0,0],'LineWidth',2,'Color','k')
%  line([0,W],[H,H],'LineWidth',2,'Color','k')
% xx1=0;
% xx2=W;
% yy1=0;
% yy2=H;
% xxx = [xx1, xx2, xx2, xx1, xx1];
% yyy = [yy1, yy1, yy2, yy2, yy1];
% %plot(xxx, yyy, 'k-', 'LineWidth', 1);
% %hold on;
% xlim([0 W]);
% ylim([0 H]);
% %set(gcf, 'Position',  [0, 500, 800, 300])
% legend( sprintfc('%g s', t),'FontSize',12 )
% 
%% serpentine_plot (only if file Serpentine_Plot.m is in folder)
if serpentine_plot==1
    Serpentine_Plot
end
vtkploti=vtkploti+1;
%% ============= Vtk Plotting =====================
% requires vtkwrite.m to be in folder

if vtk==1
    if vtkploti>vtkplot
         %filename=[num2str(t) '.vtk']
         %filename=['VTK_single/' num2str(ic) '.vtk'];
          filename=[vtkname num2str(ic) '.vtk'];
         % filename=[vtkname num2str(t) '.vtk'];
         
         x_vtk=D(:,3);
         y_vtk=D(:,4);
         z_vtk=D(:,5);
         R_vtk=D(:,7);
         DP_vtk=DP;
         u_vtk=D(:,20);
         regime_vtk=D(:,6);
         T_vtk=D(:,27);
         
         %newvtk(filename,'polydata','lines',x_vtk,y_vtk,z_vtk,R_vtk)
         %vtkwrite(filename,'structured_grid',x_vtk,y_vtk,z_vtk,'scalars','radius',R_vtk,'scalars','Pressure',DP_vtk,'scalars','velocity',u_vtk,'scalars','regime',regime_vtk) %this allows for saving of scalar data
         vtkwrite(filename,'unstructured_grid',x_vtk,y_vtk,z_vtk,'scalars','radius',R_vtk,'scalars','Pressure',DP_vtk,'scalars','velocity',u_vtk,'scalars','regime',regime_vtk,'scalars','Temperature',T_vtk) 
         sumdt=0;
     
     clear x_vtk y_vtk z_vtk R_vtk DP_vtk u_vtk regime_vtk
     vtkploti=0;
    end
end     


%% ========== Determine Droplet Regimes ============

RegimeVolume = [sum(D((D(:,6)==1),1));sum(D((D(:,6)==2),1));...
              sum(D((D(:,6)==3),1));sum(D((D(:,6)==4),1));sum(D((D(:,6)==5),1));sum(D((D(:,6)==6),1))...
              ;sum(D((D(:,6)==7),1))];
          
RegimeSaturation = RegimeVolume./sum(RegimeVolume);
RegimeSatSave(ic,:)=RegimeSaturation;


%%=========== Determine Regimes Along the Channel (INPROGRESS) ================%%
ct=ct+dt;
 if ct > write_precision
nseg=(nz)/2;                    % number of channel segments
nseg=20;                        % ''
Dseg=L/(nseg);
segment=transpose([0:Dseg:L]);

if time_value==0                %initialisation of array
    avregime_new=zeros(length(segment),8);
    plottimer=0;
    av_vol_frac_i=zeros(1,length(segment));
end

%save previous timestep average
avregime_i = avregime_new;

%determine the fraction of each regime within each segment of the channel
for is=1:1:length(segment)
segment(is,2)=sum(D((D(:,5)<segment(is,1)+Dseg) & (D(:,5)>(segment(is,1)-Dseg)) &(D(:,6)==1) & (D(:,1)>0),1));
segment(is,3)=sum(D((D(:,5)<segment(is,1)+Dseg) & (D(:,5)>(segment(is,1)-Dseg)) &(D(:,6)==2) & (D(:,1)>0),1));
segment(is,4)=sum(D((D(:,5)<segment(is,1)+Dseg) & (D(:,5)>(segment(is,1)-Dseg)) &(D(:,6)==3) & (D(:,1)>0),1));
segment(is,5)=sum(D((D(:,5)<segment(is,1)+Dseg) & (D(:,5)>(segment(is,1)-Dseg)) &(D(:,6)==4) & (D(:,1)>0),1));
segment(is,6)=sum(D((D(:,5)<segment(is,1)+Dseg) & (D(:,5)>(segment(is,1)-Dseg)) &(D(:,6)==5) & (D(:,1)>0),1));
segment(is,7)=sum(D((D(:,5)<segment(is,1)+Dseg) & (D(:,5)>(segment(is,1)-Dseg)) &(D(:,6)==6) & (D(:,1)>0),1));
segment(is,8)=sum(D((D(:,5)<segment(is,1)+Dseg) & (D(:,5)>(segment(is,1)-Dseg)) &(D(:,6)==7) & (D(:,1)>0),1));
end

 segmentRegime=segment;
 segmentRegime(:,2:end)=segment(:,2:end)./sum(segment(:,2:end),2);

 time_value=time_value+1;
 
 %%% ==========  TRACKING AVERAGE WITHOUT SAVING ALL TIMESTEPS =================
 segmentRegime(isnan(segmentRegime))=0; %remove any values which /0
 avregime_new = (((avregime_i).*(time_value-1))+segmentRegime)./((time_value)); %use new averaging procedure

%_______Plot visualisation of regimes along channel__________
%   figure(27)
%   %bar(segmentRegime(:,1),segmentRegime(:,2:end),'stacked')           %Bar representation 
%   area(segmentRegime(:,1),segmentRegime(:,2:end),'facealpha',0.3)     %Area representation
%   xlabel('Distance Along Channel (m)')
%   ylabel('Volume fraction of total water in channel')
%   legend('Emerging','Isolated','Side Wall','Corner','Capillary Brdige','film')
%   ylim([0 1])  
 
 
%% ========== calculate volume fraction over length =================

Vp = D(:,1);                        %volume of droplets
zdrop = D(:,5);                     %z position of droplets
Ap = D(:,15);                       %projected area of droplet
Lp = D(:,13);
Regime_vol = D(:,6);

% zp +- LD/2  will find the maximum range of droplet in z direction.
zdrop_min = D(:,5) - (D(:,13)./2);   %minimum z position of droplets
zdrop_max = D(:,5) + (D(:,13)./2);   %minimum z position of droplets

%assum constant variation of volume over its length in the channel
NZ = nseg+1;
DZ = L/NZ;

%define new grid for cell centers of channel
Z = (DZ/2):DZ:L-(DZ/2);

%Distance from particle p and cell k

Spk = sqrt((Z-zdrop).^2);
Lpk = (Lp/2) + (DZ/2);

%overlap length
overlap = (Lpk-Spk);
%all non-positive values are zero
overlap(overlap<0)=0;
%all overlap greater than the cell length is cell length
overlap(overlap>DZ)=DZ;
%divide overlap by the length of the droplet in length direction
overlap=overlap./Lp;
overlap(isnan(overlap)|isinf(overlap))=0;
overlap(overlap>1)=1;
water_cell_volume = overlap.*Vp;
regime_cell_volume = overlap.*Regime_vol; %%% change this maybe
vol_fraction = sum(water_cell_volume)./(DZ*W*H);

%vol_frac_time(:,time_value)=vol_fraction;

av_vol_frac = (((av_vol_frac_i).*(time_value-1))+vol_fraction)./((time_value)); 

% figure(100)
% bar(segmentRegime(:,1),avregime(:,2:end).*mean(vol_frac_time,2),'stacked')
av_vol_frac_i = av_vol_frac;
% figure(100)
% clf('reset')
% plot(av_vol_frac,'-x')

% 
 % figure(101)
  %bar(segmentRegime(:,1),avregime_new(:,2:end).*av_vol_frac','stacked')
 % bar(segmentRegime(:,1),avregime_new(:,2:end),0.9,'stacked')
 %area(segmentRegime(:,1),avregime_new(:,2:end).*av_vol_frac','facealpha',0.3)     %Area representation
 %area(segmentRegime(:,1),segmentRegime(:,2:end).*vol_fraction','facealpha',0.3)     %Area representation
% 
%  equiv_diameter = (6.*D(:,1)./pi()).^(1/3) .* 1e+6;
%  histogram(equiv_diameter)
%  xlim([0 2000])
 
 ct=0;
 end

 %End of the time loop for simulation
 %calculate percentage complete time
 pt = (t./(endtime)).*100;

end

% End of the time step loop
% Moving on to saving results for each run if it is turned on

repeat % displays the simulation iteration number for tracking progress

 if saveresults==1        
    % Automatic saving of average properties about the run in a .csv file.
    % Label file name to correspond to number of run i
    filename=[Savefilename num2str(sid) '.csv']
    
    % Create heading titles to identify what is saved.
    SaveDataTitles = ['u_air' 'u_water' 'GDL_angle' 'Wall_angle' ...
    'npores' 'Max_Pore_Radius' 'L' 'H' 'W' 'Total CPU Time (s)' 'endtime'...
    'dt' 'Saturation' 'ACR' 'Pressure Drop' 'Ndrops' 'mean vol frac emerging'...
    'mean vol frac isolated' 'mean vol frac sidewall' 'mean vol frac corner' 'mean vol frac capillary' 'mean vol frac film'];
  
   % Create list of saved variables to identify operating and average conditions.
    SaveData = [uin u_wav theta_GDL.*(180/pi()) theta_wall.*(180/pi())...
    npores Rpmax L H W sum(sum(timetest)) S(end,1) dt mean(Time(:,3)) mean(Time(:,4)) mean(S(:,2)) mean(Time(:,2)) mean(RegimeSatSave) current_density Sto];
    
   % Write the data and the headings to a .csv file
    csvwrite(filename,SaveData);

   % After the first row in the .csv, now save the transient data for the simulation at every timestep.
                    %time  , saturation, ACR,   pressure ndrops . regimes , droplet velocities av           
    dlmwrite(filename,[S(:,1) Time(:,3) Time(:,4) S(:,2) Time(:,2) RegimeSatSave S(:,4)],'delimiter',',','-append','precision',9);

    % save the file
    filename2=[Savefilename num2str(sid) '_length' '.csv'];
    % saving new data for regime distribution along the channel
    LengthData = [segmentRegime avregime_new vol_fraction' av_vol_frac'];
    csvwrite(filename2,LengthData);

    % write vtk of last time step value (to see the final distribution using paraView.

         filename = [vtkname num2str(repeat) '.vtk'];
         x_vtk = D(:,3);
         y_vtk = D(:,4);
         z_vtk = D(:,5);
         R_vtk = D(:,7);
         DP_vtk = DP;
         u_vtk = D(:,20);
         regime_vtk = D(:,6);

         vtkwrite(filename,'unstructured_grid',x_vtk,y_vtk,z_vtk,'scalars','radius',R_vtk,'scalars','Pressure',DP_vtk,'scalars','velocity',u_vtk,'scalars','regime',regime_vtk) 
         sumdt = 0;   
         clear x_vtk y_vtk z_vtk R_vtk
 end


end

% To visualise the results, the best way is to load .vtk files into paraView. Then use Glyph filter -> Spheres -> adjust radius using the scalar of
% radius. Then clip the spheres where the channel location would be, this will give you the correct distribution.






