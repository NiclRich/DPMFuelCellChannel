import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix, triu
import matplotlib.pyplot as plt
import csv


# Model initialisation
## Model Visualisation and Saving Settings
repeatmax = 1

# insert for loop which repeats the running of the model,
# later it will be extended to use several runs, but at first only one run
for repeat in range(repeatmax):
    sid = repeat

    print(f"Starting Run Number: {1}") # TODO replace the 1 with repeat later
    # Those variables are used for the control of the program, currently those are only
    # 1/0 values. TODO: Replace it with Boolean values
    # TODO: check clear variables
    vtk = 1                     # Vtk plotting of results at each time step
    vtkplot = 1                 # how many iterations to skip for vtk saving
    serpentine_plot = 0         # Vtk plotting in serpentine arrangement
    coalescence = 1             # Droplet colision and merge Algorithm

    transport = 0               # Transport Algorithm on (1 for on) (IN PROGRESS) TODO: What is this?
    evaporation = 0             # Evaporation Algorithm on (1 for on) (IN PROGRESS) TODO: What is this?
    condensation = 0            # Condensation Algorithm (IN PROGRESS) TODO: Does it need to be included?
    film_formation = 1          # Film Formation Regime (1 for on)
    plug_formation = 1          # Plug Formation Regime (1 for on)
    ct = 0                      # Iteration counter
    saveresults = 1             # Turn on to save results in .csv and .vtk of final time step
    Savefilename = "results_"   # template for the file name to save results
    timeplot = 0                # If you want to see transient Saturation, WCR, pressure
    vdclimit = 2                # After how many interation to plot the surface
    gridRes = 50                # maximum resolutoin for grid (slow plotting)

    ## Model Initial and Boundary Conditions

    # Channel Dimensions
    H = 0.001                   # channel height (m)
    W = 0.001                   # channel width (m)
    L = 0.01                    # channel length (m)
    # TODO: What means the following line in the original code?
    # L = channel_length(repeat)

    xmin = 0                    # min xaxis
    xmax = W                    # max xaxis
    ymin = 0                    # min yaxis
    ymax = H                    # max yaxis

    # GDL and Channel Wetting Properties TODO: What is GDL?
    theta_GDL = 120 * (np.pi / 180)     # GDL contact angle (radians)
    theta_wall = 70 * (np.pi / 180)     # wall contact angle (radians)
    theta_hyst = 25 * (np.pi / 180)     # contact angle hyteresis (radians)

    # Discretise channel points for transport and regime plotting along channel
    nz = 10                     # number of nodes
    dz = L / nz                 # cell size
    Cz = np.arange(0,L + dz,dz)      # location of nodes, endpoints included

    # Electrochemical Operating Parameters
    # TODO: What does the following (commented out) line in the original code:
    # Sto = OPC(repeat, 2) % Stochiometry of o2 reduction reaction (read from file)
    Sto = 2                     # Stochiometry
    Active_Area = L * W * 2     # Active Area (L*W) + 2 half channel ribs = 2 * L * W TODO: What does this mean?
    Active_Area = 1.5125e-05    # Or specify cell cross section
    # current = OPC(repeat, 1) * Active_Area # Specify current externally
    current = 1e4 * Active_Area # Operating current (must change numerical value)
    current_density = current / Active_Area # Current density based on active area

    # Operating Parameters needed for transport algorithm
    xo2 = 0.21                  # fraction of O2
    xN2 = 1 - xo2               # fraction of N2
    RHi = 1                     # Inlet RH TODO: What does this variable
    RHd = 1                     # Domain RH TODO: What does this variable?

    # Temperature are in Kelvin (K) up to 0.15K
    Tin = 70 + 273              # Inlet Air Temperature (K)
    Tdin = 70 + 273             # Droplet temperature (K)
    Tdomain = 70 + 273          # Temperature of Domain (K)
    Pin = 120000                # Inlet Pressure (Pa)
    T = Tdomain * np.ones(Cz.size)    # Nodal temperature of domain
    # define anononymous function, which is later used several times
    # Saturated vapour pressure equation (Pa), TODO: make clear what does it means
    Psat = lambda T: np.exp(34.494 - (4924.99 / ((T-273) + 237.1))) / (((T-273) + 105) ** 1.57)
    Ps = Psat(T)                # Initialising Psat in domain
    Pv = RHd * Ps               # Initialising vapour pressure in domain
    Pvin = RHi * Psat(Tin)      # Calculating water vapor pressure in inlet
    CWin = Pvin / (8.314 * Tin) # Calculating water concention inlet (mol / m3)
    Cw = Pv / (8.314 * T)       # Determine water concentration in domain (mol / m3)
    xw = Pv / Pin               # Water mol fraction
    # TODO Check the size of the expression, it's confusing if it is a scalar or an array
    u_sto = (Sto / xo2) * ((current_density * W * L) / (4 * 96487)) * ((8.314 * Tin) / (Pin - Pvin)) * (1/(H*W)) # stochiometric air velocity (m/s)
    Oin = ((xo2 * (Pin - Pvin)) / (8.314 * Tin)) # oxygen concentration at inlet
    N2 = (0.79 / 0.21) * Oin    # Nitrogen concentration
    rhow = 1e3                  # density of water (kg / m3)
    rhoa = 1                    # density of air (kg/m3)
    sigma = 0.072               # surface tension of water
    w_c = np.zeros(Cz.size)     # water in channel segments
    Sw_c = w_c                   # Saturation of water in channel segments


    ## Two phase Droplet Injection Parameters
    # Number of Injection points (i.e. number of discrete water clusters in GDL
    # TODO: Insert commented out lines
    npores = np.round(7*(W*L*1000000)) # set number of injection pores (7 clusters / mm2)
    npores = int(npores) # Python expects in integer

    # ========== Specify the pore size and location ========== 
    Rpmax = (50e-6)/2           # maximum pore radius
    #Rpmax = 100e-6             # maximum pore radius (VoF compare) commented out in the orignal
    # create a seed for replication, the original matlab code does not include a seed
    np.random.seed(42)
    xp = np.random.rand(npores) * (W-2 * Rpmax) + Rpmax # x coord of poores ensures no pore touches wall due to Rpmax
    # xp = np.random.rand(npores) * W   # x coord of pores (random)  - commented out in original file
    # xp = np.ones(npores) * (W / 2)    # x coord of pores (ordered) - commented out in original file
    yp = np.zeros(npores)       # y coord of pores (set to zero)
    zp = np.random.rand(npores) * (L / 1.0005) # z coord of pores (random along channel)
    #zp = np.ones(npores) * L / 2   # z coord of pores (ordered)
    # TODO: include the lines for reading the csv


    # ========== Preallocate of droplet data ========== 
    poresize = np.ones(npores) * Rpmax # Assign Rpmax to all the pores
    Vd = np.zeros(poresize.size)            # volume of droplets
    Rc = poresize                           # contact radius
    Rd = poresize                           # radius
    h = np.ones(poresize.size) * 1e-9       # height
    b = np.zeros(poresize.size)             # medial axis
    Lc = np.zeros(poresize.size)            # characteristic length (wall transition)
    Ld = np.zeros(poresize.size)            # length droplet in z direction
    E = np.zeros(poresize.size)             # Equation number
    v_evap = np.zeros(poresize.size)        # evaporation volume
    A_surf = np.zeros(poresize.size)        # area covering GDL surface
    A_sl = np.zeros(poresize.size)          # solid liquid area coverage
    flow_l = np.zeros(poresize.size)        # flow rate into droplet
    u_d = np.zeros(poresize.size)           # velocity of droplets
    u_i = np.zeros(poresize.size)           # velocity of air above droplets
    u_ac = np.zeros(poresize.size)          # velocity of air in channel section
    As = np.zeros(poresize.size)            # surface area of droplets
    Ac = np.zeros(poresize.size)            # cross section area droplet
    Lad = np.zeros(poresize.size)           # length of adhesion contact line
    Fdrag = np.zeros(poresize.size)         # drag force
    Fadhesion = np.zeros(poresize.size)     # adhesion force
    theta_app = np.zeros(poresize.size)     # apparent contact angle
    Rp = np.ones(poresize.size) * poresize  # radius of pore
    hlim = np.ones(poresize.size)           # height limit of slug
    Td = np.ones(poresize.size) * Tdin      # TODO: Temperature of what?


    # Determine initial adhesion force based on advancing receding contact angles
    adv = theta_GDL                         # advancing contact angle
    rec = theta_GDL - theta_hyst            # receding contact angle
    Lad = 2 * Rp                            # contact line normal to flow (adhesion force length)
    Fadhesion = (24/(np.pi**3)) * sigma * 2 * Lad * (np.cos(rec) - np.cos(adv)) # adhesion force (N)

    # ========== Water Flow Rates by current density ========== 

    #mass_flow_water = 1e-6                 # Set water injection rate by mass flow rate (kg/s)
    mass_flow_water = (current_density * W * L * 0.018)/(2*96487) # Set water injection rate by current density (m3/s)
    volume_flow_water = mass_flow_water / rhow # convert mass to volume flow (m3/s)
    #volume_flow_water = ((2400/3600)*10**-6)/1000  # set water injection rate by volume flow rate (m3s/s)
    volume_flow_water_pore = volume_flow_water / npores # volume flow per pore (m3 / s*pore)
    velocity_flow_water = volume_flow_water_pore / (np.pi * poresize**2) # velocity of water in each pore (m/s)
    u_w = velocity_flow_water               # velocity of water  in each pore (m/s)
    u_wav = u_w[0]                          # TODO: What's that?

    # ========== Water Injection by Velocity ==========
    u_wav = 0.1                             # uncomment to set water velocity per pore explicit, original file uncommented

    # Used when flow rate is specified
    volume_flow_water_pore = u_wav * (np.pi * (poresize[0]**2)) # m / s * pore
    volume_flow_water = volume_flow_water_pore * npores

    # Assign a flow rate to each pore
    flow_l = np.ones(poresize.size) * volume_flow_water_pore    # m3/ s* pore

    #========== Preallocate Droplet storage Array ========== 
    # with original pore values saved in Dorig TODO: Make the storage model in Matlab clear
    Dorig = np.asarray([Vd, flow_l, xp, yp, zp, E, Rd, Rc, Rp, h, b, Lc, Ld, As, Ac, Lad, A_surf, A_sl, v_evap, u_d, u_i, u_ac, Fdrag, Fadhesion, theta_app, hlim, Td])
    D = Dorig.transpose() # save pore injection conditions in Dorig
    Dorig = Dorig.transpose()
    # D is the main 70x27 array. Each Line is one droplet (so 70 lines) and this keep track of the changes
    # while Dorig uses the start values

    # Fuel Cell Dimensions (not used)
    GDLy = 200e-6                           # Thickness of GDL (unused)
    MPLy = 10e-6                            # Thickness of MPL (unused)
    CLy = 10e-6                             # Thickness of CL

    # Physical Properties of cell
    Do2 = 2.19e-5                           # Diffusion coefficient of oxygen in nitrogen (m2/s)
    Dh2o = 21.6*10**(-6) * (1 + (0.0071 * (T-273))) # Diffussion coeffficient of water vapour (m2 /s)
    Cpg = (1.00926*10**3)-(4.0403*10**(-2) * (T-273)) + (6.1759*(10**(-4)) * (T-273)**2) - (4.097*10**(-7) * (T-273)**3) # heat capacity of air
    Cpw = 4181                              # heat capacity water (J K-1 g-1)
    rhog = Pin / (287.058 * T)              # density of air (kg/m3)
    kg = (1.97*10**(-4)) * (T**0.858)       # thermal conductivity of air (W m-1 K-1)
    alphag = kg / (Cpg * rhog)              # thermal diffusivity of air (W m-1 K-1)
    mug = (1.691 * 10**(-5)) + (4.984*10**-8 * (T-273)) - (3.187*10**-11 * (T - 273)**2) + (1.318*10**-14 * (T-273)**3) # kg m-2 s-1
    klv = ((-2.5464 * D[:, 26]) + 2503.6) * 1000 # latent heat of vaporisation / condensation (KJ / kg) TODO This will result in a bug, since it access D
    muw = 0.0006539                         # viscosity of water (kg m-2 s-1)
    Mw = 0.018                              # molecular mass kg/mol


    # ========== Initialise Species Matrix (Transport) ========== 
    # T, O2, H20

    #Tdin = 60 + 278                        # temperature of Inlet air (K) -- commented out in original
    #Tdomain = Tin                          # Make domain same as inlet 
    Co2 = np.ones(Cz.size) * Oin            # Concentration of O2 (mol/m3)
    # Co2 = np.zeros(Cz.size) * Oin         # zeros (mol/m3)
    # T = np.ones(Cz.size) * Tdomain        # Preallocate domain T
    P = np.ones(Cz.size) * Pin              # Preallocate domain P

    # Droplet Temperature (Evaporation)
    # Td = np.ones(D[:,0].size) *Tdin
    x_o2 = Co2 * (8.314 * T) / P
    Saturation = np.zeros(Cz.size)


    # ========== Set Air Inlet Velocity ========== 
    uin = u_sto                             # air velocity based on stochiometry and current density (m/s)
    uin = 10                                # air velocity based on set value (m/s)
    # IF air velocity (uin = 0), air must flow induced by growth of droplets
    # uin = (volume_flow_water / (H*W) * dt # i.e. injection of droplets must cause air flow

    u_a = np.ones(Cz.size) * uin            # assign velocity to channel points (unused)
    V_air_in = uin * (H * W)               # volume flow of air in (unused)
    V_air = u_a * ( H * W)                  # volume flow of air in each segment (unused)

    # Find where the location is less than the channel points
    Location_Evap = D[:, 4][:, np.newaxis] > Cz  # [:, np.newaxis] reshapes Cz to align with D[:, 5] for broadcasting

    # Find the first occurrence of 1 (True) in each row
    i_evap = np.argmin(Location_Evap, axis=1)

    # Extract corresponding values from Cw and T, converting them to column arrays
    Cbulk = Cw[i_evap]
    Tbulk = T[i_evap]
    # # This is a possible bug since Cz has a different shape than D
    # # After testing: Yes, it is
    # Location_Evap = D[:, 4] > Cz            # find where the location is less than the channel points. This is an array of booleans and not 0-1 like in matlab
    # i_evap = np.argmin(Location_Evap != 1, axis=1) # find the first true value in each row (i_evap) 

    # Cbulk = Cw[i_evap].T                     # converting row to column array for droplet evaportaion assesment
    # Tbulk = T[i_evap].T                     # same ?
    # the following line is definitely a BUG with D
    Nu = 2 + 0.6 * (Cpg[0] * mug[0] / kg[0])**(1/3) * ((2 * D[:, 6] * (u_a[0] -D[:,19])* rhoa) / (mug[0]))**0.5     # nusselt number
    h_conv = Nu * kg[0] / (2 * D[:, 7]) # convective heat transfer coeffficient
    h_conv = np.nan_to_num(h_conv, posinf=0.0, neginf = 0.0) # this does not make sense NR

    ## Balacne Viscous, Intertial and Capilaary Shear Forces (Film Formation == 1)

    Ca = uin * mug[0] / sigma
    A = (rhoa * (uin**2)) / (5 * sigma * (1-np.cos(theta_wall)))
    # coefficients for equation
    a = A
    b = 1 + (3 * Ca * H /(1- np.cos(theta_wall))) - (2 * A * H)
    c = - ( 3* Ca * H / (1 - np.cos(theta_wall))) - (2 * H)
    d = H**2

    # Newton iteration to solve
    f  = lambda h: a * (h**3) + b * (h**2) + c * h + d
    df = lambda h: (3 * a * h**2) + 2 * b * h + c

    h = np.ones(np.size(uin)) * 0.0 # = 0 # really?
    holde = h
    h1 = h
    hnew = holde + 10000

    f1_error = 1e-9
    iter_i = 0
    w = 0.7

    # TODO In the Newton Iteration there is a bug, check it
    # while np.sqrt(np.abs(np.transpose(holde-hnew) @ (holde - hnew))) > 1e-12: # solving the height equation (NEwton Method)
    #     hnew = h1 - f(h1) / (w * df(h1))  #  np.linalg.solve(w * df(h1), f(h1)) not working since it is 1d
    #     holde = h1
    #     h1 = hnew
    #     iter_i += 1

    while (np.abs(holde- hnew) > 1e-12):
        hnew = h1 - f(h1) / ( w * df(h1))
        holde = h1
        h1 = hnew
        iter_i += 1

    h1 = np.minimum(h1, H)  # Cap values in h1 to H (if h1 > H, set it to H)
    h1 = np.maximum(h1, 0)  # Cap values in h1 to 0 (if h1 < 0, set it to 0)

    hlim = np.ones(poresize.size) * h1 # maximum height before deformation in capillary brdige
    D[:, 25] = hlim
    Dorig[:, 25] = hlim # height limit for each droplet based on the external forces

    ## ========== Simulation Time Integration Setting ==========
    dt = 1e-7               # timestep (may be unstable for too large flow rate)
    dtorig_i = 1            # original timestep NR: What's that?
    #dtorig_i = 0.5 * H / uin # if air velocity is high (0.5 courant number) - remove for limiting time step
    dt3 = 0.5 * H / u_wav   # if water velocity is high
    dtorig = min(dtorig_i, dt3) # finds the minimum between the two
    dt2 = dtorig

    # dtorig = 1e-3         # alternatively, set the max time step

    Evap = 0                # initialise evaporation
    V_evap = 0              # Initialise evaporation
    ic = 0                  # iteration coutner
    plotlim = 1             # plotting limit (how many timesteps to skip
    vtkploti = 0            # NR: What's that?
    plotlimcounter = 0      # plotting limit coutner
    write_count = 0         # NR: What's that?
    av_vol_fraction = 0
    maxdroplets = 10000     # determine exit condition of divergence of model

    # Determine simulation time defined by 1 Channel volume of water injected
    Channel_Volume = H * W * L                   # Channel Volume (m3) TODO: Check if unit is correct
    endtime = 3 * (Channel_Volume / volume_flow_water) # Simulation will stop at this time (now 3 time longer to be able to get better time average)
    # endtime = 100                             # alterantively, set endtime
    t = 0                                       # initial time
    time_value =0                               # time coutner
    nwrite = 2                                  # change number to change write interval
    write_precision = 1 * 10**(nwrite)
    timer_tic = 0                               # NR: What's that?
    timer_tic_coal = 0
    courant = 0.1
    visualiseDropletCounter = 0


    ## Start time Intergration
    for tt in range(10): # maximum time loop
        # Initialise simulation with droplets already present
        if t == 0:
            pass
            # Example of intoducing artificial droplets
            # Di = Dorig[1,:]
            # Di[:, 1] = 5.98e-10
            # Di[:, 3] = W/2
            # Di[:, 4] = 2 * H
            # Di[:, 5] = L / 2
            # Di[:, 7] = 0.001
            # Di[:, 6] = 5
            # Di[:, 20] = 0
            # D = np.concatenate((D, Di), axis=1)
            # D[:, 1] = 5.98e-11 # m3
        
        # adaptive time step
        adapt_t = courant * (0.001) / D[:, 20]
        adapt_t = np.nan_to_num(adapt_t, posinf = dtorig)
        dt = np.asarray([np.min(adapt_t)])
        
        # NR: What does this condition?
        if np.max(D[:, 22] / D[:, 23]) > 0.999: 
            dt2 = np.asarray([1e-3])

        dt2 = np.asarray([dt2])
        dt2[dt2 == 0] = dtorig  # This can be a BUG, since it compares directly with 0
        dt3 = 0.5 * dz / uin    # new addition for transport
        
        dt3 = 0                 # remove if transport not needed
        dt3 = np.asarray([dt3])
        if dt3 == 0:
            dt3 = dtorig
        # dt3[dt3 == 0] = dtorig
        dt = min(dt, dt2, dt3)
        dt2 = 0

        #dt = 1e-3
        #dt = 1e-4  # activate to set timestep

        t += dt
        if t > endtime:
            break
        ic += 1
        plotlimcounter += 1

        if transport == 1: # solves trasnport equations if transport == 1
            # transport of water vapour explicit
            Cwn = Cw
            Tn = T
            Co2n = Co2
            V_airn = V_air
            i = np.arange(start=1, stop=Cz.size, step=1)
            Cw[:, i] = Cwn[:, i] + dt * ((u_a[:, i-1] * Cwn[:, i-1] - u_a[:, i] * Cwn[:, i]) / dz) + dt * (Dh2o[:,i] * ((Cwn[:, i-1] - 2*Cwn[:, i] + Cwn[:, i+1]) / dz**2))

            # Transport of temperature
            # TODO

            # Transport of Oxygen
            Co2[:,i] = Co2n[:, i] + dt * ((u_a[:, i-1] * Co2n[:, i-1] - u_a[:, i] * Co2n[:, i] / dz) + dt * (Do2 * (Co2n[:, i-1] - 2 * Co2n[:, i] + Co2n[:, i+1] / dz**2))) - dt * (((current_density * W * dz) / (4 * 96487)) / (H * W *dz))

            # Droplet temperatuer update D[:, 27] = Td

            # Boundary condition update
            T[:, 0] = Tin
            # T[:,0] = Tn[:,0] + dt * ((uin * Tin - u_a[:,0] * Tn[:,0]) / dz) - alphag[:,0] * ((Tin - 2*Tn[:,0] + Tn[:, 1]) / dz**2) * dt
            # T[:,-1] = Tn[:, -1] + dt * ((u_a[:, -2] * Tn[:, -2] - u_a[:, -1] * Tn[:, -1]) / dz) - alphag[:, -1] * (( Tn[:, -2] - 2 * Tn[:, -1] + Tn[:, -1]) / dz**2) * dt
            T[:, -1] = T[:, -2]

            # water vapour
            Cw[:,0] = Cwin
            # TODO Insert intended out commented lines which do not make any sense

            Cw[:, -1] = Cwn[:, -2]
            # Cw[:, -1] = Cwn[:, -1] + dt * ((u_a[:, -2] * Cwn[:, -2] - u_a[:, -1] * Cwn[:, -1] / dz) - dt * (Dh2o[:, -1] * ((Cwn[:, -2] - 2* Cwn[:, -1] + Cwn[:, -1]) / dz**2)) 

            # oxygen concentration
            # Co2n[:, 0] + dt * ((uin * oin - u_a[:, 0] * Co2n[:,0] / dz) - dt *(Do2 * ((Oin -2 * Co2n[:,0] + Co2n[:, 1]) / dz**2) - dt * (((current_density * W * dz) / (4 * 96487) / (H * W *dz)) # missing brackets
            Co2[:,0] = oin
            # Co2[:, 0] = Co2n[:, 0] + 2 dt * ((uin * Oin - u_a[:, 0] * Co2n[:, 0]) / dz) - dt * Do2 * (( Oin - 2 * Co2n[:, 0] + Co2n[:, 1]) / dz**2) - dt * (((current_density * W * dz) / (4 * 96487)) / (H * W * dz))
            
            Co2[:, -1] = Co2n[:, -2]
            # Co2[:, -1] = Co2n[:, -1] + dt * ((u_a[:, -2] * Co2n[:, -2] - u_a[:, -1] * Co2n[:, -1] / dz) - dt * Do2 * ((Co2n[:, -2] - 2 * Co2n[:, -1] + Co2n[:, -1]) / dz**2) - dt * (((current_density * W * dz) / (4 * 96487)) / (H * W *dz)

            # volume balance
            i = np.arange(1, Cz.size, step=1)
            V_air[:,0] = V_air[:,0] + V_air_in - V_air[:, 1] - dt * (((current_density * W * dz) / (4 * 96387)) * (8.314 * T[:, 0] / P[:, 0]))
            V_air[:, i] = V_air[:, i] + V_air[:, i-1] - V_air[:, i+1] - dt * (((current_density * W * dz) / (4 * 96487)) * (8.314 * T[:, i] / P[:, i]))

            V_air[:, -1] = V_air[:, -1] + V_air[:, -2] - V_air[:, -1] - dt * (((current_density * W * dz) / (4 * 96487 )) * (8.314 * T[:, -1] / P[:, -1]))

            uaprox = V_air / (H * W)
            #Ps = Psat(T)
            RH = (Cw * 8.314 * T) / Ps
        # end of the if condition Transport == 1

        ## Droplet Transport and update
        # Search all injection points to see if droplet is being injected into
        # amend July 2021: (only search those point that are within the search radius, i.e.
        # create new search radius function
        search_radius = 4 * D[:, 6] # i.e. search radius is 4 times the droplet radius

        timer_tic += 1
        timer_tic = 1   # NR: Häää???

        if timer_tic == 1:
            if D.shape[0] > 1:
                for i in np.arange(0, D.shape[0]):
                    # improve speed algorithm
                    # TODO insert commented lines

                    distance_compare = D[i, 6] + (poresize) > np.sqrt((D[i, 4] - zp)**2 + (D[i, 2] - xp)**2 + (D[i,3] - yp)**2)
                    D[i,1] = np.sum(distance_compare) * volume_flow_water_pore

                    if D[i,5] == 1:
                        D[i,1] = volume_flow_water_pore

                    if D[i,5] == 0:
                        D[i, 1] = volume_flow_water_pore
            timer_tic = -1
        # end of the if condition timer_tic == 1

        ## ========== Injection volume balance ========== 
        D[:, 0] = D[:, 0] + ((D[:,1] - D[:, 18]) * dt)

        if evaporation == 1:
            D[:, 26] = D[:, 26] + (h_conv * D[:, 13] * (Tbulk - D[:, 26])) / (rhow * D[:, 0] * Cpw) * dt - (hlv * (D[:, 18] * rhow) / (rhow * D[:, 0] * Cpw)) * dt

        # TODO insert commented if condition with comments

        # this filters all the rows with positive entries in the first column
        D = D[ D[:, 0] > 0, :] # new addition by OA from 17/06/2020, 'may help for evaporation'

        # NR: What does the following line mean???
        # %D(:,23)=D(:,23)*0; % remove this section for code to work

        # ========== Momentum Balance (Newton Second Law) ========== 
        dvdt = np.maximum((D[:, 22] - D[:, 23]) / (D[:, 0] * rhow), 0)  # acceleration update
        # move = D[: 19] > 0
        # dcdt[move] = ((D[move, 22], D[move, 23]) / (D[move, 1] * rhow)

        D[:, 19] = D[:, 19] + (dvdt * dt)                   # velocity update
        D[D[:, 19] < 0, 19] = 0                             # ensure no negative velocity
        # TODO insert commented lines


        D[(D[:, 5] == 7), 19] = uin
        D[D[:, 19] > uin, 19] = uin
        D[:, 4] = D[:, 4] + D[:, 19] * dt                     # z position update

        ## ========== After water injection, update classifications ==========
    
        # Hier wird eine Fall Unterscheidung gemacht

        # class 0 - ghost      
        # D[D[:, 5] == 0, 6] == Rpmax
        
        # class 1 - emerging
        mask1 = D[:, 0] > 0 
        mask2 = D[:, 5] == 0
        mask = mask1 & mask2
        D[mask, 5] == 1 # finds all droplets with volume and assigns as amerging (1)
        # D[D[:,0] > 0 & D[:, 5]== 0, 5] = 1             

        mask1 = D[:, 5] == 2
        mask2 = D[:, 24] < theta_GDL
        mask = mask1 & mask2  
        D[mask, 5] = 1
        # D[D[:,5] == 2 & D[:, 24] < theta_GDL, 5] = 1
        # D[D[:,5] == 1, 6] = Rpmax

        # class 2 - isolated
        mask1 = D[:, 5] == 1
        mask2 = D[:, 24] > theta_GDL
        mask = mask1 & mask2
        D[mask, 6] = 2       # finds all emerging droplets and assigns isolated if constant angle regime (2)

        # clas 3 - side wall (west)
        mask1 = D[:, 5] < 3
        mask2 = (D[:, 2] - D[:, 6]) < 0
        r3find_minus = mask1 & mask2
        #r3find_minus = np.where((D[:, 5] < 3) & ((D[:, 2] - D[:, 6]) < 0))[0] # finds droplets in class less than 3 and that have a radius that touches the left wall
        D[r3find_minus, 5] = 3                   # assign classification 3
        D[r3find_minus, 2] = 0                   # x position is left wall
        D[r3find_minus, 3] = 0                   # y position is zero

        # class 3 - side wall (east)
        mask1 = D[:, 5] < 3
        mask2 = (D[:, 2] + D[:, 6]) > W
        r3find_plus = mask1 & mask2
        #r3find_plus = np.where(D[:, 5] < 3 & (D[:, 2] + D[:, 6]) > W)   # finds droplets in class less than 3 and that have a radius that touches the right wall
        D[r3find_plus, 5] = 3                   # assign classfification 3
        D[r3find_plus, 2] = W                   # x position is left wall
        D[r3find_plus, 3] = 0                   # y position is zero


        # class 4 - corner droplet (north)
        mask1 = D[:, 5] < 4
        mask2 = D[:, 9] > H
        r4find = mask1 & mask2
        # r4find = np.where( D[:, 5] < 4 & D[:, 9] > H) # find droplets class less than 4 and also greater than the height of channel
        D[r4find, 5] = 4
        D[r4find, 3] = H
        D[r4find, 1] = 0

        # class 5 - truncated capillary bridge
        mask1 = D[:, 5] < 5
        mask2 = D[:, 11] > W
        r5find = mask1 & mask2
        # r5find = np.where( D[:, 5] < 5 & D[:, 11] > W) # find droplets class less than 5 and where Lc (length across width) is greater than the channel width
        D[r5find, 5] = 5


        if film_formation == 1:
            mask1 = D[:, 5] == 5
            mask2 = D[:, 9] > D[:, 25]
            r6find = mask1 & mask2
            # r6find = np.where(D[:,5] == 5 & D[:, 9] > D[:, 25])
            D[r6find, 5] = 6

        if plug_formation == 1:
            # class 6 plug flow
            mask1 = D[:, 5] == 6
            mask2 = D[:, 9] > 0.99 * H
            r7find = mask1 & mask2
            # r7find = np.where(D[:, 5] == 6 & D[:, 9] > (0.99 * H))
            D[r7find, 5] = 7
            mask1 = D[:, 5] == 7
            mask2 = D[:, 9] < H
            r5find2 = mask1 & mask2
            # r5find2 = np.where(D[:,5] == 7 & D[:, 9] < H)
            D[r5find2, 5] = 5

        water_removed = np.where(D[:, 4] - D[:, 12] > L)
        water_cond = np.sum(D[water_removed], 0) # new addition, track water removed
        newmatrix = np.where(D[:, 4] - D[:, 12] < L)
        D = D[newmatrix, :]

        # timing in matlab
        # timetest[ic, 1] = toc
        # tic

        ## ========== Coalescence Algorithm ========== 

        if coalescence == 1:
            timer_tic_coal += 1
            timer_tic_coal = 1       # Das Ergibt keinen Sinn

            if timer_tic_coal == 1:
                if D.shape[0] > 1:
                    
                    # each droplet tracks all other droplets, build upper triangular matrix, not to self duplicate connection
                    G = np.ones((D.shape[0], D.shape[0]))

                    # remove self-connection from upper triangular matrix and make sparse matrix
                    K = triu(G, k=1)  # k=1 removes the diagonal (self-connections)

                    # Convert the upper triangular sparse matrix to coordinate format (coo_matrix)
                    K_sparse = coo_matrix(K)

                    # Access the row (i) and column (j) indices of the non-zero elements in the sparse matrix
                    i, j = K_sparse.row, K_sparse.col

                    # Create connectivity matrix
                    ij = np.vstack((i, j)).T

                    # Generate an array of indices for all potential connections
                    k = np.arange(ij.shape[0])

                    # The next step would be to determine the distance between centers of all droplets, depending on what your `D` array represents.

                    s = np.sqrt((D[ij[k,0], 2] - D[ij[k,1], 2])**2 + (D[ij[k,0], 4] - D[ij[k,1], 4])**2)
                    # determine sum of radii
                    DR = D[ij[k,0], 6] + D[ij[k,1], 6]
                    # special case for film flow because R is very large
                    film_index = np.where(D[ij[k,0], 5] == 6 & D[ij[k,0], 6] == 6) # Why is here a double logical?
                    DR[film_index] = (D[ij[film_index, 0], 12] / 2) + (D[ij[film_index, 1], 12]/2)

                    # mark pair for coalescence
                    coal = np.argwhere(DR > s)      # find row coordination for coalescing pair
                    d1 = ij[coal, 0]                # colescing pair d1
                    d2 = ij[coal, 2]                # coaelescing pair d2
                    ik = 0

                    # for all colliding droplets

                    Dn = np.zeros((d1.shape[0], D.shape[1]))  # Preallocate Dn matrix

                    for ik in range(d1.shape[0]):
                        # Find the largest droplet between d1 and d2
                        max_index = np.argmax([D[d1[ik], 0], D[d2[ik], 0]])

                        # Create dvolumes matrix
                        dvolumes = np.array([
                            [D[d1[ik], 2], D[d2[ik], 2]],
                            [D[d1[ik], 3], D[d2[ik], 3]],
                            [D[d1[ik], 4], D[d2[ik], 4]]
                        ])

                        # Initialize new droplet as a copy of the first droplet
                        Dn[ik, :] = D[d1[ik], :]

                        # Volume balance
                        Dn[ik, 0] = D[d1[ik], 0] + D[d2[ik], 0]

                        # x, y, z position calculations
                        Dn[ik, 2] = (D[d1[ik], 0] * rhow * D[d1[ik], 2] + D[d2[ik], 0] * rhow * D[d2[ik], 2]) / (Dn[ik, 0] * rhow)
                        Dn[ik, 3] = (D[d1[ik], 0] * rhow * D[d1[ik], 3] + D[d2[ik], 0] * rhow * D[d2[ik], 3]) / (Dn[ik, 0] * rhow)
                        Dn[ik, 4] = (D[d1[ik], 0] * rhow * D[d1[ik], 4] + D[d2[ik], 0] * rhow * D[d2[ik], 4]) / (Dn[ik, 0] * rhow)

                        # Maximum values for specific columns
                        Dn[ik, 5] = max(D[d1[ik], 5], D[d2[ik], 5])  # Classification
                        Dn[ik, 6] = max(D[d1[ik], 6], D[d2[ik], 6])
                        Dn[ik, 7] = max(D[d1[ik], 7], D[d2[ik], 7])
                        Dn[ik, 8] = max(D[d1[ik], 8], D[d2[ik], 8])  # Changed 16/06/2020
                        Dn[ik, 9] = np.max(poresize)
                        
                        # Column 19 (matlab 20) balance (similar to x, y, z position)
                        Dn[ik, 19] = (D[d1[ik], 0] * rhow * D[d1[ik], 19] + D[d2[ik], 0] * rhow * D[d2[ik], 19]) / (Dn[ik, 0] * rhow)
                        
                        # Set height
                        Dn[ik, 10] = H / 2
                        
                        # Maximum value for column 27
                        Dn[ik, 26] = max(D[d1[ik], 26], D[d2[ik], 26])

                        # Kill connecting droplets and merge the new droplet
                        D[d1[ik], :] = np.zeros(D.shape[1])  # Droplet 1 = 0
                        D[d2[ik], :] = np.zeros(D.shape[1])  # Droplet 2 = 0

                        # Replace any NaN values with 0
                        D[np.isnan(D)] = 0

                        # Append the new droplet to D
                        D = np.vstack([D, Dn[ik, :]])

                    # D now contains the updated droplet data

                    if np.min(ik.shape) == 0:
                        del Dn
                        del volumes
                        del coal
                        del d1, d2
                        del s, G, K, ij, k
                        D = D[D[:,1] != 0, :]
                    timer_tic_coal = -1
            # emergence algorithm

            #timetest(ic,2)= toc
            #ic

            i_e = 0  # Initialize emergence index

            # Loop through all points in zp
            for i_p in range(zp.shape[0]):
                # Calculate the distance between the current point and all droplets
                s2 = np.sqrt((xp[i_p] - D[:, 2])**2 + (yp[i_p] - D[:, 3])**2 + (zp[i_p] - D[:, 4])**2)
                dr2 = D[:, 6] + poresize[i_p]
                
                # Set distances to 0 based on conditions
                s2[D[:, 5] == 1] = 0  # D(:,6)==1 in MATLAB becomes D[:,5]==1 in NumPy (0-indexing)
                s2[D[:, 5] == 0] = 0  # D(:,6)==0 in MATLAB becomes D[:,5]==0 in NumPy

                # Check which pores are covered
                pore_cover = s2 < dr2
                
                # If no droplets are covering the pore, add the droplet to D_emerge
                if np.sum(pore_cover) < 1:
                    if i_e == 0:
                        D_emerge = np.zeros_like(Dorig)  # Initialize D_emerge if not initialized
                    i_e += 1
                    D_emerge = np.vstack([D_emerge, Dorig[i_p, :]])

            # Append D_emerge to D if any droplets emerged
            if i_e > 0:
                D = np.vstack([D, D_emerge])

            # Check for divergence by comparing with maxdroplets
            if D.shape[0] > maxdroplets:
                print('Model has diverged, lower time step')
                break
            
            # Replace NaNs with 0
            D = np.nan_to_num(D)
            # TODO why is D an 3 dimensional Array now with only one entry in the first dimension?
            D = D[0, :, :]

            # Remove rows where all elements are zero
            # D = D[~np.all(D == 0, axis=1)]
            
            mask = np.any(D != 0, axis = 1)
            D = D[mask, :]

            # Record the time for this step
            # timetest[ic, 2] = time.time() - start_time

            # Apply temperature phase change correction for droplets with temperature < 273
            D[D[:, 26] < 273, 26] = Dorig[0, 26]  # MATLAB indexing is adjusted for 0-based indexing in NumPy

            # Clear intermediate variables
            if "D_emerge" in globals():
                del D_emerge
            del pore_cover, dr2, s2
            
            # create Droplet Classes
            D0 = np.where(D[:, 5] == 0)         # Find all droplets that use equation 0
            D1 = np.where(D[:, 5] == 1)         # Find all droplets that use stage 1
            D2 = np.where(D[:, 5] == 2)         # Find all droplets that use stage 2
            D3 = np.where(D[:, 5] == 3)         # Find all droplets that use stage 3
            D4 = np.where(D[:, 5] == 4)         # Find all droplets that use stage 4
            D5 = np.where(D[:, 5] == 5)         # Find all droplets that use stage 5
            D6 = np.where(D[:, 5] == 6)         # Find all droplets that use stage 6
            D7 = np.where(D[:, 5] == 7)         # Find all droplets that use stage 7

            # Separate Droplets based on class
            F0 = D[D0, :]
            F1 = D[D1, :]
            F2 = D[D2, :]
            F3 = D[D3, :]
            F4 = D[D4, :]
            F5 = D[D5, :]
            F6 = D[D6, :]
            F7 = D[D7, :]

            ## Update Droplet Geometric Parameters
            # F[:, 16] = np.pi * Rpmax**2       # Surface Area GDL coverage
            ## Emerging (constant contact Radius) (CCR)
            if F1.size != 0:
                # Equation V = np.pi * h / 6 * (3 *rc**2 + h**2) ... We know rc and V -> h and theta
                # height [:, 9]
                # radius of curvature [:, 7]
                # radius of pore [:, 8]
                # theta_app [:, 24]
                # droplet radius [:, 6] used for coalescence

                holde = F1[:, 9] # h
                h1 = F1[:, 9] # h
                hnew = holde + 10000

                f1_error = 1e-9
                F1_1 = lambda rpore, h, Vd: h**3 + 3 * (rpore**2) * h - 6 * Vd / np.pi # Equation stage 1
                F1_J = lambda rpore, h: 3*h**2 + 3 * rpore**2 # Jacobian for stage 1
                F1_R = lambda Vd, h: (((3*Vd) / (np.pi * h**2)) + h) / 3 # radius as a function of height
                F1_Ac = lambda R, theta_1c: (R**2 / 2) * (theta_1c - np.sin(theta_1c)) # cross section area
                F1_As = lambda R, h: 2 * np.pi * R * h              # surface area
                F1_Ac_90 = lambda R, theta_app: np.pi * R**2 - ((R**2 / 2) * ((2 * (np.pi - theta_app)) - np.sin(2 * (np.pi - theta_app))))

                while np.sqrt(np.abs(np.transpose(holde-hnew) * (holde-hnew))) > 1e-9:
                    # solving the height equation using the newton method
                    hnew = h1 - np.linalg.solve(np.eye(len(D1)) @ F1_J(F1[:, 8], h1), F1_1(F1[:, 8], h1, F1[:, 0]))
                    holde = h1
                    h1 = hnew
                # end of while loop

                F1[:, 9] = hnew                     # update height
                F1[:, 7] = F1_R(F1[:, 0], F[:, 9])  # update radius of curvature


                # Find the indices where the radius is greater than the height
                b90 = np.where(F1[:, 9] - F1[:, 7] < 0)[0]  # Indexing adjusted for 0-based index (MATLAB uses 1-based indexing)
                a90 = np.where(F1[:, 9] - F1[:, 7] > 0)[0]

                # Compute the apparent contact angle
                F1[b90, 24] = np.real(np.arcsin(F1[b90, 8] / F1[b90, 7]))

                # Compute the cross-section area using the function F1_Ac
                F1[b90, 14] = F1_Ac(F1[b90, 7], 2 * F1[b90, 24])

                # Compute the surface area using the function F1_As
                F1[b90, 13] = F1_As(F1[b90, 7], F1[b90, 9])

                # Set the radius of curvature equal to rpore for plotting
                F1[b90, 6] = F1[b90, 8]

                F1[a90, 24] = np.pi - np.real(np.arcsin(F1[a90, 8] / F1[a90, 7]))  # apparent contact angle
                F1[a90, 14] = F1_Ac_90(F1[a90, 7], F1[a90, 24])                    # cross section area
                F1[a90, 13] = F1_As(F1[a90, 7], F1[a90, 9])                        # surface area
                F1[a90, 6] = F1[a90, 7]                                            # radius of curvature used for plotting

                r2find = np.where(F1[:, 24] > theta_GDL)[0]                        # find indices where F1(:,25) > theta_GDL
                F1[r2find, 5] = 2                                                  # set F1(:,6) to 2 at these indices
                # Outcommented line
                # F[:, 3] = F1[:, 9] - F1[:, 6]
                # Update y position (new method, updating radius)
                F1[:, 3] = F1[:, 9] - F1[:, 7]                                     

                # Medial axis y location update
                F1[:, 10] = (H - F1[:, 9]) / 2                                    

                # Length of droplet in z direction
                F1[:, 12] = 2 * F1[:, 6]                                          

                # Surface area GDL coverage
                F1[:, 16] = np.pi * (F1[:, 6])**2                                 

                # Adhesion length
                F1[:, 15] = F1[:, 8]                                              

                # Liquid-solid shear area
                F1[:, 17] = np.pi * F1[:, 8]**2                                    
                # end of the loop

            if F1.size != 0:
                # Analytical Geometry equation for R as function of V
                # (V = 4/3 * np.pi * R**3 + (2 + np.cos(theta)(1-np.cos(theta))**2
                # Finding Radius of Curvature

                # New radius of curvature found by analytical geometric equation 2
                Rnew = (F2[:, 0] / ((4/3) * np.pi - (np.pi/3) * ((2 + np.cos(np.pi - theta_GDL)) * ((1 - np.cos(np.pi - theta_GDL))**2))))**(1/3)

                F2[:, 7] = Rnew                                 # Radius of Curvature
                F2[:, 8] = F2[:, 7] * np.sin(np.pi - theta_GDL) # Rpore (also known as Rcontact)
                F2[:, 9] = F2[:, 7] * (1 - np.cos(theta_GDL))   # Height of droplet
                F2[:, 14] = F1_Ac_90(F2[:, 7], F2[:, 24])       # Cross-sectional area of droplet
                F2[:, 13] = F1_As(F2[:, 7], F2[:, 9])           # Surface area
                F2[:, 6] = F2[:, 7]                             # R for plotting and coalescence
                F2[:, 3] = F2[:, 9] - F2[:, 6]                  # y position
                F2[:, 10] = (H - F2[:, 9]) / 2                  # Medial axis b
                F2[:, 12] = 2 * F2[:, 6]                        # Length of droplet along channel z
                F2[:, 16] = (F2[:, 8])**2 * np.pi               # Area coverage of droplet (need to update this)
                F2[:, 15] = F2[:, 8]                            # Adhesion length
                F2[:, 17] = np.pi * F2[:, 8]**2                 # Liquid-solid shear area

                # Conditional update based on theta_GDL
                if theta_GDL < (np.pi / 2):
                    F2[:, 6] = F2[:, 8]
            # end of if condition

            if "Rnew" in globals():
                del Rnew

            ## Side Wall
            if F3.size != 0:
                alpha_wall = (np.pi / 2) - theta_wall
                alpha_GDL = theta_GDL - np.pi / 2
                theta_c = 2 * (np.pi - theta_GDL)

                if theta_GDL < np.pi / 2:
                    theta_c = theta_GDL

                theta_c_total = theta_wall + alpha_GDL
                a3 = 2 - (3 * np.sin(alpha_wall)) + np.sin(alpha_wall)**3
                #a3 = 2 - (3 * np.sin(theta_wall)) + np.sin(theta_wall)**3
                #b3 = 1 - (theta_c - np.sin(theta_c))/2
                b3 = 1 - (theta_c - np.sin(theta_c)) / (2 * np.pi)
                #c3 = a3 * b3
                c3 = np.pi * (a3 @ b3)

                # New Radius of curvature
                Rnew3 = (3 * F3[:,0] / c3)**(1/3)

                x31 = Rnew3 * np.sin(alpha_wall)
                x32 = Rnew3 - x31
                x33 = Rnew3 * np.sin(alpha_GDL)
                x34 = (Rnew3**2 - x31**2)**(0.5)
                x35 = (((x33 / np.tan(alpha_GDL)) - x31) / (1 - np.cos(theta_wall)))

                # Define anonymous functions
                F3_Ac = lambda R, theta_c_total, x5, x3, x4: ((R**2) / 2) * (theta_c_total - np.sin(theta_c_total)) + (x5 * (x3 + x4) / 2)
                F3_As = lambda R, theta_c, x2: 2 * np.pi * R * x2 * (1 - (theta_c - np.sin(theta_c)) / 2)

                # Update matrix F3 using the functions and other operations
                F3[:, 7] = Rnew3                                         # Radius of Curvature
                F3[:, 8] = x34                                           # Rpore (also known as Rcontact)
                F3[:, 15] = x34 + x33 + (F3[:, 7] - x31)                 # Adhesion length
                F3[:, 9] = x33 + x34                                     # Height of droplet
                F3[:, 14] = F3_Ac(Rnew3, theta_c_total, x35, x33, x34)   # Cross-sectional area of droplet
                F3[:, 13] = F3_As(Rnew3, theta_c, x32)                   # Surface area
                F3[:, 6] = F3[:, 7]                                      # R for coalescence and plotting
                F3[:, 10] = (W - x32) / 2                                # Medial axis b
                F3[:, 12] = F3[:, 7]                                     # Length along channel z

                # Ratio of contact angles to determine the bottom length
                F3[:, 16] = ((((x33 / np.tan(alpha_GDL)) - x31) / (1 - np.cos(theta_wall)))**2 / 2) * (2 * theta_wall - np.sin(2 * theta_wall))
                F3[:, 17] = np.pi * (x34**2)                             # Liquid solid shear area

                # Determine x and y location of droplets
                x3min = np.where(F3[:, 2] < (W / 2))[0]                  # For droplets at xmin wall
                x3max = np.where(F3[:, 2] > (W / 2))[0]                  # For droplets at xmax wall

                # Update x positions based on the conditions
                F3[x3min, 2] = xmin - x31[x3min]
                F3[x3max, 2] = xmax + x31[x3max]

                # Update y positions
                y3min = np.where(F3[:, 3] == ymin)[0]                   # Floating Point Ar. BUG possible due to ==
                F3[y3min, 3] = ymin + x33[y3min]
                F3[:, 3] = F3[:, 9] - x34                                # Update y position

            # end of if condition

            ## Corner
            if F4.size != 0:
                alpha_4 = 0.5 * np.pi - theta_wall
                alpha_42 = theta_wall - 0.25 * np.pi                    # used to calculate cross section
                a4 = 2 - 3 * np.sin(alpha_4) + np.sin(alpha_4)**3
                b4 = 2 * theta_wall - np.sin(2 * theta_wall)
                c4 = 2 * b4

                # Define anonymous functions
                F4_As = lambda R, alpha, theta_wall: 2 * np.pi * (R**2) * (1 - np.sin(alpha)) * (2 * theta_wall - np.sin(2 * theta_wall)) # surface area
                # need to fix this cross section equation
                # NR: WHY??? Where's the error? Why does it say on the next line to be fixed?
                F4_Ac = lambda R, x42: (x42**2) / 2 + (R**2) / 2 * (2 * alpha_42 - np.sin(2 * alpha_42))

                # Calculations
                Rnew4 = (6 * F4[:, 0] / (a4 * b4))**(1/3)
                x41 = Rnew4 * np.sin(alpha_4)
                x42 = Rnew4 - x41

                # Update F4 array
                F4[:, 7] = Rnew4                                       # Radius of Curvature
                F4[:, 8] = F4[:, 7] * (1 - np.sin(alpha_4))            # Rpore (also known as Rcontact)
                F4[:, 14] = F4_Ac(F4[:, 7], x42)                       # Cross section area
                F4[:, 13] = F4_As(F4[:, 7], alpha_4, theta_wall)       # Surface area
                F4[:, 6] = Rnew4                                       # R for coalescence and plotting
                F4[:, 11] = x42                                        # Length in x direction Lc
                F4[:, 10] = (W - x42) / 2                              # Medial axis b (updated)
                F4[:, 12] = 2 * Rnew4 * np.sin(theta_wall)             # Length in channel z direction
                F4[:, 16] = 0                                          # No GDL surface coverage for corner droplets
                F4[:, 3] = ymax + x41                                  # New y position
                F4[:, 15] = x42 * 2                                    # Adhesion length
                F4[:, 17] = 2 * np.pi * (F4[:, 7]**2) * (theta_wall / np.pi)  # May need to change

                # Logical indexing for x-position update
                x4min = np.where(F4[:, 2] < (W / 2))[0]
                x4max = np.where(F4[:, 2] > (W / 2))[0]

                # Update x positions
                F4[x4min, 2] = xmin - x41[x4min]
                F4[x4max, 2] = xmax + x41[x4max]
            
            # end of if condition

            ## Capillary bridge

            if F5.size != 0:

                del z51
                # Calculations
                V5 = F5[:, 0]
                z51 = (W / 2) / np.cos(theta_wall) * (1 - np.cos((np.pi / 2) - theta_wall))
                R51 = (W / 2) / np.cos(theta_wall)  # Equivalent to sec(theta_wall)
                theta_51 = np.pi - (2 * theta_wall)
                theta_52 = theta_wall * 2

                phi_51 = (R51**2) * (theta_51 - np.sin(theta_51)) / (2 * z51 * W)
                phi_52 = (theta_52 - np.sin(theta_52)) / (2 * np.pi)

                a5 = phi_52 * W * np.pi
                b5 = -2 * phi_52 * W * np.pi * phi_51 * z51
                c5 = phi_52 * W * np.pi * phi_51 * (z51**2) - V5

                R52 = (b5 + np.sqrt(b5**2 - 4 * a5 * c5)) / (2 * a5)
                R53 = R52 - z51
                h5 = R52 * (1 - np.cos(theta_wall))
                h53 = R53 * (1 - np.cos(theta_wall))

                # New section
                h5b = h53 - H
                h5b[h5b < 0] = 0
                c5b = 2 * R52 * np.sin(np.arccos(1 - (h5b / R52)))
                c5b[c5b < 0] = 0

                # Update F5 array
                F5[:, 7] = R52                            # Radius of Curvature
                F5[:, 6] = R52                            # R for plotting/coalescence
                F5[:, 12] = 2 * R52 * np.sin(theta_wall)  # Length along channel LD
                F5[:, 14] = R53 * (1 - np.cos(theta_wall)) * W  # Cross-sectional area
                F5[:, 14] = h5 * W                        # May affect 25/06/2020
                F5[F5[:, 14] > H * W, 14] = H * W
                F5[:, 13] = np.pi * R52 * phi_52 * W      # Surface area
                F5[:, 10] = (H - h5) / 2                  # Medial axis b
                F5[:, 16] = 0                             # GDL surface coverage
                F5[:, 16] = c5b * W                       # New estimation 11/06/2020
                F5[:, 2] = W / 2                          # x position
                F5[:, 3] = H + (R52 - h5)                 # y position
                F5[:, 8] = W + 2 * (R52 * (1 - np.cos(theta_wall)))  # Rpore
                F5[:, 15] = W                             # LAD
                F5[:, 17] = F5[:, 12] * W                 # Liquid solid area
                F5[:, 11] = h5
                F5[:, 9] = h53                            # New droplet height
                # end of if condition
            
            if "h5" in globals():
                del h5

            ## Plug # empty for now # NR: What does it mean? Is the code dead?
            if F6.size != 0:
                # New addition 16/07/2020 Inertial, viscous, capillary balance relvelocity
                rel_v6 = uin - F6[:, 19]  # Adjusted for 0-based indexing

                Ca6 = rel_v6 * mug[0] / sigma
                A6 = (rhoa * (rel_v6**2)) / (4 * sigma * (1 - np.cos(theta_wall)))

                # Coefficients for equation
                a6 = A6
                b6 = 1 + (3 * Ca6 * H / (1 - np.cos(theta_wall))) - (2 * A6 * H)
                c6 = -(3 * Ca6 * H / (1 - np.cos(theta_wall))) - (2 * H)
                d6 = H**2

                # Anonymous functions for f6 and df6
                f6 = lambda h6: a6 * (h6**3) + b6 * (h6**2) + (c6 * h6) + d6
                df6 = lambda h6: 3 * a6 * (h6**2) + 2 * b6 * h6 + c6

                # Initializations
                h6 = np.zeros_like(uin)
                holde6 = h6.copy()
                h16 = h6.copy()
                hnew6 = holde6 + 10000

                f1_error6 = 1e-9
                iter6 = 0
                w6 = 0.7

                # Newton-Raphson iterations to solve for h6
                while np.sqrt(np.abs(np.dot(holde6 - hnew6, holde6 - hnew6.T))) > 1e-12:
                    hnew6 = h16 - np.linalg.solve(w6 * f6(h16) / df6(h16))
                    holde6 = h16
                    h16 = hnew6
                    iter6 += 1

                # Post-processing: ensuring h16 remains within bounds
                h16 = np.clip(h16, 0, H)

                # Assign maximum height before deformation in capillary bridge
                F6[:, 25] = h16  # Adjusted for 0-based indexing
                            
                # Start of geometry section

                F6[:, 14] = H * W  # Cross-section area, equivalent to F6(:,15) in MATLAB (adjusted for 0-based indexing)
                h6 = F6[:, 25]     # Equivalent to F6(:,26), as MATLAB is 1-based

                theta_61 = np.pi - (2 * theta_wall)
                theta_62 = 2 * (theta_wall - (np.pi / 4))
                R61 = (W / 2) * (1 / np.cos(theta_61))
                x61 = np.tan(np.pi / 4) * h6
                x62 = x61 / np.tan(theta_62)
                R62 = F6[:, 7]     # Equivalent to F6(:,8)
                y6 = H + x62

                # Calculate L6 using the given equation
                L6 = (F6[:, 0] - (2 * h6 * x61 * W) - ((R62**2) * (theta_62 - np.sin(theta_62)) * W)) / \
                    ((h6 * W) - ((R61**2) / 2) * (theta_61 - np.sin(theta_61)))

                L6[L6 < 0] = 0  # Replace any negative values with 0
                L6 += 2 * R62

                R6plot = (L6**2) / (8 * h6) + (h6 / 2)
                y6plot = H + R6plot - h6

                # Update F6 array with calculated values
                F6[:, 7] = R62         # Radius of Curvature
                F6[:, 6] = R6plot      # R for plotting/coalescence
                F6[:, 12] = L6         # Length along channel LD

                # Cross-section area
                F6[:, 14] = np.minimum(h6 * W, H * W)  # Ensure values do not exceed H * W

                # Surface area
                F6[:, 13] = W * L6

                # Medial axis b (height)
                F6[:, 10] = (H - h6) / 2

                # x position
                F6[:, 2] = W / 2

                # y position
                F6[:, 3] = y6plot

                # Rpore
                F6[:, 8] = W + 2 * (R62 * (1 - np.cos(theta_wall)))

                # LAD
                F6[:, 15] = F6[:, 8]

                # Liquid-solid area
                F6[:, 17] = (L6 * W) + (2 * h6 * L6)

                # Other assignments
                F6[:, 11] = W
                F6[:, 9] = h6
                
                if "h6" in globals():
                    del h6
                if "R6plot" in globals():
                    del R6plot
                if  "y6plot" in globals():
                    del y6plot
                
                if "L6" in globals():
                    del L6

                if "R62" in globals():
                    del R62
                # end of if condition

            if F7.size != 0:
                # is the same algorithm as the capillary bridge rgime
                # NR: Noted for the refactoring, but the lines below changed too.

                V7 = F7[:, 0]
                # Start of calculations
                z71 = ((W / 2) / np.cos(theta_wall)) * (1 - np.cos((np.pi / 2) - theta_wall))
                R71 = (W / 2) * (1 / np.cos(theta_wall))
                theta_71 = np.pi - (2 * theta_wall)
                theta_72 = 2 * theta_wall

                phi_71 = (R71**2) * (theta_71 - np.sin(theta_71)) / (2 * z71 * W)
                phi_72 = (theta_72 - np.sin(theta_72)) / (2 * np.pi)

                a7 = phi_72 * W * np.pi
                b7 = -2 * phi_72 * W * np.pi * phi_71 * z71
                c7 = phi_72 * W * np.pi * phi_71 * (z71**2) - F7[:, 0]

                # Solve the quadratic equation for R72
                R72 = (b7 + np.sqrt(b7**2 - 4 * a7 * c7)) / (2 * a7)
                R73 = R72 - z71
                h7 = R72 * (1 - np.cos(theta_wall))
                h73 = R73 * (1 - np.cos(theta_wall))

                # New section added on 11/06/2020
                h7b = h73 - H
                h7b[h7b < 0] = 0

                c7b = 2 * R72 * np.sin(np.arccos(1 - (h7b / R72)))
                if theta_wall > (np.pi / 2):
                    c7b = 2 * R73 * np.sin(np.arccos(1 - (h7b / R73)))
                c7b[c7b < 0] = 0

                # New derivation October 2020
                L7 = F7[:, 0] / (H * W)

                # Update F7 array with the calculated values
                F7[:, 7] = R72  # Radius of Curvature for plotting
                F7[:, 6] = R72  # Radius for plotting/coalescence
                F7[:, 12] = 2 * R72 * np.sin(theta_wall)  # Length along channel LD

                # Cross-sectional area
                F7[:, 14] = np.minimum(h73 * W, H * W)  # Ensures values do not exceed H * W

                # Surface area
                F7[:, 13] = np.pi * R72 * phi_72 * W

                # Medial axis b
                F7[:, 10] = (H - h7) / 2
                F7[:, 16] = 0       # GDL surface Coverage
                # GDL surface coverage
                F7[:, 7] = L7 / 2  # Radius of Curvature
                F7[:, 6] = L7 / 2  # Radius for plotting/coalescence
                F7[:, 12] = L7  # Length along channel LD

                # x and y positions
                F7[:, 2] = W / 2  # x position
                F7[:, 3] = H / 2  # y position

                # Rpore and LAD
                F7[:, 8] = W + 2 * (R72 * (1 - np.cos(theta_wall)))  # Rpore
                F7[:, 15] = F7[:, 8]  # LAD

                # Liquid-solid area
                F7[:, 17] = F7[:, 12] * W  # Liquid-solid area

                # New droplet height
                F7[:, 11] = h7
                F7[:, 9] = h73  # New droplet height


                #end of if condition
            # Re-establish droplet matrix

            if F1.size == 0:
                F1 = np.zeros((1, D.shape[1]))

            if F2.size == 0:
                F2 = np.zeros((1, D.shape[1]))

            if F3.size == 0:
                F3 = np.zeros((1, D.shape[1]))

            if F4.size == 0:
                F4 = np.zeros((1, D.shape[1]))

            if F5.size == 0:
                F5 = np.zeros((1, D.shape[1]))

            if F6.size == 0:
                F6 = np.zeros((1, D.shape[1]))

            if F7.size == 0:
                F7 = np.zeros((1, D.shape[1]))

            if F0.size == 0:
                F0 = np.zeros((1, D.shape[1]))
            else:
                F0 = F0[0, :, :]  # BUG FIX, since the shape of 0 becomes 3 dimensional

            # Combine all matrices
            D = np.vstack([F1, F2, F3, F4, F5, F6, F7, F0])

            # Remove rows where all elements are zero
            D = D[np.any(D, axis=1)]

            # Clear variables (optional, though usually not needed in Python)
            del F1, F2, F3, F4, F5, F6, F7, F0

            # ========== Force Balance ========== 
            # First finde projected cross sectional area defined from classes
            # Update values where condition D(:,15) > (H*W)
            D[D[:, 14] > (H * W), 14] = (H * W)  

            # Assign P_A with column 14
            P_A = D[:, 14]

            # Set the air velocity in column 21 
            D[:, 21] = u_a[0]  # Assuming u_a is a single value or first entry is used

            # Update D(:,20) where it's greater than D(:,22)
            D[D[:, 19] > D[:, 21], 19] = D[D[:, 19] > D[:, 21], 21]

            # Use mass balance (continuity equation) to find velocity of air
            U_I = (D[:, 21] * H * W) / ((H * W) - P_A)

            # Determine pressure force from Bernoulli's equation (pressure to kinetic)
            DP_p = (rhoa / 2) * (((U_I - D[:, 19]) ** 2) - (D[:, 21] ** 2))

            # Pressure Force = Pressure Difference * Cross section area
            FPressure = DP_p * P_A
            FPressure = np.nan_to_num(FPressure)  # Remove NaNs, replacing with 0

            # Determine Shear Force gas-liquid on each droplet (parallel plate with b)
            DP_v = (3 * mug[0] * (U_I - D[:, 19])) / D[:, 10]
            FShear = DP_v * D[:, 13]
            FShear = np.nan_to_num(FShear)  # Remove NaNs, replacing with 0

            # Shear Force liquid-solid
            FShear_ls = ((muw * D[:, 19]) / (D[:, 6] / 2)) * D[:, 17]

            # Special case for condition D(:,6) == 6
            mask = D[:, 5] == 6
            FShear_ls[mask] = (6 * (muw * D[mask, 19]) / D[mask, 9]) * D[mask, 17]
            FShear_ls = np.nan_to_num(FShear_ls)  # Remove NaNs, replacing with 0

            # Total Drag Force
            D[:, 22] = FShear + FPressure - FShear_ls

            # Set any negative total drag forces to 0
            D[D[:, 22] < 0, 22] = 0

            # for plug flow
            mask1 = D[:, 5] == 7
            mask2 = D[:, 0] > 0
            index = mask1 & mask2
            #index = np.where(D[:, 5] == 7 & D[:, 0] > 0)
            Fshearplug = (( 6 * muw * uin ) / H) * (2 * ( H * D[index, 12]) + 2 * ( W * D[index, 12])) # N
            DPplug = Fshearplug / (H * W)
            DP3 = DPplug # this is not used anymore afterwards
            #DP3[index, 0] = DPplug

            # new oct 2020
            bp = W / 2
            cp = H / 2
            betap = (( 1/3) - (64 / (np.pi**5)) * ((cp / bp) * np.tanh(np.pi * bp / (2 * cp))))
            DPplug = (uin * muw * D[index, 12]) / (cp**2 * betap)
            DPplug = np.nan_to_num(DPplug, nan = 0.0, posinf = 0.0, neginf = 0.0)

            # Adhesion Force based on contact line normal to flow and contact angle

            # Adhesion Force based on contact line normal to flow and contact angle

            ADV = np.ones((D.shape[0],)) * theta_GDL

            ADV[D[:, 5] == 1] = theta_GDL     # emerging
            ADV[D[:, 5] == 2] = theta_GDL     # isolated
            ADV[D[:, 5] == 3] = theta_wall    # side wall
            ADV[D[:, 5] == 4] = theta_wall    # corner droplet
            ADV[D[:, 5] == 5] = theta_wall    # truncated capillary bridge
            ADV[D[:, 5] == 6] = theta_wall    # plug
            ADV[D[:, 5] == 7] = theta_wall    # plug

            # Receding angle based on constant CAH (contact angle hysteresis)
            REC = ADV - theta_hyst

            REC[D[:, 5] == 5] = theta_wall - 1  # hysteresis when not attached to GDL
            REC[D[:, 5] == 6] = theta_wall - 1
            REC[D[:, 5] == 7] = theta_wall - 1
            REC[D[:, 5] == 3] = theta_wall - 1

            # Adhesion force calculation
            D[:, 23] = sigma * D[:, 15] * (np.cos(REC) - np.cos(ADV))  # new adhesion

            # Special case for plug
            D[D[:, 5] == 6, 23] = sigma * D[D[:, 5] == 6, 15] * (np.cos(REC[D[:, 5] == 6]) - np.cos(ADV[D[:, 5] == 6]))

            # Record Pressure Drop from Droplets Alone
            DP = DP_v + DP_p
            DP[D[:, 5] == 7] = DPplug


            ## Determine Pressure Drop in Channel

            # single phase pressure Drop
            bp = W / 2
            cp = H / 2
            betap = ((1/3) - (64 / (np.pi**5)) * (cp / bp) * np.tanh((np.pi * bp) / (2 * cp)))
            DPsinglephase = (uin * mug * dz) / (cp**2 * betap)
            DP = np.nan_to_num(DP, nan=0.0, posinf=0.0, neginf=0.0)
            DP = np.abs(DP)

            TotalP = np.sum(DPsinglephase) * np.sum(DP)

            Rel_perm = uin * mug[0] / TotalP
            perm = uin * mug[0] / np.sum(DPsinglephase)

            # Dec 2020 - Determine pressure drop along the channel
            # commented lines 1296 - 1315

            # ==========  Evaporation ========== 

            if evaporation ==1:

                # Find where the location is less than the channel points
                Location_Evap = D[:, 4] > Cz  

                # Find the first True value in each row (i_evap)
                i_evap = np.argmax(Location_Evap, axis=1)

                # Get Cbulk and Tbulk
                Cbulk = Cw[i_evap]  # Converting to column array for droplet evaporation assessment
                Tbulk = T[i_evap]

                # Surface is at saturated conditions of temperature of the droplet
                Psurf = Psat(D[:, 26])  # Update Td every timestep

                # Kelvin equation
                Psurf = np.exp(2 * sigma * (18 / 1000) / (D[:, 6] * 8.314 * D[:, 26])) * Psurf

                # Concentration at surface in mol/m³ by the ideal gas law
                Csurf = Psurf / (8.314 * D[:, 26])

                # Update physical properties for boundary layer properties
                mean_temp = (D[:, 26] + Tbulk) / 2

                kgS = (1.97 * 10**-4) * (mean_temp**0.858)

                mugS = (1.691 * 10**-5) + (4.984 * 10**-8 * (mean_temp - 273)) - \
                    (3.187 * 10**-11 * (mean_temp - 273)**2) + \
                    (1.319 * 10**-14 * (mean_temp - 273)**3)

                Dh2oS = 21.6 * 10**-6 * (1 + (0.0071 * (mean_temp - 273)))

                # Dh2oS = Dh2o(1) -- Commented out as in the original code

                # Sherwood number calculation
                Sh = 2 + 0.6 * (mugS / (rhoa * Dh2oS))**(1/3) * ((2 * D[:, 6] * (U_I - D[:, 19]) * rhoa) / mugS)**(1/2)

                # Set Sh to 2 where D[:,5] == 1
                Sh[D[:, 5] == 1] = 2

                # Mass transfer coefficient
                km = Sh * Dh2oS / (2 * D[:, 6])

                # Set km to 0 where V_evap is inf or NaN
                km[np.isinf(V_evap) | np.isnan(V_evap)] = 0

                # Volume evaporation rate
                V_evap = 0.018 * km * D[:, 13] * (Csurf - Cbulk) / rhow

                # Set V_evap to 0 where it is inf or NaN
                V_evap[np.isinf(V_evap) | np.isnan(V_evap)] = 0

                # Assign the evaporation rate to column 19 of D
                D[:, 18] = V_evap

                # Latent heat of vaporization/condensation
                hlv = ((-2.5464 * D[:, 26]) + 2503.6) * 1000

                # Temperature difference factor
                BT = (Cpg[0] * (Tbulk - D[:, 26])) / hlv

                # Nusselt number calculation
                Nu = (2 + 0.6 * (Cpg[0] * mugS / kgS)**(1/3) * ((2 * D[:, 6] * (U_I - D[:, 19]) * rhoa) / mugS)**(1/2)) / (1 + BT)**0.7

                # Set Nu to 2 where D[:,5] == 1
                Nu[D[:, 5] == 1] = 2

                # Convective heat transfer coefficient
                h_conv = Nu * kgS / (2 * D[:, 6])

                # Set h_conv to 0 where it is inf or NaN
                h_conv[np.isinf(h_conv) | np.isnan(h_conv)] = 0
                
                # Find which segments the droplet is adding mass to
                w_c.fill(0)  # Reset w_c array to zeros but it is alredy 0 ???

                for ie in range(len(i_evap)):
                    Cw[i_evap[ie]] = Cw[i_evap[ie]] + ((V_evap[ie] * rhow / 0.018) * dt) / (dz * H * W)
                    # T[i_evap[ie]] = T[i_evap[ie]] - (h_conv[ie] * D[ie, 14] * (Tbulk[ie] - D[ie, 27]) * dt)
                    # +(hlv[ie] * (D[ie, 19] * rhow) / (rhog[i_evap[ie]] * Cpg[i_evap[ie]]) * (H * W * dz)) * dt

                    w_c[i_evap[ie]] = w_c[i_evap[ie]] + D[ie, 18]

                # Sw_c = w_c / (dz * W)

                # Update physical properties
                Dh2o = 21.6e-6 * (1 + (0.0071 * (T - 273)))  # Diffusion coefficient of water vapor
                Cpg = (1.00926e3) - (4.0403e-2 * (T - 273)) + (6.1759e-4 * (T - 273) ** 2) - \
                        (4.097e-7 * (T - 273) ** 3)  # Heat capacity of air
                rhog = Pin / (287.058 * T)  # Density of air
                kg = (1.97e-4) * (T ** 0.858)  # Thermal conductivity of air
                alphag = kg / (Cpg * rhog)  # Thermal diffusivity of air
                mug = (1.691e-5) + (4.984e-8 * (T - 273)) - (3.187e-11 * (T - 273) ** 2) + (1.319e-14 * (T - 273) ** 3)

                hlv = ((-2.5464 * D[:, 26]) + 2503.6) * 1000  # Latent heat of vaporization/condensation

                AverageRH[ic, :] = RH  # Update average relative humidity
            # end of if condition of evaporation == 1

            # Initialise simulation with droplets already present
            # here are some lines commented out in the original code
            RH = ( Cw * 8.314 * T ) / Ps

            # ========== Condensation (in PROGRESS) ========== 
            if condensation == 1:
            # Calculation of Cc and Vcond
                Cc = (RH - 1) * Psat(T) / (8.314 * T)
                Vcond = ((RH - 1) * Psat(T) * 0.018) / (8.314 * T * rhow)
                Vcond[Vcond < 0] = 0

                if np.any(Vcond > 0): # this can be a BUG, since the condition in matlab (1426) is ambigious
                    icond = np.arange(Cz.shape[1])
                    Di = np.zeros_like(Dorig)  # Initialize Di with the same shape as Dorig
                    
                    Di[icond, :] = Dorig[0, :]  # Copy the first row of Dorig into Di
                    Di[icond, 0] = Vcond[icond]  # Assign Vcond values to the second column (index 1)
                    Di[icond, 2] = np.round(np.random.rand(len(icond)) * W, 2)  # Random values in the third column (index 3)
                    Di[icond, 3] = H  # Assign height (H) to the fourth column (index 4)
                    Di[icond, 4] = Cz[icond]  # Assign Cz values to the fifth column (index 5)
                    
                    # Calculate the seventh column (index 7)
                    Di[icond, 6] = ((3 * Di[icond, 1]) / (4 * np.pi)) ** (1/3)
                    
                    Di[:, 5] = 4  # Set sixth column (index 6) to 4
                    D = np.vstack([D, Di])  # Append Di to D
                # timetest(ic, 4) = toc

            # ========== Saving Data ========== 

            # Store values in S matrix. BUG WARNING: Is S already declared?
            # BUG ic is not working
            S = np.zeros((ic, 4))
            S[ic-1, 0] = t
            S[ic-1, 1] = TotalP
            S[ic-1, 2] = np.sum(V_evap)
            S[ic-1, 3] = np.mean(D[D[:, 0] > 0, 19]) # index shift right?

            # Store values in Time matrix
            Time = np.zeros((ic, 4))
            Time[ic-1, 0] = 0 # toc() elapsed time measuring
            Time[ic-1, 1] = D.shape[0]
            Time[ic-1, 2] = np.sum(D[:, 0]) / (L * W * H)
            Time[ic-1, 3] = np.sum(D[:, 16]) / (L * W)

            # Save condensation data
            Water = np.zeros((ic, 1))
            Water = water_cond

            E_data = np.zeros((ic, 2))
            # Save evaporation data
            E_data[ic-1, 0] = D[0, 0]
            E_data[ic-1, 1] = D[0, 26]


            # Save relative permeability estimation
            permsave = np.zeros((ic, 1))
            permsave[ic-1, 0] = Rel_perm / perm

            # Increment visualisation droplet counter
            visualiseDropletCounter += 1

            # Check if visualisation should be triggered
            if visualiseDropletCounter > vdclimit:
                #visualiseDroplets4()  # Call to the visualisation function. BUG: Unknown function
                visualiseDropletCounter = 0
            
            # ========== Figure Plotting During Simulation ========== 


            if timeplot == 1:
                # Create a figure for saturation
                plt.figure(1)
                plt.clf()  # Clear the current figure

                # Plot Channel Saturation (Time[:,3]) and Channel Water Coverage Ratio (Time[:,4])
                plt.plot(S[:, 0], Time[:, 2], 'k', label='Channel Saturation')  # 'k' for black
                plt.plot(S[:, 0], Time[:, 3], 'b', label='Channel Water Coverage Ratio')  # 'b' for blue

                # Label y-axis (left)
                plt.ylabel('Saturation/Water coverage ratio')
                
                # Enable secondary y-axis for pressure
                ax2 = plt.gca().twinx()
                ax2.plot(S[:, 0], S[:, 1], 'r', label='Channel Air Pressure (Pa)')  # 'r' for red
                ax2.set_yscale('log')  # Log scale for pressure
                ax2.set_ylabel('Channel Pressure (Pa)')

                # Label x-axis
                plt.xlabel('Time (s)')

                # Set font size for the plot
                plt.gca().tick_params(axis='both', labelsize=16)

                # Add legend
                plt.legend(loc='upper left')
                ax2.legend(loc='upper right')

                # Show the plot
                plt.show()

                # TODO: Check if plot is correct
                # TODO: There are several plots from here on which are commented out in the original code. Insert them

            if serpentine_plot == 1:
                pass
                #Serpentine_Plot() # TODO This is a bug, since no such function
            
            vtkploti += 1

            # ========== Vtk Plotting ========== 
            # This requires vtkwrite.m
            # This is skipped, can be added later TODO


            #  ========== Determine Droplet Regimes ========== 
                    
            # Calculate the regime volume
            RegimeVolume = np.array([np.sum(D[D[:, 5] == i, 0]) for i in range(1, 8)])

            # Calculate regime saturation
            RegimeSaturation = RegimeVolume / np.sum(RegimeVolume)

            # Save the regime saturation for the current iteration
            #RegimeSatSave[ic, :] = RegimeSaturation
            RegimeSatSave = RegimeSaturation

            # ==========  Determine Regimes Along the Channel (in progress) ========== 
            ct += dt

            if ct > write_precision:
                nseg = nz / 2           # number of channel segments
                nseg = 20
                Dseg = L / nseg

                segment = np.arange(0, L, Dseg).transpose()

                if time_value == 0:
                    avregime_new = np.zeros((segment.shape[0], 8))
                    plottimer = 0
                    av_vol_frac_i = np.zeros((1, segment.shape[0]))

                # save previous timestep average
                avregime_i = avregime_new


                # Iterate through each segment to determine the fraction of each regime
                for is_idx in range(len(segment)):
                    segment_start = segment[is_idx, 0] - Dseg
                    segment_end = segment[is_idx, 0] + Dseg
                    
                    # Calculate for each regime (1 to 7)
                    for regime in range(1, 8):
                        segment[is_idx, regime + 1] = np.sum(
                            (D[:, 4] < segment_end) & 
                            (D[:, 4] > segment_start) & 
                            (D[:, 5] == regime) & 
                            (D[:, 0] > 0)
                        )

                segmentRegime = segment
                segmentRegime[:, 1:] = segment[:, 1:] / np.sum(segment[:, 1:], 1)

                time_value += 1

                # ========== Tracking average without saving all timesteps ========== 
                segmentRegime = np.nan_to_num(segmentRegime)
                avregime_new = (((avregime_i) * (time_value -1)) + segmentRegime) / time_value
                # use new average procedure

                # TODO Plot visualisation in the original file commented out

                # ========== calculate volume fraction over length ========== 

                Vp = D[:, 0]            # volume of droplets
                zdrop = D[:, 4]         # z position of droplets
                Ap = D[:, 14]           # projected area of droplet
                Lp = D[:, 12]
                Regime_vol = D[:, 5]

                # zp +- LD/2 will find the maximum range of droplet in z direction
                zdrop_min = D[:, 4] - D[:, 12] / 2      # minimum z position of droplets
                zdrop_max = D[:, 4] - D[:, 12] / 2      # maximum z posistion of droplets

                # assume constant variation of volume over its length in the channel
                NZ = nseg + 1
                DZ = L / NZ

                # Define new grid for cell centers of the channel
                Z = np.arange(DZ/2, L, DZ) # TODO check if ends are the same

                # Distance from particle p and cell k
                Spk = np.sqrt((Z - zdrop)**2)

                # Overlap length
                Lpk = (Lp/2) + (DZ/2)
                overlap = Lpk - Spk

                # Set non-positive values to zero
                overlap[overlap < 0] = 0

                # Set values greater than cell length to cell length
                overlap[overlap > DZ] = DZ

                # Divide overlap by the length of the droplet in the length direction
                overlap = overlap / Lp

                # Replace NaN or Inf values with 0
                overlap[np.isnan(overlap) | np.isinf(overlap)] = 0

                # Cap overlap at 1
                overlap[overlap > 1] = 1

                # Calculate water cell volume
                water_cell_volume = overlap * Vp

                # Calculate regime cell volume (modify if needed)
                regime_cell_volume = overlap * Regime_vol

                # Calculate volume fraction
                vol_fraction = np.sum(water_cell_volume) / (DZ * W * H)
                
                av_vol_frac = ((av_vol_frac_i) * (time_value - 1)) * vol_fraction / (time_value)
                av_vol_frac_i = av_vol_frac

                # TODO Matlab Plot code

                ct = 0
                print("Loop inner done")
            # end of the if condition with ct > write_precision
            pt = t / endtime * 100
            print("Loop done")
            
            # end of the time loop for simulation with for i in tt...
        # moving on to saving results for each run it if is turned on
        # repeat() # NR: What does this do?

        if saveresults == 1:
            # Automatic saving of average properties about the run in a .csv file
            # Label file name to correspond to number of run i
            filename = Savefilename + str(sid) + ".csv"

            # create heading titles to identify what is saved
            SaveDataTitles = ["u_air", "u_water", "GDL_angle", "Wall_angle", "npores", "Max_Pore_radius", "L", "H", "W", "Total CPU Time (s)", "endtime", \
                    "dt", "Saturation", "ACR", "Pressure Drop", "Ndrops", \
                    "mean vol frac emerging", "mean vol frac isolated", "mean vol frac sidewall",
                    "mean vol frac corner", \
                            "mean vol frac capillary", "mean vol frac film"]

            # Create list of saved variables to identify operating and average
            # conditions

            
            timetest = 100.0 # just a fake value for timing
            # TODO Check if the output is correct, done by ChatGPT
            #  SaveData array containing various parameters
            dt = dt[0] # unpack array
            SaveData = np.array([uin, u_wav, theta_GDL * (180 / np.pi), theta_wall * (180 / np.pi),
                                npores, Rpmax, L, H, W, np.sum(timetest), S[-1, 0], dt, 
                                np.mean(Time[:, 2]), np.mean(Time[:, 3]), np.mean(S[:, 1]), 
                                np.mean(Time[:, 1]), np.mean(RegimeSatSave), current_density, Sto])

            # Write the data and headings to a .csv file
            filename = 'output.csv'
            with open(filename, mode='w', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(SaveData)

            # Append transient data for the simulation at every timestep
            # Quick fix debug
            RegimeSatSave = np.ones(S.shape[0])
            transient_data = np.column_stack([S[:, 0], Time[:, 2], Time[:, 3], S[:, 1], Time[:, 1], RegimeSatSave, S[:, 3]])
            with open(filename, mode='a', newline='') as file:
                writer = csv.writer(file)
                writer.writerows(transient_data)

            # Save new data for regime distribution along the channel
            if "segmentRegime" in globals(): # BUG not definied variables
                filename2 = f'{Savefilename}{sid}_length.csv'
                LengthData = np.column_stack([segmentRegime, avregime_new, vol_fraction, av_vol_frac])
                with open(filename2, mode='w', newline='') as file:
                    writer = csv.writer(file)
                    writer.writerows(LengthData)

            # TODO VTK Plotting in matlab code
        # end of the if condition
# end of the main loop for repeat 
