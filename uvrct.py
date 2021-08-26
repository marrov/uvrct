"""
Code to compute a CSTR UV reactor with Cantera
"""
import re
import time
import pandas as pd
import numpy as np
import cantera as ct
import discretize
import matplotlib.pyplot as plt
from scipy.constants import N_A, h, c


def reactor_volume(L, r_outer, r_inner):
    """
    Computes the volume of an annular reactor (hollow cylinder).

    Parameters [units]
    ------------------
    L : int, float
        Length of the reactor/lamp [m]
    r_outer : int, float
        Outer radius [m]
    r_outer : int, float
        Inner radius [m]

    Returns [units]
    ---------------
    V_r : float
        The reactor volume [m**2]
    """
    return np.pi*L*(r_outer**2-r_inner**2)


def make_cyl_mesh(L, r, ncr, ncz):
    """
    Makes a cylindrical mesh object with one cell in the azimuthal direction.

    Parameters [units]
    ------------------
    L : int, float
        Length of the reactor/lamp [m]
    r : int, float
        Cylinder radius [m]
    ncr : int
        Number of cells in the radial direction
    ncz : int
        Number of cells in the axial direction

    Returns [units]
    ---------------
    mesh : obj
        The mesh cylinder object
    """
    dr = r/ncr              # cell width r
    dz = L/ncz              # cell width z
    hr = dr*np.ones(ncr)
    hz = dz*np.ones(ncz)
    offset = [0, 0, -L/2]   # offset to center z on half the lamp length
    return discretize.CylMesh([hr, 2*np.pi, hz], offset)


def sigma_to_eps(sigma):
    """
    Converts the molar absorption cross-section of a gas to the molar 
    absorption coefficient in S.I. units.

    Parameters [units]
    ------------------
    sigma : int, float
        Absorbtion cross section (log base e) [cm**2]

    Returns [units]
    ---------------
    eps : float
        molar absorption coefficient (log base 10) [m**2 mol**-1]
    """
    return sigma*N_A*1e-4/np.log(10)


def ppm_to_C(ppm):
    """
    Converts ppm (parts per million) concentration of a gas to
    concentration in S.I. units.

    Parameters [units]
    ------------------
    ppm : int, float
        ppm concentration of a gas [-]

    Returns [units]
    ---------------
    C : float
        concentration [mol m**-3]
    """
    rho_air = 1.2  # [kg m**-3]
    Mr_air = 28.96  # [g mol**-1] Note: this is for DRY air
    Vm_air = molar_volume(rho_air, Mr_air)  # [m**3 mol**-1]
    C = ppm/(Vm_air*1e6)
    return C


def molar_volume(rho, Mr):
    """
    Computes the molar volume of a gas.

    Parameters [units]
    ------------------
    rho : int, float
        density [kg m**-3]
    Mr : int, float
        molar mass [g mol**-1]

    Returns [units]
    ---------------
    Vm : float
        molar volume [m**3 mol**-1]
    """
    return (Mr/rho)*1e-3


def exp_medium_absorption(r, ppm=300):
    """
    Computes the absorption of the medium considered in the fluence 
    rate models. Initially, I will treat only SO2 as absorbing species.

    Parameters [units]
    ------------------
    r : int, float, array
        Radial position (i.e. of a mesh element) [m]
    ppm : int, float
        ppm concentration of a gas [-]

    Returns [units]
    ---------------
    ema : int, float, array
        Exponential medium absorption [-]
    """
    sigma = 1.551258e-19  # Value for SO2 [cm**2]
    sigma_air = sigma
    eps_air = sigma_to_eps(sigma_air)
    C_air = ppm_to_C(ppm)
    return np.exp(-np.log(10)*eps_air*C_air*r)


def coord_ls(L, nls):
    """
    Computes the coordinates of the number light source elements 
    which have been chosen to model the UV lamp.

    Parameters [units]
    ------------------
    L : int, float
        Length of the reactor/lamp [m]
    nls : int
        Number of light source elements

    Returns [units]
    ---------------
    coord_ls : array
        Coordinates of the number light source elements [m]
    """
    d_ls = L/nls
    return np.arange(-L/2, L/2+d_ls, d_ls)


def optimize_nls(case_variables, tol=1E-2):
    """
    Computes the optimum number of light sources (nls) for which to
    model a UV lamp (for the MPSS or MSSS models). This optimum is 
    the lowest nls that makes the fluence rate computed by the model
    (without medium absorption) match the fluence rate computed by 
    the LSI model for a given tolerance. Based on Powell and Lawryshyn 
    2015.
    Note: this routine takes some time, proceed with caution.

    Parameters [units]
    ------------------
    case_variables : array
        Array contining the following case variables:
        case_variables = [P, L, eff, r_outer, r_inner, ncr, ncz]    
        P       = lamp power [W]
        L       = lamp length [m]
        eff     = lamp efficiency at 254 nm wavelength
        r_outer = reactor radius [m]
        r_inner = quartz sleeve radius [m]
        ncr     = number of mesh cells in r
        ncz     = number of mesh cells in z
    tol : int, float
        Tolerance for the optimization routine

    Returns [units]
    ---------------
    nls : int
        (Optimum) number of light source elements
    """
    P, eff, L, r_outer, r_inner, ncr, ncz = case_variables
    V_r = reactor_volume(L, r_outer, r_inner)
    mesh_outer = make_cyl_mesh(L, r_outer, ncr, ncz)
    mesh_inner = make_cyl_mesh(L, r_inner, ncr, ncz)
    err, nls = 100, 23  # Initialize error and smallest possible number of light sources
    # Compute reference LSI vafr
    I_LSI = compute_vafr(case_variables, model='LSI')
    # Loop to find the optimum value for the number of light sources
    while err > tol*I_LSI:
        nls += 1
        I_i_outer = per_cell_fluence_rate_MPSS(
            mesh_outer, P, eff, L, nls=nls, no_abs=True)
        I_i_inner = per_cell_fluence_rate_MPSS(
            mesh_inner, P, eff, L, nls=nls, no_abs=True)
        I_MPSS = volume_avg_fluence_rate(
            I_i_outer, mesh_outer, I_i_inner, mesh_inner, V_r)
        err = np.abs(I_MPSS-I_LSI)
    return nls


def per_cell_fluence_rate_RAD(mesh, P, eff, L, no_abs=False, ppm=300):
    """
    Computes the fluence rate for each cell in a given mesh following
    the RAD model (Coenen 2013).

    Parameters [units]
    ------------------
    mesh : obj
        Mesh object
    P : int, float
        Lamp power [W]
    eff : float
        Lamp efficiency at 254 nm wavelength
    L : int, float
        Length of the reactor/lamp [m]
    no_abs : bool
        Flag to enable/disable medium absorption
    ppm : int, float
        ppm concentration of a gas [-]

    Returns [units]
    ---------------
    I_i : array
        Fluence rate for each cell in the mesh [W m**-2]
    """
    R = mesh.gridCC[:, 0]
    if no_abs:
        I_i = ((P*eff)/(2*np.pi*R*L))
    else:
        I_i = ((P*eff)/(2*np.pi*R*L))*exp_medium_absorption(R, ppm=ppm)
    return I_i


def per_cell_fluence_rate_LSI(mesh, P, eff, L):
    """
    Computes the fluence rate for each cell in a given mesh following
    the LSI model (Blatchley 1998). It assumes no medium absorption

    Parameters [units]
    ------------------
    mesh : obj
        Mesh object
    P : int, float
        Lamp power [W]
    eff : float
        Lamp efficiency at 254 nm wavelength
    L : int, float
        Length of the reactor/lamp [m]

    Returns [units]
    ---------------
    I_i : array
        Fluence rate for each cell in the mesh [W m**-2]
    """
    R = mesh.gridCC[:, 0]
    H = mesh.gridCC[:, 2]
    I_i = ((P*eff)/(4*np.pi*L*R))*(np.arctan((L/2+H)/R)+np.arctan((L/2-H)/R))
    return I_i


def per_cell_fluence_rate_RAD_LSI(mesh, P, eff, L):
    """
    Computes the fluence rate for each cell in a given mesh following
    the RAD-LSI model (Liu 2004). It assumes no medium absorption

    Parameters [units]
    ------------------
    mesh : obj
        Mesh object
    P : int, float
        Lamp power [W]
    eff : float
        Lamp efficiency at 254 nm wavelength
    L : int, float
        Length of the reactor/lamp [m]

    Returns [units]
    ---------------
    I_i : array
        Fluence rate for each cell in the mesh [W m**-2]
    """
    R = mesh.gridCC[:, 0]
    H = mesh.gridCC[:, 2]
    I_i_LSI = ((P*eff)/(4*np.pi*L*R)) * \
        (np.arctan((L/2+H)/R)+np.arctan((L/2-H)/R))
    I_i_RAD = ((P*eff)/(2*np.pi*R*L))
    I_i = np.minimum(I_i_RAD, I_i_LSI)
    return I_i


def per_cell_fluence_rate_MPSS(mesh, P, eff, L, nls, no_abs=False, ppm=300):
    """
    Computes the fluence rate for each cell in a given mesh following
    the MPSS model (Bolton 2000).

    Parameters [units]
    ------------------
    mesh : obj
        Mesh object
    P : int, float
        Lamp power [W]
    eff : float
        Lamp efficiency at 254 nm wavelength
    L : int, float
        Length of the reactor/lamp [m]
    nls : int
        Number of light source elements
    no_abs : bool
        Flag to enable/disable medium absorption
    ppm : int, float
        ppm concentration of a gas [-]

    Returns [units]
    ---------------
    I_i : array
        Fluence rate for each cell in the mesh [W m**-2]
    """
    x_ls = coord_ls(L, nls)
    R = mesh.gridCC[:, 0]
    H = mesh.gridCC[:, 2]
    I_i = np.zeros(len(R))
    for i in range(len(R)):
        r = R[i]*np.ones(len(x_ls))
        h = H[i]*np.ones(len(x_ls))
        rho = np.sqrt(r**2+(h-x_ls)**2)
        if no_abs:
            I_ls = (((P*eff)/nls)/(4*np.pi*rho**2))
        else:
            I_ls = (((P*eff)/nls)/(4*np.pi*rho**2)) * \
                exp_medium_absorption(r, ppm)
        I_i[i] = sum(I_ls)
    return I_i


def per_cell_fluence_rate_MSSS(mesh, P, eff, L, nls, no_abs=False, ppm=300):
    """
    Computes the fluence rate for each cell in a given mesh following
    the MSSS model (Liu 2004).

    Parameters [units]
    ------------------
    mesh : obj
        Mesh object
    P : int, float
        Lamp power [W]
    eff : float
        Lamp efficiency at 254 nm wavelength
    L : int, float
        Length of the reactor/lamp [m]
    nls : int
        Number of light source elements
    no_abs : bool
        Flag to enable/disable medium absorption
    ppm : int, float
        ppm concentration of a gas [-]

    Returns [units]
    ---------------
    I_i : array
        Fluence rate for each cell in the mesh [W m**-2]
    """
    x_ls = coord_ls(L, nls)
    R = mesh.gridCC[:, 0]
    H = mesh.gridCC[:, 2]
    I_i = np.zeros(len(R))
    for i in range(len(R)):
        r = R[i]*np.ones(len(x_ls))
        h = H[i]*np.ones(len(x_ls))
        rho = np.sqrt(r**2+(h-x_ls)**2)
        if no_abs:
            I_ls = (((P*eff)/nls)/(4*np.pi*rho**2))
        else:
            I_ls = (((P*eff)/nls)/(4*np.pi*rho**2)) * \
                exp_medium_absorption(r, ppm)
        theta = np.arctan(np.abs(h-x_ls)/r)
        I_ls = I_ls*np.cos(theta)
        I_i[i] = sum(I_ls)
    return I_i


def per_cell_fluence_rate(mesh, P, eff, L, model='RAD', nls=40, ppm=300):
    """
    Generic function that computes the fluence rate for each cell 
    in a given mesh for various model.

    Parameters [units]
    ------------------
    mesh : obj
        Mesh object
    P : int, float
        Lamp power [W]
    eff : float
        Lamp efficiency at 254 nm wavelength
    L : int, float
        Length of the reactor/lamp [m]
    model : str
        name of the model ('RAD', 'LSI', 'MPSS', 'MSSS' or 'RAD_LSI')
    nls : int
        Number of light source elements
    ppm : int, float
        ppm concentration of a gas [-]

    Returns [units]
    ---------------
    I_i : array
        Fluence rate for each cell in the mesh [W m**-2]
    """
    if isinstance(model, str) == 0:
        raise TypeError("The variable 'model' has to be a string")

    if model == 'RAD':
        I_i = per_cell_fluence_rate_RAD(mesh, P, eff, L, ppm=ppm)
    elif model == 'LSI':
        I_i = per_cell_fluence_rate_LSI(mesh, P, eff, L)
    elif model == 'MPSS':
        I_i = per_cell_fluence_rate_MPSS(mesh, P, eff, L, nls=nls, ppm=ppm)
    elif model == 'MSSS':
        I_i = per_cell_fluence_rate_MSSS(mesh, P, eff, L, nls=nls, ppm=ppm)
    elif model == 'RAD_LSI':
        I_i = per_cell_fluence_rate_RAD_LSI(mesh, P, eff, L)
    else:
        raise ValueError(
            "Only 'RAD', 'LSI', 'MPSS', 'MSSS' or 'RAD_LSI' models are implemented")
    return I_i


def volume_avg_fluence_rate(I_i_outer, mesh_outer, I_i_inner, mesh_inner, V_r):
    """
    Computes the volume averaged fluence rate (vafr) for the
    annular reactor (hollow cylinder) by subtracting the
    average of the small cylinder (_inner) from the large 
    cylinder (_outer).

    Parameters [units]
    ------------------
    I_i_outer : array
        Fluence rate for each cell in the outer cylinder mesh [W m**-2]
    mesh_outer : obj
        outer cylinder mesh object
    I_i_inner : array
        Fluence rate for each cell in the inner cylinder mesh [W m**-2]
    mesh_inner : obj
        Inner cylinder mesh object
    V_r : float
        The reactor volume [m**2]

    Returns [units]
    ---------------
    I_vol_avg : array
        Volume average fluence rate for the annular reactor [W m**-2]
    """
    # volume averaged fluence rate (vafr)
    I_outer = np.sum(I_i_outer*mesh_outer.vol)
    I_inner = np.sum(I_i_inner*mesh_inner.vol)
    return (I_outer - I_inner)/V_r


def compute_vafr(case_variables, model='RAD', nls=40, optimize=False, ppm=300):
    """
    Generic function to compute the volume average fluence rate.

    Parameters [units]
    ------------------
    case_variables : array
        Array contining the following case variables:
        case_variables = [P, L, eff, r_outer, r_inner, ncr, ncz]    
        P       = lamp power [W]
        L       = lamp length [m]
        eff     = lamp efficiency at 254 nm wavelength
        r_outer = reactor radius [m]
        r_inner = quartz sleeve radius [m]
        ncr     = number of mesh cells in r
        ncz     = number of mesh cells in z
    model : str
        name of the model ('RAD', 'LSI', 'MPSS', 'MSSS' or 'RAD_LSI')
    nls : int
        Number of light source elements
    no_abs : bool
        Flag to enable/disable nls optimization for the MPSS or MSSS models
    ppm : int, float
        ppm concentration of a gas [-]

    Returns [units]
    ---------------
    I_vol_avg : array
        Volume average fluence rate for the annular reactor [W m**-2]
    """
    # Unpack the case variables
    P, eff, L, r_outer, r_inner, ncr, ncz = case_variables

    # Compute reactor volume
    V_r = reactor_volume(L, r_outer, r_inner)

    # Make the outer and inner mesh objects
    mesh_outer = make_cyl_mesh(L, r_outer, ncr, ncz)
    mesh_inner = make_cyl_mesh(L, r_inner, ncr, ncz)

    # Find optimum nls if required
    if (model == 'MPSS' or model == 'MSSS') and optimize == True:
        nls = optimize_nls(case_variables, tol=1E-2)

    # Calculate the per-cell fluence rate for the inner and outer meshes
    I_i_outer = per_cell_fluence_rate(
        mesh_outer, P, eff, L, model=model, nls=nls, ppm=ppm)
    I_i_inner = per_cell_fluence_rate(
        mesh_inner, P, eff, L, model=model, nls=nls, ppm=ppm)

    # Compute volume average LSI fluence rate
    I_vol_avg = volume_avg_fluence_rate(
        I_i_outer, mesh_outer, I_i_inner, mesh_inner, V_r)
    return I_vol_avg


def set_gas(mechanism, T, P, C):
    """
    Define the properties of the cantera gas object.

    Parameters [units]
    ------------------
    mechanism : string
        Path to the CTI mechanism to be used
    T : int, float
        Temperature [K]
    P : int, float
        Pressure [in atm]
    C : dict
        Dictionary containing the concentrations of the different gas species [molar fraction]

    Returns [units]
    ---------------
    gas : obj
        cantera gas object
    """
    gas = ct.Solution('data/photolysis.cti')
    gas.TPX = T, P, C
    return gas


def set_reactor(gas, V_r, residence_t, p_valve_coeff):
    """
    Set-up the reactor network using Cantera functions and syntax.

    Parameters [units]
    ------------------
    gas : obj
        Cantera gas object
    V_r : float
        The reactor volume [m**2]
    residence_t : int, float
        Time that the gas will stay in the reactor [s]
    p_valve_coeff : int, float
        This is the "conductance" of the pressure valve and will 
        determine its efficiency in holding the reactor pressure 
        to the desired conditions. Set to 0.01 for default

    Returns [units]
    ---------------
    stirred_reactor : obj
        Cantera reactor object
    reactor_network : obj
        Cantera network object
    """
    fuel_air_mixture_tank = ct.Reservoir(gas)
    exhaust = ct.Reservoir(gas)
    stirred_reactor = ct.IdealGasReactor(gas, energy='off', volume=V_r)
    mass_flow_controller = ct.MassFlowController(
        upstream=fuel_air_mixture_tank, downstream=stirred_reactor, mdot=stirred_reactor.mass/residence_t)
    pressure_regulator = ct.Valve(
        upstream=stirred_reactor, downstream=exhaust, K=p_valve_coeff)
    reactor_network = ct.ReactorNet([stirred_reactor])
    return stirred_reactor, reactor_network


def intialize_results(stirred_reactor):
    """
    Initializes the DataFrame where results will be stored.

    Parameters [units]
    ------------------
    stirred_reactor : obj
        Cantera reactor object

    Returns [units]
    ---------------
    time_history : DataFrame
        Pandas DataFrame to store time evolution of results. It
        includes the gas's mass, volume, temperature and species
        concentration
    """
    column_names = [stirred_reactor.component_name(
        item) for item in range(stirred_reactor.n_vars)]
    time_history = pd.DataFrame(columns=column_names)
    return time_history


def update_results(time_history, stirred_reactor, t):
    """
    Updates the results at a given time instant.

    Parameters [units]
    ------------------
    time_history : DataFrame
        Pandas DataFrame to store time evolution of results. It
        includes the gas's mass, volume, temperature and species
        concentration
    stirred_reactor : obj
        Cantera reactor object
    t : int, float
        Current simulation time [s]

    Returns [units]
    ---------------
    time_history : DataFrame
        Pandas DataFrame to store time evolution of results. It
        includes the gas's mass, volume, temperature and species
        concentration
    """
    state = np.hstack([stirred_reactor.mass, stirred_reactor.volume,
                       stirred_reactor.T, stirred_reactor.thermo.X])
    time_history.loc[t] = state
    return time_history


def uv_update_reactions(gas, stirred_reactor, case_variables, model='RAD', nls=40):
    """
    Modifies the rate constant of the photolysis reaction of ozone. This
    has to be done for each loop iteration as the concentrations of
    the different gas species are changing, which will affect the absorption
    of the medium.

    Parameters [units]
    ------------------
    gas : obj
        Cantera gas object
    stirred_reactor : obj
        Cantera reactor object
    case_variables : array
        Array contining the following case variables:
        case_variables = [P, L, eff, r_outer, r_inner, ncr, ncz]    
        P       = lamp power [W]
        L       = lamp length [m]
        eff     = lamp efficiency at 254 nm wavelength
        r_outer = reactor radius [m]
        r_inner = quartz sleeve radius [m]
        ncr     = number of mesh cells in r
        ncz     = number of mesh cells in z
    model : str
        name of the model ('RAD', 'LSI', 'MPSS', 'MSSS' or 'RAD_LSI')
    nls : int
        Number of light source elements

    Returns [units]
    ---------------
    None
    """
    ppm = stirred_reactor.thermo.X[stirred_reactor.component_index(
        'SO2')-3]*1e6
    k_uv = (0.9) * ((254E-9)/(N_A*h*c)) * (np.log(10)) * sigma_to_eps(1.132935E-17) * \
        compute_vafr(case_variables, model=model, nls=nls, ppm=ppm)
    ID = np.size(gas.reactions())-1
    reaction = gas.reactions()[ID]
    reaction.rate = ct.Arrhenius(A=k_uv, b=0, E=0)
    gas.modify_reaction(ID, reaction)


def time_loop(gas, reactor_network, stirred_reactor, case_variables, time_history, max_simulation_t, model='RAD'):
    """
    Computes the time loop iterations for the chemical reactor
    until the maximum simulation time is reached.

    Parameters [units]
    ------------------
    gas : obj
        Cantera gas object
    stirred_reactor : obj
        Cantera reactor object
    reactor_network : obj
        Cantera network object
    case_variables : array
        Array contining the following case variables:
        case_variables = [P, L, eff, r_outer, r_inner, ncr, ncz]    
        P       = lamp power [W]
        L       = lamp length [m]
        eff     = lamp efficiency at 254 nm wavelength
        r_outer = reactor radius [m]
        r_inner = quartz sleeve radius [m]
        ncr     = number of mesh cells in r
        ncz     = number of mesh cells in z
    time_history : DataFrame
        Pandas DataFrame to store time evolution of results. It
        includes the gas's mass, volume, temperature and species
        concentration
    max_simulation_t : int
        Maximum simulation time [s]
    model : str
        name of the model ('RAD', 'LSI', 'MPSS', 'MSSS' or 'RAD_LSI')

    Returns [units]
    ---------------
    time_history : DataFrame
        Pandas DataFrame to store time evolution of results. It
        includes the gas's mass, volume, temperature and species
        concentration
    """
    counter, t = 1, 0
    while t < max_simulation_t:
        if (t == 0) or (counter % 10 == 0):  # Update every 10 iterations
            if (model == 'MPSS') or (model == 'MSSS'):
                if t == 0:
                    nls = optimize_nls(case_variables, tol=1E-2)
                uv_update_reactions(
                    gas, stirred_reactor, case_variables, model=model, nls=nls)
            else:
                uv_update_reactions(gas, stirred_reactor,
                                    case_variables, model=model)
            time_history = update_results(time_history, stirred_reactor, t)
        t = reactor_network.step()
        counter += 1
    return time_history


def run_simulation(gas, V_r, residence_t, p_valve_coeff, case_variables, max_simulation_t, model='RAD'):
    """
    Generic function to run the chemical reactor with UV light
    modeling.

    Parameters [units]
    ------------------
    gas : obj
        Cantera gas object
    V_r : float
        The reactor volume [m**2]
    residence_t : int, float
        Time that the gas will stay in the reactor [s]
    p_valve_coeff : int, float
        This is the "conductance" of the pressure valve and will 
        determine its efficiency in holding the reactor pressure 
        to the desired conditions. Set to 0.01 for default
    case_variables : array
        Array contining the following case variables:
        case_variables = [P, L, eff, r_outer, r_inner, ncr, ncz]    
        P       = lamp power [W]
        L       = lamp length [m]
        eff     = lamp efficiency at 254 nm wavelength
        r_outer = reactor radius [m]
        r_inner = quartz sleeve radius [m]
        ncr     = number of mesh cells in r
        ncz     = number of mesh cells in z
    max_simulation_t : int
        Maximum simulation time [s]
    model : str
        name of the model ('RAD', 'LSI', 'MPSS', 'MSSS' or 'RAD_LSI')

    Returns [units]
    ---------------
    time_history : DataFrame
        Pandas DataFrame to store time evolution of results. It
        includes the gas's mass, volume, temperature and species
        concentration
    """
    # Define the reactor object
    stirred_reactor, reactor_network = set_reactor(
        gas, V_r, residence_t, p_valve_coeff)
    # Initialize the results dataFrame
    time_history = intialize_results(stirred_reactor)
    # Run the time loop
    time_history = time_loop(gas, reactor_network, stirred_reactor, case_variables, time_history,
                             max_simulation_t,  model=model)
    return time_history


def plot_results(time_history, specie, log=False):
    """
    Plots the results from the reactor simulation.

    Parameters [units]
    ------------------
    time_history : DataFrame
        Pandas DataFrame to store time evolution of results. It
        includes the gas's mass, volume, temperature and species
        concentration
    specie : str
        Specie for which to plot the time evolution
    log : bool
        Flag to enable/disable the x axis as logarithmic

    Returns [units]
    ---------------
    None
    """
    plt.style.use('seaborn-pastel')
    plt.rcParams.update({'axes.labelsize': 18, 'xtick.labelsize': 14,
                         'ytick.labelsize': 14, 'legend.fontsize': 14, 'figure.autolayout': True})
    plt.figure()
    if log == True:
        plt.semilogx(time_history.index,
                     time_history[specie]*1e6, '-', linewidth=2.5)
    else:
        plt.plot(time_history.index,
                 time_history[specie]*1e6, '-', linewidth=2.5)
    plt.xlabel('Time (s)')
    plt.ylabel(re.sub(r"([0-9]+(\.[0-9]+)?)", r"_\1",
                      '$'+specie+'$' + r' Mole Fraction : $ppmv$'))
    plt.show()


if __name__ == "__main__":
    start = time.time()

    print('Defining the gas properties')
    gas = set_gas('data/photolysis.cti', 298, ct.one_atm,
                  {'N2': 0.78009, 'O2': 0.20946, 'O3': 0.00015, 'SO2': 0.00030, 'H2O': 0.01})

    print('Defining the geometry variables')
    case_variables = P, eff, L, r_outer, r_inner, ncr, ncz = [
        17, 0.33, 0.28, 0.0425, 0.0115, 15, 15]
    V_r = reactor_volume(L, r_outer, r_inner)

    print('Defining the simulation variables')
    p_valve_coeff, max_p_rise, residence_t, max_simulation_t = 0.01, 0.01, 100, 1000

    print('Running the simulation')
    time_history = run_simulation(gas, V_r, residence_t, p_valve_coeff, case_variables,
                                  max_simulation_t, model='MPSS')

    end = time.time()

    print(f'Elapsed time: {end - start:.3f}s')

    print('Plotting the simulation results')
    plot_results(time_history, specie='SO2', log=True)
