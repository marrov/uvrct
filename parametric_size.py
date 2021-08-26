import uvrct as uv

uv.plt.style.use('seaborn-pastel')
uv.plt.rcParams['axes.labelsize'] = 18
uv.plt.rcParams['legend.fontsize'] = 14
uv.plt.rcParams['xtick.labelsize'] = 14
uv.plt.rcParams['ytick.labelsize'] = 14
uv.plt.rcParams['figure.autolayout'] = True


def run_cstr_uv(mech, T, P, C, case_variables, sim_variables, model):

    gas = uv.set_gas(mech, T, P, C)
    P, eff, L, r_outer, r_inner, ncr, ncz = case_variables
    V_r = uv.reactor_volume(L, r_outer, r_inner)
    p_valve_coeff, max_p_rise, residence_t, max_simulation_t = sim_variables
    time_history = uv.run_simulation(gas, V_r, residence_t, p_valve_coeff, case_variables,
                                     max_simulation_t, model=model)
    # Calculate SO2 removal rate
    removal_so2 = (
        1-(time_history['SO2'].iloc[-1]/time_history['SO2'].iloc[0]))*100

    return removal_so2


if __name__ == "__main__":

    # Defining simulation variables
    initial_SO2_concentration = 300/1e6
    molarity_ratios = uv.np.linspace(0, 1, 21)
    mech = 'data/photolysis.cti'
    T = 298
    P = uv.ct.one_atm
    case_variables = [17, 0.33, 0.28, 0.0425, 0.0115, 15, 15]
    sim_variables = [0.01, 0.01, 100, 1000]
    model = 'LSI'
    sizes = [0.0425, 0.0675, 0.0820]
    uv.plt.figure()

    for size in sizes:
        removal_so2 = []
        for molarity_ratio in molarity_ratios:
            C = {'N2': 0.78054 - (initial_SO2_concentration * (1 + molarity_ratio)),
                 'O2' : 0.20946,
                 'O3' : initial_SO2_concentration*molarity_ratio,
                 'SO2': initial_SO2_concentration,
                 'H2O': 0.01}
            case_variables[3] = size
            removal_so2.append(run_cstr_uv(mech, T, P, C, case_variables, sim_variables, model))
        uv.plt.plot(molarity_ratios, removal_so2, '-', label=r'$r_{out}$='+str(size), linewidth=2.5)
    
    uv.plt.xlabel('Molarity ratio (O$_3$/SO$_2$)')
    uv.plt.ylabel('SO$_2$ removal [%]')
    uv.plt.legend(loc='lower right')
    uv.plt.show()
