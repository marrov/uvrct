import uvrct as uv
import time

start = time.time()

case_variables = [17, 0.33, 0.28, 0.0425, 0.0115, 15, 15]

I_LSI = uv.compute_vafr(case_variables, model='LSI')
print(f'VAFR as computed by LSI is: {I_LSI:.3f}')
tol = 1E-2 # tolerance
err = 100 # Initial error
nls = 23

P, eff, L, r_outer, r_inner, ncr, ncz = case_variables
V_r = uv.reactor_volume(L, r_outer, r_inner)
# Make the outer and inner mesh objects
mesh_outer = uv.make_cyl_mesh(L, r_outer, ncr, ncz)
mesh_inner = uv.make_cyl_mesh(L, r_inner, ncr, ncz)

while err > tol*I_LSI:
    nls += 1
    I_i_outer = uv.per_cell_fluence_rate_MPSS(mesh_outer, P, eff, L, nls, no_abs=True)
    I_i_inner = uv.per_cell_fluence_rate_MPSS(mesh_inner, P, eff, L, nls, no_abs=True)
    I_MPSS = uv.volume_avg_fluence_rate(I_i_outer, mesh_outer, I_i_inner, mesh_inner, V_r)
    err = uv.np.abs(I_MPSS-I_LSI)
    print(f'The current error with {nls} points is: {100*err/I_LSI:.3f}%')
print(f'The VAFR without absorbtion as computed by MPSS with n_opt ({nls} points) is: {I_MPSS:.3f}')
print(f'Optimum value of line sources with {100*tol} % tolerance is: {nls}')


end = time.time()
print(f'Elapsed time is to find n_opt is {end - start:.3f}s')

start = time.time()
I_MPSS_1000 = uv.compute_vafr(case_variables, model='MPSS', nls=1000, optimize=False)
end = time.time()
t_1000 = end - start

start = time.time()
I_MPSS_opt = uv.compute_vafr(case_variables, model='MPSS', nls=nls, optimize=False)
end = time.time()
t_opt = end - start

print(f'Elapsed times are: t_1000={t_1000:.3f}s, t_opt={t_opt:.3f}s')
print(f'Speed up factor is: x{t_1000/t_opt:.3f}')
