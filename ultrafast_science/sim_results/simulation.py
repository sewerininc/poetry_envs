import numpy as np
from multiprocessing import Pool
from compute import compute_values

I = 5e13
convertion_I_au = 3.51e16
I_au = I / convertion_I_au
F_0 = np.sqrt(I_au)

conversion_ev_au = 27.21
w_800 = 1.56 / conversion_ev_au
w_400 = w_800 * 800 / 400
w_s = w_400

A_0_s = F_0 / w_400

E_0_s = -13.6 / conversion_ev_au

N_s = 7
T_C_s = 2 * np.pi / w_s
T_s = N_s * T_C_s

phase_s = 0.5
epsilon_s = 0.5

length = 5
start_val = -4
end_val = 4
kxs = np.linspace(start_val, end_val, length)
kys = np.linspace(start_val, end_val, length)

pool = Pool()  # Create a pool of worker processes
results = pool.starmap(compute_values, [(kx, kys, T_s, w_s, A_0_s, E_0_s, phase_s) for kx in kxs])

kx_plot = np.repeat(kxs, length)
ky_plot = np.tile(kys, length)
values_plot = np.concatenate(results)

pool.close()  # Close the pool
pool.join()   # Join all the processes

print(kx_plot)
print(ky_plot)
print(values_plot)
