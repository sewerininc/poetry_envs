# compute.py

import numpy as np
from scipy.integrate import quad 


def phi_0(k: np.array) -> float:
    k_length = np.sqrt(k[0]**2 + k[1]**2 + k[2]**2)
    # hydrogen_100 = 1 / (np.sqrt(np.pi) * c.a0.value**(3/2)) * np.exp(-r_size / c.a0.value)
    # tot_phi_0 = 1 / (2 * np.pi)**(3/2)
    solution = 2 ** (3/2) / np.pi / (1 + k_length**2)**2
    return solution


def M_SFA_VG(t_f: float, t_0: float, k: np.array, w_sim: float,
             A_0: float, E_0: float, phase: float, epsilon: float, T_f: float, 
             #sin4cos2, sin4sin2, sin2cos, sin2sin
            ) -> float:
    """
    sin4cos2 and sin4sin2 needs to be the integral and type Piecewise 
    """
    
    i = complex(0, 1)
    
    k_squared = k[0]**2 + k[1]**2 + k[2]**2
    
    # I need to make for kx and ky
    
    test_sin4cos2 = lambda time: np.sin(np.pi * time / T_s)**4 * np.cos(w_s * time + phase_s)**2
    test_sin4sin2 = lambda time: np.sin(np.pi * time / T_s)**4 * np.sin(w_s * time + phase_s)**2
    test_sin2cos = lambda time: np.sin(np.pi * time / T_s)**2 * np.cos(w_s * time + phase_s)
    test_sin2sin = lambda time: np.sin(np.pi * time / T_s)**2 * np.sin(w_s * time + phase_s)
    
    # kx here first k[0]
    if k[0] != 0:
        i_x = lambda ts: k[0] * np.sin(np.pi * ts / T_f)**2 * np.cos(w_sim * ts + phase) * np.exp(
            (-1) * i * E_0 * ts +
                 # k^2
            i * (k_squared * ts + 
                 # k*A
                 # A_0 / np.sqrt(1 + epsilon**2) * (k[0] * sin2cos.subs({t: ts}).evalf() +  
                 #                                 k[1] * epsilon * sin2sin.subs({t: ts}).evalf()
                 #                                ) + 
                 A_0 / np.sqrt(1 + epsilon**2) * (k[0] * quad(test_sin2cos, 0, ts, epsabs=1e-3, 
                                                              full_output=False)[0] +
                                                  k[1] * epsilon * quad(test_sin2sin, 0, ts, epsabs=1e-3,
                                                                        full_output=False)[0]
                                                 ) + 
                 # A ** 2
                 # A_0**2 / (1 + epsilon**2) * (sin4cos2.subs({t: ts}).evalf() +
                 #                             epsilon**2 * sin4sin2.subs({t: ts}).evalf()
                 #                            )
                 A_0**2 / (1 + epsilon**2) * (quad(test_sin4cos2, 0, ts, epsabs=1e-3, full_output=False)[0] +
                                              epsilon**2 * quad(test_sin4sin2, 0, ts, epsabs=1e-3,
                                                                full_output=False)[0]
                                             )
                )
        )
        sol_large_x = quad(i_x, t_0, t_f, epsabs=1e-3, full_output=False, complex_func=True)
    else:
        sol_large_x = [0, 0]
    
    # ky here first k[1]
    if k[1] != 0:
        i_y = lambda ts: k[1] * np.sin(np.pi * ts / T_f)**2 * epsilon * np.sin(w_sim * ts + phase) * np.exp(
            (-1) * i * E_0 * ts +
                 # k^2
            i * (k_squared * ts + 
                 # k*A
                 # A_0 / np.sqrt(1 + epsilon**2) * (k[0] * sin2cos.subs({t: ts}).evalf() +  
                 #                                 k[1] * epsilon * sin2sin.subs({t: ts}).evalf()
                 #                                ) + 
                 A_0 / np.sqrt(1 + epsilon**2) * (k[0] * quad(test_sin2cos, 0, ts, epsabs=1e-3, 
                                                              full_output=False)[0] +
                                                  k[1] * epsilon * quad(test_sin2sin, 0, ts, epsabs=1e-3,
                                                                        full_output=False)[0]
                                                 ) + 
                 # A ** 2
                 # A_0**2 / (1 + epsilon**2) * (sin4cos2.subs({t: ts}).evalf() +
                 #                             epsilon**2 * sin4sin2.subs({t: ts}).evalf()
                 #                            )
                 A_0**2 / (1 + epsilon**2) * (quad(test_sin4cos2, 0, ts, epsabs=1e-3, full_output=False)[0] +
                                              epsilon**2 * quad(test_sin4sin2, 0, ts, epsabs=1e-3,
                                                                full_output=False)[0]
                                             )
                )
        )
        sol_large_y = quad(i_y, t_0, t_f, epsabs=1e-3, full_output=False, complex_func=True)
    else:
        sol_large_y = [0, 0]
    
    
    full_integral =  A_0 / np.sqrt(1 + epsilon**2) * (sol_large_y[0] + sol_large_x[0])
    
    M = (-1) * i * phi_0(k) * full_integral
    
    return M
    
    
def dp_dk(M: complex) -> float:
    return M.real**2 + M.imag**2



def compute_values(kx, kys, T_s, w_s, A_0_s, E_0_s, phase_s):
    return [dp_dk(M_SFA_VG(t_f=T_s, t_0=0, k=np.array([kx, ky, 0]), w_sim=w_s, A_0=A_0_s, E_0=E_0_s, 
                       phase=phase_s,
                       epsilon=0.5, T_f=T_s)) for ky in kys]
