"""
    Integrates the axon growth dynamic.

    // TODO: save the data as a .npz instead of plotting

"""

import os
from scipy.ndimage import convolve
import numpy as np
import matplotlib.pyplot as plt


kwargs = {'origin':'lower', 'interpolation':'sinc', 'cmap':'seismic',}

heaviside = np.vectorize(heaviside)
h = lambda x: x ** 2 * (3. - 2. * x)
boundary_condition = 'wrap'
dL = 1.0
stencil = (1.0 / (12.0 * dL * dL)) * np.array(
    [[0, 0, -1, 0, 0],
     [0, 0, 16, 0, 0],
     [-1, 16, -60, 16, -1],
     [0, 0, 16, 0, 0],
     [0, 0, -1, 0, 0]])

# id = '001'
# os.system('mkdir '+id)
L = np.array([40, 100])
L0 = np.array([0, 0])

neuron_position = np.array([20, 20])
neuron_radius = 15
gc_radius = 3
gc_position = np.array([neuron_position[0], neuron_position[1] + neuron_radius ])

small_box_size = 15

ngf = np.zeros(L)
phi = np.zeros(L)
psi = np.zeros(L)
mtb = np.zeros(L)
mf = np.zeros(L)
ml = np.zeros(L)

mf = init_cell(mf, neuron_position, neuron_radius-5)
ngf = init_growth_factor(ngf, neuron_position, neuron_radius )
phi = init_cell(phi, gc_position, gc_radius)
psi = init_cell(psi, neuron_position, neuron_radius)

nstep = 0  # actual time step
tstep = 25000  # total number of steps
dt = 0.001  # time interval
print_period = 5000
nprint = print_period # counter to output data


### mRNA
kd = 10
lambda_f = 0.005
lambda_l = 0.001
beta_m = 0.001
xi = 1
Mm = 1
v_m = 0.001
gamma = 1
######
dL = 1.0
lambda_ = 1.
chi = 5.
alpha_p = 1.0
D_ngf = 1.
D_mtb = 0.1

phi = smoothing(phi, 2000, stencil)  # growth cone
psi = smoothing(psi, 2000, stencil)  # neuron

for i in range(0, 100):
  mf = mf + dt * convolve(100 * mf + container(mf, psi, 0.001), stencil, mode=boundary_condition )

V_target = h(phi).sum()


r_cm = com(phi, V_target, ndim=2)
print(r_cm)
icalc = 0
ncalc = 250

while (nstep <= tstep):

    v = h(phi).sum()
    Lagrange = np.ones(L) * (V_target - v)

    il = np.int(r_cm[0]) - small_box_size
    ih = np.int(r_cm[0]) + small_box_size
    jl = np.int(r_cm[1]) - small_box_size
    jh = np.int(r_cm[1]) + small_box_size
    grad_phi = np.gradient(phi[il:ih, jl:jh], axis=1)
    grad_ngf = np.gradient(ngf)
    #grad_phi_abs = np.power(grad_phi[0], 2) + np.power(grad_phi[1], 2)
    chem = grad_phi #grad_phi[0]*grad_ngf[0] + grad_phi[1]*grad_ngf[1]

    coeff_ml = (Mm-ml)/Mm
    coeff_beta = beta_m*(1-xi)*ml/xi
    coeff_gamma = gamma*mf*xi*coeff_ml
    grad_coeff_ml = np.gradient(coeff_ml)
    grad_mtb = np.gradient(mtb)
    advection_ml = grad_coeff_ml[0]*grad_mtb[0] + grad_coeff_ml[1]*grad_mtb[1]

    chi_sum = (phi[il:ih, jl:jh]*ml[il:ih, jl:jh]).sum()
    chi_ = chi * chi_sum/(1+chi_sum)  #chi*mf/(1+mf)

    phi[il:ih, jl:jh], psi,  mtb, ngf, mf, ml = phi[il:ih, jl:jh] + (dt * (-chi_ * chem + convolve(phi[il:ih, jl:jh], stencil, mode='constant') +
                                              phi[il:ih, jl:jh] * (1.0 - phi[il:ih, jl:jh]) * ( (phi[il:ih, jl:jh] - 0.5) + Lagrange[il:ih, jl:jh]))), \
                                psi + (dt * (-convolve(8 * (psi * (1.0 - psi) * (psi - 0.5))
                                      + convolve(psi, stencil, mode=boundary_condition),
                                      stencil, mode=boundary_condition) +
                                      alpha_p * phi * ngf * psi * (1.0 - psi) )), \
                                mtb + (dt * (D_mtb * convolve(mtb + container(mtb, psi, 0.001), stencil, mode=boundary_condition))),\
                                ngf + (dt * (D_ngf * convolve(ngf, stencil, mode=boundary_condition) - psi * ngf )),\
                                mf + dt * ( convolve(kd * mf + 2*container(mf, psi, 0.001), stencil, mode=boundary_condition)
                                          - lambda_f * mf
                                          - gamma*mf*xi*coeff_ml
                                          + coeff_beta),\
                                ml + dt * ( convolve(ml + 2*container(ml, psi, 0.001), stencil, mode=boundary_condition)
                                          - lambda_l * ml
                                          + coeff_gamma
                                           - coeff_beta - advection_ml )

    if icalc >= ncalc:
        icalc = 0
        r_cm = com(phi, v, 2)
        mtb[np.int(r_cm[0]), np.int(r_cm[1])] = 10.0

    if (nprint >= print_period):
        #print(nstep, chi_)
        plt.subplot(1, 4, 1)
        plt.imshow(psi, alpha=1.0, **kwargs)
        plt.imshow(phi, alpha=0.2, **kwargs)
        plt.subplot(1, 4, 2)
        #plt.imshow(ngf, alpha=0.2, **{'interpolation':'sinc', 'cmap':'viridis'})
        plt.imshow(ml, alpha=1.0, origin='lower', interpolation='sinc')
        plt.subplot(1, 4, 3)
        plt.imshow(mf, alpha=1.0, origin='lower', interpolation='sinc')
        plt.subplot(1, 4, 4)
        plt.imshow(mtb, alpha=1.0, origin='lower', interpolation='sinc')
        plt.show()
        #plt.savefig(str(nstep)+'.png', dpi=300)

        nprint = 0
    icalc += 1
    nstep += 1
    nprint += 1

