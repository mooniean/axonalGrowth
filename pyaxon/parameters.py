import numpy as np
from functions import *

kwargs = {'origin':'lower', 'interpolation':'sinc', 'cmap':'seismic',}

alpha_ = 10e-4 # restrainement constant
boundary_condition = 'constant'
dL = 1.0  # grid separation
stencil = (1.0 / (12.0 * dL * dL)) * np.array(
    [[0, 0, -1, 0, 0],
     [0, 0, 16, 0, 0],
     [-1, 16, -60, 16, -1],
     [0, 0, 16, 0, 0],
     [0, 0, -1, 0, 0]])


L = np.array([40, 100])  # simulation size
L0 = np.array([0, 0]) 

neuron_position = np.array([20, 20])
ngf_source_position = (30, 90)
ngf_source_value = 10.0
mf_source_position = (neuron_position[0], neuron_position[1])
mf_source_value = 10.0
mtb_source = 1.0
mtb_positions = []
lambda_mtb = 0.1
neuron_radius = 15
gc_radius = 3  # growth cone radius
gc_position = np.array([neuron_position[0], neuron_position[1] + neuron_radius ])

small_box_size = 15


### mRNA
kf = 1.0  # mRNA-free diffusion coefficient 
kl = 0.5  # mRNA-linked diffusion coefficient 
lambda_f = 0.0001  # mRNA-free decay
lambda_l = 0.0001  # mRNA-linked decay
gamma = 10.0  # reaction rate mRNA free to linked
beta_everywhere = 0.001  # reaction rate mRNA linked -> free
beta_growth_cone = 10.0  # at the growth cone 
beta_m = np.zeros(L)
Mm = 1.  # microtubules saturation
chi_ml = 100.0  # mRNA-Linked transport coefficient
### 

chi = 10.  # chemotaxis
alpha_p = 10.0  # proliferation
D_ngf = 100.  # Nerve Growth Factor diffusion coefficient
D_mtb = 0.1  # microtubules diffusion coefficient


nstep = 0  # actual time step
tstep = 60000  # total number of steps
dt = 0.001  # time interval
print_period = 10000
nprint = print_period # counter to output data