"""

    this files contains the functions used in the main routine

"""

from scipy.ndimage import convolve
import numpy as np

def heaviside(x):
    if x>=0:
        return 1
    else:
        return 0
heaviside = np.vectorize(heaviside)

h = lambda x: x ** 2 * (3. - 2. * x)

def smoothing(phi, tstep, stencil):
    nstep = 0
    V_target = h(phi).sum()
    dt = 0.001
    while (nstep <= tstep):  # this loop is integrating the Allen-Cahn function in time

        v = h(phi).sum()
        Lagrange = np.ones(phi.shape) * (V_target - v)
        phi = phi + dt * (convolve(phi, stencil, mode='wrap') + phi * (1.0 - phi) * (
                phi - 0.5 + Lagrange))
        nstep += 1

    return phi


def com(phi, volume, ndim):
    r_cm = np.zeros(ndim)
    volume = phi.sum()
    for i in range(0, phi.size):
        s = np.asarray(np.unravel_index(i, phi.shape))
        r_cm = r_cm + s*phi.item(i)

    return r_cm/volume

def interface_pos(phi):
    pos = 0
    for i in range(0, phi.size):
        pos = pos + i*phi.item(i)

    return pos/phi.sum()

def init_cell(phi, cell_position, cell_radius):
    """
    this function initializes a round cell
    :param phi: 2d array
    :param cell_radius: size of the cell
    :return phi: 2d array
    """
    for i in range(0, phi.size):
        s = np.asarray(np.unravel_index(i, phi.shape))

        distance = s - cell_position
        distance = np.sqrt(np.sum(np.power(distance, 2)))
        if distance <= cell_radius:
            phi.itemset(i, 1.0)
    return phi

def init_growth_factor(ngf, neuron_position, neuron_radius):
    for i in range(0, ngf.size):
        s = np.asarray(np.unravel_index(i, ngf.shape))

        distance = s - neuron_position
        distance = np.sqrt(np.sum(np.power(distance, 2)))
        if distance > neuron_radius+2:
            ngf.itemset(i, np.random.rand()*10.)

    return ngf

def container(field_1, field_2, alpha):
    """

    :param field_1: chemical field 1
    :param field_2: chemical field 2
    :return: repulsion between field_1 and field_2
    """
    mu = np.power(1. - field_2, 2) + alpha
    return field_1*mu* np.power(mu*(np.power(field_1, 2) + alpha), -0.5)

def unit_vector(vector, norm):
    
    if norm > 0:
        return vector/norm
    else:
        return 0
unit_vector = np.vectorize(unit_vector)

def coeff_beta_(mtb):
    if mtb>10e-2:
        return  1./mtb
    else:
        return 0.
coeff_beta_ = np.vectorize(coeff_beta_)