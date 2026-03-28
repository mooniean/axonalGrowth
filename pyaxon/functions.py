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
            phi.flat[i] = 1.0
    return phi

def init_growth_factor(ngf, neuron_position, neuron_radius):
    for i in range(0, ngf.size):
        s = np.asarray(np.unravel_index(i, ngf.shape))

        distance = s - neuron_position
        distance = np.sqrt(np.sum(np.power(distance, 2)))
        if distance > neuron_radius+2:
            ngf.flat[i] = np.random.rand()*10.

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


def rasterized_line(start, end):
    """Return integer grid points along the segment connecting two points."""
    start = np.asarray(np.rint(start), dtype=int)
    end = np.asarray(np.rint(end), dtype=int)
    nsteps = int(np.max(np.abs(end - start))) + 1

    if nsteps <= 1:
        return np.asarray([start], dtype=int)

    line = np.rint(np.linspace(start, end, nsteps)).astype(int)
    keep = np.ones(line.shape[0], dtype=bool)
    keep[1:] = np.any(np.diff(line, axis=0) != 0, axis=1)
    return line[keep]


def build_transport_field(track_points, shape, support_radius=1):
    """Build a motor-velocity field aligned with an ordered microtubule track."""
    velocity = np.zeros((2,) + tuple(shape), dtype=float)
    weights = np.zeros(shape, dtype=float)

    if len(track_points) < 2:
        return velocity

    track_points = np.asarray(track_points, dtype=int)
    max_i, max_j = shape[0] - 1, shape[1] - 1

    for current, nxt in zip(track_points[:-1], track_points[1:]):
        delta = nxt - current
        norm = np.linalg.norm(delta)
        if norm == 0:
            continue

        direction = delta / norm
        for di in range(-support_radius, support_radius + 1):
            for dj in range(-support_radius, support_radius + 1):
                i = current[0] + di
                j = current[1] + dj
                if 0 <= i <= max_i and 0 <= j <= max_j:
                    velocity[:, i, j] += direction
                    weights[i, j] += 1.0

    mask = weights > 0
    velocity[0, mask] /= weights[mask]
    velocity[1, mask] /= weights[mask]
    return velocity


def conservative_upwind_advection(field, velocity, dL=1.0):
    """Return a first-order upwind approximation of ``-div(field * velocity)``."""
    advection = np.zeros_like(field, dtype=float)
    v0, v1 = velocity

    face_v0 = 0.5 * (v0[:-1, :] + v0[1:, :])
    flux0 = np.where(face_v0 >= 0.0, face_v0 * field[:-1, :], face_v0 * field[1:, :])
    advection[:-1, :] -= flux0 / dL
    advection[1:, :] += flux0 / dL

    face_v1 = 0.5 * (v1[:, :-1] + v1[:, 1:])
    flux1 = np.where(face_v1 >= 0.0, face_v1 * field[:, :-1], face_v1 * field[:, 1:])
    advection[:, :-1] -= flux1 / dL
    advection[:, 1:] += flux1 / dL
    return advection

