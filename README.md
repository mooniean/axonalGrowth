![pyaxon](https://github.com/mooniean/axonalGrowth/workflows/pyaxon/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/axonalgrowth/badge/?version=latest)](https://axonalgrowth.readthedocs.io/en/latest/?badge=latest)
[![Python 3](https://pyup.io/repos/github/mooniean/axonalGrowth/python-3-shield.svg)](https://pyup.io/repos/github/mooniean/axonalGrowth)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# axonalGrowth
Developing a phase-field model for axonal growth. We have interest on how mRNA traffic affects the growth cone navigation.


# Phase-field model for axonal growth
The growth cone is described by the order parameter $\phi$ as a droplet of radius $R_{GC}$ with velocity $\vec{v}$

![equation](https://latex.codecogs.com/gif.latex?%5Cpartial_t%20%5Cpsi%20%3D%20-M%20%5Cnabla%5E2%20%5Cleft%28%208%20%5Cpsi%20%281-%5Cpsi%29%5Cleft%28%5Cpsi%20-%20%5Cfrac%7B1%7D%7B2%7D%5Cright%29%20&plus;%20%5Cvarepsilon%5E2%20%5Cnabla%5E2%20%5Cpsi%20%5Cright%29%5C%5C%20&plus;%20%5Calpha_p%20%5Cpsi%20%5Cphi%20%281-%5Cphi%29%20%5B%5Cmathrm%7BNGF%7D%5D%20%5C%2C%20%2C)

where $\alpha_V$ is the Lagrange multiplier with volume target ![equation](https://latex.codecogs.com/gif.latex?V_%7B%5Cmathrm%7Btarget%7D%7D%20%3D%20%5Cfrac%7B4%7D%7B3%7D%20%5Cpi%20R_%7BGC%7D%5E3) is the instantaneous volume of the growth cone. The velocity is written as a function of the concentration of the Nerve Growth Factor (NGF)


![equation](https://latex.codecogs.com/gif.latex?%5Cvec%7Bv%7D_%7B%5Cmathrm%7BGC%7D%7D%20%3D%20%5Cchi%28m_f%29%20%5Cfrac%7B%5Cnabla%5B%5Cmathrm%7BNGF%7D%5D%7D%7B%7C%7C%5Cnabla%5B%5Cmathrm%7BNGF%7D%5D%7C%7C%7D%5C%2C%20%2C)

The neuron and axon are modeled by the Cahn-Hilliard equation plus a phenomenological proliferative term coupled with the growth cone movement

![equation](https://latex.codecogs.com/gif.latex?%5Cpartial_t%20%5Cpsi%20%3D%20-M%20%5Cnabla%5E2%20%5Cleft%28%208%20%5Cpsi%20%281-%5Cpsi%29%5Cleft%28%5Cpsi%20-%20%5Cfrac%7B1%7D%7B2%7D%5Cright%29%20&plus;%20%5Cvarepsilon%5E2%20%5Cnabla%5E2%20%5Cpsi%20%5Cright%29%20&plus;%20%5Calpha_p%20%5Cpsi%20%5Cphi%20%281-%5Cphi%29%20%5B%5Cmathrm%7BNGF%7D%5D%20%5C%2C%20%2C)

where the first term in the right-hand side maintains the interface and the second term is responsible by introducing material at the interface of $\phi$ and $\psi$. This proliferative term depends on the concentration of NGF.



We model the mRNA with two phenotypes: free $m_f$ and linked $m_l$ mRNA. The mRNA dynamic is described by the following reaction-diffusion equations

![equation](https://latex.codecogs.com/gif.latex?%5Cpartial_t%20m_f%20%3D%20k_f%5Cnabla%5E2m_f%20-%20%5Cgamma%20m_f%20m_t%20%5Cleft%28%5Cfrac%7BM_m-m_l%7D%7BM_m%7D%5Cright%29%20-%5Clambda_f%20m_f%20%5C%5C&plus;%20%5Cfrac%7B%5Cbeta_m%20%5Cleft%281-m_t%5Cright%29m_l%7D%7Bm_t%7D%20&plus;%20%5Cnabla%5E2%20%5CXi%5E%7B%5Cnatural%7D_C%5Cleft%5Bm_f%2C%20%5Cpsi%20%5Cright%5D%20%5C%2C%20%2C)

![equation](https://latex.codecogs.com/gif.latex?%5Cpartial_t%20m_l%20%3D%20k_l%20%5Cnabla%5E2%20m_l%20-%5Cnabla%5Ccdot%5Cleft%28%5Cleft%28%5Cfrac%7BM_m-m_l%7D%7BM_m%7D%5Cright%29m_l%20m_t%20%5Cmathbf%7Bv%7D_m%5Cright%29-%5Clambda_l%20m_l%20%5C%5C%20&plus;%20%5Cgamma%20m_f%20m_t%20%5Cleft%28%5Cfrac%7BM_m-m_l%7D%7BM_m%7D%5Cright%29%20-%20%5Cfrac%7B%5Cbeta_m%20%5Cleft%281-m_t%5Cright%29m_l%7D%7Bm_t%7D%20&plus;%20%5Cnabla%5E2%20%5CXi%5E%7B%5Cnatural%7D_C%5Cleft%5Bm_l%2C%20%5Cpsi%20%5Cright%5D%20%5C%2C%20%2C)
