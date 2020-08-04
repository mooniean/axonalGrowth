from scipy.ndimage import convolve
import numpy as np
from parameters import *
import os

prefix_ = 'deb_2_'
os.system('mkdir '+prefix_)

ngf = np.zeros(L)
phi = np.zeros(L)
psi = np.zeros(L)
mtb = np.zeros(L)
mf = np.zeros(L)
ml = np.zeros(L)
v_m = np.array([np.zeros(L),np.zeros(L)] )

mf = 1.0*init_cell(mf, neuron_position, neuron_radius-5)
ngf = init_growth_factor(ngf, neuron_position, neuron_radius )
phi = init_cell(phi, gc_position, gc_radius)
psi = init_cell(psi, neuron_position, neuron_radius)


# streamplot
y=np.linspace(0,40,40)
x=np.linspace(0,100,100)


phi = smoothing(phi, 2000, stencil)  # growth cone
psi = smoothing(psi, 2000, stencil)  # neuron

for i in range(0, 100):
    mf = mf + dt * convolve(10 * mf + container(mf, psi, 0.0001), stencil, mode=boundary_condition )
    ngf = ngf + dt * D_ngf * convolve(ngf, stencil, mode=boundary_condition )

V_target = h(phi).sum()


r_cm = com(phi, V_target, ndim=2)
icalc = 0
ncalc = 250

while (nstep <= tstep): 

    v = h(phi).sum()
    Lagrange = np.ones(L) * (V_target - v)

    il, icl = np.int(r_cm[0]) - small_box_size, np.int(r_cm[0]) - gc_radius
    ih, ich = np.int(r_cm[0]) + small_box_size, np.int(r_cm[0]) + gc_radius
    jl, jcl = np.int(r_cm[1]) - small_box_size, np.int(r_cm[1]) - gc_radius
    jh, jch = np.int(r_cm[1]) + small_box_size, np.int(r_cm[1]) + gc_radius
     
    grad_phi = np.gradient(phi[il:ih, jl:jh])
    grad_ngf = np.gradient(ngf[il:ih, jl:jh])
    grad_ngf = unit_vector(grad_ngf, norm = np.linalg.norm(grad_ngf, axis=0))
    
    coeff_ml = (Mm-ml)/Mm
    beta_m[:,:] = beta_everywhere
    beta_m[icl:ich, jcl:jch] = beta_growth_cone
    coeff_beta = coeff_beta_(mtb)*beta_m*(1-mtb)*ml
    coeff_gamma = gamma*mf*mtb*coeff_ml  
    grad_ml = np.gradient(ml)

    # change the microtubules
    advection_ml = coeff_ml*mtb*(v_m[0] * grad_ml[0] + v_m[1] * grad_ml[1])
                   # + ml*convolve(mtb, stencil, mode=boundary_condition) )

    
    grad_ngf_norm = np.linalg.norm(grad_ngf, axis=0, keepdims=True)
    unit_grad_ngf = unit_vector(grad_ngf, grad_ngf_norm)
    chem = ((grad_phi[0]*unit_grad_ngf[0] + grad_phi[1]*unit_grad_ngf[1]))
    mf_GC = mf[icl:ich, jcl:jch].sum()
    
    # troca (mf_GC/(1+mf_GC)) para apenas mf_GC e talvez ajustar chi
    chi_ = chi * (mf_GC/(1+mf_GC)) * chem
    
    mf.itemset(mf_source_position, mf_source_value ) #  source of free-mRNA at the centre of the axon
    ngf.itemset(ngf_source_position, ngf_source_value)
    #ngf[:, 80] = 10.0 # for a line of sources
    
    np.put(mtb, mtb_positions, mtb_source) # setting the center of the axon as mtb_source

  
    phi[il:ih, jl:jh], psi,  mtb, ngf, mf, ml = phi[il:ih, jl:jh] + (dt * (-chi_  + convolve(phi[il:ih, jl:jh], stencil, mode='constant') +
                phi[il:ih, jl:jh] * (1.0 - phi[il:ih, jl:jh]) * ( (phi[il:ih, jl:jh] - 0.5) + Lagrange[il:ih, jl:jh]))), \
                psi + (dt * (-convolve(8 * (psi * (1.0 - psi) * (psi - 0.5))    
                       + convolve(psi, stencil, mode=boundary_condition),
                                      stencil, mode=boundary_condition) +
                                      alpha_p * phi * ngf * psi * (1.0 - psi) )), \
                mtb + (dt * (D_mtb * convolve(mtb + 2*container(mtb, psi, alpha_), stencil, mode=boundary_condition) - lambda_mtb*mtb)),\
                ngf + (dt * (D_ngf * convolve(ngf, stencil, mode=boundary_condition) - psi * ngf )),\
                mf + dt * ( convolve(kf * mf + 2*container(mf, psi, alpha_), stencil, mode=boundary_condition) 
                            - lambda_f * mf 
                            - coeff_gamma
                            + coeff_beta),\
                ml + dt * ( convolve(kl * ml + 2*container(ml, psi, alpha_), stencil, mode=boundary_condition) 
                            - lambda_l * ml  
                            + coeff_gamma
                            - coeff_beta 
                            + chi_ml*advection_ml)
    
     
    if icalc >= ncalc:
        
        icalc = 0
        r_cm = com(phi, v, 2)
        r_cm = (np.int(r_cm[0]), np.int(r_cm[1]))
        mtb_positions.append(np.ravel_multi_index(r_cm, mtb.shape))
        positions_mtb = np.unique(np.asarray(mtb_positions), axis=0)
        mtb_diff = np.diff(np.unravel_index(positions_mtb, mtb.shape)).transpose()

        try:
            # +3 - 2
            v_m[0, icl:ich, jcl+3:jch-2] = mtb_diff[-1, 1]
            v_m[1, icl:ich, jcl+3:jch-2] = mtb_diff[-1, 0]
        except:
            pass
        
    if (nprint >= print_period):
        print(nstep)
        np.savez(prefix_+'/'+prefix_+str(nstep)+'.npz', psi=psi, mtb=mtb, v_m=v_m, ml=ml, mf=mf, phi=phi, ngf=ngf)
        nprint = 0
    icalc += 1
    nstep += 1
    nprint += 1
