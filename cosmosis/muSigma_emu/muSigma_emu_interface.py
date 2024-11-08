from cosmosis.datablock import names, option_section as opt
import cosmopower as cp
from cosmopower import cosmopower_NN
import traceback
import numpy as np
import warnings
from scipy import interpolate
import HMcode2020Emu as hmcodeemu #version 3 of the emulator from github https://github.com/MariaTsedrik/HMcode2020Emu/tree/main


from scipy.integrate import trapz,simpson, quad
emulator_nonlin = cp.cosmopower_NN(restore=True,
                        restore_filename='/home/mtsedrik/cosmosis-standard-library/structure/muSigma_emu/models/MuSigma_nonlinear_log10ps_v2',
                        )
emu_k_nonlin = emulator_nonlin.modes
emulator_lin = cp.cosmopower_NN(restore=True,
                        restore_filename='/home/mtsedrik/cosmosis-standard-library/structure/muSigma_emu/models/MuSigma_linear_log10ps',
                        )
emu_k_lin = emulator_lin.modes
print("Loaded nonlinear mu0-Sigma0 model")
cosmo = names.cosmological_parameters
hmpar = names.halo_model_parameters
kmin_lin = emu_k_lin[0]
kmin_nonlin = emu_k_nonlin[0]
kmax_nonlin = emu_k_nonlin[-1]

hmcode_emulator = hmcodeemu.Matter_powerspectrum()
emulator_bar_boost = cp.cosmopower_NN(restore=True, 
                      restore_filename='/home/mtsedrik/cosmosis-standard-library/structure/muSigma_emu/models/total_matter_bar_boost_emu',
                      )
emu_k_bar_boost = emulator_bar_boost.modes

def setup(options):
    config = {}
    config['verbose'] = options.get_bool(opt, 'verbose', default=False)
    config['save_s8'] = options.get_bool(opt, 'save_s8', default=False)
    config['baryons'] = options.get_bool(opt, 'baryons', default=False)


    config['zmin'] = options.get_double(opt, 'zmin', default=0.0)#default from emu
    config['zmax'] = options.get_double(opt, 'zmax', default=2.5)#default from emu
    config['nz'] = options.get_int(opt, 'nz', default=150)
    config['zmin_background'] = options.get_double(opt, 'zmin_background', default=config['zmin'])
    config['zmax_background'] = options.get_double(opt, 'zmax_background', default=config['zmax'])
    config['nz_background'] = options.get_int(opt, 'nz_background', default=config['nz'])

    config['kmax'] = options.get_double(opt, "kmax", default=10.)#default from emu
    config['kmin'] = options.get_double(opt, "kmin", default=kmin_lin)
    config['kmax_extrapolate'] = options.get_double(opt, "kmax_extrapolate", default=config['kmax'])
    config['kmin_extrapolate'] = options.get_double(opt, "kmin_extrapolate", default=config['kmin'])
    config['nk'] = options.get_int(opt, "nk", default=200)


    print(f"zmin={config['zmin']}")      
    if config['zmax']>2.5:
        raise ValueError("z >2.5 is out of emulator range!")

    if config['kmin']<kmin_lin:
        warnings.warn(f"k <{kmin_lin} is out of emulator range!")

    if config['kmax']>kmax_nonlin: 
        warnings.warn(f"kmax={config['kmax']} is out of nonlinear emulator range!")


    zmin = config['zmin'] 
    zmax = config['zmax']
    nz = config['nz'] 
    z_arr = np.linspace(zmin, zmax, nz)
    config['z_arr'] = z_arr

    kmax = config['kmax']
    kmin = config['kmin']
    nk = config['nk'] 
    k_arr = np.logspace(np.log10(kmin), np.log10(kmax), nk)
    config['k_arr'] = k_arr

    # Evaluate z on a different grid than the spectra, so we can easily extend it further
    config['z_background'] = np.linspace(
        config["zmin_background"], config["zmax_background"], config["nz_background"])



    return config

def check_params(params):
    ranges =  {
            'Omega_m': [0.2,     0.6], 
            'Omega_b': [0.03,    0.07], 
            #'H0': [58., 80.],
            'h': [0.58, 0.8],
            'ns': [0.93, 1.0],
            'As': [1.5e-09, 2.5e-09],
            'neutrino_mass': [0., 0.6],
            'mu0': [-0.9999, 2.9999],
            'sigma0': [-1., 3.],
            'q1': [-2., 2.]
                        }

    for key, (vmin, vmax) in ranges.items():
        if not (params[key].any() > vmin) and (params[key].any() < vmax):
            raise ValueError(f"emu: {key} must be in range [{vmin},{vmax}]")



def execute(block, config):
    try:
        emulate_power(block, config)
    except ValueError as error:
        return 1
    
    return 0    


def emulate_power(block, config):
    ################################
    # Saving power spectra
    ################################
    nz = config['nz']
    z_arr = config['z_arr']
    k_arr = config['k_arr']
    h = block.get_double(cosmo, "h0")
    Omm = block[cosmo, "omega_m"]
    Omb = block.get_double(cosmo, "omega_b")
    w = -1. #block.get_double(cosmo, "w", default=-1.)
    wa = 0.0 #block.get_double(cosmo, 'wa', default=0.0)
    As = block.get_double(cosmo, "a_s")
    ns = block.get_double(cosmo, 'n_s')
    Mnu = block.get_double(cosmo, 'mnu')
    Omnu = Mnu/(93.14*h**2)
    Omc = Omm-Omb-Omnu
    mu0 = block.get_double(cosmo, 'mu0', default=0.0)
    sigma0 = block.get_double(cosmo, 'sigma0', default=0.0)
    q1 = block.get_double(cosmo, 'q1', default=2.)
    params = {
        "Omega_m": Omm*np.ones(nz),
        "As": As*np.ones(nz),
        "Omega_b": Omb*np.ones(nz),
        "ns": ns*np.ones(nz),
        "H0": h*100.*np.ones(nz),
        "h": h*np.ones(nz),
        "neutrino_mass": Mnu*np.ones(nz),
        "Omega_nu": Omnu*np.ones(nz),
        "mu0": mu0*np.ones(nz),
        "sigma0": sigma0*np.ones(nz),
        "q1": q1*np.ones(nz),
        "z": z_arr,
    }
    # check we're within the allowed bounds for the emulator
    check_params(params)
    # This is the real trustful range for ReACT
    #zmask = z < 2.5
    #kmask = (k < 5.) & (k > 0.01)

    Pk_nl = emulator_nonlin.ten_to_predictions_np(params)
    Pk_l = emulator_lin.ten_to_predictions_np(params)


    pklin_interp = [interpolate.interp1d(emu_k_lin,
                                            pk_lin_i,
                                            kind='linear',
                                            bounds_error=False,
                                            fill_value=(pk_lin_i[0], pk_lin_i[-1]),
                                            ) for pk_lin_i in Pk_l]
    plin = np.array([pklin_interp[i](k_arr) for i in range(nz)])

    pknonlin_interp = [interpolate.interp1d(emu_k_nonlin,
                                                    pk_nonlin_i,
                                                    kind='linear',bounds_error=False,
                                                    fill_value=(pk_nonlin_i[0], pk_nonlin_i[-1]) 
                                                    ) for pk_nonlin_i in Pk_nl]
    pnonlin = np.array([pknonlin_interp[i](k_arr[k_arr>=0.01]) for i in range(nz)])

    
    # concatenation of PkLin, Pk_NL
    pl_left = plin[:, k_arr<emu_k_nonlin[0]]
    pnl  = np.concatenate((pl_left, pnonlin),axis=1)
    k_h = k_arr

    ########################################
    # Adding baryonic feedback if necessary
    ########################################    
    if config['baryons'] == True :
        params_bar_boost = {"omega_cdm": Omc*np.ones(nz),
                "As": As*np.ones(nz),
                "omega_baryon": Omb*np.ones(nz),
                "ns": ns*np.ones(nz),
                "hubble": h*np.ones(nz),
                "neutrino_mass": Mnu*np.ones(nz),
                "w0": -1.*np.ones(nz),
                "wa": np.zeros(nz),
                "log10TAGN": block.get_double(hmpar, "logt_agn", default=7.7)*np.ones(nz),
                "z": z_arr,
            }
        
        #_, Pnl = hmcode_emulator.get_nonlinear_pk(**params_bar_boost, baryonic_boost=True, no_nu=False)
        #print('Pnl', Pnl.shape)
        barboost = emulator_bar_boost.predictions_np(params_bar_boost) 
        barboost_interp = [interpolate.interp1d(emu_k_bar_boost,
                                                    barboost_i,
                                                    kind='linear',bounds_error=False,
                                                    fill_value=(barboost_i[0], barboost_i[-1]) 
                                                    ) for barboost_i in barboost]
        bar_boost = np.array([barboost_interp[i](k_arr) for i in range(nz)])
        pnl = pnl*bar_boost

    block.put_grid("matter_power_lin", "z", z_arr, "k_h", k_h, "P_k", plin)
    block.put_grid("matter_power_nl", "z",  z_arr, "k_h", k_h, "P_k", pnl)
    
    ################################
    # Saving background
    ################################
    z_background = config['z_background']
    # Write basic distances and related quantities to datablock
    block[names.distances, "nz"] = len(z_background)
    block[names.distances, "z"] = z_background
    block[names.distances, "a"] = 1/(z_background+1)

    H0 = h/2997.92458 #in 1/Mpc
    omegaL_func = lambda z: (1.-Omm) * pow(1.+z, 3.*(1.+w+wa)) * np.exp(-3.*wa*z/(1.+z))
    E_z_func = lambda z: np.sqrt(Omm*pow(1.+z, 3) + omegaL_func(z))
    E_z_grid = np.array([E_z_func(zz_i) for zz_i in z_background])
    r_z_int = lambda z: 1./np.sqrt(Omm*pow(1.+z, 3) + omegaL_func(z))
    r_z_func = lambda z_in: quad(r_z_int, 0, z_in)[0]
    r_z_grid = np.array([r_z_func(zz_i) for zz_i in z_background])/H0 #Mpc

    D_C = r_z_grid #in Mpc
    H = H0*E_z_grid #in 1/Mpc
    #D_H = 1 / H[0]
    D_M = D_C
    D_L = D_M * (1 + z_background)
    D_A = D_M / (1 + z_background)
    #D_V = ((1 + z_background)**2 * z_background * D_A**2 / H)**(1./3.)


    block[names.distances, "D_C"] = D_C
    block[names.distances, "D_M"] = D_M
    block[names.distances, "D_L"] = D_L
    block[names.distances, "D_A"] = D_A
    #block[names.distances, "D_V"] = D_V
    block[names.distances, "H"] = H

    ################################
    # Saving sigma_8 if necessary
    # (would recommend to keep at False 
    # and compute S8 from the chain 
    # after convergence) 
    ################################
    if config['save_s8'] == True :
        if z_arr[0]==0:
            params_s8 = {"omega_cdm": [Omc],
                "As": [As],
                "omega_baryon": [Omb],
                "ns": [ns],
                "hubble": [h],
                "neutrino_mass":[Mnu],
                "w0": [-1.],
                "wa": [0.],
                "z": [0.],
            }
            sigma_8_GR, _ = hmcode_emulator.get_sigma8(**params_s8)
            sigma_8_GR = sigma_8_GR[0]
            #compute linear growth from power spectra
            params = {
                'omega_cdm'     :  [Omc],
                'As'            :  [As],
                'omega_baryon'  :  [Omb],
                'ns'            :  [ns],
                'hubble'        :  [h],
                'neutrino_mass' :  [Mnu],
                'w0'            :  [-1.0],
                'wa'            :  [0.0],
                'z'             :  [0.]
            }
            _, pklin_lcdm = hmcode_emulator.get_linear_pk(**params, k=np.array([0.01]))
            D=np.sqrt(pklin_interp[0](0.01)/pklin_lcdm[0][0])
            #print('D', D)
            sigma_8 = D*sigma_8_GR
            block[cosmo, "sigma_8"] = sigma_8
            block[cosmo, "S_8"] = sigma_8*np.sqrt(Omm/0.3)
            #print("sigma_8", sigma_8)
            #print("S_8", sigma_8*np.sqrt(Omm/0.3))
        else:
            warnings.warn(f"z_min != 0 is not implemented yet")

    
    return 0     


if __name__=="__main__":
    print("Executing example case")
