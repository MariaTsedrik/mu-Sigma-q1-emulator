from cosmosis.datablock import names, option_section
cosmo = names.cosmological_parameters

"""
Skip point in parameter space if mu > (2*sigma+1) . 
"""

def setup(options):
    verbose = options.get_bool(option_section, "verbose", default=False)
    return {'verbose':verbose}

def execute(block, config):

    mu0 = block[cosmo, "mu0"]
    sigma0 = block[cosmo, "sigma0"]
    
    if mu0 > (2.*sigma0+1.):
        return 1

    return 0


def cleanup(config):
    pass