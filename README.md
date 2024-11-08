# $\mu-\Sigma$ Emulator with Screening

Emulated $\mu-\Sigma$ non-linear and linear power
spectra for fast weak lensing analysis.

The cosmosis-folder contains the required module for power spectrum computation, as well as examplary ini-files.

## Requirements
Required python packages:
* numpy
* scipy
* ...
* [HMcode2020Emu](https://github.com/MariaTsedrik/HMcode2020Emu)


## Parameter ranges
| parameter         | limits                |
| :---:             | :---:                 |
| omega_cdm         | [0.2, 0.6]            |
| omega_baryon      | [0.03, 0.07]          |
| hubble            | [0.58, 0.8]           |
| As                | [1.5e-9, 2.5e-9]      |
| ns                | [0.93, 1.]            |
| neutrino_mass [eV]| [0., 0.6]             |
| mu0               | [-1., 3.]             |
| sigma0            | [-1., 3.]             |
| q1                | [-2., 2.]             |
| log10TAGN         | [7.6, 8.3]            |
| z                 | [0.0, 2.5]            |
| :---:             | :---:                 |
| k_lin [h/Mpc]     | [1e-4, 50]            |
| k_nonlin [h/Mpc]  | [0.01, 10]            |

Note that we impose the following condition: $\mu_0 <= ( 2 \Sigma_0 + 1 )$. One could add it on the level of priors into analysis too.



## Citation

If you use this emulator at any point in your work please cite [this paper](https://arxiv.org/abs/2404.11508) where the similar [data production and emulation pipelines](https://github.com/nebblu/ReACT-emus?tab=readme-ov-file) have been used:

    @article{Tsedrik2024,
        author = Tsedrik, Maria and Bose, Benjamin and Carrilho, Pedro and Pourtsidou, Alkistis and Pamuk, Sefa and Casas, Santiago and Lesgourgues, Julien,
        title = {Stage-IV Cosmic Shear with Modified Gravity and Model-independent Screening},
        journal = {JCAP},
        year = {2024},
        month = {Oct},
        volume = {2024},
        number = {10},
        eid = {099},
        pages = {099},
        doi = {10.1088/1475-7516/2024/10/099},
        eprint = {2404.11508},
        archivePrefix = {arXiv}
    }

    


## License
[MIT](https://choosealicense.com/licenses/mit/)
