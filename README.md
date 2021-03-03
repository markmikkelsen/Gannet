# Gannet 3.1

Open-source, MATLAB-based software for automated data processing and quantification of edited magnetic resonance spectroscopy (MRS) data.

Full software documentation can be found [here](https://markmikkelsen.github.io/Gannet-docs/index.html). Visit our [blog](http://www.gabamrs.com/) for the latest news on Gannet and our developments in edited MRS methodology.

## Getting Started

### Prerequisites

Gannet runs in MATLAB (we recommend using the latest release if possible). Additionally, Gannet requires that the following MATLAB toolboxes are installed:

* Image Processing
* Optimization
* Signal Processing
* Statistics and Machine Learning

To run the voxel co-registration and structural image segmentation modules, [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) must be installed.

### Installing

Please visit the [documentation website](https://markmikkelsen.github.io/Gannet-docs/index.html) for full installation instructions.

## Getting help

If you encounter any problems, please first check the [documentation website](https://markmikkelsen.github.io/Gannet-docs/index.html) for a solution.

For bug reporting, contribution requests or queries, please contact us using our blog's [contact form](http://www.gabamrs.com/contact).

## Versioning

[Semantic versioning](https://semver.org/) is used when updates are made to Gannet using the style 'x.x.x'. Versioning is also conducted on a module-specific basis using the style 'YYMMDD'. That is, each Gannet module has its own release version.

## Developers

* **Richard Edden** (Johns Hopkins University)
* **Mark Mikkelsen** (Johns Hopkins University)
* **Georg Oeltzschner** (Johns Hopkins University)
* **Muhammad Saleh** (Johns Hopkins University)
* **C. John Evans** (Cardiff University)
* **Ashley Harris** (University of Calgary)
* **Nicolaas Puts** (King's College London)

## License and citing Gannet

This software is licensed under the open-source [BSD 3-Clause License](https://markmikkelsen.github.io/Gannet-docs/gannet-license.html). Should you disseminate material that made use of Gannet, please cite the following publications, as appropriate:

**For all work using Gannet:**

* Edden RAE, Puts NAJ, Harris AD, Barker PB, Evans CJ. [Gannet: A batch-processing tool for the quantitative analysis of gamma-aminobutyric acid-edited MR spectroscopy spectra.](https://doi.org/10.1002/jmri.24478) _Journal of Magnetic Resonance Imaging_. 2014;40(6):1445–1452

**If you perform frequency-and-phase correction (FPC) using:**

Robust spectral registration (`RobustSpecReg`):

* Mikkelsen M, Tapper S, Near J, Mostofsky SH, Puts NAJ, Edden RAE. [Correcting frequency and phase offsets in MRS data using robust spectral registration.](https://doi.org/10.1002/nbm.4368) _NMR in Biomedicine_. 2020:e4368

multi-step FPC (`SpecRegHERMES`):

* Mikkelsen M, Saleh MG, Near J, et al. [Frequency and phase correction for multiplexed edited MRS of GABA and glutathione.](https://doi.org/10.1002/mrm.27027) _Magnetic Resonance in Medicine_. 2018;80(1):21-28

or spectral registration (`SpecReg`):

* Near J, Edden R, Evans CJ, Paquin R, Harris A, Jezzard P. [Frequency and phase drift correction of magnetic resonance spectroscopy data by spectral registration in the time domain.](https://doi.org/10.1002/mrm.25094) _Magnetic Resonance in Medicine_. 2015;73(1):44-50

**If you perform tissue segmentation:**

* Ashburner J, Friston KJ. [Unified segmentation.](https://doi.org/10.1016/j.neuroimage.2005.02.018) _NeuroImage_. 2005;26(3):839–851

**If you report water-referenced measurements tissue-corrected using:**

The Harris et al. method:

* Harris AD, Puts NAJ, Edden RAE. [Tissue correction for GABA-edited MRS: Considerations of voxel composition, tissue segmentation, and tissue relaxations.](https://doi.org/10.1002/jmri.24903) _Journal of Magnetic Resonance Imaging_. 2015;42(5):1431–1440

or the Gasparovic et al. method:

* Gasparovic C, Song T, Devier D, et al. [Use of tissue water as a concentration reference for proton spectroscopic imaging.](https://doi.org/10.1002/mrm.20901) _Magnetic Resonance in Medicine_. 2006;55(6):1219–1226

## Acknowledgments

The development and dissemination of Gannet has been supported by the following NIH grants:

* R01 EB016089
* R01	EB023963
* P41 EB015909
* R01	MH106564
* R21 NS077300
* R21 MH098228
* R01 MH096263

We wish to thank the following individuals for their direct or indirect contributions:

* Yair Altman (Undocumented Matlab)
* Peter Barker (Johns Hopkins University)
* Alex Craven (University of Bergen)
* Philipp Ehses (Max Planck Institute for Biological Cybernetics)
* Robin de Graaf (Yale School of Medicine)
* Ralph Noeske (GE Healthcare)
* Jamie Near (McGill University)
* Wouter Potters (UMC Amsterdam)
