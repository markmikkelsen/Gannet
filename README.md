# Gannet

[![GitHub release (latest by date)](https://img.shields.io/github/v/release/markmikkelsen/Gannet)](https://github.com/markmikkelsen/Gannet/releases)
![GitHub commits since latest release (by date) for a branch](https://img.shields.io/github/commits-since/markmikkelsen/Gannet/latest/main)
[![Website](https://img.shields.io/website?label=website&up_message=documentation&url=https%3A%2F%2Fmarkmikkelsen.github.io%2FGannet-docs%2Findex.html)](https://markmikkelsen.github.io/Gannet-docs/index.html)
[![Website](https://img.shields.io/website?label=website&up_message=forum&url=https://forum.mrshub.org/c/mrs-software/gannet/9)](https://forum.mrshub.org/c/mrs-software/gannet/9)
[![GitHub](https://img.shields.io/github/license/markmikkelsen/Gannet)](https://github.com/markmikkelsen/Gannet/blob/main/LICENSE)
[![DOI:10.1002/jmri.24478](https://img.shields.io/badge/DOI-10.1002/jmri.24478-blue.svg)](https://doi.org/10.1002/jmri.24478)
[![Citation Badge](https://api.juleskreuer.eu/citation-badge.php?doi=10.1002/jmri.24478)](https://juleskreuer.eu/projekte/citation-badge/)

<img src="Gannet3_logo.png" alt="Gannet logo" width="150">

Open-source, MATLAB-based software for automated data processing and quantification of edited <sup>1</sup>H magnetic resonance spectroscopy (MRS) data.

Full software documentation can be found [here](https://markmikkelsen.github.io/Gannet-docs/index.html).

## Overview

Gannet is a free, open-source MATLAB-based software toolkit for analyzing edited single-voxel <sup>1</sup>H magnetic resonance spectroscopy (MRS) data. Its largely automated functions cover all the essential steps of modern MRS analysis:

- Loading raw data
- Substantial preprocessing
- Signal modeling
- Voxel co-registration with structural MR images
- Concentration estimation based on tissue composition

Several existing software packages for MRS data analysis require substantial user input or offer a wide selection of processing options. In contrast, the philosophy behind Gannet is to provide users with a complete automated pipeline without the need for significant user input.

Additionally, as open-source software, advanced users have the ability to modify the underlying routines for ad hoc purposes.

## Installation

### Prerequisites  

Gannet runs in [MATLAB](https://mathworks.com/products/matlab.html). For best performance, we recommend using the latest release if possible. Additionally, Gannet requires that the following MATLAB toolboxes are installed:

- Image Processing
- Optimization
- Signal Processing
- Statistics and Machine Learning

You can check which toolboxes you have installed by typing `ver` in the MATLAB command window. To install any missing toolboxes, please follow these [instructions](https://mathworks.com/matlabcentral/answers/101885-how-do-i-install-additional-toolboxes-into-an-existing-installation-of-matlab).

To run the voxel co-registration and structural image segmentation modules, [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) must be installed.

### Download

The simplest way to install Gannet is to download the code from the <a href="https://github.com/markmikkelsen/Gannet" target="_blank">GitHub repository</a> and move the `Gannet-main/` directory into your MATLAB directory.

Alternatively, Git users can clone the Gannet repository into a directory of their choice:

`git clone https://github.com/markmikkelsen/Gannet.git`

The development version can be downloaded from the <a href="https://github.com/markmikkelsen/Gannet/tree/dev" target="_blank">development branch</a> on GitHub or by using the following git command if the repository was cloned:

`git checkout dev`

Stable releases can be found <a href="https://github.com/markmikkelsen/Gannet/releases" target="_blank">here</a>.

### Setup

Open the *Set Path* dialog box from the MATLAB menu (or run the command `pathtool` in the Command Window), click *Add with Subfolders...*, find the downloaded Gannet directory and then select it. When done, press *Save* to permanently save the Gannet directory to MATLAB's default search path.

SPM12 can be installed in the same manner after it has been downloaded from the [SPM](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) website.

**It is highly recommended that you only add the main SPM12 directory (`spm12/`) to your search path instead of including all the subdirectories. This prevents function conflicts.**

**If you have <a href="https://github.com/schorschinho/osprey" target="_blank">Osprey</a> also installed on your computer, please ensure that you have either all the Gannet and SPM12 directories at the top of your search path or removed the Osprey directories from your search path. Gannet and Osprey share several functions that can lead to conflicts that result in Gannet not functioning correctly.**

## Compatibility

Gannet is currently being developed in MATLAB R2023b in macOS 14 Sonoma (Apple silicon). While reasonable effort is made to ensure legacy and cross-OS compatibility, an error-free user experience is not guaranteed.

## Supported file formats

At present, the following MRS data file formats are supported:

- DICOM (.dcm)
- GE P-file (.7)
- NIfTI-MRS (.nii[.gz])
- Philips .data/.list
- Philips .raw
- Philips .sdat/.spar
- Siemens DICOM (.ima)
- Siemens .rda
- Siemens TWIX (.dat)

For creating and co-registering voxel masks, structural images need to be in NIfTI format (DICOM files can also be used if processing GE P-files).

**Philips users: Do not use structural images exported using the *fsl-nifti* option as this creates problems with co-registration in Gannet.**

## Getting help

If you encounter any problems, please first check the [documentation website](https://markmikkelsen.github.io/Gannet-docs/index.html) or the [FAQ](https://markmikkelsen.github.io/Gannet-docs/faq.html) page for a solution.

Otherwise, you can post your query on the [Gannet forum](https://forum.mrshub.org/c/mrs-software/gannet/9) on the [MRSHub](https://mrshub.org/).

The Gannet team can also be <a href="mailto:mam4041@med.cornell.edu?subject=[Gannet]%20<Please enter the subject of your query here>" target="_blank">contacted</a> directly. We will do our best to work with you to solve your issue.

## Versioning

Gannet uses a form of [semantic versioning](https://semver.org/) in the style 'x.x.x' to mark code releases. Versioning is also conducted on a module-specific basis using the style 'YYMMDD'. That is, each Gannet module has its own release version. Users should note that module-specific versions sometimes are updated despite the semantic version number remaining unchanged (typically for minor updates/bug fixes).

## Developers

- Richard Edden (Johns Hopkins University) - creator
- Mark Mikkelsen (Weill Cornell Medicine) - lead developer
- Georg Oeltzschner (Johns Hopkins University) - contributor
- Muhammad Saleh (Children's Hospital of Philadelphia) - contributor
- C. John Evans (Cardiff University) - contributor
- Ashley Harris (University of Calgary) - contributor
- Nicolaas Puts (King's College London) - contributor

## License and citing Gannet

This software is licensed under the open-source [BSD-3-Clause License](https://github.com/markmikkelsen/Gannet/blob/main/LICENSE). Should you disseminate material that made use of Gannet, please cite the following publications, as appropriate:

- Edden RAE, Puts NAJ, Harris AD, Barker PB, Evans CJ. [Gannet: A batch-processing tool for the quantitative analysis of gamma-aminobutyric acid-edited MR spectroscopy spectra.](https://doi.org/10.1002/jmri.24478) *Journal of Magnetic Resonance Imaging*. 2014;40(6):1445–1452

**If you perform frequency-and-phase correction (FPC) using:**

Robust spectral registration (`RobustSpecReg`):

- Mikkelsen M, Tapper S, Near J, Mostofsky SH, Puts NAJ, Edden RAE. [Correcting frequency and phase offsets in MRS data using robust spectral registration.](https://doi.org/10.1002/nbm.4368) *NMR in Biomedicine*. 2020;33(10):e4368

multi-step FPC (`SpecRegHERMES`):

- Mikkelsen M, Saleh MG, Near J, et al. [Frequency and phase correction for multiplexed edited MRS of GABA and glutathione.](https://doi.org/10.1002/mrm.27027) *Magnetic Resonance in Medicine*. 2018;80(1):21-28

or spectral registration (`SpecReg`):

- Near J, Edden R, Evans CJ, Paquin R, Harris A, Jezzard P. [Frequency and phase drift correction of magnetic resonance spectroscopy data by spectral registration in the time domain.](https://doi.org/10.1002/mrm.25094) *Magnetic Resonance in Medicine*. 2015;73(1):44-50

**If you perform tissue segmentation:**

- Ashburner J, Friston KJ. [Unified segmentation.](https://doi.org/10.1016/j.neuroimage.2005.02.018) *NeuroImage*. 2005;26(3):839–851

**If you report water-referenced, tissue-corrected metabolite measurements using:**

The Harris et al. method:

- Harris AD, Puts NAJ, Edden RAE. [Tissue correction for GABA-edited MRS: Considerations of voxel composition, tissue segmentation, and tissue relaxations.](https://doi.org/10.1002/jmri.24903) *Journal of Magnetic Resonance Imaging*. 2015;42(5):1431–1440

or the Gasparovic et al. method:

- Gasparovic C, Song T, Devier D, et al. [Use of tissue water as a concentration reference for proton spectroscopic imaging.](https://doi.org/10.1002/mrm.20901) *Magnetic Resonance in Medicine*. 2006;55(6):1219–1226

## Acknowledgments

The development and dissemination of Gannet has been supported by the following NIH grants:

- R01 EB016089
- R01 EB023963
- P41 EB015909
- K99 EB028828
- R01 MH106564
- R21 MH098228
- R21 NS077300
- R01 MH096263

We wish to thank the following individuals for their direct or indirect contributions:

- Yair Altman (Undocumented Matlab)
- Peter Barker (Johns Hopkins University)
- Alex Craven (University of Bergen)
- Philipp Ehses (Max Planck Institute for Biological Cybernetics)
- Robin de Graaf (Yale School of Medicine)
- Xiangrui Li (Ohio State University)
- Jamie Near (McGill University)
- Ralph Noeske (GE Healthcare)
- Wouter Potters (UMC Amsterdam)
- Jan Simon (Heidelberg)
