# InSilicoHeartGen
Pipeline to automatically create large virtual patient cohort files to conduct large-scale in silico trials through cardiac electromechanical simulations
(https://doi.org/10.48550/arXiv.2503.03706)


## Overview
**InSilicoHeartGen** is a MATLAB-based computational pipeline for the fully-automated generation of simulation-ready ventricular models and associated input files for large-scale cardiac electrophysiological and electromechanical simulations. Starting from ventricular surface meshes, the pipeline standardises geometry processing, generates volumetric meshes at multiple resolutions, assigns anatomical labels, and computes the most common geometry-dependent fields required for cardiac simulations, such as fibre orientation, transmurality, and electrophysiological heterogeneity.

The pipeline is designed to support the efficient creation of large virtual patient cohorts, enabling in silico studies that would be impractical to construct manually at scale. It integrates a range of established open-source tools and algorithms within a single automated workflow, while remaining flexible with respect to mesh format, resolution, and solver requirements. The generated outputs are solver-ready and compatible with commonly used cardiac simulation environments, including both CPU- and GPU-based solvers.

## Installation

### Requirements
- **MATLAB** version **2021b** or later
- **Windows 10/11** or **Linux** operating system

- NVIDIA GPU is recommended for faster performance.

### Dependencies
To ensure full functionality, the following dependencies must be installed. Copy the required libraries into the `dependencies` folder:

#### MATLAB Toolboxes:
- [Iso2Mesh](https://github.com/fangq/iso2mesh)
- [toolbox_graph](https://github.com/gpeyre/matlab-toolboxes/tree/master/toolbox_graph)
- [gptoolbox](https://github.com/alecjacobson/gptoolbox)
- [Cobiveco](https://github.com/KIT-IBT/Cobiveco)
- [vtkToolbox](https://github.com/KIT-IBT/vtkToolbox) *(Requires additional installation steps; follow the instructions in their README.)*

### Usage

After installing the required dependencies, you can run any of the provided MATLAB scripts to launch the pipeline.

Before execution, make sure to update the scripts to reflect your specific setup:

- Update file paths and input filenames  
- Set the appropriate mesh format: `cut`, `open`, `UKBB` or `closed`  
- Adjust mesh resolution and labelling options as needed  

Common configuration options include:

- `RVseptal_threshold`: Lower this value to increase the number of detected endocardial septal RV faces  
- `Fiber_info`: Controls the fibre angle generated in the field computation  

**Supported Input Formats**

The pipeline supports common surface mesh formats (in cm or mm), including but not limited to:

- `.ply`  
- `.stl`  
- `.vtk`  
- `.vtu`  

**Output**

- VTK raw meshes (in cm)
- VTK surface mesh with labels (`labels_final.vtk`)
- `.ensi` meshes (in cm) for each resolution (`coarse`, `fine`, `hex`) with all associated field data

- `.ALG` file with the fields in the format used by MonoAlg3D solver

Representative output examples are available in the associated [Zenodo](https://doi.org/10.5281/zenodo.15372895) archive, including a description of the different required and generated fields.

## Troubleshooting and common issues

Common issues encountered by users are typically related to **input surface mesh quality**. While the pipeline performs several automatic checks and corrections, certain cases may require user intervention:

- **Atypical or highly pathological geometries** may lead to incorrect or incomplete anatomical labelling.  
  - In such cases, parameters in the labelling section of the configuration script may need to be adjusted (e.g. when the left ventricular cavity is larger than the right ventricular cavity, the corresponding thresholds should be updated accordingly).

- **Irregular face orientation near the basal plane** can affect automatic detection of basal regions and labelling.  
  - This can often be mitigated by adjusting the relevant labelling parameters in the script.

- **Incorrect label distribution due to inconsistent face normals**, for example when normals do not point outward or are irregularly distributed.  
  - Users are encouraged to visually inspect face normals in the input mesh and correct them if necessary.

- **Non-watertight surface meshes** may prevent successful volumetric mesh generation.  
  - If automatic repair routines fail, the surface mesh must be repaired manually before rerunning the pipeline.

For mesh inspection, visualisation, and manual repair, external tools may be required, such as:
- [ParaView](https://www.paraview.org) 
- [Blender](https://www.blender.org)
- [MeshLab](https://www.meshlab.net)


### License
This project is licensed under the GNU General Public License v3.
This project includes modified code from the gptoolbox (https://github.com/alecjacobson/gptoolbox) and the cobiveco (https://github.com/KIT-IBT/Cobiveco) repositories, licensed under the Apache License, Version 2.0. See the source file(s) for details. In addition, two files in this repository are based on third-party code under the BSD 2-Clause License. See individual file headers for the license text and authorship.