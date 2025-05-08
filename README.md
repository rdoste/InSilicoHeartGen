# InSilicoHeartGen
Pipeline to automatically create large virtual patient cohort files to conduct large-scale in silico trials through cardiac electromechanical simulations


## Overview
**InSilicoHeartGen** is a pipeline designed to automatically generate large virtual patient cohort files for conducting large-scale **in silico** trials through cardiac electromechanical simulations.

## Installation

### Requirements
- **MATLAB** version **2021b** or later
- **Windows 10** or **Linux** operating system

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

After installing the required dependencies, you can run any of the provided scripts to launch the pipeline.

Before execution, make sure to update the scripts to reflect your specific setup:

- Update file paths and input filenames  
- Set the appropriate mesh format: `cut`, `open`, `UKBB` or `closed`  
- Adjust mesh resolution and labelling options as needed  

Common configuration options include:

- `RVseptal_threshold`: Lower this value to increase the number of detected endocardial septal RV faces  
- `Fiber_info`: Controls the fibre angle generated in the field computation  

**Supported Input Formats**

The pipeline supports common surface mesh formats, including but not limited to:

- `.ply`  
- `.stl`  
- `.vtk`  
- `.vtu`  

**Output**

- VTK raw meshes  
- VTK surface mesh with labels (`labels_final.vtk`
- `.ensi` meshes for each resolution (`coarse`, `fine`, `hex`) with all associated field data
- `.ALG` file with the fields in the format used by MonoAlg3D solver

### License
This project is licensed under the GNU General Public License v3.
This project includes modified code from the gptoolbox (https://github.com/alecjacobson/gptoolbox) and the cobiveco (https://github.com/KIT-IBT/Cobiveco) repositories, licensed under the Apache License, Version 2.0. See the source file(s) for details. In addition, two files in this repository are based on third-party code under the BSD 2-Clause License. See individual file headers for the license text and authorship.