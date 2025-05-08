# InSilicoHeartGen
Pipeline to automatically create large virtual patient cohort files to conduct large-scale in silico trials through cardiac electromechanical simulations


## Overview
**InSilicoHeartGen** is a pipeline designed to automatically generate large virtual patient cohort files for conducting large-scale **in silico** trials through cardiac electromechanical simulations.

## Installation

### Requirements
- **MATLAB** version **2021b** or later
- **NVIDIA GPU**
- **Windows 10** or **Linux** operating system

### Dependencies
To ensure full functionality, the following dependencies must be installed. Copy the required libraries into the `dependencies` folder:

#### MATLAB Toolboxes:
- [Iso2Mesh](https://github.com/fangq/iso2mesh)
- [toolbox_graph](https://github.com/gpeyre/matlab-toolboxes/tree/master/toolbox_graph)
- [gptoolbox](https://github.com/alecjacobson/gptoolbox)
- [Cobiveco](https://github.com/KIT-IBT/Cobiveco)
- [vtkToolbox](https://github.com/KIT-IBT/vtkToolbox) *(Requires additional installation steps; follow the instructions in their README.)*

### Usage
Once the dependencies are installed, execute any provided script within the repository to initiate the pipeline.

### License
This project is licensed under the GNU General Public License v3.
This project includes modified code from the gptoolbox (https://github.com/alecjacobson/gptoolbox) and the cobiveco (https://github.com/KIT-IBT/Cobiveco) repositories, licensed under the Apache License, Version 2.0. See the source file(s) for details. In addition, one or more files in this repository are based on third-party code under the BSD 2-Clause License. See individual file headers for the license text and authorship.