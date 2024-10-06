# NeuronQuickStart
A wrapper to common functionality to work with Neuron package (https://www.neuron.yale.edu/neuron/)

## Overview
### Installation:
* git clone <current repository path>
* conda create --name  neuron  python=3.8
* conda activate neuron
* pip install -r requirements.txt
* pip install neuron==8.0.2  # in linux. maybe different in windows

* For windows, may requrire to manually install neuron separately: https://www.neuron.yale.edu/neuron/download 
* Can updgrade for newer python & neuron version if needed (not validated against the code here)

### Compile mod (mechanism) files:
* cd L5PC_Mouse_model_example
* nrnivmodl mods  # Important: after conda activate neuron!

### Usage example: 
python main_example.py # or via IDE run the python code

* main_example.py / ipynb for a short overview of functionality
* Neuron's python tutorials: https://nrn.readthedocs.io/en/8.2.6/tutorials/index.html
  
## Repository structure:
* Python files:
  * main_example.py & main_example.ipynb- show functionality in current repository
  * utils/neuron_model.py - NeuronCell wraps neuron h loading (common functionality and replacement for some hoc code) 
  * utils/factory_neuron_model.py - wrap & refactor converters (swc, asc, NeuronCell)
  * utils/neuron_viewer.py - plot cell (with extensions of markers, electrodes and scalebars)
* Neuron models:
  * h01_Human_swc_example - contains downloaded examples from h01 data
  * L5PC_Mouse_model_example - commonly used neuron model (Etay Hay, 2011). 
  Contains:
    * mods - folder with needed mechanisms (mod) files. 
      Compile with: "nrnivmodl mods" (mods is the mechanisms' folder name). 
      Will create output file within the folder you run the script from.
    * morphologies - *.ASC or *.SWC ascii (textual) description files with morphology description.
      ASC - tree-like structure.
      SWC - table-like structure with columns for: (id, type, x, y , z, diam, parent-id).
      Can be viewed with the following tool: https://neuroinformatics.nl/HBP/morphology-viewer/
    * hoc files - templates used when loading the model (can be a single generic template instead)
       * generic_template.hoc - common/basic hoc file for swc/asc neuron loading
       * hoc files can be re-written to use the python API (as in this repository)
    (NeuronCell complete additional functionality converted from hoc to python code)
