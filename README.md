# NeuronQuickStart
A wrapper to common functionality to work with Neuron package (https://www.neuron.yale.edu/neuron/)

Installation:
* pip install -r requirements.txt
* Manually install neuron separately: https://www.neuron.yale.edu/neuron/download 

Usage example: look (& run) main_example.py for a short overview of functionality

Repository structure:
* Python files:
  * main_example.py - show functionality in current repository
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
    * morphologies - *.ASC or *.SWC files with morphology description.
      ASC - ascii (textual) description.
      SWC - binary description.
      Can be viewed with the following tool: https://neuroinformatics.nl/HBP/morphology-viewer/
    * hoc files - templates used when loading the model (can be a single generic template instead)
       * generic_template.hoc - common/basic hoc file for swc/asc neuron loading 
    (NeuronCell complete additional functionality converted from hoc to python code)
