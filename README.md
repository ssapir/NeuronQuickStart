# NeuronQuickStart
A wrapper to common functionality to work with Neuron package (https://www.neuron.yale.edu/neuron/)

* generic_template.hoc - common/basic hoc file for swc/asc neuron loading 
  (NeuronCell complete additional functionality converted from hoc to python code)
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
* Python files:
  * main_example.py - show common functionality in current repository
  * neuron_model.py - NeuronCell wraps neuron h loading (common functionality and replacement for some hoc code) 
  * factory_neuron_model.py - wrap & refactor converters (swc, asc, NeuronCell)