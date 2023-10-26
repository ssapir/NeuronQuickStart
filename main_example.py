import os
import logging
import sys

from neuron_model import NeuronCell


def use_neurom_visualization(curr_path):
    # import within only if needed
    from neurom import viewer, load_morphology

    nrn = load_morphology(curr_path)
    viewer.draw(nrn)                    # 2d plot
    viewer.draw(nrn, mode='3d')         # 3d plot
    viewer.draw(nrn.neurites[0])        # 2d plot of neurite tree
    viewer.draw(nrn, mode='dendrogram') # dendrogram plot


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    use_neurom = True  # if installed (pip install neurom)
    curr_path = os.path.join("L5PC_Mouse_model_example", "morphologies", "cell1.ASC")

    cell1 = NeuronCell(use_cvode=True,
                       model_path=os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])),
                                                              "L5PC_Mouse_model_example"))
    # visualize etc?
    # add clamp and play and plot

    if use_neurom:
        use_neurom_visualization()
