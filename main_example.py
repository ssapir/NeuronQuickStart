import os
import logging
import sys

from neuron_model import NeuronCell


def use_neurom_visualization(curr_path):
    # import within only if needed
    import matplotlib.pyplot as plt
    from neurom import load_morphology
    from neurom.view import matplotlib_impl as mpl_n

    nrn = load_morphology(curr_path)
    plt.close("all")
    mpl_n.plot_morph(nrn)                    # 2d plot
    plt.tight_layout()
    mpl_n.plot_morph3d(nrn)         # 3d plot
    plt.tight_layout()
    mpl_n.plot_tree(nrn.neurites[0])        # 2d plot of neurite tree
    plt.tight_layout()
    mpl_n.plot_dendrogram(nrn)  # dendrogram plot
    plt.tight_layout()
    plt.show()


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
        use_neurom_visualization(curr_path)
