import os
import logging
import sys

from factory_neuron_model import FactoryNeuronCell, Skeleton
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


def compare_swc_without_header(swc_str1: str, swc_str2: str):
    swc_str1_no_header = "\n".join([l for l in swc_str1.split("\n") if not l.startswith("#")])
    swc_str2_no_header = "\n".join([l for l in swc_str2.split("\n") if not l.startswith("#")])
    print(swc_str1_no_header == swc_str2_no_header)
    if swc_str1_no_header != swc_str2_no_header:
        ls1, ls2 = swc_str1_no_header.split("\n"), swc_str2_no_header.split("\n")
        if len(ls1) != len(ls2):
            print(f"Wrong length: {len(ls1)} vs {len(ls2)}")
            return
        for ind, (l1, l2) in enumerate(zip(ls1, ls2)):
            if l1 != l2:
                print(ind, l1, l2)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    # swc converters
    with open(os.path.join("h01_Human_swc_examples", "4559469707.swc"), "r") as f:
        swc_content = f.read()
        skeleton: Skeleton = FactoryNeuronCell.swc_to_skeleton(swc_content=swc_content)
        print(skeleton)
        swc_str: str = FactoryNeuronCell.skeleton_to_swc(skeleton=skeleton)
        compare_swc_without_header(swc_str1=swc_content, swc_str2=swc_str)

    # load neuron as cell (neuron h wrapper)
    use_neurom = True  # if installed (pip install neurom)
    curr_path = os.path.join("L5PC_Mouse_model_example", "morphologies", "cell1.ASC")

    cell1 = NeuronCell(use_cvode=True,
                       model_path=os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])),
                                               "L5PC_Mouse_model_example"))
    # visualize etc?
    # add clamp and play and plot

    if use_neurom:
        use_neurom_visualization(curr_path)

