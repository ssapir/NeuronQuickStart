import os
import logging
import sys

from matplotlib import pyplot as plt

from factory_neuron_model import FactoryNeuronCell, Skeleton
from neuron_model import NeuronCell
from neuron_viewer import plot_morphology_from_cell


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


def compare_swc_without_header(swc_str1: str, swc_str2: str, ignore_ids=True):
    """Compare 2 strings (test swc converters, back and forth)

    :param ignore_ids: id can change when creating the skeleton... todo ignore here or
    :return:
    """
    swc_str1_no_header = "\n".join([l for l in swc_str1.split("\n") if not l.startswith("#") and l != ""])
    swc_str2_no_header = "\n".join([l for l in swc_str2.split("\n") if not l.startswith("#") and l != ""])
    print(swc_str1_no_header == swc_str2_no_header)
    if swc_str1_no_header != swc_str2_no_header:
        ls1, ls2 = swc_str1_no_header.split("\n"), swc_str2_no_header.split("\n")
        print(len(ls1), len(ls2))
        if len(ls1) != len(ls2):
            print(f"Wrong length: {len(ls1)} vs {len(ls2)}")
        for ind, (l1, l2) in enumerate(zip(ls1, ls2)):
            if l1 != l2 and ind < 100:
                pass
                # print(l1, " || ", l2)
            else:
                pass
                # print(ind, l1)


def load():  # todo currently fails because the morph has no types (must have soma, apical & basal. it's undefined for h01 swc files)
    cell1 = NeuronCell(use_cvode=True, morphologyFilename=os.path.join("morphologies", "4559469707.swc"),
                       biophysicalModelTemplateFilename="generic_template.hoc", templateName="Celltemplate",
                       biophysicalModelFilename=None,  # inside template already
                       model_path=os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])),
                                               "h01_Human_swc_examples"))
    print(cell1.soma, cell1.apical_sections)
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    plot_morphology_from_cell(ax, cell1, color_by_type=True, fontsize=8)
    ax.axis('equal')
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    # swc converters
    if False:
        with open(os.path.join("h01_Human_swc_examples", "morphologies", "4559469707.swc"), "r") as f:
            swc_content = f.read()
            skeleton: Skeleton = FactoryNeuronCell.swc_to_skeleton(swc_content=swc_content)
            print(skeleton)
            print(skeleton.vertex_types)
            swc_str: str = FactoryNeuronCell.skeleton_to_swc(skeleton=skeleton)
            compare_swc_without_header(swc_str1=swc_content, swc_str2=swc_str)
            # with open(os.path.join("h01_Human_swc_examples", "morphologies", "4559469707_rewritten.swc"), "w") as f:
            #     f.write(swc_str)

    # load neuron as cell (neuron h wrapper)
    use_neurom = True  # if installed (pip install neurom)
    curr_path = os.path.join("L5PC_Mouse_model_example", "morphologies", "cell1.swc")  # converted. Can also the asc
    cell1 = NeuronCell(use_cvode=True,
                       morphologyFilename=os.path.join("morphologies", "cell1.swc"),
                       model_path=os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])),
                                               "L5PC_Mouse_model_example"))
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    plot_morphology_from_cell(ax, cell1, color_by_type=True, fontsize=8)
    ax.axis('equal')
    plt.tight_layout()
    plt.show()

    # add clamp and play and plot

    if use_neurom:
        use_neurom_visualization(curr_path)
