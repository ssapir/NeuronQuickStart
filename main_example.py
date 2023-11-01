import os
import logging
import sys

import numpy as np
from matplotlib import pyplot as plt

from utils.neuron_model import NeuronCell
from utils.neuron_viewer import plot_morphology_from_cell


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


def simulate(cell: NeuronCell, dt=0.1, delay_ms=200, dur_from_delay_ms=400, initial_voltage=-70):
    def get_recording_vectors():
        rec_time = h.Vector()
        rec_time.record(h._ref_t)
        # record soma voltage
        rec_voltage_soma = h.Vector()
        rec_voltage_soma.record(cell.soma(0.5)._ref_v)
        rec_voltage_all_segments, all_segments = [], []
        for indsec, section in enumerate(cell1.apical_sections + cell.basal_sections):
            for indseg, segment in enumerate(section):
                rec_voltage_seg = h.Vector()
                rec_voltage_seg.record(segment._ref_v)
                rec_voltage_all_segments.append(rec_voltage_seg)
                all_segments.append(segment)  # make sure we know what was recorded (keep order)
        return dict(time=rec_time, soma_v=rec_voltage_soma, recorded_segments=all_segments,
                    segments_v=rec_voltage_all_segments)

    from neuron import h
    from neuron.units import mV, ms

    h.dt = dt
    recorded_data = get_recording_vectors()
    h.finitialize(initial_voltage * mV)
    h.continuerun((delay_ms + dur_from_delay_ms) * ms)  # todo before interp this is 70 points for small dt
    for k, v in recorded_data.items():
        if k == "recorded_segments":
            recorded_data[k] = v
        else:
            if isinstance(v, list):
                recorded_data[k] = [np.array(iv.to_python()) for iv in v]
            else:
                recorded_data[k] = np.array(v.to_python())
    return recorded_data


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


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    # load neuron as cell (neuron h wrapper)
    use_neurom = True  # if installed (pip install neurom)
    curr_path = os.path.join("L5PC_Mouse_model_example", "morphologies", "cell1.swc")  # converted. Can also the asc
    cell1 = NeuronCell(use_cvode=True,
                       morphologyFilename=os.path.join("morphologies", "cell1.ASC"),
                       model_path=os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])),
                                               "L5PC_Mouse_model_example"))
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    plot_morphology_from_cell(ax, cell1, color_by_type=True, fontsize=8)
    ax.axis('equal')
    plt.tight_layout()
    plt.show()

    # add clamp and play and plot
    cell1.init_passive_params()
    cell1.change_passive_params(CM=1.2, RA=200, RM=22000.0, E_PAS=-70, SPINE_START=60, F_factor=1.9)
    delay_ms = 50
    dur_from_delay_ms = 200
    cell1.add_alpha_current_stim(seg=cell1.soma(0.5), delay_ms=delay_ms, dur_from_delay_ms=dur_from_delay_ms)
    recorded_data_traces = simulate(cell=cell1, dt=0.1, initial_voltage=-70,
                                    delay_ms=delay_ms, dur_from_delay_ms=dur_from_delay_ms)
    plt.figure()
    plt.plot(recorded_data_traces["time"], recorded_data_traces["soma_v"], "--k")
    for seg, seg_v in zip(recorded_data_traces["recorded_segments"], recorded_data_traces["segments_v"]):
        plt.plot(recorded_data_traces["time"], recorded_data_traces["soma_v"],
                 label="{0}-{1}".format(seg.sec.name().split(".")[1], seg.x))
    plt.tight_layout()
    # plt.legend()
    plt.show()

    if use_neurom:
        use_neurom_visualization(curr_path)

    # swc converters - will work with proper types / after fixing soma/apical/basal in downloaded swc
    if False:
        with open(os.path.join("h01_Human_swc_examples", "morphologies", "4559469707.swc"), "r") as f:
            swc_content = f.read()
            skeleton: Skeleton = FactoryNeuronCell.swc_to_skeleton(swc_content=swc_content)
            print(skeleton)
            print(skeleton.vertex_types)
            swc_str: str = FactoryNeuronCell.skeleton_to_swc(skeleton=skeleton)
            compare_swc_without_header(swc_str1=swc_content, swc_str2=swc_str)
            with open(os.path.join("h01_Human_swc_examples", "morphologies", "4559469707_rewritten.swc"), "w") as f:
                f.write(swc_str)

        # todo currently fails because the morph has no types (must have soma, apical & basal. it's undefined for h01 swc files)
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
