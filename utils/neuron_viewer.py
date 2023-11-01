import math

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.linalg import norm
import matplotlib
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, Circle, Polygon
from matplotlib_scalebar.scalebar import ScaleBar

from utils.neuron_model import NeuronCell


def generate_cylindrical_points(start, end, start_radius, end_radius,
                                linspace_count=300):
    """Generate a 3d mesh of a cylinder with start and end points, and varying radius.
    Based on: http://stackoverflow.com/a/32383775
    """

    def _get_normals(v):
        """Get two vectors that form a basis w/ v.
        Note: returned vectors are unit
        """
        not_v = np.array([1, 0, 0])
        if np.all(np.abs(v) == not_v):
            not_v = np.array([0, 1, 0])
        n1 = np.cross(v, not_v)
        n1 /= norm(n1)
        n2 = np.cross(v, n1)
        return n1, n2

    v = end - start
    length = norm(v)
    v = v / length
    n1, n2 = _get_normals(v)

    # pylint: disable=unbalanced-tuple-unpacking
    l, theta = np.meshgrid(np.linspace(0, length, linspace_count),
                           np.linspace(0, 2 * np.pi, linspace_count))

    radii = np.linspace(start_radius, end_radius, linspace_count)
    rsin = np.multiply(radii, np.sin(theta))
    rcos = np.multiply(radii, np.cos(theta))

    return np.array([start[i] +
                     v[i] * l +
                     n1[i] * rsin + n2[i] * rcos
                     for i in range(3)])


def get_cylinder_hull(plane, start, end, start_radius, end_radius):
    from scipy.spatial import ConvexHull
    points = generate_cylindrical_points(start, end, start_radius, end_radius, 100)
    pts = np.array([points[0].ravel(), points[1].ravel()]).T
    return ConvexHull(pts), pts


def map_cell_to_xyzd(cell: NeuronCell, times_dict={}):
    from neuron import h

    all_segment_coords, all_section_coords = {}, {}
    soma = {"x": np.mean([h.x3d(i, sec=cell.soma) for i in range(int(h.n3d(sec=cell.soma)))]),
            "y": np.mean([h.y3d(i, sec=cell.soma) for i in range(int(h.n3d(sec=cell.soma)))]),
            "z": np.mean([h.z3d(i, sec=cell.soma) for i in range(int(h.n3d(sec=cell.soma)))]),
            "d": np.mean([h.diam3d(i, sec=cell.soma) for i in range(int(h.n3d(sec=cell.soma)))])}
    # for what, sections_list in zip(["apical", "basal", "soma", "axon"],
    #                               [cell.ApicalSectionsList, cell.BasalSectionsList,
    #                                cell.SomaSectionsList, cell.AxonalSectionsList]):
    for what, sections_list in zip(["all"], [cell.all]):
        for sec_ind, sec in enumerate(sections_list):
            what = "basal" if "dend" in str(sec) else "all"
            what = "apical" if "apic" in str(sec) else what
            what = "soma" if "soma" in str(sec) else what
            what = "axon" if "axon" in str(sec) else what
            x_path = [h.x3d(i, sec=sec) for i in range(int(h.n3d(sec=sec)))]
            y_path = [h.y3d(i, sec=sec) for i in range(int(h.n3d(sec=sec)))]
            z_path = [h.z3d(i, sec=sec) for i in range(int(h.n3d(sec=sec)))]
            d_path = [h.diam3d(i, sec=sec) for i in range(int(h.n3d(sec=sec)))]
            sec_name = sec.name().split('.')[-1].replace('[', '_')[:-1]

            sec_type = "trunk" if sec in cell.trunk_sections else "oblique" if sec in cell.oblique_sections else what
            ks = [k for k in times_dict.keys() if sec.name().split(".")[-1] in k]
            color_ = None
            if len(times_dict.keys()) > 0:
                color_ = times_dict[ks[0]] if len(ks) == 1 else np.nan
            all_section_coords[(what, sec_ind, 'all')] = {'sec name': sec_name, 'seg index': 0, 'what': sec_type,
                                                          "color": color_,
                                                          'x': x_path, 'y': y_path, 'z': z_path, 'd': d_path}
        # for seg_ind in range(sec.nseg):  # this needs to be calculated
        #     all_segment_coords[(sec_ind, seg_ind)] = {}
        #     all_segment_coords[(sec_ind, seg_ind)]['sec name'] = sec_name
        #     all_segment_coords[(sec_ind, seg_ind)]['seg index'] = seg_ind
        #     all_segment_coords[(sec_ind, seg_ind)]['x'] = x_path  # [x_path[k] for k in curr_seg_inds]
        #     all_segment_coords[(sec_ind, seg_ind)]['y'] = y_path  # [y_path[k] for k in curr_seg_inds]
        #     all_segment_coords[(sec_ind, seg_ind)]['z'] = z_path  # [z_path[k] for k in curr_seg_inds]
        #     all_segment_coords[(sec_ind, seg_ind)]['d'] = d_path  #[d_path[k] for k in curr_seg_inds]
    return soma, all_section_coords, all_segment_coords


def plot_electrode(ax, x, y, color="gray", inner_color="silver", bigger=False,
                   shift_bottom=False, red_x=0, red_y=0, is_filled=True):
    sft_x, sft_y = 0, 0
    if bigger:
        if shift_bottom:
            sft_x, sft_y = 15, -30
        up_add_x, up_add_y = 60, 180
        ellipse_x, ellipse_y = 70, 170
        down_add_x, down_add_y = 80, 160
        ellipse_d_x, ellipse_d_y = 25, 7
    else:
        if shift_bottom:
            sft_x, sft_y = 15, -5
        up_add_x, up_add_y = 42, 90
        ellipse_x, ellipse_y = 46, 85
        down_add_x, down_add_y = 50, 80
        ellipse_d_x, ellipse_d_y = 10, 7
    ax.plot([x, x + sft_x + up_add_x], [y, y + sft_y + up_add_y], color=color, linewidth=.4)
    ax.plot([x, x + sft_x + down_add_x + red_x], [y, y + sft_y + down_add_y + red_y], color=color, linewidth=.4)
    ax.add_patch(
        Ellipse([x + sft_x + ellipse_x, y + sft_y + ellipse_y], ellipse_d_x, ellipse_d_y, angle=-45, fill=False,
                color=color, linewidth=.4, zorder=2))
    if is_filled:
        ax.fill([x, x + sft_x + up_add_x, x + sft_x + down_add_x, x],
                [y, y + sft_y + up_add_y, y + sft_y + down_add_y, y],
                color=inner_color, alpha=alpha)
        ax.add_patch(Ellipse([x + sft_x + ellipse_x, y + sft_y + ellipse_y], ellipse_d_x, ellipse_d_y, angle=-45,
                             fill=True, alpha=alpha + 0.3, color=color, linewidth=0, zorder=2))
        print("sapsap Ellipse ", color, [x + sft_x + ellipse_x, y + sft_y + ellipse_y])


def plot_morphology_from_cell(ax, cell: NeuronCell, segment_colors=None, width_mult_factors=None, colors_dict={},
                              fontsize=3, plot_per_segment=False, color_by_type=False, soma_as_cylinder=False,
                              with_legend=False, with_text=False, is_scalebar=False, text_dict={}, times_dict={},
                              only_soma_marker=False, alpha=0.2,
                              with_markers=True, only_apical=False, only_basal=False, with_soma=True, with_lims=True,
                              shift_x=0, shift_y=0, seg_colors_cmap=plt.cm.jet, is_electrode=False, is_marker=True,
                              seg_path=None, seg_path_color="lime", seg_electrode=None, soma_elect_color="purple",
                              width_fraction=0.01 / 5, fixed_value=100, no_scalebar_no_ax=False):
    soma_mean, all_section_coords, all_segment_coords = map_cell_to_xyzd(cell, times_dict=times_dict)
    data_dict = all_segment_coords if plot_per_segment else all_section_coords

    if segment_colors is None:
        segment_colors = np.arange(len(data_dict.keys()))
        segment_colors = segment_colors / segment_colors.max()
    if width_mult_factors is None:
        width_mult_factors = 1.2 * np.ones(segment_colors.shape)
    colors = seg_colors_cmap(segment_colors)
    colors[np.isnan(segment_colors)] = (0, 0, 0, 1)
    c_to_t = {"axon": "gray", "apical": plt.cm.Blues, "basal": plt.cm.Reds}  # default jet
    map_name_c_to_t = {"axon": "gray", "apical": "blue", "basal": "red"}
    existing_colors = []

    if color_by_type:
        curr_types = [curr_data['what'] for curr_data in data_dict.values()]
        for curr in np.unique(curr_types):
            locs = np.where(curr == np.array(curr_types))[0]
            if len(colors_dict.keys()) > 0:
                cmap = ListedColormap([colors_dict.get(curr)] * len(locs))
                existing_colors.append(colors_dict.get(curr))
            else:
                cmap = c_to_t[curr] if curr in c_to_t.keys() else plt.cm.jet
                get_c = lambda v, n: v if isinstance(cmap, str) else map_name_c_to_t[n]
                existing_colors.append(get_c(c_to_t[curr], n=curr) if curr in c_to_t.keys() else None)
            from_ind = 1 if len(locs) == 1 else len(locs) // 3
            segment_colors = np.arange(from_ind, from_ind + len(locs))
            colors[locs] = cmap(segment_colors / segment_colors.max()) if not isinstance(cmap, str) else plt.cm.gray(
                0.5)

    data_x, data_y = [], []
    colors_added = []
    to_vec = lambda curr_data, k: curr_data[k] - soma_mean[k]  # normalized to 0 of soma
    get_loc = lambda curr_data, k: to_vec(curr_data, k)[
        np.max(np.argmax(np.abs([to_vec(curr_data, 'x'), to_vec(curr_data, 'y')]), axis=1))]  # max in 2d values
    for ind, key in enumerate(data_dict.keys()):
        curr_data = data_dict[key]
        if len(curr_data['x']) == 0 or len(curr_data['y']) == 0 or len(curr_data['d']) == 0:
            print("Error. Lacking {0} {1}".format(ind, key))
            continue
        seg_line_width = width_mult_factors[ind] * np.array(curr_data['d']).mean()
        if curr_data['what'] == "soma":
            continue
        if not only_basal and only_apical and curr_data['what'] not in ["apical", "oblique", "trunk"]:
            continue
        if not only_apical and only_basal and curr_data['what'] != "basal":
            continue
        part_name = ((" " + curr_data['what'][:1] if curr_data['what'] != "axon" else " axon") if not color_by_type else "")
        add_me = "" if not any([_ in curr_data['sec name'].split("_")[1] for _ in ["36", "70"]]) else "(1)"  # (1) is the end/max
        if with_text and norm([get_loc(curr_data, 'x'), get_loc(curr_data, 'y')]) > 40:
            if len(text_dict.keys()) == 0 or text_dict.get(curr_data['what'], None) is not None:
                ax.text(get_loc(curr_data, 'x') - shift_x, get_loc(curr_data, 'y') - shift_y,
                        curr_data['sec name'].split("_")[1] + add_me + part_name, size=fontsize)
        cl = colors[ind] if curr_data["color"] is None else \
            seg_colors_cmap(curr_data["color"]) if not np.isnan(curr_data["color"]) else "k"  # colors[ind]
        order_ = 1
        if seg_path is not None:
            if curr_data['sec name'] in seg_path:
                cl = seg_path_color
                order_ = 2

        # if "k" != cl and curr_data["color"] is not None:
        #     colors_added.append(curr_data["color"])
        #     print("------", curr_data['what'], cl, curr_data["color"])
        ax.plot(to_vec(curr_data, 'x'), to_vec(curr_data, 'y'), lw=seg_line_width, color=cl, zorder=order_)
        data_x.extend(to_vec(curr_data, 'x'))
        data_y.extend(to_vec(curr_data, 'y'))

    if with_lims:
        ax.set_ylim([min(data_y), max(data_y)])
        ax.set_xlim([min(data_x), max(data_x)])
    if with_soma:
        for ind, key in enumerate(data_dict.keys()):
            curr_data = data_dict[key]
            seg_line_width = width_mult_factors[ind] * np.array(curr_data['d']).mean()
            if curr_data['what'] == "soma":
                points = np.array([curr_data['x'] - soma_mean['x'], curr_data['y'] - soma_mean['y'],
                                   curr_data['z'] - soma_mean['z']]).T
                if soma_as_cylinder:
                    hull, pts = get_cylinder_hull(points, points[0, :], points[-1, :],
                                                  curr_data['d'][0] / 2, curr_data['d'][-1] / 2)
                    ax.add_patch(Polygon(pts[hull.vertices, :], fill=True, color='k', linewidth=1, zorder=2))
                else:
                    dist = lambda p0, p1: np.dot(np.subtract(p0, p1), np.subtract(p0, p1))
                    area = lambda p0, p1, r0, r1: math.pi * (r0 + r1) * math.sqrt((r0 - r1) ** 2 + dist(p0, p1))
                    soma_areas = [area(p0, p1, d0 / 2, d1 / 2) for (p0, p1, d0, d1)
                                  in zip(points, points[1:], curr_data['d'], curr_data['d'][1:])]
                    soma_radius = math.sqrt(sum(soma_areas) / (4. * math.pi))
                    center = np.mean(points, axis=0)
                    ax.add_patch(Circle(center, soma_radius, fill=True, color='k', linewidth=1, zorder=2))
                    print("Soma ", list(center), soma_radius)

    if len(cell.markers) > 0 and with_markers and (is_electrode or is_marker):
        circle_patch = [v for v in cell.markers if "circle" in v['label'].lower()]
        arrow_patch = [v for v in cell.markers if "arrow" in v['label'].lower()]
        # for c in circle_patch:
        #     curr_data = {'x': c['points'][0][0], 'y': c['points'][0][1]}
        #     x, y = to_vec(curr_data, 'x'), to_vec(curr_data, 'y')
        #     ax.add_patch(Circle([x, y], 1, fill=True, color='m', linewidth=1, zorder=2))
        matplotlib.style.use("seaborn")
        for a, c in zip(circle_patch, arrow_patch):
            curr_data = {'x': c['points'][0][0], 'y': c['points'][0][1]}
            ar_curr_data = {'x': a['points'][0][0], 'y': a['points'][0][1]}
            x, y = to_vec(curr_data, 'x'), to_vec(curr_data, 'y')
            ar_x, ar_y = to_vec(ar_curr_data, 'x'), to_vec(ar_curr_data, 'y')
            # ax.add_patch(Arrow(x, y, abs(x - ar_x), abs(y - ar_y), fill=True, color='k', linewidth=2, zorder=2))
            if not is_electrode:
                ax.add_patch(Circle([x, y], 5, fill=True, color='m', linewidth=3, zorder=2))
                ax.add_patch(Circle([ar_x, ar_y], 5, fill=True, color='orange', linewidth=1, zorder=2))
            else:
                if not only_soma_marker:
                    plot_electrode(x, y, ax=ax, color="c", inner_color="c")
                if seg_electrode is not None:
                    from neuron import h
                    x_path = [h.x3d(i, sec=seg_electrode.sec) for i in range(int(h.n3d(sec=seg_electrode.sec)))]
                    y_path = [h.y3d(i, sec=seg_electrode.sec) for i in range(int(h.n3d(sec=seg_electrode.sec)))]
                    plot_electrode(x_path[-1], y_path[-1], ax=ax, color="c", inner_color="c")
        if is_electrode:
            plot_electrode(0, 0, ax=ax, color=soma_elect_color, inner_color=soma_elect_color, shift_bottom=True,
                           red_y=-0, red_x=0)  # 60, -90 | 80, -60 | 100, -30

    if with_legend and color_by_type:
        legend_elements = []
        for curr_color in list(set(existing_colors)):
            if curr_color in colors_dict.values():
                key = [k for k, v in colors_dict.items() if v == curr_color][0]
                legend_elements.append(Line2D([0], [0], color=colors_dict[key], lw=2, label=key))
            elif curr_color in map_name_c_to_t.values():
                key = [k for k, v in map_name_c_to_t.items() if v == curr_color][0]
                legend_elements.append(Line2D([0], [0], color=map_name_c_to_t[key], lw=2, label=key))
        ax.legend(handles=legend_elements, loc="upper left")

    if not is_scalebar and not no_scalebar_no_ax:
        ax.set_xlabel("um")
        ax.set_ylabel("um")
    else:
        ax.axis('off')
        if not no_scalebar_no_ax:
            ax.add_artist(
                ScaleBar(1, "um", fixed_value=fixed_value, location="lower left", width_fraction=width_fraction,
                         pad=0, frameon=False, border_pad=0, sep=5, ))
