# Based on: https://github.com/seung-lab/cloud-volume/wiki/Advanced-Topic%3A-Skeletons
import datetime
import enum
from collections import defaultdict
import re

import numpy as np


class SwcSectionType(enum.Enum):
    """Note: undefined is not read by morphio & neurom (parsing error)
    From: https://swc-specification.readthedocs.io/en/latest/swc.html
    """
    undefined = 0
    soma = 1
    axon = 2
    basal_dendrite = 3
    apical_dendrite = 4
    point = 5
    end_point = 6
    # custom = 5

    def __str__(self):
        return str(self.value)


class Skeleton:
    IDENTITY = np.array([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
    ], dtype=np.float32)

    def __init__(self, vertices=[], edges=[], radii=[], vertex_types=[], transform_matrix=None):
        """

        :param vertices: [[x,y,z], ...] Nx3 float32
        :param edges: [[v1,v2], ...] Nx2 uint32 (ID of connected vertices)
        :param radii: [r1,r2,...] Nx1 float32 distance from vertex to nearest boudary
        :param vertex_types: uint8 SWC vertex types
        :param transform_matrix:  3x4 scaling and translation matrix (ie homogenous coordinates)
        that represents the transformation from voxel to physical coordinates.
        """
        self.vertices = np.array(vertices, dtype=np.float32)
        self.edges = np.array(edges, dtype=np.uint32)
        self.radii = np.array(radii, dtype=np.float32)
        # patch for an error-  will crush with undefined type
        vertex_types = [v if SwcSectionType(v) != SwcSectionType.undefined else SwcSectionType.apical_dendrite for v in vertex_types]
        self.vertex_types = np.array(vertex_types, dtype=SwcSectionType)
        self.transform_matrix = np.copy(Skeleton.IDENTITY) if transform_matrix is None \
            else np.array(transform_matrix).reshape((3, 4))

    def __str__(self):
        return f"{self.vertices.shape} vertices, {self.edges.shape} edges, {self.radii.shape} radius," + \
               f"{self.vertex_types.shape} vertex_types, transformation: {self.transform_matrix}"


class FactoryNeuronCell:
    """Wraps converters that create NeuronCell object:
    - EM coordinates to swc file / cell instance

    """
    @staticmethod
    def em_coords_to_swc_file():  # todo what are coordinates?
        pass

    @staticmethod
    def swc_to_skeleton(swc_content: str):
        lines = swc_content.split("\n")
        header = [l for l in lines if re.match(r'[#\s]', l)]
        lines = [l for l in lines if l.replace(r"\s", '') != '' and l != '' and l not in header]

        sx, sy, sz = 1, 1, 1
        scale_header = [l.replace("#", "") for l in header if "SCALE" in l]
        if len(scale_header) == 1:
            sx, sy, sz = [float(_) for _ in scale_header[0].replace("SCALE", "").split()]
        transform_matrix = np.array([
            [sx, 0, 0, 0],
            [0, sy, 0, 0],
            [0, 0, sz, 0],
        ], dtype=np.float32)

        if len(lines) == 0:
            return Skeleton(transform_matrix=transform_matrix)

        edges, radii, vertex_types, vertices = \
            FactoryNeuronCell.__generate_skeleton_components_from_swc(lines)

        return Skeleton(vertices, edges, radii, vertex_types, transform_matrix=transform_matrix)


    @staticmethod
    def write_cell_builder_hoc_to_swc(secs: list):
        """Convert sections list (created via Cell Builder for example) to valid swc file to test format
        Note: should contain some (as 1 section!)
        :param secs:
        :return: swc_str
        """
        from neuron import h
        def sec_to_str(c_sec, offset, p_ind, curr_type: SwcSectionType):
            res = ""
            n_points = c_sec.n3d()
            map_parent_to_id[c_sec.name()] = (offset + 1, offset + 1 + n_points - 1)
            for i in range(n_points):
                # index     type         X            Y            Z       radius       parent
                res += "{n} {T} {x:0.6f} {y:0.6f} {z:0.6f} {R:0.6f} {P} # {N}{PN}\n".format(
                    n=(i + 1 + offset), T=int(curr_type.value), N=c_sec.name(),
                    x=h.x3d(i, sec=c_sec), y=h.y3d(i, sec=c_sec), z=h.z3d(i, sec=c_sec), R=h.diam3d(i, sec=c_sec)/2,
                    P=p_ind if i == 0 else (i + offset),
                    PN=f" son of {c_sec.parentseg().sec.name()}" if p_ind != -1 else "")
            return res, offset + n_points

        swc_str = """# swc file format
# Labels: {labels}
""".format(labels=" ,".join([f"{v.name} - {v.value}" for v in SwcSectionType]))
        somas = [s for s in secs if "soma" in s.name()]
        secs_no_soma = [s for s in secs if "soma" not in s.name()]
        if len(somas) != 1:
            print("Error. not single soma, ", somas)
            return
        # add soma
        map_parent_to_id = {"soma": (0, int(h.n3d(somas[0])) - 1)}  # name: (min, max)
        offset = 0
        res, offset = sec_to_str(somas[0], offset=offset, p_ind=-1, curr_type=SwcSectionType.soma)
        swc_str += res
        for sec in secs_no_soma:
            p_ind = map_parent_to_id[sec.parentseg().sec.name()][1]
            res, offset = sec_to_str(sec, offset=offset, p_ind=p_ind,
                                     curr_type=SwcSectionType.apical_dendrite if "apic" in sec.name() else
                                     SwcSectionType.basal_dendrite if "dend" in sec.name() else
                                     SwcSectionType.axon if "axon" in sec.name() else
                                     SwcSectionType.basal_dendrite)
            swc_str += res
        return swc_str

    @staticmethod
    def skeleton_to_swc(skeleton: Skeleton, contributors="", version=""):
        """
        Returns: swc as a string
        """
        swc = FactoryNeuronCell.__swc_header(transform_matrix=skeleton.transform_matrix,
                                             contributors=contributors, version=version) + "\n"
        swc += FactoryNeuronCell.__generate_single_swc_from_skeleton(skeleton, offset=0) + "\n"
        # offset = 0
        # for skel in skeleton.components():  # todo needed? several trees in single coords
        #     swc += FactoryNeuronCell.__generate_single_swc_from_skeleton(skel, offset) + "\n"
        #     offset += skel.vertices.shape[0]
        return swc

    @staticmethod
    def __generate_skeleton_components_from_swc(lines):
        """Refactored from https://github.com/seung-lab/cloud-volume skeleton.py
        :return: edges, radii, vertex_types, vertices #todo add transform from header?
        """
        def str_to_float(val):
            try:
                return float(val)
            except ValueError:
                return -1  # e.g. NA or N/A

        vertices, edges, radii, vertex_types = [], [], [], []
        label_index = {}
        for N, line in enumerate(lines):
            (vid, v_type, x, y, z, radius, parent_id) = line.split(" ")

            vertices.append(tuple([float(_) for _ in (x, y, z)]))
            vertex_types.append(int(v_type))
            radii.append(str_to_float(radius))

            vid = int(vid)
            parent_id = int(parent_id)
            label_index[vid] = N

            if parent_id >= 0:
                edge = [vid, parent_id] if vid < parent_id else [parent_id, vid]
                edges.append(edge)
        for edge in edges:
            edge[0] = label_index[edge[0]]
            edge[1] = label_index[edge[1]]
        return edges, radii, vertex_types, vertices

    @staticmethod
    def __generate_single_swc_from_skeleton(skeleton: Skeleton, offset: int):
        """Refactored from https://github.com/seung-lab/cloud-volume skeleton.py
        :return: swc string
        """
        swc = ""
        if skeleton.edges.size == 0:
            return ""

        index = defaultdict(set)
        visited = defaultdict(bool)
        for e1, e2 in skeleton.edges:
            index[e1].add(e2)
            index[e2].add(e1)

        stack, parents = [skeleton.edges[0, 0]], [-1]
        while stack:
            node = stack.pop()
            parent = parents.pop()

            if visited[node]:
                continue

            swc += "{n} {T} {x:0.6f} {y:0.6f} {z:0.6f} {R:0.6f} {P}\n".format(
                n=(node + 1 + offset),
                T=skeleton.vertex_types[node],
                x=skeleton.vertices[node][0],
                y=skeleton.vertices[node][1],
                z=skeleton.vertices[node][2],
                R=skeleton.radii[node],
                P=parent if parent == -1 else (parent + 1 + offset))
            visited[node] = True

            for child in index[node]:
                stack.append(child)
                parents.append(node)

        return swc

    @staticmethod
    def __swc_header(transform_matrix, contributors="", version=""):
        """Refactored from Taken from https://github.com/seung-lab/cloud-volume skeleton.py

        :param transform_matrix:
        :param contributors:
        :param version:
        :return:
        """
        sx, sy, sz = np.diag(transform_matrix)[:3]
        return f"""#
# CREATURE 
# REGION
# FIELD/LAYER
# TYPE
# CONTRIBUTOR {contributors}
# REFERENCE
# RAW 
# EXTRAS 
# SOMA_AREA
# SHINKAGE_CORRECTION 
# VERSION_NUMBER {version}
# VERSION_DATE {datetime.datetime.utcnow().isoformat()}
# SCALE {sx:.6f} {sy:.6f} {sz:.6f}
"""
