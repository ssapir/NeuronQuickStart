import numpy as np
from neuron import h
import os
from morphio import Morphology
import logging


class NeuronCell:
    """Wraps the basic neuron loading and model creation
    """
    @property
    def soma(self):
        soma_sections = [self.L5PC.soma[x] for x in range(len(self.L5PC.soma))]
        assert(len(soma_sections) == 1)
        return soma_sections[0]

    @property
    def all(self):
        return self.basal_sections + self.apical_sections + [self.soma] + self.axonal_sections

    @property
    def basal_sections(self):
        return [self.L5PC.dend[x] for x in range(len(self.L5PC.dend))
                if self.L5PC.dend[x].name() not in self.deleted_secs]

    @property
    def apical_sections(self):
        return [self.L5PC.apic[x] for x in range(len(self.L5PC.apic))
                if self.L5PC.apic[x].name() not in self.deleted_secs]

    @property
    def axonal_sections(self):
        return [self.L5PC.axon[x] for x in range(len(self.L5PC.axon))
                if self.L5PC.axon[x].name() not in self.deleted_secs]

    @property
    def oblique_sections(self):
        if hasattr(self, "oblique_list"):
            return self.oblique_list
        return []

    @property
    def trunk_sections(self):
        if hasattr(self, "trunk_list"):
            return self.trunk_list
        return []

    @property
    def tuft_sections(self):
        if hasattr(self, "tuft_list"):
            return self.tuft_list
        return []

    def __init__(self,
                 is_init_trunk_oblique=False,
                 is_delete_axon=False,
                 use_cvode=True, model_path="./L5PC_Mouse_model_example",
                 templateName='L5PCtemplate',
                 name="L5PC_Mouse_model_example",
                 morphologyFilename=os.path.join("morphologies", "cell1.ASC"),
                 biophysicalModelFilename="L5PCbiophys5b.hoc",
                 biophysicalModelTemplateFilename="L5PCtemplate_2.hoc"):

        h.load_file('nrngui.hoc')
        h.load_file("import3d.hoc")
        h.load_file("stdrun.hoc")

        self.deleted_secs = []
        self.default_v = -81.0

        self.L5PC = None  # model from template

        # Current clamp, vectors & parameters
        self.icl = None
        self.alpha_current_vec = None
        self.alpha_time_vec = None
        self.alpha_params = None

        self.__load_mechanisms(model_path)

        morphologyFilename = self.__load_morphology(biophysicalModelFilename,
                                                    biophysicalModelTemplateFilename, model_path,
                                                    morphologyFilename, templateName)
        self.name = name

        if is_delete_axon:  # can cause internal neuron failure in some cases without deletion
            self.L5PC.delete_axon()
            logging.info("Done delete axon")

        self.markers = self.__get_markers(os.path.join(model_path, morphologyFilename))
        logging.info(f"Done get_markers: {np.unique([c.keys() for c in self.markers])}")

        if is_init_trunk_oblique:
            # default. Should match the cell morph from trial and error
            diam_diff_threshold_um = 0.5
            stop_trunk_at = None
            logging.info("Set trunk oblique")
            self.__set_trunk_oblique_from_morph(diam_diff_threshold_um=diam_diff_threshold_um, stop_trunk_at=stop_trunk_at)

        if use_cvode:
            h.CVode().active(1)

        # todo add passive/active parameters?

        print("Done __init__")

    def nearest_segment_to_marker_coords(self, marker_dict):
        closest_sec, sec_x, dist = None, None, np.inf
        x, y, z = marker_dict["points"][0]
        for sec_ind, sec in enumerate(self.all):
            x_path = [h.x3d(i, sec=sec) for i in range(int(h.n3d(sec=sec)))]
            y_path = [h.y3d(i, sec=sec) for i in range(int(h.n3d(sec=sec)))]
            z_path = [h.z3d(i, sec=sec) for i in range(int(h.n3d(sec=sec)))]
            distances = np.sqrt((np.array(x_path) - x) ** 2 + (np.array(y_path) - y) ** 2 + (np.array(z_path) - z) ** 2)
            if distances.min() < dist:
                dist = distances.min()
                closest_sec = sec
                sec_x = distances.argmin() / len(x_path)
        return closest_sec(sec_x), dist

    def add_current_stim(self, seg, delay_ms=200, dur_from_delay_ms=400, amp_ns=0.1):
        """

        :param seg: segment to stimulation location
        :param delay_ms: let model stabilize before simulating (at least 100ms)
        :param dur_from_delay_ms:
        :param amp_ns:
        :return:
        """
        if self.icl is not None:
            del self.icl
        self.icl = h.IClamp(seg.x, sec=seg.sec)
        self.icl.delay = delay_ms  # ms - let model stabilize before changing
        self.icl.dur = dur_from_delay_ms  # ms
        self.icl.amp = amp_ns  # nS

    @staticmethod
    def __alpha_stim(delay_ms=200, dur_from_delay_ms=400, amp_ns=0.1, tau0=10, tau1=15, dt=0.1):
        """ Create time & current vectors of alpha shape with following parameters

        :param delay_ms: start shape from certain delay
        :param dur_from_delay_ms:
        :param amp_ns:
        :param tau0:
        :param tau1:
        :param dt:
        :return:
        """
        time = np.arange(0, delay_ms + dur_from_delay_ms, dt)
        time_for_exp = np.arange(0, dur_from_delay_ms, dt)
        current_vec = np.zeros(time.shape)
        from_t = int(np.round(delay_ms/dt))
        current_vec[from_t:] = amp_ns * ((1 - np.exp(-time_for_exp / tau0)) - (1 - np.exp(-time_for_exp / tau1)))
        return time, current_vec

    def add_alpha_current_stim(self, seg, delay_ms=200, dur_from_delay_ms=400, amp_ns=0.1, tau0=5, tau1=8):
        if self.icl is not None:
            del self.icl

        time, current = self.__alpha_stim(delay_ms=delay_ms, dur_from_delay_ms=dur_from_delay_ms, amp_ns=amp_ns,
                                          tau0=tau0, tau1=tau1, dt=h.dt)
        self.icl = h.IClamp(seg.x, sec=seg.sec)
        self.icl.dur = 1e4  # ms
        self.alpha_current_vec = h.Vector(current)
        self.alpha_time_vec = h.Vector(time)
        self.alpha_current_vec.play(self.icl._ref_amp, self.alpha_time_vec)
        self.alpha_params = dict(delay_ms=delay_ms, dur_from_delay_ms=dur_from_delay_ms, amp_ns=amp_ns)

    def init_passive_params(self):
        for sec in self.all:
            sec.insert('pas')
            for seg in sec:
                seg.pas.e = 0

    @staticmethod
    def sec_to_type_name(sec, map_name=None):
        if map_name is None:
            map_name = {'apic': 'apical', 'dend': 'basal', 'axon': 'axonal', 'soma': 'soma'}
        curr_name = sec.name().split(".")[1].split("[")[0].replace("]", "")
        return map_name.get(curr_name, curr_name)

    @staticmethod
    def calculate_impedance_Rin_Rtr(x_sec, input_electrode_sec, freq_hz=0):
        imp = h.Impedance()
        imp.loc(input_electrode_sec.x, sec=input_electrode_sec.sec)  # location of input electrode (middle of section).
        imp.compute(freq_hz)
        return {"Rin": imp.input(x_sec.x, sec=x_sec.sec), "Rtr_M_ohm": imp.transfer(x_sec.x, sec=x_sec.sec)}

    def change_passive_params(self, CM=1, RA=250, RM=20000.0, E_PAS=-70, SPINE_START=60, F_factor=1.9):
        h.distance(0, 0.5, sec=self.soma)

        for sec in self.all:
            sec.Ra = RA  # Ohm-cm
            sec.cm = CM  # uF/cm^2
            sec.g_pas = 1.0 / RM  # RM: Ohm-cm^2
            sec.e_pas = E_PAS

        for sec in self.apical_sections + self.basal_sections:
            for seg in sec:
                if self._get_distance_between_segments(origin_segment=self.soma(0.5), to_segment=seg) > SPINE_START:
                    seg.cm *= F_factor
                    seg.g_pas *= F_factor

    @staticmethod
    def _get_distance_between_segments(origin_segment, to_segment):
        """ (sec contains segments, equivalent to RC-circuits)
        Example: soma(0.5) origin, and a segment to connect to
        :param origin_segment:
        :param to_segment:
        :return: distance in um
        """
        h.distance(0, origin_segment.x, sec=origin_segment.sec)
        return h.distance(to_segment.x, sec=to_segment.sec)

    def __load_morphology(self, biophysicalModelFilename, biophysicalModelTemplateFilename, model_path,
                          morphologyFilename, templateName):
        if biophysicalModelFilename is not None:
            if os.path.isfile(os.path.join(model_path, biophysicalModelFilename)):
                h.load_file(os.path.join(model_path + "/", biophysicalModelFilename))
            else:
                logging.error("Missing biophysics file {0}".format(os.path.join(model_path, biophysicalModelFilename)))
        if not os.path.isfile(os.path.join(model_path, biophysicalModelTemplateFilename)):
            raise Exception(
                "Missing model template {0}".format(os.path.join(model_path, biophysicalModelTemplateFilename)))
        if hasattr(h, templateName):
            return
            # delattr(h, templateName)
        if morphologyFilename.startswith("/"):
            morphologyFilename = morphologyFilename[1:]

        if morphologyFilename.endswith(".swc"):  # patch. Can load swc but for some reason it fails now
            # Import2 = h.Import3d_SWC_read()
            # # Import2.input(morphologyFilename)  # should work but template fails so just convert
            from morph_tool import convert
            convert(os.path.join(model_path, morphologyFilename),
                    os.path.join(model_path, morphologyFilename.replace(".swc", "_converted.ASC")))
            morphologyFilename = morphologyFilename.replace(".swc", "_converted.ASC")
        h.load_file(os.path.join(model_path + "/", biophysicalModelTemplateFilename))

        logging.info("Loading {0}".format(os.path.join(model_path, morphologyFilename)))
        cell_temp_func = getattr(h, templateName)  # todo can be automatic in the template as cell?
        self.L5PC = cell_temp_func(os.path.join(model_path, morphologyFilename))
        logging.info("Loaded {0}".format(morphologyFilename))
        return morphologyFilename

    def __load_mechanisms(self, model_path):
        self.mods_dll = None
        if os.path.isfile(os.path.join(model_path, "x86_64", "libnrnmech.so")):
            self.mods_dll = os.path.join(model_path, "x86_64", "libnrnmech.so")
        elif os.path.isfile(os.path.join(model_path, "./nrnmech.dll")):
            self.mods_dll = os.path.join(model_path, "./nrnmech.dll")
        else:
            logging.error(f"Missing dll/so mechanisms in {model_path}")
        if "CaDynamics_E2" not in dir(h) and self.mods_dll is not None:
            h.nrn_load_dll(self.mods_dll)
            logging.info(f"Loaded {self.mods_dll}")

    @staticmethod
    def __get_markers(morph_full_path):
        """Use blueBrain morphio package to parse asc (tbd not sure how to do the same in neuron loading)

        :param morph_full_path:
        :return: dict with markers data
        """
        def to_dict(cls):
            return dict([(k, getattr(cls, k)) for k in dir(cls) if not k.startswith("__")])
        curr_cell_data = Morphology(morph_full_path)
        return [to_dict(marker) for marker in curr_cell_data.markers]

    def __set_trunk_oblique_from_morph(self, diam_diff_threshold_um=0.2, stop_trunk_at=None):
        """For clear tree like shape we can use morph to split apical
        :return:
        """
        apical_children_of_soma = [c for c in self.soma.children() if "apic" in c.name()]
        apical_children_of_soma = [s for s in apical_children_of_soma if len(s.children()) > 0]
        if len(apical_children_of_soma) != 1:
            logging.error("Can't find single apical trunk in soma children: {0}".format(self.soma.children()))
            return
        child = apical_children_of_soma[0]
        trunk = [child]
        while max(abs(np.diff([c.diam for c in child.children()]))) > diam_diff_threshold_um:
            child = max(child.children(), key=lambda c: c.diam)
            trunk.append(child)
            if len(child.children()) <= 1 or (stop_trunk_at is not None and stop_trunk_at in child.name()):
                break
        obliques = []  # choose those who are on the path but not on trunk - todo prone to errors in trunk detection
        for sec in trunk[:-1]:
            added_obliques = [c for c in sec.children() if c not in trunk]
            obliques.extend(added_obliques)
            for o in added_obliques:
                obliques.extend(o.subtree())
        self.trunk_list = trunk
        self.oblique_list = list(set(obliques))
        self.tuft_list = [tuft for tuft in trunk[-1].subtree() if tuft not in trunk]  # remove trunk[-1]

    @staticmethod
    def __add_few_spines(sref_list, x_vec, neck_diam, neck_len, spine_head_area, ra, cm, rm, e_pas,
                         head_psd_percentage=0.1):
        def create(num, is_head=True, is_psd=False):
            if is_head and not is_psd:
                sp = h.Section(name=f"spine_head{num}")
                sp.L = (1 - head_psd_percentage) * L_head
                sp.diam = diam_head
            elif is_head and is_psd:
                sp = h.Section(name=f"spine_head_psd{num}")
                sp.L = head_psd_percentage * L_head
                sp.diam = diam_head
            else:
                sp = h.Section(name=f"spine_neck{num}")
                sp.L = neck_len
                sp.diam = neck_diam
            sp.insert('pas')
            sp.g_pas = 1 / rm  # 1/Rm - Rm ohm*cm^2
            sp.e_pas = e_pas
            sp.cm = cm
            sp.Ra = ra  # ohm*cm
            return sp

        L_head = 2 * np.sqrt(spine_head_area / 4 / np.pi)  # sphere has the same surface area as cylinder with L=diam
        diam_head = L_head
        spines = []
        spine_psds = []

        for i, sec in enumerate(sref_list):
            for j, shaft_x in enumerate(x_vec[i]):
                sp_head = create(j, is_head=True, is_psd=False)
                sp_head_psd = create(j, is_head=True, is_psd=True)
                sp_neck = create(j, is_head=False)
                spine_psds.append(sp_head_psd)
                spines.append(sp_neck)  # 2j
                spines.append(sp_head)  # 2j + 1
                sp_head_psd.connect(sp_head(1), 0)
                sp_head.connect(sp_neck(1), 0)
                sp_neck.connect(sec(shaft_x), 0)
                print(sp_head(1), " connect to begin of ", sp_head_psd, " with diam ", sp_head_psd.diam, " length ",
                      sp_head_psd.L)
                print(sp_neck(1), " connect to begin of ", sp_head, " with diam ", sp_head.diam, " length ",
                      sp_head.L)
                print(sec(shaft_x), " connect to begin of ", sp_neck, " with diam ", sp_neck.diam, " length ",
                      sp_neck.L)
        return spines, spine_psds
