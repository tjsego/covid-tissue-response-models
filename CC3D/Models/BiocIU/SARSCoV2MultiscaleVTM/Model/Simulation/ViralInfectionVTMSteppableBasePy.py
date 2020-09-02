# This is a library of steppable classes for the viral infection modeling project using CompuCell3D
# by the Biocomplexity Institute at Indiana University using CompuCell3D

import os
import sys

# Import project libraries
sys.path.append(os.path.dirname(__file__))
import ViralInfectionVTMLib

# Import toolkit
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from nCoVToolkit.nCoVSteppableBase import nCoVSteppableBase
from nCoVToolkit import nCoVUtils


class ViralInfectionVTMSteppableBasePy(nCoVSteppableBase):

    vr_model_name = ViralInfectionVTMLib.vr_model_name
    _viral_replication_model_string_gen = ViralInfectionVTMLib.viral_replication_model_string

    def __init__(self, frequency=1):
        nCoVSteppableBase.__init__(self, frequency)

    def start(self):
        pass

    def step(self, mcs):
        pass

    def finish(self):
        pass

    def viral_replication_model_string(self, *args, **kwargs):
        """
        Antimony model string generator for viral replication model; can be set with set_viral_replication_model(), and
        subclasses can override
        :param args:
        :param kwargs:
        :return {str}: Antimony model string for viral replication model
        """
        return ViralInfectionVTMSteppableBasePy._viral_replication_model_string_gen(*args, **kwargs)

    @staticmethod
    def set_viral_replication_model(_fnc, _name: str = ViralInfectionVTMLib.vr_model_name):
        """
        Sets the viral replication model
        :param _fnc: Antimony model string generator
        :param _name: name of Antimony model
        :return: None
        """
        ViralInfectionVTMSteppableBasePy._viral_replication_model_string_gen = _fnc
        ViralInfectionVTMSteppableBasePy.vr_model_name = _name

    def load_viral_replication_model(self, *args, **kwargs):
        """
        Loads viral replication model for a cell; subclasses can override
        :param args:
        :param kwargs:
        :return: None
        """
        assert 'cell' in kwargs, 'Specify viral replication model cell with keyword cell'
        assert 'vr_step_size' in kwargs, 'Specify viral replication model step size with keyword vr_step_size'
        return self._load_viral_replication_model(**kwargs)

    def _load_viral_replication_model(self, cell, vr_step_size, unpacking_rate=0, replicating_rate=0, r_half=0,
                                      translating_rate=0, packing_rate=0, secretion_rate=0):
        """
        Loads viral replication model for a cell
        :param cell: cell for which the viral replication model is loaded
        :param vr_step_size: Antimony/SBML model step size
        :param unpacking_rate: model unpacking rate
        :param replicating_rate: model replicating rate
        :param r_half: Value of R at which the replication rate is half max
        :param translating_rate: model translating rate
        :param packing_rate: model packing rate
        :param secretion_rate: model secretion rate
        :return: None
        """
        if cell.dict[ViralInfectionVTMLib.vrl_key]:
            self.delete_sbml_from_cell(ViralInfectionVTMSteppableBasePy.vr_model_name, cell)

        # Generate Antimony model string
        model_string = self.viral_replication_model_string(
            unpacking_rate, replicating_rate, r_half, translating_rate, packing_rate, secretion_rate)
        self.add_antimony_to_cell(model_string=model_string,
                                  model_name=ViralInfectionVTMSteppableBasePy.vr_model_name,
                                  cell=cell,
                                  step_size=vr_step_size)
        cell.dict[ViralInfectionVTMLib.vrl_key] = True
        ViralInfectionVTMLib.enable_viral_secretion(cell=cell, enable=cell.type == self.VIRUSRELEASING)

    def new_cell_in_time(self, cell_type, mcs=None):
        """
        Add cell and record MCS
        :param cell_type: type id of cell (e.g., for cell.type)
        :param mcs: step when cell is created; defaults from steppable mcs attribute
        :return: new cell instance
        """
        cell = self.new_cell(cell_type)
        if mcs is None:
            if self.mcs < 0:
                mcs = 0
            else:
                mcs = self.mcs

        cell.dict[ViralInfectionVTMLib.new_cell_mcs_key] = mcs
        return cell

    def new_uninfected_cell_in_time(self, mcs=None):
        """
        Add an uninfected cell with default initial configuration and record MCS
        :param mcs: step when cell is created; defaults from steppable mcs attribute
        :return: new cell instance of immune cell type
        """
        cell = self.new_cell_in_time(self.UNINFECTED, mcs)
        cell.dict[ViralInfectionVTMLib.vrl_key] = False
        cell.dict['Survived'] = False
        return cell

    def new_immune_cell_in_time(self, ck_production, ck_consumption, mcs=None, activated=False):
        """
        Add an immune cell with default initial configuration and record MCS
        :param ck_production: cytokine production rate
        :param ck_consumption: cytokine consumption rate
        :param mcs: step when cell is created; defaults from steppable mcs attribute
        :param activated: flag for immune cell being naive or activated
        :return: new cell instance of immune cell type
        """
        cell = self.new_cell_in_time(self.IMMUNECELL, mcs)
        # cyttokine params
        cell.dict['ck_production'] = ck_production  # TODO: replace secretion by hill
        cell.dict['ck_consumption'] = ck_consumption  # TODO: replace by hill
        cell.dict['activated'] = activated
        cell.dict['tot_ck_upt'] = 0
        return cell

    def total_seen_field(self, field, cell, estimate=True):
        """
        Calculates total value of field in the cell.
        :param field: the field to be looked
        :param cell: the cell
        :param estimate: when true assumes homogeneous field. false is slower
        :return: calculated total field value
        """
        if estimate:
            tot_field = field[cell.xCOM, cell.yCOM, cell.zCOM] * cell.volume
        else:
            tot_field = 0
            for ptd in self.get_cell_pixel_list(cell):
                tot_field += field[ptd.pixel.x, ptd.pixel.y, ptd.pixel.z]

        return tot_field

    def kill_cell(self, cell):
        """
        Model-specific cell death routines
        :param cell: cell to kill
        :return: None
        """
        cell.type = self.DYING

        # Remove viral replication model: no model for dead cell type
        self.remove_viral_replication_model(cell=cell)

    def remove_viral_replication_model(self, cell):
        """
        Removes viral replication model for a cell
        :param cell: cell for which to remove the viral replication model
        :return: None
        """
        if cell.dict[ViralInfectionVTMLib.vrl_key]:
            self.delete_sbml_from_cell(ViralInfectionVTMSteppableBasePy.vr_model_name, cell)
            cell.dict[ViralInfectionVTMLib.vrl_key] = False
