import os
import sys

# Import simulation stuff
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), f'Simulation'))
from ViralInfectionVTMModelInputs import *
from ViralInfectionVTMSteppableBasePy import ViralInfectionVTMSteppableBasePy
import ViralInfectionVTMLib
from nCoVToolkit import SigGens
from Dosing import DosingInputs

if DosingInputs.use_sawtooth:
    r_max_oscillator = SigGens.SawtoothOscillator()
elif DosingInputs.use_exp_decay:
    r_max_oscillator = SigGens.PeriodicExpDecay()
else:
    raise EnvironmentError('Please select an oscillator in DosingInputs.py')


def update_replicating_rate(cell, replicating_rate):
    getattr(cell.sbml, ViralInfectionVTMLib.vr_model_name)['replicating_rate'] = replicating_rate


class DosingSteppableSawtooth(ViralInfectionVTMSteppableBasePy):
    """
    Implements randomly varying viral susceptibility in space by setting surface receptors to zero
    """
    def __init__(self, frequency=1):
        ViralInfectionVTMSteppableBasePy.__init__(self, frequency)

        # Initialize oscillator
        r_max_oscillator.set_on_time(DosingInputs.on_time)
        r_max_oscillator.set_off_time(DosingInputs.off_time)
        r_max_oscillator.set_init_state(DosingInputs.init_state)
        r_max_oscillator.set_init_time(DosingInputs.init_time)

        # replicating rate values are stored here and tested for if updates to cell models need to be made
        self._last_replicating_rate = self.get_replicating_rate(0)

    def step(self, mcs):
        # Check for update to replicating rate and update if necessary
        current_replicating_rate = self.get_replicating_rate(mcs * s_to_mcs)
        if current_replicating_rate != self._last_replicating_rate:
            for cell in self.cell_list_by_type(self.UNINFECTED, self.INFECTED, self.VIRUSRELEASING):
                update_replicating_rate(cell, current_replicating_rate)
        self._last_replicating_rate = current_replicating_rate

    def get_replicating_rate(self, _time):
        return replicating_rate * r_max_oscillator.get_signal(_time)


class DosingSteppableExpDecay(ViralInfectionVTMSteppableBasePy):
    """
    Implements randomly varying viral susceptibility in space by setting surface receptors to zero
    """
    def __init__(self, frequency=1):
        ViralInfectionVTMSteppableBasePy.__init__(self, frequency)

        # Initialize oscillator
        r_max_oscillator.set_on_time(DosingInputs.on_time)
        r_max_oscillator.set_init_time(DosingInputs.init_time)
        r_max_oscillator.set_decay_rate(DosingInputs.decay_rate)

    def step(self, mcs):
        current_replicating_rate = self.get_replicating_rate(mcs * s_to_mcs)
        for cell in self.cell_list_by_type(self.INFECTED, self.VIRUSRELEASING):
            update_replicating_rate(cell, current_replicating_rate)

    def get_replicating_rate(self, _time):
        return replicating_rate * (1 - r_max_oscillator.get_signal(_time))


if DosingInputs.use_sawtooth:
    class DosingSteppable(DosingSteppableSawtooth):
        def __init__(self, frequency=1):
            DosingSteppableSawtooth.__init__(self, frequency)
elif DosingInputs.use_exp_decay:
    class DosingSteppable(DosingSteppableExpDecay):
        def __init__(self, frequency=1):
            DosingSteppableExpDecay.__init__(self, frequency)
