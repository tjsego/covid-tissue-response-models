import math


class _TwoStateOscillator:
    def __init__(self, on_time=1, off_time=None, init_state=True, init_time=None):
        """
        Base class for two-state periodic signal in time
        :param on_time: period of high state
        :param off_time: period of low state; defaults to *on_time*
        :param init_state: initial state; defaults to high
        :param init_time: period of initial state; defaults to *on_time*
        """
        self._on_time = on_time
        if off_time is not None:
            self._off_time = off_time
        else:
            self._off_time = self._on_time

        if init_time is not None:
            self._init_time = init_time
        else:
            self._init_time = self._on_time

        self._flip_time = self._init_time
        self._init_state = init_state
        self._state = self._init_state
        self.initial_period = True

    def set_on_time(self, _on_time):
        self._on_time = _on_time

    def set_off_time(self, _off_time):
        self._off_time = _off_time

    def set_init_state(self, _state: bool):
        self._init_state = _state
        if self.initial_period:
            self.set_state(_state)

    def set_init_time(self, _init_time):
        self._init_time = _init_time

    def set_state(self, _state: bool):
        self._state = _state

    def _update_state(self, _time):
        if self.initial_period:
            if _time < self._init_time:
                self._state = self._init_state
            else:
                self._state = not self._init_state
                self.initial_period = False
                self._flip_time = self._init_time
        else:
            if self._state:
                if self._flip_time + self._on_time < _time:
                    self._state = False
                    self._flip_time += self._on_time
            else:
                if self._flip_time + self._off_time < _time:
                    self._state = True
                    self._flip_time += self._off_time

    def get_signal(self, _time):
        raise NotImplementedError


class SawtoothOscillator(_TwoStateOscillator):
    def __init__(self, on_time=1, off_time=None, init_state=True, init_time=None):
        """
        Implements a sawtooth periodic signal in time
        :param on_time: period of high state
        :param off_time: period of low state; defaults to *on_time*
        :param init_state: initial state; defaults to high
        :param init_time: period of initial state; defaults to *on_time*
        """
        super().__init__(on_time=on_time, off_time=off_time, init_state=init_state, init_time=init_time)

    def get_signal(self, _time):
        self._update_state(_time)
        if self._state:
            return 1.0
        return 0.0


class PeriodicExpDecay(_TwoStateOscillator):
    def __init__(self, on_time=1, init_time=None, decay_rate=1):
        """
        Implements a periodic exponential decay
        :param on_time: period of signal
        :param init_time: period before first instance of signal
        :param decay_rate: decay rate of signal
        """
        super().__init__(on_time=on_time, init_time=init_time)
        self._decay_rate = decay_rate

    def set_decay_rate(self, _rate):
        self._decay_rate = _rate

    def get_signal(self, _time):
        self._update_state(_time)
        if self.initial_period:
            return 0.0
        return math.exp(-self._decay_rate * (_time - self._flip_time))
