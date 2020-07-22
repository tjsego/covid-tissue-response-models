# When using Sawtooth
use_sawtooth = False  # Set to True to use sawtooth
if use_sawtooth:
    on_time = 24 * 60 * 60  # period of high state, in seconds
    off_time = 8 * 60 * 60  # period of low state, in seconds
    init_state = True  # initial state; defaults to high
    init_time = 24 * 60 * 60  # period of initial state, in seconds

# When using exponential decay
use_exp_decay = True  # Set to True to use exponential decay
if use_exp_decay:
    on_time = 8 * 60 * 60  # period of dose, in seconds
    init_time = 24 * 60 * 60  # time before first dose, in seconds
    decay_rate = 1 / (4 * 60 * 60)  # Decay rate of dose, in 1/seconds
