#!/usr/bin/env python3
"""Simulated experiment: triangle test for synchronization holonomy."""
import numpy as np
from tep_model import save

def run():
    print("Building simulated experiment for Triangle Test...")
    
    # Simulate a network of 3 clocks: Ground 1, Ground 2, Satellite.
    # We simulate the time-of-flight of photons between them.
    
    # Parameters
    # Loop duration: ~12 hours for a MEO orbit.
    # Let's just generate a synthetic time series of the holonomy signal.
    # The signal is expected to be a diurnal/orbital modulation of amplitude ~ 0.1 fs.
    
    t_hours = np.linspace(0, 24, 100) # 24 hours of data
    
    # True TEP signal: H_resid ~ H0 * cos(omega_orbit * t + phase)
    amplitude_fs = 0.1
    omega_orbit = 2 * np.pi / 12.0 # 12-hour orbit
    
    # Add some GR "noise" (known effects that need to be subtracted, e.g. Sagnac, Shapiro)
    # GR effects are typically on the order of nanoseconds to microseconds!
    # Sagnac effect for MEO is ~ 100 ns.
    # We assume the pipeline subtracts this perfectly, leaving only clock noise + TEP signal.
    
    # Clock noise: white noise + random walk
    # Best optical clocks have fractional stability ~ 1e-18 at 10^4 s.
    # Over 12 hours (~4e4 s), timing error ~ 4e4 * 1e-18 = 4e-14 s = 40 fs.
    # Wait, the expected signal is 0.1 fs!
    # So single-shot measurement is heavily noise-dominated.
    # This requires integration over many orbits.
    
    noise_amplitude = 40.0 # fs
    noise = np.random.normal(0, noise_amplitude, len(t_hours))
    
    signal = amplitude_fs * np.cos(omega_orbit * t_hours)
    
    measured = signal + noise
    
    result = {
        't_hours': t_hours.tolist(),
        'signal_fs': signal.tolist(),
        'noise_fs': noise.tolist(),
        'measured_fs': measured.tolist(),
        'conclusion': 'Synthetic data generated. Single-orbit signal (0.1 fs) is buried in optical clock noise (40 fs). Detection requires phase-locked integration over ~1.6e5 orbits, or ~200 years with a single triangle, OR improvement in clock stability to 1e-19, OR a larger disformal parameter. Alternatively, interplanetary baselines (AU) increase the signal to picoseconds, making it immediately detectable.'
    }
    
    save('step_15_simulated_experiment.json', result)
    print(result['conclusion'])
    return result

if __name__ == '__main__':
    run()
