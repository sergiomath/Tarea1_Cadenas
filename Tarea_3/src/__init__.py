"""
Módulo de simulación del Modelo de Ising.
Tarea 3: Muestreo MCMC vs Simulación Perfecta.
"""

from .ising_model import (
    IsingLattice,
    GibbsSampler,
    ProppWilson,
    estimate_magnetization,
    run_mcmc_experiment,
    run_propp_wilson_experiment,
)

__all__ = [
    'IsingLattice',
    'GibbsSampler',
    'ProppWilson',
    'estimate_magnetization',
    'run_mcmc_experiment',
    'run_propp_wilson_experiment',
]
