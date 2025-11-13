"""
Módulo de conteo aproximado con MCMC para Tarea 2.

Contiene implementaciones de:
- Q-coloraciones en grafos reticulares
- Modelo Hard-Core
- Algoritmos basados en Teorema 9.1 (FPRAS)
"""

from .mcmc_counting import (
    MCMCConfig,
    LatticeGraph,
    QColoringMCMC,
    QColoringApproximation,
    HardCoreMCMC,
    HardCoreApproximation,
    exact_q_colorings_small,
    exact_chromatic_polynomial,
    exact_hardcore_small,
    exact_hardcore_count,
    run_q_coloring_experiment,
    run_hardcore_experiment,
)

__all__ = [
    'MCMCConfig',
    'LatticeGraph',
    'QColoringMCMC',
    'QColoringApproximation',
    'HardCoreMCMC',
    'HardCoreApproximation',
    'exact_q_colorings_small',
    'exact_chromatic_polynomial',
    'exact_hardcore_small',
    'exact_hardcore_count',
    'run_q_coloring_experiment',
    'run_hardcore_experiment',
]
