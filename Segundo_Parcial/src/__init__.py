"""
Módulo para Recocido Simulado aplicado a diseño de materiales magnéticos.
Segundo Parcial - Cadenas de Markov Aplicadas
"""

from .morse_potential import (
    calcular_potencial_morse,
    calcular_energia_total,
    calcular_energia_total_eficiente,
    calcular_delta_energia_swap
)

from .simulated_annealing import (
    recocido_simulado,
    esquema_enfriamiento_exponencial,
    esquema_enfriamiento_lineal,
    esquema_enfriamiento_logaritmico
)

from .lattice_2d import (
    crear_red_2d_4x4,
    crear_red_2d_10x10,
    visualizar_red_2d
)

from .lattice_3d import (
    cargar_posiciones_3d,
    visualizar_red_3d
)

__all__ = [
    'calcular_potencial_morse',
    'calcular_energia_total',
    'calcular_energia_total_eficiente',
    'calcular_delta_energia_swap',
    'recocido_simulado',
    'esquema_enfriamiento_exponencial',
    'esquema_enfriamiento_lineal',
    'esquema_enfriamiento_logaritmico',
    'crear_red_2d_4x4',
    'crear_red_2d_10x10',
    'visualizar_red_2d',
    'cargar_posiciones_3d',
    'visualizar_red_3d'
]