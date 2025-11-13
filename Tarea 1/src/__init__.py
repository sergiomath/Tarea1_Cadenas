"""
Módulos para Tarea 1 - Muestreo de Gibbs en Modelos Estocásticos

Módulos disponibles:
- hard_core: Implementación del Gibbs Sampler para modelo Hard-Core
- q_coloraciones: Implementación del Gibbs Sampler para q-coloraciones
- estadisticas: Funciones para análisis estadístico
- visualizacion: Funciones para visualización de resultados
"""

from .hard_core import (
    gibbs_sampler_hard_core,
    es_configuracion_factible,
    obtener_vecinos,
    contar_particulas
)

from .q_coloraciones import (
    gibbs_sampler_q_coloraciones,
    es_coloracion_propia,
    contar_colores
)

from .estadisticas import (
    calcular_estadisticas,
    analizar_multiple_K,
    crear_tabla_estadisticas
)

from .visualizacion import (
    visualizar_configuracion,
    graficar_histograma,
    graficar_escalamiento,
    graficar_distribucion_colores
)

__all__ = [
    'gibbs_sampler_hard_core',
    'es_configuracion_factible',
    'obtener_vecinos',
    'contar_particulas',
    'gibbs_sampler_q_coloraciones',
    'es_coloracion_propia',
    'contar_colores',
    'calcular_estadisticas',
    'analizar_multiple_K',
    'crear_tabla_estadisticas',
    'visualizar_configuracion',
    'graficar_histograma',
    'graficar_escalamiento',
    'graficar_distribucion_colores'
]
