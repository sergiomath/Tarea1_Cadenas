"""
Funciones para análisis estadístico de los modelos
"""
import numpy as np
from .hard_core import gibbs_sampler_hard_core, contar_particulas

def calcular_estadisticas(muestras):
    """
    Calcula estadísticas descriptivas de un conjunto de muestras

    Parámetros:
    - muestras: array con valores numéricos

    Retorna:
    - dict con media, mediana, desviación estándar, min, max
    """
    return {
        'media': np.mean(muestras),
        'mediana': np.median(muestras),
        'std': np.std(muestras),
        'min': np.min(muestras),
        'max': np.max(muestras)
    }

def analizar_multiple_K(K_valores, T, n_muestras=100):
    """
    Analiza el modelo Hard-Core para múltiples tamaños de rejilla

    Parámetros:
    - K_valores: lista de tamaños de rejilla a analizar
    - T: número de iteraciones del Gibbs Sampler
    - n_muestras: número de muestras independientes por cada K

    Retorna:
    - dict con resultados por cada K
    """
    resultados = {}

    for K in K_valores:
        particulas = []
        for i in range(n_muestras):
            config = gibbs_sampler_hard_core(K, T, semilla=i)
            n_particulas = contar_particulas(config)
            particulas.append(n_particulas)

        resultados[K] = {
            'particulas': particulas,
            'estadisticas': calcular_estadisticas(particulas)
        }

    return resultados

def crear_tabla_estadisticas(resultados):
    """
    Crea una tabla formateada con las estadísticas por tamaño K

    Parámetros:
    - resultados: diccionario retornado por analizar_multiple_K

    Retorna:
    - string con tabla formateada
    """
    lineas = []
    lineas.append("K\tMedia\tMediana\tStd\tMin\tMax\tDensidad")
    lineas.append("-" * 60)

    for K in sorted(resultados.keys()):
        stats = resultados[K]['estadisticas']
        densidad = stats['media'] / (K * K)
        lineas.append(f"{K}\t{stats['media']:.2f}\t{stats['mediana']:.2f}\t"
                     f"{stats['std']:.2f}\t{stats['min']}\t{stats['max']}\t{densidad:.4f}")

    return "\n".join(lineas)
