"""
Implementación del Gibbs Sampler para el modelo de q-coloraciones
"""
import numpy as np

try:
    from .hard_core import obtener_vecinos
except ImportError:
    from hard_core import obtener_vecinos

def es_coloracion_propia(config, q):
    """Verifica si una configuración es una q-coloración propia"""
    K = config.shape[0]
    for i in range(K):
        for j in range(K):
            color_actual = config[i, j]
            if not (0 <= color_actual < q):
                return False
            for ni, nj in obtener_vecinos(i, j, K):
                if config[ni, nj] == color_actual:
                    return False
    return True

def contar_colores(config, q):
    """Cuenta cuántas celdas tienen cada color"""
    conteo = np.zeros(q, dtype=int)
    for color in range(q):
        conteo[color] = np.sum(config == color)
    return conteo

def gibbs_sampler_q_coloraciones(K, q, T, semilla=None):
    """
    Gibbs Sampler para el modelo de q-coloraciones

    Parámetros:
    - K: tamaño de la rejilla (K×K)
    - q: número de colores disponibles
    - T: número de iteraciones
    - semilla: semilla para reproducibilidad

    Retorna:
    - Configuración final después de T iteraciones
    """
    if semilla is not None:
        np.random.seed(semilla)

    # Inicialización: coloración tipo tablero de ajedrez para q=2
    # Para q>2, inicialización aleatoria válida
    config = np.zeros((K, K), dtype=int)
    for i in range(K):
        for j in range(K):
            config[i, j] = (i + j) % min(q, 2)

    # Iteraciones del Gibbs Sampler
    for _ in range(T):
        # Seleccionar sitio aleatorio
        i, j = np.random.randint(0, K, size=2)

        # Obtener colores de vecinos
        vecinos = obtener_vecinos(i, j, K)
        colores_vecinos = set(config[ni, nj] for ni, nj in vecinos)

        # Colores disponibles (no usados por vecinos)
        colores_disponibles = [c for c in range(q) if c not in colores_vecinos]

        # Seleccionar color uniformemente de los disponibles
        if colores_disponibles:
            config[i, j] = np.random.choice(colores_disponibles)

    return config
