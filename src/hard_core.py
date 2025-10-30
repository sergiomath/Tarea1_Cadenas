"""
Implementación del Gibbs Sampler para el modelo Hard-Core
"""
import numpy as np

def obtener_vecinos(i, j, K):
    """Retorna las coordenadas de los vecinos de la celda (i,j) en una rejilla K×K"""
    vecinos = []
    for di, dj in [(-1,0), (1,0), (0,-1), (0,1)]:
        ni, nj = i + di, j + dj
        if 0 <= ni < K and 0 <= nj < K:
            vecinos.append((ni, nj))
    return vecinos

def es_configuracion_factible(config):
    """Verifica si una configuración satisface la restricción Hard-Core"""
    K = config.shape[0]
    for i in range(K):
        for j in range(K):
            if config[i,j] == 1:
                for ni, nj in obtener_vecinos(i, j, K):
                    if config[ni, nj] == 1:
                        return False
    return True

def contar_particulas(config):
    """Cuenta el número total de partículas en la configuración"""
    return np.sum(config)

def gibbs_sampler_hard_core(K, T, semilla=None):
    """
    Gibbs Sampler para el modelo Hard-Core

    Parámetros:
    - K: tamaño de la rejilla (K×K)
    - T: número de iteraciones
    - semilla: semilla para reproducibilidad

    Retorna:
    - Configuración final después de T iteraciones
    """
    if semilla is not None:
        np.random.seed(semilla)

    # Inicialización: configuración vacía
    config = np.zeros((K, K), dtype=int)

    # Iteraciones del Gibbs Sampler
    for _ in range(T):
        # Seleccionar sitio aleatorio
        i, j = np.random.randint(0, K, size=2)

        # Verificar vecinos
        vecinos = obtener_vecinos(i, j, K)
        tiene_vecino_ocupado = any(config[ni, nj] == 1 for ni, nj in vecinos)

        # Actualizar sitio
        if tiene_vecino_ocupado:
            config[i, j] = 0  # Debe estar vacío
        else:
            config[i, j] = np.random.choice([0, 1])  # Uniforme entre {0,1}

    return config
