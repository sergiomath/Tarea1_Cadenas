"""
Implementación del Potencial de Morse para cálculo de energías en sistemas atómicos.
"""

import numpy as np
from typing import Dict, List, Tuple, Optional

# Intentar importar numba, si no está disponible, crear un decorador dummy
try:
    from numba import jit
except ImportError:
    def jit(nopython=True):
        def decorator(func):
            return func
        return decorator


# Parámetros del Potencial de Morse
PARAMETROS_MORSE = {
    ('Fe', 'Fe'): {'D0': 0.764, 'alpha': 1.5995, 'r0': 2.7361},
    ('Fe', 'Nd'): {'D0': 0.6036, 'alpha': 1.6458, 'r0': 3.188},
    ('Nd', 'Fe'): {'D0': 0.6036, 'alpha': 1.6458, 'r0': 3.188},  # Simétrico
    ('Nd', 'Nd'): {'D0': 0.312, 'alpha': 0.945, 'r0': 4.092},
    ('Fe', 'Ti'): {'D0': 0.8162, 'alpha': 1.448, 'r0': 2.914},
    ('Ti', 'Fe'): {'D0': 0.8162, 'alpha': 1.448, 'r0': 2.914},  # Simétrico
    ('Nd', 'Ti'): {'D0': 0.4964, 'alpha': 1.4401, 'r0': 3.4309},
    ('Ti', 'Nd'): {'D0': 0.4964, 'alpha': 1.4401, 'r0': 3.4309},  # Simétrico
    ('Ti', 'Ti'): {'D0': 0.6540, 'alpha': 1.2118, 'r0': 3.3476},
}


def calcular_distancia(pos1: np.ndarray, pos2: np.ndarray) -> float:
    """
    Calcula la distancia euclidiana entre dos posiciones.

    Args:
        pos1: Posición del primer átomo [x, y] o [x, y, z]
        pos2: Posición del segundo átomo [x, y] o [x, y, z]

    Returns:
        Distancia euclidiana
    """
    return np.sqrt(np.sum((pos1 - pos2) ** 2))


def calcular_potencial_morse(r: float, tipo1: str, tipo2: str) -> float:
    """
    Calcula el potencial de Morse entre dos átomos.

    U(r) = D0 * [exp(-2α(r-r0)) - 2*exp(-α(r-r0))]

    Args:
        r: Distancia entre los átomos
        tipo1: Tipo del primer átomo ('Fe', 'Nd', 'Ti')
        tipo2: Tipo del segundo átomo ('Fe', 'Nd', 'Ti')

    Returns:
        Energía de interacción
    """
    params = PARAMETROS_MORSE[(tipo1, tipo2)]
    D0 = params['D0']
    alpha = params['alpha']
    r0 = params['r0']

    # Potencial de Morse: U(r) = D0 * [e^(-2α(r-r0)) - 2*e^(-α(r-r0))]
    exp_single = np.exp(-alpha * (r - r0))
    exp_double = np.exp(-2 * alpha * (r - r0))
    return D0 * (exp_double - 2 * exp_single)


def calcular_energia_total(
    posiciones: Dict[str, np.ndarray],
    tipos: Dict[str, str]
) -> float:
    """
    Calcula la energía total del sistema sumando todas las interacciones de pares.

    Args:
        posiciones: Diccionario {id_atomo: posicion_np_array}
        tipos: Diccionario {id_atomo: tipo_str}

    Returns:
        Energía total del sistema
    """
    energia_total = 0.0
    atomos = list(posiciones.keys())

    for i, atomo_i in enumerate(atomos):
        for atomo_j in atomos[i+1:]:
            r = calcular_distancia(posiciones[atomo_i], posiciones[atomo_j])
            energia_total += calcular_potencial_morse(
                r, tipos[atomo_i], tipos[atomo_j]
            )

    return energia_total


def calcular_energia_total_eficiente(
    posiciones_array: np.ndarray,
    tipos_array: np.ndarray,
    parametros_matriz: np.ndarray
) -> float:
    """
    Versión optimizada del cálculo de energía usando arrays de NumPy.

    Args:
        posiciones_array: Array (n_atomos, dims) con posiciones
        tipos_array: Array (n_atomos,) con índices de tipos (0=Fe, 1=Nd, 2=Ti)
        parametros_matriz: Matriz (3, 3, 3) con parámetros [D0, alpha, r0]

    Returns:
        Energía total del sistema
    """
    n_atomos = len(posiciones_array)
    energia_total = 0.0

    for i in range(n_atomos):
        for j in range(i + 1, n_atomos):
            r = np.linalg.norm(posiciones_array[i] - posiciones_array[j])
            tipo_i = tipos_array[i]
            tipo_j = tipos_array[j]

            D0 = parametros_matriz[tipo_i, tipo_j, 0]
            alpha = parametros_matriz[tipo_i, tipo_j, 1]
            r0 = parametros_matriz[tipo_i, tipo_j, 2]

            exp_term = np.exp(-alpha * (r - r0))
            energia_total += D0 * (exp_term * exp_term - 2 * exp_term)

    return energia_total


def calcular_delta_energia_swap(
    posiciones: Dict[str, np.ndarray],
    tipos: Dict[str, str],
    atomo_ti: str,
    atomo_fe: str
) -> float:
    """
    Calcula el cambio en energía al intercambiar un átomo de Ti con uno de Fe.
    Método eficiente que solo recalcula las interacciones afectadas.

    Args:
        posiciones: Diccionario de posiciones
        tipos: Diccionario de tipos
        atomo_ti: ID del átomo de Ti a intercambiar
        atomo_fe: ID del átomo de Fe a intercambiar

    Returns:
        Cambio en energía (ΔE = E_nuevo - E_actual)
    """
    # Energía antes del intercambio (solo interacciones afectadas)
    energia_antes = 0.0
    # Energía después del intercambio
    energia_despues = 0.0

    pos_ti = posiciones[atomo_ti]
    pos_fe = posiciones[atomo_fe]

    # Para todos los otros átomos, calcular cambio en interacciones
    for atomo, pos in posiciones.items():
        if atomo == atomo_ti or atomo == atomo_fe:
            continue

        tipo_atomo = tipos[atomo]

        # Interacción con Ti en su posición actual
        r_ti = calcular_distancia(pos_ti, pos)
        energia_antes += calcular_potencial_morse(r_ti, 'Ti', tipo_atomo)
        # Después del swap, Fe estará en esta posición
        energia_despues += calcular_potencial_morse(r_ti, 'Fe', tipo_atomo)

        # Interacción con Fe en su posición actual
        r_fe = calcular_distancia(pos_fe, pos)
        energia_antes += calcular_potencial_morse(r_fe, 'Fe', tipo_atomo)
        # Después del swap, Ti estará en esta posición
        energia_despues += calcular_potencial_morse(r_fe, 'Ti', tipo_atomo)

    # Interacción directa Ti-Fe antes y después
    r_ti_fe = calcular_distancia(pos_ti, pos_fe)
    energia_antes += calcular_potencial_morse(r_ti_fe, 'Ti', 'Fe')
    energia_despues += calcular_potencial_morse(r_ti_fe, 'Fe', 'Ti')  # Mismo resultado

    return energia_despues - energia_antes


def preparar_parametros_matriz() -> np.ndarray:
    """
    Prepara una matriz de parámetros para cálculos eficientes.

    Returns:
        Matriz (3, 3, 3) con parámetros [D0, alpha, r0]
        Índices: 0=Fe, 1=Nd, 2=Ti
    """
    matriz = np.zeros((3, 3, 3))
    tipo_a_indice = {'Fe': 0, 'Nd': 1, 'Ti': 2}

    for (tipo1, tipo2), params in PARAMETROS_MORSE.items():
        i = tipo_a_indice[tipo1]
        j = tipo_a_indice[tipo2]
        matriz[i, j, 0] = params['D0']
        matriz[i, j, 1] = params['alpha']
        matriz[i, j, 2] = params['r0']

    return matriz