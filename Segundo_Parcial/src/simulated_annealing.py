"""
Implementación del algoritmo de Recocido Simulado (Simulated Annealing).
"""

import numpy as np
from typing import Callable, List, Tuple, Dict, Optional
import copy
from dataclasses import dataclass
from tqdm import tqdm


@dataclass
class ResultadoSA:
    """Resultado del algoritmo de recocido simulado."""
    mejor_estado: any
    mejor_energia: float
    historia_energia: List[float]
    historia_temperatura: List[float]
    historia_aceptacion: List[float]
    n_iteraciones: int
    n_aceptadas: int


def recocido_simulado(
    estado_inicial: any,
    funcion_energia: Callable,
    generar_vecino: Callable,
    temperatura_inicial: float,
    esquema_enfriamiento: Callable,
    n_iteraciones_max: int = 10000,
    n_iteraciones_por_T: int = 100,
    temperatura_min: float = 1e-6,
    verbose: bool = True
) -> ResultadoSA:
    """
    Algoritmo de Recocido Simulado para optimización.

    Args:
        estado_inicial: Estado inicial del sistema
        funcion_energia: Función que calcula la energía de un estado
        generar_vecino: Función que genera un vecino del estado actual
        temperatura_inicial: Temperatura inicial T0
        esquema_enfriamiento: Función que actualiza la temperatura
        n_iteraciones_max: Número máximo de iteraciones totales
        n_iteraciones_por_T: Iteraciones a temperatura constante
        temperatura_min: Temperatura mínima para detener
        verbose: Mostrar progreso

    Returns:
        ResultadoSA con el mejor estado encontrado y estadísticas
    """
    # Inicialización
    estado_actual = copy.deepcopy(estado_inicial)
    energia_actual = funcion_energia(estado_actual)

    mejor_estado = copy.deepcopy(estado_actual)
    mejor_energia = energia_actual

    temperatura = temperatura_inicial

    # Historiales
    historia_energia = [energia_actual]
    historia_temperatura = [temperatura]
    historia_aceptacion = []

    n_iteraciones = 0
    n_aceptadas = 0
    n_rechazadas = 0

    # Barra de progreso
    if verbose:
        pbar = tqdm(total=n_iteraciones_max, desc="Recocido Simulado")

    while n_iteraciones < n_iteraciones_max and temperatura > temperatura_min:
        aceptadas_en_T = 0

        for _ in range(n_iteraciones_por_T):
            if n_iteraciones >= n_iteraciones_max:
                break

            # Generar vecino
            estado_vecino, movimiento_info = generar_vecino(estado_actual)

            # Calcular energía del vecino
            # Si la función devuelve delta directamente (más eficiente)
            if movimiento_info and 'delta_energia' in movimiento_info:
                delta_energia = movimiento_info['delta_energia']
                energia_vecino = energia_actual + delta_energia
            else:
                energia_vecino = funcion_energia(estado_vecino)
                delta_energia = energia_vecino - energia_actual

            # Criterio de Metropolis
            if delta_energia < 0:
                # Siempre aceptar si mejora
                aceptar = True
            else:
                # Aceptar con probabilidad exp(-ΔE/T)
                probabilidad_aceptacion = np.exp(-delta_energia / temperatura)
                aceptar = np.random.random() < probabilidad_aceptacion

            if aceptar:
                estado_actual = estado_vecino
                energia_actual = energia_vecino
                n_aceptadas += 1
                aceptadas_en_T += 1

                # Actualizar mejor solución
                if energia_actual < mejor_energia:
                    mejor_estado = copy.deepcopy(estado_actual)
                    mejor_energia = energia_actual
            else:
                n_rechazadas += 1

            historia_energia.append(energia_actual)
            n_iteraciones += 1

            if verbose:
                pbar.update(1)
                pbar.set_postfix({
                    'T': f'{temperatura:.4f}',
                    'E_mejor': f'{mejor_energia:.4f}',
                    'E_actual': f'{energia_actual:.4f}'
                })

        # Guardar tasa de aceptación en esta temperatura
        tasa_aceptacion = aceptadas_en_T / n_iteraciones_por_T if n_iteraciones_por_T > 0 else 0
        historia_aceptacion.append(tasa_aceptacion)

        # Actualizar temperatura
        temperatura = esquema_enfriamiento(temperatura, n_iteraciones)
        historia_temperatura.append(temperatura)

    if verbose:
        pbar.close()
        print(f"\nRecocido Simulado completado:")
        print(f"  Iteraciones: {n_iteraciones}")
        print(f"  Mejor energía: {mejor_energia:.6f}")
        print(f"  Tasa aceptación global: {n_aceptadas/n_iteraciones:.2%}")

    return ResultadoSA(
        mejor_estado=mejor_estado,
        mejor_energia=mejor_energia,
        historia_energia=historia_energia,
        historia_temperatura=historia_temperatura,
        historia_aceptacion=historia_aceptacion,
        n_iteraciones=n_iteraciones,
        n_aceptadas=n_aceptadas
    )


def esquema_enfriamiento_exponencial(alpha: float = 0.95) -> Callable:
    """
    Esquema de enfriamiento exponencial: T(k+1) = α * T(k)

    Args:
        alpha: Factor de enfriamiento (0 < α < 1)

    Returns:
        Función de enfriamiento
    """
    def enfriamiento(T_actual: float, iteracion: int) -> float:
        return alpha * T_actual
    return enfriamiento


def esquema_enfriamiento_lineal(delta_T: float = 0.01) -> Callable:
    """
    Esquema de enfriamiento lineal: T(k+1) = T(k) - δT

    Args:
        delta_T: Decremento de temperatura

    Returns:
        Función de enfriamiento
    """
    def enfriamiento(T_actual: float, iteracion: int) -> float:
        return max(0, T_actual - delta_T)
    return enfriamiento


def esquema_enfriamiento_logaritmico(c: float = 1.0) -> Callable:
    """
    Esquema de enfriamiento logarítmico: T(k) = c / log(k + 2)

    Args:
        c: Constante de escala

    Returns:
        Función de enfriamiento
    """
    def enfriamiento(T_actual: float, iteracion: int) -> float:
        return c / np.log(iteracion + 2)
    return enfriamiento


def calcular_temperatura_inicial(
    estado_inicial: any,
    funcion_energia: Callable,
    generar_vecino: Callable,
    probabilidad_aceptacion_inicial: float = 0.8,
    n_muestras: int = 100
) -> float:
    """
    Estima una temperatura inicial adecuada mediante muestreo.

    Args:
        estado_inicial: Estado inicial
        funcion_energia: Función de energía
        generar_vecino: Función generadora de vecinos
        probabilidad_aceptacion_inicial: Probabilidad deseada de aceptación
        n_muestras: Número de muestras para estimar

    Returns:
        Temperatura inicial estimada
    """
    estado = copy.deepcopy(estado_inicial)
    energia_actual = funcion_energia(estado)

    deltas_positivos = []

    for _ in range(n_muestras):
        estado_vecino, info = generar_vecino(estado)
        energia_vecino = funcion_energia(estado_vecino)
        delta = energia_vecino - energia_actual

        if delta > 0:
            deltas_positivos.append(delta)

        # Actualizar aleatoriamente para explorar
        if np.random.random() < 0.5:
            estado = estado_vecino
            energia_actual = energia_vecino

    if deltas_positivos:
        delta_promedio = np.mean(deltas_positivos)
        # De exp(-ΔE/T) = p, despejamos T = -ΔE/ln(p)
        T0 = -delta_promedio / np.log(probabilidad_aceptacion_inicial)
        return T0
    else:
        # Si no hay deltas positivos, usar un valor por defecto
        return 1.0


def recocido_simulado_red(
    posiciones: Dict[str, np.ndarray],
    tipos: Dict[str, str],
    num_ti: int,
    temp_inicial: float = 5.0,
    alpha: float = 0.95,
    iteraciones_temp: int = 50,
    max_iter: int = 10000,
    verbose: bool = False
) -> ResultadoSA:
    """
    Función simplificada de recocido simulado para redes de átomos.

    Args:
        posiciones: Diccionario {id_atomo: np.array([x, y, ...])}
        tipos: Diccionario {id_atomo: tipo_str} con tipos 'Fe', 'Nd', etc.
        num_ti: Número de átomos de Ti a colocar
        temp_inicial: Temperatura inicial
        alpha: Factor de enfriamiento exponencial
        iteraciones_temp: Iteraciones por nivel de temperatura
        max_iter: Máximo de iteraciones totales
        verbose: Mostrar progreso

    Returns:
        ResultadoSA con el mejor estado encontrado
    """
    from .morse_potential import calcular_energia_total

    # Crear estado inicial colocando Ti aleatoriamente
    estado = copy.deepcopy(tipos)
    atomos_fe = [id_a for id_a, t in estado.items() if t == 'Fe']

    # Colocar num_ti átomos de Ti aleatoriamente
    if num_ti > 0 and num_ti <= len(atomos_fe):
        indices_ti = np.random.choice(len(atomos_fe), size=num_ti, replace=False)
        for idx in indices_ti:
            estado[atomos_fe[idx]] = 'Ti'

    # Función de energía
    def funcion_energia(est):
        return calcular_energia_total(posiciones, est)

    # Función generadora de vecinos (swap Ti-Fe)
    def generar_vecino(est):
        nuevo_estado = copy.deepcopy(est)
        atomos_ti = [id_a for id_a, t in nuevo_estado.items() if t == 'Ti']
        atomos_fe = [id_a for id_a, t in nuevo_estado.items() if t == 'Fe']

        if atomos_ti and atomos_fe:
            ti_sel = np.random.choice(atomos_ti)
            fe_sel = np.random.choice(atomos_fe)
            nuevo_estado[ti_sel] = 'Fe'
            nuevo_estado[fe_sel] = 'Ti'
            return nuevo_estado, {'ti_swap': ti_sel, 'fe_swap': fe_sel}
        return nuevo_estado, {}

    # Esquema de enfriamiento
    esquema = esquema_enfriamiento_exponencial(alpha)

    # Ejecutar recocido simulado
    return recocido_simulado(
        estado_inicial=estado,
        funcion_energia=funcion_energia,
        generar_vecino=generar_vecino,
        temperatura_inicial=temp_inicial,
        esquema_enfriamiento=esquema,
        n_iteraciones_max=max_iter,
        n_iteraciones_por_T=iteraciones_temp,
        verbose=verbose
    )