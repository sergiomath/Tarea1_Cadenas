"""
Funciones para crear y manipular redes cristalinas 2D.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple, Optional
import copy


def crear_red_2d_4x4() -> Tuple[Dict[str, np.ndarray], Dict[str, str]]:
    """
    Crea una red 2D de 4x4 con 4 átomos de Nd (R) y 12 de Fe.

    Returns:
        Tupla (posiciones, tipos) donde:
        - posiciones: {id_atomo: np.array([x, y])}
        - tipos: {id_atomo: tipo_str}
    """
    posiciones = {}
    tipos = {}

    # Átomos de Nd (R) en las esquinas interiores
    posiciones_nd = [(1, 1), (1, 2), (2, 1), (2, 2)]
    for i, (x, y) in enumerate(posiciones_nd):
        id_atomo = f"Nd_{i}"
        posiciones[id_atomo] = np.array([x, y], dtype=float)
        tipos[id_atomo] = 'Nd'

    # Átomos de Fe en las posiciones restantes
    contador_fe = 0
    for x in range(4):
        for y in range(4):
            if (x, y) not in posiciones_nd:
                id_atomo = f"Fe_{contador_fe}"
                posiciones[id_atomo] = np.array([x, y], dtype=float)
                tipos[id_atomo] = 'Fe'
                contador_fe += 1

    return posiciones, tipos


def crear_red_2d_10x10() -> Tuple[Dict[str, np.ndarray], Dict[str, str]]:
    """
    Crea una red 2D de 10x10 con 16 átomos de Nd en el centro y 84 de Fe.

    Returns:
        Tupla (posiciones, tipos)
    """
    posiciones = {}
    tipos = {}

    # Átomos de Nd en el cuadrado central (3,3) a (6,6)
    contador_nd = 0
    for x in range(3, 7):
        for y in range(3, 7):
            id_atomo = f"Nd_{contador_nd}"
            posiciones[id_atomo] = np.array([x, y], dtype=float)
            tipos[id_atomo] = 'Nd'
            contador_nd += 1

    # Átomos de Fe en las posiciones restantes
    contador_fe = 0
    for x in range(10):
        for y in range(10):
            if not (3 <= x <= 6 and 3 <= y <= 6):
                id_atomo = f"Fe_{contador_fe}"
                posiciones[id_atomo] = np.array([x, y], dtype=float)
                tipos[id_atomo] = 'Fe'
                contador_fe += 1

    return posiciones, tipos


def obtener_atomos_fe(tipos: Dict[str, str]) -> List[str]:
    """Obtiene lista de IDs de átomos de Fe."""
    return [id_atomo for id_atomo, tipo in tipos.items() if tipo == 'Fe']


def obtener_atomos_ti(tipos: Dict[str, str]) -> List[str]:
    """Obtiene lista de IDs de átomos de Ti."""
    return [id_atomo for id_atomo, tipo in tipos.items() if tipo == 'Ti']


def sustituir_fe_por_ti(
    posiciones: Dict[str, np.ndarray],
    tipos: Dict[str, str],
    indice_fe: int
) -> Tuple[Dict[str, np.ndarray], Dict[str, str]]:
    """
    Sustituye un átomo de Fe por Ti en la posición indicada.

    Args:
        posiciones: Diccionario de posiciones
        tipos: Diccionario de tipos
        indice_fe: Índice del átomo de Fe a sustituir (0-based)

    Returns:
        Nuevas copias de (posiciones, tipos) con la sustitución
    """
    posiciones_copia = copy.deepcopy(posiciones)
    tipos_copia = copy.deepcopy(tipos)

    atomos_fe = obtener_atomos_fe(tipos)
    if indice_fe >= len(atomos_fe):
        raise ValueError(f"Índice {indice_fe} fuera de rango. Hay {len(atomos_fe)} átomos de Fe.")

    id_fe = atomos_fe[indice_fe]
    # Cambiar tipo de Fe a Ti
    tipos_copia[id_fe] = 'Ti'

    return posiciones_copia, tipos_copia


def colocar_multiples_ti(
    posiciones: Dict[str, np.ndarray],
    tipos: Dict[str, str],
    indices_ti: List[int]
) -> Tuple[Dict[str, np.ndarray], Dict[str, str]]:
    """
    Coloca múltiples átomos de Ti en las posiciones indicadas.

    Args:
        posiciones: Diccionario de posiciones
        tipos: Diccionario de tipos
        indices_ti: Lista de índices de Fe a sustituir por Ti

    Returns:
        Nuevas copias de (posiciones, tipos) con las sustituciones
    """
    posiciones_copia = copy.deepcopy(posiciones)
    tipos_copia = copy.deepcopy(tipos)

    atomos_fe = obtener_atomos_fe(tipos)

    for idx in indices_ti:
        if idx >= len(atomos_fe):
            raise ValueError(f"Índice {idx} fuera de rango")
        id_fe = atomos_fe[idx]
        tipos_copia[id_fe] = 'Ti'

    return posiciones_copia, tipos_copia


def generar_vecino_swap(
    posiciones: Dict[str, np.ndarray],
    tipos: Dict[str, str]
) -> Tuple[Tuple[Dict, Dict], Dict]:
    """
    Genera un vecino intercambiando un átomo de Ti con uno de Fe.

    Args:
        posiciones: Diccionario de posiciones
        tipos: Diccionario de tipos

    Returns:
        Tupla ((nuevas_posiciones, nuevos_tipos), info_movimiento)
    """
    atomos_ti = obtener_atomos_ti(tipos)
    atomos_fe = obtener_atomos_fe(tipos)

    if not atomos_ti or not atomos_fe:
        # No hay átomos para intercambiar
        return (posiciones, tipos), {}

    # Seleccionar aleatoriamente
    ti_seleccionado = np.random.choice(atomos_ti)
    fe_seleccionado = np.random.choice(atomos_fe)

    # Hacer copia y realizar swap
    nuevos_tipos = copy.deepcopy(tipos)
    nuevos_tipos[ti_seleccionado] = 'Fe'
    nuevos_tipos[fe_seleccionado] = 'Ti'

    info = {
        'ti_swap': ti_seleccionado,
        'fe_swap': fe_seleccionado
    }

    return (posiciones, nuevos_tipos), info


def visualizar_red_2d(
    posiciones: Dict[str, np.ndarray],
    tipos: Dict[str, str],
    titulo: str = "Red 2D",
    guardar: Optional[str] = None,
    mostrar_etiquetas: bool = False
) -> None:
    """
    Visualiza la red 2D con diferentes colores para cada tipo de átomo.

    Args:
        posiciones: Diccionario de posiciones
        tipos: Diccionario de tipos
        titulo: Título del gráfico
        guardar: Path para guardar la figura (opcional)
        mostrar_etiquetas: Si mostrar etiquetas de átomos
    """
    fig, ax = plt.subplots(figsize=(8, 8))

    # Colores para cada tipo
    colores = {'Fe': 'blue', 'Nd': 'red', 'Ti': 'green'}
    tamaños = {'Fe': 200, 'Nd': 300, 'Ti': 250}

    # Agrupar por tipo
    for tipo_atomo in ['Nd', 'Fe', 'Ti']:
        x_coords = []
        y_coords = []
        for id_atomo, tipo in tipos.items():
            if tipo == tipo_atomo:
                pos = posiciones[id_atomo]
                x_coords.append(pos[0])
                y_coords.append(pos[1])

        if x_coords:
            ax.scatter(x_coords, y_coords,
                      c=colores[tipo_atomo],
                      s=tamaños[tipo_atomo],
                      label=f'{tipo_atomo} ({len(x_coords)})',
                      edgecolors='black',
                      linewidths=1.5,
                      alpha=0.8)

            if mostrar_etiquetas:
                for i, (x, y) in enumerate(zip(x_coords, y_coords)):
                    ax.annotate(f'{tipo_atomo[0]}', (x, y),
                              ha='center', va='center',
                              fontsize=10, fontweight='bold')

    ax.set_xlabel('X', fontsize=12)
    ax.set_ylabel('Y', fontsize=12)
    ax.set_title(titulo, fontsize=14, fontweight='bold')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')

    # Ajustar límites
    todas_x = [p[0] for p in posiciones.values()]
    todas_y = [p[1] for p in posiciones.values()]
    margin = 0.5
    ax.set_xlim(min(todas_x) - margin, max(todas_x) + margin)
    ax.set_ylim(min(todas_y) - margin, max(todas_y) + margin)

    if guardar:
        plt.savefig(guardar, dpi=150, bbox_inches='tight')
        print(f"Figura guardada en: {guardar}")

    plt.show()


def calcular_distancia_promedio_ti_nd(
    posiciones: Dict[str, np.ndarray],
    tipos: Dict[str, str]
) -> Tuple[float, List[float]]:
    """
    Calcula la distancia promedio entre cada átomo de Ti y su vecino Nd más cercano.

    Returns:
        Tupla (promedio, lista_distancias_minimas)
    """
    atomos_ti = [id_a for id_a, t in tipos.items() if t == 'Ti']
    atomos_nd = [id_a for id_a, t in tipos.items() if t == 'Nd']

    if not atomos_ti or not atomos_nd:
        return 0.0, []

    distancias_minimas = []

    for ti_id in atomos_ti:
        pos_ti = posiciones[ti_id]
        dist_min = float('inf')

        for nd_id in atomos_nd:
            pos_nd = posiciones[nd_id]
            dist = np.linalg.norm(pos_ti - pos_nd)
            if dist < dist_min:
                dist_min = dist

        distancias_minimas.append(dist_min)

    return np.mean(distancias_minimas), distancias_minimas


def calcular_distancia_promedio_ti_ti(
    posiciones: Dict[str, np.ndarray],
    tipos: Dict[str, str]
) -> Tuple[float, List[float]]:
    """
    Calcula la distancia promedio entre cada átomo de Ti y su vecino Ti más cercano.

    Returns:
        Tupla (promedio, lista_distancias_minimas)
    """
    atomos_ti = [id_a for id_a, t in tipos.items() if t == 'Ti']

    if len(atomos_ti) < 2:
        return 0.0, []

    distancias_minimas = []

    for i, ti_id in enumerate(atomos_ti):
        pos_ti = posiciones[ti_id]
        dist_min = float('inf')

        for j, otro_ti_id in enumerate(atomos_ti):
            if i == j:
                continue
            pos_otro = posiciones[otro_ti_id]
            dist = np.linalg.norm(pos_ti - pos_otro)
            if dist < dist_min:
                dist_min = dist

        distancias_minimas.append(dist_min)

    return np.mean(distancias_minimas), distancias_minimas