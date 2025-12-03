"""
Funciones para crear y manipular redes cristalinas 3D.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from typing import Dict, List, Tuple, Optional
import copy


def cargar_posiciones_3d(
    archivo_nd: str,
    archivo_fe: str
) -> Tuple[Dict[str, np.ndarray], Dict[str, str]]:
    """
    Carga las posiciones 3D desde archivos de texto.

    Args:
        archivo_nd: Path al archivo con posiciones de Nd
        archivo_fe: Path al archivo con posiciones de Fe

    Returns:
        Tupla (posiciones, tipos)
    """
    posiciones = {}
    tipos = {}

    # Cargar posiciones de Nd
    try:
        datos_nd = np.loadtxt(archivo_nd)
        if datos_nd.ndim == 1:
            datos_nd = datos_nd.reshape(1, -1)
        for i, pos in enumerate(datos_nd):
            id_atomo = f"Nd_{i}"
            posiciones[id_atomo] = pos
            tipos[id_atomo] = 'Nd'
    except Exception as e:
        print(f"Error cargando {archivo_nd}: {e}")

    # Cargar posiciones de Fe
    try:
        datos_fe = np.loadtxt(archivo_fe)
        if datos_fe.ndim == 1:
            datos_fe = datos_fe.reshape(1, -1)
        for i, pos in enumerate(datos_fe):
            id_atomo = f"Fe_{i}"
            posiciones[id_atomo] = pos
            tipos[id_atomo] = 'Fe'
    except Exception as e:
        print(f"Error cargando {archivo_fe}: {e}")

    return posiciones, tipos


def generar_supercelda_ndfe12() -> Tuple[Dict[str, np.ndarray], Dict[str, str]]:
    """
    Genera una supercelda 2x2x1 de NdFe12 con posiciones sintéticas.
    Se usa cuando no están disponibles los archivos de datos reales.

    Returns:
        Tupla (posiciones, tipos)
    """
    posiciones = {}
    tipos = {}

    # Parámetros de red (aproximados para NdFe12)
    a, b, c = 8.5, 8.5, 4.8  # Angstroms

    # Posiciones de Nd en la celda unitaria (4 átomos)
    nd_unitaria = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5]
    ]

    # Generar supercelda 2x2x1 para Nd (16 átomos)
    contador_nd = 0
    for i in range(2):
        for j in range(2):
            for pos in nd_unitaria:
                id_atomo = f"Nd_{contador_nd}"
                x = pos[0] * a + i * a
                y = pos[1] * b + j * b
                z = pos[2] * c
                posiciones[id_atomo] = np.array([x, y, z])
                tipos[id_atomo] = 'Nd'
                contador_nd += 1

    # Posiciones de Fe (96 átomos en supercelda 2x2x1)
    # Generamos posiciones en sitios específicos del grupo espacial
    contador_fe = 0

    # Sitios 8i, 8j, 8f para Fe en la estructura tetragonal
    for i in range(2):
        for j in range(2):
            # Generar 24 átomos de Fe por celda unitaria
            for k in range(24):
                id_atomo = f"Fe_{contador_fe}"
                # Distribuir Fe en posiciones características
                theta = 2 * np.pi * k / 24
                r = 0.3 + 0.2 * (k % 3) / 3
                x = (0.5 + r * np.cos(theta)) * a + i * a
                y = (0.5 + r * np.sin(theta)) * b + j * b
                z = (0.25 + 0.5 * (k % 2)) * c
                posiciones[id_atomo] = np.array([x, y, z])
                tipos[id_atomo] = 'Fe'
                contador_fe += 1

    return posiciones, tipos


def visualizar_red_3d(
    posiciones: Dict[str, np.ndarray],
    tipos: Dict[str, str],
    titulo: str = "Red 3D",
    guardar: Optional[str] = None,
    vista: Tuple[float, float] = (30, 45)
) -> None:
    """
    Visualiza la red 3D con diferentes colores para cada tipo de átomo.

    Args:
        posiciones: Diccionario de posiciones
        tipos: Diccionario de tipos
        titulo: Título del gráfico
        guardar: Path para guardar la figura
        vista: Ángulos de elevación y azimut para la vista (elev, azim)
    """
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Colores y tamaños
    colores = {'Fe': 'blue', 'Nd': 'red', 'Ti': 'green'}
    tamaños = {'Fe': 50, 'Nd': 100, 'Ti': 75}

    # Graficar por tipo
    for tipo_atomo in ['Nd', 'Fe', 'Ti']:
        coords = []
        for id_atomo, tipo in tipos.items():
            if tipo == tipo_atomo:
                coords.append(posiciones[id_atomo])

        if coords:
            coords = np.array(coords)
            ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2],
                      c=colores[tipo_atomo],
                      s=tamaños[tipo_atomo],
                      label=f'{tipo_atomo} ({len(coords)})',
                      alpha=0.7,
                      edgecolors='black',
                      linewidths=0.5)

    ax.set_xlabel('X (Å)', fontsize=10)
    ax.set_ylabel('Y (Å)', fontsize=10)
    ax.set_zlabel('Z (Å)', fontsize=10)
    ax.set_title(titulo, fontsize=14, fontweight='bold')
    ax.legend(loc='upper right')
    ax.view_init(elev=vista[0], azim=vista[1])

    # Ajustar límites
    todas_coords = np.array(list(posiciones.values()))
    for i, axis in enumerate(['x', 'y', 'z']):
        min_val = todas_coords[:, i].min()
        max_val = todas_coords[:, i].max()
        margin = (max_val - min_val) * 0.1
        if axis == 'x':
            ax.set_xlim(min_val - margin, max_val + margin)
        elif axis == 'y':
            ax.set_ylim(min_val - margin, max_val + margin)
        else:
            ax.set_zlim(min_val - margin, max_val + margin)

    if guardar:
        plt.savefig(guardar, dpi=150, bbox_inches='tight')
        print(f"Figura 3D guardada en: {guardar}")

    plt.show()


def proyeccion_2d(
    posiciones: Dict[str, np.ndarray],
    tipos: Dict[str, str],
    plano: str = 'xy',
    titulo: str = "Proyección 2D",
    guardar: Optional[str] = None
) -> None:
    """
    Proyecta la estructura 3D en un plano 2D.

    Args:
        posiciones: Diccionario de posiciones
        tipos: Diccionario de tipos
        plano: Plano de proyección ('xy', 'xz', 'yz')
        titulo: Título del gráfico
        guardar: Path para guardar
    """
    fig, ax = plt.subplots(figsize=(10, 10))

    # Índices según el plano
    if plano == 'xy':
        idx1, idx2 = 0, 1
        xlabel, ylabel = 'X (Å)', 'Y (Å)'
    elif plano == 'xz':
        idx1, idx2 = 0, 2
        xlabel, ylabel = 'X (Å)', 'Z (Å)'
    elif plano == 'yz':
        idx1, idx2 = 1, 2
        xlabel, ylabel = 'Y (Å)', 'Z (Å)'
    else:
        raise ValueError(f"Plano inválido: {plano}")

    # Colores y tamaños
    colores = {'Fe': 'blue', 'Nd': 'red', 'Ti': 'green'}
    tamaños = {'Fe': 100, 'Nd': 200, 'Ti': 150}

    # Graficar por tipo
    for tipo_atomo in ['Nd', 'Fe', 'Ti']:
        coords = []
        for id_atomo, tipo in tipos.items():
            if tipo == tipo_atomo:
                pos = posiciones[id_atomo]
                coords.append([pos[idx1], pos[idx2]])

        if coords:
            coords = np.array(coords)
            ax.scatter(coords[:, 0], coords[:, 1],
                      c=colores[tipo_atomo],
                      s=tamaños[tipo_atomo],
                      label=f'{tipo_atomo} ({len(coords)})',
                      alpha=0.7,
                      edgecolors='black',
                      linewidths=1)

    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(f"{titulo} - Plano {plano.upper()}", fontsize=14, fontweight='bold')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')

    if guardar:
        plt.savefig(guardar, dpi=150, bbox_inches='tight')
        print(f"Proyección guardada en: {guardar}")

    plt.show()


def analizar_distribucion_ti(
    posiciones: Dict[str, np.ndarray],
    tipos: Dict[str, str]
) -> Dict[str, float]:
    """
    Analiza la distribución de átomos de Ti en la red.

    Returns:
        Diccionario con métricas de análisis
    """
    metricas = {}

    atomos_ti = [id_a for id_a, t in tipos.items() if t == 'Ti']
    atomos_nd = [id_a for id_a, t in tipos.items() if t == 'Nd']
    atomos_fe = [id_a for id_a, t in tipos.items() if t == 'Fe']

    if not atomos_ti:
        return metricas

    # Distancia promedio Ti-Nd (al más cercano)
    dist_ti_nd = []
    for ti_id in atomos_ti:
        pos_ti = posiciones[ti_id]
        if atomos_nd:
            distancias = [np.linalg.norm(pos_ti - posiciones[nd_id])
                         for nd_id in atomos_nd]
            dist_ti_nd.append(min(distancias))

    if dist_ti_nd:
        metricas['dist_promedio_ti_nd_cercano'] = np.mean(dist_ti_nd)
        metricas['dist_min_ti_nd'] = min(dist_ti_nd)
        metricas['dist_max_ti_nd'] = max(dist_ti_nd)

    # Distancia promedio Ti-Ti (al más cercano)
    if len(atomos_ti) > 1:
        dist_ti_ti = []
        for i, ti_id in enumerate(atomos_ti):
            pos_ti = posiciones[ti_id]
            distancias = []
            for j, otro_ti in enumerate(atomos_ti):
                if i != j:
                    distancias.append(np.linalg.norm(pos_ti - posiciones[otro_ti]))
            if distancias:
                dist_ti_ti.append(min(distancias))

        if dist_ti_ti:
            metricas['dist_promedio_ti_ti_cercano'] = np.mean(dist_ti_ti)
            metricas['dist_min_ti_ti'] = min(dist_ti_ti)
            metricas['dist_max_ti_ti'] = max(dist_ti_ti)

    # Número de coordinación promedio de Ti (vecinos cercanos)
    radio_corte = 3.5  # Ångstroms
    coordinacion_ti = []
    for ti_id in atomos_ti:
        pos_ti = posiciones[ti_id]
        vecinos = 0
        for id_atomo, pos in posiciones.items():
            if id_atomo != ti_id:
                if np.linalg.norm(pos_ti - pos) < radio_corte:
                    vecinos += 1
        coordinacion_ti.append(vecinos)

    if coordinacion_ti:
        metricas['coordinacion_promedio_ti'] = np.mean(coordinacion_ti)

    return metricas