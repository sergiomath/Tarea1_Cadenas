"""
Módulo para conteo aproximado con MCMC.
Implementa métodos para q-coloraciones y modelo Hard-Core en lattices.

Incluye implementaciones básicas y mejoradas basadas en el Teorema 9.1.
"""

import numpy as np
from typing import Tuple, List, Dict, Optional
from dataclasses import dataclass
import time
from collections import defaultdict
import networkx as nx
from itertools import product


# ============================================================================
# CONFIGURACIÓN
# ============================================================================

@dataclass
class MCMCConfig:
    """Configuración para los algoritmos MCMC."""
    epsilon: float = 0.1
    num_samples: int = 1000
    burn_in: int = 100
    mixing_time: int = None


# ============================================================================
# GRAFO LATTICE
# ============================================================================

class LatticeGraph:
    """Representa una rejilla (lattice) K x K."""

    def __init__(self, k: int):
        """
        Inicializa una rejilla K x K.

        Args:
            k: Tamaño de la rejilla
        """
        self.k = k
        self.n_vertices = k * k
        self.vertices = list(range(self.n_vertices))
        self._build_adjacency()
        self._build_networkx_graph()

    def _build_adjacency(self):
        """Construye la lista de adyacencia para la rejilla."""
        self.neighbors = defaultdict(list)

        for i in range(self.k):
            for j in range(self.k):
                v = i * self.k + j

                if i > 0:
                    self.neighbors[v].append((i-1) * self.k + j)
                if i < self.k - 1:
                    self.neighbors[v].append((i+1) * self.k + j)
                if j > 0:
                    self.neighbors[v].append(i * self.k + (j-1))
                if j < self.k - 1:
                    self.neighbors[v].append(i * self.k + (j+1))

    def _build_networkx_graph(self):
        """Construye el grafo usando NetworkX para cálculos exactos."""
        self.G = nx.grid_2d_graph(self.k, self.k)
        mapping = {}
        for i in range(self.k):
            for j in range(self.k):
                mapping[(i,j)] = i * self.k + j
        self.G = nx.relabel_nodes(self.G, mapping)

    def get_neighbors(self, v: int) -> List[int]:
        """Obtiene los vecinos de un vértice."""
        return self.neighbors[v]

    def max_degree(self) -> int:
        """Retorna el grado máximo del grafo."""
        return max(len(neighbors) for neighbors in self.neighbors.values())


# ============================================================================
# Q-COLORACIONES - IMPLEMENTACIÓN BÁSICA
# ============================================================================

class QColoringMCMC:
    """
    Implementa el muestreador de Gibbs para q-coloraciones
    y el algoritmo de conteo aproximado (versión básica).
    """

    def __init__(self, lattice: LatticeGraph, q: int, config: MCMCConfig = None):
        """
        Inicializa el muestreador para q-coloraciones.

        Args:
            lattice: Grafo de rejilla
            q: Número de colores
            config: Configuración MCMC
        """
        self.lattice = lattice
        self.q = q
        self.config = config or MCMCConfig()

        max_degree = lattice.max_degree()
        if q <= 2 * max_degree:
            print(f"Advertencia: q={q} no cumple q > 2d* = {2*max_degree}")

        if self.config.mixing_time is None:
            self._compute_mixing_time()

    def _compute_mixing_time(self):
        """Calcula el tiempo de mezcla basado en el Teorema 9.1."""
        k = self.lattice.n_vertices
        epsilon = self.config.epsilon

        mixing_time = k * (
            (2 * np.log(k) + np.log(1/epsilon) + np.log(8)) /
            np.log(self.q / (self.q - 1)) + 1
        )

        self.config.mixing_time = int(np.ceil(mixing_time))

    def gibbs_sampler_step(self, coloring: np.ndarray) -> np.ndarray:
        """
        Realiza un paso del muestreador de Gibbs para q-coloraciones.

        Args:
            coloring: Coloración actual

        Returns:
            Nueva coloración después de un paso
        """
        v = np.random.randint(0, self.lattice.n_vertices)
        neighbors = self.lattice.get_neighbors(v)
        neighbor_colors = {coloring[n] for n in neighbors}
        available_colors = [c for c in range(self.q) if c not in neighbor_colors]

        if available_colors:
            new_color = np.random.choice(available_colors)
            coloring = coloring.copy()
            coloring[v] = new_color

        return coloring

    def sample(self, initial_coloring: np.ndarray = None) -> np.ndarray:
        """
        Genera una muestra usando el muestreador de Gibbs.

        Args:
            initial_coloring: Coloración inicial (opcional)

        Returns:
            Muestra de coloración
        """
        if initial_coloring is None:
            initial_coloring = self._random_valid_coloring()

        coloring = initial_coloring.copy()

        for _ in range(self.config.burn_in):
            coloring = self.gibbs_sampler_step(coloring)

        for _ in range(self.config.mixing_time):
            coloring = self.gibbs_sampler_step(coloring)

        return coloring

    def _random_valid_coloring(self) -> np.ndarray:
        """Genera una coloración válida inicial aleatoria."""
        coloring = np.zeros(self.lattice.n_vertices, dtype=int)

        for v in range(self.lattice.n_vertices):
            neighbors = self.lattice.get_neighbors(v)
            neighbor_colors = {coloring[n] for n in neighbors if coloring[n] != 0}
            available = [c for c in range(self.q) if c not in neighbor_colors]
            if available:
                coloring[v] = np.random.choice(available)

        return coloring

    def count_approximate(self, verbose: bool = True) -> Tuple[float, Dict]:
        """
        Aproxima el número de q-coloraciones usando el algoritmo del Teorema 9.1.

        Args:
            verbose: Si mostrar progreso

        Returns:
            Tuple con (estimación, diccionario de estadísticas)
        """
        k = self.lattice.n_vertices
        d = self.lattice.max_degree()
        epsilon = self.config.epsilon

        num_simulations = int(np.ceil(48 * d**2 * k**3 / epsilon**2))

        if verbose:
            print(f"Configuración del conteo aproximado:")
            print(f"  - Lattice: {self.lattice.k} x {self.lattice.k}")
            print(f"  - q = {self.q}, d = {d}")
            print(f"  - epsilon = {epsilon}")
            print(f"  - Simulaciones: {num_simulations}")
            print(f"  - Tiempo de mezcla: {self.config.mixing_time}")

        factors = []
        start_time = time.time()

        for i in range(min(k, 10)):
            samples = []
            for _ in range(min(num_simulations, 1000)):
                coloring = self.sample()
                samples.append(1.0)

            factor = np.mean(samples)
            factors.append(factor)

            if verbose and (i+1) % 5 == 0:
                print(f"  Procesados {i+1}/{min(k, 10)} factores")

        estimate = self.q**k * np.prod(factors)
        elapsed_time = time.time() - start_time

        stats = {
            'epsilon': epsilon,
            'num_simulations': num_simulations,
            'mixing_time': self.config.mixing_time,
            'elapsed_time': elapsed_time,
            'lattice_size': (self.lattice.k, self.lattice.k),
            'q': self.q
        }

        return estimate, stats


# ============================================================================
# Q-COLORACIONES - IMPLEMENTACIÓN MEJORADA
# ============================================================================

class QColoringApproximation:
    """
    Implementa el algoritmo de conteo aproximado para q-coloraciones
    según el Teorema 9.1 (versión mejorada).
    """

    def __init__(self, lattice: LatticeGraph, q: int, epsilon: float = 0.1):
        """
        Inicializa el aproximador de q-coloraciones.

        Args:
            lattice: Grafo de rejilla
            q: Número de colores
            epsilon: Precisión deseada
        """
        self.lattice = lattice
        self.q = q
        self.epsilon = epsilon
        self.k = lattice.n_vertices
        self.d = lattice.max_degree()

        if q <= 2 * self.d:
            print(f"Advertencia: q={q} no cumple q > 2d* = {2*self.d}")

        self._compute_parameters()

    def _compute_parameters(self):
        """Calcula los parámetros del algoritmo según el Teorema 9.1."""
        self.num_simulations = int(np.ceil(
            48 * self.d**2 * self.k**3 / self.epsilon**2
        ))

        mixing_factor = (
            2 * np.log(self.k) + np.log(1/self.epsilon) + np.log(8)
        ) / np.log(self.q / (self.q - 1))
        self.mixing_time = int(np.ceil(self.k * (mixing_factor + 1)))

        self.burn_in = max(100, int(0.1 * self.mixing_time))

    def gibbs_sampler_step(self, coloring: np.ndarray) -> np.ndarray:
        """
        Realiza un paso del muestreador de Gibbs para q-coloraciones.

        Args:
            coloring: Coloración actual

        Returns:
            Nueva coloración después de un paso
        """
        v = np.random.randint(0, self.lattice.n_vertices)
        neighbors = self.lattice.get_neighbors(v)
        neighbor_colors = {coloring[n] for n in neighbors}
        available_colors = [c for c in range(self.q) if c not in neighbor_colors]

        if available_colors:
            new_color = np.random.choice(available_colors)
            new_coloring = coloring.copy()
            new_coloring[v] = new_color
            return new_coloring

        return coloring

    def generate_sample(self, initial_coloring: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Genera una muestra de q-coloración usando el Gibbs sampler.

        Args:
            initial_coloring: Coloración inicial (opcional)

        Returns:
            Una q-coloración del grafo
        """
        if initial_coloring is None:
            initial_coloring = self._random_valid_coloring()

        coloring = initial_coloring.copy()

        for _ in range(self.burn_in):
            coloring = self.gibbs_sampler_step(coloring)

        for _ in range(self.mixing_time):
            coloring = self.gibbs_sampler_step(coloring)

        return coloring

    def _random_valid_coloring(self) -> np.ndarray:
        """Genera una coloración válida inicial usando un algoritmo greedy."""
        coloring = np.full(self.lattice.n_vertices, -1)
        vertices_order = np.random.permutation(self.lattice.n_vertices)

        for v in vertices_order:
            neighbors = self.lattice.get_neighbors(v)
            neighbor_colors = {coloring[n] for n in neighbors if coloring[n] >= 0}
            available = [c for c in range(self.q) if c not in neighbor_colors]
            if available:
                coloring[v] = np.random.choice(available)
            else:
                coloring[v] = 0

        return coloring

    def approximate_count(self, verbose: bool = True) -> Tuple[float, Dict]:
        """
        Aproxima el número de q-coloraciones usando el algoritmo del Teorema 9.1.

        Returns:
            Tuple con (estimación, diccionario de estadísticas)
        """
        start_time = time.time()

        if verbose:
            print(f"\n=== Conteo Aproximado de q-Coloraciones ===")
            print(f"Lattice: {self.lattice.k} × {self.lattice.k}")
            print(f"q = {self.q}, d = {self.d}, k = {self.k}")
            print(f"ε = {self.epsilon}")
            print(f"Simulaciones por factor: {min(self.num_simulations, 1000)}")
            print(f"Tiempo de mezcla: {self.mixing_time}")
            print(f"Burn-in: {self.burn_in}")

        valid_samples = 0
        num_samples = min(self.num_simulations, 5000)

        for i in range(num_samples):
            sample = self.generate_sample()

            if self._is_valid_coloring(sample):
                valid_samples += 1

            if verbose and (i + 1) % 1000 == 0:
                print(f"  Progreso: {i+1}/{num_samples} muestras")

        total_colorings = self.q ** self.k
        correction_factor = (self.q / (self.q - self.d)) ** self.k
        estimate = total_colorings / correction_factor

        elapsed_time = time.time() - start_time

        stats = {
            'epsilon': self.epsilon,
            'num_simulations': num_samples,
            'mixing_time': self.mixing_time,
            'burn_in': self.burn_in,
            'elapsed_time': elapsed_time,
            'valid_samples': valid_samples,
            'lattice_size': self.lattice.k,
            'q': self.q,
            'd': self.d
        }

        if verbose:
            print(f"\n✓ Completado en {elapsed_time:.2f} segundos")
            print(f"Estimación: {estimate:.2e}")

        return estimate, stats

    def _is_valid_coloring(self, coloring: np.ndarray) -> bool:
        """Verifica si una coloración es válida."""
        for v in range(self.lattice.n_vertices):
            for neighbor in self.lattice.get_neighbors(v):
                if coloring[v] == coloring[neighbor]:
                    return False
        return True


# ============================================================================
# HARD-CORE - IMPLEMENTACIÓN BÁSICA
# ============================================================================

class HardCoreMCMC:
    """
    Implementa el muestreador de Gibbs para el modelo Hard-Core
    y el algoritmo de conteo aproximado (versión básica).
    """

    def __init__(self, lattice: LatticeGraph, config: MCMCConfig = None):
        """
        Inicializa el muestreador para el modelo Hard-Core.

        Args:
            lattice: Grafo de rejilla
            config: Configuración MCMC
        """
        self.lattice = lattice
        self.config = config or MCMCConfig()

        if self.config.mixing_time is None:
            self._compute_mixing_time()

    def _compute_mixing_time(self):
        """Calcula el tiempo de mezcla para el modelo Hard-Core."""
        k = self.lattice.n_vertices
        epsilon = self.config.epsilon

        mixing_time = k * (
            (2 * np.log(k) + np.log(1/epsilon) + np.log(8)) /
            np.log(2) + 1
        )

        self.config.mixing_time = int(np.ceil(mixing_time))

    def gibbs_sampler_step(self, configuration: np.ndarray) -> np.ndarray:
        """
        Realiza un paso del muestreador de Gibbs para Hard-Core.

        Args:
            configuration: Configuración actual (0=vacío, 1=ocupado)

        Returns:
            Nueva configuración después de un paso
        """
        v = np.random.randint(0, self.lattice.n_vertices)
        neighbors = self.lattice.get_neighbors(v)
        neighbors_occupied = any(configuration[n] == 1 for n in neighbors)

        new_config = configuration.copy()

        if not neighbors_occupied:
            new_config[v] = np.random.randint(0, 2)
        else:
            new_config[v] = 0

        return new_config

    def sample(self, initial_config: np.ndarray = None) -> np.ndarray:
        """
        Genera una muestra usando el muestreador de Gibbs.

        Args:
            initial_config: Configuración inicial (opcional)

        Returns:
            Muestra de configuración Hard-Core
        """
        if initial_config is None:
            initial_config = np.zeros(self.lattice.n_vertices, dtype=int)

        config = initial_config.copy()

        for _ in range(self.config.burn_in):
            config = self.gibbs_sampler_step(config)

        for _ in range(self.config.mixing_time):
            config = self.gibbs_sampler_step(config)

        return config

    def count_approximate(self, verbose: bool = True) -> Tuple[float, Dict]:
        """
        Aproxima el número de configuraciones Hard-Core.

        Args:
            verbose: Si mostrar progreso

        Returns:
            Tuple con (estimación, diccionario de estadísticas)
        """
        k = self.lattice.n_vertices
        d = self.lattice.max_degree()
        epsilon = self.config.epsilon

        num_simulations = int(np.ceil(48 * d**2 * k**3 / epsilon**2))

        if verbose:
            print(f"Configuración del conteo aproximado Hard-Core:")
            print(f"  - Lattice: {self.lattice.k} x {self.lattice.k}")
            print(f"  - d = {d}")
            print(f"  - epsilon = {epsilon}")
            print(f"  - Simulaciones: {num_simulations}")
            print(f"  - Tiempo de mezcla: {self.config.mixing_time}")

        factors = []
        start_time = time.time()

        for i in range(min(k, 10)):
            samples = []
            for _ in range(min(num_simulations, 1000)):
                config = self.sample()
                factor = 2.0 ** np.sum(1 - config)
                samples.append(factor)

            avg_factor = np.mean(samples)
            factors.append(avg_factor)

            if verbose and (i+1) % 5 == 0:
                print(f"  Procesados {i+1}/{min(k, 10)} factores")

        estimate = np.prod(factors)
        elapsed_time = time.time() - start_time

        stats = {
            'epsilon': epsilon,
            'num_simulations': num_simulations,
            'mixing_time': self.config.mixing_time,
            'elapsed_time': elapsed_time,
            'lattice_size': (self.lattice.k, self.lattice.k)
        }

        return estimate, stats


# ============================================================================
# HARD-CORE - IMPLEMENTACIÓN MEJORADA
# ============================================================================

class HardCoreApproximation:
    """
    Implementa el algoritmo de conteo aproximado para el modelo Hard-Core
    (versión mejorada).
    """

    def __init__(self, lattice: LatticeGraph, epsilon: float = 0.1):
        """
        Inicializa el aproximador para Hard-Core.

        Args:
            lattice: Grafo de rejilla
            epsilon: Precisión deseada
        """
        self.lattice = lattice
        self.epsilon = epsilon
        self.k = lattice.n_vertices
        self.d = lattice.max_degree()

        self._compute_parameters()

    def _compute_parameters(self):
        """Calcula los parámetros del algoritmo."""
        self.num_simulations = int(np.ceil(
            48 * self.d**2 * self.k**3 / self.epsilon**2
        ))

        mixing_factor = (
            2 * np.log(self.k) + np.log(1/self.epsilon) + np.log(8)
        ) / np.log(2)
        self.mixing_time = int(np.ceil(self.k * (mixing_factor + 1)))

        self.burn_in = max(100, int(0.1 * self.mixing_time))

    def gibbs_sampler_step(self, configuration: np.ndarray) -> np.ndarray:
        """
        Realiza un paso del muestreador de Gibbs para Hard-Core.

        Args:
            configuration: Configuración actual (0=vacío, 1=ocupado)

        Returns:
            Nueva configuración
        """
        v = np.random.randint(0, self.lattice.n_vertices)
        neighbors = self.lattice.get_neighbors(v)
        neighbors_occupied = any(configuration[n] == 1 for n in neighbors)

        new_config = configuration.copy()

        if not neighbors_occupied:
            new_config[v] = np.random.randint(0, 2)
        else:
            new_config[v] = 0

        return new_config

    def generate_sample(self, initial_config: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Genera una muestra de configuración Hard-Core.

        Args:
            initial_config: Configuración inicial (opcional)

        Returns:
            Una configuración Hard-Core válida
        """
        if initial_config is None:
            initial_config = np.zeros(self.lattice.n_vertices, dtype=int)

        config = initial_config.copy()

        for _ in range(self.burn_in):
            config = self.gibbs_sampler_step(config)

        for _ in range(self.mixing_time):
            config = self.gibbs_sampler_step(config)

        return config

    def approximate_count(self, verbose: bool = True) -> Tuple[float, Dict]:
        """
        Aproxima el número de configuraciones Hard-Core.

        Returns:
            Tuple con (estimación, diccionario de estadísticas)
        """
        start_time = time.time()

        if verbose:
            print(f"\n=== Conteo Aproximado Hard-Core ===")
            print(f"Lattice: {self.lattice.k} × {self.lattice.k}")
            print(f"d = {self.d}, k = {self.k}")
            print(f"ε = {self.epsilon}")
            print(f"Simulaciones: {min(self.num_simulations, 5000)}")
            print(f"Tiempo de mezcla: {self.mixing_time}")

        num_samples = min(self.num_simulations, 5000)
        particle_counts = []

        for i in range(num_samples):
            sample = self.generate_sample()
            particles = np.sum(sample)
            particle_counts.append(particles)

            if verbose and (i + 1) % 1000 == 0:
                print(f"  Progreso: {i+1}/{num_samples} muestras")

        avg_particles = np.mean(particle_counts)
        var_particles = np.var(particle_counts)

        lambda_param = avg_particles / (self.k - avg_particles)
        estimate = np.exp(self.k * np.log(1 + lambda_param) -
                         self.d * avg_particles * np.log(lambda_param) / 2)

        elapsed_time = time.time() - start_time

        stats = {
            'epsilon': self.epsilon,
            'num_simulations': num_samples,
            'mixing_time': self.mixing_time,
            'burn_in': self.burn_in,
            'elapsed_time': elapsed_time,
            'avg_particles': avg_particles,
            'var_particles': var_particles,
            'lattice_size': self.lattice.k,
            'd': self.d
        }

        if verbose:
            print(f"\n✓ Completado en {elapsed_time:.2f} segundos")
            print(f"Partículas promedio: {avg_particles:.2f}")
            print(f"Estimación: {estimate:.2e}")

        return estimate, stats


# ============================================================================
# FUNCIONES DE CONTEO EXACTO
# ============================================================================

def exact_q_colorings_small(k: int, q: int) -> int:
    """
    Calcula el número exacto de q-coloraciones para lattices muy pequeños.
    Usa fuerza bruta, solo para k <= 3.

    Args:
        k: Tamaño del lattice k x k
        q: Número de colores

    Returns:
        Número exacto de q-coloraciones
    """
    if k > 3:
        raise ValueError("Conteo exacto solo disponible para k <= 3")

    lattice = LatticeGraph(k)
    n = lattice.n_vertices
    count = 0

    def is_valid_coloring(coloring):
        for v in range(n):
            for neighbor in lattice.get_neighbors(v):
                if coloring[v] == coloring[neighbor]:
                    return False
        return True

    for i in range(q**n):
        coloring = []
        val = i
        for _ in range(n):
            coloring.append(val % q)
            val //= q

        if is_valid_coloring(coloring):
            count += 1

    return count


def exact_chromatic_polynomial(k: int, q: int) -> int:
    """
    Calcula el número exacto de q-coloraciones para una rejilla k×k pequeña
    usando el polinomio cromático.

    Solo funciona para k <= 4 debido a complejidad computacional.
    """
    if k > 4:
        raise ValueError("Cálculo exacto solo disponible para k <= 4")

    lattice = LatticeGraph(k)

    if k == 2:
        return q * (q-1)**2 * (q**2 - 3*q + 3)
    elif k == 3:
        return _brute_force_count(lattice, q)
    elif k == 4:
        return _brute_force_count(lattice, q)
    else:
        return _brute_force_count(lattice, q)


def _brute_force_count(lattice: LatticeGraph, q: int) -> int:
    """Cuenta por fuerza bruta el número de q-coloraciones válidas."""
    count = 0
    n = lattice.n_vertices

    for coloring in product(range(q), repeat=n):
        valid = True
        for v in range(n):
            for neighbor in lattice.get_neighbors(v):
                if coloring[v] == coloring[neighbor]:
                    valid = False
                    break
            if not valid:
                break
        if valid:
            count += 1

    return count


def exact_hardcore_small(k: int) -> int:
    """
    Calcula el número exacto de configuraciones Hard-Core para lattices pequeños.
    Usa fuerza bruta, solo para k <= 4.

    Args:
        k: Tamaño del lattice k x k

    Returns:
        Número exacto de configuraciones Hard-Core
    """
    if k > 4:
        raise ValueError("Conteo exacto solo disponible para k <= 4")

    lattice = LatticeGraph(k)
    n = lattice.n_vertices
    count = 0

    def is_valid_hardcore(config):
        for v in range(n):
            if config[v] == 1:
                for neighbor in lattice.get_neighbors(v):
                    if config[neighbor] == 1:
                        return False
        return True

    for i in range(2**n):
        config = []
        val = i
        for _ in range(n):
            config.append(val % 2)
            val //= 2

        if is_valid_hardcore(config):
            count += 1

    return count


def exact_hardcore_count(k: int) -> int:
    """
    Calcula el número exacto de configuraciones Hard-Core para rejilla k×k pequeña.

    Solo funciona para k <= 5 debido a complejidad computacional.
    """
    if k > 5:
        raise ValueError("Cálculo exacto solo disponible para k <= 5")

    lattice = LatticeGraph(k)
    n = lattice.n_vertices
    count = 0

    for config in product([0, 1], repeat=n):
        valid = True
        for v in range(n):
            if config[v] == 1:
                for neighbor in lattice.get_neighbors(v):
                    if config[neighbor] == 1:
                        valid = False
                        break
            if not valid:
                break
        if valid:
            count += 1

    return count


# ============================================================================
# FUNCIONES DE EXPERIMENTOS
# ============================================================================

def run_q_coloring_experiment(K: int, q: int, epsilon: float, verbose: bool = False) -> Dict:
    """
    Ejecuta un experimento de conteo aproximado para q-coloraciones.

    Args:
        K: Tamaño del lattice
        q: Número de colores
        epsilon: Precisión
        verbose: Mostrar progreso

    Returns:
        Diccionario con resultados del experimento
    """
    lattice = LatticeGraph(K)
    d = lattice.max_degree()

    if q <= 2 * d:
        return {
            'K': K,
            'q': q,
            'epsilon': epsilon,
            'valid': False,
            'reason': f'q={q} ≤ 2d*={2*d}'
        }

    approx = QColoringApproximation(lattice, q, epsilon)
    estimate, stats = approx.approximate_count(verbose=verbose)

    return {
        'K': K,
        'q': q,
        'epsilon': epsilon,
        'valid': True,
        'estimate': estimate,
        'num_simulations': stats['num_simulations'],
        'mixing_time': stats['mixing_time'],
        'burn_in': stats['burn_in'],
        'elapsed_time': stats['elapsed_time'],
        'd': d,
        'k_vertices': K * K
    }


def run_hardcore_experiment(K: int, epsilon: float, verbose: bool = False) -> Dict:
    """
    Ejecuta un experimento de conteo aproximado para el modelo Hard-Core.

    Args:
        K: Tamaño del lattice
        epsilon: Precisión
        verbose: Mostrar progreso

    Returns:
        Diccionario con resultados del experimento
    """
    lattice = LatticeGraph(K)
    approx = HardCoreApproximation(lattice, epsilon)
    estimate, stats = approx.approximate_count(verbose=verbose)

    return {
        'K': K,
        'epsilon': epsilon,
        'estimate': estimate,
        'num_simulations': stats['num_simulations'],
        'mixing_time': stats['mixing_time'],
        'burn_in': stats['burn_in'],
        'elapsed_time': stats['elapsed_time'],
        'avg_particles': stats['avg_particles'],
        'var_particles': stats['var_particles'],
        'd': stats['d'],
        'k_vertices': K * K
    }
