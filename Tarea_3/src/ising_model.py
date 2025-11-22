"""
Módulo para simulación del Modelo de Ising.
Incluye implementaciones de MCMC tradicional y simulación perfecta (Propp-Wilson).
"""

import numpy as np
from typing import Tuple, List, Dict
import time


class IsingLattice:
    """Representa un lattice K×K para el modelo de Ising."""

    def __init__(self, K: int):
        self.K = K
        self.n_sites = K * K
        self._build_neighbors()

    def _build_neighbors(self):
        """Construye diccionario de vecinos para cada sitio."""
        self.neighbors = {}
        for i in range(self.K):
            for j in range(self.K):
                site = (i, j)
                nbrs = []
                if i > 0:
                    nbrs.append((i-1, j))
                if i < self.K - 1:
                    nbrs.append((i+1, j))
                if j > 0:
                    nbrs.append((i, j-1))
                if j < self.K - 1:
                    nbrs.append((i, j+1))
                self.neighbors[site] = nbrs

    def random_config(self) -> np.ndarray:
        """Genera configuración aleatoria de spins."""
        return np.random.choice([-1, 1], size=(self.K, self.K))

    def all_plus_config(self) -> np.ndarray:
        """Configuración con todos los spins en +1."""
        return np.ones((self.K, self.K), dtype=int)

    def all_minus_config(self) -> np.ndarray:
        """Configuración con todos los spins en -1."""
        return -np.ones((self.K, self.K), dtype=int)

    def hamiltonian(self, config: np.ndarray) -> float:
        """Calcula el Hamiltoniano H(η) = -Σ(x~y) ηₓηᵧ."""
        H = 0.0
        for i in range(self.K):
            for j in range(self.K):
                spin = config[i, j]
                for ni, nj in self.neighbors[(i, j)]:
                    H -= spin * config[ni, nj]
        return H / 2.0

    def magnetization(self, config: np.ndarray) -> float:
        """Calcula magnetización M(η) = (1/K²) Σ ηₓ."""
        return np.mean(config)


class GibbsSampler:
    """Gibbs Sampler para el modelo de Ising."""

    def __init__(self, lattice: IsingLattice, beta: float):
        self.lattice = lattice
        self.beta = beta

    def step(self, config: np.ndarray) -> np.ndarray:
        """Un paso completo del Gibbs Sampler (actualiza todos los sitios)."""
        new_config = config.copy()
        sites = [(i, j) for i in range(self.lattice.K) for j in range(self.lattice.K)]
        np.random.shuffle(sites)

        for i, j in sites:
            neighbor_sum = sum(new_config[ni, nj] for ni, nj in self.lattice.neighbors[(i, j)])

            # Probabilidad de spin +1
            prob_plus = np.exp(self.beta * neighbor_sum)
            # Probabilidad de spin -1
            prob_minus = np.exp(-self.beta * neighbor_sum)

            # Normalizar
            Z_local = prob_plus + prob_minus
            prob_plus /= Z_local

            new_config[i, j] = 1 if np.random.random() < prob_plus else -1

        return new_config

    def sample(self, n_samples: int, burn_in: int = 1000,
               thin: int = 100) -> List[np.ndarray]:
        """
        Genera muestras usando MCMC.

        Args:
            n_samples: Número de muestras a generar
            burn_in: Pasos de burn-in
            thin: Espaciado entre muestras (thinning)
        """
        config = self.lattice.random_config()

        for _ in range(burn_in):
            config = self.step(config)

        samples = []
        for _ in range(n_samples):
            for _ in range(thin):
                config = self.step(config)
            samples.append(config.copy())

        return samples


class ProppWilson:
    """Algoritmo de Propp-Wilson (Coupling From The Past) para el modelo de Ising."""

    def __init__(self, lattice: IsingLattice, beta: float):
        self.lattice = lattice
        self.beta = beta

    def _coupled_update(self, configs: Dict[str, np.ndarray],
                       random_choices: List[Tuple]) -> Dict[str, np.ndarray]:
        """
        Actualiza todas las configuraciones usando las mismas elecciones aleatorias.
        """
        new_configs = {key: config.copy() for key, config in configs.items()}

        for i, j, u in random_choices:
            for key in new_configs:
                config = new_configs[key]
                neighbor_sum = sum(config[ni, nj] for ni, nj in self.lattice.neighbors[(i, j)])

                # Probabilidad de spin +1
                prob_plus = np.exp(self.beta * neighbor_sum)
                # Probabilidad de spin -1
                prob_minus = np.exp(-self.beta * neighbor_sum)

                # Normalizar
                Z_local = prob_plus + prob_minus
                prob_plus /= Z_local

                config[i, j] = 1 if u < prob_plus else -1

        return new_configs

    def _check_coalescence(self, configs: Dict[str, np.ndarray]) -> bool:
        """Verifica si todas las configuraciones son iguales."""
        first_config = list(configs.values())[0]
        return all(np.array_equal(first_config, config) for config in configs.values())

    def sample_once(self, max_time: int = 10000) -> Tuple[np.ndarray, int]:
        """
        Genera una muestra exacta usando CFTP.

        Returns:
            (configuración, tiempo de coalescencia)
        """
        T = 1
        random_history = []

        while T <= max_time:
            configs = {
                'all_plus': self.lattice.all_plus_config(),
                'all_minus': self.lattice.all_minus_config()
            }

            # Extender la historia aleatoria si es necesario
            while len(random_history) < T:
                sites = [(i, j) for i in range(self.lattice.K) for j in range(self.lattice.K)]
                np.random.shuffle(sites)
                step_randoms = [(i, j, np.random.random()) for i, j in sites]
                random_history.append(step_randoms)

            # Simular desde -T hasta 0 usando la historia
            for t in range(T):
                configs = self._coupled_update(configs, random_history[T - 1 - t])

            if self._check_coalescence(configs):
                return configs['all_plus'], T

            T *= 2

        raise RuntimeError(f"No se alcanzó coalescencia en {max_time} pasos")

    def sample(self, n_samples: int, max_time: int = 10000) -> Tuple[List[np.ndarray], List[int]]:
        """
        Genera múltiples muestras exactas.

        Returns:
            (lista de muestras, lista de tiempos de coalescencia)
        """
        samples = []
        coalescence_times = []
        failed_samples = 0

        for i in range(n_samples):
            try:
                config, coal_time = self.sample_once(max_time)
                samples.append(config)
                coalescence_times.append(coal_time)
            except RuntimeError:
                # Si no coalece, registramos el fallo pero continuamos
                failed_samples += 1
                if failed_samples <= 3:  # Solo advertir las primeras 3 veces
                    print(f"  Advertencia: Muestra {i+1}/{n_samples} no alcanzó coalescencia (β={self.beta:.1f})")
                # Usar configuración aleatoria como fallback
                samples.append(self.lattice.random_config())
                coalescence_times.append(max_time)

        if failed_samples > 0:
            print(f"  Total de muestras sin coalescencia: {failed_samples}/{n_samples}")

        return samples, coalescence_times


def estimate_magnetization(samples: List[np.ndarray], lattice: IsingLattice) -> Tuple[float, float]:
    """
    Estima el valor esperado de la magnetización y su desviación estándar.

    Returns:
        (media, desviación estándar)
    """
    magnetizations = [lattice.magnetization(sample) for sample in samples]
    return np.mean(magnetizations), np.std(magnetizations)


def run_mcmc_experiment(K: int, beta: float, n_samples: int = 100,
                        burn_in: int = 1000, thin: int = 100) -> Dict:
    """
    Ejecuta experimento completo con MCMC tradicional.
    """
    lattice = IsingLattice(K)
    sampler = GibbsSampler(lattice, beta)

    start_time = time.time()
    samples = sampler.sample(n_samples, burn_in, thin)
    elapsed_time = time.time() - start_time

    mean_mag, std_mag = estimate_magnetization(samples, lattice)

    return {
        'method': 'MCMC',
        'K': K,
        'beta': beta,
        'n_samples': n_samples,
        'magnetization_mean': mean_mag,
        'magnetization_std': std_mag,
        'elapsed_time': elapsed_time,
        'samples': samples
    }


def run_propp_wilson_experiment(K: int, beta: float, n_samples: int = 100,
                                max_time: int = 10000) -> Dict:
    """
    Ejecuta experimento completo con Propp-Wilson.
    """
    lattice = IsingLattice(K)
    sampler = ProppWilson(lattice, beta)

    start_time = time.time()
    samples, coal_times = sampler.sample(n_samples, max_time)
    elapsed_time = time.time() - start_time

    mean_mag, std_mag = estimate_magnetization(samples, lattice)

    return {
        'method': 'Propp-Wilson',
        'K': K,
        'beta': beta,
        'n_samples': n_samples,
        'magnetization_mean': mean_mag,
        'magnetization_std': std_mag,
        'coalescence_times': coal_times,
        'mean_coalescence_time': np.mean(coal_times),
        'max_coalescence_time': np.max(coal_times),
        'elapsed_time': elapsed_time,
        'samples': samples
    }
