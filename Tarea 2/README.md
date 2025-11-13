# Tarea 2: Conteo Aproximado con MCMC

Implementación de algoritmos de conteo aproximado basados en Cadenas de Markov Monte Carlo (MCMC) para modelos Hard-Core y q-coloraciones.

## Estructura del Proyecto

```
Tarea_2/
├── src/
│   ├── __init__.py           # Paquete Python
│   ├── mcmc_counting.py      # Implementaciones principales (913 líneas)
│   └── mcmc_improved.py      # Implementación alternativa (legacy)
├── notebooks/
│   ├── tarea2_conteo_aproximado.ipynb    # Implementación básica
│   └── tarea2_completa.ipynb             # Implementación mejorada (Teorema 9.1)
├── resultados/
│   ├── q_coloraciones_analisis.png
│   ├── hardcore_analisis.png
│   └── *.csv                 # Resultados de experimentos
├── informe_taller_2/         # Informe LaTeX
│   └── main.tex
├── docs/
│   └── CMA_conteo aproximado.pdf  # Enunciado original
└── test_notebooks.py         # Script de verificación
```

## Módulo Principal: `src/mcmc_counting.py`

### Clases Implementadas

#### Configuración
- `MCMCConfig`: Configuración de parámetros (epsilon, num_samples, etc.)

#### Grafos
- `LatticeGraph`: Retículas K×K con construcción de adyacencia

#### Q-Coloraciones
- `QColoringMCMC`: Implementación básica con Gibbs sampler
- `QColoringApproximation`: Implementación mejorada basada en Teorema 9.1

#### Hard-Core
- `HardCoreMCMC`: Implementación básica
- `HardCoreApproximation`: Implementación mejorada

### Funciones de Conteo Exacto
- `exact_q_colorings_small(k, q)`: Fuerza bruta para k ≤ 3
- `exact_chromatic_polynomial(k, q)`: Polinomio cromático para k ≤ 4
- `exact_hardcore_small(k)`: Fuerza bruta Hard-Core para k ≤ 4
- `exact_hardcore_count(k)`: Conteo exacto para k ≤ 5

### Funciones de Experimentos
- `run_q_coloring_experiment(K, q, epsilon, verbose)`
- `run_hardcore_experiment(K, epsilon, verbose)`

## Notebooks

### `tarea2_conteo_aproximado.ipynb`
Implementación básica con:
- Experimentos para q-coloraciones (3 ≤ K ≤ 20, 2 ≤ q ≤ 15)
- Experimentos para Hard-Core
- Comparación con valores exactos
- Visualizaciones completas

### `tarea2_completa.ipynb`
Implementación mejorada basada en Teorema 9.1:
- Cálculo de parámetros según teoría (m, τ)
- Análisis de escalamiento
- Comparación con polinomio cromático
- Reportes detallados de parámetros

## Instalación

```bash
# Desde la raíz del proyecto
cd Tareas/Tarea_2

# Activar entorno virtual
source ../../.venv/bin/activate  # Linux/macOS
../../.venv/Scripts/activate     # Windows

# Verificar instalación
python test_notebooks.py
```

## Uso

### Desde Python

```python
import sys
sys.path.insert(0, '.')

from src.mcmc_counting import (
    LatticeGraph,
    QColoringApproximation,
    run_q_coloring_experiment
)

# Crear lattice
lattice = LatticeGraph(K=5)

# Ejecutar experimento
result = run_q_coloring_experiment(K=5, q=10, epsilon=0.1, verbose=True)
print(f"Estimación: {result['estimate']}")
print(f"Simulaciones: {result['num_simulations']}")
print(f"Tiempo de mezcla: {result['mixing_time']}")
```

### Desde Notebooks

```bash
# Jupyter Notebook
jupyter notebook notebooks/

# JupyterLab
jupyter lab notebooks/

# O abrir directamente en VSCode
code notebooks/tarea2_completa.ipynb
```

## Parámetros del Teorema 9.1

Para grafos reticulares con grado máximo d=4:

### Número de Simulaciones
```
m = ⌈48d²k³/ε²⌉
```
Donde:
- k = K² (número de vértices)
- ε = precisión deseada

### Tiempo de Mezcla
```
τ = ⌈k((2log k + log(1/ε) + log 8) / log(q/(q-1)) + 1)⌉
```
Donde:
- q = número de colores
- Condición: q > 2d = 8

## Resultados Principales

### Q-Coloraciones
- Rango explorado: 3 ≤ K ≤ 20, 2 ≤ q ≤ 15
- Precisiones: ε ∈ {0.5, 0.2, 0.1, 0.05}
- Error típico con ε=0.1: < 5%

### Hard-Core
- Rango explorado: 3 ≤ K ≤ 20
- Densidad crítica observada: ρ∞ ≈ 0.236
- Escalamiento: Z_HC(L_K) ~ exp(αK²) con α ≈ 0.89

## Verificación

El script `test_notebooks.py` verifica:
1. ✓ Importaciones del módulo
2. ✓ Creación de objetos básicos
3. ✓ Experimentos de Q-Coloring
4. ✓ Experimentos de Hard-Core
5. ✓ Conteo exacto para validación

```bash
python test_notebooks.py
```

## Dependencias

- numpy
- pandas
- matplotlib
- seaborn
- tqdm
- networkx
- scipy

## Autores

- Sergio Andrés Díaz Vera (seadiazve@unal.edu.co)
- Julián Mateo Espinosa Ospina (juespinosao@unal.edu.co)

## Referencias

- Levin, D.A., Peres, Y. (2017). *Markov Chains and Mixing Times*. AMS.
- Teorema 9.1 sobre esquemas de aproximación polinomial aleatorizada (FPRAS)
