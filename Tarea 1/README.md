# Tarea 1: Muestreo de Gibbs - Entregables

**Estudiantes:**
- Sergio Andrés Díaz Vera (seadiazve@unal.edu.co)
- Julián Mateo Espinosa Ospina (juespinosao@unal.edu.co)

**Curso:** Cadenas de Markov y Aplicaciones (2025-II)
**Profesor:** Freddy Hernández-Romero
**Departamento:** Matemáticas - Universidad Nacional de Colombia

---
## 🚀 Instalación Rápida

### Prerequisitos
- **Python 3.8 o superior**
- **Git** (para clonar el repositorio)

### Opción 1: Con UV (Recomendado)

```bash
# 1. Instalar UV (si no lo tienes)
curl -LsSf https://astral.sh/uv/install.sh | sh  # Linux/macOS
# O con pip: pip install uv

# 2. Navegar al directorio
cd Tarea_1_Entregables

# 3. Instalar dependencias automáticamente
uv sync

# 4. Verificar instalación
uv run python verificar_instalacion.py

# 5. Ejecutar notebooks
uv run jupyter notebook notebooks/
```

### Opción 2: Con pip y venv

```bash
# 1. Navegar al directorio
cd Tarea_1_Entregables

# 2. Crear ambiente virtual
python3 -m venv .venv
source .venv/bin/activate  # Linux/macOS
# O en Windows: .venv\Scripts\activate

# 3. Instalar dependencias
pip install numpy>=1.21.0 matplotlib>=3.4.0 jupyter>=1.0.0 ipykernel>=6.0.0

# 4. Verificar instalación
python verificar_instalacion.py

# 5. Ejecutar notebooks
jupyter notebook notebooks/
```

### Opción 3: Con Conda

```bash
conda create -n tarea1 python=3.8 numpy matplotlib jupyter ipykernel -y
conda activate tarea1
jupyter notebook notebooks/
```

**📖 Para instrucciones detalladas:** Ver [INSTALACION.md](INSTALACION.md)

---

## 📦 Librerías Necesarias

| Librería | Versión | Propósito |
|----------|---------|-----------|
| **numpy** | ≥1.21.0 | Cálculos numéricos y matrices |
| **matplotlib** | ≥3.4.0 | Visualización de configuraciones |
| **jupyter** | ≥1.0.0 | Ejecución de notebooks |
| **ipykernel** | ≥6.0.0 | Kernel de Python para Jupyter |

---

## 🧪 Verificar Instalación

```bash
# Con UV
uv run python verificar_instalacion.py

# Con venv activado
python verificar_instalacion.py
```

**Salida esperada:** Mensaje confirmando que todas las librerías y módulos funcionan correctamente.

---

## 📂 Estructura del Proyecto

```
Tarea_1_Entregables/
├── notebooks/                 # Notebooks de análisis
│   ├── ejercicio_1a.ipynb    # Modelo Hard-Core
│   ├── ejercicio_1b.ipynb    # q-Coloraciones
│   └── ejercicio_2.ipynb     # Análisis comparativo
├── src/                       # Código fuente
│   ├── hard_core.py          # Gibbs Sampler Hard-Core
│   ├── q_coloraciones.py     # Gibbs Sampler q-Coloraciones
│   ├── visualizacion.py      # Funciones de gráficos
│   └── estadisticas.py       # Análisis estadístico
├── pyproject.toml            # Configuración del proyecto (UV)
├── verificar_instalacion.py  # Script de verificación
├── INSTALACION.md            # Guía detallada
└── README.md                 # Este archivo
```

### Orden de Ejecución Recomendado

1. **Ejercicio 1a** (`ejercicio_1a.ipynb`)
   - Implementación del Gibbs Sampler
   - Visualización de configuraciones
   - Trayectorias de la cadena

2. **Ejercicio 1b** (`ejercicio_1b.ipynb`)
   - Estimación del número de partículas
   - Histogramas de frecuencias
   - Verificación en diferentes tiempos

3. **Ejercicio 2** (`ejercicio_2.ipynb`)
   - Generalización a q-coloraciones
   - Distribución de colores
   - Comparación entre diferentes valores de q

---

## 📊 Resultados Principales

### Ejercicio 1: Modelo Hard-Core
- Rejillas de tamaño K×K (3 ≤ K ≤ 20)
- Tiempo de convergencia: 10,000 - 100,000 iteraciones
- Distribución uniforme sobre configuraciones factibles
- Número típico de partículas estimado

### Ejercicio 2: q-Coloraciones
- Número de colores: 2 ≤ q ≤ 10
- Coloraciones propias de la rejilla
- Distribución uniforme entre colores
- Análisis de escalamiento

---

## 🔍 Verificación de Funcionamiento

Para verificar que todo funciona correctamente:

```bash
# 1. Instalar dependencias (si no están instaladas)
pip install -r requirements.txt

# 2. Ejecutar script de verificación automática
python verificar_entregables.py

# 3. Si todo está OK, ejecutar notebooks
jupyter notebook notebooks/
```

