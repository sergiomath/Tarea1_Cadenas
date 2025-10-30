# Guía de Instalación - Tarea 1: Gibbs Sampler

**Estudiantes:**
- Sergio Andrés Díaz Vera (seadiazve@unal.edu.co)
- Julián Mateo Espinosa Ospina (juespinosao@unal.edu.co)


## Opción 1: Instalación con UV (Recomendado)

UV es un manejador de paquetes moderno y rápido para Python.

### Paso 1: Instalar UV

**Linux/macOS:**
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

**Windows (PowerShell como Administrador):**
```powershell
powershell -c "irm https://astral.sh/uv/install.ps1 | iex"
```

**Alternativa con pip:**
```bash
pip install uv
```

**Verificar instalación:**
```bash
uv --version
```

### Paso 2: Navegar al directorio del proyecto

```bash
cd Tarea_1_Entregables
```

### Paso 3: Crear el ambiente e instalar dependencias

```bash
uv sync
```

Este comando:
- Creará automáticamente un entorno virtual en `.venv/`
- Instalará todas las dependencias especificadas en `pyproject.toml`
- Configurará el proyecto para uso inmediato

### Paso 4: Verificar instalación

```bash
uv run python -c "
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('.')
from src.hard_core import gibbs_sampler_hard_core, contar_particulas
print('✅ Todas las librerías instaladas correctamente')
"
```
