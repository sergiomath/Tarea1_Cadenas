#!/usr/bin/env python
"""
Script de verificación de instalación para Tarea 1
Ejecutar con: python verificar_instalacion.py
"""
import sys

print("=" * 60)
print("VERIFICACIÓN DE INSTALACIÓN - TAREA 1")
print("=" * 60)
print()

# Test 1: Importar librerías básicas
print("📦 Verificando librerías básicas...")
try:
    import numpy as np
    print(f"  ✅ NumPy {np.__version__}")
except ImportError as e:
    print(f"  ❌ NumPy no disponible: {e}")
    sys.exit(1)

try:
    import matplotlib
    print(f"  ✅ Matplotlib {matplotlib.__version__}")
except ImportError as e:
    print(f"  ❌ Matplotlib no disponible: {e}")
    sys.exit(1)

try:
    import IPython
    print(f"  ✅ IPython {IPython.__version__}")
except ImportError as e:
    print(f"  ❌ IPython no disponible: {e}")
    sys.exit(1)

print()

# Test 2: Importar módulos del proyecto
print("🔧 Verificando módulos del proyecto...")
sys.path.insert(0, '.')

try:
    from src.hard_core import (
        gibbs_sampler_hard_core,
        contar_particulas,
        es_configuracion_factible
    )
    print("  ✅ Módulo hard_core importado")
except ImportError as e:
    print(f"  ❌ Error importando hard_core: {e}")
    sys.exit(1)

try:
    from src.q_coloraciones import gibbs_sampler_q_coloraciones
    print("  ✅ Módulo q_coloraciones importado")
except ImportError as e:
    print(f"  ❌ Error importando q_coloraciones: {e}")
    sys.exit(1)

try:
    from src.visualizacion import visualizar_configuracion
    print("  ✅ Módulo visualizacion importado")
except ImportError as e:
    print(f"  ❌ Error importando visualizacion: {e}")
    sys.exit(1)

try:
    from src.estadisticas import calcular_estadisticas
    print("  ✅ Módulo estadisticas importado")
except ImportError as e:
    print(f"  ❌ Error importando estadisticas: {e}")
    sys.exit(1)

print()

# Test 3: Ejecutar funciones básicas
print("🧪 Verificando funcionalidad...")

try:
    config = gibbs_sampler_hard_core(K=5, T=100, semilla=42)
    n_part = contar_particulas(config)
    factible = es_configuracion_factible(config)
    print(f"  ✅ Hard-Core: Configuración 5×5 con {n_part} partículas")

    if not factible:
        print("  ⚠️  Advertencia: Configuración no es factible")

except Exception as e:
    print(f"  ❌ Error ejecutando Hard-Core: {e}")
    sys.exit(1)

try:
    # Grafo simple: triángulo
    edges = [(0, 1), (1, 2), (2, 0)]
    config_q = gibbs_sampler_q_coloraciones(n=3, edges=edges, q=3, T=100, semilla=42)
    print(f"  ✅ q-Coloraciones: Configuración de {len(config_q)} nodos")
except Exception as e:
    print(f"  ❌ Error ejecutando q-Coloraciones: {e}")
    sys.exit(1)

try:
    muestras = [10, 15, 12, 18, 14]
    stats = calcular_estadisticas(muestras)
    print(f"  ✅ Estadísticas: Media = {stats['media']:.2f}")
except Exception as e:
    print(f"  ❌ Error calculando estadísticas: {e}")
    sys.exit(1)

print()
print("=" * 60)
print("✅ INSTALACIÓN COMPLETA Y FUNCIONAL")
print("=" * 60)
print()
print("Puedes ejecutar los notebooks con:")
print("  • uv run jupyter notebook notebooks/")
print("  • jupyter notebook notebooks/  (con ambiente activado)")
print()
