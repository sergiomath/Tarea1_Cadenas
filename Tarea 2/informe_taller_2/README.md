# Informe Tarea 2 - Conteo Aproximado con MCMC

Informe profesional de la Tarea 2 del curso Cadenas de Markov y Aplicaciones.

## Contenido del Informe

**Tema**: Implementación de algoritmos de conteo aproximado basados en el Teorema 9.1 para:
- q-Coloraciones en grafos reticulares
- Modelo Hard-Core (conjuntos independientes)

## Archivos Principales

- **`main.tex`** - Documento principal LaTeX
- **`Informe_Tarea2_Final.pdf`** - PDF compilado final (730 KB)
- **`plantilla_src/`** - Sistema de plantilla LaTeX completo
- **`resultados/`** - Imágenes de resultados experimentales
- **`img/`** - Recursos gráficos (logos UNAL)

## Estructura del Documento

1. **Introducción**
   - Motivación y contexto
   - Objetivos del trabajo

2. **Marco Teórico**
   - Teorema 9.1 (FPRAS)
   - Muestreador de Gibbs
   - Definiciones matemáticas

3. **Implementación**
   - Arquitectura del software
   - Clases principales
   - Parámetros del algoritmo

4. **Resultados**
   - Experimentos con q-coloraciones
   - Análisis del modelo Hard-Core
   - Validación con conteo exacto

5. **Visualizaciones**
   - Gráficas de análisis completo
   - Figuras de resultados

6. **Conclusiones**
   - Logros principales
   - Observaciones clave
   - Trabajo futuro

## Compilación

### Método 1: Usando pdflatex
```bash
pdflatex main.tex
pdflatex main.tex  # Segunda vez para referencias
```

### Método 2: Usando latexmk (recomendado)
```bash
latexmk -pdf main.tex
```

## Autores

- Sergio Andrés Díaz Vera (seadiazve@unal.edu.co)
- Julián Mateo Espinosa Ospina (juespinosao@unal.edu.co)

## Profesor

Freddy Hernández-Romero

## Universidad

Universidad Nacional de Colombia
Facultad de Ciencias
Departamento de Matemáticas
2025-II

## Resultados Principales

- ✅ Implementación exitosa del algoritmo FPRAS
- ✅ Validación con valores exactos (error < 5%)
- ✅ Complejidad O(k⁴) confirmada experimentalmente
- ✅ Densidad Hard-Core: ρ ≈ 0.236

## Requisitos de Compilación

- LaTeX distribution (TeXLive o MikTeX)
- Paquetes: amsmath, amsthm, graphicx, hyperref, listings
- pdfLaTeX como compilador

## Notas

- El informe usa la plantilla institucional de la Universidad Nacional
- Las imágenes están en `resultados/`
- El código fuente está en `../src/`
- Los experimentos están en `../notebooks/`
