"""
Funciones para visualización de configuraciones y resultados
"""
import matplotlib.pyplot as plt
import numpy as np

def visualizar_configuracion(config, titulo="Configuración", colors=None):
    """Visualiza una configuración de la rejilla"""
    plt.figure(figsize=(8, 8))
    if colors is None:
        plt.imshow(config, cmap='binary', interpolation='nearest')
    else:
        plt.imshow(config, cmap=colors, interpolation='nearest', vmin=0)
    plt.title(titulo)
    plt.colorbar()
    plt.grid(True, alpha=0.3)
    return plt.gcf()

def graficar_histograma(datos, titulo="Histograma", xlabel="Valor"):
    """Crea un histograma de los datos"""
    plt.figure(figsize=(10, 6))
    plt.hist(datos, bins=30, edgecolor='black', alpha=0.7)
    plt.title(titulo)
    plt.xlabel(xlabel)
    plt.ylabel("Frecuencia")
    plt.grid(True, alpha=0.3)
    return plt.gcf()

def graficar_escalamiento(K_valores, medias):
    """Grafica el escalamiento del número de partículas vs K²"""
    plt.figure(figsize=(10, 6))
    K_cuadrado = [K**2 for K in K_valores]
    plt.scatter(K_cuadrado, medias, s=100, alpha=0.6)
    plt.plot(K_cuadrado, medias, 'b--', alpha=0.5)
    plt.xlabel("K² (Tamaño de la rejilla)")
    plt.ylabel("Número medio de partículas")
    plt.title("Escalamiento del modelo Hard-Core")
    plt.grid(True, alpha=0.3)
    return plt.gcf()

def graficar_distribucion_colores(conteos_colores, q):
    """Grafica la distribución de colores en una q-coloración"""
    plt.figure(figsize=(10, 6))
    colores = list(range(q))
    plt.bar(colores, conteos_colores, edgecolor='black', alpha=0.7)
    plt.xlabel("Color")
    plt.ylabel("Número de celdas")
    plt.title(f"Distribución de colores (q={q})")
    plt.xticks(colores)
    plt.grid(True, alpha=0.3, axis='y')
    return plt.gcf()
