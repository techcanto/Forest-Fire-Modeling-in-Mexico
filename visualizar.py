import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# --- Gráfico de resumen temporal ---
archivo_resumen = "resumen_incendio.csv"
datos = pd.read_csv(archivo_resumen)

fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(datos["tiempo_s"], datos["frios"], label="Fríos", color="blue")
ax.plot(datos["tiempo_s"], datos["calentando"], label="Calentándose", color="orange")
ax.plot(datos["tiempo_s"], datos["ardiendo"], label="Ardiendo", color="red")
ax.set_xlabel("Tiempo (s)")
ax.set_ylabel("Número de árboles")
ax.set_title("Evolución del estado del bosque")
ax.grid(True)
ax.legend()
fig.tight_layout()
plt.savefig("evolucion_estados.png", dpi=150)
plt.show()

# --- Mapas de calor de temperatura ---
archivo_datos = "resultados_incendio.csv"
datos_incendio = pd.read_csv(archivo_datos)

# Buscamos los tiempos que realmente existen en el CSV (múltiplos de 20)
# Ajustamos los instantes para que coincidan con los datos guardados
tiempos_disponibles = sorted(datos_incendio["tiempo_s"].unique())
# Seleccionamos 5 instantes representativos
instantes = [tiempos_disponibles[0], tiempos_disponibles[len(tiempos_disponibles)//4], 
             tiempos_disponibles[len(tiempos_disponibles)//2], 
             tiempos_disponibles[3*len(tiempos_disponibles)//4], 
             tiempos_disponibles[-1]]

fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.flatten()

for idx, t in enumerate(instantes):
    datos_t = datos_incendio[datos_incendio["tiempo_s"] == t]
    if datos_t.empty:
        continue

    matriz_temp = np.full((100, 100), np.nan)
    for _, row in datos_t.iterrows():
        if row["estado"] != 0:
            matriz_temp[int(row["fila"]), int(row["columna"])] = row["temperatura_c"]

    ax = axes[idx]
    im = ax.imshow(matriz_temp, cmap="hot", origin="lower", vmin=25, vmax=400)
    ax.set_title(f"t = {t} s")
    ax.set_xlabel("Columna")
    ax.set_ylabel("Fila")
    plt.colorbar(im, ax=ax, label="Temperatura (°C)")

if len(instantes) < len(axes):
    axes[-1].axis("off")

fig.suptitle("Mapas de calor de la temperatura del incendio", fontsize=14)
fig.tight_layout()
plt.savefig("mapas_calor.png", dpi=150)
plt.show()

print("Visualizaciones generadas:")
print("  - evolucion_estados.png")
print("  - mapas_calor.png")