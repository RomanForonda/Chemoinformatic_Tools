import matplotlib.pyplot as plt

x=[1, 2, 3, 4, 5]
y=[2, 4, 6, 8, 10]

# Crear gráfico de puntos
plt.scatter(x, y, color='black', label='Puntos')

# Crear gráfico de línea
plt.plot(x, y, color='red', linestyle='--', label='Línea')

# Añadir título y etiquetas
plt.title("Mi primer gráfico con matplotlib")
plt.xlabel("Eje X")
plt.ylabel("Eje Y")
plt.legend()  # muestra la leyenda
plt.grid(True) # añade rejilla

# Mostrar el gráfico
plt.show()