import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
import sys

# Функция по варианту
def f1(x, y):
	return 3 * np.exp(-((x - 2) / 1) ** 2 - ((y - 3) / 2) ** 2) + 1 * np.exp(-((x - 1) / 2) ** 2 - ((y - 1) / 1) ** 2)

# Квадратичная функция
def f2(x, y):
	return 100 * (y - x) * (y - x) + (1 - x) * (1 - x)

# Функция Розенброка
def f3(x, y):
	return 100 * (y - x ** 2) ** 2 + (1 - x) ** 2

# Читаем точки и название метода и функции
def input(filename):
	x = []
	y = []

	with open (filename) as file:
		for line in file:
			x_i, y_i = map(float, line.split(" "))
			x.append(x_i)
			y.append(y_i)

	return x, y

# Строим сетку
def build_mesh(x_vals, y_vals, func):
	x = np.arange(min(x_vals) - 5, max(x_vals) + 5, 0.05)
	y = np.arange(min(y_vals) - 5, max(y_vals) + 5, 0.05)

	X, Y = np.meshgrid(x, y)

	Z = func(X, Y)

	return X, Y, Z

# Отрисовка линии
def draw_line(x, y, name):
	plt.title(name)
	plt.scatter(x, y, s = 20, color = "blue")
	plt.plot(x, y, color = "cyan")
	plt.xlabel("X")
	plt.ylabel("Y")

def main():
	name = sys.argv[1]
	filename = sys.argv[2]
	x, y = input(filename)

	X, Y, Z = build_mesh(x, y, f1)

	residual = pl.contourf(X, Y, Z, levels = 50)
	plt.colorbar(residual, location = "right")

	plt.contour(X, Y, Z)
	draw_line(x, y, name)
	plt.savefig("images\\" + name + ".png")

if __name__ == "__main__":
	main()