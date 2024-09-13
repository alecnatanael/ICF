import numpy as np

velocidades = open('velocidades.txt', 'r')
forcavstempo = open('forcavstempo.txt', 'w')
k = 0.03

matriz = np.loadtxt(velocidades, dtype = 'float')
print(matriz)

for i in range(0, 101):
    forca = -k*(matriz[i][1])**2
    forcavstempo.write(f'{matriz[i][1]:.5f}  {forca:.5f}\n')



