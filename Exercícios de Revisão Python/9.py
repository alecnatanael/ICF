import math
import numpy as np

x = float(input('Coordenada x do vetor: '))
y = float(input('Coordenada y do vetor: '))
z = float(input('Coordenada z do vetor: '))
vetor = np.array([[x],
                  [y],
                  [z]])

angulo = float(input('Angulo de rotação no eixo z: '))
t = math.radians(angulo)

R = np.array([[math.cos(t), -math.sin(t), 0],
              [math.sin(t), math.cos(t), 0],
              [0, 0, 1]])

matriz = [[f'{math.cos(t):5.2f}', f'{-math.sin(t):5.2f}', 0],
     [f'{math.sin(t):5.2f}', f'{math.cos(t):5.2f}', 0],
     [f'{0:5.0f}', f'{0:5.0f}', 1]]

matrizrotacao = open('matrizrotacao.dat', 'w')
for linha in matriz:
    linha_texto = ', '.join(map(str, linha))
    matrizrotacao.write(f'[{linha_texto}]\n')

novovetor = R @ vetor
print(novovetor)