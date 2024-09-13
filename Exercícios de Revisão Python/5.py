h0 = float(input('Informe a altura inicial (metros): '))
t = float(input('Informe o tempo em segundos: '))
g = 9.81

h = h0 - (g*(t**2)/2)

print(f'A altura da bola após {t:.2f} segundos é de {h:.2f} metros.')
