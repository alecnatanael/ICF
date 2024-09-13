import math

soma1 = soma2 = 0
for n in range(1, 1001):
    soma1 += 1/(n**4)

print(soma1)

for n in range(1, 2001):
    soma2 += 1/(n**4)

print(soma2)

valor = (math.pi**4)/90
print(valor)