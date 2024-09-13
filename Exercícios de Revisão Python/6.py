import math

r = float(input('Raio da circunferência: '))
angle = float(input('Ângulo de giro (graus): '))
theta = math.radians(angle)

if angle == 90:
    x = 0.0
    y = 1.
x = r * math.cos(theta)
y = r * math.sin(theta)

print(f'{x}, {y}')