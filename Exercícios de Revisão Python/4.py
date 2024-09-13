primos = [2]
for i in range(3, 101):
    statusprimo = True
    for j in range(2, i):
        if (i % j) == 0:
            statusprimo = False
            break
    if statusprimo == True:
        primos.append(i)

print(primos)

arquivo = open("texto.txt", "w")
for num in primos:
    arquivo.write(f'{num} ')