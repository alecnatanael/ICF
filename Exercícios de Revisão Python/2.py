from sys import exit

def maiornumero (numeros):
    maior = numeros[0]
    for num in numeros:
        if num > maior:
            maior = num

    return maior


first = float(input('Primeiro número real: '))
secnd = float(input('Segundo número real: '))
third = float(input('Terceiro número real: '))

if first == secnd or third == secnd or first == third:
    exit()

numeros = [first, secnd, third]

maior = maiornumero(numeros)
print(maior)



