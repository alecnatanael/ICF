dicionario = {'esfera':240, 'cilindro':100, 'paralelepipedo':20}
g = 9.81

for nome, massa in dicionario.items():
    peso = massa*g
    print(f'O peso exercido sobre o(a) {nome} é de {peso:.2f} Newtons.')