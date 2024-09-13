import numpy as np
from numpy.linalg import qr


### Zero de Funções #################################################################

### Bissecção:

def bissec(f, a, b, erro, zero):
    k = 1

    while (b-a) > erro:
        E = (b+a)/2
        if (f(a) < zero and f(b) > zero and f(E) > zero) or (f(a) > zero and f(b) < zero and f(E) < zero):
            b = E
            print(f'Iteração {k} ; x = {E} ; f(x) = {f(E)}')
            k += 1
        elif (f(a) < zero and f(b) > zero and f(E) < zero) or (f(a) > zero and f(b) < zero and f(E) > zero):
            a = E
            print(f'Iteração {k} ; x = {E} ; f(x) = {f(E)}')
            k += 1
        else:
            break
    
    return E


### Newton-Raphson

def newton_raphson(f, x0, erro, phif):
    k = 0

    if abs(f(x0)) < erro:
        k += 1
        print(f'Iteração {k} ; x = {x0} ; f(x) = {f(x0)}')
    else:
        k += 1
        while True:
            x = phif(x0)
            print(f'Iteração {k} ; x = {x} ; f(x) = {f(x)}')
            
            if abs(f(x)) < erro or abs(x-x0) < erro:
                break

            x0 = x
            k += 1
    
    return x

    
# Método das Secantes:

def secantes(f, a, b, erro, phisf):
    x0 = a
    x = b
    k = 0

    if abs(f(x0)) < erro:
        k += 1
        print(f'Iteração {k} ; x = {x0} ; f(x) = {f(x0)}')
    else:
        k += 1
        while True:
            x = phisf(x0, x)
            print(f'Iteração {k} ; x = {x} ; f(x) = {f(x)}')

            if abs(f(x)) < erro or abs(x-x0) < erro:
                break

            m = (x + x0)/2

            if (f(x0) < 0 and f(x) > 0 and f(m) > 0) or (f(x0) > 0 and f(x) < 0 and f(m) < 0):
                x = m
            elif (f(x0) < 0 and f(x) > 0 and f(m) < 0) or (f(x0) > 0 and f(x) < 0 and f(m) > 0):
                x0 = m

            k += 1
    
    return x


#####################################################################################



### Derivação Numérica ##############################################################

def foward_diff(h, x, f): # Argumentos: (precisão, ponto, função)
    fl = (1/h)*(f(x + h) - f(x))
    return fl

def back_diff(h, x, f):
    fl = (1/h)*(f(x) - f(x - h))
    return fl

def central_diff(h, x, f):
    fl = (1/(2*h))*(f(x + h) - f(x - h))
    return fl

def central_diff_5(h, x, f):
    fl5p1 = f(x - 2*h) - 8*f(x-h)
    fl5p2 = 8*f(x+h) - f(x + 2*h)
    fl = (1/(12*h))*(fl5p1 + fl5p2)
    return fl

def diff_2(h, x, f):
    fll = (1/(h**2))*(f(x+h)-(2*f(x))+f(x-h))
    return fll


### Para pontos discretos:

def discret_central_diff(x, xf, xb, fxf, fxb): # "x" que você quer calcular; "xf" posterior ao que você quer calcular; "xb", anterior ao que você quer calcular; "fxf" posterior ao f(x); "fxb" anterior ao f(x)
    h = ((xf - x) + (x - xb))/2
    fl = (1/(2*h))*(fxf - fxb)
    return fl

def discret_foward_diff(x, xf, fx, fxf):
    h = xf - x
    fl = (1/h)*(fxf - fx)
    return fl

def discret_back_diff(x, xb, fx, fxb):
    h = x - xb
    fl = (1/h)*(fx - fxb)
    return fl

def discret_central_diff_5(x, xf, xb, fxf1, fxf2, fxb1, fxb2):
    h = ((xf - x) + (x - xb))/2
    fl5p1 = fxb2 - 8*fxb1
    fl5p2 = 8*fxf1 - fxf2
    fl = (1/(12*h))*(fl5p1 + fl5p2)
    return fl

def discret_diff_2(x, xf, xb, fx, fxf, fxb): ## RUIM!
    h = ((xf - x) + (x - xb))/2
    fll = (1/(h**2))*(fxf-(2*fx)+fxb)
    return fll

#####################################################################################



### Integração Numérica #############################################################


def regra_trapezio(a, b, f): # Argumentos: (início, fim, função)
    h = b-a
    return (h/2)*(f(a)+f(b))

def regra_simpson(a, b, f):
    h = (b-a)/2
    return (h/3)*(f(a)+4*f(a+h)+f(b))


# Argumentos: (início, fim, número de pontos)
def pontos(a, b, m): # ATENÇÃO: m deve ser par para utilizar o método de Simpson
    x = []
    tam = np.abs(b-a)
    h = tam/(m-1)


    for i in range (0, m):
        x.append(a + i*h)
    
    return x


# Argumentos (vetor de pontos, função)
def regra_trapezio_repetida(x, f): 

    h = x[1] - x[0]

    I = f(x[0]) + f(x[len(x)-1]) # Início e fim

    #Integrando os pontos do meio
    for i in range(1, len(x)-1):
        I += 2*f(x[i])

    I = (h/2)*I
    
    return I

def regra_simpson_repetida(x, f):
    h = x[1] - x[0]

    I = f(x[0]) + f(x[len(x)-1]) # Início e fim

    #Integrando os pontos do meio
    for i in range(1, len(x)-1):
        if (i % 2) == 0:
            I += 2*f(x[i])
        else:
            I += 4*f(x[i])

    I = (h/3)*I
    
    return I

#####################################################################################



### Sistemas Lineares ###############################################################


## Métodos Iterativos (Aproximações) ##

def gauss_jacobi(A, B, chute, eps): # Argumentos: (Matriz coeficientes, Matriz resultatos, matriz chute, precisão)
    it = 1
    i = 0
    n = len(A)
    
    x = chute.copy()
    xk = chute.copy()

    print('xk = ', xk)
    print('x = ', x)

    while True:

        for i in range(0, n):   
            soma = 0
            for j in range(0, n):
                if j != i:
                    soma += (A[i][j])*(x[j])
                    
            xk[i] = (1/(A[i][i]))*(B[i] - soma)


        print('Iteracao = ', it)
        print('xk = ', xk)
        print('x = ', x)
        it += 1

        ds = []
        for m in range(0, n):
            ds.append(abs(xk[m]-x[m]))
        
        if (max(ds) < eps):
            break

        x = xk.copy()
        

    return xk


def gauss_seidel(A, B, chute, eps):
    it = 1
    i = 0
    n = len(A)

    
    x = chute.copy()
    xk = chute.copy()

    print('xk = ', xk)
    print('x = ', x)

    while True:

        for i in range(0, n):   
            soma = 0
            
            for j in range(0, n):
                if i > j:
                    soma += (A[i][j])*(xk[j])
                elif i < j:
                    soma += (A[i][j]*x[j])
            xk[i] = (1/(A[i][i]))*(B[i] - soma)


        print('Iteracao = ', it)
        print('xk = ', xk)
        print('x = ', x)
        it += 1

        ds = []
        
        for m in range(0, n):
            ds.append(abs(xk[m]-x[m]))
        
        if (max(ds) < eps):
            break

        x = xk.copy()
        

    return xk

#######################################


## Métodos diretos ####################

# Argumentos: (Matriz coeficientes, Matriz resultatos)
def subs_retro(A, B): ## Resolve um sistema triangular
    n = len(A)
    x = []

    for i in range(0, n):   
        x.append(0)

    x[n-1] = B[n-1]/A[n-1][n-1]

    c = 1
    
    for i in range(n-2, -1, -1):
        soma = 0
        c += 1
        for j in range(n-1, i, -1):
            soma += (A[i][j]*x[j])
        x[n-c] = (B[i] - soma)/(A[i][i])

    return x

def eliminacao_Gauss(A, B): ## Transforma um sistema em triangular
    n = len(A)
    for k in range(0, n-1):
        for i in range(k+1, n):
            m = (A[i][k])/(A[k][k])

            for j in range(k, n):
                A[i][j] = A[i][j] - m*A[k][j]
            B[i] = B[i] - m*B[k]
    
    return A, B

def solucao_direta(A, B): ## Solução direta de sistemas lineares
    n = len(A)

    # Triangulando com a Eliminação de Gauss
    for k in range(0, n-1):
        for i in range(k+1, n):
            m = (A[i][k])/(A[k][k])

            for j in range(k, n):
                A[i][j] = A[i][j] - m*A[k][j]
            B[i] = B[i] - m*B[k]


    # Resolvendo o sistema triangular
    x = []

    for i in range(0, n):
        x.append(0)

    x[n-1] = B[n-1]/A[n-1][n-1]

    c = 1
    
    for i in range(n-2, -1, -1):
        soma = 0
        c += 1
        for j in range(n-1, i, -1):
            soma += (A[i][j]*x[j])
        x[n-c] = (B[i] - soma)/(A[i][i])

    return x


#####################################################################################


### Mínimos Quadrados ###############################################################

### Gs úteis:

def funcgs2(xc): ## para g1 = 1; g2 = x;

    m = len(xc)

    g1 = np.array([1.0]*m)
    g2 = np.array(xc)

    gm = [g1, g2]

    return gm

def funcgs3(xc): ## para g1 = 1; g2 = x; g3 = x^2

    m = len(xc)

    g1 = np.array([1.0]*m)
    g2 = np.array(xc)
    g3 = np.array(xc)**2

    gm = [g1, g2, g3]

    return gm

def funcsen(xc): ## para g1 = 1; g2 = sin(pi*x/2)

    n = len(xc)

    g1 = np.array([1.0]*n)
    
    g2 = []
    for i in range(0, n):
        g2.append(np.sin((np.pi*xc[i])/2))

    gm = [g1, g2]

    return gm



## Resolução:

def montarMatriz(gs, xc, fx): # Monta o sistema linear a ser resolvido
    n = len(gs)
    m = len(xc)
    A = np.zeros((n, n))
    B = np.array([0.0]*n) ## é realmente importante especificar que é um float com 0.0 !

    for i in range(0, n):
        for j in range(0, n):
            for k in range(0, m):
                A[i][j] += gs[j][k]*gs[i][k]
    
    for i in range(0, n):
        soma = 0
        for k in range(0, m):
            soma += fx[k]*gs[i][k]
        B[i] = soma

    return A, B


def minimos_quadrados(xc, fx, gm): ## Função geral para obter os coeficientes alphas pelo método mínimos quadrados

    # Montando a matriz:
    A, B = montarMatriz(gm, xc, fx)

    # Resolvendo o sistema linear, encontrando os alphas:
    alphas = solucao_direta(A, B)

    return alphas


def erro_min_quadrados(xc, fx, gm, alpha):

    n = len(gm)
    m = len(xc)
    
    phis = []
    for k in range(0, m):
        phix = 0.0
        for i in range(0, n):
            phix += alpha[i]*gm[i][k]
        phis.append(phix)


    d = 0
    for k in range(0, m):
        d += (fx[k] - phis[k])**2
    
    return d

#####################################################################################


### Interpolação ####################################################################

def interpolacao(x, fx): ### Retorna os coeficientes "a" do polinômio Pn(x)
    n = len(x)
    A = np.zeros((n, n))
    for i in range(0, n):
        for j in range(0, n):
            A[i][j] = x[i]**j

    B = np.copy(fx)

    coef = solucao_direta(A, B)


    return coef

def polinomio(x, a): ## Retorna um valor correspondente a Pn(x) para um polinômio com os coeficientes do vetor "a"
    n = len(a)


    Pnx = 0
    for j in range(0, n):
        Pnx += a[j]*(x**j)
        
    
    return Pnx


### Outros Métodos:

def Polinomio_Lagrange(x, xref, yref):
    n = len(xref)
    P = 0

 
    for k in range(0, n):

        # Produtório do Numerador:
        produto1 = 1
        for j in range(0, n):
            if (j != k):
                produto1 *= (x - xref[j])

        # Produtório do Denominador:
        produto2 = 1
        for j in range(0, n):
            if (j != k):
                produto2 *= (xref[k] - xref[j])
        
        Lk = produto1/produto2
        
        P += yref[k]*Lk
    
    return P ### Retorna um ponto de P(x), polinômio formado com xref e yref


def Polinomio_Newton(x, xc, f):
    n = len(xc)

    # Obtendo operadores de diferença
    opn = np.zeros([n, n])
    opn[0] = f
    
    # Obtendo Operadores de diferença dividida ordem 1
    for k in range(0, n-1):
        opn[1][k] = ((f[k]-f[k+1])/(xc[k]-xc[k+1]))

    # Obtendo matriz de operadores até a ordem n:
    for i in range(2, n):
        for k in range(0, n-i):
            opn[i][k] = (opn[i-1][k+1]-opn[i-1][k])/(xc[k+i]-xc[k])
    # Cada linha é uma lista de operadores de ordem correspondente a da linha

    # Resolvendo o polinômio com determinado valor x:

    Pn = f[0]
    prodx = (x-xc[0])
    for i in range(1, n):
        Pn += prodx*opn[i][0]
        prodx *= (x-xc[i])

    return Pn


### Splines:

def Coef_Splines(x, fx):

    a = np.copy(fx) # Obtendo 'a'
    n = len(x) 
    h = [] # Obtendo 'h'
    for j in range(0, n-1):
        h.append(x[j+1] - x[j])


    # Obtendo a matriz A:
    A = np.zeros((n, n))
    for i in range(0, n):
        for j in range(0, n):
            if ((i == 0) and (j == 0)) or ((i == n-1) and (j == n-1)):
                A[i][j] = 1
            elif i==j:
                A[i][j] = 2*(h[i-1]+h[j])
                A[i][j-1] = h[i-1]
                A[i][j+1] = h[j]

    # Obtendo a matriz B:
    B = [0]*(n)

    for j in range(1, n-1):
        B[j] = (3/h[j])*(a[j+1]-a[j]) - (3/h[j-1])*(a[j]-a[j-1])


    # Resolvendo para obter c:

    c = np.linalg.solve(A,B)

    # Obtendo b e d:
    d = []
    b = []
    for j in range(0, n-1):
        b.append(((a[j+1]-a[j])/h[j]) - (h[j]/3)*(2*c[j]+c[j+1]))
        d.append((c[j+1]-c[j])/(3*h[j]))

    return a, b, c, d # Obtem os vetores a, b, c e d


def Ponto_Splines(a, b, c, d, xc, x): # Retorna um valor correspondende a determinado x que é parte de determinada spline

    xd = np.array(xc)
    n = len(xd)

    ind = 0

    ### Obtendo indice da spline que sera usada
    for j in range(0, n-1, 1):
        if ((xc[j] <= x) and (x <= xc[j+1])):
            ind = int(j)

    ### Calculo do valor de S(x)
    S = a[ind]+(b[ind]*(x-xc[ind]))+(c[ind]*pow(x-xc[ind],2))+(d[ind]*pow(x-xc[ind],3)) 

    return S


def Obter_Splines(a, b, c, d, x, precisao): # Retorna uma matriz onde cada linha é o conjunto de cordenadas y associdados a pontos de uma Spline (cada linha é uma Spline)
    n = len(x)

    S = np.zeros((n-1, precisao))
    for j in range(0, n-1):
        xm = np.linspace(x[j], x[j+1], precisao)
        for i in range(0, precisao):
            S[j][i] = a[j] + b[j]*(xm[i]-x[j]) + c[j]*((xm[i]-x[j])**2) + d[j]*((xm[i]-x[j])**3)
            
    return S


#####################################################################################



### Autovalores ####################################################################

def m_potencia(A, y0, iteracoes):
    n = len(A)

    z = A @ y0 # Produto de Matrizes
    for i in range(1, iteracoes):
        alpha = z[0]
        for j in range(0, n):  # Máxima componente
            if z[j] > alpha: alpha = z[j]

        y = z/alpha
        z = A @ y

    return alpha


def m_potencia_inv(A, y0, k):
    A_inv = np.linalg.inv(A)
    n = len(A)

    z = A_inv @ y0 # Produto de Matrizes
    for i in range(1, k):
        alpha = z[0]
        for j in range(0, n):  # Máxima componente
            if z[j] > alpha: alpha = z[j]

        y = z/alpha
        z = A_inv @ y

    return 1/alpha


def autovalores_qr(MatrizA, Q, R, precisao):
    
    An = np.copy(MatrizA)

    # Iterando até obter a diagonal com autovalores
    maior = precisao + 1.0
    while (maior > precisao):
        # Obtendo o maior valor absoluto de A
        maior = abs(An[1][0])  
        for i in range(0, len(An)):
            for j in range(0, len(An)):
                if (i > j):
                    if (abs(An[i][j]) > maior):
                        maior = abs(An[i][j])
        
        An = R @ Q  # Próximo A
        Q, R = qr(An) # Novos Q, R

    
    # Obtendo os autovetores
    Autovalores = []
    for i in range(0, len(An)):
        for j in range(0, len(An)):
            if (i == j):
                Autovalores.append(An[i][j])
    

    return Autovalores


#####################################################################################


### Outros ##########################################################################

# Obtem pontos pro gráfico baseado no polinômio de coeficientes "a"
def obtempontosP(x, a, prec): # Argumentos: pontos de referência "x", vetor de coeficientes "a", precisão (sugerido: 100)
    xm = np.linspace(x[0], x[-1], prec) ## cria uma "malha fina" de "prec" pontos do primeiro x ao último

    fi = []

    for i in range(0, len(xm)):
        soma = 0
        for j in range(0, len(a)):
            soma += a[j]*(xm[i]**j)
        fi.append(soma)
    
    return xm, fi
