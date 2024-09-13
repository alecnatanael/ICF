import numpy as np



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

#####################################################################################



### Integração Numérica #############################################################

# Argumentos: (início, fim, função)
def regra_trapezio(a, b, f): 
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

def solucao_direta(A, B): ## Junta os dois passos anteriores em uma função só
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

def funcgs3(xc): ## para g1 = 1; g2 = x; g3 = x^2

    m = len(xc)

    g1 = np.array([1.0]*m)
    g2 = np.array(xc)
    g3 = np.array(xc)**2

    gm = [g1, g2, g3]

    return gm

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


def minimos_quadrados_3gs(xc, fx, gm): ## Função geral para obter os coeficientes alphas pelo método mínimos quadrados

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
