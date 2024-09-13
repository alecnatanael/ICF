import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt

def cubic_splineN(xc,fc):
    
    xd = np.array(xc)
    fd = np.array(fc)
    n = len(xd)
    
    ### Declaracao matriz A e vetor b
    Ac = np.zeros([n,n])
    bc = np.zeros([n,1])
    
    ### Primeiro e ultimo elemento da diagonal
    Ac[0,0]= 1.0
    Ac[-1,-1] = 1.0
    
    ###
    
    for i in range(1,n-1):
        
     Ac[i,i-1] = xd[i]-xd[i-1]
     Ac[i,i+1] = xd[i+1]-xd[i]
     ### elementos diagonal
     Ac[i,i] = 2.0*((xd[i]-xd[i-1])+(xd[i+1]-xd[i]))
     ## Vetor b
     bc[i] = 3.0*(((fd[i+1]-fd[i])/(xd[i+1]-xd[i]))-((fd[i]-fd[i-1])/(xd[i]-xd[i-1])))

    print(' Ac = ', Ac)
    
    print('bc = ', bc)
    
    ### Precisamos resolver Ac=b
    cf = np.linalg.solve(Ac,bc)

    af = np.array(fd)
    
    bf = np.zeros([n-1,1])
    df = np.zeros([n-1,1])
    
    for i in range(0,n-1):
        df[i]=(cf[i+1]-cf[i])/(3.0*(xd[i+1]-xd[i]))
        bf[i]=((fd[i+1]-fd[i])/(xd[i+1]-xd[i])) - (((xd[i+1]-xd[i])/3.0)*(2.0*cf[i]+cf[i+1]))
        
    return af,bf,cf,df

def funcinter(ac,bc,cc,dc,xc,xp):
    
    xd = np.array(xc)
    n = len(xd)
    h = [0.0]*(n-1)
    for i in range(0, n-1, 1):
        h[i] = xd[i+1]-xd[i]
        
    ind = 0
    
    ### Obtendo indice da spline que sera usada
    for j in range(0,n-1,1):
        if ((xc[j] <= xp) and (xp <= xc[j+1])):
            ind = int(j)
    
    ### Calculo do valor de S(x)
    spl = ac[ind]+(bc[0][ind]*(xp-xc[ind]))+(cc[0][ind]*pow(xp-xc[ind],2))+(dc[0][ind]*pow(xp-xc[ind],3))
    
    return spl


#### Pontos
x = np.linspace(-1.0,1.0,15)
f = 1.0/(1.0 + 25.0*(x**2))

print('x = ', x)
print('f = ', f)


#### Obtencao dos coefientes dos Splines
a, b, c, d = cubic_splineN(x,f)

bff = np.reshape(b,(1,14))
cff = np.reshape(c,(1,15))
dff = np.reshape(d,(1,14))


print(' a= ', a)
print(' b= ', bff)
print(' c= ', cff)
print(' d = ', dff)

### Intervalo entre pontos dados
xmsh = np.linspace(-1.0,1.0, 500)
funcpl = []

for i in range(0,len(xmsh)):
    funcpl.append(funcinter(a, bff, cff, dff, x, xmsh[i]))
    
plt.scatter(x,f, color='red', label='f(x)')
plt.plot(xmsh, funcpl, color='blue', label= 'S(x)')
plt.grid()
plt.legend()
plt.show()




