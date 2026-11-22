import numpy as np
import matplotlib.pyplot as plt

#PRIMER APARTAT

# 1. FUNCIO GAUSS-SEIDEL

def gauss_seidel(A, b, tol=1e-8, max_iter=10000):

    n = len(b) #num columnes matriu A
    x = np.zeros(n)

    for k in range(max_iter):

        x_old = x.copy()

        for i in range(n):

            suma = 0

            for j in range(n):

                if j != i: #si j és diferent de i
                    suma += A[i,j]*x[j] #sumem contribucions

            x[i] = (b[i]-suma)/A[i,i]

        error = np.max(np.abs(x-x_old))

        if error < tol:
            break

    return x


# 2. DADES DEL PROBLEMA

N = 101
L = 1
ta = 0.025

# CONSTANT DE NORMALITZACIO

dV = 40
sigma = 0.472
kappa = 0.56

C = 2*kappa/(sigma*dV**2)

Tc_no = 36.5
Tc = (Tc_no + 273.15)*C

# 3. DISCRETITZACIO ESPACIAL I TEMPORAL

dx = L/(N-1)
dt = dx**2

Nt = int(ta/dt)

x = np.linspace(0, L, N)


# 4. PARAMETRE r

r = dt/dx**2 #r=1

# 5. CONSTRUCCIO DE LA MATRIU A

M = N-2 #matriu és MxM

A = (1+2*r)*np.identity(M) #1+2r=3 #matriu diagonal 3
A += -r*np.diag(np.ones(M-1), 1) #diagonal superior
A += -r*np.diag(np.ones(M-1), -1) #diagonal inferior

# 6. CONDICIO INICIAL

T = np.ones(N)*Tc

# 7. EVOLUCIO TEMPORAL: EULER IMPLICIT

for n in range(Nt):

    b = np.zeros(N-2)

    for i in range(N-2):
        b[i] = T[i+1] + dt

    # Afegim les condicions de contorn
    b[0] += r*Tc
    b[-1] += r*Tc

    T_nou = gauss_seidel(A, b)

    for i in range(N-2):
        T[i+1] = T_nou[i]

# DESFER LA NORMALITZACIO

T_real = T/C - 273.15

#REPRESENTACIÓ GRÀFICA
plt.plot(2*x, T_real, color='forestgreen')
plt.xlabel("Posició (cm)")
plt.ylabel("Temperatura (ºC)")
plt.axhline(50, color='brown', label='50 ºC')
plt.axvspan(0.75, 1.25, color='skyblue', label='Zona malalta')
plt.ylim(35,55)
plt.tick_params(axis='both',which='both',top=True,right=True,direction='in')
plt.legend()
plt.show()
