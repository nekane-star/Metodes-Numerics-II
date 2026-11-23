import numpy as np
import matplotlib.pyplot as plt

#PRIMER APARTAT

# 1. FUNCIÓ GAUSS-SEIDEL

def gauss_seidel(A, b, x0, tol=1e-8, max_iter=10000): #Definim una funció que resol el sistema d'equacions. Considerem un error màxim de 1e-8 per a determinar la convergència de Gauss_Seidel.

    n = len(b) #num columnes matriu A
    x = x0.copy()

    for k in range(max_iter):

        x_old = x.copy() #Guardem una còpia del vector de la iteració anterior per després poder comparar-les i veure si convergeix.

        for i in range(n): #Creem nous valors de x.

            suma = 0

            for j in range(n):

                if j != i: #si j és diferent de i (és a dir, no agafem el terme diagonal de la matriu, perquè aquest terme el deixarem aïllat en calcular x_i).
                    suma += A[i,j]*x[j] #Anem acumulant els productes dels coeficients de la matriu pels valors de les incògnites.

            x[i] = (b[i]-suma)/A[i,i] #Aïllem x_i de l'equació principal del mètode.

        error = np.max(np.abs(x-x_old)) #Calculem la diferència màxima entre els valors de les incògnites de dues iteracions consecutives.

        if error < tol: #Si es compleix, el mètode ha convergit i podem aturar el bucle d'iteracions.
            break

    return x


# 2. DADES DEL PROBLEMA

N = 101 #Nombre de punts del mallat espaial
L = 1 #Longitud del sistema normalitzat
ta = 0.025 

# CONSTANT DE NORMALITZACIÓ + DADES DEL PROBLEMA

dV = 40
sigma = 0.472
kappa = 0.56

C = 2*kappa/(sigma*dV**2) #Normalitzem temperatura

Tc_no = 36.5
Tc = (Tc_no + 273.15)*C

# 3. DISCRETITZACIÓ ESPACIAL I TEMPORAL

dx = L/(N-1) #Distància entre 2 punts del mallat espaial
dt = dx**2

Nt = int(ta/dt) #nombre de passos temporals que necessitem per arribar al temps final

x = np.linspace(0, L, N) #Vector de posicions espaials


# 4. PARÀMETRE r
#r=gamma
r = dt/dx**2 #r=1

# 5. CONSTRUCCIÓ DE LA MATRIU A

M = N-2 #La matriu té dimensions MxM=99x99: les temperatures de les fronteres ja venen dictades per les condicions de contorn.

A = (1+2*r)*np.identity(M) #1+2r=3 #Matriu amb 3 a tota la diagonal
A += -r*np.diag(np.ones(M-1), 1) #Afegim diagonal superior i multipliquem per -r (ve de l'equació d'Euler implícit).
A += -r*np.diag(np.ones(M-1), -1) #Afegim diagonal inferior i multipliquem per -r 
#print(A)

# 6. CONDICIÓ INICIAL

T = np.ones(N)*Tc #Vector de 101 components, totes igual a T_c.

# 7. EVOLUCIÓ TEMPORAL: EULER IMPLÍCIT

for n in range(Nt):

    b = np.zeros(N-2)

    for i in range(N-2):
        b[i] = T[i+1] + dt #Actualitzem el vector de temperatures amb cada iteració

    # Afegim les condicions de contorn al vector de termes independents
    b[0] += r*Tc
    b[-1] += r*Tc

    T_nou = gauss_seidel(A, b, T[1:-1]) #La funció Gauss-Seidel ens retorna el vector de temperatures del nou instant temporal.

    for i in range(N-2): #Actualitzem el vector de temperatures amb els valors que acabem de calcular: anem recorrent els 99 punts interiors del mallat i substituïm les temperatures antigues per les noves.
        T[i+1] = T_nou[i]

# DESFER LA NORMALITZACIÓ (per a fer les gràfiques)

T_real = T/C - 273.15

#REPRESENTACIÓ GRÀFICA
plt.plot(2*x, T_real, color='forestgreen'label='Euler amb gamma=1',)
plt.xlabel("Posició (cm)")
plt.ylabel("Temperatura (ºC)")
plt.axhline(50, color='brown', label='50 ºC')
plt.axvspan(0.75, 1.25, color='skyblue', label='Zona malalta')
plt.ylim(35,55)
plt.tick_params(axis='both',which='both',top=True,right=True,direction='in')
plt.legend()
plt.show()

#ERRORS
#Solució analítica
def sol_analitica(k,q,xn,t_final):
    return k+2*q/(np.pi**3)*sum(((1-(-1)**n)/n**3)*np.sin(n*np.pi*xn)*(1-np.exp(-n**2 *np.pi**2 * t_final)) for n in range(1,100) )

T_final1 = T.copy()
Errors1 = []
T_norm = sol_analitica(Tc, 1, x, ta)
Err = np.abs(T_norm - T_final1)
Errors1.append(Err)

#SEGON APARTAT: dt = 0.5*dx**2
import numpy as np
import matplotlib.pyplot as plt

# 1. FUNCIÓ GAUSS-SEIDEL

def gauss_seidel(A, b, x0, tol=1e-8, max_iter=10000): #Definim una funció que resol el sistema d'equacions. Considerem un error màxim de 1e-8 per a determinar la convergència de Gauss_Seidel.

    n = len(b) #num columnes matriu A
    x = x0.copy()

    for k in range(max_iter):

        x_old = x.copy() #Guardem una còpia del vector de la iteració anterior per després poder comparar-les i veure si convergeix.

        for i in range(n): #Creem nous valors de x.

            suma = 0

            for j in range(n):

                if j != i: #si j és diferent de i (és a dir, no agafem el terme diagonal de la matriu, perquè aquest terme el deixarem aïllat en calcular x_i).
                    suma += A[i,j]*x[j] #Anem acumulant els productes dels coeficients de la matriu pels valors de les incògnites.

            x[i] = (b[i]-suma)/A[i,i] #Aïllem x_i de l'equació principal del mètode.

        error = np.max(np.abs(x-x_old)) #Calculem la diferència màxima entre els valors de les incògnites de dues iteracions consecutives.

        if error < tol: #Si es compleix, el mètode ha convergit i podem aturar el bucle d'iteracions.
            break

    return x


# 2. DADES DEL PROBLEMA

N = 101 #Nombre de punts del mallat espaial
L = 1 #Longitud del sistema normalitzat
ta = 0.025 

# CONSTANT DE NORMALITZACIÓ + DADES DEL PROBLEMA

dV = 40
sigma = 0.472
kappa = 0.56

C = 2*kappa/(sigma*dV**2) #Normalitzem temperatura

Tc_no = 36.5
Tc = (Tc_no + 273.15)*C

# 3. DISCRETITZACIÓ ESPACIAL I TEMPORAL

dx = L/(N-1) #Distància entre 2 punts del mallat espaial
dt = 0.5*dx**2

Nt = int(ta/dt) #nombre de passos temporals que necessitem per arribar al temps final

x = np.linspace(0, L, N) #Vector de posicions espaials


# 4. PARÀMETRE r

r = dt/dx**2 

# 5. CONSTRUCCIÓ DE LA MATRIU A

M = N-2 #La matriu té dimensions MxM=99x99: les temperatures de les fronteres ja venen dictades per les condicions de contorn.

A = (1+2*r)*np.identity(M) #1+2r=3 #Matriu amb 3 a tota la diagonal
A += -r*np.diag(np.ones(M-1), 1) #Afegim diagonal superior i multipliquem per -r (ve de l'equació d'Euler implícit).
A += -r*np.diag(np.ones(M-1), -1) #Afegim diagonal inferior i multipliquem per -r 
#print(A)

# 6. CONDICIÓ INICIAL

T = np.ones(N)*Tc #Vector de 101 components, totes igual a T_c.

# 7. EVOLUCIÓ TEMPORAL: EULER IMPLÍCIT

for n in range(Nt):

    b = np.zeros(N-2)

    for i in range(N-2):
        b[i] = T[i+1] + dt #Actualitzem el vector de temperatures amb cada iteració

    # Afegim les condicions de contorn al vector de termes independents
    b[0] += r*Tc
    b[-1] += r*Tc

    T_nou = gauss_seidel(A, b, T[1:-1]) #La funció Gauss-Seidel ens retorna el vector de temperatures del nou instant temporal.

    for i in range(N-2): #Actualitzem el vector de temperatures amb els valors que acabem de calcular: anem recorrent els 99 punts interiors del mallat i substituïm les temperatures antigues per les noves.
        T[i+1] = T_nou[i]

T_final2 = T.copy()
# DESFER LA NORMALITZACIÓ (per a fer les gràfiques)

T_real = T/C - 273.15

#REPRESENTACIÓ GRÀFICA
plt.plot(2*x, T_real, color='forestgreen',label='Euler amb gamma=0.5')
plt.xlabel("Posició (cm)")
plt.ylabel("Temperatura (ºC)")
plt.axhline(50, color='brown', label='50 ºC')
plt.axvspan(0.75, 1.25, color='skyblue', label='Zona malalta')
plt.ylim(35,55)
plt.tick_params(axis='both',which='both',top=True,right=True,direction='in')
plt.legend()
plt.show()

#ERRORS

#Solució analítica
def sol_analitica(k,q,xn,t_final):
    return k+2*q/(np.pi**3)*sum(((1-(-1)**n)/n**3)*np.sin(n*np.pi*xn)*(1-np.exp(-n**2 *np.pi**2 * t_final)) for n in range(1,100) )


Errors2 = []
Err = np.abs(T_norm - T_final2)
Errors2.append(Err)

# REPRESENTACIÓ DELS ERRORS

plt.plot(2*x, Errors1[0]/C, color='forestgreen',
         label=r'$\Delta t=(\Delta x)^2$')

plt.plot(2*x, Errors2[0]/C, color='royalblue',
         label=r'$\Delta t=0.5(\Delta x)^2$')

plt.xlabel("Posició (cm)")
plt.ylabel("Error absolut (ºC)")

plt.tick_params(axis='both', which='both',
                top=True, right=True, direction='in')

plt.legend()
plt.show()
print(np.max(np.abs(T_final1-T_final2))) #Diferència entre les dues solucions


