import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
k=0.459227
t_final=0.025
x_disc=np.linspace(0,1,101)
# a és la relació entre delta_t i (delta_x)^2
def calcul (a,N):
        
    # Discretitzacions
    x_disc=np.linspace(0,1,N)
    t_disc=x_disc**2
    delta_x= 1/(N-1)
    delta_t= a*(delta_x)**2
    N_t = int(t_final / delta_t)
    #Definim gamma
    g=delta_t/(2*(delta_x)**2)

    # Matriu de temperatures
    T= [[sp.Symbol(f"x_{s}_{d}") for d in range(N)]
    for s in range(N_t +1)]

    # Condicions de contorn: per a tots els temps
    for r in range(N_t + 1):
        T[r][0] = k
        T[r][N-1] = k

    # Condició inicial: per a totes les posicions
    for r in range(N):
        T[0][r] = k
    
    #Creem la matriu A i el vector b (primer les omplim totes de 0)
    A=sp.zeros(N-2,N-2)
    b=sp.zeros(N-2,1)

    #Creem la llista d, és el vector que correpson a les temperatures de la fila n-1 quan estiguem calculant les temperatures de la fila n
    #Inicialment d està ple de k (condicio de contorn) però l'anirem acualitzant cada cop q calculem les files de temperatures
    d=[k]*N

    #Comencem el procés iteratiu, ho fem N-1 vegades perquè ja coneixem la primera fila de temperatures
    for n in range(N_t):
        #Els primers valors de temperatura el calcularem manualment (pq la fòrmula és diferent)
        A[0,0]=1+2*g
        A[0,1]=-g
        b[0,0]=(1-2*g)*d[1]+g*(d[2]+d[0]+d[0]) + delta_t

        #Els valors seguents els calculem amb un bucle (pq tenen la mateixa formula)
        for j in range(N-4):
            A[j+1,j]=-g
            A[j+1,j+1]=1+2*g
            A[j+1,j+2]=-g
            b[j+1,0]=d[j+2]*(1-2*g)+g*(d[j+3]+d[j+1]) + delta_t

        #Els ultims valors els calculem manualment (pq la fòrmula és diferent)
        A[-1,-2]=-g
        A[-1,-1]=1+2*g
        b[N-3,0]=(1-2*g)*d[N-2]+g*(d[N-1]+d[N-1]+d[N-3]) + delta_t

        #Ara tenim la matriu A i el vector b completes, per tant resolem per Gauss-seidel
        llista_inicial=[0]*(N-2)
        
        for j in range(10):
            for i in range(N-2):
                x=(b[i,0]-sum(A[i,m]*llista_inicial[m] for m in range(N-2))+A[i,i]*llista_inicial[i])/A[i,i]
                llista_inicial[i]=x

        #Posem  els valors de T trobats a la matriu T
        for i in range(N-2):
            T[n+1][i+1]=llista_inicial[i]
        #Substituim la llista d
        d=llista_inicial
        d = [k] + llista_inicial + [k]
    #Passem les temperatures normalitzades a temperatures reals
    for j in range(N_t+1):
        for i in range(N):
            T[j][i]=(4720*T[j][i])/7 -273.15
    return(T)
T_1 = calcul(1, 101)
T_05 = calcul(0.5, 101)


plt.plot(x_disc, T_1[-1], label="a = 1")
plt.plot(x_disc, T_05[-1], label="a = 0.5")
plt.plot(x_disc, T_5[-1], label='a=5')
plt.xlabel("x")
plt.ylabel("Temperatura")
plt.title("Temperatura a t = 0.025")
plt.legend()
plt.grid()
plt.show()

        
            


