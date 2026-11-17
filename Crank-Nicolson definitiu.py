#LLibreries utilitzades
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

#Condició de contorn i inicial (En Celsius)
T_contorn=36.5

#Condició de contorn normalitzada
k=(T_contorn+273.15)*(7/4720)

#Definim el temps final, el nombre de punts i l'espai x (normalitzat)

t_final=0.025
N=101
x=np.linspace(0,1,N)

# a és la relació entre delta_t i (delta_x)^2

#Programem la funció que té com a output la matriu de temperatures a qualsevol posició x i temps t i té com a input el paràmetre que relaciona les dues discretitzacions normalitzades (de x i de t)
def calcul (a):
       
    # Discretitzacions
    x_disc=np.linspace(0,1,N)
    t_disc=x_disc**2
    delta_x= 1/(N-1)
    delta_t= a*(delta_x)**2

    #Nombre de punts pel temps
    N_t = int(t_final / delta_t)

    #Definim gamma
    g=delta_t/(2*(delta_x)**2)

    #Matriu de temperatures
    T= [[sp.Symbol(f"x_{s}_{d}") for d in range(N)]
    for s in range(N_t +1)]

    #Condicions de contorn per a tots els temps
    for r in range(N_t + 1):
        T[r][0] = k
        T[r][N-1] = k

    #Condició inicial per a totes les posicions
    for r in range(N):
        T[0][r] = k
   
    #Creem la matriu A i el vector b (primer les omplim totes de 0). Aquests elements ens serviran per iterar el mètode de Gauss_seidel
    A=sp.zeros(N-2,N-2)
    b=sp.zeros(N-2,1)

    #Creem la llista d, és el vector que correpson a les temperatures de la fila n-1 quan estiguem calculant les temperatures de la fila n
    #Inicialment d està ple de k (condicio de contorn i inicial) però l'anirem acualitzant cada cop q calculem les files de temperatures
    d=[k]*N

    #Comencem el procés iteratiu, ho fem N-1 vegades perquè ja coneixem la primera fila de temperatures (CI)
    for n in range(N_t):
        #Els primers valors de temperatura el calcularem manualment (perquè la fòrmula és diferent)
        A[0,0]=1+2*g
        A[0,1]=-g
        b[0,0]=(1-2*g)*d[1]+g*(d[2]+d[0]+d[0]) + delta_t

        #Els valors seguents els calculem amb un bucle (perquè tenen la mateixa formula)
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

        #Creem una llista que consisteix en agafar tots els valors de d exceptuant el primer i el ultim
        llista_inicial=d[1:-1]

        #Iterem 20 cops el mètode (Hem considerat que iterar-lo 20 cops és suficient, a més, amb 20 iteracions el codi ja triga molt en donar les gràfiques)
        #Cal recalcar que hem fet varies proves considerant un nombre d'iteracions més gran i el resultat ha sortit idèntic al d'utilitzar solament 20 iteracions
        for j in range(20):

            #Apliquem les fòrmules d'iteració del mètode de Gauss_Seidel
            for i in range(N-2):
                x=(b[i,0]-sum(A[i,m]*llista_inicial[m] for m in range(N-2))+A[i,i]*llista_inicial[i])/A[i,i]
                llista_inicial[i]=x
           
        #Posem  els valors de T trobats a la matriu T
        for i in range(N-2):
            T[n+1][i+1]=llista_inicial[i]

        #Actualitzem la llista d per tal que per la pròxima iteració el codi utilitzi els valors trobats a l'anterior iteració per calcular els nous valors
        d=llista_inicial

        #Afegim dos valors a l'inici i al final de la llista perquè així quan el codi faci la següent iteració agafi els valors de llista_inicial
        d = [k] + llista_inicial + [k]

    #Passem les temperatures normalitzades a temperatures reals
    for j in range(N_t+1):
        for i in range(N):
            T[j][i]=(4720*T[j][i])/7 -273.15
    return(T)

#Utilitzem la funció per definir les matrius de temperatura corresponents als dos valors de a
T_1 = calcul(1)
T_05 = calcul(0.5)

#SOLUCIÓ ANALÍTICA#

#Definim q
q=1

#Utilitzem l'expressió de l'annex i substituim el valor de q i de k per calcular la temperatura (normalitzada) en funció de x (normalitzat) a t=0.025 (normalitzat)
T_norm=k+2*q/(np.pi**3)*sum(((1-(-1)**n)/n**3)*np.sin(n*np.pi*x)*(1-np.exp(-n**2 *np.pi**2 * t_final)) for n in range(1,100) )

#Utilitzem la fòrmula de la conversió entre temperatura normalitzada i temperatura real
T_real=(4720/7)*T_norm -273.15

#Representem gràficament la solució analítica i les aproximacions
plt.plot(x, T_1[-1], label="a = 1")
plt.plot(x, T_05[-1], label="a = 0.5")
plt.plot(x,T_real,color='black',label='solució analítica')
plt.xlabel("x")
plt.ylabel("Temperatura")
plt.title("Temperatura a t = 0.025")
plt.legend()
plt.grid()
plt.show()

#ERRORS NUMÈRICS#

#Definim els errors numèrics
e_05=np.abs(T_real-T_05[-1])
e_1=np.abs(T_real-T_1[-1])

#Representem gràficament els errors numèrics
plt.plot(x,e_05,label='error amb a=0.5',color='red')
plt.plot(x,e_1,label='error amb a=1',color='blue')
plt.legend()
plt.grid()
plt.xlabel('x discretitzat')
plt.ylabel('Error numèric')
plt.title('Errors numèrics amb a=0.5 i a=1')
plt.show()

