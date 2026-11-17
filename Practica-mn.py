import matplotlib.pyplot as plt
import numpy as np

'Aqui es resoldra el problema fent servir el metode de Euler explicit'

#Definim algunes de les constants necessaries per aplicar el metode
N=101  #Nombre de punts del maiat espaial
x0=0
dx=1/(N-1) #Delta de x (separació entre els punts del maiat espaial)
t0=0
dt1=0.51*dx**2 #Delta de t donats a l'enunciat
dt2=0.49*dx**2
dt3=0.25*dx**2
tmax=0.025

#Definirem la constant de normalització (per normalitzar multiplicarem i per desfer la normalització dividirem)
dV=40
sigma=0.472
kappa=0.56
C=2*kappa/(sigma*dV**2)

#Definim els punts del maiat espaial(per poder fer els grafics despres)
x=[]
for u in range(N):
    xin=x0+u*dx #Aquesta esta normalitzada
    xi=xin*2 #Aquesta esta sense normalitzar
    x.append(xi)

#Condicions de contorn
Tc_no=36.5 #En ºC (aquesta és la no normalitzada)
Tc= (Tc_no+273.15)*C #Aqui tenim la normalitzada en K per poder fer els calculs

#Farem les llistes dels dt
dt1=0.51*dx**2 #Deltes donandes a l'enunciat
dt2=0.49*dx**2
dt3=0.25*dx**2
Dt=[dt1,dt2,dt3]

#Per poder fer més d'una grafica (ens interesen la grafica que compara els tres metodes i la que compara cada metode amb l'analitica, ja que un metode esta molt allunyat i no es pot apreciar),farem una llista on guardem les ultimes fileres de cada matriu
T_imp=[]

#Farem una llista de colors pel grafic per poder-lo personalitzar
colors=['mediumspringgreen','violet','deepskyblue']

for dt,color in zip(Dt,colors):#el zip ens permet que hi hagi dues dependencies a la vegada al bucle, i aixi podrem personalitzar els colors del grafic
    g=dt/dx**2 #Calculem la gamma
    Nt=int(tmax/dt)+1 #Calculem el nombre de files que haura de tenir la matriu per asegurar-nos que arribem al temps desitjat
    
    #Ara definirem una matriu de tot 0 
    m=[]
    for k in range(Nt):    
        k = []           
        for i in range(N):  
            k.append(0.0)
        m.append(k)
             
    #Afegim les condicions de contorn (aixi si canviesim Tc es canviarien directament)
    for i in range(N):
     m[0][i]=Tc
    for v in range(Nt):
     m[v][0]=Tc
     m[v][N-1]=Tc

    #Apliquem la formula d'Euler explicit a la matriu
    for n in range(Nt-1): #n és el nombre de fileres, i correspon al subindex de les t
        for j in range(1,N-1): #j és el nombre de columnes i correspon al subindex de les x
            m[n+1][j]=m[n][j]+g*(m[n][j+1]-2*m[n][j]+m[n][j-1])+dt

    #D'aqui nomes ens interesarà l'ultima filera, que correspon a T(x,ta), aixi que definirem
    T_final_n=m[-1][:]

    #Desfem la normalitzacio
    T_final=[]
    for t in T_final_n:
       tf0=t/C #Aquest pas desfa la normalitació
       tf=tf0-273.15 #Aquest pas passa de Kelvin a Celsius per poder comparar millor amb les dades de l'enunciat
       T_final.append(tf)

    T_imp.append(T_final) #Afegim la filera a la llista per despres poder fer diferents grafics si fos necessari

    #Fem els grafics
    plt.plot(x,T_final,color=color,alpha=1,label=f'gamma={g:.2f}')

plt.tick_params(axis='both',which='both',top=True,right=True,direction='in')
plt.legend(fancybox=False,shadow=False,fontsize=10,edgecolor='black')
plt.xlabel('Posició en cm')
plt.ylabel('Temperatura en ºC')
plt.show()

plt.plot(x,T_imp[-1],color='slateblue',alpha=0.5,label='gamma 0.25')
plt.plot(x,T_imp[2],color='blueviolet',alpha=0.5,label='gamma 0.45')
plt.axhline(50, color='indianred')
plt.axvspan(0.75,1.25, color='lavender',label='Zona malalta')
plt.tick_params(axis='both',which='both',top=True,right=True,direction='in')
plt.xlabel('Posició en cm')
plt.ylabel('Temperatura en ºC')
plt.legend(fancybox=False,shadow=False,fontsize=10,edgecolor='black')
plt.show()

    