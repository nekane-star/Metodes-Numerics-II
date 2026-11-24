import matplotlib.pyplot as plt
import numpy as np

"""Aqui es resoldra el problema fent servir el metode de Euler explicit"""

'Preparació'

#Definim algunes de les constants necessaries per aplicar el metode
N=101  #Nombre de punts del maiat espaial
x0=0
dx=1/(N-1) #Delta de x (separació entre els punts del maiat espaial)
t0=0
tmax=0.025

#Definirem la constant de normalització (per normalitzar multiplicarem i per desfer la normalització dividirem)
dV=40
sigma=0.472
kappa=0.56
C=2*kappa/(sigma*dV**2)

#Definim els punts del maiat espaial(per poder fer els grafics despres)
x=[] #Aquesta serà la llista de x no normalitzades (la que farem servir en general)
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

#Pels errors tambe farem una llista
Errors=[]

#Llista de ordres per les grafiques (aixi ordenem si estan per sobre o per sota de la solució analitica perque es vegin millor)
ordre=[1,3,3]

'Solució analitica'
#Aquesta part la necessitarem per calcular l'error i per comparar els resultats
#Farem servir la mateixa funció que vam defirnir al codi de Solució Analitica
def T(k,q,xn,t_final):
    return k+2*q/(np.pi**3)*sum(((1-(-1)**n)/n**3)*np.sin(n*np.pi*xn)*(1-np.exp(-n**2 *np.pi**2 * t_final)) for n in range(1,100) )

Xn=np.linspace(0,1,N)
Q=1
T_real=(4720/7)*T(Tc,Q,Xn,tmax) -273.15 #Sense normalitzar i en ºC 


'Metode Euler explicit'
for (dt,z) in zip(Dt,ordre):
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

    #Calculem l'error de les funcions i fem el grafic de l'error
    Err=np.abs(T_real - np.array(T_final)) #El np.array es perque es faci l'operació a cada numero de la llista perque sino python no ho enten (tambe es podria fer amb un bucle com la desnosmalització de T_final)
    Errors.append(Err)

    #Fem els grafics  (la idea principal era fer un unic grafic amb totes tres però no s'aprecia be, de manera que es faran 3 grafics)
    plt.plot(x,T_final,color='forestgreen',linestyle='-' ,zorder=z,alpha=1,label=rf'Euler amb $\gamma={g:.2f}$')
    plt.plot(x,T_real,color='slateblue',linestyle='-',label='Solució analitíca')
    plt.tick_params(axis='both',which='both',top=True,right=True,direction='in')
    plt.axhline(50, color='brown',label='Temperatura 50ºC',zorder=1)
    plt.axvspan(0.75,1.25, color='skyblue',alpha=0.75,label='Zona malalta',zorder=0)
    plt.xlabel('Posició (cm)')
    plt.ylabel('Temperatura (ºC)')
    plt.legend(fancybox=False,shadow=False,fontsize=10,edgecolor='black',loc='upper right')
    plt.show()

#Farem una llista de colors per poder personalitzar el grafic
colors=['mediumseagreen','darkviolet','darkcyan']

#Grafiquem els errors
for er,dt,color in zip(Errors,Dt,colors): #el zip ens permet que hi hagi dues dependencies a la vegada al bucle, i aixi podrem personalitzar els colors del grafic
   g=dt/dx**2 #Aixo ho farem per poder posar els labels
   plt.plot(x,er,color=color,label=rf'Euler amb $\gamma={g:.2f}$')
plt.xlabel('Posició (cm)')
plt.ylabel('Error (ºC)')
plt.tick_params(axis='both',which='both',top=True,right=True,direction='in')
plt.legend(fancybox=False,shadow=False,fontsize=10,edgecolor='black',loc='upper right')
plt.show()
#Tambe farem una grafica sense l'error de dt1 per poder comparar millor els altres
plt.plot(x,Errors[-1],color=colors[-1] ,label=r'Euler amb $\gamma=0.25$')
plt.plot(x,Errors[1],color=colors[1], label=r'Euler amb $\gamma=0.49$')
plt.xlabel('Posició (cm)')
plt.ylabel('Error absolut (ºC)')
plt.tick_params(axis='both',which='both',top=True,right=True,direction='in')
plt.legend(fancybox=False,shadow=False,fontsize=10,edgecolor='black')
plt.show()    
