import matplotlib.pyplot as plt
import numpy as np

#Per resoldre el problema farem el mateix que amb el metode de Euler explicit pero afegirem una comprovació al bucle per asegurar-nos que mai es crema la part sana

'Dades'
N=101 #Nombre de punts del maiat espaial
x0=0
t0=0
tmax=0.025 #Sabem que aquest tmax és més gran que el temps necessari de manera que el farem servir per definir la matriu, i aixi ens asegurarem que sigui lo suficientment gran com per fer arribar a la solució

#Diferencies del maiat i gamma
dx=1/(N-1)
dt=0.25*dx**2
g=dt/dx**2

#Definim els punts del maiat espaial(per poder fer els grafics despres)
x=[] #Aquesta serà la llista de x no normalitzades (la que farem servir en general)
xn=[] #Aquesta serà la llista de x normalitzades 
for u in range(N):
    xin=x0+u*dx #Aquesta esta normalitzada
    xi=xin*2 #Aquesta esta sense normalitzar
    x.append(xi)
    xn.append(xin)

#Constants de normalització 
dV=40
sigma=0.472
kappa=0.56
cv=3698
ro=1081
C=2*kappa/(sigma*dV**2) #Normalització de la temperatura
K= 0.02**2*ro*cv/kappa #Normalització del temps

#Condicions de contorn
Tc_no=36.5 #En ºC (aquesta és la no normalitzada)
Tc= (Tc_no+273.15)*C #Aqui tenim la normalitzada en K per poder fer els calculs

#Condició per la comprovació
TM_no=50 #En ºC i no normalitzada
TM=(TM_no+273.15)*C #En K i normalitzada

'Preparacio--Matriu'
Nt=int(tmax/dt)+1 #Nombre de files de la matriu

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

'Preparació--Condició'
#Sabem que els extrems de la regio malalta son x=0,75 i x=1,25. Hem de trobar a quina columna pertanyen per poder posar la condició

#Cal que trobem en quin punt de la llista es troben les fronteres (per aixo farem un bucle de comprovació)
for i in range(len(x)):
    if x[i]<=0.75001: #Fiquem el 001 perque a vegades python pot escriure els numeros com 0,0000000001 i aquest al calcul l'hauria d'interpretar com 0
        ind_front_1=i #Aixo ens donarà l'index del punt del maiat més proper a 0,75 però sense pasar-lo

for i in range(len(x)):
    if x[i]>=1.2499: #Igual que abans fem el marge de seguretat
        ind_front_2=i
        break #Aixi fem que pari exactament al punt més proper a 1,25 per la dreta
#Aquests ind_front els posarem a la condició perque el codi s'aturi un cop arribin a la temperatura maxima


'Solució'
#Apliquem la formula d'Euler explicit a la matriu
for n in range(Nt-1):
    for j in range(1,N-1):
        m[n+1][j]=m[n][j]+g*(m[n][j+1]-2*m[n][j]+m[n][j-1])+dt
       
    T_final_n=m[n+1][:] #Això ens redefinira T_final_n cada cop com l'ultima fila que ha completat. Quan pari el codi T_final_n serà la ultima fila amb numeros de la matriu, i la fila que ens interesa
    #Calculem en quin temps ens trobem
    t_final_n=t0+n*dt
    #Afegim la comprovació
    if m[n+1][ind_front_1]>=TM:
        break
    if m[n+1][ind_front_2]>=TM:
        break



#Desfem la normalitzacio
T_final=[]
for t in T_final_n:
    tf0=t/C #Aquest pas desfa la normalitació
    tf=tf0-273.15 #Aquest pas passa de Kelvin a Celsius per poder comparar millor amb les dades de l'enunciat
    T_final.append(tf)

'Gràfica'
plt.plot(x,T_final,color='forestgreen',label=r'Solució per Euler $\gamma = 0.25$',zorder=3)
plt.tick_params(axis='both',which='both',top=True,right=True,direction='in')
plt.axhline(50, color='brown',label='Temperatura 50ºC')
plt.axvspan(0.75,1.25, color='skyblue',alpha=0.75,label='Zona malalta')
plt.xlabel('Posició (cm)')
plt.ylabel('Temperatura (ºC)')
plt.legend(fancybox=False,shadow=False,fontsize=10,edgecolor='black')
plt.show()

#Desfem la normalització del temps
t_final=t_final_n*K
print(f'{t_final:.6f}')
        