import matplotlib.pyplot as plt
import numpy as np

'Aqui es resoldra el problema fent servir el metode de Euler explicit'

#Definim el maiat espaial
N=101  #Nombre de punts del maiat espaial  
x0=0
dx=1/(N-1) #Delta de x (separació entre els punts del maiat espaial)

#Definim el maiat temporal
t0=0
dt1=0.51*dx**2 #Delta de t donats a l'enunciat
dt2=0.49*dx**2
dt3=0.25*dx**2
tmax=0.025

#Definim les gamma per cada dt
g1=dt1/dx**2
g2=dt2/dx**2
g3=dt3/dx**2

#Definim les condicions de contorn (Cada T_n_j t de la malla on n es reereix al temps i j a l'espai)
Tc=36.5
#Definim una matriu de tot 0
m=[]
for j in range(N):
    j=[]
    for i in range(N):
        j.append(0)
    m.append(j)

#Afegim les condicions de contorn a la matriu
for i in range(N):
    m[0][i]=Tc
    m[i][0]=Tc
    m[i][N-1]=Tc

#Definim la nostra formula
for n in range(1,N-2):
    t=t0+n*dt1
    for j in range(1,N-2):
      m[n+1][j]=m[n][j]+g1*(m[n][j+1]-2*m[n][j]+m[n][j-1])+dt1
      x=x0+j*dx
    if t>tmax:
        break


'Repetim ultim pas per una matriu per despres fer el bucle'
#Definim el numero de files que haura de tenir (numero de pases per arribar a tmax)
Nt=int(tmax/dt1)+1
#Definim una matriu de tot 0
m=[]
for j in range(Nt):
    j=[]
    for i in range(N):
        j.append(0)
    m.append(j)

#Afegim les condicions de contorn a la matriu
for i in range(N):
    m[0][i]=Tc
    m[i][0]=Tc
    m[i][N-1]=Tc


#Ara apliquem Euler explicit
for n in range(Nt-1):
    for j in range(1,N-2):
        m[n+1][j]=m[n][j]+g1*(m[n][j+1]-2*m[n][j]+m[n][j-1])+dt1
        
#Fem la llista dels punts de x per poder fer el grafic (aquest es comu a les 3 dt)
x=[]
for j in range(N):
    xi=x0+j*dx
    x.append(xi)

 #Aquesta llista es relacionara amb
 T_final=m[-1] #aixi podrem fer la grafica


    
