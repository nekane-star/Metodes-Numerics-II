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
for n in [0,N-2]:
    t=t0+n*dt1
    for j in [1,N-2]:
      m[n+1][j]=m[n][j]+g1*(m[n][j+1]-2*m[n][j]+m[n][j-1])+dt1
      x=x0+j*dx
      print(m[n+1][j],"es el nostre nou punt i es troba a:",n+1,j)
    if t>tmax:
        break


    
