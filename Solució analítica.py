import numpy as np
import matplotlib.pyplot as plt

#SOLUCIÓ ANALÍTICA#

#Utilitzem l'expressió de l'annex i substituim el valor de q i de k per calcular la temperatura (normalitzada) en funció de x (normalitzat) a t=0.025 (normalitzat) per definir una funció
def T(k,q):
    return k+2*q/(np.pi**3)*sum(((1-(-1)**n)/n**3)*np.sin(n*np.pi*x)*(1-np.exp(-n**2 *np.pi**2 * t_final)) for n in range(1,100) )

#Condició de contorn i inicial (En Celsius)
T_contorn=36.5

#Condició de contorn normalitzada
tc=(T_contorn+273.15)*(7/4720)

#Altres paràmetres
N=101
Q=1
t_final=0.025
x=np.linspace(0,1,N)

T_norm=T(tc,Q)
#Utilitzem la fòrmula de la conversió entre temperatura normalitzada i temperatura real
T_real=(4720/7)*T_norm -273.15
print(T_real)

#Representem gràficament la solució analítica i les aproximacions
plt.plot(x,T_real,color='black',label='solució analítica')
plt.xlabel("x")
plt.ylabel("Temperatura")
plt.title("Temperatura a t = 0.025")
plt.legend()
plt.show()
