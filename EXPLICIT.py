import numpy as np
import matplotlib.pyplot as plt

# initial condition
from numpy.core._multiarray_umath import ndarray

Pi = int(input('Masukkan nilai Pi ='))
Pleft = int(input('Masukkan nilai Pleft ='))
Pright = int(input('Masukkan nilai Pright ='))
Alpha = int(input('Masukkan nilai Alpha ='))
dt = float(input('Masukkan nilai dt ='))
dx = float(input('Masukkan nilai dx ='))
beta: float = (Alpha * dx**2) / dt
y = int(input('Masukkan nilai index time ='))
z = int(input('Masukkan jumlah titik ='))
P: ndarray = np.random.rand(y,z)
for i in range (z):
      P[0,i] = Pi

# boundary condition
print('Batas Kiri-Kanan')
print()
print('Jika Dirichlet-Dirichlet, masukkan angka 1; Jika Dirichlet-Neumann, masukkan angka 2,'
      ' Jika Neumann-Dirichlet, masukkan angka 3; Jika Neumann-Neumann, masukkan angka 4')
x = int(input('Batas Kiri-Kanan ='))
t = np.arange(0.0,z,1)


# dirichlet-dirichlet
if x == 1:
      for i in range (0,y-1):
            P[i+1,0] = Pleft
            for j in range (0,z-2):
                  P[i+1,z-1] = Pright
                  P[i+1,j+1] = (P[i,j]+(-2+beta)*P[i,j+1]+P[i,j+2])/beta
            plt.plot(t,P[i,])
      plt.show()
print(P)


# dirichlet-neumann
if x == 2:
      for i in range (0,y-1):
            P[i+1,0] = Pleft
            for j in range (0,z-2):
                  P[i+1,z-1] = ((-2+beta)*P[i,j]+2*P[i,j-1]-Pright*2*dx)/beta
                  P[i+1,j+1] = (P[i,j]+((-2+beta)*P[i,j+1])+P[i,j+2])/beta
            plt.plot(t, P[i,])
      plt.show()
print(P)

# neumann-dirichlet
if x == 3:
      for i in range (0,y-1):
            P[i+1,0] = ((-2+beta)*P[i,0]+2*P[i,1]-Pleft*2*dx)/beta
            for j in range (0,z-2):
                  P[i+1,z-1] = Pright
                  P[i+1,j+1] = (P[i,j]+(-2+beta)*P[i,j+1]+P[i,j+2])/beta
            plt.plot(t,P[i,])
      plt.show()
print(P)

# neumann-neumann
if x == 4:
      for i in range (0,y-1):
            for j in range (0,z-2):
                  P[i+1,0] = ((-2+beta)*P[i,j+1]+2*P[i,j+1]-(Pleft*2*dx))/beta
                  P[i+1,z-1] = ((-2+beta)*P[i,j]+2*P[i,j-1]-Pright*2*dx)/beta
                  P[i+1,j+1] = (P[i,j]+(-2+beta)*P[i,j+1]+P[i,j+2])/beta
            plt.plot(t,P[i,])
      plt.show()
print(P)

