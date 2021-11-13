#Hazman Hatadi_101317004
#Hana Khayyirah H_101317024
#Innocentia Ukkana_101317070
#Intan Umull Magfira_101317104

import numpy as np
from numpy.linalg import inv

import matplotlib.pyplot as plt

#Input Data
PP = float(input("Masukkan Nilai Panjang Pipa: "))
dx = float(input("Masukkan Nilai dx:  "))
dt = float(input ("Masukkan Nilai dt: " ))
Al = float(input("Masukkan Nilai Aplha: "))
TF = float(input("Masukkan Nilai Time Final: "))
JN = int(PP/dx)+1                                  #Jumlah_N
JTS = int(TF/dt)+1                                 #Jumlah_Time_Step
b = Al*dx**2 /dt


#Boundary Condition
PR = float(input("Masukkan Nilai Pressure Right: "))
PL = float(input("Masukkan Nilai Pressure Left: "))
Pi = float(input("Masukkan Nilai Pressure Intial:"))

#IMPLISIT

KS = float(0.5*Al*(dx**2))
print('Kriteria Stabilitas = ', KS)

if dt > KS:
    print('Data Salah! Masukkan ulang data')

# Looping for P(x,t>0)
if dt <= KS:

    #Boundary Condition

    print('Batas Kiri-Kanan')
    print()
    print('Jika Dirichlet-Dirichlet, masukkan angka 1; Jika Dirichlet-Neumann, masukkan angka 2,'
          'Jika Neumann-Dirichlet, masukkan angka 3; Jika Neumann-Neumann, masukkan angka 4')
    BKK = int(input('Batas Kiri-Kanan ='))
    X = np.random.rand(JTS, JN)

    for i in range(JN):
        X[0, i] = Pi

    A = np.random.rand(JN,)
    W = np.arange(0.0, JN, 1)

      # Bagian Sisi Kiri
      # Base Matrix Inverse
    D = np.identity(JN, dtype=float)
    D[0, 0] = 1
    D[JN - 1, JN - 1] = 1

    for i in range(JN - 2):
        D[i + 1, i] = 1
        D[i + 1, i + 1] = -2 - b
        D[i + 1, i + 2] = 1
    D = np.linalg.inv((D))

    #Bagian Sisi Kanan
    if BKK == 1:
        for i in range(JTS - 1):
            X[i + 1, 0] = PL
            X[i + 1, JN - 1] = PR

        for i in range(JTS - 1):
            A[0] = PL
            A[JN - 1] = PR
            for n in range( JN - 2):
                A[n + 1] = X[i, n + 1] * -b
            A = np.dot(D, A)
            for n in range(JN - 2):
                X[i + 1, n + 1] = A[n + 1]
            plt.plot(W, X[i,])
        plt.show()

    elif BKK == 2:

        for i in range(JTS - 1):
            X[i + 1, 0] = PL
            X[i + 1, JN - 1] = X[i + 1, 0] + ((2 * X[i + 1, 0] - 2*X[i + 1, 0] - PL*dx*2)/b)

        for i in range(JTS - 1):
            A[0] = PL
            A[JN - 1] = X[i + 1, 0] + ((2 * X[i + 1, 0] - 2*X[i + 1, 0] - PL*dx*2)/b)
            for n in range(JN - 2):
                A[n + 1] = X[i, n + 1] * -b
            A = np.dot(D, A)
            for n in range(JN - 2):
                X[i + 1, n + 1] = A[n + 1]
            plt.plot(W, X[i,])
        plt.show()

    elif BKK == 3:

        for i in range(JTS - 1):
            X[i + 1, 0] = X[i + 1, 0] + ((2 * X[i + 1, 0] - 2*X[i + 1, 0] - PL*dx*2)/b)
            X[i + 1, JN - 1] = PR
        for i in range(JTS - 1):
            A[0] = X[i + 1, 0] + ((2 * X[i + 1, 0] - 2*X[i + 1, 0] - PL*dx*2)/b)
            X[i + 1, JN - 1] = PR
            for n in range(JN - 2):
                A[n + 1] = X[i, n + 1] * -b
            A = np.dot(D, A)
            for n in range(JN - 2):
                X[i + 1, n + 1] = A[n + 1]
            plt.plot(W, X[i,])
        plt.show()

    elif BKK == 4:

        for i in range(JTS - 1):
            X[i + 1, 0] = X[i + 1, 0] + ((2 * X[i + 1, 0] - 2*X[i + 1, 0] - PL*dx*2)/b)
            X[i + 1, JN - 1] = X[i + 1, JN - 1] + (((2*X[i + 1, JN - 1] - 2*X[i + 1, JN - 1] + dx*2)/b))

        for i in range(JTS - 1):
            A[0] = X[i + 1, 0] + ((2 * X[i + 1, 0] - 2*X[i + 1, 0] - PL*dx*2)/b)
            A[JN - 1] = X[i + 1, JN - 1] + (((2*X[i + 1, JN - 1] - 2*X[i + 1, JN - 1] + dx*2)/b))
            for n in range(JN - 2):
                A[n + 1] = X[i, n + 1] * -b
            A = np.dot(D, A)
            for n in range(JN - 2):
                X[i + 1, n + 1] = A[n + 1]
            plt.plot(W, X[i,])
        plt.show()
