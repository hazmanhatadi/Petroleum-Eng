##################################################################
# Kelompok 19 #
# Hazman Hatadi - 101317004 #
# Hana Khayyirah H - 101317024 #
# Innocentia Ukkana - 101317070 #
# Intan Umull M. Sary - 101317104 #

from __future__ import division

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
from scipy.sparse.linalg import gmres
# np.set_printoptions(threshold=sys.maxsize)


##################################################################
#                       PENGELOMPOKAN CLASS                      #
##################################################################
class prop_rock(object):
# Class untuk menyimpan nilai dari properti batuan seperti permeabilitas, porositas, dan kompresibilitas #
    def __init__(self, kx=0, ky=0, por=0, cr=0):
        self.kx = kx
        self.ky = ky
        self.por = por
        self.cr = cr

class prop_fluid(object):
# Class untuk menyimpan nilai dari properti fluida seperti kompresibilitas, viskositas, densitas, dan FVF(b) yang merupakan fungsi tekanan #
    def __init__(self, c_o=0, mu_o=0, rho_o=0, b_o=0):
        self.c_o = c_o
        self.mu_o = mu_o
        self.rho_o = rho_o
        self.b_o = b_o

    def calc_b(self, p):
        return 1 / (1 + self.c_o * (p - 14.7))

class prop_grid(object):
# Class untuk menyimpan ukuran grid and posisi #
    def __init__(self, Nx=0, Ny=0, Nz=0):
        self.Nx = Nx
        self.Ny = Ny
        self.Nz = Nz

class prop_res(object):
# Class untuk menyimpan ukuran reservoir dan tekanan awal #
    def __init__(self, Lx=0, Ly=0, Lz=0, p_init=0):
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        self.p_init = p_init

class prop_well(object):
# Class untuk menyimpan grid sumur dengan laju alir tertentu dan konversi dari lokasi sumur ke grid #
    def __init__(self, loc=0, q=0):
        self.loc = loc
        self.q = q

    def index_to_grid(self, Nx):
        return self.loc[1] * Nx + self.loc[0]

class prop_time(object):
# Class untuk menyimpan nilai dari timestep dan interval waktu #
    def __init__(self, tstep=0, timeint=0):
        self.tstep = tstep
        self.timeint = timeint

############################################################
#                   INISIALISASI FUNGSI                    #
############################################################
def calc_transmissibility_x(k_x, mu, B_o, props, i, j):
# Menghitung transmisibilitas di arah x #
    dx = props['res'].Lx / props['grid'].Nx
    dy = props['res'].Ly / props['grid'].Ny
    dz = props['res'].Lz / props['grid'].Nz
    k_x = (k_x[j, i] + k_x[j, i + 1]) / 2
    mu = (mu[j, i] + mu[j, i + 1]) / 2
    B_o = (B_o[j, i] + B_o[j, i + 1]) / 2
    return k_x * dy * dz / (mu * B_o * dx)

def calc_transmissibility_y(k_y, mu, B_o, props, i, j):
# Menghitung transmisibilitas di arah y #
    dx = props['res'].Lx / props['grid'].Nx
    dy = props['res'].Ly / props['grid'].Ny
    dz = props['res'].Lz / props['grid'].Nz
    k_y = (k_y[j, i] + k_y[j + 1, i]) / 2
    mu = (mu[j, i] + mu[j + 1, i]) / 2
    B_o = (B_o[j, i] + B_o[j + 1, i]) / 2
    return k_y * dx * dz / (mu * B_o * dy)

def ij_to_grid(i, j, Nx):
# Mengubah nilai grid ke angka sebagai panjang reservoir #
    return (i) + Nx * j

def construct_T(mat, params, props):
# Membuat matrix transmisibilitas dengan memasukkan properti fluida ke grid #
    k_x = params['k_x']
    k_y = params['k_y']
    B_o = params['B_o']
    mu = params['mu']

    m = mat.shape[0]
    n = mat.shape[1]
    A = np.zeros((m*n , m*n))

####################
##### IMPLISIT #####
####################
    for j in range(m):              #looping transmisibilitas arah y
        for i in range(n):          #looping transmisibilitas arah x
            # 2 tetangga di arah x
            if i < n - 1:           #kalkulasi kondisional, menghitung dua kotak di arah x
                A[mat[j, i] - 1, mat[j, i + 1] - 1] = calc_transmissibility_x(k_x, mu, B_o, props, i, j)
                A[mat[j, i + 1] - 1, mat[j, i] - 1] = A[mat[j, i] - 1, mat[j, i + 1] - 1]

            # 2 tetangga di arah y
            if j < m - 1:           #kalkulasi kondisional, menghitung dua kotak di arah y
                A[mat[j, i] - 1, mat[j + 1, i] - 1] = calc_transmissibility_y(k_y, mu, B_o, props, i, j)
                A[mat[j + 1, i] - 1, mat[j, i] - 1] = A[mat[j, i] - 1, mat[j + 1, i] - 1]

    for k in range(A.shape[0]):
        p = np.sum(A[k, :]) *-1     #menjumlahkan nilai tiap baris dikali -1
        A[k, k] = p                 #mengembalikan nilai array selain point of interest ke nilai semula pada matriks
    return A /1000


def run_simulation(props):
    rock = props['rock']
    fluid = props['fluid']
    grid = props['grid']
    res = props['res']
    wells = props['well']
    sim_time = props['time']

    # Memasukkan nilai properties awal kedalam matrix
    k_x = np.full((grid.Ny, grid.Nx), rock.kx)              #membuat matrix sepanjang Ny, Nx dengan isi kx
    k_y = np.full((grid.Ny, grid.Nx), rock.ky)              #membuat matrix sepanjang Ny, Nx dengan isi ky
    B_o = np.full((grid.Ny, grid.Nx), fluid.calc_b(res.p_init))
    mu = np.full((grid.Ny, grid.Nx), fluid.mu_o)
    p_grids = np.full((grid.Ny * grid.Nx, 1), res.p_init)   #membuat matrix sepanjang Ny*Nx dan sumbu x=1
    params = {'k_x': k_x, 'k_y': k_y, 'B_o': B_o, 'mu': mu} #memasukkan parameter tersebut ke dictionary params

    # Membuat transmisibilitas matrix T
    mat = np.reshape(np.arange(1, grid.Ny * grid.Nx + 1), (grid.Ny, grid.Nx))   #membuat matrix berukuran Ny*Nx
    T = construct_T(mat, params, props)

    # Membuat matrix A = transmisibilitas matrix - akumulasi matrix
    dx = res.Lx / grid.Nx
    dy = res.Lx / grid.Nx
    dz = res.Lx / grid.Nx
    V = dx * dy * dz                                #Mencari Volume
    accumulation = V * rock.por * fluid.c_o / (5.615 * (sim_time.tstep))        #rumus akumulasi volume per timestep, dikonversi ke bbl
    A = T - np.eye(T.shape[0]) * accumulation       #menambah pengaruh akumulasi di tiap arah diagonal dari T

    #memberi pengaruh laju alir ke matrix
    Q = np.zeros((T.shape[0], 1))   #membuat matrix dengan sumbu Y sebesar Y dari matrix T, matrix dengan sumbu X sebesar 1
    for well in wells:
        Q[well.index_to_grid(grid.Nx)] = -well.q    #menambah adanya laju alir sehingga ada pengurangan volume di lokasi grid sumur

    # Kalkulasi Right Hand Side (RHS)
    p_RHS = np.full((grid.Ny * grid.Nx, 1), -accumulation * res.p_init) #membuat matrix ukuran Nx*Ny, dengan sumbu x=1 dengan isi -akumulasi*p_init
    b = p_RHS - Q   #mengurangi volume yang ada di matrix p_RHS dengan pengaruh laju alir, sehingga hanya cell terakhir yang berkurang

    # Time-loop
    for t in sim_time.timeint:
        print('evaluating t = %1.1f (days)' % t)

        # Menghitung tekanan di waktu n+1
        p_grids = (gmres(A, b))[0]  #mencari nilai p_grids dengan melakukan penyelesaian x = gmres(A,b), dimana A*x = b for x
        p_grids = np.reshape(p_grids, (len(p_grids), 1))    #melakukan penyusunan kembali array p_grids dengan x=1

        # Update Bo, b, dan transmisibilitas matrix
        for i in range(grid.Nx):
            for j in range(grid.Ny):
                B_o[i, j] = fluid.calc_b(p_grids[ij_to_grid(i, j, grid.Nx)])
        params['B_o'] = B_o
        A = construct_T(mat, params, props)
        A = A - np.eye(A.shape[0]) * accumulation

        b = -accumulation * p_grids - Q
    return p_grids

##########################################################
#                PROGRAM INPUT DATA                      #
##########################################################
def main():
    # INISIALISASI DATA
    tstep = 1
    endtime = int(input('Masukkan jumlah hari simulasi (hari) ='))
    timeint = np.arange(0, endtime, tstep)
    sim_time = prop_time(tstep=tstep,
                         timeint=timeint)
    rock = prop_rock(kx=int(input('Masukkan nilai permeabilitas arah x (mD) =')),
                     ky=int(input('Masukkan nilai permeabilitas arah y (mD) =')),
                     por=float(input('Masukkan nilai porositas (fraksi) =')),
                     cr=float(input('Masukkan nilai kompresibilitas batuan (1/Psi) =')) )
    fluid = prop_fluid(c_o=float(input('Masukkan nilai kompresibilitas minyak (1/Psi) =')),
                       mu_o=float(input('Masukkan nilai viskositas minyak (cP) =')),
                       rho_o=float(input('Masukkan nilai densitas minyak (lbm/ft3) =')) )
    grid = prop_grid(Nx=int(input('Masukkan jumlah grid arah x =')),
                     Ny=int(input('Masukkan jumlah grid arah y =')),
                     Nz=int(input('Masukkan jumlah grid arah z =')))
    res = prop_res(Lx=int(input('Masukkan panjang reservoir (ft) =')),
                   Ly=int(input('Masukkan lebar reservoir (ft) =')),
                   Lz=int(input('Masukkan tebal reservoir (ft) =')),
                   p_init=int(input('Masukkan nilai tekanan awal (Psi) =')))
    prod_x = int(input('Masukkan koordinat sumur produksi arah x (harus lebih kecil dari jumlah grid arah x) ='))
    prod_y = int(input('Masukkan koordinat sumur produksi arah y (harus lebih kecil dari jumlah grid arah y) ='))
    well_prod = prop_well(loc=(prod_y, prod_x),
                      q=int(input('Masukkan nilai laju alir sumur (STB/D) =')))

    props = {'rock': rock, 'fluid': fluid, 'grid': grid, 'res': res, 'well': [well_prod], 'time': sim_time}
    if well_prod.loc[0] > grid.Ny or well_prod.loc[1]>grid.Nx:
        print('Masukkan data lokasi sumur dengan benar!')
        sys.exit()


## Case: 1 producer 1 injector
    # Define 1 injector
    inj_x = int(input('Masukkan koordinat sumur injeksi arah x (harus lebih kecil dari jumlah grid arah x) ='))
    inj_y = int(input('Masukkan koordinat sumur injeksi arah y (harus lebih kecil dari jumlah grid arah y) ='))
    well_inj = prop_well(loc=(inj_y, inj_x), q=-int(input('Masukkan nilai laju alir sumur (STB/D) =')))
    props['well'] = [well_prod, well_inj]

    # Run simulasi
    print('Running Case: 1 producer 1 injector')
    p_grids = run_simulation(props)

    # Plot
    p_2D = np.reshape(p_grids, (grid.Nx, grid.Ny))
    plt.show()
    plt.plot(p_2D)

    plt.matshow(p_2D)
    plt.colorbar()
    plt.xlabel('grid arah x')
    plt.ylabel('grid arah y')
    plt.draw()
    plt.title('CASE = 1 PROD + 1 INJ')
    plt.show()
    plt.show(block=True)

#######################################################################################

if __name__ == '__main__':
    main()
