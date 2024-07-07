import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

N = 20
L = 10.0
h = 0.002
epsilon = 1
sigma = 1
m = 1
v_0 = 1
k=1

Time = 50
time = np.arange(0, Time, h)
f = open(f"pruebas/datos_{v_0}.txt", "w")
f_energia = open(f"pruebas/energia_{v_0}.txt", "w")
f_velocidad = open(f"pruebas/velocidad_{v_0}.txt", "w")

def init_cond():
    r = np.zeros((N, 2))
    v =np.zeros((N, 2))
    for i in range(N):
        r[i] = np.array([i%(L/2)*2+1, i%4*2 +1 ]) + 1* np.random.rand(1, 2)[0]
        # r[i] = np.array([i%(L/2)*2+1, i%4*2.5 +1 ])
        theta = np.random.rand()*2*np.pi
        v[i] = v_0 * np.array([np.sin(theta), np.cos(theta)])
    return r, v

def cond_contorno(r):
    if(r[0] > L):
        r[0] = r[0] - L
    if(r[0] < 0):
        r[0] = r[0] + L
    if(r[1] > L):
        r[1] = r[1] - L
    if(r[1] < 0):
        r[1] = r[1] + L
    return r

def cond_contorno_distancia(r):
    if(r[0] > L/2):
        r[0] = r[0] - L
    elif(r[0] < -L/2):
        r[0] = r[0] + L
    if(r[1] > L/2):
        r[1] = r[1] - L
    elif(r[1] < -L/2):
        r[1] = r[1] + L
    return r

def lennard_jones(r):
    R = compute_distance(r)
    acc = np.zeros((N, 2))
    for i in range(N):
        for j in range(N):
            if(i!=j):
                norm  = np.linalg.norm(R[i, j])
                if (norm < 3):
                    acc[i] = acc[i] + 4*R[i, j]* epsilon * (6*np.power((sigma/norm), 5) - 12*np.power((sigma/norm), 11))/(norm*m)
    return acc, R

def compute_distance(r):
    R = np.zeros((N, N, 2))
    for i in range(0, N-1):
        for j in range(i+1, N):
            R[i, j] = r[j]- r[i]
            R[i, j] = cond_contorno_distancia(R[i, j])
            R[j, i] = -R[i, j]
    return R

def verlet_algorithm(r, v, a):
    w = np.zeros((N, 2))
    for i in range(N):
        r[i] = r[i] + h*v[i] + h*h*a[i]/2
        r[i] = cond_contorno(r[i])
        w[i] = v[i] + h*a[i]/2
    a, R = lennard_jones(r)
    for i in range(N):
        v[i] = w[i] + h*a[i]/2
    return r, v, a, R

def compute_energy(v, R):
    T = 0
    V = 0
    for i in range(N):
        T = T + 0.5*m*np.linalg.norm(v[i])**2
        for j in range(N):
            if(i!=j):
                norm  = np.linalg.norm(R[i, j])
                V = V + 4*epsilon * (np.power((sigma/norm), 12) - np.power((sigma/norm), 6))
    return T, V

def compute_average_speed(v):
    v_prom = 0
    for i in range(N):
        v_prom = v_prom + np.linalg.norm(v[i])
    return v_prom/N

def write_vector(r, f):
    for i in range(N):
        f.write(f"{r[i, 0]}, {r[i, 1]}\n")
    f.write(f"\n")

def write_velocity(v, f):
    for i in range(N):
        f.write(f"{v[i, 0]}, {v[i, 1]}, {np.linalg.norm(v[i])}\n")

T, V = np.zeros(len(time)), np.zeros(len(time))

r, v = init_cond()
v_0 = v
write_vector(r, f)
a, R = lennard_jones(r)
T[0], V[0] = compute_energy(v, R)
f_energia.write(f"{T[0]}, {V[0]}\n")

for i in range(1, len(time)):
    r, v, a, R = verlet_algorithm(r, v, a)
    T[i], V[i] = compute_energy(v, R)
    write_vector(r, f)
    if(i >= int(20/h)):
        write_velocity(v, f_velocidad)
    f_energia.write(f"{T[i]}, {V[i]}\n")

f.close()
f_velocidad.close()
f_energia.close()


velocidad = pd.read_csv(f'pruebas/velocidad_{v_0}.txt', delimiter=',', index_col=False, header=0 ,names=['vx', 'vy', 'v'])

Temp = m/(2*k)*np.mean(velocidad.v**2)

def distribution_module (v):
    return (m/(k*Temp))*v*np.exp(-m*v*v/(2*k*Temp))

def distribution_velocity (v):
    return np.sqrt(2*m/(np.pi*k*Temp))*np.exp(-m*v*v/(2*k*Temp))

fig=plt.figure(figsize=(10,7)) 
ax=fig.add_subplot(111)
v = np.linspace(0, np.max(velocidad.v), 100)
plt.grid(zorder=0)
plt.hist(velocidad.v, bins=100 , rwidth=0.7, zorder = 3, density=True, label='Histograma de velocidades')
plt.plot(v, distribution_module(v), label=f'Distribución de velocidades ', zorder = 4)
plt.ylabel('nº de cuenta ')
plt.xlabel('velocidad')
plt.legend()
plt.title(f'Histograma del módulo de la velocidad, con temperatura T = {Temp:.2f}')
plt.savefig(f"latex/plots/velocidad_{v_0}.png")

fig=plt.figure(figsize=(10,7)) 
ax=fig.add_subplot(111)
vx = np.linspace(0, np.max(velocidad.vx), 100)
plt.grid(zorder=0)
plt.hist(np.abs(velocidad.vx), bins=100 , rwidth=0.7, zorder = 3, density=True, label='Histograma de velocidades')
plt.plot(vx, distribution_velocity(vx), label=f'Distribución de velocidades ', zorder = 4)
plt.ylabel('nº de cuenta ')
plt.xlabel('velocidad')
plt.legend()
plt.title(f'Histograma del valor absoluto de la velocidad en el eje x, con temperatura T = {Temp:.2f}')
plt.savefig(f"latex/plots/velocidad_x_{v_0}.png")

fig=plt.figure(figsize=(10,7)) 
ax=fig.add_subplot(111)
vy = np.linspace(0, np.max(velocidad.vy), 100)
plt.grid(zorder=0)
plt.hist(np.abs(velocidad.vy), bins=100 , rwidth=0.7, zorder = 3, density=True, label='Histograma de velocidades')
plt.plot(vy, distribution_velocity(vy), label=f'Distribución de velocidades', zorder = 4)
plt.ylabel('nº de cuenta ')
plt.xlabel('velocidad')
plt.legend()
plt.title(f'Histograma del valor absoluto de la velocidad en el eje y, con temperatura T = {Temp:.2f}')
plt.savefig(f"latex/plots/velocidad_y_{v_0}.png")


fig=plt.figure(figsize=(10,7)) 
ax=fig.add_subplot(111)
v = np.linspace(0, np.max(velocidad.v), 100)
plt.grid(zorder=0)
plt.hist(velocidad.v, bins=100 , rwidth=0.7, zorder = 3, density=True, label='Histograma para el modulo de la velocidad', alpha=0.7, color='tab:blue')
plt.hist(np.abs(velocidad.vx), bins=100 , rwidth=0.7, zorder = 4, density=True, label='Histograma para el eje x', alpha=0.7, color='tab:purple')
plt.hist(np.abs(velocidad.vy), bins=100 , rwidth=0.7, zorder = 5, density=True, label='Histograma para el eje y', alpha=0.7, color='tab:red')

plt.plot(v, distribution_module(v), label=f'Distribución del módulo de velocidades ', zorder = 7, color='tab:blue')
plt.plot(v, distribution_velocity(v), label=f'Distribución de la componente x', zorder = 8, color='tab:purple')
plt.plot(v, distribution_velocity(v), label=f'Distribución de la componente y', zorder = 9, color='tab:red')

plt.ylabel('nº de cuenta (normalizado)')
plt.xlabel('velocidad')
plt.legend()
plt.title(f'Histograma del valor absoluto de la velocidad, con temperatura T = {Temp:.2f}')
plt.savefig(f"latex/plots/velocidad_todos_{v_0}.png")

plt.show()
