import numpy as np
import matplotlib.pyplot as plt



def initialize_grid(size, initial):
    
    if initial == -1:
        return np.full((size, size, size), -1)
    if initial == 1:
        return np.full((size, size, size), 1)
    else:
        return np.full((size, size, size), np.random.choice([-1,1], size = (size, size, size)))

               
def transitionvalues(t, h):
    
    possible_values_sum = np.arange(-6, 8, 2)
    possible_spins = np.arange(-1, 3, 2)
    
    num = [[j + i * h for i in possible_spins] for j in possible_values_sum]
    transvalue = [[1 if elem <= 0 else np.exp( - elem / t) for elem in row] for row in num]
   
    return np.array(transvalue)
    
    
             
def neighborstable(size):
    
    
    positive = [i + 1 for i in range(size)]
    positive[-1] = 0
    negative = [i - 1 for i in range(size)]
    
       
    return np.array([positive, negative])


def w(sigma, sigmasum, transvalue):
    
    
    i = int(sigmasum / 2 + 3)
    j = int(sigma / 2 + 1 / 2)
    
    return transvalue[i, j]



 
def calc_sus(size, order, t):
    """
    ...
    """
    
    return order.var() * size ** 3 / t
    
   

def calc_cap(size, e, t):
    """
    ...
    """

    return e.var() / (size ** 3 * t ** 2)
    
def cycle(grid, size, neighbors, transvalue, h):
    
    energy = 0
    moment = 0
    for i in range(size):
        for j in range(size):
            for k in range(size):
            
                sigma = grid[i, j, k] 
            
                sigmasum = grid[neighbors[0, i], j, k] + grid[neighbors[1, i], j, k] + grid[i, neighbors[0, j], k] + grid[i, neighbors[1, j], k] + grid[i, j, neighbors[0, k]] + grid[i, j, neighbors[1, k]]
                sigmasum = sigmasum * sigma
                moment += abs(sigmasum)
                
                
                energy_point = -0.5 * sigmasum - sigma * h 
                
               
                p = np.random.random()
            
                if p < w(sigma, sigmasum, transvalue):
                    grid[i, j, k] = -sigma 
                    energy_point = -energy_point
                    
                energy += energy_point
                    
    
   
    return grid, energy

def simulation(size, initial, cycles, t, h):
    
    
    grid = initialize_grid(size, initial)
    
    transvalue = transitionvalues(t, h)
    
    neighbors = neighborstable(size)
    
   
    order = np.zeros(cycles)
    total_energy = np.zeros(cycles)
    
    
    for i in range(cycles):
        
        grelha, energy = cycle(grid, size, neighbors, transvalue, h)
        
       
        total_energy[i] = energy
        order[i] = 2 * grid[grid == 1].shape[0] - size ** 3
        
    
    order /= size**3
    total_energy /= size **3
    
    
                
    return grid, order, total_energy


def ferro_graft(m, sus, e, c, t):
    """
    ...
    """
    
    fig = plt.figure()
    axm = fig.add_subplot(2, 2, 1)
    axm.plot(t, m)
    axm.set_title('m vs t')
    axm.set_ylabel(r'm')
    axm.set_xlabel(r't')
    
    axs = fig.add_subplot(2, 2, 2)
    axs.plot(t, sus)
    axs.set_title('$\chi$  vs  t')
    axs.set_ylabel(r'$\chi$')
    axs.set_xlabel(r't')
    
    axe = fig.add_subplot(2, 2, 3)
    axe.plot(t, e)
    axe.set_title('e vs t')
    axe.set_ylabel(r'e')
    axe.set_xlabel(r't')
    
    axc = fig.add_subplot(2, 2, 4)
    axc.plot(t, c)
    axc.set_title('C vs t')
    axc.set_ylabel(r'C')
    axc.set_xlabel(r't')
    axc.ticklabel_format(useMathText=True, style='sci')
    
    fig.set_size_inches(12, 12)
    
    
    
    
    return fig

def ferro_grafh(m, sus, e, c, h):
    """
    ...
    """
    
    fig = plt.figure()
    axm = fig.add_subplot(2, 2, 1)
    axm.plot(h, m)
    axm.set_title('m vs h')
    axm.set_ylabel(r'm')
    axm.set_xlabel(r'h')
    
    axs = fig.add_subplot(2, 2, 2)
    axs.plot(h, sus)
    axs.set_title('$\chi$  vs  h')
    axs.set_ylabel(r'$\chi$')
    axs.set_xlabel(r'h')
    
    axe = fig.add_subplot(2, 2, 3)
    axe.plot(h, e)
    axe.set_title('e vs h')
    axe.set_ylabel(r'e')
    axe.set_xlabel(r'h')
    
    axc = fig.add_subplot(2, 2, 4)
    axc.plot(h, c)
    axc.set_title('C vs h')
    axc.set_ylabel(r'C')
    axc.set_xlabel(r'h')
    axc.ticklabel_format(useMathText=True, style='sci')
    
    fig.set_size_inches(12, 12)
    
    
    
    
    return fig


def hyst_graph(m, hs, temps):
    
    plt.figure(figsize=(12, 12))
    for i, temp in enumerate(temps):
        plt.plot(hs, m[i], label=f'Temperature: {temp}')
    plt.title('m vs h with different temperatures')
    plt.xlabel('h')
    plt.ylabel('m')
    plt.legend()
    
def tails_graph(m, hs, initials):
    plt.figure(figsize=(12, 12))
    for i, s in enumerate(initials):
        plt.plot(hs, m[i], label=f'Spin: {s}')
    plt.title('m vs h with different initial conditions')
    plt.xlabel('h')
    plt.ylabel('m')
    plt.legend()
    
    
    
def simulacao_temp(temps, size, initial, nmax, h):
    """
    ...
    """
    
    n_pontos = temps.size
    # Inicializa os vetores para recolher os valores
    m = np.zeros(n_pontos)
    sus = np.zeros(n_pontos)
    e = np.zeros(n_pontos)
    c = np.zeros(n_pontos)
    
    start_indx = int(nmax / 100)
    
    indx = 0
    for t in temps:
        
        _, order, en = simulation(size, initial, nmax, t, h)
        
        order = np.abs(order[start_indx:])
        en = en[start_indx:]
        m[indx] = order.mean()
        sus[indx] = calc_sus(size, order, t)
        e[indx] = en.mean()
        c[indx] = calc_cap(size, en, t)
        indx += 1
    
    return m, sus, e, c

def simul_h(hs, size, initial, cycles, t):
    
   n_pontos = hs.size
   # Inicializa os vetores para recolher os valores
   m = np.zeros(n_pontos)
   sus = np.zeros(n_pontos)
   e = np.zeros(n_pontos)
   c = np.zeros(n_pontos)
   
   start_indx = int(cycles / 10)
   
   indx = 0
   for h in hs:
       
       _, order, en = simulation(size, initial, cycles, t, h)
       
       order = order[start_indx:]
       en = en[start_indx:]
       m[indx] = order.mean()
       sus[indx] = calc_sus(size, order, t)
       e[indx] = en.mean()
       c[indx] = calc_cap(size, en, t)
       indx += 1
   
   return m, sus, e, c
    


def hysteresis(hs, temps, size, initial, cycles):
    
    n_hs = hs.size
    n_temps = temps.size
    
    # Initialize arrays to store values
    m = np.zeros((n_temps, n_hs))
    
    
    start_indx = int(cycles / 10)
    
    for i, t in enumerate(temps):
        for j, h in enumerate(hs):
            _, order, _ = simulation(size, initial, cycles, t, h)
            order = order[start_indx:]
            m[i, j] = order.mean()
           
    
    return m

def tails(hs, temp, size, initials, cycles):
    
    n_hs = hs.size
    n_initials = initials.size
    
    # Initialize arrays to store values
    m = np.zeros((n_initials, n_hs))
    
    
    start_indx = int(cycles / 10)
    
    for i, s in enumerate(initials):
        for j, h in enumerate(hs):
            _, order, _ = simulation(size, s, cycles, temp, h)
            order = order[start_indx:]
            m[i, j] = order.mean()
            
    return m
# Define magnetic field strengths and temperatures
hs = np.arange(-5, 5.5, .5)
temps = np.arange(0.5, 5.5, .5)  
initials = np.arange(-1, 3, 2)

m0 = hysteresis(hs, temps, 10, 2, 100) 
fig0 = hyst_graph(m0, hs, temps)


m, sus, e, c = simulacao_temp(temps, 10, 3, 100, 1)
fig = ferro_graft(m, sus, e, c, temps)


m1, sus1, e1, c1 = simul_h(hs, 10, 2, 100,1)
fig1 = ferro_grafh(m1, sus1, e1, c1, hs)

m2 = tails(hs, 1, 10, initials, 100)
fig2 = tails_graph(m2, hs, initials)
plt.show()







   
  


  
