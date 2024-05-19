import multiprocessing
import time
import numpy as np
import matplotlib.pyplot as plt



def initialize_grid(size, initial):
    """
    Inicializa uma rede cúbica com um determinado tamanho.
    
    size - número de elementos numa aresta do cubo
    initial - valor inicial de todos os elementos na rede: 1 ou -1.
    Se o atributo for diferente desses dois valores, a rede será
    gerada aleatoriamente.
    
    return: rede cúbica

    """
    
    if initial == -1:
        return np.full((size, size, size), -1)
    if initial == 1:
        return np.full((size, size, size), 1)
    else:
        return np.full((size, size, size), np.random.choice([-1,1], size = (size, size, size)))

               
def transitionvalues(t, h):
    """
    Gera uma tabela com todos os valores da função de transição possíveis.
    
    t - temperatura reduzida
    h - campo magnético externo reduzido
    
    return: tabela com os valores da função transição.
    
    """
    possible_values_sum = np.arange(-6, 8, 2)
    possible_spins = np.arange(-1, 3, 2)
    
    num = [[j + i * h for i in possible_spins] for j in possible_values_sum]
    transvalue = [[1 if elem <= 0 else np.exp( - elem / t) for elem in row] for row in num]
   
    return np.array(transvalue)
    
    
             
def neighborstable(size):
    """
    Cria uma tabela com os indices dos vizinhos mais próximos
    para uma rede.
    
    size - aresta da rede cúbica
    
    returns: tabela com os indices do vizinhos mais próximos na
    direção positiva e direção negativa.
    
    """
    
    
    positive = [i + 1 for i in range(size)] #Vizinho na direção positiva
    positive[-1] = 0
    negative = [i - 1 for i in range(size)] #Vizinho na direção negativa
    
       
    return np.array([positive, negative])


def w(sigma, sigmasum, transvalue):
    """
    Permite localizar o valor da função de transição para um elemento
    com um dado sigma e soma dos sigmas dos vizinhos.
    
    sigma - valor do spin do elemento da rede (-1 ou 1)
    sigasum - valor da soma dos spins dos vizinhos mais próximos
    nas 3 direções.
    transvalue - tabela com todos os valores possiveís para a funcão
    de transição
    
    returns: valor de função de transição
    
    """
    
    
    
    i = int(sigmasum / 2 + 3)
    j = int(sigma / 2 + 1 / 2)
    
    return transvalue[i, j]



 
def calc_sus(size, moment, t):
    """
    Calcula o valor da susceptibilidade magnética da rede cúbica.
    
    size - aresta da rede
    moment - momento magnético
    t- temperatura reduzida
    
    returns: susceptibilidade magnética.
    
    """
    
    return moment.var() * size ** 3 / t
    
   

def calc_cap(size, e, t):
    """
    Calcula o valor da capacidade térmica da rede cúbica.
    
    size- aresta da rede
    e - energia
    t - temperatura reduzida
    
    """

    return e.var() / (size ** 3 * t ** 2)
    
def cycle(grid, size, neighbors, transvalue, h):
    """
    Corre um ciclo de Monte Carlo para a simulação de ferromagnetismo.
    
    grid - rede cúbica
    size - aresta da rede
    neighbors - tabela com os indices dos vizinhos mais próximos
    transvalue - valores possíveis de função de transição
    h - campo magnético externo reduzido
    
    returns:
    grid - rede cúbica após um ciclo
    energy - energia total da rede após um ciclo
    
    """
    
    energy = 0
    
    for i in range(size):
        for j in range(size):
            for k in range(size):
            
                sigma = grid[i, j, k] 
            
                sigmasum = grid[neighbors[0, i], j, k] + grid[neighbors[1, i], j, k] + grid[i, neighbors[0, j], k] + grid[i, neighbors[1, j], k] + grid[i, j, neighbors[0, k]] + grid[i, j, neighbors[1, k]]
                sigmasum = sigmasum * sigma
               
                
                
                energy_point = -0.5 * sigmasum - sigma * h #energia por ponto da rede
                
               
                p = np.random.random() #aleatório entre 0 e 1
            
                if p < w(sigma, sigmasum, transvalue): #inversão de spin
                    grid[i, j, k] = -sigma 
                    energy_point = -energy_point
                    
                energy += energy_point
                    
    
   
    return grid, energy

def simulation(size, initial, cycles, t, h):
    """
    Corre uma simulação de ferromagnetismo com um dado número de
    ciclos de Monte Carlo
    
    size - aresta da rede
    initial - variável de inicialização da rede
    cycles - número de ciclos na simulação
    t - temperatura reduzida
    h - campo magnético externo reduzido
    
    returns:
    grid - rede cúbica após a simulação
    moment - lista com os momentos magnéticos da rede após cada ciclo
    total_energy - lista com as energias totais da rede após cada ciclo
    
    """
    
    grid = initialize_grid(size, initial)
    
    transvalue = transitionvalues(t, h)
    
    neighbors = neighborstable(size)
    
   
    moment = np.zeros(cycles)
    total_energy = np.zeros(cycles)
    
    
    for i in range(cycles):
        
        grid, energy = cycle(grid, size, neighbors, transvalue, h)
        
       
        total_energy[i] = energy
        moment[i] = 2 * grid[grid == 1].shape[0] - size ** 3 #momento magnético
        
        
    
    moment /= size**3
    total_energy /= size **3
    
    
                
    return grid, moment, total_energy


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
    
    plt.show()

    
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
    
    plt.show()


def hyst_graph(m, hs, temps):
    
    plt.figure(figsize=(12, 12))
    for i, temp in enumerate(temps):
        plt.plot(hs, m[i], label=f'Temperature: {temp}')
    plt.title('m vs h with different temperatures')
    plt.xlabel('h')
    plt.ylabel('m')
    plt.legend()
    plt.show()

    
def tails_graph(m, hs, initials):
    plt.figure(figsize=(12, 12))
    for i, s in enumerate(initials):
        plt.plot(hs, m[i], label=f'Spin: {s}')
    plt.title('m vs h with different initial conditions')
    plt.xlabel('h')
    plt.ylabel('m')
    plt.legend()
    plt.show()
    
    
    
def simulacao_temp(temps, size, initial, cycles, h):
    """
    Corre simulação de ferromagnetismo com um campo magnétco externo fixo
    e para várias temperaturas reduzidas
    
    temps - lista com as temperaturas reduzidas a usar
    size - aresta da rede
    initial - variável de inicialização
    cycles - número de ciclos
    h - campo magnético externo reduzido
    
    returns:
    m - lista com momentos magnéticos médios após cada simulação
    sus - lista com susceptibilidades magnéticas médias após cada simulação
    e - lista com energias totais médias após cada simulação
    c - lista com capacidades térmicas médias após cada simulação
    """
    
    n_pontos = temps.size
    
    # armazenar valores
    m = np.zeros(n_pontos)
    sus = np.zeros(n_pontos)
    e = np.zeros(n_pontos)
    c = np.zeros(n_pontos)
    
    start_indx = int(cycles / 100)
    
    indx = 0
    for t in temps:
        
        _, moment, en = simulation(size, initial, cycles, t, h)
        
        moment = np.abs(moment[start_indx:])
        en = en[start_indx:]
        m[indx] = moment.mean()
        sus[indx] = calc_sus(size, moment, t)
        e[indx] = en.mean()
        c[indx] = calc_cap(size, en, t)
        indx += 1
    
    return m, sus, e, c

def simul_h(hs, size, initial, cycles, t):
    """
    Corre simulação de ferromagnetismo com uma temperatura reduzida fixa
    e para vários campos magnéticos externos reduzidos
    
    hss - lista com os campos magnéticos reduzidos a usar
    size - aresta da rede
    initial - variável de inicialização
    cycles - número de ciclos
    t - temperatura reduzida
    
    returns:
    m - lista com momentos magnéticos médios após cada simulação
    sus - lista com susceptibilidades magnéticas médias após cada simulação
    e - lista com energias totais médias após cada simulação
    c - lista com capacidades térmicas médias após cada simulação
    """
    
    n_pontos = hs.size
    
    # armazenar valores
    m = np.zeros(n_pontos)
    sus = np.zeros(n_pontos)
    e = np.zeros(n_pontos)
    c = np.zeros(n_pontos)
   
    start_indx = int(cycles / 10)
   
    indx = 0
    for h in hs:
       
        _, moment, en = simulation(size, initial, cycles, t, h)
       
        moment = moment[start_indx:]
        en = en[start_indx:]
        m[indx] = moment.mean()
        sus[indx] = calc_sus(size, moment, t)
        e[indx] = en.mean()
        c[indx] = calc_cap(size, en, t)
        indx += 1
   
    return m, sus, e, c
    


def hysteresis(hs, temps, size, initial, cycles):
    """
    Corre a simulação de ferromagnetismo para vários valores de campo magnético
    reduzido e temperatura reduzida.
    
    hs - lista com os valores de campo magnético reduzido a usar
    temps - lista com os valores de temperatura reduzida a usar
    size - aresta da rede
    initial - variável de inicialização
    cycles - número de ciclos
    
    returns: tabela com os valores de momento magnético em função do campo
    magnético externo reduzido para várias temperaturas reduzidas
    """
    
    n_hs = hs.size
    n_temps = temps.size
    
   
    m = np.zeros((n_temps, n_hs)) # armazenar valores
    
    
    start_indx = int(cycles / 10)
    
    for i, t in enumerate(temps):
        for j, h in enumerate(hs):
            _, moment, _ = simulation(size, initial, cycles, t, h)
            moment = moment[start_indx:]
            m[i, j] = moment.mean()
           
    
    return m

def tails(hs, temp, size, initials, cycles):
    """
    Corre a simulação de ferromagnetismo com temperatura reduzida fixa para
    ambas as condições de inicialização e rede.
    
    hs - valores de campo magnético reduzido a usar
    temp - temperatura reduzida
    size - aresta da rede
    initials - valores de inicialização de rede (-1 e 1)
    cycles - número de ciclos
    
    returns: tabela com os valores de momento magnético em função do campo
    magnético externo reduzido para ambas as condições de inicialização.
    
    """
    
    n_hs = hs.size
    n_initials = initials.size
    
    
    m = np.zeros((n_initials, n_hs)) # armazenar valores
    
    
    start_indx = int(cycles / 10)
    
    for i, s in enumerate(initials):
        for j, h in enumerate(hs):
            _, moment, _ = simulation(size, s, cycles, temp, h)
            moment = moment[start_indx:]
            m[i, j] = moment.mean()
            
    return m


def curieT(sus, temps):
    """
    Determina a temperatura de Curie
    
    sus - lista com susceptibilidades magnéticas médias
    temps - lista de temperaturas correspondentes às susceptibilidades magnéticas
    
    returns: Temperatura de Curie
    
    """
        
    
    
    maxsus = max(sus)
    curieindex = np.where(sus == maxsus)
    curieT = temps[curieindex]
    return curieT

# Define magnetic field strengths and temperatures
hs = np.arange(-5, 5.5, .5)
temps = np.arange(0.5, 5.5, .5)  
initials = np.arange(-1, 3, 2)
cycles = 10
size = 10

#Define Functions
def fun1 ():
   print("Runnning f1")
   m0 = hysteresis(hs, temps, size, 2, cycles) 
   hyst_graph(m0, hs, temps)
   print("Finished f1")

def fun2():
   print("Runnning f2")
   m, sus, e, c = simulacao_temp(temps, size, 3, cycles, 1)
   ferro_graft(m, sus, e, c, temps)
   curie = curieT(sus, temps)
   print(curie)
   print("Finished f2")


def fun3():
   print("Runnning f3")
   m1, sus1, e1, c1 = simul_h(hs, size, 2, cycles,1)
   ferro_grafh(m1, sus1, e1, c1, hs)
   print("Finished f3")

 

def fun4():
   print("Runnning f4")
   m2 = tails(hs, 1, size, initials, cycles)
   tails_graph(m2, hs, initials)
   print("Finished f4")



if __name__ == '__main__':
   
   start_time = time.time()
   # Create Processes
   process1 = multiprocessing.Process(target=fun1)
   process2 = multiprocessing.Process(target=fun2)
   process3 = multiprocessing.Process(target=fun3)
   process4 = multiprocessing.Process(target=fun4)

   # Start Processes
   process1.start()
   process2.start()
   process3.start()
   process4.start()

   # Join Processes
   process1.join()
   process2.join()
   process3.join()
   process4.join()

   end_time = time.time()
   print("Total execution time with multiprocessing: {:.2f} seconds".format(end_time - start_time))

