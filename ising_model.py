import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# Function to initialize a cubic grid
def initialize_grid(N):
    return np.full((N, N, N), 1)

# Function to plot the cubic grid
def plot_grid(grid):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    N = grid.shape[0]
    x, y, z = np.indices((N, N, N))  # Create 3D grid indices

    # Plot points for spins with value 1
    ax.scatter(x[grid == 1], y[grid == 1], z[grid == 1], color='r', marker='o', label='Spin Up')
    # Plot points for spins with value -1
    ax.scatter(x[grid == -1], y[grid == -1], z[grid == -1], color='b', marker='o', label='Spin Down')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Cubic Grid Visualization')
    ax.legend()

    plt.show()


# Function to initialize spins randomly
def initialize_spins(N):
    spins = np.random.choice([-1, 1], size=(N, N, N))
    return spins

# Function to calculate energy of a spin configuration
def calculate_energy(spins, N, J, H):
    energy = 0
    for i in range(N):
        for j in range(N):
            for k in range(N):
                energy += -J * spins[i, j, k] * (spins[(i+1)%N, j, k] + spins[(i-1)%N, j, k] +
                                                  spins[i, (j+1)%N, k] + spins[i, (j-1)%N, k] +
                                                  spins[i, j, (k+1)%N] + spins[i, j, (k-1)%N]) - H * spins[i, j, k]
    return energy


def compute_energy_change(spins, i, j, k, N, J, H):
    current_spin = spins[i, j, k]
    
    # Neighboring spins
    neighbors_sum = spins[(i+1)%N, j, k] + spins[(i-1)%N, j, k] + \
                    spins[i, (j+1)%N, k] + spins[i, (j-1)%N, k] + \
                    spins[i, j, (k+1)%N] + spins[i, j, (k-1)%N]
    
    # Contribution to energy from neighboring spins
    energy_from_neighbors = -J * current_spin * neighbors_sum
    
    # Contribution to energy from external magnetic field
    energy_from_field = -H * current_spin
    
    # Total energy change
    delta_E = 2 * (energy_from_neighbors + energy_from_field)
    
    return delta_E

# Function to perform a Monte Carlo step
def monte_carlo_step(spins, J, H, T):
    N = spins.shape[0]
    i = np.random.randint(N)
    j = np.random.randint(N)
    k = np.random.randint(N)
    dE = 2 * J * spins[i, j, k] * (spins[(i+1)%N, j, k] + spins[(i-1)%N, j, k] +
                                   spins[i, (j+1)%N, k] + spins[i, (j-1)%N, k] +
                                   spins[i, j, (k+1)%N] + spins[i, j, (k-1)%N]) + 2 * H * spins[i, j, k]
    if dE <= 0 or np.random.rand() < np.exp(-dE / T):
        spins[i, j, k] *= -1
    return spins

# Function to run the simulation
def simulate(N, J, H, T, steps):
    spins = initialize_grid(N)
    energies = np.zeros(steps)
    for step in range(steps):
        spins = monte_carlo_step(spins, J, H, T)
        energies[step] = calculate_energy(spins, J, H)
    return spins, energies

# Main function
def main():
    initial_grid = initialize_grid(10)
    # plot_grid(initial_grid)
    N = 10  # Size of the grid
    J = 1   # Interaction strength
    H = 0   # External magnetic field
    T = 1   # Temperature
    steps = 1000  # Number of Monte Carlo steps

    # spins, energies = simulate(N, J, H, T, steps)

    print("Energy of the spin configuration of 3 , 3 , 3: ", compute_energy_change(initial_grid, 3, 3, 3, N, J, H))

    # Plot energy vs. step
    # plt.plot(energies)
    # plt.xlabel('Step')
    # plt.ylabel('Energy')
    # plt.title('Energy vs. Step')
    # plt.show()

if __name__ == "__main__":
    main()
