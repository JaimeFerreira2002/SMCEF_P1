import numpy as np
import matplotlib.pyplot as plt

# Function to initialize spins randomly
def initialize_spins(N):
    spins = np.random.choice([-1, 1], size=(N, N, N))
    return spins

# Function to calculate energy of a spin configuration
def calculate_energy(spins, J, H):
    N = spins.shape[0]
    energy = 0
    for i in range(N):
        for j in range(N):
            for k in range(N):
                energy += -J * spins[i, j, k] * (spins[(i+1)%N, j, k] + spins[(i-1)%N, j, k] +
                                                  spins[i, (j+1)%N, k] + spins[i, (j-1)%N, k] +
                                                  spins[i, j, (k+1)%N] + spins[i, j, (k-1)%N]) - H * spins[i, j, k]
    return energy

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
    spins = initialize_spins(N)
    energies = np.zeros(steps)
    for step in range(steps):
        spins = monte_carlo_step(spins, J, H, T)
        energies[step] = calculate_energy(spins, J, H)
    return spins, energies

# Main function
def main():
    N = 10  # Size of the grid
    J = 1   # Interaction strength
    H = 0   # External magnetic field
    T = 1   # Temperature
    steps = 1000  # Number of Monte Carlo steps

    spins, energies = simulate(N, J, H, T, steps)

    # Plot energy vs. step
    plt.plot(energies)
    plt.xlabel('Step')
    plt.ylabel('Energy')
    plt.title('Energy vs. Step')
    plt.show()

if __name__ == "__main__":
    main()
