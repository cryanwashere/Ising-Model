#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define L 20    // the size of the lattice
#define J 1     // the coupling constant
#define kB 1    // the Boltzmann constant
#define T 10     // the temperature
#define steps 10000 // the number of Monte Carlo steps



// calculate the energy of a given spin configuration
double energy(int spins[L][L]) 
{
    double energy = 0;
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            int spin = spins[i][j];
            int up_spin = spins[i][(j+1)%L];
            int down_spin = spins[i][(j-1+L)%L];
            int right_spin = spins[(i+1)%L][j];
            int left_spin = spins[(i-1+L)%L][j];
            energy -= J * spin * (up_spin + down_spin + right_spin + left_spin);
        }
    }
    return energy;
}

// calculate the magnetization of the given spin configuration
double magnetization(int spins[L][L]) 
{
    double magnetization = 0.0;
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            magnetization += spins[i][j];
        }
    }
    return magnetization;
}


// flip a random spin in the lattice
void flip_spin(int spins[L][L], double * energy, double * magnetization) 
{
    int i = rand() % L;
    int j = rand() % L;

    int spin = spins[i][j];

    double delta_E = 2.0 * J * spin * (spins[(i+1)%L][j] + spins[i][(j+1)%L] + spins[(i-1+L)%L][j] + spins[i][(j-1+L)%L]);
    if (delta_E <= 0.0 || exp(-delta_E / (kB * T)) > ((double) rand() / (double) RAND_MAX)) {
        spins[i][j] = -spin;
        * energy += delta_E;
        * magnetization += 2.0 * spin;
    }
}

// Function to print the current spin configuration
void print_grid(int spins[L][L],)
{
    // Move cursor to top-left corner of console
    printf("\033[1;1H");

    // Move cursor down one line from where the command was entered
    //printf("\033[1E");

   

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            if (spins[i][j] == 1) {
                printf("\033[31m+\033[0m ");
            } else {
                printf("\033[34m-\033[0m ");
            }
        }
        printf("\n");
    }
    
    fflush(stdout); // Flush output buffer to ensure characters are printed immediately
}

int main() 
{
    //initialize the random number generator
    srand(time(NULL));

    // initialize the spin grid
    int spins[L][L];

    // the energy of the spin grid
    double E = 0.0;

    // the magnetization of the spin grid
    double M = 0.0;

    // initialize the spin grid randomly and 
    // calculate the inital magnetization
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            spins[i][j] = ((rand() % 2) * 2) - 1;
            M += spins[i][j];
        }
    }
    E = energy(spins);


    // Clear the screen and move the cursor to the top-left corner of the console
    printf("\033[2J");

    // run the Monte Carlo simulation
    for (int step = 0; step < steps; step ++) {
        flip_spin(spins, &E, &M);
        print_grid(spins);
    }
    //printf("\n");

    
    return 0;
}