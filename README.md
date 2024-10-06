# 2D Advection Simulation

This project simulates the advection of material from a chimney, carried by wind in the planetary boundary layer, using a Gaussian distribution. The simulation is implemented in C and optimised for high-performance computing using OpenMP for parallelisation.

### Specification

The program models a 2D advection scenario in which material is moved by wind across a grid. It uses a Gaussian function to initialise the material concentration and updates the state over time by solving advection equations. The time step is calculated using the Courant–Friedrichs–Lewy (CFL) condition, and the advection occurs with vertical shear effects factored in. The grid is evolved for a fixed number of time steps, and the output consists of initial and final material distributions.

### Prerequisites

To compile and run the simulation, the following tools are required:

- **Linux** or another POSIX-compliant operating system
- **GCC** with OpenMP support
- **Gnuplot** (for visualising the results)

### Compilation Instructions

To compile the serial version of the program, use the following command:

```bash
gcc -o advection2D -std=c99 advection2D.c -lm
