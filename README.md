# 2D Advection Simulation

This project simulates the advection of material from a chimney, carried by wind in the planetary boundary layer, using a Gaussian distribution. The simulation is implemented in C and optimised for high-performance computing using OpenMP for parallelisation.

### Specification

The program models a 2D advection scenario in which material is moved by wind across a grid. It uses a Gaussian function to initialise the material concentration and updates the state over time by solving advection equations. The time step is calculated using the Courant–Friedrichs–Lewy (CFL) condition, and the advection occurs with vertical shear effects factored in. The grid is evolved for a fixed number of time steps, and the output consists of initial and final material distributions.

### Prerequisites

To compile and run the simulation, the following tools are required:

- **Linux** or another POSIX-compliant operating system
- **GCC** with OpenMP support
- **gnuplot** (for visualising the results)

### Compilation Instructions

To compile the serial version of the program, use the following command:

```bash
gcc -o advection2D -std=c99 advection2D.c -lm
```
For the parallelised version (requires OpenMP support), use:

```bash
gcc -fopenmp -o advection2D -std=c99 advection2D.c -lm
```

### Running the Program

After compiling the program, you can run it with:

```bash
./advection2D
```

This will generate three output data files:
- **initial.dat** - Initial values of the concentration **u(x,y)**
- **final.dat** - Final values after the simulation
- **vert.dat** - Vertically averaged distribution at each x-position

### Visualisation
To visualise the results, you can use the provided gnuplot scripts. Make sure you have already generated the .dat files by running the program. Then, execute the following commands to create the plots:

```bash
gnuplot plot_initial
gnuplot plot_final
gnuplot plot_vert
```

Alternatively, you can use:

```bash
gnuplot plot_*
```
These commands will generate PNG images of the initial, final, and vertically averaged concentration distributions.

### Output Files
- **initial.png**: Shows the initial concentration distribution.
- **final.png**: Displays the final concentration distribution after advection.
If you encounter any issues generating the images, pre-generated results are available in the current directory.

### Program Details
- The grid has a size of 1000x1000 points, and boundary conditions are applied to all edges of the grid.
- The initial material concentration is modelled as a Gaussian function centred at a specific point in the grid.
- The simulation runs for 800 time steps, with the time step size calculated based on the CFL condition.
- Vertical shear is introduced to simulate realistic wind conditions.
### Notes
- Only certain loops have been parallelised to avoid race conditions, especially during file writing.
- The OpenMP directives are used to improve performance on multi-core machines by parallelising computations.
