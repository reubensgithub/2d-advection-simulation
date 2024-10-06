/*******************************************************************************
2D advection example program which advects a Gaussian u(x,y) at a fixed velocity



Outputs: initial.dat - inital values of u(x,y)
         final.dat   - final values of u(x,y)

         The output files have three columns: x, y, u

         Compile with: gcc -o advection2D -std=c99 advection2D.c -lm

Notes: The time step is calculated using the CFL condition

********************************************************************************/

/*********************************************************************
                     Include header files
**********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

/*********************************************************************
                      Main function
**********************************************************************/

int main(){

  /* Grid properties */
  const int NX=1000;    // Number of x points
  const int NY=1000;    // Number of y points
  const float xmin=0.0; // Minimum x value
  const float xmax=30.0; // Maximum x value
  const float ymin=0.0; // Minimum y value
  const float ymax=30.0; // Maximum y value

  /* Parameters for the Gaussian initial conditions */
  const float x0=3.0;                    // Centre(x)
  const float y0=15.0;                    // Centre(y)
  const float sigmax=1.0;               // Width(x)
  const float sigmay=5.0;               // Width(y)
  const float sigmax2 = sigmax * sigmax; // Width(x) squared
  const float sigmay2 = sigmay * sigmay; // Width(y) squared
  
  /* Adding vertical shear (Question 2.3) */
  const float ustar = 0.2;
  const float z0 = 1.0;
  const float vk = 0.41;

  /* Boundary conditions */
  const float bval_left=0.0;    // Left boudnary value
  const float bval_right=0.0;   // Right boundary value
  const float bval_lower=0.0;   // Lower boundary
  const float bval_upper=0.0;   // Upper bounary

  /* Time stepping parameters */
  const float CFL=0.9;   // CFL number
  const int nsteps=800; // Number of time steps

  /* Velocity */
  float velx=1.0; // Velocity in x direction
  const float vely=0.0; // Velocity in y direction

  /* Arrays to store variables. These have NX+2 elements
     to allow boundary values to be stored at both ends */
  float x[NX+2];          // x-axis values
  float y[NX+2];          // y-axis values
  float u[NX+2][NY+2];    // Array of u values
  float dudt[NX+2][NY+2]; // Rate of change of u

  float x2;   // x squared (used to calculate iniital conditions)
  float y2;   // y squared (used to calculate iniital conditions)

  /* Calculate distance between points */
  float dx = (xmax-xmin) / ( (float) NX);
  float dy = (ymax-ymin) / ( (float) NY);

  /* Calculate time step using the CFL condition */
  /* The fabs function gives the absolute value in case the velocity is -ve */
  float dt = CFL / ( (fabs(velx) / dx) + (fabs(vely) / dy) );

  /*** Report information about the calculation ***/
  printf("Grid spacing dx     = %g\n", dx);
  printf("Grid spacing dy     = %g\n", dy);
  printf("CFL number          = %g\n", CFL);
  printf("Time step           = %g\n", dt);
  printf("No. of time steps   = %d\n", nsteps);
  printf("End time            = %g\n", dt*(float) nsteps);
  printf("Distance advected x = %g\n", velx*dt*(float) nsteps);
  printf("Distance advected y = %g\n", vely*dt*(float) nsteps);

  /*** Place x points in the middle of the cell ***/
  /* LOOP 1 */
  // This loop (loop 1) can be parallelised.
  #pragma omp parallel for shared(u,x,y)
  for (int i=0; i<NX+2; i++){
    x[i] = ( (float) i - 0.5) * dx;
  }

  /*** Place y points in the middle of the cell ***/
  /* LOOP 2 */
  // This loop (loop 2) can be parallelised.
  #pragma omp for
  for (int j=0; j<NY+2; j++){
    y[j] = ( (float) j - 0.5) * dy;
  }

  /*** Set up Gaussian initial conditions ***/
  /* LOOP 3 */
  // This loop (loop 3) can be parallelised.
  #pragma omp for collapse(2)
  for (int i=0; i<NX+2; i++){
    for (int j=0; j<NY+2; j++){
      x2      = (x[i]-x0) * (x[i]-x0);
      y2      = (y[j]-y0) * (y[j]-y0);
      u[i][j] = exp( -1.0 * ( (x2/(2.0*sigmax2)) + (y2/(2.0*sigmay2)) ) );
    }
  }

  /*** Write array of initial u values out to file ***/
  FILE *initialfile;
  initialfile = fopen("initial.dat", "w");
  /* LOOP 4 */
  // This loop cannot be parallelised. Due to the file writing to initial.dat
  // being sequential, parallelising this loop can cause race conditions and
  // the file may write in an incorrect order.
  for (int i=0; i<NX+2; i++){
    for (int j=0; j<NY+2; j++){
      fprintf(initialfile, "%g %g %g\n", x[i], y[j], u[i][j]);
    }
  }
  fclose(initialfile);

  /*** Update solution by looping over time steps ***/
  /* LOOP 5 */
  // This loop (loop 5) can be parallelised.
  #pragma omp parallel shared(u,x,y,dudt)
  for (int m=0; m<nsteps; m++){

    /*** Apply boundary conditions at u[0][:] and u[NX+1][:] ***/
    /* LOOP 6 */
    // This loop (loop 6) can be parallelised.
    #pragma omp single
    for (int j=0; j<NY+2; j++){
      u[0][j]    = bval_left;
      u[NX+1][j] = bval_right;
    }

    /*** Apply boundary conditions at u[:][0] and u[:][NY+1] ***/
    /* LOOP 7 */
    // This loop (loop 7) can be parallelised.
    #pragma omp single
    for (int i=0; i<NX+2; i++){
      u[i][0]    = bval_lower;
      u[i][NY+1] = bval_upper;
    }

    /*** Calculate rate of change of u using leftward difference ***/
    /* Loop over points in the domain but not boundary values */
    /* LOOP 8 */
    // This loop (loop 8) can be parallelised.
    #pragma omp for collapse(2) private(velx)
    for (int i=1; i<NX+1; i++){
      for (int j=1; j<NY+1; j++){
		if (y[j] > z0){
	        velx = (ustar/vk) * log(y[j]/z0);
	    } else {
            velx = 0.0;
	    }
        dudt[i][j] = -velx * (u[i][j] - u[i-1][j]) / dx
                    - vely * (u[i][j] - u[i][j-1]) / dy;
		dudt[i][j] *= 0.6;
      }
    }

    /*** Update u from t to t+dt ***/
    /* Loop over points in the domain but not boundary values */
    /* LOOP 9 */
    // This loop (loop 9) can be parallelised.
    #pragma omp for collapse(2)
    for (int i=1; i<NX+1; i++){
      for (int j=1; j<NY+1; j++){
		u[i][j] = u[i][j] + dudt[i][j] * dt;
      }
    }

  } // time loop

  /*** Write array of final u values out to file ***/
  FILE *finalfile;
  finalfile = fopen("final.dat", "w");
  /* LOOP 10 */
  // This loop cannot be parallelised. Due to the same reasons as loop 4.
  // As the file writing to final.dat is sequential, parallelising this process
  // may not be ideal due to race conditions and may produce the file in an
  // incorrect order.
  for (int i=0; i<NX+2; i++){
    for (int j=0; j<NY+2; j++){
      fprintf(finalfile, "%g %g %g\n", x[i], y[j], u[i][j]);
    }
  }
  fclose(finalfile);
  
  
  /*** 2.4 (Task 4) Calculate vertically averaged distribution ***/
  /* Calculate sum of all u(x,y) values at each y position at that x
     position then divide by the total number of y values considered.
     Do not include boundary values.
  */
  
  float* total_sum = (float*)malloc((NX+2) * sizeof(float));

  for (int i=1; i<NX+1; i++){
	  int num_rows = NY;
	  float sum = 0.0;
	  for (int j=1; j<NY+1; j++){
		  sum += u[i][j];
	  }
	  float average_u_at_current_x = sum / num_rows;
	  total_sum[i-1] = average_u_at_current_x;
  }
  
  /*** Write array of total_sum values out to file ***/
  FILE *vertfile;
  vertfile = fopen("vert.dat", "w");
  /* Loop over total_sum array */
  for (int i=0; i<=NX; i++) {
	  float x_val = i * dx;
      fprintf(vertfile, "%f %f\n", x_val, total_sum[i]);
  }
  fclose(vertfile);
  free(total_sum);

  return 0;
}

/* End of file ******************************************************/