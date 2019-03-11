/*
Simulating equilibrium magnetization in a 2-dimensional lattice using the
Ising model and the Metropolis Monte Carlo algorithm.
B = 0, J = 1, k = 1.
*/

#include <math.h>  /* exp */
#include <stdio.h>  /* input/output */
#include <stdlib.h>  /* exit, random, srandom */
#include <string.h>  /* memset */

long random_at_most(long max) {
  /* Generates random integer number within [0, max]
     from almost uniform distribution.
     Taken from Ryan Reich's answer on
     https://stackoverflow.com/questions/2509679/how-to-generate-a-random-
     integer-number-from-within-a-range. */
  unsigned long num_bins = (unsigned long) max + 1;
  unsigned long num_rand = (unsigned long) RAND_MAX + 1;
  unsigned long bin_size = num_rand / num_bins;
  unsigned long defect = num_rand % num_bins;  /* (remainder) */

  /* Create random number x from [0, RAND_MAX], "normalize" by bin_size.
     Prevent overflow by checking that x is small enough. */
  long x;
  do {
    x = random();
  } while (num_rand - defect <= (unsigned long) x);

  return x / bin_size;
}

long sum_nn(int lattice, long L, long idx) {
  /* Performs sum over nearest-neighbor states in 2-dim lattice,
     applying helical boundary conditions. */

	long Lsquare = L * L;
  int snn = lattice[(idx + 1) % Lsquare] +
            lattice[(idx - 1 + Lsquare) % Lsquare] +
            lattice[(idx + L) % Lsquare] +
            lattice[(idx - L + Lsquare) % Lsquare];
  return snn;
}

int main() {
  /* set seed */
  unsigned int seed = 1;
  void srandom(unsigned int seed);

  /* prompt for parameters  */
  /* initial temperature? */
  double Ti;
  printf("Specify initial temperature: \n");
  scanf("Ti = %f", &Ti);
  printf("\nFirst simulation runs at temperature %f.\n", Ti);

  /* temperature step? */
  double dT;
  printf("\nSpecify temperature step: \n");
  scanf("dT = %f", &dT);

  /* final temperature? */
  double Tf;
  printf("\nSpecify final temperature: \n");
  scanf("Tf = %f", &Tf);

  /* error message if sign(dT) does not match sign(Tf - Ti) */
  if ((Tf < Ti && dT > 0) || (Tf > Ti && dT < 0)) {
    fprintf(stderr, "\nCheck that the sign of dT is consistent with Ti and Tf! Exiting.\n");
    exit(1);
  }

  /* determine temperature range */
  double len_float = (Tf - Ti) / dT;
  int len_int = ((int) len_float) + 1;
  double Trange[len_int] = {Ti};  /* first entry is Ti */

  if (Tf == Ti) {
    printf("\nSimulation runs for temperature %f.\n", Ti);
  } else {
    for (int i = 1; i < len_int; i++) {
      Trange[i] = Trange[i - 1] + dT;
    }
    printf("\nSimulations run for temperatures %f through %f at steps of %f.\n",
           Ti, Trange[len_int - 1], dT);
  }

  /* pre-calculate Boltzmann factors for later,
     only cases Edif = 4 and Edif = 8 for each T are interesting (see below) */
  double b[len_int][2];
  for (int i = 0; i < len_int; i++) {
    b[i][0] = exp(-4.0 / Trange[i]);  /* if Edif = E_new - E_old = 4 */
    b[i][1] = exp(-8.0 / Trange[i]);  /* if Edif = 8 */
  }

  /* number of iterations? */
  int N;
  printf("\nSpecify number of iterations for each temperature: \n");
  scanf("N = %d", &N);
  printf("\nEach simulation has %d iterations.\n", N);

  /* lattice size? */
  long L;
  printf("Specify width of the square lattice: \n");
  scanf("L = %d", &L);  /* read in value of L */
  printf("\nProceeding with lattices of size %d x %d.", L, L);
  long Lsquare = L * L;

  /* make two initial lattice states for two parallel simulations,
     starting at T = 0 */
  int s1[Lsquare];
  int s2[Lsquare];
  memset(s1, -1, Lsquare);  /* T = 0, first possibility */
  memset(s2, 1, Lsquare);  /* T = 0, second possibility */

  /* Metropolis algorithm */
  for (int i = 0; i < len_int; i++) {
    /* iterate through Trange */

		/* open file for writing */
    sprintf(path1, "psm1_L%d_N%d_Ti%f_Tf%f_dT%f_Tact%f.csv",
						L, N, Ti, Trange[len_int - 1], dT, Trange[i]);
    sprintf(path2, "psm2_L%d_N%d_Ti%f_Tf%f_dT%f_Tact%f.csv",
						L, N, Ti, Trange[len_int - 1], dT, Trange[i]);
    FILE *fp1;
		FILE *fp2;
    fp1 = fopen(path1, "w");
    fp2 = fopen(path2, "w");

    if (fp1 == NULL) {
      fprintf(stderr, "Cannot open output file for simulation 1! Exiting.\n");
      exit(1);
    }
    if (fp2 == NULL) {
      fprintf(stderr, "Cannot open output file for simulation 2! Exiting.\n");
      exit(1);
    }

    /* intialize per-spin magnetization */
    double psm1[N + 1];  /* should include the initial state */
    double psm2[N + 1];
    for (int k = 0; k < Lsquare; k++) {
      psm1[0] = psm1[0] + s1[k];
      psm2[0] = psm2[0] + s2[k];
    }

    /* write initial values to files */
    /* output format: iteration,psm */
    fprintf(fp1, "0,%f\n", psm1[0]);
    fprintf(fp2, "0,%f\n", psm2[0]);

    for (int j = 1; j < N; j++) {
      /* Metropolis iterations at fixed T */

      /* choose random lattice site */
      long rand_idx1 = random_at_most(Lsquare - 1);
      long rand_idx2 = random_at_most(Lsquare - 1);

      /* create new state by flipping spin at chosen site,
         calculate energy difference (new - old),
         can takes values in [-8, -4, 0, 4, 8] */
      int Edif1 = 2 * s1[rand_idx1] * sum_nn(s1, L, rand_idx1);
      int Edif2 = 2 * s2[rand_idx2] * sum_nn(s2, L, rand_idx2);

      /* decide whether to flip spin permanently */
      /* first simulation */
      if (Edif1 <= 0) {
        /* flip spin */
        s1[rand_idx1] = -s1[rand_idx1];
      } else {
        /* flip spin with probability exp(-beta*Edif);
           since Edif > 0, only consider Edif = 4 or Edif = 8 */
        long rand01 = (double) random() / (double) ((unsigned) RAND_MAX + 1);
        if (rand01 < b[i][Edif1 / 4 - 1]) {
          /* random number in [0, 1) is smaller than acceptance ratio */
          s1[rand_idx1] = -s1[rand_idx1];
        }
      }

      /* second simulation */
      if (Edif2 <= 0) {
        /* flip spin */
        s2[rand_idx2] = -s2[rand_idx2];
      } else {
        /* flip spin with probability exp(-beta*Edif) */
        long rand01 = (double) random() / (double) ((unsigned) RAND_MAX + 1);
        if (rand01 < b[i][Edif2 / 4 - 1]) {
          /* random number in [0, 1) is smaller than acceptance ratio -> flip */
          s2[rand_idx2] = -s2[rand_idx2];
        }
      }

      /* compute per-spin magnetization at iteration j at fixed T */
      psm1[j + 1] = psm1[j] + 2 * s1[rand_idx1];
      psm2[j + 1] = psm2[j] + 2 * s2[rand_idx2];

      /* write to file */
      fprintf(fp1, "%d,%f\n", j + 1, psm1[j + 1]);
      fprintf(fp2, "%d,%f\n", j + 1, psm2[j + 1]);
    }
		fclose(fp1);
		fclose(fp2);
  }
  return 0;
}
