// Integration of the Lorenz system
// using Runge-Kutta (4) integration
//
// Lorenz system:
// dot{x} = 10(y - x)
// dot{y} = 28x - y - xz
// dot{z} = -(8/3)z + xy

// Problem with the xdiff calculation (difference between undisturbed & disturbed trajetory):
// choose initial value for xprime that is not a disturbance of x0 but of a value inside the attractor!

#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

// define functions dot{x}, dot{y} and dot{z} according to the 3-dim dynamical system above

double fx(double x, double y) {
  // this is dot{x}(x, y)
	return 10.0*(y - x);
}

double fy(double x, double y, double z) {
	// this is dot{y}(x, y, z)
	return 28.0*x - y - x*z;
}

double fz(double x, double y, double z) {
	// this is dot{z}(x, y, z)
	return -8.0*z/3.0 + x*y;
}

// we need a function for integrating
// three functions (dot{x}, dot{y}) at the same time

void integrate_xyz(double *px, double *py, double *pz, double dt) {
  // Runge-Kutta (4) coefficients
	double kx1 = fx(*px, *py);
	double ky1 = fy(*px, *py, *pz);
	double kz1 = fz(*px, *py, *pz);
	double kx2 = fx(*px + kx1*dt/2, *py + ky1*dt/2);  // for VDP osci, we used e.g. *py + kx1*dt/2, where comes the difference from?
	double ky2 = fy(*px + kx1*dt/2, *py + ky1*dt/2, *pz + kz1*dt/2);
	double kz2 = fz(*px + kx1*dt/2, *py + ky1*dt/2, *pz + kz1*dt/2);
	double kx3 = fx(*px + kx2*dt/2, *py + ky2*dt/2);
  double ky3 = fy(*px + kx2*dt/2, *py + ky2*dt/2, *pz + kz2*dt/2);
	double kz3 = fz(*px + kx2*dt/2, *py + ky2*dt/2, *pz + kz2*dt/2);
	double kx4 = fx(*px + kx3*dt, *py + ky3*dt);
	double ky4 = fy(*px + kx3*dt, *py + ky3*dt, *pz + kz3*dt);
	double kz4 = fz(*px + kx3*dt, *py + ky3*dt, *pz + kz3*dt);

	// integrate for one time step
	*px = *px + (kx1 + 2*kx2 + 2*kx3 + kx4)/6.0*dt;
	*py = *py + (ky1 + 2*ky2 + 2*ky3 + ky4)/6.0*dt;
	*pz = *pz + (kz1 + 2*kz2 + 2*kz3 + kz4)/6.0*dt;
}

int main() {
  // parameters
	double x0 = 1.0;
	double xprime0 = 1.01;
	double y0 = 1.0;
	double yprime0 = 1.0;
	double z0 = 1.0;
	double zprime0 = 1.0;
	double t0 = 0.0;
	double tmax = 100.0;
  double dt = 0.01;
  
  // output files
	ofstream fout_x("lorenz_x.dat");
  ofstream fout_xxprime("lorenz_xxprime.dat");
	ofstream fout_xdiff("lorenz_xdiff.dat");
	ofstream fout_y("lorenz_y.dat");
	ofstream fout_z("lorenz_z.dat");
	ofstream fout_xy("lorenz_xy.dat");
	ofstream fout_xyz("lorenz_xyz.dat");
  
  // initialize
	double x = x0;
	double xprime = xprime0;
	double y = y0;
	double yprime = yprime0;
	double z = z0;
	double zprime = zprime0;
  
  // integration over [t0, tmax]
	for (double t = t0; t < tmax; t = t + dt) {
    integrate_xyz(&x, &y, &z, dt);
		integrate_xyz(&xprime, &yprime, &zprime, dt);
		fout_x << t << ", " << x << endl;
		fout_xxprime << t << ", " << x << ", " << xprime << endl;
		fout_xdiff << t << ", " << log(x - xprime) << endl;
		fout_y << t << ", " << y << endl;
		fout_z << t << ", " << z << endl;
		fout_xy << x << ", " << y << endl;
		fout_xyz << x << ", " << y << ", " << z << endl;
	}
	
  return 0;
}
