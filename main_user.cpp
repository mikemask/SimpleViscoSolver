#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <list>
#include "modulus.h"
#include "part.h"
#include "distance.h"
#include "prop_bond.h"
#include "sc.h"
#include "fcc.h"
#include "forces_bond_undamped.h"
#include "forces_bond_kelvin.h"
#include "forces_bond_maxwell.h"
#include "integrate.h"
#include "stress_relax.h"
#include "creep.h"
#include "time_sweep.h"


int main()
{
    double rho, d, E, nu, dt, lambda, kn, ks, mun, mus, eps;
    int type, model, n_part;

    std::cout << "Insert particles density [kg/m3]: ";
    std::cin >> rho;

    std::cout << "Insert particles diameter [m]: ";
    std::cin >> d;

    std::cout << "Insert particles Young's modulus [Pa]: ";
    std::cin >> E;

    std::cout << "Insert particles Poisson ratio: ";
    std::cin >> nu;

    std::cout << "Insert lattice type, sc=1 fcc=2: ";
    std::cin >> type;

    int np_x, np_y, np_z;
    int n_col, n_row, n_lay;

    if (type==1)
    {

        std::cout << "Insert number of particles in x-direction: ";
        std::cin >> np_x;

        std::cout << "Insert number of particles in y-direction: ";
        std::cin >> np_y;

        std::cout << "Insert number of particles in z-direction: ";
        std::cin >> np_z;

        std::cout << "Insert minimum distance between particles [m]: ";
        std::cin >> eps;

        n_part = np_x*np_y*np_z;
    }

    else
    {
        int nc_ext, nr_ext, nl_ext;

        std::cout << "Insert number of columns in x-direction: ";
        std::cin >> n_col;

        std::cout << "Insert number of rows in y-direction: ";
        std::cin >> n_row;

        std::cout << "Insert number of layers in z-direction: ";
        std::cin >> n_lay;

        std::cout << "Insert minimum distance between particles [m]: ";
        std::cin >> eps;

        nc_ext = ceil(n_col/2.);
        nr_ext = ceil(n_row/2.);
        nl_ext = ceil(n_lay/2.);

        n_part = nc_ext*nr_ext*nl_ext + (n_lay-nl_ext)*nc_ext + (n_col-nc_ext)*nl_ext + (n_lay-nl_ext)*(n_col-nc_ext)*nr_ext;
    }

    Part vec_part[n_part];

    switch (type)
    {
        case 1:

            std_cube(d, eps, np_x, np_y, np_z, vec_part, rho);

            break;

        case 2:

            fc_cube(d, eps, n_col, n_row, n_lay, vec_part, rho);

            break;
    }


    std::cout << "Insert visco-elastic model, undamped=0 Kelvin=1 Maxwell=2: ";
    std::cin >> model;

    std::cout << "Insert spring constant [Pa/m]: ";
    std::cin >> kn;

    ks = kn;

    if (model != 0)
    {
        std::cout << "Insert damper constant [Ns/m]: ";
        std::cin >> mun;

        mus = mun;
    }

    int n_bond = 0;
    int exp = 0;

    /*Bond initialization*/

    std::cout << "Insert bond radius multiplier: ";
    std::cin >> lambda;

    int pp = 1;
    int *vec_ids_1 = new int[n_part*10];
    int *vec_ids_2 = new int[n_part*10];
    double *vec_rad = new double[n_part*10];

    for (int i=0; i<n_part; i++)
    {
        for (int j=pp; j<n_part; j++)
        {
            if (i != j)
            {
                double *posi, *posj, radi, radj, radi_cr, radj_cr, dist;

                posi = vec_part[i].getPos();
                radi = vec_part[i].getDiameter()/2.;
                radi_cr = radi*lambda;

                posj = vec_part[j].getPos();
                radj = vec_part[j].getDiameter()/2.;
                radj_cr = radj*lambda;

                dist = distance(posi, posj);

                if (dist < (radi_cr+radj_cr))
                {
                    vec_ids_1[n_bond] = i;
                    vec_ids_2[n_bond] = j;
                    vec_rad[n_bond] = lambda*std::min(radi,radj);
                    n_bond++;
                }

            }
        }

        pp++;
    }

    Bond vec_bond[n_bond];

    for (int i=0; i<n_bond; i++)
    {
        double radius, inertia[2];
        int IDS[2];

        radius = vec_rad[i];

        IDS[0] = vec_ids_1[i];
        IDS[1] = vec_ids_2[i];

        inertia[0] = 0.25*3.14159265*pow(radius,4.);
        inertia[1] = 2.*inertia[0];

        vec_bond[i].setRad(radius);
        vec_bond[i].setIds(IDS);
        vec_bond[i].setInertia(inertia);

        std::cout << "bond: " << i << "\n"
                  << "ids: " << IDS[0] << "," << IDS[1] << "\n";
    }

    delete [] vec_ids_1;
    delete [] vec_ids_2;
    delete [] vec_rad;

    std::cout << "Which experiment do you want to perform? " << "\n"
              << "1 --> Stress relaxation" << "\n"
              << "2 --> Creep test" << "\n"
              << "3 --> Time sweep" << "\n"
              << "4 --> Frequency sweep" << "\n";

    std::cin >> exp;

    switch (exp)
    {
        case 1:
            stress_relax(model, n_part, n_bond, vec_part, vec_bond, kn, ks, mun, mus);
            break;

        case 2:
            creep(model, n_part, n_bond, vec_part, vec_bond, kn, ks, mun, mus);
            break;

        case 3:
            time_sweep(model, n_part, n_bond, vec_part, vec_bond, kn, ks, mun, mus);
            break;

        case 4:
            freq_sweep(model, n_part, n_bond, vec_part, vec_bond, kn, ks, mun, mus);
            break;
    }


    return 0;
}
