void time_sweep(int model, int n_part, int n_bond, Part *vec_part, Bond *vec_bond, double kn, double ks, double mun, double mus)
{
    double ampl[3], omega, t_tot, n_double, t_out, n_out_double, dt;
    int n, n_out, count=0;

    std::cout << "Insert oscillation amplitude vector applied to top layer of particles [m]:  ";
    std::cin >> ampl[0] >> ampl[1] >> ampl[2];

    std::cout << "Insert oscillation frequency applied to top layer of particles [rad/s]:  ";
    std::cin >> omega;

    std::cout <<"Insert total time of simulation [s]: ";
    std::cin >> t_tot;

    std::cout << "Insert simulation time-step [s]: ";
    std::cin >> dt;

    std::cout << "Insert output time interval [s]: ";
    std::cin >> t_out;

    n_double = t_tot/dt;
    n = (int) n_double;

    n_out_double = t_out/dt;
    n_out = (int) n_out_double;

    std::ofstream out_file;

    if (model==0)
        out_file.open("time_sweep_undamped.csv");

    else if (model==1)
        out_file.open("time_sweep_kelvin.csv");

    else
        out_file.open("time_sweep_maxwell.csv");

    out_file << "timestep,ID,x,y,z,u,v,w,fx,fy,fz,omx,omy,omz,tx,ty,tz\n";
    bool flag = false;

    double *xi, *xtop, *xbot;
    double d;
    int ind_top, ind_bot;

    xtop = vec_part[n_part-1].getPos();
    xbot = vec_part[0].getPos();
    d = vec_part[0].getDiameter();

    for (int i=0; i<n_part; i++)
    {
        xi = vec_part[i].getPos();

        if (fabs(xi[2] - xtop[2]) < 1.e-15)
        {
            ind_top = i;
            break;
        }
    }

    for (int i=0; i<n_part; i++)
    {
        xi = vec_part[i].getPos();

        if (fabs(xi[2] - xbot[2]) > d/2.)
        {
            ind_bot = i-1;
            break;
        }
    }

    std::cout << ind_top << "," << ind_bot << "\n";
    std::cout << n_part << "\n";

    for (int i=0; i<n; i++)
    {
        for (int j=ind_bot+1; j<n_part; j++)
        {
            if (j<ind_top)
            {
                Part *particle;
                particle = &vec_part[j];

                integrate(particle, dt, i);
            }
            else
            {
               double *xin, *x, *v;

               xin = vec_part[j].getPosIn();
               x = vec_part[j].getPos();
               v = vec_part[j].getVel();

               x[0] = xin[0] + ampl[0]*sin(omega*dt*i);
               x[1] = xin[1] + ampl[1]*sin(omega*dt*i);
               x[2] = xin[2] + ampl[2]*sin(omega*dt*i);

               v[0] = omega*ampl[0]*cos(omega*dt*i);
               v[1] = omega*ampl[1]*cos(omega*dt*i);
               v[2] = omega*ampl[2]*cos(omega*dt*i);
            }
        }

        for (int j=0; j<n_bond; j++)
        {
            int *IDS;
            Bond *bondino;

            IDS = vec_bond[j].getIds();
            bondino = &vec_bond[j];

            if (model==0)
                forces_bond_undamped(bondino, vec_part[IDS[0]], vec_part[IDS[1]], kn, ks, dt, i);

            else if (model==1)
                forces_bond_kelvin(bondino, vec_part[IDS[0]], vec_part[IDS[1]], kn, ks, mun, mus, dt, i);

            else
                forces_bond_maxwell(bondino, vec_part[IDS[0]], vec_part[IDS[1]], kn, ks, mun, mus, dt, i);
        }

        /*Forces on particles*/

        for (int j=0; j<n_part; j++)
        {
            double F_part[3]={}, T_part[3]={};

            /*Bonds contribution*/

            for (int k=0; k<n_bond; k++)
            {
                int *IDS;

                IDS = vec_bond[k].getIds();

                if (j == IDS[0])
                {
                    F_part[0] += -(vec_bond[k].getForce_n()[0] + vec_bond[k].getForce_s()[0]);
                    F_part[1] += -(vec_bond[k].getForce_n()[1] + vec_bond[k].getForce_s()[1]);
                    F_part[2] += -(vec_bond[k].getForce_n()[2] + vec_bond[k].getForce_s()[2]);

                    T_part[0] += -(vec_bond[k].getTorque_n()[0] + vec_bond[k].getTorque_s()[0]);
                    T_part[1] += -(vec_bond[k].getTorque_n()[1] + vec_bond[k].getTorque_s()[1]);
                    T_part[2] += -(vec_bond[k].getTorque_n()[2] + vec_bond[k].getTorque_s()[2]);
                }

                else if (j == IDS[1])
                {
                    F_part[0] += (vec_bond[k].getForce_n()[0] + vec_bond[k].getForce_s()[0]);
                    F_part[1] += (vec_bond[k].getForce_n()[1] + vec_bond[k].getForce_s()[1]);
                    F_part[2] += (vec_bond[k].getForce_n()[2] + vec_bond[k].getForce_s()[2]);

                    T_part[0] += (vec_bond[k].getTorque_n()[0] + vec_bond[k].getTorque_s()[0]);
                    T_part[1] += (vec_bond[k].getTorque_n()[1] + vec_bond[k].getTorque_s()[1]);
                    T_part[2] += (vec_bond[k].getTorque_n()[2] + vec_bond[k].getTorque_s()[2]);
                }
            }

            vec_part[j].setForce(F_part);
            vec_part[j].setTorque(T_part);
        }

        for(int j=0; j<n_part; j++)
        {
            if (i == (count*n_out))
            {
                out_file << i << "," << j << ","
                         << vec_part[j].getPos()[0] << "," << vec_part[j].getPos()[1] << ","
                         << vec_part[j].getPos()[2] << "," << vec_part[j].getVel()[0] << ","
                         << vec_part[j].getVel()[1] << "," << vec_part[j].getVel()[2] << ","
                         << vec_part[j].getForce()[0] << "," << vec_part[j].getForce()[1] << ","
                         << vec_part[j].getForce()[2] << "," << vec_part[j].getRot()[0] << ","
                         << vec_part[j].getRot()[1] << "," << vec_part[j].getRot()[2] << ","
                         << vec_part[j].getTorque()[0] << "," << vec_part[j].getTorque()[1] << ","
                         << vec_part[j].getTorque()[2] << "," << "\n";
                flag = true;
            }
        }

        if (flag)
        {
            count++;
            flag = false;
            std::cout << count << "\n";
        }
    }
}

