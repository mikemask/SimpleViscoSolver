void forces_bond_maxwell(Bond *bond, Part parti, Part partj, double kn, double ks, double mun, double mus, double dt, int i)
{
    double *omi, *omj;
    double *posi, *posj, *veli, *velj;
    double *posi_in, *posj_in;
    double radius, *inertia, area;
    double *r_old, *fn, *fs, *tn, *ts;

    /*Getting bond values*/

    double lambdan, lambdas;

    r_old = bond -> getRold();

    fn = bond -> getForce_n();
    fs = bond -> getForce_s();
    tn = bond -> getTorque_n();
    ts = bond -> getTorque_s();

    inertia = bond -> getInertia();

    radius = bond -> getRad();
    area = M_PI*radius*radius;

    lambdan = mun/(kn*area);
    lambdas = mus/(ks*area);

    /*Getting particles positions, velocities and torques*/

    posi = parti.getPos();
    posj = partj.getPos();

    posi_in = parti.getPosIn();
    posj_in = partj.getPosIn();

    veli = parti.getVel();
    velj = partj.getVel();

    omi = parti.getRot();
    omj = partj.getRot();

    /*Computing relative motion and contact plane versors*/

    double r_f[3], r_0[3], v_rel[3], v_c[3];

    r_f[0] = posi[0] - posj[0];
    r_f[1] = posi[1] - posj[1];
    r_f[2] = posi[2] - posj[2];

    r_0[0] = posi_in[0] - posj_in[0];
    r_0[1] = posi_in[1] - posj_in[1];
    r_0[2] = posi_in[2] - posj_in[2];

    v_rel[0] = veli[0] - velj[0];
    v_rel[1] = veli[1] - velj[1];
    v_rel[2] = veli[2] - velj[2];

    double x_cont[3];

    x_cont[0] = 0.5*(posi[0] + posj[0]);
    x_cont[1] = 0.5*(posi[1] + posj[1]);
    x_cont[2] = 0.5*(posi[2] + posj[2]);

    int e;

    for (int i=0; i<3; i++)
    {
        double perm = 0.;

        for (int j=0; j<3; j++)
        {
            for (int k=0; k<3; k++)
            {
                if (i==j || i==k || j==k || i==j==k)
                    e = 0;
                else if ((i==1 && j==2 && k==3) || (i==2 && j==3 && k==1) || (i==3 && j==1 && k==2))
                    e = 1;
                else
                    e = -1;

                perm += e*(omi[j]*(x_cont[k] - posi[k]) - omj[j]*(x_cont[k] - posj[k]));
            }
        }

        v_c[i] = v_rel[i] + perm;
    }

    double mod_rf, mod_r0, mod_rold, beta;
    double n[3], S[3], s[3];

    mod_rf = modulus(r_f);
    mod_r0 = modulus(r_0);
    mod_rold = modulus(r_old);

    n[0] = r_f[0]/mod_rf;
    n[1] = r_f[1]/mod_rf;
    n[2] = r_f[2]/mod_rf;

    S[0] = r_f[1]*(r_f[0]*r_0[1]-r_0[0]*r_f[1]) - r_f[2]*(r_0[0]*r_f[2]-r_f[0]*r_0[2]);
    S[1] = r_f[2]*(r_f[1]*r_0[2]-r_0[1]*r_f[2]) - r_f[0]*(r_0[1]*r_f[0]-r_f[1]*r_0[0]);
    S[2] = r_f[0]*(r_f[2]*r_0[0]-r_0[2]*r_f[0]) - r_f[1]*(r_0[2]*r_f[1]-r_f[2]*r_0[1]);

    double vc_s;
    double mod_S = modulus(S);

    if (mod_S < 1.e-30)
    {
        s[0] = 0.;
        s[1] = 0.;
        s[2] = 0.;
    }

    else
    {
        s[0] = S[0]/mod_S;
        s[1] = S[1]/mod_S;
        s[2] = S[2]/mod_S;
    }

    vc_s = v_c[0]*s[0] + v_c[1]*s[1] + v_c[2]*s[2];

    beta = acos((r_f[0]*r_0[0] + r_f[1]*r_0[1] + r_f[2]*r_0[2])/(mod_rf*mod_r0));

    if (std::isnan(beta))
        beta = 0.;

    double delta_theta[3];

    delta_theta[0] = (omi[0] - omj[0])*dt;
    delta_theta[1] = (omi[1] - omj[1])*dt;
    delta_theta[2] = (omi[2] - omj[2])*dt;

    double delta_theta_n, delta_theta_s;

    delta_theta_n = delta_theta[0]*n[0] + delta_theta[1]*n[1] + delta_theta[2]*n[2];
    delta_theta_s = delta_theta[0]*s[0] + delta_theta[1]*s[1] + delta_theta[2]*s[2];

    /*Computing forces and torques*/

    if (i==0)
    {
        double En, Es;

        En = 2.*kn*area/(M_PI*radius);
        Es = 2.*ks*area/(M_PI*radius);

        std::cout << "En: " << En << "\n";

        fn[0] = En*area*(mod_rf - mod_r0)/mod_r0*n[0];
        fn[1] = En*area*(mod_rf - mod_r0)/mod_r0*n[1];
        fn[2] = En*area*(mod_rf - mod_r0)/mod_r0*n[2];

        fs[0] = Es*area*beta*s[0];
        fs[1] = Es*area*beta*s[1];
        fs[2] = Es*area*beta*s[2];

    }

    else
    {
        double fn_old, fs_old;

        fn_old = fn[0]*n[0] + fn[1]*n[1] + fn[2]*n[2];
        fs_old = fs[0]*s[0] + fs[1]*s[1] + fs[2]*s[2];

        fn[0] = (fn_old/(1. + dt/lambdan) + (mod_rf - mod_rold)/(1./(kn*area) + dt/mun))*n[0];
        fn[1] = (fn_old/(1. + dt/lambdan) + (mod_rf - mod_rold)/(1./(kn*area) + dt/mun))*n[1];
        fn[2] = (fn_old/(1. + dt/lambdan) + (mod_rf - mod_rold)/(1./(kn*area) + dt/mun))*n[2];

        fs[0] = (fs_old/(1. + dt/lambdas) + vc_s*dt/(1./(ks*area) + dt/mus))*s[0];
        fs[1] = (fs_old/(1. + dt/lambdas) + vc_s*dt/(1./(ks*area) + dt/mus))*s[1];
        fs[2] = (fs_old/(1. + dt/lambdas) + vc_s*dt/(1./(ks*area) + dt/mus))*s[2];
    }

    tn[0] += ks*inertia[1]*delta_theta_n*n[0];
    tn[1] += ks*inertia[1]*delta_theta_n*n[1];
    tn[2] += ks*inertia[1]*delta_theta_n*n[2];

    ts[0] += kn*inertia[0]*delta_theta_s*s[0];
    ts[1] += kn*inertia[0]*delta_theta_s*s[1];
    ts[2] += kn*inertia[0]*delta_theta_s*s[2];

    bond -> setRold(r_f);

    if (std::isnan(fs[0]) || std::isnan(fs[1]) || std::isnan(fs[2]))
    {
        std::cout << "error at iteration: " << i << "\n";
        std::cout << "normal force: " << fn[0] << "," << fn[1] << "," << fn[2] << "\n";
        std::cout << "shear versor: " << s[0] << "," << s[1] << "," << s[2] << "\n";
        std::cout << "beta: " << beta << "\n";
        exit(0);
    }
}
