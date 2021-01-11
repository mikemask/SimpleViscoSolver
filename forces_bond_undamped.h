void forces_bond_undamped(Bond *bond, Part parti, Part partj, double kn, double ks, double dt, int i)
{
    double *omi, *omj, *veli, *velj;
    double *posi, *posj;
    double *posi_in, *posj_in;
    double radius, *inertia, area;
    double *fn, *fs, *tn, *ts;

    /*Getting bond values*/

    fn = bond -> getForce_n();
    fs = bond -> getForce_s();
    tn = bond -> getTorque_n();
    ts = bond -> getTorque_s();

    inertia = bond -> getInertia();

    radius = bond -> getRad();
    area = pi*radius*radius;

    /*Getting particles positions, velocities and torques*/

    posi = parti.getPos();
    posj = partj.getPos();

    posi_in = parti.getPosIn();
    posj_in = partj.getPosIn();

    omi = parti.getRot();
    omj = partj.getRot();

    veli = parti.getVel();
    velj = partj.getVel();

    /*Computing relative motion and contact plane versors*/

    double r_f[3], r_0[3], r_c[3], v_rel[3], v_c[3];

    r_f[0] = posi[0] - posj[0];
    r_f[1] = posi[1] - posj[1];
    r_f[2] = posi[2] - posj[2];

    r_0[0] = posi_in[0] - posj_in[0];
    r_0[1] = posi_in[1] - posj_in[1];
    r_0[2] = posi_in[2] - posj_in[2];

    v_rel[0] = veli[0] - velj[0];
    v_rel[1] = veli[1] - velj[1];
    v_rel[2] = veli[2] - velj[2];

    double mod_r0, mod_rf, beta;
    double n[3], S[3], s[3];

    mod_r0 = modulus(r_0);
    mod_rf = modulus(r_f);

    n[0] = r_f[0]/mod_rf;
    n[1] = r_f[1]/mod_rf;
    n[2] = r_f[2]/mod_rf;

    beta = acos((r_0[0]*r_f[0]+r_0[1]*r_f[1]+r_0[2]*r_f[2])/(mod_rf*mod_r0));

    if (std::isnan(beta))
        beta = 0.;

    S[0] = r_f[1]*(r_f[0]*r_0[1]-r_0[0]*r_f[1]) - r_f[2]*(r_0[0]*r_f[2]-r_f[0]*r_0[2]);
    S[1] = r_f[2]*(r_f[1]*r_0[2]-r_0[1]*r_f[2]) - r_f[0]*(r_0[1]*r_f[0]-r_f[1]*r_0[0]);
    S[2] = r_f[0]*(r_f[2]*r_0[0]-r_0[2]*r_f[0]) - r_f[1]*(r_0[2]*r_f[1]-r_f[2]*r_0[1]);

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

    double delta_theta[3];

    delta_theta[0] = (omi[0] - omj[0])*dt;
    delta_theta[1] = (omi[1] - omj[1])*dt;
    delta_theta[2] = (omi[2] - omj[2])*dt;

    double delta_theta_n, delta_theta_s;

    delta_theta_n = delta_theta[0]*n[0] + delta_theta[1]*n[1] + delta_theta[2]*n[2];
    delta_theta_s = delta_theta[0]*s[0] + delta_theta[1]*s[1] + delta_theta[2]*s[2];

    /*Computing forces and torques*/

    fn[0] = kn*area*(mod_rf - mod_r0)*n[0];
    fn[1] = kn*area*(mod_rf - mod_r0)*n[1];
    fn[2] = kn*area*(mod_rf - mod_r0)*n[2];

    fs[0] = ks*area*mod_r0*sin(beta)*s[0];
    fs[1] = ks*area*mod_r0*sin(beta)*s[1];
    fs[2] = ks*area*mod_r0*sin(beta)*s[2];

    tn[0] += ks*inertia[1]*delta_theta_n*n[0];
    tn[1] += ks*inertia[1]*delta_theta_n*n[1];
    tn[2] += ks*inertia[1]*delta_theta_n*n[2];

    ts[0] += kn*inertia[0]*delta_theta_s*s[0];
    ts[1] += kn*inertia[0]*delta_theta_s*s[1];
    ts[2] += kn*inertia[0]*delta_theta_s*s[2];

    if (std::isnan(fs[0]) || std::isnan(fs[1]) || std::isnan(fs[2]))
    {
        std::cout << "error at iteration: " << i << "\n";
        std::cout << "normal force: " << fn[0] << "," << fn[1] << "," << fn[2] << "\n";
        std::cout << "shear versor: " << s[0] << "," << s[1] << "," << s[2] << "\n";
        std::cout << "beta: " << beta << "\n";
        exit(0);
    }
}
