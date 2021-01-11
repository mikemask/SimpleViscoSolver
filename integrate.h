void integrate(Part *part, double dt, int i)
{
    double *x, *v, *f, *om, *t;
    double *v_old, *f_old, *t_old;
    double x_act[3], v_act[3];
    double m, d, inertia;

    x = part -> getPos();
    v = part -> getVel();
    f = part -> getForce();
    om = part -> getRot();
    t = part -> getTorque();

    v_old = part -> getVelOld();
    f_old = part -> getForceOld();
    t_old = part -> getTorqueOld();

    m = part -> getMass();
    d = part -> getDiameter();
    inertia = (2./5.)*m*pow((d/2.),2.);

    if(i==1)
    {
        for (int k=0; k<3; k++)
        {
            v_act[k] = v[k];

            v[k] += (dt/m)*f[k];
            x[k] += v[k]*dt;
            om[k] += (dt/inertia)*t[k];
        }
    }

    else
    {
        for (int k=0; k<3; k++)
        {
            v_act[k] = v[k];

            v[k] += 0.5*(dt/m)*(3.*f[k] - f_old[k]);
            x[k] += 0.5*dt*(3.*v[k] - v_old[k]);
            om[k] += 0.5*(dt/inertia)*(3.*t[k] - t_old[k]);
        }
    }

    part -> setVelOld(v_act);
    part -> setForceOld(f);
    part -> setTorqueOld(t);

}
