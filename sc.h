void std_cube(double d, double eps, int np_x, int np_y, int np_z, Part *vec_part, double rho)
{
    double pos[3] = {};
    int id = 0;

    double const m = rho*(4./3.)*pi*pow((d/2.),3.);

    for (int i=0; i<np_z; i++)
    {
        for (int j=0; j<np_y; j++)
        {
            for (int k=0; k<np_x; k++)
            {
                pos[0] = k*(d+eps);
                pos[1] = j*(d+eps);
                pos[2] = i*(d+eps);

                vec_part[id].setPos(pos);
                vec_part[id].setPosIn(pos);

                vec_part[id].setMass(m);
                vec_part[id].setDiameter(d);
                id++;
            }
        }
    }
}
