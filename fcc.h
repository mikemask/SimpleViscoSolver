int fc_cube(double d, double eps, int n_cols, int n_rows, int n_lays, Part *vec_part, double rho)
{
    double pos[3]={};
    int id = 0;

    double const a = (d+eps)*sqrt(2.);
    double const m = rho*(4./3.)*pi*pow((d/2.),3.);

    double iterd = floor(n_cols/2.);
    int iter = int(iterd);

    for (int i=0; i<n_lays; i++)
    {
        for (int j=0; j<n_rows; j++)
        {
            if ((i==0 || i%2==0) && (j==0 || j%2==0))
            {
                for (int k=0; k<=iter; k++)
                {
                    pos[0] = k*a;
                    pos[1] = j*0.5*a;
                    pos[2] = i*0.5*a;

                    vec_part[id].setPos(pos);
                    vec_part[id].setPosIn(pos);
                    vec_part[id].setMass(m);
                    vec_part[id].setDiameter(d);
                    id++;
                }
            }

            else if ((i==0 || i%2==0) && (j%2!=0))
            {
                for (int k=1; k<=iter; k++)
                {
                    pos[0] = k*0.5*a;
                    pos[1] = j*0.5*a;
                    pos[2] = i*0.5*a;

                    vec_part[id].setPos(pos);
                    vec_part[id].setPosIn(pos);
                    vec_part[id].setMass(m);
                    vec_part[id].setDiameter(d);
                    id++;
                }
            }

            else if ((i%2!=0) && (j==0 || j%2==0))
            {
                for (int k=1; k<=iter; k++)
                {
                    pos[0] = k*0.5*a;
                    pos[1] = j*0.5*a;
                    pos[2] = i*0.5*a;

                    vec_part[id].setPos(pos);
                    vec_part[id].setPosIn(pos);
                    vec_part[id].setMass(m);
                    vec_part[id].setDiameter(d);
                    id++;
                }
            }

            else if ((i%2!=0) && (j%2!=0))
            {
                for (int k=0; k<=iter; k++)
                {
                    pos[0] = k*a;
                    pos[1] = j*0.5*a;
                    pos[2] = i*0.5*a;

                    vec_part[id].setPos(pos);
                    vec_part[id].setPosIn(pos);
                    vec_part[id].setMass(m);
                    vec_part[id].setDiameter(d);
                    id++;
                }
            }
        }
    }
}
