#include<iostream>
#include<fstream>
#include<cmath>
#include<cassert>



int main()
{
    /*data */
    double q_ext,sigma_t,X,q,del_x;

    q_ext=3;
    //therefore
    q=3/2;
    sigma_t=1;
    //problem domain
    X=1;


    const int tot_cell=20, tot_quad_point=16;
    del_x=X/20;

    double angl[tot_quad_point][2]={
        {-0.9894009349916499, 0.0271524594117541},
        {-0.9445750230732326, 0.0622535239386479},
        {-0.8656312023878318, 0.0951585116824928},
        {-0.7554044083550030, 0.1246289712555339},
        {-0.6178762444026438, 0.1495959888165767},
        {-0.4580167776572274, 0.1691565193950025},
        {-0.2816035507792589, 0.1826034150449236},
        {-0.0950125098376374, 0.1894506104550685},
        { 0.0950125098376374, 0.1894506104550685},
        { 0.2816035507792589, 0.1826034150449236},
        { 0.4580167776572274, 0.1691565193950025},
        { 0.6178762444026438, 0.1495959888165767},
        { 0.7554044083550030, 0.1246289712555339},
        { 0.8656312023878318, 0.0951585116824928},
        { 0.9445750230732326, 0.0622535239386479},
        { 0.9894009349916499, 0.0271524594117541}
    };


    double psi_avg[tot_cell][tot_quad_point],psi_in[tot_cell][tot_quad_point],psi_out[tot_cell][tot_quad_point],scalar_flux[tot_cell];

    for (int i=0;i<tot_cell;i++)
    {
        scalar_flux[i]=0;
    }

    int c_count, ang_count;

    //forward sweep
    for (ang_count=tot_quad_point-1;ang_count>=(tot_quad_point/2);ang_count--)
    {
        for (c_count=0;c_count<tot_cell;c_count++)
            {
                if(c_count==0)
                {
                    psi_in[0][ang_count]=0;
                    psi_out[0][ang_count]=(2*q*del_x-psi_in[0][ang_count]*(sigma_t*del_x-2*angl[ang_count][0]))/(sigma_t*del_x+2*angl[ang_count][0]);
                    psi_avg[c_count][ang_count]=(psi_in[c_count][ang_count]+psi_out[c_count][ang_count])/2;
                    std::cout<<"Cell:"<<c_count<<"\n"<<"Angle:"<<ang_count<<"\n"<<psi_in[c_count][ang_count]<<"\n"<<psi_in[c_count][ang_count]<<"\n"<<psi_avg[c_count][ang_count]<<"\n";
                }
                else
                {
                    psi_in[c_count][ang_count]=psi_out[c_count-1][ang_count];
                    psi_out[c_count][ang_count]=(2*q*del_x-psi_in[c_count][ang_count]*(sigma_t*del_x-2*angl[ang_count][0]))/(sigma_t*del_x+2*angl[ang_count][0]);
                    psi_avg[c_count][ang_count]=(psi_in[c_count][ang_count]+psi_out[c_count][ang_count])/2;
                    std::cout<<"Cell:"<<c_count<<"\n"<<"Angle:"<<ang_count<<"\n"<<psi_in[c_count][ang_count]<<"\n"<<psi_out[c_count][ang_count]<<"\n"<<psi_avg[c_count][ang_count]<<"\n";
                }
                scalar_flux[c_count]=scalar_flux[c_count]+psi_avg[c_count][ang_count]*angl[ang_count][1];
            }

    }

    //backward sweep
    for (ang_count=0;ang_count<(tot_quad_point/2);ang_count++)
    {
        angl[ang_count][0]=fabs(angl[ang_count][0]);
        for (c_count=tot_cell-1;c_count>=0;c_count--)
            {
                if(c_count==tot_cell-1)
                {
                    psi_in[c_count][ang_count]=0;
                    psi_out[c_count][ang_count]=(2*q*del_x-psi_in[c_count][ang_count]*(sigma_t*del_x-2*angl[ang_count][0]))/(sigma_t*del_x+2*angl[ang_count][0]);
                    psi_avg[c_count][ang_count]=(psi_in[c_count][ang_count]+psi_out[c_count][ang_count])/2;
                    std::cout<<"Cell:"<<c_count<<"\n"<<"Angle:"<<ang_count<<"\n"<<psi_in[c_count][ang_count]<<"\n"<<psi_out[c_count][ang_count]<<"\n"<<psi_avg[c_count][ang_count]<<"\n";
                }
                else
                {
                    psi_in[c_count][ang_count]=psi_out[c_count+1][ang_count];
                    psi_out[c_count][ang_count]=(2*q*del_x-psi_in[c_count][ang_count]*(sigma_t*del_x-2*angl[ang_count][0]))/(sigma_t*del_x+2*angl[ang_count][0]);
                    psi_avg[c_count][ang_count]=(psi_in[c_count][ang_count]+psi_out[c_count][ang_count])/2;
                    std::cout<<"Cell:"<<c_count<<"\n"<<"Angle:"<<ang_count<<"\n"<<psi_in[c_count][ang_count]<<"\n"<<psi_out[c_count][ang_count]<<"\n"<<psi_avg[c_count][ang_count]<<"\n";
                }
                scalar_flux[c_count]=scalar_flux[c_count]+psi_avg[c_count][ang_count]*angl[ang_count][1];
            }


    }
    std::ofstream write_output_scalar_flux("Output_scalar_flux.txt");
    assert(write_output_scalar_flux.is_open());
    for (int j=0;j<tot_cell;j++)
        {
            write_output_scalar_flux<<scalar_flux[j]<<"\n";
        }

    write_output_scalar_flux.close();


    return 0;
}
