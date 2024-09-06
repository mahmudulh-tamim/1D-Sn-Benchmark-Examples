//
// Created by mhtam on 9/6/2024.
//
#include<iostream>
using namespace std;
#include<cmath>
#include<fstream>
#include<cassert>


int main() {

    /*data */
    double q_ext,sigma_s,q_s,sigma_t,X,del_x,psi_left,psi_right,a_number_for_which_I_have_to_calculate_e_base_log;

    q_ext=3;
    //therefore
    q_s=q_ext/2;
    sigma_t=1;
    a_number_for_which_I_have_to_calculate_e_base_log=2;
    sigma_s=log(a_number_for_which_I_have_to_calculate_e_base_log);
    //problem domain
    X=1;
    psi_left=0;
    psi_right=0;


    const int tot_cell=20, tot_quad_point=16;
    del_x=X/tot_cell;

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


    double psi_avg,psi_in,psi_out,scalar_flux[tot_cell],scalar_flux_temp[tot_cell],q[tot_cell];

    for (int i=0;i<tot_cell;i++)
    {
        scalar_flux[i]=0;
        scalar_flux_temp[i]=0;
    }

    int c_count, ang_count,itr_count=0;
    double tol=0.001, diff_2=1,square_diff_2;
    int test;
    while (diff_2>tol)
    {

        for (int i=0;i<tot_cell;i++)
        {
            q[i]=q_s+(sigma_s/2.0)*scalar_flux[i];
            scalar_flux_temp[i]=scalar_flux[i];
            scalar_flux[i]=0;
        }
        for (ang_count=tot_quad_point-1;ang_count>=0;ang_count--)
        {
            if (angl[ang_count][0]>0) {
                //forward sweep
                angl[ang_count][0]=fabs(angl[ang_count][0]);
                for (c_count=0;c_count<tot_cell;c_count++) {
                    if(c_count==0) {
                        psi_in=psi_left;

                    }
                    else {
                        psi_in=psi_out;
                    }

                    psi_out=(2*q[c_count]*del_x-psi_in*(sigma_t*del_x-2*angl[ang_count][0]))/(sigma_t*del_x+2*angl[ang_count][0]);
                    psi_avg=(psi_in+psi_out)/2;
                    scalar_flux[c_count]=scalar_flux[c_count]+psi_avg*angl[ang_count][1];
                }
            }
            else {
                //backward sweep
                angl[ang_count][0]=fabs(angl[ang_count][0]);
                for (c_count=tot_cell-1;c_count>=0;c_count--){
                    if(c_count==tot_cell-1){
                        psi_in=psi_right;
                    }
                    else{
                        psi_in=psi_out;
                    }


                    psi_out=(2*q[c_count]*del_x-psi_in*(sigma_t*del_x-2*angl[ang_count][0]))/(sigma_t*del_x+2*angl[ang_count][0]);
                    psi_avg=(psi_in+psi_out)/2;
                    scalar_flux[c_count]=scalar_flux[c_count]+psi_avg*angl[ang_count][1];
                }
            }



        }
        square_diff_2=0;
        for (int i=0;i<tot_cell;i++) {
            square_diff_2=square_diff_2+(scalar_flux[i]-scalar_flux_temp[i])*(scalar_flux[i]-scalar_flux_temp[i]);
        }
        diff_2=sqrt(square_diff_2);
        itr_count+=1;
    }
    cout<<"Total Number of Iterations:"<<itr_count<<endl;
    std::ofstream write_output_scalar_flux("Output_scalar_flux.txt");
    assert(write_output_scalar_flux.is_open());
    cout<<"Writing scalar fluxes to Output_scalar_flux.txt"<<endl;
    for (int j=0;j<tot_cell;j++)
        {
            write_output_scalar_flux<<scalar_flux[j]<<"\n";
            cout<<scalar_flux[j]<<endl;
        }

    write_output_scalar_flux.close();


    return 0;


}