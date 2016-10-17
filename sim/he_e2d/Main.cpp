#include <string>
#include <vector>
#include "constant.h"
#include "Cartesian4D.h"
#include "ParameterMap.h"
#include "InputStuff.h" 
#include <fstream>
#include "Laser.h"
#include "mpi.h"
#include <time.h>
using namespace std;
using namespace Cartesian_4D;


#define DeBugOut(var) do{cout<< var <<endl;}while(0);                
/*****************************************************************************/
// grid parameters

Parameter<vector<int> > n_points( "n_points", vector<int>(4,50), "number of points for (X1,X2,X3,X4)" );
Parameter<vector<int> > n_points_read( "n_points_read", vector<int>(4,50), "number of points READED for (X1,X2,X3,X4)" );
Parameter<vector<double> > spatial_step( "spatial_step", vector<double>(4,0.3), "spatial_steps for (X1,X2,X3,X4)" );


// Gaussian Initial guess
Parameter<vector<double> > sigma("sigma", vector<double>(4,1.), "sigma");
Parameter<vector<double> > gauss_center("gauss_center", vector<double>(4,0.), "Center of Gaussian");
Parameter<vector<double> > initial_momentum("initial_momentum", vector<double>(4,0.), "Initial momentum of wavepacket");
Parameter<double> GroundStateEnergyThreshold("GroundStateEnergyThreshold",1.,"the value for telling imaginary time.");

// The Name of the job and wavefunction
Parameter<string> jobname("jobname", "Job", "The Name of the job");
Parameter<string> outputFolder("outputFolder", "out", "The basename of the output directory");  // just a default, is set in Input.cpp from jobname()
Parameter<string> wavefunctionFolder("wavefunctionFolder", "/data/npsf/lib/states/tmp", "The wavefunction archive");  // just a default, is set in Input.cpp from jobname()
Parameter<string> GroundStateName("GroundStateName","ground", "label of the waveFunction to load"); 
Parameter<string> NewStateName("NewStateName","ground", "label of the waveFunction to load"); 

// Laser parameters

Parameter<double> wavelength_nm_atto("wavelength_nm_atto",1., "wavelenth of the Laser(nm) ");
Parameter<double> intensity_wcm2_atto("intensity_wcm2_atto",3.5e16, "Intensity of the Laser (W/cm2)");
Parameter<double> fwhm_cycles_atto("fwhm_cycles_atto",2, "Number of cycles in the laser");
Parameter<double> N_cycles_atto("N_cycles_atto",2, "Number of cycles in the laser");
Parameter<double> cep_atto("cep_atto",0., "fraction of Pi for the initial phase");
Parameter<double> epsilon_atto("epsilon_atto",0., "the parameter for the ellipticity");
Parameter<double> pol_angle_atto("pol_angle_atto",0., "the parameter for the polarization.");

Parameter<double> wavelength_nm("wavelength_nm",1., "wavelenth of the Laser(nm) ");
Parameter<double> intensity_wcm2("intensity_wcm2",3.5e16, "Intensity of the Laser (W/cm2)");
Parameter<double> fwhm_cycles("fwhm_cycles",2, "Number of cycles in the laser");
Parameter<double> N_cycles("N_cycles",2, "Number of cycles in the laser");
Parameter<double> N_rise("N_rise",2, "Number of cycles in the laser");
Parameter<double> N_fall("N_fall",2, "Number of cycles in the laser");
Parameter<double> cep("cep",0., "fraction of Pi for the initial phase");
Parameter<double> epsilon("epsilon",0., "the parameter for the ellipticity");
Parameter<double> pol_angle("pol_angle",0., "the parameter for the polarization.");

Parameter<double> dt("dt",0.1, "Time step");
Parameter<double> delay("delay",0., "time delay between two laser pulses");

// target parameters

Parameter<double> soft_en("soft_en",1.,"for H2");
Parameter<double> soft_ee("soft_ee",1.,"for H2");
Parameter<double> soft_nn("soft_nn",1.,"for H2");
Parameter<double> r0("r0",1.4,"for H2");
Parameter<double> R_fix("R_fix",7,"for H2");
Parameter<double> R_fix1("R_fix1",7,"for H2");
Parameter<double> R_fix2("R_fix2",7,"for H2");
Parameter<double> R_fix3("R_fix3",7,"for H2");

Parameter<double> mass("mass",918.,"for h2+");
Parameter<double> mass_2("mass_2",1836.,"for h2+");
Parameter<double> mass_3("mass_3",1836.,"for h2+");
Parameter<double> charge_1("charge_1",1,"for h2+");
Parameter<double> charge_2("charge_2",1,"for h2+");
Parameter<double> charge_3("charge_3",1,"for h2+");
Parameter<double> femto("femto",1,"starting");

// save parameters

Parameter<int> store_sampler("store_sampler",1, "Store the wavefunction n times per cycle");
Parameter<double> pro_again("pro_again", 0., "for calculation");
Parameter<int> nthreads("nthreads",16, "the number of threads in parallel compution");

/*****************************************************************************/


int main( int argc , char *argv[] )
{
  
  parseInput(argc, argv); 
  int nlocal0, nlocal;
  int myid,numprocs,left,right;    
  double elapsed_time;
  MPI_Status status;
 

  MPI_Init (&argc, &argv);
  parseFile(argc, argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &myid);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  

  if ((n_points[3] % numprocs!=0)||(n_points_read[3]%numprocs !=0))
    {
      cerr<<"Data decomposition failed ! Please reset nodes ! " <<endl;     
      exit(1);
    }

  cout << "myid is : "  << myid << endl;
  cout << "numprocs is : "  << numprocs << endl;

  nlocal0=n_points_read[3]/numprocs;
  nlocal=n_points[3]/numprocs;

  // ////********************************************* start to comment
  
  wavefunction phi0;
  wavefunction phi;
  wavefunction phi_mom;


  Initialize( phi0 , n_points_read() , spatial_step(), 0, 0, 0,0);
  Initialize( phi , n_points() , spatial_step(), 0, 0, 0,0);
  Guess_Function_Gauss( phi0, sigma() , gauss_center() , initial_momentum() );
  

  //Potential_H_like_e3minus2D pot0(phi0);
  //Potential_H_like_e3minus2D pot(phi);

  Potential_He_4D pot0(phi0);
  Potential_He_4D pot(phi);

  // // Potential_Helium_plus_e3minus2D pot0(phi0);
  // // Potential_Helium_plus_e3minus2D pot(phi);

  Observable<double,wavefunction, ABVparam> energy0( "energy0" , 1 , Obs_Energy_4D , phi0 , pot0 );
  Observable<vector<double>,wavefunction> wi_x4( "wi_x4", 1 ,Distribution_X4, phi0);
  Observable<vector<double>,wavefunction> wi_x1x2( "wi_x1x2", 1 , Obs_Projection_X3X4, phi0);
  Observable<vector<double>,wavefunction> wi_x3x4( "wi_x3x4", 1 , Obs_Projection_X1X2, phi0);
 
  Observable<vector<double>,wavefunction> wf_x4( "wf_x4", 1 ,Distribution_X4, phi);
  Observable<vector<double>,wavefunction> wf_x1x2( "wf_x1x2", 1 , Obs_Projection_X3X4, phi);
  Observable<vector<double>,wavefunction> wf_x1x3( "wf_x1x3", 1 , Obs_Projection_X2X4, phi);
  Observable<vector<double>,wavefunction> wf_x2x4( "wf_x2x4", 1 , Obs_Projection_X1X3, phi);
  Observable<vector<double>,wavefunction> wf_x3x4( "wf_x3x4", 1 , Obs_Projection_X1X2, phi);


  Hamiltonian H0(phi0);
  
  
  double tmp_eng=0.;
  double eng_new=100., eng_old=0.;
  double time=0.;
  
  do
    {
      //elapsed_time=MPI_Wtime();
      
      time=time+1;
      eng_old=tmp_eng;
      H0.X1( -I*dt()*0.5, phi0 ,pot0, nthreads());
      H0.X2( -I*dt()*0.5, phi0 ,pot0, nthreads());
      H0.X3( -I*dt()*0.5, phi0 ,pot0, nthreads());
      H0.X4_MPI( -I*dt(), phi0 ,pot0, nlocal0, nthreads());
      H0.X3( -I*dt()*0.5, phi0 ,pot0, nthreads());
      H0.X2( -I*dt()*0.5, phi0 ,pot0, nthreads());
      H0.X1( -I*dt()*0.5, phi0 ,pot0, nthreads());
      Normalize(phi0);
     
      tmp_eng = energy0.measure();
      eng_new=tmp_eng;
      cout<<"Energy "<<tmp_eng<<endl;
      
      //elapsed_time=MPI_Wtime()-elapsed_time;
      //cout << " elapsed_time is : " <<  elapsed_time << endl;
    }
   while( abs(eng_new-eng_old) > GroundStateEnergyThreshold() );
   phi0.save(GroundStateName);

   wi_x4.record(time);
   wi_x1x2.record(time);
   wi_x3x4.record(time);

   MPI_Finalize();
   exit(1);

  // //////////************************************* laser begin  
   Hamiltonian H(phi); 
   phi0.load(GroundStateName);
   PlaceWaveFunction(phi,phi0); 


  // ////********************************************* laser parameters
   Sin2Laser laser1( wavelength_nm_atto(), intensity_wcm2_atto(), dt(), N_cycles_atto(), cep_atto(), epsilon_atto(), pol_angle_atto()*180./pi);
  Sin2Laser laser2( wavelength_nm(), intensity_wcm2(), dt(), N_cycles(),cep(), epsilon(), pol_angle()*180./pi );
  //GaussLaser laser1( wavelength_nm_atto(), intensity_wcm2_atto(), dt(), fwhm_cycles_atto(), cep_atto(), epsilon(), 0.);   
  //GaussLaser laser2( wavelength_nm(), intensity_wcm2(), dt(), fwhm_cycles(),cep(), epsilon(), 0.);

  Field field( laser1, 0., 0.);
  //Field field( laser1, laser2, delay(), 0., 0.);
  //Field field( laser2, 0., 0.);

  ofstream out5((outputFolder()+"/laserx_E.txt").c_str());
  ofstream out6((outputFolder()+"/laserx_A.txt").c_str());
  //ofstream out_time((outputFolder()+"/time_elapsed.txt").c_str());

  for (int it=0; it<field.t.size(); it++)
    {   
      out5<< it*dt() << "\t" << field.Ex[it] << "\t" << field.Ey[it] << endl;
      out6<< it*dt() << "\t" << field.Ax[it] << "\t" << field.Ay[it] << endl;
    }
  
  // //////////********************************************* end of  laser parameters 


   
   ////////////************************************* Halmitonian with laser now

   for(int time_index = 0; time_index < int(field.t.size()+pro_again/dt() ); time_index++ )
     {
       time=time_index*dt();
       if (time_index>=0 & time_index < field.t.size() )
       	 {
	   H.X1_Laser( dt()*0.5, phi , pot, field.Ax[time_index], velocitygauge, nthreads());
	   H.X2( dt()*0.5, phi ,pot, nthreads());
	   H.X3_Laser( dt()*0.5, phi , pot, field.Ax[time_index], velocitygauge, nthreads());
	   H.X4_MPI( dt(), phi ,pot, nlocal, nthreads());
	   //H.X4_Laser_MPI( dt(), phi, pot, field.Ax[time_index], velocitygauge, nlocal, nthreads());
	   H.X3_Laser( dt()*0.5, phi , pot, field.Ax[time_index], velocitygauge, nthreads());
	   H.X2( dt()*0.5, phi ,pot, nthreads());
	   H.X1_Laser( dt()*0.5, phi , pot, field.Ax[time_index], velocitygauge, nthreads());
       	 }
       else
       	 {
 	   H.X1( dt()*0.5, phi ,pot, nthreads());
	   H.X2( dt()*0.5, phi ,pot, nthreads());
	   H.X3( dt()*0.5, phi ,pot, nthreads());
	   H.X4_MPI( dt(), phi ,pot, nlocal, nthreads());
	   H.X3( dt()*0.5, phi ,pot, nthreads());
	   H.X2( dt()*0.5, phi ,pot, nthreads());
	   H.X1( dt()*0.5, phi ,pot, nthreads());
       	 }
       if(time_index % 5 == 0)
	 {
	   Mask_Function(phi, 0.1, 0.1, 0.1, 0.8, 0.1, 0.1, 0.1, 0.8, 0.167, nthreads());
	 }
      
   //     ////-----------------------------------------
       if(time_index % store_sampler == 0)
       	 {
	   wf_x4.record(time);
	   wf_x1x2.record(time);
	   wf_x3x4.record(time);
	   wf_x1x3.record(time);
	   wf_x2x4.record(time);	   
   	 }
     }
 

  MPI_Finalize();

} //end of main
   
