// $Id:$
#ifndef CARTESIAN4D_H
#define CARTESIAN4D_H
#include "constant.h"
#include "wavefunction.h"
#include "potentials.h"
#include "Laser.h"
#include <fftw3-mpi.h>

using namespace std; // to save the std:: in front of vector and cout

namespace Cartesian_4D{
  
  /**
   * Declare everything what is needed for the Crank-Nicholson-Scheme.  
   * independent of geometry namespace
   */

  struct temp_data_t {
    /**
     * temporary storage for Tridag
     */
    vector<complex> gam;
    
    /**
     * the 3 diagonals
     */
    vector<complex> tridag_upp;
    vector<complex> tridag_mid;
    vector<complex> tridag_low;
    
    /**
     * the potential with respect to the coordinate
     */
    vector<complex> v_1D;
    vector<complex> wf_1D;
    vector<complex> wf_1D_rightside;
    vector<complex> wf_1D_solution;
    
  };
  
  /**
   *  The namespace Cartesian deals only Hamiltonians H=Dxx + V.
   */
  
  
  /*--------------------------  Hamiltonian class  ----------------------------*/
  
  class Hamiltonian {
  public:
    
  /**
   * The crucial grid parameters
   * (introduced here to check if wf dimensions match grid dimensions)
   */
    double dx[DIM+1];/**< Spatial Step x[j]( i ) */  
    int n[DIM+1];/**< Number of points in the N[j] dimension */  

    temp_data_t temp_data[DIM+1]; /* provide workspace for DIM dimensions */
    
    /**
     * Constructor just Initializes the workspace
     */
    Hamiltonian( const wavefunction &wf );
    
    /**
     * Propagation Operators
     */
    
    void X1( const complex time_step, wavefunction &wf, const ABVparam &p, const int nthreads );
    void X2( const complex time_step, wavefunction &wf, const ABVparam &p, const int nthreads );
    void X3( const complex time_step, wavefunction &wf, const ABVparam &p, const int nthreads );
    void X4_MPI( const complex time_step , wavefunction &wf , const ABVparam &p, int nlocal, const int nthreads );
      
    void X1_Laser( const complex time_step, wavefunction &wf, const ABVparam &p,
		   const double field, const gauge_t gauge, const int nthreads );
    void X2_Laser( const complex time_step, wavefunction &wf, const ABVparam &p,
		   const double field, const gauge_t gauge, const int nthreads );
    void X3_Laser( const complex time_step, wavefunction &wf, const ABVparam &p,
		   const double field, const gauge_t gauge, const int nthreads );
    void X4_Laser_MPI( const complex time_step , wavefunction &wf , const ABVparam &p , const double field , const gauge_t gauge, int nlocal, const int nthreads );
    
    vector<complex> Inverse_X1_tau(const complex time_step , wavefunction &wf , const ABVparam &p);
    vector<complex> Inverse_X1_alpha(const complex time_step , wavefunction &wf , const ABVparam &p);
    vector<complex> Inverse_X1_gama(const complex time_step , wavefunction &wf , const ABVparam &p);

    void X1_inve( const complex time_step , wavefunction &wf , const ABVparam &p, vector<complex> &tau_load, vector<complex> &alpha_load, vector<complex> &gama_load, const int nthreads );

  }; // end of class Hamiltonian
  
  
 
  /*------------------	Wavefunction related functions	-----------------------*/
  
  // Attention !!!!!
  
	
  /**
   * A 3D initializer (void).  
   * Simply allocate the arrays and instanciate the number of points and the spatial step.   
   * Allocate the wavefunction (wave) as an array of size (n1+2)*(n2+2)*(n3+2)*(n4+2).   
   * Allocate the potential (v) as an array of size (n1+2)*(n2+2)*(n3+2)*(n4+2).	 
   * Allocate the x1 (x1) as an array of (n1+2).	 
   */
  
  void Initialize( wavefunction &wf , vector<int> n_points , vector<double> spatial_steps , int symmetry_x1, int symmetry_x2,int symmetry_x3, int symmetry_x4 );
  void Initialize_Momentum( wavefunction &wf , wavefunction &wf_mom);
  void Guess_Function_Gauss( wavefunction &wf , vector<double> sigmas , vector<double> gauss_centers , vector<double> initial_momenta );
  void Normalize( wavefunction &wf );
  
  
  void PlaceWaveFunction( wavefunction &wf , wavefunction &small_wf ); 
  complex Project( wavefunction &wf , wavefunction &wf2 );
  double Overlap( wavefunction &wf , wavefunction &wf2);
  void ProjectOUT_Diff_Sizes( wavefunction &wf , wavefunction &small_wf );
  
  
  void dig_h2p_e3d_r1d(wavefunction &wf, double radius, wavefunction &wf_new);
  void dig_h2p_e2d_r2d(wavefunction &wf, double radius,wavefunction &wf_new);
  void dig_h2p_e2d_r1d_ionz(wavefunction &wf, double radius,wavefunction &wf_new);
  void Rotate(wavefunction &tmp1, wavefunction &tmp2, wavefunction &xg, wavefunction &xu, vector<double> cp, const double dt, const double field );



 
  /**----------------   * Calculation of Observables:--------------*/
  double Obs_Norm( wavefunction &wf );
  double Obs_Norm_Squared( wavefunction &wf );
  double Obs_Energy_1D( wavefunction &wf, ABVparam &p );
  double Obs_Energy_2D( wavefunction &wf, ABVparam &p);
  double Obs_Energy_3D( wavefunction &wf, ABVparam &p);
  double Obs_Energy_4D( wavefunction &wf, ABVparam &p);


  double Obs_Expectation_Value_X1( wavefunction &wf );
  double Obs_Expectation_Value_Abs_X1( wavefunction &wf );
  double Obs_Expectation_Value_X2( wavefunction &wf );
  double Obs_Expectation_Value_X3( wavefunction &wf );
  double Obs_Expectation_Value_X4( wavefunction &wf );

  double Obs_Alignment( wavefunction &wf );

  vector<double> Obs_Fix_X1 (wavefunction &wf, double x1);
  vector<double> Obs_Fix_X2 (wavefunction &wf, double x2);
  vector<double> Obs_Fix_X3 (wavefunction &wf, double x3);
  vector<double> Obs_Projection_X1 (wavefunction &wf);
  vector<double> Obs_Projection_X2 (wavefunction &wf);
  vector<double> Obs_Projection_X3 (wavefunction &wf);
  vector<double> Obs_Projection_X4 (wavefunction &wf);
  vector<double> Obs_Projection_X1X2 (wavefunction &wf);
  vector<double> Obs_Projection_X1X3 (wavefunction &wf);
  vector<double> Obs_Projection_X1X4 (wavefunction &wf);
  vector<double> Obs_Projection_X2X3 (wavefunction &wf);
  vector<double> Obs_Projection_X2X4 (wavefunction &wf);
  vector<double> Obs_Projection_X3X4 (wavefunction &wf);

  vector<double> Distribution_X1( wavefunction &wf );
  vector<double> Distribution_X2( wavefunction &wf );
  vector<double> Distribution_X3( wavefunction &wf );
  vector<double> Distribution_X4( wavefunction &wf );




  
  /**----------   * Calculation of Derivatives: ----------------------------*/
  complex First_Derivative_X1( wavefunction &wf , int l, int k , int j , int i );
  complex First_Derivative_X2( wavefunction &wf , int l, int k , int j , int i );
  complex First_Derivative_X3( wavefunction &wf , int l, int k , int j , int i );
  complex Second_Derivative_X1( wavefunction &wf , int l, int k , int j , int i );
  complex Second_Derivative_X2( wavefunction &wf , int l, int k , int j , int i );
  complex Second_Derivative_X3( wavefunction &wf , int l, int k , int j , int i );
  complex Second_Derivative_X4_MPI( wavefunction &wf , int l, int k , int j , int i );

/*------------------  Stuff needed for FFT and Masks  -----------------------*/


  void Mask_Function(wavefunction &wf, double frac_n1_right, double frac_n2_right, double frac_n3_right, double frac_n4_right, double frac_n1_left, double frac_n2_left, double frac_n3_left, double frac_n4_left, double exponent, const int nthreads );

 
  void FFT_n1( wavefunction &wf , wavefunction &wf_transformed);
  void FFT_n1_shift( wavefunction &wf , wavefunction &wf_transformed);
  void FFT_n2( wavefunction &wf , wavefunction &wf_transformed);
  void FFT_n2_shift( wavefunction &wf , wavefunction &wf_transformed);
  void FFT_n3( wavefunction &wf , wavefunction &wf_transformed);
  void FFT_n3_shift( wavefunction &wf , wavefunction &wf_transformed);
  void FFT_n3n4( wavefunction &wf , wavefunction &wf_transformed);

  void fft_correlate_shift(vector<complex> &auto_relate, string outputFolder, double dt);





/*----------------  Classical Nuclear Motion functions  ---------------------*/


  
 }

#endif	/* CARTESIAN4D_H */

