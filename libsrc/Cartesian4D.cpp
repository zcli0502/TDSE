#include <fftw3-mpi.h>
#include "constant.h"
#include "Cartesian4D.h"
#include "utils.h"
#include "potentials.h"
#include <math.h>
#include <complex>
#define complex complex<double>
#include "ParameterMap.h" 
#include <omp.h>


#include <time.h>
using namespace std;
namespace Cartesian_4D{

/*--------------------------  Hamiltonian class  ----------------------------*/

Hamiltonian::Hamiltonian(const wavefunction &wf) {
  /**
   *  Note the grid dimensions that were used to initialize the
   * ABV and temp data vectors
   */
  for (int coord_index=1; coord_index<=4; ++coord_index) {
    n[coord_index]  = wf.n[coord_index];
    dx[coord_index] = wf.x[coord_index].delta;
  }

  /**
   * Allocate the temporary storage for Tridag
   */
  for ( int coord_index=1; coord_index <= 4; ++coord_index ) {
    int num = wf.n[coord_index] + 2;
    temp_data_t &d=temp_data[coord_index];

    d.gam.resize(num, complex(0., 0.) );
    d.tridag_low.resize(num, complex(0., 0.) );
    d.tridag_mid.resize(num, complex(0., 0.) );
    d.tridag_upp.resize(num, complex(0., 0.) );
    d.v_1D.resize(num, complex(0., 0.) );
    d.wf_1D.resize(num, complex(0., 0.) );
    d.wf_1D_rightside.resize(num, complex(0., 0.) );
    d.wf_1D_solution.resize(num, complex(0., 0.) );
    
  }
  
 
 

};

  /**
   * One time step in Real time (NO LASER), the default value for a in a*(d/dx^2) is a=1/2
   * IMPORTANT, the time step used comes from the Cranck-Nichols equation
   * The auxiliary vector are allocated and initialized inside the function, it uses the tridag (www.nr.com) function implemented in utils.h 
   */  




  void X1_kernel(const complex time_step , wavefunction &wf , const ABVparam &p, temp_data_t d, int l, int k, int j)
  {
    const complex idt=I*time_step;
    const complex arg_A = ( idt*.5*wf.one_by_dx1sqr )*p.ABV_A_x1;
    const complex arg_B = 0.;//complex(0.,0.);//( idt*.25*wf.one_by_dx1   )*p.ABV_B_x1;
    complex arg_V;
    complex tridag_low_Fast =arg_A-arg_B;
    complex tridag_upp_Fast =arg_A+arg_B;

    
    for( int i = 0 ; i < wf.n1+2 ; ++i )//Copy the wavefunction
      {
	/**
	 * Creates a 1D potential and a 1D wavefunction out of the 3D ones
	 */
        const int index=wf.in4(l, k , j , i );
	d.v_1D[ i ] = p.ABV_V[ index ];
	d.wf_1D[ i ] = wf.wave[ index ];
      }
    
    
    for( int i = 1 ; i <= wf.n1 ; ++i )
      {
	/**
	 * Definition of the arguments in Eq. (7-9) in doc/ABV/abv.tex
	 */
	
	arg_V = ( idt*.5 )*d.v_1D[ i ];
		
	/**
	 * Define the 3 diagonals X of the left side of Eq. (6)  in doc/ABV/abv.tex => look at Eqs. (7-9) in doc/ABV/abv.tex
	 */
	d.tridag_mid[ i ] =1.-2.*arg_A+arg_V;
	
	d.wf_1D_rightside[ i ]=-tridag_low_Fast*d.wf_1D[ i-1 ]+( 1.+2.*arg_A-arg_V )*d.wf_1D[ i ] -tridag_upp_Fast*d.wf_1D[ i+1 ];//general case, symmetry==1
	
	/**
	 * for (anti-)symmetric wavefunctions one has to take care of the right boundary conditions (phi(0)=+/-phi(1). That leads to modified values.
	 */
	if (i==1) {
	  if (wf.symmetry_x1==1)
	    {
	      d.tridag_mid[ i ] =1.-arg_A-arg_B+arg_V;
	      d.wf_1D_rightside[ i ]=-tridag_low_Fast*d.wf_1D[ i-1 ] +( 1.+arg_A+arg_B-arg_V )*d.wf_1D[ i ]-tridag_upp_Fast*d.wf_1D[ i+1 ];
	    }
	  
	  if (wf.symmetry_x1==-1)
	    {
		    d.tridag_mid[ i ] =1.-3.*arg_A+arg_B+arg_V;
		    d.wf_1D_rightside[ i ]=-tridag_low_Fast*d.wf_1D[ i-1 ] +( 1.+3.*arg_A-arg_B-arg_V )*d.wf_1D[ i ]-tridag_upp_Fast*d.wf_1D[ i+1 ];
	    }
	}
	
      }
    
    /**
	     * Solve the system of linear equation of Eq. (6) in doc/ABV/abv.tex
	     */
    Tridag_Fast( tridag_low_Fast , d.tridag_mid, tridag_upp_Fast , d.wf_1D_rightside, d.wf_1D_solution, d.gam );
    
    /**
     * Copy back the 1D wavefunction to the 3D one
     */
    for( int i = 1 ; i <= wf.n1 ; i++ )
      {
	wf.wave[ wf.in4(l, k , j , i ) ] = d.wf_1D_solution[ i ];
      }
    
    
  }
  
  void X2_kernel(const complex time_step , wavefunction &wf , const ABVparam &p, temp_data_t d, int l, int k, int i)
  {
    const complex idt=I*time_step;
    const complex arg_A = idt*.5*wf.one_by_dx2sqr  *p.ABV_A_x2;
    const complex arg_B = 0;//( idt*.25*wf.one_by_dx2 )*p.ABV_B_x2;
    complex arg_V;
    complex tridag_upp_Fast;
    complex tridag_low_Fast;
    for( int j = 0 ; j < wf.n2+2 ; j++ )
      {
	const int index=wf.in4(l, k , j , i );
	d.v_1D[ j ]  = p.ABV_V[ index ];
	d.wf_1D[ j ] = wf.wave[ index ];
      }
    
    for( int j = 1 ; j <= wf.n2 ; j++ )
      {
	
	arg_V = ( idt*.5  )*d.v_1D[ j ];
	tridag_low_Fast    = arg_A-arg_B;	     
	d.tridag_mid[ j ] = 1.-2.*arg_A+arg_V;
	tridag_upp_Fast    = arg_A+arg_B;
	
	d.wf_1D_rightside[ j ]= -tridag_low_Fast*d.wf_1D[ j-1 ] + ( 1.+2.*arg_A-arg_V )*d.wf_1D[ j ] - tridag_upp_Fast*d.wf_1D[ j+1 ];//sym==1
	
	
	if (j==1) {
	  if (wf.symmetry_x2==1)
	    {
	      d.tridag_mid[ j ] = 1.-arg_A-arg_B+arg_V;
	      d.wf_1D_rightside[ j ]=-tridag_low_Fast*d.wf_1D[ j-1 ] +( 1.+arg_A+arg_B-arg_V )*d.wf_1D[ j ]-tridag_upp_Fast*d.wf_1D[ j+1 ];
	      
	    }
	  
	  if (wf.symmetry_x2==-1)
	    {
	      d.tridag_mid[ j ] = 1.-3.*arg_A+arg_B+arg_V;
	      d.wf_1D_rightside[ j ]=-tridag_low_Fast*d.wf_1D[ j-1 ] +( 1.+3.*arg_A-arg_B-arg_V )*d.wf_1D[ j ]-tridag_upp_Fast*d.wf_1D[ j+1 ];
	      
	    }
	}
	
	
      }
    
    Tridag_Fast( tridag_low_Fast , d.tridag_mid, tridag_upp_Fast , d.wf_1D_rightside, d.wf_1D_solution, d.gam );
    for( int j = 1 ; j <= wf.n2 ; j++ )
      {
	wf.wave[ wf.in4( l, k , j , i ) ] = d.wf_1D_solution[ j ];
      }
    
  }



  void X3_kernel(const complex time_step , wavefunction &wf , const ABVparam &p, temp_data_t d, int l, int j, int i)
  {
    const complex idt=I*time_step;
    const complex arg_A = idt*.5*wf.one_by_dx3sqr    *p.ABV_A_x3;
    const complex arg_B =0;// idt*.25*wf.one_by_dx3      *p.ABV_B_x3;
    complex arg_V;
    complex tridag_upp_Fast;
    complex tridag_low_Fast;
    
    for( int k = 0 ; k < wf.n3+2 ; k++ )
      {
  	const int index=wf.in4(l, k , j , i );
  	d.v_1D[ k ]  = p.ABV_V[ index ];
  	d.wf_1D[ k ] = wf.wave[ index ];
      }
    
    for( int k = 1 ; k <= wf.n3 ; k++ )
      {		
  	arg_V = idt*.5                   *d.v_1D[ k ];
	
  	tridag_low_Fast    = arg_A-arg_B;	     
  	d.tridag_mid[ k ] = 1.-2.*arg_A+arg_V;
  	tridag_upp_Fast    = arg_A+arg_B;	     
	
  	d.wf_1D_rightside[ k ]= -tridag_low_Fast*d.wf_1D[ k-1 ] +( 1.+2.*arg_A-arg_V )*d.wf_1D[ k ] -tridag_upp_Fast*d.wf_1D[ k+1 ];//sym==1
	
  	if (k==1) {
  	  if (wf.symmetry_x3==1)
  	    {
  	      d.tridag_mid[ k ] = 1.-arg_A-arg_B+arg_V;
  	      d.wf_1D_rightside[ k ]= -tridag_low_Fast*d.wf_1D[ k-1 ] +( 1.+arg_A+arg_B-arg_V )*d.wf_1D[ k ] -tridag_upp_Fast*d.wf_1D[ k+1 ];
  	    }
	  
  	  if (wf.symmetry_x3==-1)
  	    {
  	      d.tridag_mid[ k ] = 1.-3.*arg_A+arg_B+arg_V;
  	      d.wf_1D_rightside[ k ]= -tridag_low_Fast*d.wf_1D[ k-1 ] +( 1.+3.*arg_A-arg_B-arg_V )*d.wf_1D[ k ] -tridag_upp_Fast*d.wf_1D[ k+1 ];
  	    }
  	}
	
      }
    
    Tridag_Fast( tridag_low_Fast , d.tridag_mid, tridag_upp_Fast , d.wf_1D_rightside, d.wf_1D_solution, d.gam );
    
    
    for( int k = 1 ; k <= wf.n3 ; k++ )
      {
    	wf.wave[ wf.in4(l, k , j , i ) ] = d.wf_1D_solution[ k ];
      }
  }
   
  

  void X1_Laser_kernel(const complex time_step , wavefunction &wf, const ABVparam &p, const double field, const gauge_t gauge,  temp_data_t d, int l, int k, int j )
  {
    const complex idt=I*time_step;
    const complex arg_A = idt*.5*wf.one_by_dx1sqr *p.ABV_A_x1;
    const complex f_arg_B=time_step*.25*wf.one_by_dx1; 
    complex arg_B;
    complex arg_V;
    complex tridag_upp_Fast;
    complex tridag_low_Fast;

    for( int i = 0 ; i < wf.n1+2 ; i++ )
      {
	const int index=wf.in4( l, k , j , i );
	d.v_1D[ i ]  = p.ABV_V[ index ];
	d.wf_1D[ i ] = wf.wave[ index ];
      }
    
    for( int i = 1 ; i <= wf.n1 ; i++ )
      {
	
	if ( gauge == lengthgauge )
	  {
	    arg_B = 0;//f_arg_B       *p.ABV_B_x2;
	    arg_V = idt*.5                  *( p.ABV_B_x1/p.ABV_A_x1*0.5*wf.x1[ i ]*field + d.v_1D[ i ] ); //p.ABV_A_x2*0.5 is coming from /doc/ABV/factors
	  }
	else if ( gauge == velocitygauge )
	  {
	    arg_B = f_arg_B*( -p.ABV_B_x1*field*one_by_lightC_au );
	    arg_V = idt*.5                  *d.v_1D[ i ];
	  }
	//	else
	// {
	//   arg_B = 0.;
	//    arg_V = idt*.5                  *d.v_1D[ j ];
	//  }
	
	tridag_low_Fast    = arg_A-arg_B;	     
	d.tridag_mid[ i ] = 1.-2.*arg_A+arg_V;
	tridag_upp_Fast    = arg_A+arg_B;	     
		
	d.wf_1D_rightside[ i ]= -tridag_low_Fast*d.wf_1D[ i-1 ] +( 1.+2.*arg_A-arg_V )*d.wf_1D[ i ] -tridag_upp_Fast*d.wf_1D[ i+1 ];//sym==1
	
	if (i==1) {
	  if (wf.symmetry_x1==1)
	    {
	      d.tridag_mid[ i ] = 1.-arg_A-arg_B+arg_V;
	      d.wf_1D_rightside[ i ]= -tridag_low_Fast*d.wf_1D[ i-1 ] +( 1.+arg_A+arg_B-arg_V )*d.wf_1D[ i ] -tridag_upp_Fast*d.wf_1D[ i+1 ];
	    }
	  
	  if (wf.symmetry_x1==-1)
	    {
	      d.tridag_mid[ i ] = 1.-3.*arg_A+arg_B+arg_V;
	      d.wf_1D_rightside[ i ]= -tridag_low_Fast*d.wf_1D[ i-1 ] +( 1.+3.*arg_A-arg_B-arg_V )*d.wf_1D[ i ] -tridag_upp_Fast*d.wf_1D[ i+1 ];
	    }
	}
	
      }
    
    Tridag_Fast( tridag_low_Fast , d.tridag_mid, tridag_upp_Fast , d.wf_1D_rightside, d.wf_1D_solution, d.gam );
    
	    
    for( int i = 1 ; i <= wf.n1 ; i++ )
      {
	wf.wave[ wf.in4(l, k , j , i ) ] = d.wf_1D_solution[ i ];
      } 
  }


  void X2_Laser_kernel(const complex time_step , wavefunction &wf, const ABVparam &p, const double field, const gauge_t gauge,  temp_data_t d, int l, int k, int i )
 {
   
   const complex idt=I*time_step;
   const complex arg_A = idt*.5*wf.one_by_dx2sqr *p.ABV_A_x2;
   const complex f_arg_B=time_step*.25*wf.one_by_dx2; 
   complex arg_B;
   complex arg_V;
   complex tridag_upp_Fast;
   complex tridag_low_Fast;
   
   for( int j = 0 ; j < wf.n2+2 ; j++ )
     {
       const int index=wf.in4(l, k , j , i );
       d.v_1D[ j ]  = p.ABV_V[ index ];
       d.wf_1D[ j ] = wf.wave[ index ];
     }
   
   for( int j = 1 ; j <= wf.n2 ; j++ )
     {
       
       if ( gauge == lengthgauge )
	 {
	   arg_B = 0;//f_arg_B       *p.ABV_B_x2;
	   arg_V = idt*.5                  *( p.ABV_B_x2/p.ABV_A_x2*0.5*wf.x2[ j ]*field + d.v_1D[ j ] ); //p.ABV_A_x2*0.5 is coming from /doc/ABV/factors
	 }
       else if ( gauge == velocitygauge )
	 {
	   arg_B = f_arg_B*( -p.ABV_B_x2*field*one_by_lightC_au );
	   arg_V = idt*.5                  *d.v_1D[ j ];
	 }
       //	else
       // {
       //   arg_B = 0.;
       //    arg_V = idt*.5                  *d.v_1D[ j ];
       //  }
       
       tridag_low_Fast    = arg_A-arg_B;	     
       d.tridag_mid[ j ] = 1.-2.*arg_A+arg_V;
       tridag_upp_Fast    = arg_A+arg_B;	     
       
       d.wf_1D_rightside[ j ]= -tridag_low_Fast*d.wf_1D[ j-1 ] +( 1.+2.*arg_A-arg_V )*d.wf_1D[ j ] -tridag_upp_Fast*d.wf_1D[ j+1 ];//sym==1

       if (j==1) {
	 if (wf.symmetry_x2==1)
	   {
	     d.tridag_mid[ j ] = 1.-arg_A-arg_B+arg_V;
	     d.wf_1D_rightside[ j ]= -tridag_low_Fast*d.wf_1D[ j-1 ] +( 1.+arg_A+arg_B-arg_V )*d.wf_1D[ j ] -tridag_upp_Fast*d.wf_1D[ j+1 ];
	   }
	 
	 if (wf.symmetry_x2==-1)
	   {
	     d.tridag_mid[ j ] = 1.-3.*arg_A+arg_B+arg_V;
	     d.wf_1D_rightside[ j ]= -tridag_low_Fast*d.wf_1D[ j-1 ] +( 1.+3.*arg_A-arg_B-arg_V )*d.wf_1D[ j ] -tridag_upp_Fast*d.wf_1D[ j+1 ];
	   }
       }
		
     }
   
   Tridag_Fast( tridag_low_Fast , d.tridag_mid, tridag_upp_Fast , d.wf_1D_rightside, d.wf_1D_solution, d.gam );
   
   
   for( int j = 1 ; j <= wf.n2 ; j++ )
     {
       wf.wave[ wf.in4(l, k , j , i ) ] = d.wf_1D_solution[ j ];
     }   
 }


  void X3_Laser_kernel(const complex time_step , wavefunction &wf, const ABVparam &p, const double field, const gauge_t gauge,  temp_data_t d, int l, int j, int i )
 {
    const complex idt=I*time_step;
    const complex arg_A = idt*.5*wf.one_by_dx3sqr *p.ABV_A_x3;
    const complex f_arg_B=time_step*.25*wf.one_by_dx3;
    complex arg_B;
    complex arg_V;
    complex tridag_upp_Fast;
    complex tridag_low_Fast;

    for( int k = 0 ; k < wf.n3+2 ; k++ )
      {
	const int index=wf.in4(l, k , j , i );
	d.v_1D[ k ]  = p.ABV_V[ index ];
	d.wf_1D[ k ] = wf.wave[ index ];
      }
    for( int k = 1 ; k <= wf.n3 ; k++ )
      {
	
	if ( gauge == lengthgauge )
	  {
	    arg_B = 0.;//f_arg_B       *p.ABV_B_x3;
	    arg_V = idt*.5                  *( p.ABV_B_x3/p.ABV_A_x3*0.5*wf.x3[ k ]*field + d.v_1D[ k ] ); //p.ABV_A_x3*0.5 is coming from /doc/ABV/factors
	  }
	else if ( gauge == velocitygauge )
	  {
	    arg_B = f_arg_B       *( -p.ABV_B_x3*field*one_by_lightC_au );
	    arg_V = idt*.5                  *d.v_1D[ k ];
	  }
	else
	  {
	    arg_B = 0.;
	    arg_V = idt*.5                  *d.v_1D[ k ];
	  }
	
	tridag_low_Fast    = arg_A-arg_B;	     
	d.tridag_mid[ k ] = 1.-2.*arg_A+arg_V;
	tridag_upp_Fast    = arg_A+arg_B;	     
	
	d.wf_1D_rightside[ k ]= -tridag_low_Fast*d.wf_1D[ k-1 ] +( 1.+2.*arg_A-arg_V )*d.wf_1D[ k ] -tridag_upp_Fast*d.wf_1D[ k+1 ];//sym==1
	
	if (k==1) {
	  if (wf.symmetry_x3==1)
	    {
	      d.tridag_mid[ k ] = 1.-arg_A-arg_B+arg_V;
	      d.wf_1D_rightside[ k ]= -tridag_low_Fast*d.wf_1D[ k-1 ] +( 1.+ arg_A+arg_B-arg_V )*d.wf_1D[ k ] -tridag_upp_Fast*d.wf_1D[ k+1 ];
	    }
	  
	  if (wf.symmetry_x3==-1)
	    {
	      d.tridag_mid[ k ] = 1.-3.*arg_A+arg_B+arg_V;
	      d.wf_1D_rightside[ k ]= -tridag_low_Fast*d.wf_1D[ k-1 ] +( 1.+ 3.*arg_A-arg_B-arg_V )*d.wf_1D[ k ] -tridag_upp_Fast*d.wf_1D[ k+1 ];
	    }
	}
	
      }
    
    Tridag_Fast( tridag_low_Fast , d.tridag_mid, tridag_upp_Fast , d.wf_1D_rightside, d.wf_1D_solution, d.gam );
    
    
    for( int k = 1 ; k <= wf.n3 ; k++ )
      {
	wf.wave[ wf.in4( l,  k , j , i ) ] = d.wf_1D_solution[ k ];
      }
 }
  



  void Hamiltonian::X1( const complex time_step , wavefunction &wf , const ABVparam &p,const int nthreads )
  {
    temp_data_t &d=temp_data[1]; 
    if(n[1] != wf.n1 || dx[1] != wf.dx1) {
      cerr<<"wavefunction has wrong size for Hamiltonian\n";
      exit(1);
    }
    
    // int nthreads=16;
    // double wall_timer=omp_get_wtime();
    // clock_t clock_timer=clock();
    //int k,j;
    #pragma omp parallel for  num_threads(nthreads) collapse(3)   
    for(int l = 1 ; l <= wf.n4 ; l++ )
      for(int k = 1 ; k <= wf.n3 ; k++ )
	{   
	  for(int j = 1 ; j <= wf.n2 ; j++ )
	    {
	   
	      X1_kernel(time_step, wf , p, d, l, k, j);
	    }
	}
    //std::cout<< "threads:"<< nthreads <<"time on clock():"<<(double) (clock() -clock_timer)/ CLOCKS_PER_SEC<<" time on wall :"<<omp_get_wtime()-wall_timer<<endl;
  } // end of X1_3D 

  void Hamiltonian::X2(const complex time_step , wavefunction &wf , const ABVparam &p,const int nthreads )
  {
    temp_data_t &d=temp_data[2]; 
    if(n[2] != wf.n2 || dx[2] != wf.dx2) {
      cerr<<"wavefunction has wrong size for Hamiltonian\n";
      exit(1);
    }


    //  int nthreads=16;
    //int k,i;
#pragma omp parallel for  num_threads(nthreads) collapse(3)
    for(int l = 1 ; l <= wf.n4 ; l++ )
      for(int k = 1 ; k <= wf.n3 ; k++ )
	{
	  for(int i = 1 ; i <= wf.n1 ; i++ )
	    {
	      X2_kernel( time_step, wf, p, d, l, k, i);
	    }
	}
  } // end of X2_3D



  void Hamiltonian::X3( const complex time_step , wavefunction &wf , const ABVparam &p,const int nthreads )
  {
    temp_data_t &d=temp_data[3]; 

    if(n[3] != wf.n3 || dx[3] != wf.dx3) {
      cerr<<"wavefunction has wrong size for Hamiltonian\n";
      exit(1);
    }
    //  int nthreads=16;
   #pragma omp parallel for  num_threads(nthreads) collapse(3)    
    for( int l = 1 ; l <= wf.n4 ; l++ )
      for( int j = 1 ; j <= wf.n2 ; j++ )
	{
	  for( int i = 1 ; i <= wf.n1 ; i++ )
	    {
	      X3_kernel( time_step , wf , p, d, l, j, i);
	    }
	}
  } // end of X3_3D




   void Hamiltonian::X4_MPI( const complex time_step , wavefunction &wf , const ABVparam &p, const int nlocal, const int nthreads)
   {
     if(n[4] != wf.n4 || dx[4] != wf.dx4) {
       cerr<<"wavefunction has wrong size for Hamiltonian\n";
       exit(1);
     }
  
     const complex idt=I*time_step;
     const complex arg_A = idt*.5*wf.one_by_dx4sqr    *p.ABV_A_x4;
     const complex arg_B =0.;// idt*.25*wf.one_by_dx3      *p.ABV_B_x3;
     complex arg_V;
     complex tridag_upp_Fast=arg_A+arg_B;
     complex tridag_low_Fast=arg_A-arg_B;
     int tag1=3, tag2=4, tag3=5, tag4=6;
     int left, right, myid, numprocs;
     int num_comm=(wf.n3+2)*(wf.n2+2)*(wf.n1+2);
     int num_data=(wf.n1+2)*(wf.n2+2)*(wf.n3+2)*(wf.n4+2);
     MPI_Status status;
     MPI_Comm_rank (MPI_COMM_WORLD, &myid);
     MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
     vector<complex> tridag_mid(num_data,0.);     
     vector<complex> wf_rightside(num_data,0.);     

    if ( myid > 0 )
      left=myid-1;
    else 
      left=MPI_PROC_NULL;
    if ( myid < numprocs-1  )
      right=myid+1;
    else
      right=MPI_PROC_NULL;  //// define MPI_NULL buffer
    vector<complex> wf_left_send(num_comm,0.), wf_right_send(num_comm,0.);
    vector<complex> wf_left_recv(num_comm,0.), wf_right_recv(num_comm,0.);

#pragma omp parallel for  num_threads(nthreads)
     for( int k = 1 ; k <= wf.n3 ; k++ )
     for( int j = 1 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1 ; i++ )
   	  {
   	    wf_left_send[wf.in3(k,j,i)]=wf.wave[wf.in4(1,k,j,i)];
            wf_right_send[wf.in3(k,j,i)]=wf.wave[wf.in4(nlocal,k,j,i)];
   	  }
       }
     //cout<<"num2d"<<num_2D<<" num3d"<<num_3D<<"nlocal"<<nlocal<<endl;
     //MPI_Send(&wf_left_send[0],num_2D,MPI_DOUBLE_COMPLEX,left,tag1,MPI_COMM_WORLD);
     //cout<<"send left"<<myid<<endl;
     //MPI_Recv(&wf_right_recv[0],num_2D,MPI_DOUBLE_COMPLEX,right,tag1,MPI_COMM_WORLD, &status);
     MPI_Sendrecv(&wf_left_send[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag1,&wf_right_recv[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag1,MPI_COMM_WORLD, &status);
     MPI_Sendrecv(&wf_right_send[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag2,&wf_left_recv[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag2,MPI_COMM_WORLD, &status);

     //cout<<"Send_wf_right"<<endl;

#pragma omp parallel for  num_threads(nthreads) private(arg_V) collapse(3)
     for( int k = 1 ; k <= wf.n3 ; k++ )
     for( int j = 1 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1 ; i++ )
   	  {
            temp_data_t d=temp_data[4]; 
   	    for( int l = 0 ; l < wf.n4+2 ; l++ )
   	      {
                 const int index=wf.in4( l, k , j , i );
   		d.v_1D[ l ]  = p.ABV_V[ index ];
   		d.wf_1D[ l ] = wf.wave[ index ];
   	      }
            d.wf_1D[0]=wf_left_recv[wf.in3(k,j,i)];
	    d.wf_1D[nlocal+1]=wf_right_recv[wf.in3(k,j,i)]; 
   	    for( int l = 1 ; l <= wf.n4 ; l++ )
   	      {		
   		arg_V = idt*.5*d.v_1D[ l ];
   		tridag_mid[wf.in4(l,k,j,i)] = 1.-2.*arg_A+arg_V;
      	
   		wf_rightside[wf.in4(l,k,j,i)]= -tridag_low_Fast*d.wf_1D[ l-1 ] +( 1.+2.*arg_A-arg_V )*d.wf_1D[ l ] -tridag_upp_Fast*d.wf_1D[ l+1 ];//sym==1

		if ( (l==1) && (myid==0) ) {
   		if (wf.symmetry_x4==1)
   		  {
   		    tridag_mid[wf.in4(l,k,j,i)] = 1.-arg_A-arg_B+arg_V;
   		    wf_rightside[wf.in4(l,k,j,i)]= -tridag_low_Fast*d.wf_1D[ l-1 ] +( 1.+arg_A+arg_B-arg_V )*d.wf_1D[ l ] -tridag_upp_Fast*d.wf_1D[ l+1 ];
   		  }
      	
   		if (wf.symmetry_x4==-1)
   		  {
   		    tridag_mid[wf.in4(l,k,j,i)] = 1.-3.*arg_A+arg_B+arg_V;
   		    wf_rightside[wf.in4(l,k,j,i)] = -tridag_low_Fast*d.wf_1D[ l-1 ] +( 1.+3.*arg_A-arg_B-arg_V )*d.wf_1D[ l ] -tridag_upp_Fast*d.wf_1D[ l+1 ];
   		  }
                 }
   	      }}}
     //cout<<"tridag_mpi_start"<<endl;       
          
     //Tridag_X4_MPI( wf, tridag_low_Fast , tridag_mid, tridag_upp_Fast , wf_rightside, nthreads);          
     Tridag_X4_MPI_opti( wf, tridag_low_Fast , tridag_mid, tridag_upp_Fast , wf_rightside, nthreads);  // the number for every dimention should be even.
     //cout<<"wf[wf.in4(1,1,1,1)]"<< wf.wave[wf.in4(1,1,1,1)]<<endl;         
   } // end of X4_4D_MPI









  void Hamiltonian::X1_Laser(const complex time_step , wavefunction &wf, const ABVparam &p, const double field, const gauge_t gauge, const int nthreads )
  {
    temp_data_t &d=temp_data[1]; 
    if(n[1] != wf.n1 || dx[1] != wf.dx1) {
      cerr<<"wavefunction has wrong size for Hamiltonian\n";
      exit(1);
    }
    // int nth;
    // nth=omp_get_num_threads();
    //   cout << "we get the num of threads: " << nth << endl;
    //  int nthreads=16;
   #pragma omp parallel for  num_threads(nthreads) collapse(2) 
    for( int l = 1 ; l <= wf.n4 ; l++ )
      for( int j = 1 ; j <= wf.n2 ; j++ )
	{
	  for( int k = 1 ; k <= wf.n3 ; k++ )
	    {
	      X1_Laser_kernel( time_step , wf, p, field, gauge, d,l, k, j );
	    }
	}
    
  } // end of X1_3D_Laser
  

  


  void Hamiltonian::X2_Laser( const complex time_step , wavefunction &wf , const ABVparam &p , const double field , const gauge_t gauge, const int nthreads)
  {
    temp_data_t &d=temp_data[2]; 
    if(n[2] != wf.n2 || dx[2] != wf.dx2) {
      cerr<<"wavefunction has wrong size for Hamiltonian\n";
      exit(1);
    }
    
    // int nthreads=16;
  #pragma omp parallel for  num_threads(nthreads) collapse(2)
    for( int l = 1 ; l <= wf.n4 ; l++ )
      for( int k = 1 ; k <= wf.n3 ; k++ )
	{
	  for( int i = 1 ; i <= wf.n1 ; i++ )
	    {
	      X2_Laser_kernel( time_step , wf, p, field, gauge, d, l, k, i );
	    }
	}
  } // end of X2_3D_Laser

  void Hamiltonian::X3_Laser( const complex time_step , wavefunction &wf , const ABVparam &p , const double field , const gauge_t gauge, const int nthreads )
  {
   temp_data_t &d=temp_data[3]; 
    if(n[3] != wf.n3 || dx[3] != wf.dx3) {
      cerr<<"wavefunction has wrong size for Hamiltonian\n";
      exit(1);
    }

    // int nthreads=16;
  #pragma omp parallel for  num_threads(nthreads) collapse(2)
    for( int l = 1 ; l <= wf.n4 ; l++ )
      for( int j = 1 ; j <= wf.n2 ; j++ )
	{
	  for( int i = 1 ; i <= wf.n1 ; i++ )
	    {
	      X3_Laser_kernel( time_step , wf, p, field, gauge, d,l, j, i );
	    }
	}
  } // end of X3_3D_Laser




   void Hamiltonian::X4_Laser_MPI( const complex time_step , wavefunction &wf , const ABVparam &p, const double field , const gauge_t gauge, const int nlocal, const int nthreads)
   {
     if(n[4] != wf.n4 || dx[4] != wf.dx4) {
       cerr<<"wavefunction has wrong size for Hamiltonian\n";
       exit(1);
     }
  
     const complex idt=I*time_step;
     const complex arg_A = idt*.5*wf.one_by_dx4sqr    *p.ABV_A_x4;
     const complex f_arg_B =time_step*.25*wf.one_by_dx4;
     complex arg_B, arg_V;
     complex tridag_upp_Fast;
     complex tridag_low_Fast;
     int tag1=3, tag2=4, tag3=5, tag4=6;
     int left, right, myid, numprocs;
     int num_comm=(wf.n3+2)*(wf.n2+2)*(wf.n1+2);
     int num_data=(wf.n1+2)*(wf.n2+2)*(wf.n3+2)*(wf.n4+2);
     MPI_Status status;
     MPI_Comm_rank (MPI_COMM_WORLD, &myid);
     MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
     vector<complex> tridag_mid(num_data,0);     
     vector<complex> wf_rightside(num_data,0);     

    if ( myid > 0 )
      left=myid-1;
    else 
      left=MPI_PROC_NULL;
    if ( myid < numprocs-1  )
      right=myid+1;
    else
      right=MPI_PROC_NULL;  //// define MPI_NULL buffer
    vector<complex> wf_left_send(num_comm,0.), wf_right_send(num_comm,0.);
    vector<complex> wf_left_recv(num_comm,0.), wf_right_recv(num_comm,0.);

#pragma omp parallel for  num_threads(nthreads)
     for( int k = 1 ; k <= wf.n3 ; k++ )
       for( int j = 1 ; j <= wf.n2 ; j++ )
	 {
	   for( int i = 1 ; i <= wf.n1 ; i++ )
	     {
	       wf_left_send[wf.in3(k,j,i)]=wf.wave[wf.in4(1,k,j,i)];
	       wf_right_send[wf.in3(k,j,i)]=wf.wave[wf.in4(nlocal,k,j,i)];
	     }
	 }

     MPI_Sendrecv(&wf_left_send[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag1,&wf_right_recv[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag1,MPI_COMM_WORLD, &status);
     MPI_Sendrecv(&wf_right_send[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag2,&wf_left_recv[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag2,MPI_COMM_WORLD, &status);


     if ( gauge == lengthgauge )
       {
	 arg_B = 0.;
       }
     else if ( gauge == velocitygauge )
       {
	 arg_B = f_arg_B       *( -p.ABV_B_x4*field*one_by_lightC_au );
       }
     else
       {
	 arg_B = 0.;
       }
     tridag_low_Fast    = arg_A-arg_B;	     
     tridag_upp_Fast    = arg_A+arg_B;
     

#pragma omp parallel for  num_threads(nthreads) private( arg_V ) collapse(3)
     for( int k = 1 ; k <= wf.n3 ; k++ )
     for( int j = 1 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1 ; i++ )
   	  {
            temp_data_t d=temp_data[4]; 
   	    for( int l = 0 ; l < wf.n4+2 ; l++ )
   	      {
		const int index=wf.in4( l, k , j , i );
   		d.v_1D[ l ]  = p.ABV_V[ index ];
   		d.wf_1D[ l ] = wf.wave[ index ];
   	      }
            d.wf_1D[0]=wf_left_recv[wf.in3(k,j,i)];
	    d.wf_1D[nlocal+1]=wf_right_recv[wf.in3(k,j,i)]; 
   	    for( int l = 1 ; l <= wf.n4 ; l++ )
   	      {		
   		if ( gauge == lengthgauge )
	        {
	         arg_V = idt*.5 *( p.ABV_B_x4/p.ABV_A_x4*0.5*wf.x4[ l ]*field + d.v_1D[ l ] ); 
	        }
	       else if ( gauge == velocitygauge )
	        {
	         arg_V = idt*.5                  *d.v_1D[ l ];
	        }
	       else
	        {
	         arg_V = idt*.5                  *d.v_1D[ l ];
	        }

   		tridag_mid[wf.in4(l,k,j,i)] = 1.-2.*arg_A+arg_V;
   		wf_rightside[wf.in4(l,k,j,i)]= -tridag_low_Fast*d.wf_1D[ l-1 ] +( 1.+2.*arg_A-arg_V )*d.wf_1D[ l ] -tridag_upp_Fast*d.wf_1D[ l+1 ];//sym==1

		if ( (l==1) && (myid==0) ) {
   		if (wf.symmetry_x4==1)
   		  {
   		    tridag_mid[wf.in4(l,k,j,i)] = 1.-arg_A-arg_B+arg_V;
   		    wf_rightside[wf.in4(l,k,j,i)]= -tridag_low_Fast*d.wf_1D[ l-1 ] +( 1.+arg_A+arg_B-arg_V )*d.wf_1D[ l ] -tridag_upp_Fast*d.wf_1D[ l+1 ];
   		  }
      	
   		if (wf.symmetry_x4==-1)
   		  {
   		    tridag_mid[wf.in4(l,k,j,i)] = 1.-3.*arg_A+arg_B+arg_V;
   		    wf_rightside[wf.in4(l,k,j,i)] = -tridag_low_Fast*d.wf_1D[ l-1 ] +( 1.+3.*arg_A-arg_B-arg_V )*d.wf_1D[ l ] -tridag_upp_Fast*d.wf_1D[ l+1 ];
   		  }
                 }
   	      }}}
          
     //Tridag_X4_MPI( wf, tridag_low_Fast , tridag_mid, tridag_upp_Fast , wf_rightside, nthreads);          
     Tridag_X4_MPI_opti( wf, tridag_low_Fast , tridag_mid, tridag_upp_Fast , wf_rightside, nthreads); // the number for every dimension should be even.
          
   } // end of X4_Laser_MPI



  vector<complex> Inverse_X1(const complex time_step , wavefunction &wf , const ABVparam &p)
  {
    const complex idt=I*time_step;
    const complex arg_A = ( idt*.5*wf.one_by_dx1sqr )*p.ABV_A_x1;
    const complex arg_B = 0.;
    complex arg_V;
    complex tridag_low_Fast =arg_A-arg_B;
    complex tridag_upp_Fast =arg_A+arg_B;
    vector<complex> tridag_mid(wf.n1+2,0.);
    vector<complex> tau(wf.n1+2,0.);
    vector<complex> alpha(wf.n1+2,0.);
    vector<complex> gama(wf.n1+2,0.);
    vector<complex> result((wf.n4+2)*(wf.n3+2)*(wf.n2+2)*(wf.n1+2)*(wf.n1+2),0.);
    int index;
    
    for(int l=1; l<=wf.n4; ++l)
      for(int k=1; k<=wf.n3; ++k)
	for(int j=1; j<=wf.n2; ++j)
	  {
	    for(int i=1; i<=wf.n1; ++i)
	      {
		arg_V=(idt*0.5)*p.ABV_V[wf.in4(l,k,j,i)];
		tridag_mid[i]=1.-2.*arg_A+arg_V;
		if (i==1) {
		  if (wf.symmetry_x1==1)
		    {
		      tridag_mid[ i ] =1.-arg_A-arg_B+arg_V;
		    }
		  if (wf.symmetry_x1==-1)
		    {
		      tridag_mid[ i ] =1.-3.*arg_A+arg_B+arg_V;
		    }
		  alpha[i]=tridag_mid[i];
		}	
		
	      }
	    for(int i=2; i<=wf.n1; ++i)
	      {
		tau[i-1]=tridag_upp_Fast/alpha[i-1];
		alpha[i]=tridag_mid[i]-tridag_low_Fast*tau[i-1];
		gama[i]=tridag_low_Fast/alpha[i-1];
	      }
	    index=l*(wf.n3+2)*(wf.n2+2)*(wf.n1+2)*(wf.n1+2)+k*(wf.n2+2)*(wf.n1+2)*(wf.n1+2)+j*(wf.n1+2)*(wf.n1+2);
	    result[index+wf.n1*(wf.n1+2)+wf.n1]=1./alpha[wf.n1];
	    for(int i=wf.n1-1; i>=1; --i)
	      {
		result[index+i*(wf.n1+2)+i]=1./alpha[i]+tau[i]*gama[i+1]*result[index+(i+1)*(wf.n1+2)+i+1];
	      }
	    for(int ir=wf.n1; ir>=2; --ir)
	      for(int ic=ir-1; ic>=1; --ic)
		{
		  result[index+ir*(wf.n1+2)+ic]=-gama[ic+1]*result[index+ir*(wf.n1+2)+ic+1];
		  result[index+ic*(wf.n1+2)+ir]=-tau[ic]*result[index+(ic+1)*(wf.n1+2)+ir];
		}
	  }
    return result;
  }


  vector<complex> Hamiltonian::Inverse_X1_tau(const complex time_step , wavefunction &wf , const ABVparam &p)
  {
    const complex idt=I*time_step;
    const complex arg_A = ( idt*.5*wf.one_by_dx1sqr )*p.ABV_A_x1;
    const complex arg_B = 0.;
    complex arg_V;
    complex tridag_low_Fast =arg_A-arg_B;
    complex tridag_upp_Fast =arg_A+arg_B;
    vector<complex> tridag_mid(wf.n1+2,0.);
    vector<complex> tau(wf.n1+2,0.);
    vector<complex> alpha(wf.n1+2,0.);
    //vector<complex> gama(wf.n1+2,0.);
    vector<complex> result((wf.n4+2)*(wf.n3+2)*(wf.n2+2)*(wf.n1+2),0.);
    int index;
    
    for(int l=1; l<=wf.n4; ++l)
      for(int k=1; k<=wf.n3; ++k)
	for(int j=1; j<=wf.n2; ++j)
	  {
	    for(int i=1; i<=wf.n1; ++i)
	      {
		arg_V=(idt*0.5)*p.ABV_V[wf.in4(l,k,j,i)];
		tridag_mid[i]=1.-2.*arg_A+arg_V;
		if (i==1) {
		  if (wf.symmetry_x1==1)
		    {
		      tridag_mid[ i ] =1.-arg_A-arg_B+arg_V;
		    }
		  if (wf.symmetry_x1==-1)
		    {
		      tridag_mid[ i ] =1.-3.*arg_A+arg_B+arg_V;
		    }
		  alpha[i]=tridag_mid[i];
		}	
		
	      }
	    for(int i=2; i<=wf.n1; ++i)
	      {
		tau[i-1]=tridag_upp_Fast/alpha[i-1];
		alpha[i]=tridag_mid[i]-tridag_low_Fast*tau[i-1];
		//gama[i]=tridag_low_Fast/alpha[i-1];
	      
		result[wf.in4(l,k,j,i-1)]=tau[i-1];
	      }
	  }
    return result;
  }


  vector<complex> Hamiltonian::Inverse_X1_alpha(const complex time_step , wavefunction &wf , const ABVparam &p)
  {
    const complex idt=I*time_step;
    const complex arg_A = ( idt*.5*wf.one_by_dx1sqr )*p.ABV_A_x1;
    const complex arg_B = 0.;
    complex arg_V;
    complex tridag_low_Fast =arg_A-arg_B;
    complex tridag_upp_Fast =arg_A+arg_B;
    vector<complex> tridag_mid(wf.n1+2,0.);
    vector<complex> tau(wf.n1+2,0.);
    vector<complex> alpha(wf.n1+2,0.);
    //vector<complex> gama(wf.n1+2,0.);
    vector<complex> result((wf.n4+2)*(wf.n3+2)*(wf.n2+2)*(wf.n1+2),0.);
    int index;
    
    for(int l=1; l<=wf.n4; ++l)
      for(int k=1; k<=wf.n3; ++k)
	for(int j=1; j<=wf.n2; ++j)
	  {
	    for(int i=1; i<=wf.n1; ++i)
	      {
		arg_V=(idt*0.5)*p.ABV_V[wf.in4(l,k,j,i)];
		tridag_mid[i]=1.-2.*arg_A+arg_V;
		if (i==1) {
		  if (wf.symmetry_x1==1)
		    {
		      tridag_mid[ i ] =1.-arg_A-arg_B+arg_V;
		    }
		  if (wf.symmetry_x1==-1)
		    {
		      tridag_mid[ i ] =1.-3.*arg_A+arg_B+arg_V;
		    }
		  alpha[i]=tridag_mid[i];
		}	
		
	      }
	    for(int i=2; i<=wf.n1; ++i)
	      {
		tau[i-1]=tridag_upp_Fast/alpha[i-1];
		alpha[i]=tridag_mid[i]-tridag_low_Fast*tau[i-1];
		//gama[i]=tridag_low_Fast/alpha[i-1];
	      
		result[wf.in4(l,k,j,i)]=alpha[i];
	      }
	  }
    return result;
  }


  vector<complex> Hamiltonian::Inverse_X1_gama(const complex time_step , wavefunction &wf , const ABVparam &p)
  {
    const complex idt=I*time_step;
    const complex arg_A = ( idt*.5*wf.one_by_dx1sqr )*p.ABV_A_x1;
    const complex arg_B = 0.;
    complex arg_V;
    complex tridag_low_Fast =arg_A-arg_B;
    complex tridag_upp_Fast =arg_A+arg_B;
    vector<complex> tridag_mid(wf.n1+2,0.);
    vector<complex> tau(wf.n1+2,0.);
    vector<complex> alpha(wf.n1+2,0.);
    vector<complex> gama(wf.n1+2,0.);
    vector<complex> result((wf.n4+2)*(wf.n3+2)*(wf.n2+2)*(wf.n1+2),0.);
    int index;
    
    for(int l=1; l<=wf.n4; ++l)
      for(int k=1; k<=wf.n3; ++k)
	for(int j=1; j<=wf.n2; ++j)
	  {
	    for(int i=1; i<=wf.n1; ++i)
	      {
		arg_V=(idt*0.5)*p.ABV_V[wf.in4(l,k,j,i)];
		tridag_mid[i]=1.-2.*arg_A+arg_V;
		if (i==1) {
		  if (wf.symmetry_x1==1)
		    {
		      tridag_mid[ i ] =1.-arg_A-arg_B+arg_V;
		    }
		  if (wf.symmetry_x1==-1)
		    {
		      tridag_mid[ i ] =1.-3.*arg_A+arg_B+arg_V;
		    }
		  alpha[i]=tridag_mid[i];
		}	
		
	      }
	    for(int i=2; i<=wf.n1; ++i)
	      {
		tau[i-1]=tridag_upp_Fast/alpha[i-1];
		alpha[i]=tridag_mid[i]-tridag_low_Fast*tau[i-1];
		gama[i]=tridag_low_Fast/alpha[i-1];
	      
		result[wf.in4(l,k,j,i)]=gama[i];
	      }
	  }
    return result;
  }





  void X1_kernel_inve(const complex time_step , wavefunction &wf , const ABVparam &p, temp_data_t d, int l, int k, int j, vector<complex> &tau, vector<complex> &alpha, vector<complex> &gama)
  {
    const complex idt=I*time_step;
    const complex arg_A = ( idt*.5*wf.one_by_dx1sqr )*p.ABV_A_x1;
    const complex arg_B = 0.;//complex(0.,0.);//( idt*.25*wf.one_by_dx1   )*p.ABV_B_x1;
    complex arg_V;
    complex tridag_low_Fast =arg_A-arg_B;
    complex tridag_upp_Fast =arg_A+arg_B;
    vector<complex> factor_inverse( (wf.n1+2)*(wf.n1+2),0.);

    
    for( int i = 0 ; i < wf.n1+2 ; ++i )//Copy the wavefunction
      {
	/**
	 * Creates a 1D potential and a 1D wavefunction out of the 3D ones
	 */
        const int index=wf.in4(l, k , j , i );
	d.v_1D[ i ] = p.ABV_V[ index ];
	d.wf_1D[ i ] = wf.wave[ index ];
      }
    
    
    for( int i = 1 ; i <= wf.n1 ; ++i )
      {
	/**
	 * Definition of the arguments in Eq. (7-9) in doc/ABV/abv.tex
	 */	
	arg_V = ( idt*.5 )*d.v_1D[ i ];
	
	/**
	 * Define the 3 diagonals X of the left side of Eq. (6)  in doc/ABV/abv.tex => look at Eqs. (7-9) in doc/ABV/abv.tex
	 */
	d.tridag_mid[ i ] =1.-2.*arg_A+arg_V;	
	d.wf_1D_rightside[ i ]=-tridag_low_Fast*d.wf_1D[ i-1 ]+( 1.+2.*arg_A-arg_V )*d.wf_1D[ i ] -tridag_upp_Fast*d.wf_1D[ i+1 ];//general case, symmetry==1
	
	/**
	 * for (anti-)symmetric wavefunctions one has to take care of the right boundary conditions (phi(0)=+/-phi(1). That leads to modified values.
	 */
	if (i==1) {
	  if (wf.symmetry_x1==1)
	    {
	      d.tridag_mid[ i ] =1.-arg_A-arg_B+arg_V;
	      d.wf_1D_rightside[ i ]=-tridag_low_Fast*d.wf_1D[ i-1 ] +( 1.+arg_A+arg_B-arg_V )*d.wf_1D[ i ]-tridag_upp_Fast*d.wf_1D[ i+1 ];
	    }
	  
	  if (wf.symmetry_x1==-1)
	    {
	      d.tridag_mid[ i ] =1.-3.*arg_A+arg_B+arg_V;
	      d.wf_1D_rightside[ i ]=-tridag_low_Fast*d.wf_1D[ i-1 ] +( 1.+3.*arg_A-arg_B-arg_V )*d.wf_1D[ i ]-tridag_upp_Fast*d.wf_1D[ i+1 ];
	    }
	}
	
      }


    //------------ assign the value of inverse matrix---------//
    factor_inverse[wf.in2(wf.n1,wf.n1)]=1./alpha[wf.n1];
    for(int i=wf.n1-1; i>=1; --i)
      {
	factor_inverse[i*(wf.n1+2)+i]=1./alpha[i]+tau[i]*gama[i+1]*factor_inverse[(i+1)*(wf.n1+2)+i+1];
      }
			
    for(int ir=wf.n1; ir>=2; --ir)
      for(int ic=ir-1; ic>=1; --ic)
	{
	  factor_inverse[ir*(wf.n1+2)+ic]=-gama[ic+1]*factor_inverse[ir*(wf.n1+2)+ic+1];
	  factor_inverse[ic*(wf.n1+2)+ir]=-tau[ic]*factor_inverse[(ic+1)*(wf.n1+2)+ir];
	}

    //-------- new wf calculation---------------//
    
    
    for( int i = 1 ; i <= wf.n1 ; i++ )
      {
	for (int ic = 1; ic <= wf.n1; ++ic)
	  {
	    d.wf_1D_solution[ i ] += factor_inverse[i*(wf.n1+2)+ic]*d.wf_1D_rightside[ ic ];
	  }
	wf.wave[ wf.in4(l, k , j , i ) ] = d.wf_1D_solution[ i ];
      }
    
    
  }


  void Hamiltonian::X1_inve( const complex time_step , wavefunction &wf , const ABVparam &p, vector<complex> &tau_load, vector<complex> &alpha_load, vector<complex> &gama_load, const int nthreads )
  {
    temp_data_t &d=temp_data[1]; 
    if(n[1] != wf.n1 || dx[1] != wf.dx1) {
      cerr<<"wavefunction has wrong size for Hamiltonian\n";
      exit(1);
    }

    vector<complex> tau(wf.n1+2,0.);
    vector<complex> alpha(wf.n1+2,0.);
    vector<complex> gama(wf.n1+2,0.);
    
    // int nthreads=16;
    // double wall_timer=omp_get_wtime();
    // clock_t clock_timer=clock();

    #pragma omp parallel for  num_threads(nthreads) collapse(3)   
    for(int l = 1 ; l <= wf.n4 ; l++ )
      for(int k = 1 ; k <= wf.n3 ; k++ )
	{   
	  for(int j = 1 ; j <= wf.n2 ; j++ )
	    {
	      for(int i=1; i<=wf.n1; i++)
		{
		  tau[i]=tau_load[wf.in4(l,k,j,i)];
		  alpha[i]=alpha_load[wf.in4(l,k,j,i)];
		  gama[i]=gama_load[wf.in4(l,k,j,i)];
		}	   
	      X1_kernel_inve(time_step, wf , p, d, l, k, j, tau, alpha, gama);
	    }
	}

  } // end of X1_3D 







/*------------------  Wavefunction related functions  -----------------------*/



   /**
   * Initialize wavefunction
   */
  void Initialize( wavefunction &wf , vector<int> n_points , vector<double> spatial_steps , int symmetry_x1, int symmetry_x2,int symmetry_x3,int symmetry_x4 )
  {
    int numprocs, nlocal, myid;
    MPI_Comm_rank (MPI_COMM_WORLD, &myid);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    nlocal=n_points[ 3 ]/numprocs;

    wf.verbose=false;
    /**
     * Initialize symmetries of the wavefunction.  
     */
    wf.symmetry_x1=symmetry_x1;
    wf.symmetry_x2=symmetry_x2;
    wf.symmetry_x3=symmetry_x3;
    wf.symmetry_x4=symmetry_x4;
   
    
    /**
     * Initialize the number of grid points.
     */
    
   
    wf.n1 = n_points[ 0 ];
    wf.n2 = n_points[ 1 ];
    wf.n3 = n_points[ 2 ];
    wf.n4 = nlocal;

   
    
    /**
     * Allocate the wavefunction and the potential
     */
    wf.wave.resize(( wf.n1+2 )*( wf.n2+2 )*( wf.n3+2 )*(wf.n4+2), 0.);
    for (int l=1;l<=wf.n4;l++)
      for (int k=1;k<=wf.n3;k++)
	for (int j=1;j<=wf.n2;j++)
	  for (int i=1;i<=wf.n1;i++)
	        {
		  wf.wave[ wf.in4( l, k , j , i ) ] =0.;
	        }


    /**
     * Allocate the grid
     */
    if (abs(symmetry_x1) == 1) {
      wf.x1.type=HalfAxis;
      wf.p1.type=FFTWHalfAxis;
    }
    else {
      wf.x1.type=FullAxis;
      wf.p1.type=FFTWAxis;
    }
    if (abs(symmetry_x2) == 1) {
      wf.x2.type=HalfAxis;
      wf.p2.type=FFTWHalfAxis;
    }
    else {
      wf.x2.type=FullAxis;
      wf.p2.type=FFTWAxis;
    }
    if (abs(symmetry_x3) == 1) {
      wf.x3.type=HalfAxis;
      wf.p3.type=FFTWHalfAxis;
    }
    else {
      wf.x3.type=FullAxis;
      wf.p3.type=FFTWAxis;
    }

    if (abs(symmetry_x4) == 1) {
      wf.x4.type=HalfAxis;
      wf.p4.type=FFTWHalfAxis;
    }
    else {
      wf.x4.type=FullAxis;
      wf.p4.type=FFTWAxis;
    }

    wf.x1.Init("x1", n_points[0], spatial_steps[0], wf.x1.type);
    wf.x2.Init("y1", n_points[1], spatial_steps[1], wf.x2.type);
    wf.x3.Init("x2", n_points[2], spatial_steps[2], wf.x3.type);

    wf.x4.Init_MPI("y2", n_points[3], nlocal, myid, spatial_steps[3], wf.x4.type);

    wf.dp1 = 2*pi/(wf.dx1*wf.n1);
    wf.dp2 = 2*pi/(wf.dx2*wf.n2);
    wf.dp3 = 2*pi/(wf.dx3*wf.n3);
    wf.dp4 = 2*pi/(wf.dx4*wf.n4*numprocs);

    wf.p1.Init("p_x1", n_points[0], wf.dp1, wf.p1.type);
    wf.p2.Init("p_y1", n_points[1], wf.dp2, wf.p2.type);
    wf.p3.Init("p_x2", n_points[2], wf.dp3, wf.p3.type);

    wf.p4.Init_MPI("p_y2",n_points[3], nlocal, myid, wf.dp4, wf.p4.type);
    
    /**
     * Initialize the spatial grid helper constants.
     */
    wf.one_by_dx1=1. / wf.dx1;     // wf.dx1=x[1].delta (which itself is set in Axis.Init)
    wf.one_by_dx2=1. / wf.dx2;
    wf.one_by_dx3=1. / wf.dx3;
    wf.one_by_dx4=1. / wf.dx4;
    wf.one_by_dx1sqr=1. / ( wf.dx1*wf.dx1 );
    wf.one_by_dx2sqr=1. / ( wf.dx2*wf.dx2 );
    wf.one_by_dx3sqr=1. / ( wf.dx3*wf.dx3 );
    wf.one_by_dx4sqr=1. / ( wf.dx4*wf.dx4 );
 
    //}
  } // Initialize_wavefunction



  
   void Initialize_Momentum( wavefunction &wf , wavefunction &wf_mom)
  {
    
    int myid,numprocs;
    MPI_Comm_rank (MPI_COMM_WORLD, &myid);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);

    wf_mom.verbose=false;
    
    /**
     * Initialize symmetries of the wavefunction.  
     */
    wf_mom.symmetry_x1=wf.symmetry_x1;
    wf_mom.symmetry_x2=wf.symmetry_x2;
    wf_mom.symmetry_x3=wf.symmetry_x3;
    wf_mom.symmetry_x4=wf.symmetry_x4;
    //symmetry=0: X starts from -X_max, no symmetry or antisymmetry for wavefunction necessary
    //symmetry=1: X starts from 0, wavefunction is symmetric
    //symmetry=-1: X starts from 0, wavefunction is antisymmetric
    wf_mom.n1 = wf.n1;
    wf_mom.n2 = wf.n2;
    wf_mom.n3 = wf.n3;
    wf_mom.n4 = wf.n4;

    
    /**
     * Allocate the wavefunction and the potential
     */
    wf_mom.wave.resize(( wf_mom.n1+2 )*( wf_mom.n2+2 )*( wf_mom.n3+2 )*(wf_mom.n4+2), 0.);
    for (int l=1;l<=wf_mom.n4;l++)
      for (int k=1;k<=wf_mom.n3;k++)
	for (int j=1;j<=wf_mom.n2;j++)
	  for (int i=1;i<=wf_mom.n1;i++)
	  {
	    wf_mom.wave[ wf_mom.in4(l, k , j , i ) ] =0.;
	  }

    
    /**
     * Initialize the spatial grid steps and helper constants
     */
    wf_mom.dx1 = 2*pi/(wf.dx1*wf_mom.n1);
    wf_mom.dx2 = 2*pi/(wf.dx2*wf_mom.n2);
    wf_mom.dx3 = 2*pi/(wf.dx3*wf_mom.n3);
    wf_mom.dx4 = 2*pi/(wf.dx4*wf_mom.n4*numprocs);
    wf_mom.one_by_dx1=1. / wf_mom.dx1;     // constants
    wf_mom.one_by_dx2=1. / wf_mom.dx2;
    wf_mom.one_by_dx3=1. / wf_mom.dx3;
    wf_mom.one_by_dx4=1. / wf_mom.dx4;
    wf_mom.one_by_dx1sqr=1. / ( wf_mom.dx1*wf_mom.dx1 );
    wf_mom.one_by_dx2sqr=1. / ( wf_mom.dx2*wf_mom.dx2 );
    wf_mom.one_by_dx3sqr=1. / ( wf_mom.dx3*wf_mom.dx3 );
    wf_mom.one_by_dx4sqr=1. / ( wf_mom.dx4*wf_mom.dx4 );
    
    /**
     * Allocate the grid
     */
    // !!!!!!!!!!! only symmetry =0 !!!!!!!!!!!!
    // if (abs(wf_mom.symmetry_x1) == 1) wf_mom.x1.type=FFTWHalfAxis;
    // else wf_mom.x1.type=FFTWAxis;
    // if (abs(wf_mom.symmetry_x2) == 1) wf_mom.x2.type=FFTWHalfAxis;
    // else wf_mom.x2.type=FFTWAxis;
    // if (abs(wf_mom.symmetry_x3) == 1) wf_mom.x3.type=FFTWHalfAxis;
    // else wf_mom.x3.type=FFTWAxis;
    // if (abs(wf_mom.symmetry_x4) == 1) wf_mom.x4.type=FFTWHalfAxis;
    // else wf_mom.x4.type=FFTWAxis;
    wf_mom.x1.type=FFTWAxis;
    wf_mom.x2.type=FFTWAxis;
    wf_mom.x3.type=FFTWAxis;
    wf_mom.x4.type=FFTWAxis;
    wf_mom.x1.Init("p_x1", wf_mom.n1, wf_mom.dx1,  wf_mom.x1.type);
    wf_mom.x2.Init("p_y1", wf_mom.n2, wf_mom.dx2,  wf_mom.x2.type);
    wf_mom.x3.Init("p_x2", wf_mom.n3, wf_mom.dx3,  wf_mom.x3.type);
    wf_mom.x4.Init_MPI("p_y2", wf_mom.n4*numprocs, wf_mom.n4, myid, wf_mom.dx4, wf_mom.x4.type);
    
  } // Initialize_momentum
 





    /**
     * Initialize a Gaussian function as a initial wavepacket for imaginary time propagation
     */ 
  void Guess_Function_Gauss( wavefunction &wf , vector<double> sigmas , vector<double> gauss_centers , vector<double> initial_momenta )
  {
    
    for (int l=1;l<= wf.n4;l++)
      for( int k = 1 ; k <= wf.n3 ; k++ )
	{
	  for( int j = 1 ; j <= wf.n2 ; j++ )
	    {
	      for( int i = 1 ; i <= wf.n1 ; i++ )
		{
		double arg_x1 = ( wf.x1[ i ]-gauss_centers[ 0 ] )/( 2.*sigmas[ 0 ] );
		double arg_x2 = ( wf.x2[ j ]-gauss_centers[ 1 ] )/( 2.*sigmas[ 1 ] );
		double arg_x3 = ( wf.x3[ k ]-gauss_centers[ 2 ] )/( 2.*sigmas[ 2 ] );
		double arg_x4 = ( wf.x4[ l ]-gauss_centers[ 3 ] )/( 2.*sigmas[ 3 ] );
		complex arg_p1 = I*initial_momenta[ 0 ]*wf.x1[ i ];
		complex arg_p2 = I*initial_momenta[ 1 ]*wf.x2[ j ];
		complex arg_p3 = I*initial_momenta[ 2 ]*wf.x3[ k ];
		complex arg_p4 = I*initial_momenta[ 3 ]*wf.x4[ l ];
		wf.wave[ wf.in4(l, k , j , i ) ] = exp( -arg_x1*arg_x1-arg_x2*arg_x2-arg_x3*arg_x3 - arg_x4*arg_x4 )*exp( arg_p1+arg_p2+arg_p3+arg_p4 );
		}
	    }
	}   
  } // end of Guess_Function_Gauss

  
    /**
     * Normalize of the wavefunction phi' = phi/sqrt( <phi|phi> )
     */
   void Normalize( wavefunction &wf )
  {
    double norm = Obs_Norm_Squared( wf );
    if(wf.verbose==true)
      cout << "Normalizing: Input Norm = " << norm << "\n";

    double one_by_norm = 1./sqrt(norm);
   
    
    vector<complex>::const_iterator off=wf.wave.end();
   

    for( vector<complex>::iterator element_iterator=wf.wave.begin() + 1 ;
	 element_iterator != off ;
	 ++element_iterator ) 
      {
	(*element_iterator) *= one_by_norm; 
      }
    
    norm = Obs_Norm_Squared( wf );
    
    if(wf.verbose==true)
      cout << "Normalizing: Output Norm = " << norm << "\n";
    
  } // end of Normalize

  

  void PlaceWaveFunction( wavefunction &wf , wavefunction &small_wf)
  {
    
    for (int i=1; i<=3; i++) {
      if(small_wf.n[i] > wf.n[i]) {
	cerr<<"cannot load bigger wavefunction into small grid"<<endl;
	exit(1);
      }
    }
    if(small_wf.n[4] != wf.n[4])
      {
	cerr<<"cannot load wf which has diff points in X4"<<endl;
	exit(1);    
      }
    
    int shifter_n1=0;
    int shifter_n2=0;
    int shifter_n3=0;
    //int shifter_n4=0;
    
    if (wf.symmetry_x1==0)
      {
	shifter_n1=(wf.n1-small_wf.n1)/2;
      }	  
    if (wf.symmetry_x2==0)
      {
	shifter_n2=(wf.n2-small_wf.n2)/2;
      }
    if (wf.symmetry_x3==0)
      {
	shifter_n3=(wf.n3-small_wf.n3)/2;
      }
    // if (wf.symmetry_x4==0)
    //	{
    //  shifter_n4=(wf.n4-small_wf.n4)/2;
    //	}
    
    for(int l=0;l<wf.n4+2;l++)
      for(int k=0;k<wf.n3+2;k++)
	for(int j=0;j<wf.n2+2;j++)
	  for(int i=0;i<wf.n1+2;i++)
	    {
	      wf.wave[wf.in4(l,k,j,i)]=0.;
	    }
    
    for(int l=1;l<=small_wf.n4;l++)
      for(int k=1;k<=small_wf.n3;k++)
	for(int j=1;j<=small_wf.n2;j++)
	  for(int i=1;i<=small_wf.n1;i++)
	    {
	      wf.wave[wf.in4(l,shifter_n3+k,shifter_n2+j,shifter_n1+i)]=small_wf.wave[small_wf.in4(l,k,j,i)];
	    }
    
  }
  
	
  complex Project( wavefunction &wf , wavefunction &wf2)
  {
    complex result=0.;
    complex projection=complex(0.,0.);
    for(int l=1;l<=wf.n4;l++)
      for(int k=1;k<=wf.n3;k++)
	for(int j=1;j<=wf.n2;j++)
	  for(int i=1;i<=wf.n1;i++)
	    {
	      projection+= conj(wf.wave[wf.in4(l,k,j,i)])*wf2.wave[wf2.in4(l,k,j,i)];
	    }
    projection*=wf.dx1*wf.dx2*wf.dx3*wf.dx4;
    MPI_Allreduce(&projection, &result, 1, MPI_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
    return result;
  }
	
	
	
	
  double Overlap( wavefunction &wf , wavefunction &wf2)
  {
    double overlap_;
    complex proj = Project ( wf , wf2 );
    overlap_ = ( real( proj*conj( proj ) ) );
    return overlap_;
  }
  

  void ProjectOUT_Diff_Sizes( wavefunction &wf , wavefunction &small_wf)
  {
    for (int i=1; i<=3; i++) {
      if(small_wf.n[i] > wf.n[i]) {
	cerr<<"cannot project bigger wavefunction from small grid"<<endl;
	exit(1);
      }
    }
    if(small_wf.n[4] != wf.n[4])
      {
	cerr<<"cannot project wf which has diff points in X4"<<endl;
	exit(1);    
      }



    int shifter_n1=0;
    int shifter_n2=0;
    int shifter_n3=0;//if sym==0 or 2, shift_n==0;
    
    if (wf.symmetry_x1==0)
      {
	shifter_n1=(wf.n1-small_wf.n1)/2;
      }
    if (wf.symmetry_x2==0)
      {
	shifter_n2=(wf.n2-small_wf.n2)/2;
      }
    if (wf.symmetry_x3==0)
      {
	shifter_n3=(wf.n3-small_wf.n3)/2;
      }
    
    int nthreads=16;
    complex projection=complex(0.,0.);	
    
    for(int l=1;l<=small_wf.n4;l++)
      for(int k=1;k<=small_wf.n3;k++)
	for(int j=1;j<=small_wf.n2;j++)
	  for(int i=1;i<=small_wf.n1;i++)
	    {
	      projection+= conj(small_wf.wave[small_wf.in4(l,k,j,i)])
		*wf.wave[wf.in4(l,shifter_n3+k,shifter_n2+j,shifter_n1+i)];
	    }
    
    projection*=small_wf.dx1*small_wf.dx2*small_wf.dx3*small_wf.dx4;
    
    for(int l=1;l<=small_wf.n4;l++)
      for(int k=1;k<=small_wf.n3;k++)
	for(int j=1;j<=small_wf.n2;j++)
	  for(int i=1;i<=small_wf.n1;i++)
	    {
	      wf.wave[wf.in4(l,shifter_n3+k,shifter_n2+j,shifter_n1+i)]-=projection*small_wf.wave[small_wf.in4(l,k,j,i)];
	    }
    
    //	  return projection;
  }



	
  void dig_h2p_e3d_r1d(wavefunction &wf, double radius,wavefunction &wf_new)
  {	/// to be rewrite
    double e1r,e2r;
    for (int l=1; l<=wf.n4; l++)
      for (int k=1; k<=wf.n3; k++)
	for (int j=1; j<=wf.n2; j++)
	  for (int i=1; i<=wf.n1; i++)
	    {
	      e1r= ( wf.x2[ j ]+ 0.5*wf.x1[i])*( wf.x2[ j ]+ 0.5*wf.x1[i]) + wf.x3[k]*wf.x3[k]+wf.x4[l]*wf.x4[l] ;
	      e2r= ( wf.x2[ j ]- 0.5*wf.x1[i])*( wf.x2[ j ]- 0.5*wf.x1[i]) + wf.x3[k]*wf.x3[k]+wf.x4[l]*wf.x4[l] ;
	      if (sqrt(e1r)<radius || sqrt(e2r)<radius  )	
		wf_new.wave[wf_new.in4(l,k,j,i)]= 0.;
	      else
		wf_new.wave[wf_new.in4(l,k,j,i)]= wf.wave[wf.in4(l,k,j,i)];
	    }
  }



  void dig_h2p_e2d_r2d(wavefunction &wf, double radius,wavefunction &wf_new)
  {	/// to be rewrite
    double e1r,e2r;
    for (int l=1; l<=wf.n4; l++)
      for (int k=1; k<=wf.n3; k++)
	for (int j=1; j<=wf.n2; j++)
	  for (int i=1; i<=wf.n1; i++)
	    {
	      e1r= sqrt( (wf.x1[i]/2.+wf.x3[k])*(wf.x1[i]/2.+wf.x3[k])+(wf.x2[j]/2.+wf.x4[l])*(wf.x2[j]/2.+wf.x4[l]) ) ;
	      e2r= sqrt( (wf.x1[i]/2.-wf.x3[k])*(wf.x1[i]/2.-wf.x3[k])+(wf.x2[j]/2.-wf.x4[l])*(wf.x2[j]/2.-wf.x4[l]) ) ;
	      if ( e1r<radius || e2r<radius  )	
		wf_new.wave[wf_new.in4(l,k,j,i)]= 0.;
	      else
		wf_new.wave[wf_new.in4(l,k,j,i)]= wf.wave[wf.in4(l,k,j,i)];
	    }
  }

    void dig_h2p_e2d_r1d_ionz(wavefunction &wf, double radius,wavefunction &wf_new)
  {	/// to be rewrite
    double e1r,e2r;
    for (int l=1; l<=wf.n4; l++)
      for (int k=1; k<=wf.n3; k++)
	for (int j=1; j<=wf.n2; j++)
	  for (int i=1; i<=wf.n1; i++)
	    {
	      e1r= sqrt( (wf.x2[j]/2.+wf.x3[k])*(wf.x2[j]/2.+wf.x3[k])+wf.x4[l]*wf.x4[l] ) ;
	      e2r= sqrt( (wf.x2[j]/2.-wf.x3[k])*(wf.x2[j]/2.-wf.x3[k])+wf.x4[l]*wf.x4[l] ) ;
	      if ( e1r<radius || e2r<radius  )	
		wf_new.wave[wf_new.in4(l,k,j,i)]= wf.wave[wf.in4(l,k,j,i)]*exp(-10.*abs(radius-e1r)*abs(radius-e2r));
	      else
		wf_new.wave[wf_new.in4(l,k,j,i)]= wf.wave[wf.in4(l,k,j,i)];
	    }
  }





  void Rotate(wavefunction &tmp1, wavefunction &tmp2, wavefunction &xg, wavefunction &xu, vector<double> cp, const double dt, const double field )
  {
    for ( int l=1; l<=xg.n4; l++)
    for ( int j=1; j<=xg.n2; j++)
      for ( int k=1; k<=xg.n3; k++)
        for ( int i=1; i<=xg.n1; i++)
          {
            tmp1.wave[xg.in4(l,k,j,i )]=xg.wave[xg.in4(l,k,j,i )]*cos(cp[i]*field*dt)-I*sin(cp[i]*field*dt)*xu.wave[xg.in4(l,k,j,i )];
            tmp2.wave[xg.in4(l,k,j,i )]=xu.wave[xu.in4(l,k,j,i )]*cos(cp[i]*field*dt)-I*sin(cp[i]*field*dt)*xg.wave[xg.in4(l,k,j,i )];
          }
    for ( int l=1; l<=xg.n4; l++)
    for ( int j=1; j<=xg.n2; j++)
      for ( int k=1; k<=xg.n3; k++)
        for ( int i=1; i<=xg.n1; i++)
          {
            xg.wave[xg.in4(l,k,j,i )] = tmp1.wave[xg.in4(l,k,j,i )];
            xu.wave[xu.in4(l,k,j,i )] = tmp2.wave[xg.in4(l,k,j,i )];
          }
  }// the rotating in the two state





/*-------------------Calculation of Observables:--------------*/

  
  
  /**
   * Calculate the norm ( = sqrt( <phi|phi> )
   */
  double Obs_Norm( wavefunction &wf )
  {
    double obs=0.;
    double result=0.;
    
    for( int l = 1 ; l <= wf.n4 ; l++ )
      for( int k = 1 ; k <= wf.n3 ; k++ )
	for( int j = 1 ; j <= wf.n2 ; j++ )
	  for( int i = 1 ; i <= wf.n1 ; i++ )
	    {
	      obs += norm( wf.wave[ wf.in4( l, k , j , i ) ] );
	    }

    MPI_Allreduce(&obs, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    result *= wf.dx1*wf.dx2*wf.dx3*wf.dx4;
    result = sqrt( result );
    return result;

  } // end of Obs_Norm

 

   double Obs_Norm_Squared( wavefunction &wf )
  {

    double obs=0.;
    double result;
    for( int l = 1 ; l <= wf.n4 ; l++ )
      {
	for( int k = 1 ; k <= wf.n3 ; k++ )
	  {
	    for( int j = 1 ; j <= wf.n2 ; j++ )
	      {
		for( int i = 1 ; i <= wf.n1 ; i++ )
		  {
		    obs += norm( wf.wave[ wf.in4(l, k , j , i ) ] );
		  }
	      }
	  }
      }
    obs *= wf.dx1*wf.dx2*wf.dx3*wf.dx4;
 
    MPI_Allreduce(&obs, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   
    return result;

  } // end of Obs_Norm_Squared

 
 double Obs_Energy_1D( wavefunction &wf, ABVparam &p)
 {
   double obs=0.;
   double obs_all=0.;

  
	// //////////////**************** calculate the observable
   for( int l = 1 ; l <= wf.n4 ; l++ )
     for( int k = 1 ; k <= wf.n3 ; k++ )
       {
	 for( int j = 1 ; j <= wf.n2 ; j++ )
	   {
	     for( int i = 1 ; i <= wf.n1 ; i++ )
	       {
		 obs += real( conj( wf.wave[ wf.in4(l, k , j , i ) ] )*(
									wf.wave[ wf.in4(l, k , j , i ) ]*( p.ABV_V[ wf.in4(l, k , j , i ) ] )
									// +Second_Derivative_X1( wf , k , j , i )*p.ABV_A_x1
									//+Second_Derivative_X2( wf , k , j , i )*p.ABV_A_x2
									+Second_Derivative_X4_MPI( wf ,l, k , j , i )*p.ABV_A_x4
									) );
	       }
	   }
       }
   
   
   obs *= wf.dx1*wf.dx2*wf.dx3*wf.dx4;
   
 
   MPI_Allreduce(&obs, &obs_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   
 
   
   return obs_all;
 
 } // end of Obs_Energy_1D




  
   double Obs_Energy_2D( wavefunction &wf, ABVparam &p)
  {
    double obs=0.;
    double obs_all=0.;
    for( int l = 1 ; l <= wf.n4 ; l++ )
    for( int k = 1 ; k <= wf.n3 ; k++ )
      {
  	for( int j = 1 ; j <= wf.n2 ; j++ )
  	  {
  	    for( int i = 1 ; i <= wf.n1 ; i++ )
  	      {
  		obs += real( conj( wf.wave[ wf.in4( l, k , j , i ) ] )*(
									wf.wave[ wf.in4( l, k , j , i ) ]*( p.ABV_V[ wf.in4( l, k , j , i ) ]*2.  )
									+Second_Derivative_X3( wf , l, k , j , i )*p.ABV_A_x3
									+Second_Derivative_X4_MPI( wf , l, k , j , i )*p.ABV_A_x4
  								     ) );
  	      }
  	  }
      }
    obs *= wf.dx1*wf.dx2*wf.dx3*wf.dx4;
    MPI_Allreduce(&obs, &obs_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return obs_all;
   
  } // end of Obs_Energy_2D

   double Obs_Energy_3D( wavefunction &wf, ABVparam &p)
  {
    double obs=0.;
    double obs_all=0.;
    for( int l = 1 ; l <= wf.n4 ; l++ )
    for( int k = 1 ; k <= wf.n3 ; k++ )
      {
	for( int j = 1 ; j <= wf.n2 ; j++ )
	  {
	    for( int i = 1 ; i <= wf.n1 ; i++ )
	      {
		obs += real( conj( wf.wave[ wf.in4( l, k , j , i ) ] )*(
									wf.wave[ wf.in4(l, k , j , i ) ]*( p.ABV_V[ wf.in4(l, k , j , i ) ]*3.  )
									//+Second_Derivative_X1( wf ,l, k , j , i )*p.ABV_A_x1
									+Second_Derivative_X2( wf ,l, k , j , i )*p.ABV_A_x2
									+Second_Derivative_X3( wf ,l, k , j , i )*p.ABV_A_x3
									+Second_Derivative_X4_MPI( wf ,l, k , j , i )*p.ABV_A_x4
								     ) );

	
	      }
	  }
      }
    obs *= wf.dx1*wf.dx2*wf.dx3*wf.dx4;
    MPI_Allreduce(&obs, &obs_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return obs_all;
   

  } // end of Obs_Energy_3D


   double Obs_Energy_4D( wavefunction &wf, ABVparam &p)
  {
    double obs=0.;
    double obs_all=0.;
    for( int l = 1 ; l <= wf.n4 ; l++ )
    for( int k = 1 ; k <= wf.n3 ; k++ )
      {
	for( int j = 1 ; j <= wf.n2 ; j++ )
	  {
	    for( int i = 1 ; i <= wf.n1 ; i++ )
	      {
		obs += real( conj( wf.wave[ wf.in4( l, k , j , i ) ] )*(
									wf.wave[ wf.in4(l, k , j , i ) ]*( p.ABV_V[ wf.in4(l, k , j , i ) ]*4.  )
									+Second_Derivative_X1( wf ,l, k , j , i )*p.ABV_A_x1
									+Second_Derivative_X2( wf ,l, k , j , i )*p.ABV_A_x2
									+Second_Derivative_X3( wf ,l, k , j , i )*p.ABV_A_x3
									+Second_Derivative_X4_MPI( wf ,l, k , j , i )*p.ABV_A_x4
								     ) );

	
	      }
	  }
      }
    obs *= wf.dx1*wf.dx2*wf.dx3*wf.dx4;
    MPI_Allreduce(&obs, &obs_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return obs_all;
   

  } // end of Obs_Energy_4D


 


 
    /**
     * Calculate the expectation value <x1> = <phi|x1|phi>
     */
  double Obs_Expectation_Value_X1( wavefunction &wf )
  {
    double obs = 0.;
    double obs_all = 0.;
    for( int l = 1 ; l <= wf.n4 ; l++ )
    for( int k = 1 ; k <= wf.n3 ; k++ )
      {
	for( int j = 1 ; j <= wf.n2 ; j++ )
	  {
	    for( int i = 1 ; i <= wf.n1 ; i++ )
	      {
		obs += wf.x1[ i ] *norm( wf.wave[ wf.in4( l, k , j , i ) ] );
	      }
	  }
      }
    obs *= wf.dx1*wf.dx2*wf.dx3*wf.dx4;

    MPI_Allreduce(&obs, &obs_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return obs_all;
    
    
  } // end of Obs_Expectation_Value_X1



    /**
     * Calculate the expectation value <|x1|> = <phi||x1||phi>
     */
  double Obs_Expectation_Value_Abs_X1( wavefunction &wf )
  {
    double obs = 0.;
    double obs_all = 0.;
    for( int l = 1 ; l <= wf.n4 ; l++ )
    for( int k = 1 ; k <= wf.n3 ; k++ )
      {
	for( int j = 1 ; j <= wf.n2 ; j++ )
	  {
	    for( int i = 1 ; i <= wf.n1 ; i++ )
	      {
		obs += fabs(wf.x1[ i ]) *norm( wf.wave[ wf.in4( l, k , j , i ) ] );
	      }
	  }
      }
    obs *= wf.dx1*wf.dx2*wf.dx3*wf.dx4;

    MPI_Allreduce(&obs, &obs_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return obs_all; 

    
    
  } // end of Obs_Expectation_Value_AbsX1



    /**
     * Calculate the expectation value <x2> = <phi|x2|phi>
     */
  double Obs_Expectation_Value_X2( wavefunction &wf )
  {
    double obs = 0.;
    double obs_all = 0.;
    for( int l = 1 ; l <= wf.n4 ; l++ )
    for( int k = 1 ; k <= wf.n3 ; k++ )
      {
	for( int j = 1 ; j <= wf.n2 ; j++ )
	  {
	    for( int i = 1 ; i <= wf.n1 ; i++ )
	      {
		obs += wf.x2[ j ] *norm( wf.wave[ wf.in4( l, k , j , i ) ] );
	      }
	  }
      }
    obs *= wf.dx1*wf.dx2*wf.dx3*wf.dx4;

    MPI_Allreduce(&obs, &obs_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return obs_all;


    
  } // end of Obs_Expectation_Value_X2



    /**
     * Calculate the expectation value <x3> = <phi|x3|phi>
     */
  double Obs_Expectation_Value_X3( wavefunction &wf )
  {
    double obs = 0.;
    double obs_all = 0.;
    for( int l = 1 ; l <= wf.n4 ; l++ )
    for( int k = 1 ; k <= wf.n3 ; k++ )
      {
	for( int j = 1 ; j <= wf.n2 ; j++ )
	  {
	    for( int i = 1 ; i <= wf.n1 ; i++ )
	      {
		obs += wf.x3[ k ] *norm( wf.wave[ wf.in4( l, k , j , i ) ] );
	      }
	  }
      }
    obs *= wf.dx1*wf.dx2*wf.dx3*wf.dx4;

    MPI_Allreduce(&obs, &obs_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return obs_all;

    
  } // end of Obs_Expectation_Value_X3


  double Obs_Expectation_Value_X4( wavefunction &wf )
  {
    double obs = 0.;
    double obs_all = 0.;
    for( int l = 1 ; l <= wf.n4 ; l++ )
    for( int k = 1 ; k <= wf.n3 ; k++ )
      {
	for( int j = 1 ; j <= wf.n2 ; j++ )
	  {
	    for( int i = 1 ; i <= wf.n1 ; i++ )
	      {
		obs += wf.x4[ l ] *norm( wf.wave[ wf.in4( l, k , j , i ) ] );
	      }
	  }
      }
    obs *= wf.dx1*wf.dx2*wf.dx3*wf.dx4;

    MPI_Allreduce(&obs, &obs_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return obs_all;

    
  } // end of Obs_Expectation_Value_X4


  double Obs_Alignment( wavefunction &wf )
  {
    double obs=0.;
    double obs_all=0.;
    for( int l = 1 ; l <= wf.n4 ; l++ )
    for( int k = 1 ; k <= wf.n3 ; k++ )
      {
	for( int j = 1 ; j <= wf.n2 ; j++ )
	  {
	    for( int i = 1 ; i <= wf.n1 ; i++ )
	      {
		obs += wf.x1[i]*wf.x1[i]/(wf.x1[i]*wf.x1[i]+wf.x2[j]*wf.x2[j])*norm( wf.wave[ wf.in4( l, k , j , i ) ] );
		
	      }
	  }
      }
    obs *= wf.dx1*wf.dx2*wf.dx3*wf.dx4;
    obs = obs/Obs_Norm_Squared(wf);
    MPI_Allreduce(&obs, &obs_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return obs_all;

  } // end of Obs_Alignment between x1 and molecular axis






  vector<double> Obs_Fix_X1 (wavefunction &wf, double x1)
  {
    vector<double> result;
    result.resize( wf.n4*wf.n3*wf.n2, 0.);

    for( int i = 1 ; i <= wf.n1 ; i++ )
      {if(abs(x1-wf.x1[i]) < (wf.dx1) )
	  {
	    for( int l = 1 ; l <= wf.n4 ; l++ )
	      for( int k = 1 ; k <= wf.n3 ; k++ )
		for( int j = 1 ; j <= wf.n2 ; j++)
		  {
		    result[(l-1)*wf.n3*wf.n2+(k-1)*wf.n2+j-1]= norm( wf.wave[ wf.in4( l, k , j , i ) ] );	
		  }
	  }
      }
    
    return result;

  } // end of Obs_Fix_X1

  vector<double> Obs_Fix_X2 (wavefunction &wf, double x2)
  {
    vector<double> result;
    result.resize( wf.n4*wf.n3*wf.n1, 0.);

    for( int j = 1 ; j <= wf.n2 ; j++ )
      {if(abs(x2-wf.x2[j]) < (wf.dx2) )
	  {
	    for( int l = 1 ; l <= wf.n4 ; l++ )
	      for( int k = 1 ; k <= wf.n3 ; k++ )
		for( int i = 1 ; i <= wf.n1 ; i++)
		  {
		    result[(l-1)*wf.n3*wf.n1+(k-1)*wf.n1+i-1]= norm( wf.wave[ wf.in4( l, k , j , i ) ] );	
		  }
	  }
      }

    
    return result;

  } // end of Obs_Fix_X2

  vector<double> Obs_Fix_X3 (wavefunction &wf, double x3)
  {
    vector<double> result;
    result.resize( wf.n4*wf.n2*wf.n1, 0.);

    for( int k = 1 ; k <= wf.n3 ; k++ )
      {if(abs(x3-wf.x3[k]) < (wf.dx3) )
	  {
	    for( int l = 1 ; l <= wf.n4 ; l++ )
	      for( int j = 1 ; j <= wf.n2 ; j++ )
		for( int i = 1 ; i <= wf.n1 ; i++)
		  {
		    result[(l-1)*wf.n2*wf.n1+(j-1)*wf.n1+i-1]= norm( wf.wave[ wf.in4( l, k , j , i ) ] );	
		  }
	  }
      }
    
    return result;

  } // end of Obs_Fix_X3

  
  


  vector<double> Obs_Projection_X1 (wavefunction &wf)
  {
    vector<double> result;    
    result.resize(wf.n4*wf.n3*wf.n2);
    
    for( int l = 1 ; l <= wf.n4 ; l++ )
      for( int k = 1 ; k <= wf.n3 ; k++ )
	for( int j = 1 ; j <= wf.n2 ; j++)
	  {
	    double sum = 0.;
	    for( int i = 1 ; i <= wf.n1 ; i++ )
	      sum += norm( wf.wave[ wf.in4( l, k , j , i ) ] )*wf.dx1;
	    result[(l-1)*wf.n3*wf.n2+(k-1)*wf.n2+j-1]=sum;
	  }
    
    return result;
    
  } // end of Obs_Projection_X1


 vector<double> Obs_Projection_X2 (wavefunction &wf)
  {
    vector<double> result;
    result.resize(wf.n4*wf.n3*wf.n1);
    
    for( int l = 1 ; l <= wf.n4 ; l++ )
      for( int k = 1 ; k <= wf.n3 ; k++ )
	for( int i = 1 ; i <= wf.n1 ; i++)
	  {
	    double sum = 0.;
	    for( int j = 1 ; j <= wf.n2 ; j++ )
	      sum += norm( wf.wave[ wf.in4( l, k , j , i ) ] )*wf.dx2;
	    result[(l-1)*wf.n3*wf.n1+(k-1)*wf.n1+i-1]=sum;
	  }
       
    return result;

  } // end of Obs_Projection_X2



 vector<double> Obs_Projection_X3 (wavefunction &wf)
  {
    vector<double> result;
    result.resize(wf.n4*wf.n2*wf.n1);
    
    for( int l = 1 ; l <= wf.n4 ; l++ )
      for( int j = 1 ; j <= wf.n2 ; j++ )
	for( int i = 1 ; i <= wf.n1 ; i++)
	  {
	    double sum = 0.;
	    for( int k = 1 ; k <= wf.n3 ; k++ )
	      sum += norm( wf.wave[ wf.in4( l, k , j , i ) ] )*wf.dx3;
	    result[(l-1)*wf.n2*wf.n1+(j-1)*wf.n1+i-1]=sum;
	  }

    return result;

  } // end of Obs_Projection_X3

 vector<double> Obs_Projection_X4 (wavefunction &wf)
  {
    vector<double> result;
    result.resize(wf.n3*wf.n2*wf.n1);
    
    for( int k = 1 ; k <= wf.n3 ; k++ )
      for( int j = 1 ; j <= wf.n2 ; j++ )
	for( int i = 1 ; i <= wf.n1 ; i++)
	  {
	    double sum = 0.;
	    for( int l = 1 ; l <= wf.n4 ; l++ )
	      sum += norm( wf.wave[ wf.in4( l, k , j , i ) ] )*wf.dx4;
	    result[(k-1)*wf.n2*wf.n1+(j-1)*wf.n1+i-1]=sum;
	  }

    return result;

  } // end of Obs_Projection_X4


 vector<double> Obs_Projection_X1X2 (wavefunction &wf)
  {
    vector<double> result;
    
    result.resize((wf.n4)*(wf.n3));
    for( int l = 1 ; l <= wf.n4 ; l++ )
      for( int k = 1 ; k <= wf.n3 ; k++ )
        {
          double sum = 0.;
	  for( int j = 1 ; j <= wf.n2 ; j++ )
	    for( int i = 1 ; i <= wf.n1 ; i++)
	      {
		sum += norm( wf.wave[ wf.in4( l, k , j , i ) ] )*wf.dx1*wf.dx2;
	      }
          result[(l-1)*wf.n3+k-1]=sum;
        }

    return result;

  } // end of Obs_Projection_X1_X2

 vector<double> Obs_Projection_X1X3 (wavefunction &wf)
  {
    vector<double> result;   
    result.resize((wf.n4)*(wf.n2));

    for( int l = 1 ; l <= wf.n4 ; l++ )
      for( int j = 1 ; j <= wf.n2 ; j++ )
        {
          double sum = 0.;
	  for( int k = 1 ; k <= wf.n3 ; k++ )
	    for( int i = 1 ; i <= wf.n1 ; i++)
	      {
		sum += norm( wf.wave[ wf.in4( l, k , j , i ) ] )*wf.dx1*wf.dx3;
	      }
          result[(l-1)*wf.n2+j-1]=sum;
        }

    return result;

  } // end of Obs_Projection_X1_X3

 vector<double> Obs_Projection_X1X4 (wavefunction &wf)
  {
    vector<double> result;   
    result.resize((wf.n3)*(wf.n2));
    
    for( int k = 1 ; k <= wf.n3 ; k++ )    
      for( int j = 1 ; j <= wf.n2 ; j++ )
        {
          double sum = 0.;
	  for( int l = 1 ; l <= wf.n4 ; l++ )
	    for( int i = 1 ; i <= wf.n1 ; i++)
	      {
		sum += norm( wf.wave[ wf.in4( l, k , j , i ) ] )*wf.dx1*wf.dx4;
	      }
          result[(k-1)*wf.n2+j-1]=sum;
        }

    return result;

  } // end of Obs_Projection_X1_X4

 vector<double> Obs_Projection_X2X3 (wavefunction &wf)
  {
    vector<double> result;   
    result.resize((wf.n4)*(wf.n1));
 
    for( int l = 1 ; l <= wf.n4 ; l++ )
      for( int i = 1 ; i <= wf.n1 ; i++)
        {
          double sum = 0.;
	  for( int k = 1 ; k <= wf.n3 ; k++ )
	    for( int j = 1 ; j <= wf.n2 ; j++ )
	      {
		sum += norm( wf.wave[ wf.in4( l, k , j , i ) ] )*wf.dx2*wf.dx3;
	      }
          result[(l-1)*wf.n1+i-1]=sum;
        }

    return result;

  } // end of Obs_Projection_X2_X3

 vector<double> Obs_Projection_X2X4 (wavefunction &wf)
  {
    vector<double> result;   
    result.resize((wf.n3)*(wf.n1));

    for( int k = 1 ; k <= wf.n3 ; k++ )
      for( int i = 1 ; i <= wf.n1 ; i++)
        {
          double sum = 0.;
	  for( int l = 1 ; l <= wf.n4 ; l++ )
	    for( int j = 1 ; j <= wf.n2 ; j++ )
	      {
		sum += norm( wf.wave[ wf.in4( l, k , j , i ) ] )*wf.dx2*wf.dx4;
	      }
          result[(k-1)*wf.n1+i-1]=sum;
        }

    return result;

  } // end of Obs_Projection_X2_X4

 vector<double> Obs_Projection_X3X4 (wavefunction &wf)
  {
    vector<double> result;

    result.resize((wf.n2)*(wf.n1));
    for( int j = 1 ; j <= wf.n2 ; j++ )
      for( int i = 1 ; i <= wf.n1 ; i++)
        {
          double sum = 0.;
	  for( int l = 1 ; l <= wf.n4 ; l++ )
	    for( int k = 1 ; k <= wf.n3 ; k++ )
	      {
		sum += norm( wf.wave[ wf.in4( l, k , j , i ) ] )*wf.dx3*wf.dx4;
	      }
          result[(j-1)*wf.n1+i-1]=sum;
        }

    return result;

  } // end of Obs_Projection_X3_X4



  vector<double> Distribution_X1( wavefunction &wf )
  {
    vector<double> dist(wf.n1,0.);
    double sum;
    for( int i = 1 ; i <= wf.n1 ; i++) 
      {
	sum = 0.;
	for( int l = 1 ; l <= wf.n4 ; l++)
	  for( int k = 1 ; k <= wf.n3 ; k++ ) 
	    { 
	      for( int j = 1 ; j <= wf.n2 ; j++ ) 
		{ 
		  sum += norm(wf.wave[ wf.in4( l, k , j , i ) ]);
		}
	    }
	dist[i-1] = wf.dx2*wf.dx3*wf.dx4*sum;
      }
    return dist;
    
  }
  
	
	
	
  
  vector<double> Distribution_X2( wavefunction &wf )
  {
    vector<double> dist(wf.n2,0.);
    double sum;
    for( int j = 1 ; j <= wf.n2 ; j++) 
      {
	sum = 0.;
	for( int l = 1 ; l <= wf.n4 ; l++)
	  for( int k = 1 ; k <= wf.n3 ; k++ ) 
	    { 
	      for( int i = 1 ; i <= wf.n1 ; i++ ) 
		{ 
		  sum += norm(wf.wave[ wf.in4( l, k , j , i ) ]);
		}
	    }
	dist[j-1] = wf.dx1*wf.dx3*wf.dx4*sum;
      }
    return dist;
    
  }
  
  
  vector<double> Distribution_X3( wavefunction &wf )
  {
    vector<double> dist(wf.n3,0.);
    double sum;
    for( int k = 1 ; k <= wf.n3 ; k++) 
      {
	sum = 0.;
	for( int l = 1 ; l <= wf.n4 ; l++)
	  for( int j = 1 ; j <= wf.n2 ; j++ ) 
	    { 
	      for( int i = 1 ; i <= wf.n1 ; i++ ) 
		{ 
		  sum += norm(wf.wave[ wf.in4( l, k , j , i ) ]);
		}
	    }
	dist[k-1] = wf.dx1*wf.dx2*wf.dx4*sum;
      }
    return dist;
    
  }
  
	
  vector<double> Distribution_X4( wavefunction &wf )
  {
    vector<double> dist(wf.n4,0.);
    double sum;
    for( int l = 1 ; l <= wf.n4 ; l++) 
      {
	sum = 0.;
	for( int k = 1 ; k <= wf.n3 ; k++ ) 
	  for( int j = 1 ; j <= wf.n2 ; j++ ) 
	    { 
	      for( int i = 1 ; i <= wf.n1 ; i++ ) 
		{ 
		  sum += norm(wf.wave[ wf.in4(l, k , j , i ) ]);
		}
	    }
	dist[l-1] = wf.dx1*wf.dx2*wf.dx3*sum;
      }
    return dist;
    
  }



 
 /**----------	 * Calculation of Derivatives: ----------------------------*/
   




    /**
     * Calculate first derivative in x1-direction - look at Eq. (5) in doc/ABV/abv.tex
     * note the boundary condition phi( dx1/2 ) = phi( -dx1/2 ) which leads to phi( index = 0 ) = phi( index = 1 )
     */
  complex First_Derivative_X1( wavefunction &wf , int l, int k , int j , int i )
  {
    complex deriv;

    if (wf.symmetry_x1==1 && i ==1 )//symmetry
      {
	
	deriv = ( wf.wave[ wf.in4( l, k , j , i+1 ) ]-wf.wave[ wf.in4( l, k , j , i ) ] ) / ( 2.0*wf.dx1 );
	
      }
    if (wf.symmetry_x1==-1 && i ==1 )//antisymmetry
      {
	
	deriv = ( wf.wave[ wf.in4( l, k , j , i+1 ) ]+wf.wave[ wf.in4( l, k , j , i ) ] ) / ( 2.0*wf.dx1 );
	
      }
    else
      {
	deriv = ( wf.wave[ wf.in4( l, k , j , i+1 ) ]-wf.wave[ wf.in4( l, k , j , i-1 ) ] ) / ( 2.0*wf.dx1 );
      }//sym==1
    return deriv;

  } // end of First_Derivative_X1
  


    /**
     * Calculate first derivative in x2-direction - look at Eq. (5) in doc/ABV/abv.tex
     */
  complex First_Derivative_X2( wavefunction &wf , int l, int k , int j , int i )
  {
    complex deriv;

    if (wf.symmetry_x2==1 && j==1)//symmetry
      {
	deriv = ( wf.wave[ wf.in4( l, k , j+1 , i ) ]-wf.wave[ wf.in4( l, k , j , i ) ] ) / ( 2.0*wf.dx2 );
	
      }
    if (wf.symmetry_x2==-1 && j==1)//antisymmetry
      {
	deriv = ( wf.wave[ wf.in4( l, k , j+1 , i ) ]+wf.wave[ wf.in4( l, k , j , i ) ] ) / ( 2.0*wf.dx2 );
	
      }
    else
      {
	deriv = ( wf.wave[ wf.in4( l, k , j+1 , i ) ]-wf.wave[ wf.in4( l, k , j-1 , i ) ] ) / ( 2.0*wf.dx2 );
      }

    return deriv;
  } // end of First_Derivative_X2
  


    /**
     * Calculate first derivative in x3-direction - look at Eq. (5) in doc/ABV/abv.tex
     */
  complex First_Derivative_X3( wavefunction &wf , int l, int k , int j , int i )
  {
    complex deriv;

    if (wf.symmetry_x3==1 && k==1)//symmetry
      {
	deriv = ( wf.wave[ wf.in4( l, k+1 , j , i ) ]-wf.wave[ wf.in4( l, k , j , i ) ] ) / ( 2.0*wf.dx3 );
	
      }
    if (wf.symmetry_x3==-1 && k==1)//antisymmetry
      {
	deriv = ( wf.wave[ wf.in4( l, k+1 , j , i ) ]+wf.wave[ wf.in4( l, k , j , i ) ] ) / ( 2.0*wf.dx3 );
	
      }
    else
      {
	deriv = ( wf.wave[ wf.in4( l, k+1 , j , i ) ]-wf.wave[ wf.in4( l, k-1 , j , i ) ] ) / ( 2.0*wf.dx3 );
      }//sym==1

    return deriv;
  } // end of First_Derivative_X3
  


    /**
     * Calculate second derivative in x1-direction - look at Eq. (5) in doc/ABV/abv.tex
     * note the boundary condition phi( dx1/2 ) = phi( -dx1/2 ) which leads to phi( index = 0 ) = phi( index = 1 )
     */
  complex Second_Derivative_X1( wavefunction &wf , int l, int k , int j , int i )
  {

    if (i==1) {
    if (wf.symmetry_x1==1)//symmetry
      {
	return ( wf.wave[ wf.in4( l, k , j , i+1 ) ]-wf.wave[ wf.in4( l, k , j , i ) ] ) * wf.one_by_dx1sqr;
      }
    if (wf.symmetry_x1==-1)//antisymmetry
      {
	return ( wf.wave[ wf.in4( l, k , j , i+1 ) ]-3.*wf.wave[ wf.in4( l, k , j , i ) ] ) * wf.one_by_dx1sqr;
      }
    //if (wf.symmetry_x1==1)
    } else {

      return ( wf.wave[ wf.in4( l, k , j , i+1 ) ]-2.*wf.wave[ wf.in4( l, k , j , i ) ] + wf.wave[ wf.in4( l, k , j , i-1 ) ] ) * wf.one_by_dx1sqr;//sym==1

    }
      
  } // end of Second_Derivative_X1
  


    /**
     * Calculate second derivative in x2-direction - look at Eq. (5) in doc/ABV/abv.tex
     */
  complex Second_Derivative_X2( wavefunction &wf , int l, int k , int j , int i )
  {

    if (j==1) {
    if (wf.symmetry_x2==1)//symmetry
      {
	return ( wf.wave[ wf.in4( l, k , j+1 , i ) ]-wf.wave[ wf.in4( l, k , j , i ) ] ) * wf.one_by_dx2sqr;
	
      }
    if (wf.symmetry_x2==-1)//antisymmetry
      {
	return ( wf.wave[ wf.in4( l, k , j+1 , i ) ]-3.*wf.wave[ wf.in4( l, k , j , i ) ] ) * wf.one_by_dx2sqr;
	
      }
    // if (wf.symmetry_x2==1)
    } else {

      return ( wf.wave[ wf.in4( l, k , j+1 , i ) ]-2.*wf.wave[ wf.in4( l, k , j , i ) ] + wf.wave[ wf.in4( l, k , j-1 , i ) ] ) * wf.one_by_dx2sqr;  //sym==1

    }

  } // end of Second_Derivative_X2
  


    /**
     * Calculate second derivative in x3-direction - look at Eq. (5) in doc/ABV/abv.tex
     */
  complex Second_Derivative_X3( wavefunction &wf , int l, int k , int j , int i )
  {

    if (k==1) {
    if (wf.symmetry_x3==1)//symmetry
      {
	return ( wf.wave[ wf.in4( l, k+1 , j , i ) ]-wf.wave[ wf.in4( l, k , j , i ) ] ) * wf.one_by_dx3sqr;
	
      }
    if (wf.symmetry_x3==-1)//antisymmetry
      {
	return ( wf.wave[ wf.in4( l, k+1 , j , i ) ]-3.*wf.wave[ wf.in4( l, k , j , i ) ] ) * wf.one_by_dx3sqr;
	
      }
    // if (wf.symmetry_x3==1)
    } else {

      return ( wf.wave[ wf.in4( l, k+1 , j , i ) ]-2.*wf.wave[ wf.in4( l, k , j , i ) ] + wf.wave[ wf.in4( l, k-1 , j , i ) ] ) * wf.one_by_dx3sqr;  //sym==1

    }

  } // end of Second_Derivative_X3

  

  

  complex Second_Derivative_X4_MPI( wavefunction &wf , int l, int k , int j , int i )
  {


    int tag3=3, tag4=4;
    int left, right, myid, numprocs;
    int send_num=1;
    MPI_Status status;
   
    int nlocal= wf.n4;
    MPI_Comm_rank (MPI_COMM_WORLD, &myid);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);

    if ( myid > 0 )
      left=myid-1;
    else 
      left=MPI_PROC_NULL;
    if ( myid < numprocs-1  )
      right=myid+1;
    else
      right=MPI_PROC_NULL;  //// define MPI_NULL buffer

    // //////////////**************** arrange psi
    
   
    MPI_Sendrecv(&wf.wave[wf.in4(nlocal,k,j,i)],send_num,MPI_DOUBLE_COMPLEX,right,tag3,&wf.wave[wf.in4(0,k,j,i)],send_num,MPI_DOUBLE_COMPLEX,left,tag3,MPI_COMM_WORLD, &status);
    MPI_Sendrecv(&wf.wave[wf.in4(1,k,j,i)],send_num,MPI_DOUBLE_COMPLEX,left,tag4,&wf.wave[wf.in4(nlocal+1,k,j,i)],send_num,MPI_DOUBLE_COMPLEX,right,tag4,MPI_COMM_WORLD, &status);
   
       
    MPI_Barrier(MPI_COMM_WORLD);

    if (myid==0&&l==1) 
      {
	if (wf.symmetry_x4==1)//symmetry
	  {
	    return ( wf.wave[ wf.in4( l+1 , k , j , i ) ]-wf.wave[ wf.in4( l , k , j , i ) ] ) * wf.one_by_dx4sqr;
	
	  }
	if (wf.symmetry_x4==-1)//antisymmetry
	  {
	    return ( wf.wave[ wf.in4( l+1 , k , j , i ) ]-3.*wf.wave[ wf.in4( l , k , j , i ) ] ) * wf.one_by_dx4sqr;
	
	  }
    // if (wf.symmetry_x3==1)
      } 
    else
      {

	return ( wf.wave[ wf.in4( l+1 , k , j , i ) ]-2.*wf.wave[ wf.in4( l , k , j , i ) ] + wf.wave[ wf.in4( l-1 , k , j , i ) ] ) * wf.one_by_dx4sqr;  //sym==1
      }
  
    // cout<<"l"<<l<<"second_derive "<<( wf.wave[ wf.in4( l+1 , k , j , i ) ]-2.*wf.wave[ wf.in4( l , k , j , i ) ] + wf.wave[ wf.in4( l-1 , k , j , i ) ] ) * wf.one_by_dx4sqr<<endl;
  } // end of Second_Derivative_X4_MPI


  




  /*------------------	Stuff needed for FFT and Masks	-----------------------*/











  /**
   * mask function for x1, x2, and x3. frac_n*_* are the fraction for absorbing in three corordinates.
   */


  void mask_kernel(wavefunction &wf, int i, int j, int k, int l, double frac_n1_right, double frac_n2_right, double frac_n3_right, double frac_n4_right, double frac_n1_left, double frac_n2_left, double frac_n3_left, double frac_n4_left, double exponent )
  {
    int myid, numprocs;
    MPI_Comm_rank (MPI_COMM_WORLD, &myid);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);

    double mask_start_x1_right=wf.x1[int((wf.n1)*(1.-frac_n1_right))];
    double mask_start_x2_right=wf.x2[int((wf.n2)*(1.-frac_n2_right))];
    double mask_start_x3_right=wf.x3[int((wf.n3)*(1.-frac_n3_right))];
    double mask_start_x4_right=wf.x4[int((wf.n4)*(1.-frac_n4_right))];
    double mask_start_x1_left=wf.x1[int((wf.n1)*frac_n1_left)+1];
    double mask_start_x2_left=wf.x2[int((wf.n2)*frac_n2_left)+1];
    double mask_start_x3_left=wf.x3[int((wf.n3)*frac_n3_left)+1];
    double mask_start_x4_left=wf.x4[int((wf.n4)*frac_n4_left)+1];

    double argument_x1_right;
    double argument_x2_right;
    double argument_x3_right;
    double argument_x4_right;
    double argument_x1_left;
    double argument_x2_left;
    double argument_x3_left;
    double argument_x4_left;

    double mask_x1_right;
    double mask_x2_right;
    double mask_x3_right;
    double mask_x4_right;
    double mask_x1_left;
    double mask_x2_left;
    double mask_x3_left;
    double mask_x4_left;

    argument_x1_right=(pi/2.)*(wf.x1[i]-mask_start_x1_right)/(wf.x1[wf.n1]-mask_start_x1_right+1.e-20);
    argument_x2_right=(pi/2.)*(wf.x2[j]-mask_start_x2_right)/(wf.x2[wf.n2]-mask_start_x2_right+1.e-20);
    argument_x3_right=(pi/2.)*(wf.x3[k]-mask_start_x3_right)/(wf.x3[wf.n3]-mask_start_x3_right+1.e-20);
    argument_x4_right=(pi/2.)*(wf.x4[l]-mask_start_x4_right)/(wf.x4[wf.n4]-mask_start_x4_right+1.e-20);
    
    argument_x1_left=(pi/2.)*(wf.x1[i]-mask_start_x1_left)/(wf.x1[1]-mask_start_x1_left+1.e-20);
    argument_x2_left=(pi/2.)*(wf.x2[j]-mask_start_x2_left)/(wf.x2[1]-mask_start_x2_left+1.e-20);
    argument_x3_left=(pi/2.)*(wf.x3[k]-mask_start_x3_left)/(wf.x3[1]-mask_start_x3_left+1.e-20);
    argument_x4_left=(pi/2.)*(wf.x4[l]-mask_start_x4_left)/(wf.x4[1]-mask_start_x4_left+1.e-20);

    mask_x1_right=pow(fabs(cos(argument_x1_right)),exponent);
    mask_x2_right=pow(fabs(cos(argument_x2_right)),exponent);
    mask_x3_right=pow(fabs(cos(argument_x3_right)),exponent);
    mask_x4_right=pow(fabs(cos(argument_x4_right)),exponent);
    
    mask_x1_left=pow(fabs(cos(argument_x1_left)),exponent);
    mask_x2_left=pow(fabs(cos(argument_x2_left)),exponent);
    mask_x3_left=pow(fabs(cos(argument_x3_left)),exponent);
    mask_x4_left=pow(fabs(cos(argument_x4_left)),exponent);
    
    
    if (wf.symmetry_x1==0)
      {
	if (i< int(wf.n1*frac_n1_left))		  
	  wf.wave[wf.in4(l,k,j,i)]*=mask_x1_left;		  
      }//if sym==1 or sym==-1, we needn't this part.

    if (i> int(wf.n1*(1.-frac_n1_right)))
      {
	wf.wave[wf.in4(l,k,j,i)]*=mask_x1_right;
      }
    
    if (wf.symmetry_x2==0)
      {
	if (j< int(wf.n2*frac_n2_left))
	  wf.wave[wf.in4(l,k,j,i)]*=mask_x2_left;
      }//if sym==1 or sym==-1, we needn't this part.
    
    if (j> int(wf.n2*(1.-frac_n2_right)))
      wf.wave[wf.in4(l,k,j,i)]*=mask_x2_right;
    
    if (wf.symmetry_x3==0)
      {
	if (k< int(wf.n3*frac_n3_left))
	  wf.wave[wf.in4(l,k,j,i)]*=mask_x3_left;
      }//if sym==1 or sym==-1, we needn't this part.
    
    
    if (k> int(wf.n3*(1.-frac_n3_right)))
      wf.wave[wf.in4(l,k,j,i)]*=mask_x3_right;

    if ((wf.symmetry_x4==0)&&(myid==0))
      {
	if (l< int(wf.n4*frac_n4_left))
	  wf.wave[wf.in4(l,k,j,i)]*=mask_x4_left;
      }//if sym==1 or sym==-1, we needn't this part.
    
    
    if ( (l> int(wf.n4*(1.-frac_n4_right)))&&(myid==numprocs-1) )
      wf.wave[wf.in4(l,k,j,i)]*=mask_x4_right;

  }


  void Mask_Function(wavefunction &wf, double frac_n1_right, double frac_n2_right, double frac_n3_right, double frac_n4_right, double frac_n1_left, double frac_n2_left, double frac_n3_left, double frac_n4_left,double exponent, const int nthreads )
  {
    if ( frac_n1_right>1. || frac_n1_right<0. )
      {
	cout<< "frac_n1_right is wrong";
	exit(1);
      }
    
    if ( frac_n2_right>1. || frac_n2_right<0. )
      {
	cout<< "frac_n2_right is wrong";
	exit(1);
      }

    if ( frac_n3_right>1. || frac_n3_right<0. )
      {
	cout<< "frac_n3_right is wrong";
	exit(1);
      }
    if ( frac_n4_right>1. || frac_n4_right<0. )
      {
	cout<< "frac_n4_right is wrong";
	exit(1);
      }
    if ( frac_n1_left>1. || frac_n1_left<0. )
      {
	cout<< "frac_n1_left is wrong";
	exit(1);
      }
    if ( frac_n2_left>1. || frac_n2_left<0. )
      {
	cout<< "frac_n2_left is wrong";
	exit(1);
      }

    if ( frac_n3_left>1. || frac_n3_left<0. )
      {
	cout<< "frac_n3_left is wrong";
	exit(1);
      }
    if ( frac_n4_left>1. || frac_n4_left<0. )
      {
	cout<< "frac_n4_left is wrong";
	exit(1);
      }
   
   
    // int nthreads=16;
    //double wall_timer=omp_get_wtime();
    //clock_t clock_timer=clock();
#pragma omp parallel for num_threads(nthreads) schedule(static)  collapse(4)
    for(int l=1;l<=wf.n4;l++)
      for(int k=1;k<=wf.n3;k++)
	for(int j=1;j<=wf.n2;j++)
	  for(int i=1;i<=wf.n1;i++)
	    mask_kernel( wf, i, j, k, l, frac_n1_right, frac_n2_right, frac_n3_right, frac_n4_right, frac_n1_left, frac_n2_left, frac_n3_left, frac_n4_left, exponent );
  }



 

 







 

  void FFT_n1( wavefunction &wf , wavefunction &wf_transformed )
  {
  
    double trafo_factor1=sqrt(2.0*pi)/(wf_transformed.n1*wf_transformed.dx1);

  
    complex tmp;
    fftw_complex *in=new fftw_complex[wf.n1];
    fftw_plan p1 = fftw_plan_dft_1d(wf.n1,&in[0], &in[0],FFTW_FORWARD, FFTW_ESTIMATE);
    for(int l=1;l<=wf.n4;l++)
      for(int j=1;j<=wf.n2;j++)
	for(int k=1;k<=wf.n3;k++)
	  {	    
	    for(int i=1;i<=wf.n1;i++) {
	      int fft_index1= i-1;
	      complex element=wf.wave[wf.in4(l,k,j,i)];
	      in[fft_index1][0]=real(element);
	      in[fft_index1][1]=imag(element);
	    }
	    
	    fftw_execute(p1);
	    
	    for(int i=1;i<=wf.n1;i++) {
	      int fft_index1= i-1;
	      int wf_index=wf_transformed.in4(l,k,j,i);
	      wf_transformed.wave[wf_index] = in[fft_index1][0]+I*in[fft_index1][1];
	      wf_transformed.wave[wf_index]*= trafo_factor1;
	    }	    
	  }
    fftw_destroy_plan(p1);
    delete [] in;
  }
  
  void FFT_n1_shift( wavefunction &wf , wavefunction &wf_transformed )
  {
    
    double trafo_factor1=sqrt(2.0*pi)/(wf_transformed.n1*wf_transformed.dx1);   /// original

    complex tmp;
    fftw_complex *in=new fftw_complex[wf.n1];
    fftw_plan p1 = fftw_plan_dft_1d(wf.n1,&in[0], &in[0],FFTW_FORWARD, FFTW_ESTIMATE);
    for(int l=1;l<=wf.n4;l++)
    for(int j=1;j<=wf.n2;j++)
      for(int k=1;k<=wf.n3;k++)
	{
        for(int i=1;i<=wf.n1;i++) {
          int fft_index1= i-1;
          complex element=wf.wave[wf.in4(l,k,j,i)];
          in[fft_index1][0]=real(element);
          in[fft_index1][1]=imag(element);
        }

        fftw_execute(p1);

        for(int i=1;i<=wf.n1;i++) {
          int fft_index1= i-1;
          int wf_index=wf_transformed.in4(l,k,j,i);
          wf_transformed.wave[wf_index] = in[fft_index1][0]+I*in[fft_index1][1];
          wf_transformed.wave[wf_index]*= trafo_factor1;
        }

        for(int i=1;i<=wf.n1/2;i++)
          {
	    tmp=wf_transformed.wave[wf_transformed.in4(l,k,j,i+wf.n1/2)];
	    wf_transformed.wave[wf_transformed.in4(l,k,j,i+wf.n1/2)]=wf_transformed.wave[wf_transformed.in4(l,k,j,i)];;
	    wf_transformed.wave[wf_transformed.in4(l,k,j,i)]=tmp;
          }
        
      }

    fftw_destroy_plan(p1);
    delete [] in;
  }


void FFT_n2( wavefunction &wf , wavefunction &wf_transformed )
  {

    double trafo_factor1=sqrt(2.0*pi)/(wf_transformed.n2*wf_transformed.dx2);
    fftw_complex *in=new fftw_complex[wf.n2];

    complex tmp;
    fftw_plan p1 = fftw_plan_dft_1d(wf.n2,&in[0], &in[0],FFTW_FORWARD, FFTW_ESTIMATE);

    for(int l=1;l<=wf.n4;l++)
      for(int k=1;k<=wf.n3;k++)
	for(int i=1;i<=wf.n1;i++)
	  {	   
	    for(int j=1;j<=wf.n2;j++) {
	      int fft_index1= j-1;
	      complex element=wf.wave[wf.in4(l,k,j,i)];
	      in[fft_index1][0]=real(element);
	      in[fft_index1][1]=imag(element);
	    }
	    fftw_execute(p1);
	    for(int j=1;j<=wf.n2;j++) {
	      int fft_index1= j-1;
	      int wf_index=wf_transformed.in4(l,k,j,i);
	      wf_transformed.wave[wf_index] = in[fft_index1][0]+I*in[fft_index1][1];
	      wf_transformed.wave[wf_index]*= trafo_factor1;
	    }
	    
	  }

    fftw_destroy_plan(p1);
    delete [] in;
  }




void FFT_n2_shift( wavefunction &wf , wavefunction &wf_transformed )
  {
 
    double trafo_factor1=sqrt(2.0*pi)/(wf_transformed.n2*wf_transformed.dx2);
    fftw_complex *in=new fftw_complex[wf.n2];

    complex tmp;
    fftw_plan p1= fftw_plan_dft_1d(wf.n2,&in[0], &in[0],FFTW_FORWARD, FFTW_ESTIMATE);

    for(int l=1;l<=wf.n4;l++)
    for(int k=1;k<=wf.n3;k++)
      for(int i=1;i<=wf.n1;i++)
        {        
          for(int j=1;j<=wf.n2;j++) {
            int fft_index1= j-1;
            complex element=wf.wave[wf.in4(l,k,j,i)];
	    in[fft_index1][0]=real(element);
	    in[fft_index1][1]=imag(element);
          }
          fftw_execute(p1);
          for(int j=1;j<=wf.n2;j++) {
            int fft_index1= j-1;
            int wf_index=wf_transformed.in4(l,k,j,i);
            wf_transformed.wave[wf_index] = in[fft_index1][0]+I*in[fft_index1][1];
            wf_transformed.wave[wf_index]*= trafo_factor1;
          }

          for(int j=1;j<=wf.n2/2;j++)
            {
              tmp=wf_transformed.wave[wf_transformed.in4(l,k,j+wf.n2/2,i)];
              wf_transformed.wave[wf_transformed.in4(l,k,j+wf.n2/2,i)]=wf_transformed.wave[wf_transformed.in4(l,k,j,i)];;
              wf_transformed.wave[wf_transformed.in4(l,k,j,i)]=tmp;
            }
	 
        }
    
    fftw_destroy_plan(p1);
    delete [] in;
  }





void FFT_n3( wavefunction &wf , wavefunction &wf_transformed )
{

    double trafo_factor1=sqrt(2.0*pi)/(wf_transformed.n3*wf_transformed.dx3);

    complex tmp;
    fftw_complex *in=new fftw_complex[wf.n3];
    fftw_plan p1 = fftw_plan_dft_1d(wf.n3,&in[0], &in[0],FFTW_FORWARD, FFTW_ESTIMATE);

    for(int l=1;l<=wf.n4;l++)
    for(int j=1;j<=wf.n2;j++)
      for(int i=1;i<=wf.n1;i++)
	{
        for(int k=1;k<=wf.n3;k++) {
          int fft_index1= k-1;
          complex element=wf.wave[wf.in4(l,k,j,i)];
          in[fft_index1][0]=real(element);
          in[fft_index1][1]=imag(element);
        }

        fftw_execute(p1);

        for(int k=1;k<=wf.n3;k++) {
          int fft_index1= k-1;
          int wf_index=wf_transformed.in4(l,k,j,i);
          wf_transformed.wave[wf_index] = in[fft_index1][0]+I*in[fft_index1][1];
          wf_transformed.wave[wf_index]*= trafo_factor1;
        }

        for(int k=1;k<=wf.n3/2;k++)
          {
	    tmp=wf_transformed.wave[wf_transformed.in4(l,k+wf.n3/2,j,i)];
	    wf_transformed.wave[wf_transformed.in4(l,k+wf.n3/2,j,i)]=wf_transformed.wave[wf_transformed.in4(l,k,j,i)];;
	    wf_transformed.wave[wf_transformed.in4(l,k,j,i)]=tmp;
          } 
       
      }
    fftw_destroy_plan(p1);
    delete [] in;
  }

  void FFT_n3_shift( wavefunction &wf , wavefunction &wf_transformed )
  {
   
    double trafo_factor1=sqrt(2.0*pi)/(wf_transformed.n3*wf_transformed.dx3); // original
    complex tmp;
    fftw_complex *in=new fftw_complex[wf.n3];
    fftw_plan p1 = fftw_plan_dft_1d(wf.n3,&in[0], &in[0],FFTW_FORWARD, FFTW_ESTIMATE);
    
    for(int l=1;l<=wf.n4;l++)
    for(int j=1;j<=wf.n2;j++)
      for(int i=1;i<=wf.n1;i++)
	{ 
	  for(int k=1;k<=wf.n3;k++) {
          int fft_index1= k-1;
          complex element=wf.wave[wf.in4(l,k,j,i)];
          in[fft_index1][0]=real(element);
          in[fft_index1][1]=imag(element);
        }
	
        fftw_execute(p1);
	
        for(int k=1;k<=wf.n3;k++) {
          int fft_index1= k-1;
          int wf_index=wf_transformed.in4(l,k,j,i);
          wf_transformed.wave[wf_index] = in[fft_index1][0]+I*in[fft_index1][1];
          wf_transformed.wave[wf_index]*= trafo_factor1;
        }
	
	for(int k=1;k<=wf.n3/2;k++)
          {
	    tmp=wf_transformed.wave[wf_transformed.in4(l,k+wf.n3/2,j,i)];
	    wf_transformed.wave[wf_transformed.in4(l,k+wf.n3/2,j,i)]=wf_transformed.wave[wf_transformed.in4(l,k,j,i)];
	    wf_transformed.wave[wf_transformed.in4(l,k,j,i)]=tmp;
          }
        
	}

    fftw_destroy_plan(p1);
    delete [] in;
  }

void FFT_n3n4( wavefunction &wf , wavefunction &wf_transformed )
{
  

  int myid, numprocs;
  fftw_mpi_init();
  MPI_Comm_rank (MPI_COMM_WORLD, &myid);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  
  double trafo_factor=(2.0*pi)/(wf_transformed.n3*wf_transformed.n4*wf_transformed.dx3*wf_transformed.dx4*numprocs);

  const ptrdiff_t N0 = wf.n4*numprocs, N1 = wf.n3;
  fftw_plan plan;
  fftw_complex *data;
  ptrdiff_t alloc_local, local_n0, local_0_start;
 



  local_n0= wf.n4;
  local_0_start=myid*wf.n4;
  /* get local data size and allocate */
  alloc_local = fftw_mpi_local_size_2d(N0, N1, MPI_COMM_WORLD,
                                              &local_n0, &local_0_start);
  data = fftw_alloc_complex(alloc_local);
  
  /* create plan for in-place forward DFT */
  plan = fftw_mpi_plan_dft_2d(N0, N1, data, data, MPI_COMM_WORLD,
			      FFTW_FORWARD, FFTW_ESTIMATE);

  /* initialize data  */
  for(int j=1;j<=wf.n2;j++)
    for(int i=1;i<=wf.n1;i++)
      {
	for (int l = 1; l <= local_n0; ++l) 
	  for (int k = 1; k <= N1; ++k)
	    { 
	      data[(l-1)*N1 + k-1][0] = real(wf.wave[wf.in4(l,k,j,i)]);
	      data[(l-1)*N1 + k-1][1] = imag(wf.wave[wf.in4(l,k,j,i)]);
	    }
	/* compute transforms, in-place, as many times as desired */
	fftw_execute(plan);
	for (int l = 1; l <= local_n0; ++l) 
	  for (int k = 1; k <= N1; ++k)
	    {
	      wf_transformed.wave[wf.in4(l,k,j,i)]=data[(l-1)*N1 + k-1][0]+I*data[(l-1)*N1 + k-1][1];
	      wf_transformed.wave[wf.in4(l,k,j,i)]*=trafo_factor;
	    }
      }

  fftw_destroy_plan(plan);
  delete [] data;
}

  




  void fft_correlate_shift(vector<complex> &auto_relate, string outputFolder, double dt)
  {
    ofstream out((outputFolder+"/auto_eng.txt").c_str());
  
    complex tmp;
    fftw_complex *in=new fftw_complex[auto_relate.size() ];
    double trafo_factor = 1./auto_relate.size();

    fftw_plan p1 = fftw_plan_dft_1d(auto_relate.size(), &in[0], &in[0],FFTW_FORWARD, FFTW_ESTIMATE);
    for (int i=1; i<=auto_relate.size(); i++)
      { 
	int fft_index=i-1;
	in[fft_index][0]=auto_relate[fft_index].real();
	in[fft_index][1]=auto_relate[fft_index].imag();
      }
    fftw_execute(p1);
    fftw_destroy_plan(p1);

    for(int j=1;j<=auto_relate.size();j++) 
      {
	int fft_index1= j-1;
	auto_relate[fft_index1] = in[fft_index1][0]+I*in[fft_index1][1];
	auto_relate[fft_index1]*= trafo_factor;
      }
    
    for(int j=1;j<=auto_relate.size()/2;j++) 
      {
	int fft_index1= j-1;
	tmp = auto_relate[fft_index1];
	auto_relate[fft_index1] = auto_relate[fft_index1+auto_relate.size()/2];
	auto_relate[fft_index1+auto_relate.size()/2]= tmp;
      }
    double dw=2.*pi/auto_relate.size()/dt;
    
    for( size_t step = 0 ; step < auto_relate.size(); step++ ) //plot the energy spectra
      { 
	out << ( step-auto_relate.size()/2. )*dw <<"\t"<< abs(auto_relate[step])
	    <<"\t" <<arg(auto_relate[step]) << endl;
      }

    delete [] in;
  }



   
} 
