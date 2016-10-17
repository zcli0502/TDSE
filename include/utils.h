#ifndef UTILS_H
#define UTILS_H
#include <math.h>
#include <stdio.h>
#include <vector>
#include <complex>
#include "wavefunction.h"
#include <omp.h>
#include "mpi.h"
//#include <fstream>
//#include <iostream>
#define complex complex<double>

using namespace std;

inline void nrerror ( char error_text[] )
{
    fprintf( stderr, "numerical recipes runtime error..\n" );
    fprintf( stderr, "%s\n", error_text );
    fprintf( stderr,"... now exiting system..\n" );
}

/**
 *  The tridag routine to solve the threediagonal system of linear equation in Eq. (6) in doc/ABV/abv.tex
 */
inline void Tridag( const vector<complex> &tridag_lower_diagonal , const vector<complex> &tridag_middle_diagonal , const vector<complex>&tridag_upper_diagonal, const vector<complex> &wf_1D_rightside ,
	     vector<complex> &wf_1D_tmp2 , vector<complex> &gam )
{   
  const unsigned int length=gam.size() -2 ;

/**
 *  forth substitution
 */

  complex bet=tridag_middle_diagonal[ 1 ];
  wf_1D_tmp2[ 1 ] = wf_1D_rightside[ 1 ]/bet;
  for( unsigned int index = 2 ; index <= length; ++index )
    {
      gam[ index ] = tridag_upper_diagonal[ index-1 ]/bet;
      bet  = tridag_middle_diagonal[ index ];
      bet -= tridag_lower_diagonal[ index ]*gam[ index ];
      /*      if ( bet == 0.0 )
	{
	 cout<<"Error: bet==0 in tridag! at "<<index<<"\n";
         cout<<"tridag_middle_diagonal="<<tridag_middle_diagonal[ index ]<<"\n";
         cout<<"tridag_lower_diagonal="<<tridag_lower_diagonal[ index ]<<"\n";
         cout<<"gam="<<gam[ index ]<<"\n";
	 } */

      wf_1D_tmp2[ index ] = ( wf_1D_rightside[ index ]-tridag_lower_diagonal[ index ]*wf_1D_tmp2[ index-1 ] )/bet;
    }
  
/**
 *  back substitution
 */
  for( unsigned int index = length-1 ; index >= 1 ; --index )
    {
      wf_1D_tmp2[ index ] -= gam[ index+1 ]*wf_1D_tmp2[ index+1 ];
    }
  
  //  free( gam );
}// end of tridag


inline void Tridag_Fast( const complex &tridag_lower_diagonal , const vector<complex> &tridag_middle_diagonal , const complex &tridag_upper_diagonal, const vector<complex> &wf_1D_rightside ,
		  vector<complex> &wf_1D_tmp2 , vector<complex> &gam )
{   
  const unsigned int length=gam.size() - 2;

/*  if ( tridag_middle_diagonal[ 1 ] == 0. )
    {
      nerror( "Error: tridag_middle_diagonal[1]==0 in tridag!" );
    }
*/

  complex bet=tridag_middle_diagonal[1];
  wf_1D_tmp2[ 1 ] = wf_1D_rightside[ 1 ]/bet;

  
  for( unsigned int index = 2 ; index <= length ; ++index )
    {
      gam[ index ] = tridag_upper_diagonal/bet;
      bet  = tridag_middle_diagonal[ index ];
      bet -= tridag_lower_diagonal*gam[ index ];
/*      if (bet==0.0)
	{
	nerror("Error: bet==0 in tridag!");
	}
*/
      wf_1D_tmp2[ index ] = ( wf_1D_rightside[ index ]-tridag_lower_diagonal*wf_1D_tmp2[ index-1 ] )/bet;
    }
  
  for( unsigned int index = length-1 ; index >= 1 ; --index )
    {
      wf_1D_tmp2[ index ] -= gam[ index+1 ]*wf_1D_tmp2[ index+1 ];
    }
  
 
}// end of tridag


inline void Tridag_Fast_MPI( const complex &tridag_lower_diagonal , const vector<complex> &tridag_middle_diagonal , const complex &tridag_upper_diagonal, const vector<complex> &wf_1D_rightside ,
			     vector<complex> &wf_1D_tmp2 , vector<complex> &gam, int nlocal )
{   
  const unsigned int length=gam.size() - 2;

/*  if ( tridag_middle_diagonal[ 1 ] == 0. )
    {
      nerror( "Error: tridag_middle_diagonal[1]==0 in tridag!" );
    }
*/

  int tag1=3, tag2=4, tag3=5, tag4=6;
  int left, right, myid, numprocs;
  int send_num=1;
  MPI_Status status;
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


  complex u_from_left, bet_from_left,u_to_right, bet_to_right;
  complex bet;
  complex gam_from_right, u_from_right, gam_to_left, u_to_left;

  if (myid==0)
    {
      bet=tridag_middle_diagonal[1];
      wf_1D_tmp2[ 1 ] = wf_1D_rightside[ 1 ]/bet;
    }
  else
    {
      MPI_Recv(&u_from_left,1,MPI_DOUBLE_COMPLEX,left,tag1,MPI_COMM_WORLD, &status);
      MPI_Recv(&bet_from_left,1,MPI_DOUBLE_COMPLEX,left,tag2,MPI_COMM_WORLD, &status);
      
      gam[1]=tridag_upper_diagonal/bet_from_left;
      bet  = tridag_middle_diagonal[ 1 ];
      bet -= tridag_lower_diagonal*gam[ 1 ];
      wf_1D_tmp2[ 1 ] = ( wf_1D_rightside[ 1 ]-tridag_lower_diagonal*u_from_left )/bet;
          
    }
 
  for( unsigned int index = 2 ; index <= length ; ++index )
    {
      gam[ index ] = tridag_upper_diagonal/bet;
      bet  = tridag_middle_diagonal[ index ];
      bet -= tridag_lower_diagonal*gam[ index ];
/*      if (bet==0.0)
	{
	nerror("Error: bet==0 in tridag!");
	}
*/
      wf_1D_tmp2[ index ] = ( wf_1D_rightside[ index ]-tridag_lower_diagonal*wf_1D_tmp2[ index-1 ] )/bet;
    }
  

  if (myid < numprocs-1 )
    {
      u_to_right=wf_1D_tmp2[nlocal];
      bet_to_right=bet;
      MPI_Send(&u_to_right,1,MPI_DOUBLE_COMPLEX,right,tag1,MPI_COMM_WORLD);
      MPI_Send(&bet_to_right,1,MPI_DOUBLE_COMPLEX,right,tag2,MPI_COMM_WORLD);
    }

  ///////////////////////////////////  Backward substitution below

 
  if(myid != numprocs-1)
    {
      MPI_Recv(&gam_from_right,1,MPI_DOUBLE_COMPLEX,right,tag3,MPI_COMM_WORLD, &status);
      MPI_Recv(&u_from_right,1,MPI_DOUBLE_COMPLEX,right,tag4,MPI_COMM_WORLD, &status);
      wf_1D_tmp2[ nlocal ] -= gam_from_right*u_from_right;    
    }
  
  for( int index = nlocal-1 ; index >= 1 ; index-- )
    {
      wf_1D_tmp2[ index ] -= gam[ index+1 ]*wf_1D_tmp2[ index+1 ];
    }

  if (myid != 0)
    {
      gam_to_left = gam[1];
      u_to_left = wf_1D_tmp2[ 1 ];

      MPI_Send(&gam_to_left,1,MPI_DOUBLE_COMPLEX,left,tag3,MPI_COMM_WORLD);
      MPI_Send(&u_to_left,1,MPI_DOUBLE_COMPLEX,left,tag4,MPI_COMM_WORLD);
    }
  
 
}// end of tridag_MPI

inline void Tridag_X3_MPI(wavefunction &wf,  const complex &tridag_low, const vector<complex> &tridag_mid, const complex &tridag_upp, const vector<complex> &wf_rightside, const int nthreads )
{   

/*  if ( tridag_middle_diagonal[ 1 ] == 0. )
    {
      nerror( "Error: tridag_middle_diagonal[1]==0 in tridag!" );
    }
*/

  int tag1=3, tag2=4, tag3=5, tag4=6;
  int left, right, myid, numprocs;
  int num_comm=(wf.n1+2)*(wf.n2+2)*(wf.n4+2);
  int num_data=(wf.n1+2)*(wf.n2+2)*(wf.n3+2)*(wf.n4+2);
  MPI_Status status;
  MPI_Comm_rank (MPI_COMM_WORLD, &myid);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);

  vector<complex> u_from_left(num_comm), bet_from_left(num_comm),u_to_right(num_comm), bet_to_right(num_comm);
  vector<complex> bet(num_comm);
  vector<complex> gam_from_right(num_comm), u_from_right(num_comm), gam_to_left(num_comm), u_to_left(num_comm);
  vector<complex> wf_tmp(num_data);
  vector<complex> gam(num_data);
  if ( myid > 0 )
    left=myid-1;
  else 
    left=MPI_PROC_NULL;
  if ( myid < numprocs-1  )
    right=myid+1;
  else
    right=MPI_PROC_NULL;  //// define MPI_NULL buffer

  if (myid==0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int l = 1 ; l <= wf.n4 ; l++ )
     for( int j = 1 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1 ; i++ )
   	  {
	   bet[wf.in3(l,j,i)]=tridag_mid[wf.in4(l,1,j,i)];
      	   wf_tmp[wf.in4(l,1,j,i) ] = wf_rightside[wf.in4(l,1,j,i)]/bet[wf.in3(l,j,i)];
	  }}}
  else
    {
      MPI_Recv(&u_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag1,MPI_COMM_WORLD, &status);
      MPI_Recv(&bet_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag2,MPI_COMM_WORLD, &status); 
   #pragma omp parallel for  num_threads(nthreads)  collapse(3)
     for( int l = 1 ; l <= wf.n4 ; l++ )
     for( int j = 1 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1 ; i++ )
   	  {
   	   gam[wf.in4(l,1,j,i)]=tridag_upp/bet_from_left[wf.in3(l,j,i)];
           bet[wf.in3(l,j,i)]  = tridag_mid[wf.in4(l,1,j,i) ];
           bet[wf.in3(l,j,i)] -= tridag_low*gam[wf.in4(l,1,j,i)];
           wf_tmp[wf.in4(l,1,j,i)] = ( wf_rightside[wf.in4(l,1,j,i) ]-tridag_low*u_from_left[wf.in3(l,j,i)] )/bet[wf.in3(l,j,i)];
     }}     
    }
 #pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int l = 1 ; l <= wf.n4 ; l++ )
     for( int j = 1 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1 ; i++ )
   	  {
   	for( int k = 2 ; k <= wf.n3 ; k++ )
   	  {
      	   gam[wf.in4(l,k,j,i)] = tridag_upp/bet[wf.in3(l,j,i)];
           bet[wf.in3(l,j,i)]  = tridag_mid[wf.in4(l,k,j,i)];
      	   bet[wf.in3(l,j,i)] -= tridag_low*gam[wf.in4(l,k,j,i)];
           wf_tmp[wf.in4(l,k,j,i)] = ( wf_rightside[wf.in4(l,k,j,i)]-tridag_low*wf_tmp[wf.in4(l,k-1,j,i)] )/bet[wf.in3(l,j,i)];
    }}}
  

  if (myid < numprocs-1 )
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int l = 1 ; l <= wf.n4 ; l++ )
     for( int j = 1 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1 ; i++ )
   	  {
           u_to_right[wf.in3(l,j,i)]=wf_tmp[wf.in4(l,wf.n3,j,i)];
           bet_to_right[wf.in3(l,j,i)]=bet[wf.in3(l,j,i)];
          }
       }
      MPI_Send(&u_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag1,MPI_COMM_WORLD);
      MPI_Send(&bet_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag2,MPI_COMM_WORLD);
    }


  if(myid != numprocs-1)
    {
      MPI_Recv(&gam_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag3,MPI_COMM_WORLD, &status);
      MPI_Recv(&u_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag4,MPI_COMM_WORLD, &status);
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int l = 1 ; l <= wf.n4 ; l++ )
     for( int j = 1 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1 ; i++ )
   	  {
           wf_tmp[wf.in4(l,wf.n3,j,i)] -= gam_from_right[wf.in3(l,j,i)]*u_from_right[wf.in3(l,j,i)];    
          }}
    }
  #pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int l = 1 ; l <= wf.n4 ; l++ )
     for( int j = 1 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1 ; i++ )
   	  {
           for( int k = wf.n3-1 ; k >= 1 ; k-- )
             {
              wf_tmp[ wf.in4(l,k,j,i)] -= gam[ wf.in4(l,k+1,j,i) ]*wf_tmp[wf.in4(l,k+1,j,i)];
     }}}

  if (myid != 0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int l = 1 ; l <= wf.n4 ; l++ )
     for( int j = 1 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1 ; i++ )
   	  {
            gam_to_left[wf.in3(l,j,i)] = gam[wf.in4(l,1,j,i)];
      	    u_to_left[wf.in3(l,j,i)] = wf_tmp[wf.in4(l,1,j,i)];
          }
       }
      MPI_Send(&gam_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag3,MPI_COMM_WORLD);
      MPI_Send(&u_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag4,MPI_COMM_WORLD);
    }
#pragma omp parallel for  num_threads(nthreads) collapse(4)
   for( int l = 1 ; l <= wf.n4 ; l++ )
   for( int k = 1 ; k <= wf.n3 ; k++ )
     {
     for( int j = 1 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1 ; i++ )
   	  {
           wf.wave[ wf.in4( l, k , j , i ) ] = wf_tmp[wf.in4(l,k,j,i)];
   	  }}} 
}// end of tridag_X3_MPI


inline void Tridag_X3_MPI_opti(wavefunction &wf,  const complex &tridag_low, const vector<complex> &tridag_mid, const complex &tridag_upp, const vector<complex> &wf_rightside, const int nthreads )
{   

/*  if ( tridag_middle_diagonal[ 1 ] == 0. )
    {
      nerror( "Error: tridag_middle_diagonal[1]==0 in tridag!" );
    }
*/

  int tag1_1=11, tag2_1=21, tag3_1=31, tag4_1=41;
  int tag1_2=12, tag2_2=22, tag3_2=32, tag4_2=42;
  int tag1_3=13, tag2_3=23, tag3_3=33, tag4_3=43;
  int tag1_4=14, tag2_4=24, tag3_4=34, tag4_4=44;
  int tag1_5=15, tag2_5=25, tag3_5=35, tag4_5=45;
  int tag1_6=16, tag2_6=26, tag3_6=36, tag4_6=46;
  int tag1_7=17, tag2_7=27, tag3_7=37, tag4_7=47;
  int tag1_8=18, tag2_8=28, tag3_8=38, tag4_8=48;
  int left, right, myid, numprocs;
  int num_3D=(wf.n1+2)*(wf.n2+2)*(wf.n4+2);
  int num_4D=(wf.n1+2)*(wf.n2+2)*(wf.n3+2)*(wf.n4+2);
  int num_comm=wf.n1*wf.n2*wf.n4/8;
  int index;
  MPI_Status status;
  MPI_Comm_rank (MPI_COMM_WORLD, &myid);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);

  vector<complex> u_from_left(num_comm), bet_from_left(num_comm),u_to_right(num_comm), bet_to_right(num_comm);
  vector<complex> gam_from_right(num_comm), u_from_right(num_comm), gam_to_left(num_comm), u_to_left(num_comm);
  vector<complex> bet(num_3D);
  vector<complex> wf_tmp(num_4D);
  vector<complex> gam(num_4D);
  if ( myid > 0 )
    left=myid-1;
  else 
    left=MPI_PROC_NULL;
  if ( myid < numprocs-1  )
    right=myid+1;
  else
    right=MPI_PROC_NULL;  //// define MPI_NULL buffer

//----------------------area 1 in Low part
  if (myid==0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int l = 1 ; l <= wf.n4/2 ; l++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
	   bet[wf.in3(l,j,i)]=tridag_mid[wf.in4(l,1,j,i)];
      	   wf_tmp[wf.in4(l,1,j,i) ] = wf_rightside[wf.in4(l,1,j,i)]/bet[wf.in3(l,j,i)];
	  }}}
  else
    {
      MPI_Recv(&u_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag1_1,MPI_COMM_WORLD, &status);
      MPI_Recv(&bet_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag2_1,MPI_COMM_WORLD, &status); 
     #pragma omp parallel for  num_threads(nthreads) collapse(3) 
     for( int l = 1 ; l <= wf.n4/2 ; l++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
	   index=(l-1)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+i-1;
   	   gam[wf.in4(l,1,j,i)]=tridag_upp/bet_from_left[index];
           bet[wf.in3(l,j,i)]  = tridag_mid[wf.in4(l,1,j,i) ];
           bet[wf.in3(l,j,i)] -= tridag_low*gam[wf.in4(l,1,j,i)];
           wf_tmp[wf.in4(l,1,j,i)] = ( wf_rightside[wf.in4(l,1,j,i) ]-tridag_low*u_from_left[index] )/bet[wf.in3(l,j,i)];
     }}     
    }
 #pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int l = 1 ; l <= wf.n4/2 ; l++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
   	for( int k = 2 ; k <= wf.n3 ; k++ )
   	  {
      	   gam[wf.in4(l,k,j,i)] = tridag_upp/bet[wf.in3(l,j,i)];
           bet[wf.in3(l,j,i)]  = tridag_mid[wf.in4(l,k,j,i)];
      	   bet[wf.in3(l,j,i)] -= tridag_low*gam[wf.in4(l,k,j,i)];
           wf_tmp[wf.in4(l,k,j,i)] = ( wf_rightside[wf.in4(l,k,j,i)]-tridag_low*wf_tmp[wf.in4(l,k-1,j,i)] )/bet[wf.in3(l,j,i)];
    }}}
  

  if (myid < numprocs-1 )
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int l = 1 ; l <= wf.n4/2 ; l++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
	   index=(l-1)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+i-1;
           u_to_right[index]=wf_tmp[wf.in4(l,wf.n3,j,i)];
           bet_to_right[index]=bet[wf.in3(l,j,i)];
          }
       }
      MPI_Send(&u_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag1_1,MPI_COMM_WORLD);
      MPI_Send(&bet_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag2_1,MPI_COMM_WORLD);
    }

//----------------------area 2 in Low part
  if (myid==0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int l = 1 ; l <= wf.n4/2 ; l++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
	   bet[wf.in3(l,j,i)]=tridag_mid[wf.in4(l,1,j,i)];
      	   wf_tmp[wf.in4(l,1,j,i) ] = wf_rightside[wf.in4(l,1,j,i)]/bet[wf.in3(l,j,i)];
	  }}}
  else
    {
      MPI_Recv(&u_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag1_2,MPI_COMM_WORLD, &status);
      MPI_Recv(&bet_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag2_2,MPI_COMM_WORLD, &status); 
     #pragma omp parallel for  num_threads(nthreads)  collapse(3)
     for( int l = 1 ; l <= wf.n4/2 ; l++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
	   index=(l-1)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+(i-wf.n1/2)-1;
   	   gam[wf.in4(l,1,j,i)]=tridag_upp/bet_from_left[index];
           bet[wf.in3(l,j,i)]  = tridag_mid[wf.in4(l,1,j,i) ];
           bet[wf.in3(l,j,i)] -= tridag_low*gam[wf.in4(l,1,j,i)];
           wf_tmp[wf.in4(l,1,j,i)] = ( wf_rightside[wf.in4(l,1,j,i) ]-tridag_low*u_from_left[index] )/bet[wf.in3(l,j,i)];
     }}     
    }
 #pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int l = 1 ; l <= wf.n4/2 ; l++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
   	for( int k = 2 ; k <= wf.n3 ; k++ )
   	  {
      	   gam[wf.in4(l,k,j,i)] = tridag_upp/bet[wf.in3(l,j,i)];
           bet[wf.in3(l,j,i)]  = tridag_mid[wf.in4(l,k,j,i)];
      	   bet[wf.in3(l,j,i)] -= tridag_low*gam[wf.in4(l,k,j,i)];
           wf_tmp[wf.in4(l,k,j,i)] = ( wf_rightside[wf.in4(l,k,j,i)]-tridag_low*wf_tmp[wf.in4(l,k-1,j,i)] )/bet[wf.in3(l,j,i)];
    }}}
  

  if (myid < numprocs-1 )
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int l = 1 ; l <= wf.n4/2 ; l++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
	   index=(l-1)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+(i-wf.n1/2)-1;
           u_to_right[index]=wf_tmp[wf.in4(l,wf.n3,j,i)];
           bet_to_right[index]=bet[wf.in3(l,j,i)];
          }
       }
      MPI_Send(&u_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag1_2,MPI_COMM_WORLD);
      MPI_Send(&bet_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag2_2,MPI_COMM_WORLD);
    }

//----------------------area 3 in Low part
  if (myid==0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int l = 1 ; l <= wf.n4/2 ; l++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
	   bet[wf.in3(l,j,i)]=tridag_mid[wf.in4(l,1,j,i)];
      	   wf_tmp[wf.in4(l,1,j,i) ] = wf_rightside[wf.in4(l,1,j,i)]/bet[wf.in3(l,j,i)];
	  }}}
  else
    {
      MPI_Recv(&u_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag1_3,MPI_COMM_WORLD, &status);
      MPI_Recv(&bet_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag2_3,MPI_COMM_WORLD, &status); 
     #pragma omp parallel for  num_threads(nthreads)  collapse(3)
     for( int l = 1 ; l <= wf.n4/2 ; l++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
	   index=(l-1)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-1;
   	   gam[wf.in4(l,1,j,i)]=tridag_upp/bet_from_left[index];
           bet[wf.in3(l,j,i)]  = tridag_mid[wf.in4(l,1,j,i) ];
           bet[wf.in3(l,j,i)] -= tridag_low*gam[wf.in4(l,1,j,i)];
           wf_tmp[wf.in4(l,1,j,i)] = ( wf_rightside[wf.in4(l,1,j,i) ]-tridag_low*u_from_left[index] )/bet[wf.in3(l,j,i)];
     }}     
    }
 #pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int l = 1 ; l <= wf.n4/2 ; l++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
   	for( int k = 2 ; k <= wf.n3 ; k++ )
   	  {
      	   gam[wf.in4(l,k,j,i)] = tridag_upp/bet[wf.in3(l,j,i)];
           bet[wf.in3(l,j,i)]  = tridag_mid[wf.in4(l,k,j,i)];
      	   bet[wf.in3(l,j,i)] -= tridag_low*gam[wf.in4(l,k,j,i)];
           wf_tmp[wf.in4(l,k,j,i)] = ( wf_rightside[wf.in4(l,k,j,i)]-tridag_low*wf_tmp[wf.in4(l,k-1,j,i)] )/bet[wf.in3(l,j,i)];
    }}}
  

  if (myid < numprocs-1 )
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int l = 1 ; l <= wf.n4/2 ; l++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
	   index=(l-1)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-1;
           u_to_right[index]=wf_tmp[wf.in4(l,wf.n3,j,i)];
           bet_to_right[index]=bet[wf.in3(l,j,i)];
          }
       }
      MPI_Send(&u_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag1_3,MPI_COMM_WORLD);
      MPI_Send(&bet_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag2_3,MPI_COMM_WORLD);
    }

//----------------------area 4 in Low part
  if (myid==0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int l = 1 ; l <= wf.n4/2 ; l++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
	   bet[wf.in3(l,j,i)]=tridag_mid[wf.in4(l,1,j,i)];
      	   wf_tmp[wf.in4(l,1,j,i) ] = wf_rightside[wf.in4(l,1,j,i)]/bet[wf.in3(l,j,i)];
	  }}}
  else
    {
      MPI_Recv(&u_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag1_4,MPI_COMM_WORLD, &status);
      MPI_Recv(&bet_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag2_4,MPI_COMM_WORLD, &status); 
     #pragma omp parallel for  num_threads(nthreads)  collapse(3)
	 for( int l = 1 ; l <= wf.n4/2 ; l++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
	   index=(l-1)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-wf.n1/2-1;
   	   gam[wf.in4(l,1,j,i)]=tridag_upp/bet_from_left[index];
       bet[wf.in3(l,j,i)]  = tridag_mid[wf.in4(l,1,j,i) ];
       bet[wf.in3(l,j,i)] -= tridag_low*gam[wf.in4(l,1,j,i)];
       wf_tmp[wf.in4(l,1,j,i)] = ( wf_rightside[wf.in4(l,1,j,i) ]-tridag_low*u_from_left[index] )/bet[wf.in3(l,j,i)];
     }}     
    }
 #pragma omp parallel for  num_threads(nthreads) collapse(3)
   for( int l = 1 ; l <= wf.n4/2 ; l++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
   	for( int k = 2 ; k <= wf.n3 ; k++ )
   	  {
      	   gam[wf.in4(l,k,j,i)] = tridag_upp/bet[wf.in3(l,j,i)];
           bet[wf.in3(l,j,i)]  = tridag_mid[wf.in4(l,k,j,i)];
      	   bet[wf.in3(l,j,i)] -= tridag_low*gam[wf.in4(l,k,j,i)];
           wf_tmp[wf.in4(l,k,j,i)] = ( wf_rightside[wf.in4(l,k,j,i)]-tridag_low*wf_tmp[wf.in4(l,k-1,j,i)] )/bet[wf.in3(l,j,i)];
    }}}
  

  if (myid < numprocs-1 )
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int l = 1 ; l <= wf.n4/2 ; l++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
	   index=(l-1)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-wf.n1/2-1;
           u_to_right[index]=wf_tmp[wf.in4(l,wf.n3,j,i)];
           bet_to_right[index]=bet[wf.in3(l,j,i)];
          }
       }
      MPI_Send(&u_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag1_4,MPI_COMM_WORLD);
      MPI_Send(&bet_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag2_4,MPI_COMM_WORLD);
    }

//----------------------area 5 in Low part
  if (myid==0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int l = 1+wf.n4/2 ; l <= wf.n4 ; l++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
	   bet[wf.in3(l,j,i)]=tridag_mid[wf.in4(l,1,j,i)];
      	   wf_tmp[wf.in4(l,1,j,i) ] = wf_rightside[wf.in4(l,1,j,i)]/bet[wf.in3(l,j,i)];
	  }}}
  else
    {
      MPI_Recv(&u_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag1_5,MPI_COMM_WORLD, &status);
      MPI_Recv(&bet_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag2_5,MPI_COMM_WORLD, &status); 
     #pragma omp parallel for  num_threads(nthreads) collapse(3) 
     for( int l = 1+wf.n4/2 ; l <= wf.n4 ; l++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
	   index=(l-1-wf.n4/2)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+i-1;
   	   gam[wf.in4(l,1,j,i)]=tridag_upp/bet_from_left[index];
           bet[wf.in3(l,j,i)]  = tridag_mid[wf.in4(l,1,j,i) ];
           bet[wf.in3(l,j,i)] -= tridag_low*gam[wf.in4(l,1,j,i)];
           wf_tmp[wf.in4(l,1,j,i)] = ( wf_rightside[wf.in4(l,1,j,i) ]-tridag_low*u_from_left[index] )/bet[wf.in3(l,j,i)];
     }}     
    }
 #pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int l = 1+wf.n4/2 ; l <= wf.n4 ; l++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
   	for( int k = 2 ; k <= wf.n3 ; k++ )
   	  {
      	   gam[wf.in4(l,k,j,i)] = tridag_upp/bet[wf.in3(l,j,i)];
           bet[wf.in3(l,j,i)]  = tridag_mid[wf.in4(l,k,j,i)];
      	   bet[wf.in3(l,j,i)] -= tridag_low*gam[wf.in4(l,k,j,i)];
           wf_tmp[wf.in4(l,k,j,i)] = ( wf_rightside[wf.in4(l,k,j,i)]-tridag_low*wf_tmp[wf.in4(l,k-1,j,i)] )/bet[wf.in3(l,j,i)];
    }}}
  

  if (myid < numprocs-1 )
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int l = 1+wf.n4/2 ; l <= wf.n4 ; l++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
	   index=(l-1-wf.n4/2)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+i-1;
           u_to_right[index]=wf_tmp[wf.in4(l,wf.n3,j,i)];
           bet_to_right[index]=bet[wf.in3(l,j,i)];
          }
       }
      MPI_Send(&u_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag1_5,MPI_COMM_WORLD);
      MPI_Send(&bet_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag2_5,MPI_COMM_WORLD);
    }

//----------------------area 6 in Low part
  if (myid==0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int l = 1+wf.n4/2 ; l <= wf.n4 ; l++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
	   bet[wf.in3(l,j,i)]=tridag_mid[wf.in4(l,1,j,i)];
      	   wf_tmp[wf.in4(l,1,j,i) ] = wf_rightside[wf.in4(l,1,j,i)]/bet[wf.in3(l,j,i)];
	  }}}
  else
    {
      MPI_Recv(&u_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag1_6,MPI_COMM_WORLD, &status);
      MPI_Recv(&bet_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag2_6,MPI_COMM_WORLD, &status); 
     #pragma omp parallel for  num_threads(nthreads)  collapse(3)
     for( int l = 1+wf.n4/2 ; l <= wf.n4 ; l++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
	   index=(l-1-wf.n4/2)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+(i-wf.n1/2)-1;
   	   gam[wf.in4(l,1,j,i)]=tridag_upp/bet_from_left[index];
           bet[wf.in3(l,j,i)]  = tridag_mid[wf.in4(l,1,j,i) ];
           bet[wf.in3(l,j,i)] -= tridag_low*gam[wf.in4(l,1,j,i)];
           wf_tmp[wf.in4(l,1,j,i)] = ( wf_rightside[wf.in4(l,1,j,i) ]-tridag_low*u_from_left[index] )/bet[wf.in3(l,j,i)];
     }}     
    }
 #pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int l = 1+wf.n4/2 ; l <= wf.n4 ; l++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
   	for( int k = 2 ; k <= wf.n3 ; k++ )
   	  {
      	   gam[wf.in4(l,k,j,i)] = tridag_upp/bet[wf.in3(l,j,i)];
           bet[wf.in3(l,j,i)]  = tridag_mid[wf.in4(l,k,j,i)];
      	   bet[wf.in3(l,j,i)] -= tridag_low*gam[wf.in4(l,k,j,i)];
           wf_tmp[wf.in4(l,k,j,i)] = ( wf_rightside[wf.in4(l,k,j,i)]-tridag_low*wf_tmp[wf.in4(l,k-1,j,i)] )/bet[wf.in3(l,j,i)];
    }}}
  

  if (myid < numprocs-1 )
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int l = 1+wf.n4/2 ; l <= wf.n4 ; l++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
	   index=(l-1-wf.n4/2)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+(i-wf.n1/2)-1;
           u_to_right[index]=wf_tmp[wf.in4(l,wf.n3,j,i)];
           bet_to_right[index]=bet[wf.in3(l,j,i)];
          }
       }
      MPI_Send(&u_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag1_6,MPI_COMM_WORLD);
      MPI_Send(&bet_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag2_6,MPI_COMM_WORLD);
    }

//----------------------area 7 in Low part
  if (myid==0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int l = 1+wf.n4/2 ; l <= wf.n4 ; l++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
	   bet[wf.in3(l,j,i)]=tridag_mid[wf.in4(l,1,j,i)];
      	   wf_tmp[wf.in4(l,1,j,i) ] = wf_rightside[wf.in4(l,1,j,i)]/bet[wf.in3(l,j,i)];
	  }}}
  else
    {
      MPI_Recv(&u_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag1_7,MPI_COMM_WORLD, &status);
      MPI_Recv(&bet_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag2_7,MPI_COMM_WORLD, &status); 
     #pragma omp parallel for  num_threads(nthreads)  collapse(3)
     for( int l = 1+wf.n4/2 ; l <= wf.n4 ; l++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
	   index=(l-1-wf.n4/2)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-1;
   	   gam[wf.in4(l,1,j,i)]=tridag_upp/bet_from_left[index];
           bet[wf.in3(l,j,i)]  = tridag_mid[wf.in4(l,1,j,i) ];
           bet[wf.in3(l,j,i)] -= tridag_low*gam[wf.in4(l,1,j,i)];
           wf_tmp[wf.in4(l,1,j,i)] = ( wf_rightside[wf.in4(l,1,j,i) ]-tridag_low*u_from_left[index] )/bet[wf.in3(l,j,i)];
     }}     
    }
 #pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int l = 1+wf.n4/2 ; l <= wf.n4 ; l++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
   	for( int k = 2 ; k <= wf.n3 ; k++ )
   	  {
      	   gam[wf.in4(l,k,j,i)] = tridag_upp/bet[wf.in3(l,j,i)];
           bet[wf.in3(l,j,i)]  = tridag_mid[wf.in4(l,k,j,i)];
      	   bet[wf.in3(l,j,i)] -= tridag_low*gam[wf.in4(l,k,j,i)];
           wf_tmp[wf.in4(l,k,j,i)] = ( wf_rightside[wf.in4(l,k,j,i)]-tridag_low*wf_tmp[wf.in4(l,k-1,j,i)] )/bet[wf.in3(l,j,i)];
    }}}
  

  if (myid < numprocs-1 )
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int l = 1+wf.n4/2 ; l <= wf.n4 ; l++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
	   index=(l-1-wf.n4/2)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-1;
           u_to_right[index]=wf_tmp[wf.in4(l,wf.n3,j,i)];
           bet_to_right[index]=bet[wf.in3(l,j,i)];
          }
       }
      MPI_Send(&u_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag1_7,MPI_COMM_WORLD);
      MPI_Send(&bet_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag2_7,MPI_COMM_WORLD);
    }

//----------------------area 8 in Low part
  if (myid==0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int l = 1+wf.n4/2 ; l <= wf.n4 ; l++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
	   bet[wf.in3(l,j,i)]=tridag_mid[wf.in4(l,1,j,i)];
      	   wf_tmp[wf.in4(l,1,j,i) ] = wf_rightside[wf.in4(l,1,j,i)]/bet[wf.in3(l,j,i)];
	  }}}
  else
    {
      MPI_Recv(&u_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag1_8,MPI_COMM_WORLD, &status);
      MPI_Recv(&bet_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag2_8,MPI_COMM_WORLD, &status); 
     #pragma omp parallel for  num_threads(nthreads)  collapse(3)
	 for( int l = 1+wf.n4/2 ; l <= wf.n4 ; l++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
	   index=(l-1-wf.n4/2)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-wf.n1/2-1;
   	   gam[wf.in4(l,1,j,i)]=tridag_upp/bet_from_left[index];
       bet[wf.in3(l,j,i)]  = tridag_mid[wf.in4(l,1,j,i) ];
       bet[wf.in3(l,j,i)] -= tridag_low*gam[wf.in4(l,1,j,i)];
       wf_tmp[wf.in4(l,1,j,i)] = ( wf_rightside[wf.in4(l,1,j,i) ]-tridag_low*u_from_left[index] )/bet[wf.in3(l,j,i)];
     }}     
    }
 #pragma omp parallel for  num_threads(nthreads) collapse(3)
   for( int l = 1+wf.n4/2 ; l <= wf.n4 ; l++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
   	for( int k = 2 ; k <= wf.n3 ; k++ )
   	  {
      	   gam[wf.in4(l,k,j,i)] = tridag_upp/bet[wf.in3(l,j,i)];
           bet[wf.in3(l,j,i)]  = tridag_mid[wf.in4(l,k,j,i)];
      	   bet[wf.in3(l,j,i)] -= tridag_low*gam[wf.in4(l,k,j,i)];
           wf_tmp[wf.in4(l,k,j,i)] = ( wf_rightside[wf.in4(l,k,j,i)]-tridag_low*wf_tmp[wf.in4(l,k-1,j,i)] )/bet[wf.in3(l,j,i)];
    }}}
  

  if (myid < numprocs-1 )
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int l = 1+wf.n4/2 ; l <= wf.n4 ; l++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
	   index=(l-1-wf.n4/2)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-wf.n1/2-1;
           u_to_right[index]=wf_tmp[wf.in4(l,wf.n3,j,i)];
           bet_to_right[index]=bet[wf.in3(l,j,i)];
          }
       }
      MPI_Send(&u_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag1_8,MPI_COMM_WORLD);
      MPI_Send(&bet_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag2_8,MPI_COMM_WORLD);
    }



//----------------------area 1 in up part
  if(myid != numprocs-1)
    {
      MPI_Recv(&gam_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag3_1,MPI_COMM_WORLD, &status);
      MPI_Recv(&u_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag4_1,MPI_COMM_WORLD, &status);
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	for( int l = 1 ; l <= wf.n4/2 ; l++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
           index=(l-1)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+i-1;
           wf_tmp[wf.in4(l,wf.n3,j,i)] -= gam_from_right[index]*u_from_right[index]; 
          }}
    }
  #pragma omp parallel for  num_threads(nthreads) collapse(3)
  	for( int l = 1 ; l <= wf.n4/2 ; l++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
           for( int k = wf.n3-1 ; k >= 1 ; k-- )
             {
              wf_tmp[ wf.in4(l,k,j,i)] -= gam[ wf.in4(l,k+1,j,i) ]*wf_tmp[wf.in4(l,k+1,j,i)];
     }}}

  if (myid != 0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int l = 1 ; l <= wf.n4/2 ; l++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
            index=(l-1)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+i-1;
            gam_to_left[index] = gam[wf.in4(l,1,j,i)];
      	    u_to_left[index] = wf_tmp[wf.in4(l,1,j,i)];
          }
       }
      MPI_Send(&gam_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag3_1,MPI_COMM_WORLD);
      MPI_Send(&u_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag4_1,MPI_COMM_WORLD);
    }

//----------------------area 2 in up part
  if(myid != numprocs-1)
    {
      MPI_Recv(&gam_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag3_2,MPI_COMM_WORLD, &status);
      MPI_Recv(&u_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag4_2,MPI_COMM_WORLD, &status);
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int l = 1 ; l <= wf.n4/2 ; l++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
           index=(l-1)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+i-wf.n1/2-1;
           wf_tmp[wf.in4(l,wf.n3,j,i)] -= gam_from_right[index]*u_from_right[index]; 
          }}
    }
  #pragma omp parallel for  num_threads(nthreads) collapse(3)
  	for( int l = 1 ; l <= wf.n4/2 ; l++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
           for( int k = wf.n3-1 ; k >= 1 ; k-- )
             {
              wf_tmp[ wf.in4(l,k,j,i)] -= gam[ wf.in4(l,k+1,j,i) ]*wf_tmp[wf.in4(l,k+1,j,i)];
     }}}

  if (myid != 0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int l = 1 ; l <= wf.n4/2 ; l++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
            index=(l-1)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+i-wf.n1/2-1;
            gam_to_left[index] = gam[wf.in4(l,1,j,i)];
      	    u_to_left[index] = wf_tmp[wf.in4(l,1,j,i)];
          }
       }
      MPI_Send(&gam_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag3_2,MPI_COMM_WORLD);
      MPI_Send(&u_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag4_2,MPI_COMM_WORLD);
    }

//----------------------area 3 in up part
  if(myid != numprocs-1)
    {
      MPI_Recv(&gam_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag3_3,MPI_COMM_WORLD, &status);
      MPI_Recv(&u_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag4_3,MPI_COMM_WORLD, &status);
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int l = 1 ; l <= wf.n4/2 ; l++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
           index=(l-1)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-1;
           wf_tmp[wf.in4(l,wf.n3,j,i)] -= gam_from_right[index]*u_from_right[index]; 
          }}
    }
  #pragma omp parallel for  num_threads(nthreads) collapse(3)
  	for( int l = 1 ; l <= wf.n4/2 ; l++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
           for( int k = wf.n3-1 ; k >= 1 ; k-- )
             {
              wf_tmp[ wf.in4(l,k,j,i)] -= gam[ wf.in4(l,k+1,j,i) ]*wf_tmp[wf.in4(l,k+1,j,i)];
     }}}

  if (myid != 0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int l = 1 ; l <= wf.n4/2 ; l++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
            index=(l-1)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-1;
            gam_to_left[index] = gam[wf.in4(l,1,j,i)];
      	    u_to_left[index] = wf_tmp[wf.in4(l,1,j,i)];
          }
       }
      MPI_Send(&gam_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag3_3,MPI_COMM_WORLD);
      MPI_Send(&u_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag4_3,MPI_COMM_WORLD);
    }

//----------------------area 4 in up part
  if(myid != numprocs-1)
    {
      MPI_Recv(&gam_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag3_4,MPI_COMM_WORLD, &status);
      MPI_Recv(&u_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag4_4,MPI_COMM_WORLD, &status);
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int l = 1 ; l <= wf.n4/2 ; l++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
           index=(l-1)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-wf.n1/2-1;
           wf_tmp[wf.in4(l,wf.n3,j,i)] -= gam_from_right[index]*u_from_right[index]; 
          }}
    }
  #pragma omp parallel for  num_threads(nthreads) collapse(3)
  	for( int l = 1 ; l <= wf.n4/2 ; l++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
           for( int k = wf.n3-1 ; k >= 1 ; k-- )
             {
              wf_tmp[ wf.in4(l,k,j,i)] -= gam[ wf.in4(l,k+1,j,i) ]*wf_tmp[wf.in4(l,k+1,j,i)];
     }}}

  if (myid != 0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int l = 1 ; l <= wf.n4/2 ; l++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
            index=(l-1)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-wf.n1/2-1;
            gam_to_left[index] = gam[wf.in4(l,1,j,i)];
      	    u_to_left[index] = wf_tmp[wf.in4(l,1,j,i)];
          }
       }
      MPI_Send(&gam_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag3_4,MPI_COMM_WORLD);
      MPI_Send(&u_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag4_4,MPI_COMM_WORLD);
    }


//----------------------area 5 in up part
  if(myid != numprocs-1)
    {
      MPI_Recv(&gam_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag3_5,MPI_COMM_WORLD, &status);
      MPI_Recv(&u_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag4_5,MPI_COMM_WORLD, &status);
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	for( int l = 1+wf.n4/2 ; l <= wf.n4 ; l++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
           index=(l-1-wf.n4/2)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+i-1;
           wf_tmp[wf.in4(l,wf.n3,j,i)] -= gam_from_right[index]*u_from_right[index]; 
          }}
    }
  #pragma omp parallel for  num_threads(nthreads) collapse(3)
  	for( int l = 1+wf.n4/2 ; l <= wf.n4 ; l++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
           for( int k = wf.n3-1 ; k >= 1 ; k-- )
             {
              wf_tmp[ wf.in4(l,k,j,i)] -= gam[ wf.in4(l,k+1,j,i) ]*wf_tmp[wf.in4(l,k+1,j,i)];
     }}}

  if (myid != 0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int l = 1+wf.n4/2 ; l <= wf.n4 ; l++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
            index=(l-1-wf.n4/2)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+i-1;
            gam_to_left[index] = gam[wf.in4(l,1,j,i)];
      	    u_to_left[index] = wf_tmp[wf.in4(l,1,j,i)];
          }
       }
      MPI_Send(&gam_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag3_5,MPI_COMM_WORLD);
      MPI_Send(&u_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag4_5,MPI_COMM_WORLD);
    }

//----------------------area 6 in up part
  if(myid != numprocs-1)
    {
      MPI_Recv(&gam_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag3_6,MPI_COMM_WORLD, &status);
      MPI_Recv(&u_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag4_6,MPI_COMM_WORLD, &status);
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int l = 1+wf.n4/2 ; l <= wf.n4 ; l++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
           index=(l-1-wf.n4/2)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+i-wf.n1/2-1;
           wf_tmp[wf.in4(l,wf.n3,j,i)] -= gam_from_right[index]*u_from_right[index]; 
          }}
    }
  #pragma omp parallel for  num_threads(nthreads) collapse(3)
  	for( int l = 1+wf.n4/2 ; l <= wf.n4 ; l++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
           for( int k = wf.n3-1 ; k >= 1 ; k-- )
             {
              wf_tmp[ wf.in4(l,k,j,i)] -= gam[ wf.in4(l,k+1,j,i) ]*wf_tmp[wf.in4(l,k+1,j,i)];
     }}}

  if (myid != 0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int l = 1+wf.n4/2 ; l <= wf.n4 ; l++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
            index=(l-1-wf.n4/2)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+i-wf.n1/2-1;
            gam_to_left[index] = gam[wf.in4(l,1,j,i)];
      	    u_to_left[index] = wf_tmp[wf.in4(l,1,j,i)];
          }
       }
      MPI_Send(&gam_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag3_6,MPI_COMM_WORLD);
      MPI_Send(&u_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag4_6,MPI_COMM_WORLD);
    }

//----------------------area 7 in up part
  if(myid != numprocs-1)
    {
      MPI_Recv(&gam_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag3_7,MPI_COMM_WORLD, &status);
      MPI_Recv(&u_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag4_7,MPI_COMM_WORLD, &status);
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int l = 1+wf.n4/2 ; l <= wf.n4 ; l++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
           index=(l-1-wf.n4/2)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-1;
           wf_tmp[wf.in4(l,wf.n3,j,i)] -= gam_from_right[index]*u_from_right[index]; 
          }}
    }
  #pragma omp parallel for  num_threads(nthreads) collapse(3)
  	for( int l = 1+wf.n4/2 ; l <= wf.n4 ; l++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
           for( int k = wf.n3-1 ; k >= 1 ; k-- )
             {
              wf_tmp[ wf.in4(l,k,j,i)] -= gam[ wf.in4(l,k+1,j,i) ]*wf_tmp[wf.in4(l,k+1,j,i)];
     }}}

  if (myid != 0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int l = 1+wf.n4/2 ; l <= wf.n4 ; l++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
            index=(l-1-wf.n4/2)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-1;
            gam_to_left[index] = gam[wf.in4(l,1,j,i)];
      	    u_to_left[index] = wf_tmp[wf.in4(l,1,j,i)];
          }
       }
      MPI_Send(&gam_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag3_7,MPI_COMM_WORLD);
      MPI_Send(&u_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag4_7,MPI_COMM_WORLD);
    }

//----------------------area 8 in up part
  if(myid != numprocs-1)
    {
      MPI_Recv(&gam_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag3_8,MPI_COMM_WORLD, &status);
      MPI_Recv(&u_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag4_8,MPI_COMM_WORLD, &status);
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int l = 1+wf.n4/2 ; l <= wf.n4 ; l++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
           index=(l-1-wf.n4/2)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-wf.n1/2-1;
           wf_tmp[wf.in4(l,wf.n3,j,i)] -= gam_from_right[index]*u_from_right[index]; 
          }}
    }
  #pragma omp parallel for  num_threads(nthreads) collapse(3)
  	for( int l = 1+wf.n4/2 ; l <= wf.n4 ; l++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
           for( int k = wf.n3-1 ; k >= 1 ; k-- )
             {
              wf_tmp[ wf.in4(l,k,j,i)] -= gam[ wf.in4(l,k+1,j,i) ]*wf_tmp[wf.in4(l,k+1,j,i)];
     }}}

  if (myid != 0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int l = 1+wf.n4/2 ; l <= wf.n4 ; l++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
            index=(l-1-wf.n4/2)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-wf.n1/2-1;
            gam_to_left[index] = gam[wf.in4(l,1,j,i)];
      	    u_to_left[index] = wf_tmp[wf.in4(l,1,j,i)];
          }
       }
      MPI_Send(&gam_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag3_8,MPI_COMM_WORLD);
      MPI_Send(&u_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag4_8,MPI_COMM_WORLD);
    }



//-------------------final step
#pragma omp parallel for  num_threads(nthreads) collapse(4)
   for( int l = 1 ; l <= wf.n4 ; l++ )
   for( int k = 1 ; k <= wf.n3 ; k++ )
     {
     for( int j = 1 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1 ; i++ )
   	  {
           wf.wave[ wf.in4(l, k , j , i ) ] = wf_tmp[wf.in4(l,k,j,i)];
   	  }}} 
}// end of tridag_MPI



inline void Tridag_X4_MPI(wavefunction &wf,  const complex &tridag_low, const vector<complex> &tridag_mid, const complex &tridag_upp, const vector<complex> &wf_rightside, const int nthreads )
{   

  int tag1=3, tag2=4, tag3=5, tag4=6;
  int left, right, myid, numprocs;
  int num_comm=(wf.n1+2)*(wf.n2+2)*(wf.n3+2);
  int num_data=(wf.n1+2)*(wf.n2+2)*(wf.n3+2)*(wf.n4+2);
  MPI_Status status;
  MPI_Comm_rank (MPI_COMM_WORLD, &myid);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);

  vector<complex> u_from_left(num_comm), bet_from_left(num_comm),u_to_right(num_comm), bet_to_right(num_comm);
  vector<complex> bet(num_comm);
  vector<complex> gam_from_right(num_comm), u_from_right(num_comm), gam_to_left(num_comm), u_to_left(num_comm);
  vector<complex> wf_tmp(num_data);
  vector<complex> gam(num_data);
  if ( myid > 0 )
    left=myid-1;
  else 
    left=MPI_PROC_NULL;
  if ( myid < numprocs-1  )
    right=myid+1;
  else
    right=MPI_PROC_NULL;  //// define MPI_NULL buffer

  if (myid==0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int k = 1 ; k <= wf.n3 ; k++ )
     for( int j = 1 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1 ; i++ )
   	  {
	   bet[wf.in3(k,j,i)]=tridag_mid[wf.in4(1,k,j,i)];
      	   wf_tmp[wf.in4(1,k,j,i) ] = wf_rightside[wf.in4(1,k,j,i)]/bet[wf.in3(k,j,i)];
	  }}}
  else
    {
      MPI_Recv(&u_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag1,MPI_COMM_WORLD, &status);
      MPI_Recv(&bet_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag2,MPI_COMM_WORLD, &status); 
   #pragma omp parallel for  num_threads(nthreads)  collapse(3)
     for( int k = 1 ; k <= wf.n3 ; k++ )
     for( int j = 1 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1 ; i++ )
   	  {
   	   gam[wf.in4(1,k,j,i)]=tridag_upp/bet_from_left[wf.in3(k,j,i)];
           bet[wf.in3(k,j,i)]  = tridag_mid[wf.in4(1,k,j,i) ];
           bet[wf.in3(k,j,i)] -= tridag_low*gam[wf.in4(1,k,j,i)];
           wf_tmp[wf.in4(1,k,j,i)] = ( wf_rightside[wf.in4(1,k,j,i) ]-tridag_low*u_from_left[wf.in3(k,j,i)] )/bet[wf.in3(k,j,i)];
     }}     
    }
 #pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int k = 1 ; k <= wf.n3 ; k++ )
     for( int j = 1 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1 ; i++ )
   	  {
   	for( int l = 2 ; l <= wf.n4 ; l++ )
   	  {
      	   gam[wf.in4(l,k,j,i)] = tridag_upp/bet[wf.in3(k,j,i)];
           bet[wf.in3(k,j,i)]  = tridag_mid[wf.in4(l,k,j,i)];
      	   bet[wf.in3(k,j,i)] -= tridag_low*gam[wf.in4(l,k,j,i)];
           wf_tmp[wf.in4(l,k,j,i)] = ( wf_rightside[wf.in4(l,k,j,i)]-tridag_low*wf_tmp[wf.in4(l-1,k,j,i)] )/bet[wf.in3(k,j,i)];
    }}}
  

  if (myid < numprocs-1 )
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int k = 1 ; k <= wf.n3 ; k++ )
     for( int j = 1 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1 ; i++ )
   	  {
           u_to_right[wf.in3(k,j,i)]=wf_tmp[wf.in4(wf.n4,k,j,i)];
           bet_to_right[wf.in3(k,j,i)]=bet[wf.in3(k,j,i)];
          }
       }
      MPI_Send(&u_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag1,MPI_COMM_WORLD);
      MPI_Send(&bet_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag2,MPI_COMM_WORLD);
    }


  if(myid != numprocs-1)
    {
      MPI_Recv(&gam_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag3,MPI_COMM_WORLD, &status);
      MPI_Recv(&u_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag4,MPI_COMM_WORLD, &status);
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int k = 1 ; k <= wf.n3 ; k++ )
     for( int j = 1 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1 ; i++ )
   	  {
           wf_tmp[wf.in4(wf.n4,k,j,i)] -= gam_from_right[wf.in3(k,j,i)]*u_from_right[wf.in3(k,j,i)];    
          }}
    }
  #pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int k = 1 ; k <= wf.n3 ; k++ )
     for( int j = 1 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1 ; i++ )
   	  {
           for( int l = wf.n4-1 ; l >= 1 ; l-- )
             {
              wf_tmp[ wf.in4(l,k,j,i)] -= gam[ wf.in4(l+1,k,j,i) ]*wf_tmp[wf.in4(l+1,k,j,i)];
     }}}

  if (myid != 0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int k = 1 ; k <= wf.n3 ; k++ )
     for( int j = 1 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1 ; i++ )
   	  {
            gam_to_left[wf.in3(k,j,i)] = gam[wf.in4(1,k,j,i)];
      	    u_to_left[wf.in3(k,j,i)] = wf_tmp[wf.in4(1,k,j,i)];
          }
       }
      MPI_Send(&gam_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag3,MPI_COMM_WORLD);
      MPI_Send(&u_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag4,MPI_COMM_WORLD);
    }
#pragma omp parallel for  num_threads(nthreads) collapse(4)
   for( int l = 1 ; l <= wf.n4 ; l++ )
   for( int k = 1 ; k <= wf.n3 ; k++ )
     {
     for( int j = 1 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1 ; i++ )
   	  {
           wf.wave[ wf.in4( l, k , j , i ) ] = wf_tmp[wf.in4(l,k,j,i)];
   	  }}} 
}// end of tridag_X4_MPI




inline void Tridag_X4_MPI_opti(wavefunction &wf,  const complex &tridag_low, const vector<complex> &tridag_mid, const complex &tridag_upp, const vector<complex> &wf_rightside, const int nthreads )
{   

/*  if ( tridag_middle_diagonal[ 1 ] == 0. )
    {
      nerror( "Error: tridag_middle_diagonal[1]==0 in tridag!" );
    }
*/

  if ( (wf.n1 % 2!=0)||(wf.n2 % 2!=0)||(wf.n3 % 2!=0)||(wf.n4 % 2!=0) )
    {
      cerr<<"Please use Tridag_X4_MPI in x4 propagation " <<endl;     
      exit(1);
    }

  int tag1_1=11, tag2_1=21, tag3_1=31, tag4_1=41;
  int tag1_2=12, tag2_2=22, tag3_2=32, tag4_2=42;
  int tag1_3=13, tag2_3=23, tag3_3=33, tag4_3=43;
  int tag1_4=14, tag2_4=24, tag3_4=34, tag4_4=44;
  int tag1_5=15, tag2_5=25, tag3_5=35, tag4_5=45;
  int tag1_6=16, tag2_6=26, tag3_6=36, tag4_6=46;
  int tag1_7=17, tag2_7=27, tag3_7=37, tag4_7=47;
  int tag1_8=18, tag2_8=28, tag3_8=38, tag4_8=48;
  int left, right, myid, numprocs;
  int num_3D=(wf.n1+2)*(wf.n2+2)*(wf.n3+2);
  int num_4D=(wf.n1+2)*(wf.n2+2)*(wf.n3+2)*(wf.n4+2);
  int num_comm=wf.n1*wf.n2*wf.n3/8;
  int index;
  MPI_Status status;
  MPI_Comm_rank (MPI_COMM_WORLD, &myid);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);

  vector<complex> u_from_left(num_comm), bet_from_left(num_comm),u_to_right(num_comm), bet_to_right(num_comm);
  vector<complex> gam_from_right(num_comm), u_from_right(num_comm), gam_to_left(num_comm), u_to_left(num_comm);
  vector<complex> bet(num_3D);
  vector<complex> wf_tmp(num_4D);
  vector<complex> gam(num_4D);
  if ( myid > 0 )
    left=myid-1;
  else 
    left=MPI_PROC_NULL;
  if ( myid < numprocs-1  )
    right=myid+1;
  else
    right=MPI_PROC_NULL;  //// define MPI_NULL buffer

//----------------------area 1 in Low part
  if (myid==0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int k = 1 ; k <= wf.n3/2 ; k++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
	   bet[wf.in3(k,j,i)]=tridag_mid[wf.in4(1,k,j,i)];
      	   wf_tmp[wf.in4(1,k,j,i) ] = wf_rightside[wf.in4(1,k,j,i)]/bet[wf.in3(k,j,i)];
	  }}}
  else
    {
      MPI_Recv(&u_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag1_1,MPI_COMM_WORLD, &status);
      MPI_Recv(&bet_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag2_1,MPI_COMM_WORLD, &status); 
     #pragma omp parallel for  num_threads(nthreads) collapse(3) 
     for( int k = 1 ; k <= wf.n3/2 ; k++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
	   index=(k-1)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+i-1;
   	   gam[wf.in4(1,k,j,i)]=tridag_upp/bet_from_left[index];
           bet[wf.in3(k,j,i)]  = tridag_mid[wf.in4(1,k,j,i) ];
           bet[wf.in3(k,j,i)] -= tridag_low*gam[wf.in4(1,k,j,i)];
           wf_tmp[wf.in4(1,k,j,i)] = ( wf_rightside[wf.in4(1,k,j,i) ]-tridag_low*u_from_left[index] )/bet[wf.in3(k,j,i)];
     }}     
    }
 #pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int k = 1 ; k <= wf.n3/2 ; k++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
   	for( int l = 2 ; l <= wf.n4 ; l++ )
   	  {
      	   gam[wf.in4(l,k,j,i)] = tridag_upp/bet[wf.in3(k,j,i)];
           bet[wf.in3(k,j,i)]  = tridag_mid[wf.in4(l,k,j,i)];
      	   bet[wf.in3(k,j,i)] -= tridag_low*gam[wf.in4(l,k,j,i)];
           wf_tmp[wf.in4(l,k,j,i)] = ( wf_rightside[wf.in4(l,k,j,i)]-tridag_low*wf_tmp[wf.in4(l-1,k,j,i)] )/bet[wf.in3(k,j,i)];
    }}}
  

  if (myid < numprocs-1 )
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int k = 1 ; k <= wf.n3/2 ; k++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
	   index=(k-1)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+i-1;
           u_to_right[index]=wf_tmp[wf.in4(wf.n4,k,j,i)];
           bet_to_right[index]=bet[wf.in3(k,j,i)];
          }
       }
      MPI_Send(&u_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag1_1,MPI_COMM_WORLD);
      MPI_Send(&bet_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag2_1,MPI_COMM_WORLD);
    }

//----------------------area 2 in Low part
  if (myid==0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int k = 1 ; k <= wf.n3/2 ; k++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
	   bet[wf.in3(k,j,i)]=tridag_mid[wf.in4(1,k,j,i)];
      	   wf_tmp[wf.in4(1,k,j,i) ] = wf_rightside[wf.in4(1,k,j,i)]/bet[wf.in3(k,j,i)];
	  }}}
  else
    {
      MPI_Recv(&u_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag1_2,MPI_COMM_WORLD, &status);
      MPI_Recv(&bet_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag2_2,MPI_COMM_WORLD, &status); 
     #pragma omp parallel for  num_threads(nthreads)  collapse(3)
     for( int k = 1 ; k <= wf.n3/2 ; k++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
	   index=(k-1)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+(i-wf.n1/2)-1;
   	   gam[wf.in4(1,k,j,i)]=tridag_upp/bet_from_left[index];
           bet[wf.in3(k,j,i)]  = tridag_mid[wf.in4(1,k,j,i) ];
           bet[wf.in3(k,j,i)] -= tridag_low*gam[wf.in4(1,k,j,i)];
           wf_tmp[wf.in4(1,k,j,i)] = ( wf_rightside[wf.in4(1,k,j,i) ]-tridag_low*u_from_left[index] )/bet[wf.in3(k,j,i)];
     }}     
    }
 #pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int k = 1 ; k <= wf.n3/2 ; k++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
   	for( int l = 2 ; l <= wf.n4 ; l++ )
   	  {
      	   gam[wf.in4(l,k,j,i)] = tridag_upp/bet[wf.in3(k,j,i)];
           bet[wf.in3(k,j,i)]  = tridag_mid[wf.in4(l,k,j,i)];
      	   bet[wf.in3(k,j,i)] -= tridag_low*gam[wf.in4(l,k,j,i)];
           wf_tmp[wf.in4(l,k,j,i)] = ( wf_rightside[wf.in4(l,k,j,i)]-tridag_low*wf_tmp[wf.in4(l-1,k,j,i)] )/bet[wf.in3(k,j,i)];
    }}}
  

  if (myid < numprocs-1 )
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int k = 1 ; k <= wf.n3/2 ; k++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
	   index=(k-1)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+(i-wf.n1/2)-1;
           u_to_right[index]=wf_tmp[wf.in4(wf.n4,k,j,i)];
           bet_to_right[index]=bet[wf.in3(k,j,i)];
          }
       }
      MPI_Send(&u_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag1_2,MPI_COMM_WORLD);
      MPI_Send(&bet_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag2_2,MPI_COMM_WORLD);
    }

//----------------------area 3 in Low part
  if (myid==0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int k = 1 ; k <= wf.n3/2 ; k++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
	   bet[wf.in3(k,j,i)]=tridag_mid[wf.in4(1,k,j,i)];
      	   wf_tmp[wf.in4(1,k,j,i) ] = wf_rightside[wf.in4(1,k,j,i)]/bet[wf.in3(k,j,i)];
	  }}}
  else
    {
      MPI_Recv(&u_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag1_3,MPI_COMM_WORLD, &status);
      MPI_Recv(&bet_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag2_3,MPI_COMM_WORLD, &status); 
     #pragma omp parallel for  num_threads(nthreads)  collapse(3)
     for( int k = 1 ; k <= wf.n3/2 ; k++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
	   index=(k-1)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-1;
   	   gam[wf.in4(1,k,j,i)]=tridag_upp/bet_from_left[index];
           bet[wf.in3(k,j,i)]  = tridag_mid[wf.in4(1,k,j,i) ];
           bet[wf.in3(k,j,i)] -= tridag_low*gam[wf.in4(1,k,j,i)];
           wf_tmp[wf.in4(1,k,j,i)] = ( wf_rightside[wf.in4(1,k,j,i) ]-tridag_low*u_from_left[index] )/bet[wf.in3(k,j,i)];
     }}     
    }
 #pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int k = 1 ; k <= wf.n3/2 ; k++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
   	for( int l = 2 ; l <= wf.n4 ; l++ )
   	  {
      	   gam[wf.in4(l,k,j,i)] = tridag_upp/bet[wf.in3(k,j,i)];
           bet[wf.in3(k,j,i)]  = tridag_mid[wf.in4(l,k,j,i)];
      	   bet[wf.in3(k,j,i)] -= tridag_low*gam[wf.in4(l,k,j,i)];
           wf_tmp[wf.in4(l,k,j,i)] = ( wf_rightside[wf.in4(l,k,j,i)]-tridag_low*wf_tmp[wf.in4(l-1,k,j,i)] )/bet[wf.in3(k,j,i)];
    }}}
  

  if (myid < numprocs-1 )
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int k = 1 ; k <= wf.n3/2 ; k++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
	   index=(k-1)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-1;
           u_to_right[index]=wf_tmp[wf.in4(wf.n4,k,j,i)];
           bet_to_right[index]=bet[wf.in3(k,j,i)];
          }
       }
      MPI_Send(&u_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag1_3,MPI_COMM_WORLD);
      MPI_Send(&bet_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag2_3,MPI_COMM_WORLD);
    }

//----------------------area 4 in Low part
  if (myid==0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int k = 1 ; k <= wf.n3/2 ; k++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
	   bet[wf.in3(k,j,i)]=tridag_mid[wf.in4(1,k,j,i)];
      	   wf_tmp[wf.in4(1,k,j,i) ] = wf_rightside[wf.in4(1,k,j,i)]/bet[wf.in3(k,j,i)];
	  }}}
  else
    {
      MPI_Recv(&u_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag1_4,MPI_COMM_WORLD, &status);
      MPI_Recv(&bet_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag2_4,MPI_COMM_WORLD, &status); 
     #pragma omp parallel for  num_threads(nthreads)  collapse(3)
	 for( int k = 1 ; k <= wf.n3/2 ; k++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
	   index=(k-1)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-wf.n1/2-1;
   	   gam[wf.in4(1,k,j,i)]=tridag_upp/bet_from_left[index];
       bet[wf.in3(k,j,i)]  = tridag_mid[wf.in4(1,k,j,i) ];
       bet[wf.in3(k,j,i)] -= tridag_low*gam[wf.in4(1,k,j,i)];
       wf_tmp[wf.in4(1,k,j,i)] = ( wf_rightside[wf.in4(1,k,j,i) ]-tridag_low*u_from_left[index] )/bet[wf.in3(k,j,i)];
     }}     
    }
 #pragma omp parallel for  num_threads(nthreads) collapse(3)
   for( int k = 1 ; k <= wf.n3/2 ; k++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
   	for( int l = 2 ; l <= wf.n4 ; l++ )
   	  {
      	   gam[wf.in4(l,k,j,i)] = tridag_upp/bet[wf.in3(k,j,i)];
           bet[wf.in3(k,j,i)]  = tridag_mid[wf.in4(l,k,j,i)];
      	   bet[wf.in3(k,j,i)] -= tridag_low*gam[wf.in4(l,k,j,i)];
           wf_tmp[wf.in4(l,k,j,i)] = ( wf_rightside[wf.in4(l,k,j,i)]-tridag_low*wf_tmp[wf.in4(l-1,k,j,i)] )/bet[wf.in3(k,j,i)];
    }}}
  

  if (myid < numprocs-1 )
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int k = 1 ; k <= wf.n3/2 ; k++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
	   index=(k-1)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-wf.n1/2-1;
           u_to_right[index]=wf_tmp[wf.in4(wf.n4,k,j,i)];
           bet_to_right[index]=bet[wf.in3(k,j,i)];
          }
       }
      MPI_Send(&u_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag1_4,MPI_COMM_WORLD);
      MPI_Send(&bet_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag2_4,MPI_COMM_WORLD);
    }

//----------------------area 5 in Low part
  if (myid==0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int k = 1+wf.n3/2 ; k <= wf.n3 ; k++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
	   bet[wf.in3(k,j,i)]=tridag_mid[wf.in4(1,k,j,i)];
      	   wf_tmp[wf.in4(1,k,j,i) ] = wf_rightside[wf.in4(1,k,j,i)]/bet[wf.in3(k,j,i)];
	  }}}
  else
    {
      MPI_Recv(&u_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag1_5,MPI_COMM_WORLD, &status);
      MPI_Recv(&bet_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag2_5,MPI_COMM_WORLD, &status); 
     #pragma omp parallel for  num_threads(nthreads) collapse(3) 
     for( int k = 1+wf.n3/2 ; k <= wf.n3 ; k++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
	   index=(k-1-wf.n3/2)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+i-1;
   	   gam[wf.in4(1,k,j,i)]=tridag_upp/bet_from_left[index];
           bet[wf.in3(k,j,i)]  = tridag_mid[wf.in4(1,k,j,i) ];
           bet[wf.in3(k,j,i)] -= tridag_low*gam[wf.in4(1,k,j,i)];
           wf_tmp[wf.in4(1,k,j,i)] = ( wf_rightside[wf.in4(1,k,j,i) ]-tridag_low*u_from_left[index] )/bet[wf.in3(k,j,i)];
     }}     
    }
 #pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int k = 1+wf.n3/2 ; k <= wf.n3 ; k++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
   	for( int l = 2 ; l <= wf.n4 ; l++ )
   	  {
      	   gam[wf.in4(l,k,j,i)] = tridag_upp/bet[wf.in3(k,j,i)];
           bet[wf.in3(k,j,i)]  = tridag_mid[wf.in4(l,k,j,i)];
      	   bet[wf.in3(k,j,i)] -= tridag_low*gam[wf.in4(l,k,j,i)];
           wf_tmp[wf.in4(l,k,j,i)] = ( wf_rightside[wf.in4(l,k,j,i)]-tridag_low*wf_tmp[wf.in4(l-1,k,j,i)] )/bet[wf.in3(k,j,i)];
    }}}
  

  if (myid < numprocs-1 )
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int k = 1+wf.n3/2 ; k <= wf.n3 ; k++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
	   index=(k-1-wf.n3/2)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+i-1;
           u_to_right[index]=wf_tmp[wf.in4(wf.n4,k,j,i)];
           bet_to_right[index]=bet[wf.in3(k,j,i)];
          }
       }
      MPI_Send(&u_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag1_5,MPI_COMM_WORLD);
      MPI_Send(&bet_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag2_5,MPI_COMM_WORLD);
    }

//----------------------area 6 in Low part
  if (myid==0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int k = 1+wf.n3/2 ; k <= wf.n3 ; k++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
	   bet[wf.in3(k,j,i)]=tridag_mid[wf.in4(1,k,j,i)];
      	   wf_tmp[wf.in4(1,k,j,i) ] = wf_rightside[wf.in4(1,k,j,i)]/bet[wf.in3(k,j,i)];
	  }}}
  else
    {
      MPI_Recv(&u_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag1_6,MPI_COMM_WORLD, &status);
      MPI_Recv(&bet_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag2_6,MPI_COMM_WORLD, &status); 
     #pragma omp parallel for  num_threads(nthreads)  collapse(3)
     for( int k = 1+wf.n3/2 ; k <= wf.n3 ; k++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
	   index=(k-1-wf.n3/2)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+(i-wf.n1/2)-1;
   	   gam[wf.in4(1,k,j,i)]=tridag_upp/bet_from_left[index];
           bet[wf.in3(k,j,i)]  = tridag_mid[wf.in4(1,k,j,i) ];
           bet[wf.in3(k,j,i)] -= tridag_low*gam[wf.in4(1,k,j,i)];
           wf_tmp[wf.in4(1,k,j,i)] = ( wf_rightside[wf.in4(1,k,j,i) ]-tridag_low*u_from_left[index] )/bet[wf.in3(k,j,i)];
     }}     
    }
 #pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int k = 1+wf.n3/2 ; k <= wf.n3 ; k++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
   	for( int l = 2 ; l <= wf.n4 ; l++ )
   	  {
      	   gam[wf.in4(l,k,j,i)] = tridag_upp/bet[wf.in3(k,j,i)];
           bet[wf.in3(k,j,i)]  = tridag_mid[wf.in4(l,k,j,i)];
      	   bet[wf.in3(k,j,i)] -= tridag_low*gam[wf.in4(l,k,j,i)];
           wf_tmp[wf.in4(l,k,j,i)] = ( wf_rightside[wf.in4(l,k,j,i)]-tridag_low*wf_tmp[wf.in4(l-1,k,j,i)] )/bet[wf.in3(k,j,i)];
    }}}
  

  if (myid < numprocs-1 )
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int k = 1+wf.n3/2 ; k <= wf.n3 ; k++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
	   index=(k-1-wf.n3/2)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+(i-wf.n1/2)-1;
           u_to_right[index]=wf_tmp[wf.in4(wf.n4,k,j,i)];
           bet_to_right[index]=bet[wf.in3(k,j,i)];
          }
       }
      MPI_Send(&u_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag1_6,MPI_COMM_WORLD);
      MPI_Send(&bet_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag2_6,MPI_COMM_WORLD);
    }

//----------------------area 7 in Low part
  if (myid==0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int k = 1+wf.n3/2 ; k <= wf.n3 ; k++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
	   bet[wf.in3(k,j,i)]=tridag_mid[wf.in4(1,k,j,i)];
      	   wf_tmp[wf.in4(1,k,j,i) ] = wf_rightside[wf.in4(1,k,j,i)]/bet[wf.in3(k,j,i)];
	  }}}
  else
    {
      MPI_Recv(&u_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag1_7,MPI_COMM_WORLD, &status);
      MPI_Recv(&bet_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag2_7,MPI_COMM_WORLD, &status); 
     #pragma omp parallel for  num_threads(nthreads)  collapse(3)
     for( int k = 1+wf.n3/2 ; k <= wf.n3 ; k++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
	   index=(k-1-wf.n3/2)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-1;
   	   gam[wf.in4(1,k,j,i)]=tridag_upp/bet_from_left[index];
           bet[wf.in3(k,j,i)]  = tridag_mid[wf.in4(1,k,j,i) ];
           bet[wf.in3(k,j,i)] -= tridag_low*gam[wf.in4(1,k,j,i)];
           wf_tmp[wf.in4(1,k,j,i)] = ( wf_rightside[wf.in4(1,k,j,i) ]-tridag_low*u_from_left[index] )/bet[wf.in3(k,j,i)];
     }}     
    }
 #pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int k = 1+wf.n3/2 ; k <= wf.n3 ; k++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
   	for( int l = 2 ; l <= wf.n4 ; l++ )
   	  {
      	   gam[wf.in4(l,k,j,i)] = tridag_upp/bet[wf.in3(k,j,i)];
           bet[wf.in3(k,j,i)]  = tridag_mid[wf.in4(l,k,j,i)];
      	   bet[wf.in3(k,j,i)] -= tridag_low*gam[wf.in4(l,k,j,i)];
           wf_tmp[wf.in4(l,k,j,i)] = ( wf_rightside[wf.in4(l,k,j,i)]-tridag_low*wf_tmp[wf.in4(l-1,k,j,i)] )/bet[wf.in3(k,j,i)];
    }}}
  

  if (myid < numprocs-1 )
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int k = 1+wf.n3/2 ; k <= wf.n3 ; k++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
	   index=(k-1-wf.n3/2)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-1;
           u_to_right[index]=wf_tmp[wf.in4(wf.n4,k,j,i)];
           bet_to_right[index]=bet[wf.in3(k,j,i)];
          }
       }
      MPI_Send(&u_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag1_7,MPI_COMM_WORLD);
      MPI_Send(&bet_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag2_7,MPI_COMM_WORLD);
    }

//----------------------area 8 in Low part
  if (myid==0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
     for( int k = 1+wf.n3/2 ; k <= wf.n3 ; k++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
	   bet[wf.in3(k,j,i)]=tridag_mid[wf.in4(1,k,j,i)];
      	   wf_tmp[wf.in4(1,k,j,i) ] = wf_rightside[wf.in4(1,k,j,i)]/bet[wf.in3(k,j,i)];
	  }}}
  else
    {
      MPI_Recv(&u_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag1_8,MPI_COMM_WORLD, &status);
      MPI_Recv(&bet_from_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag2_8,MPI_COMM_WORLD, &status); 
     #pragma omp parallel for  num_threads(nthreads)  collapse(3)
	 for( int k = 1+wf.n3/2 ; k <= wf.n3 ; k++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
	   index=(k-1-wf.n3/2)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-wf.n1/2-1;
   	   gam[wf.in4(1,k,j,i)]=tridag_upp/bet_from_left[index];
       bet[wf.in3(k,j,i)]  = tridag_mid[wf.in4(1,k,j,i) ];
       bet[wf.in3(k,j,i)] -= tridag_low*gam[wf.in4(1,k,j,i)];
       wf_tmp[wf.in4(1,k,j,i)] = ( wf_rightside[wf.in4(1,k,j,i) ]-tridag_low*u_from_left[index] )/bet[wf.in3(k,j,i)];
     }}     
    }
 #pragma omp parallel for  num_threads(nthreads) collapse(3)
   for( int k = 1+wf.n3/2 ; k <= wf.n3 ; k++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
   	for( int l = 2 ; l <= wf.n4 ; l++ )
   	  {
      	   gam[wf.in4(l,k,j,i)] = tridag_upp/bet[wf.in3(k,j,i)];
           bet[wf.in3(k,j,i)]  = tridag_mid[wf.in4(l,k,j,i)];
      	   bet[wf.in3(k,j,i)] -= tridag_low*gam[wf.in4(l,k,j,i)];
           wf_tmp[wf.in4(l,k,j,i)] = ( wf_rightside[wf.in4(l,k,j,i)]-tridag_low*wf_tmp[wf.in4(l-1,k,j,i)] )/bet[wf.in3(k,j,i)];
    }}}
  

  if (myid < numprocs-1 )
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int k = 1+wf.n3/2 ; k <= wf.n3 ; k++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
	   index=(k-1-wf.n3/2)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-wf.n1/2-1;
           u_to_right[index]=wf_tmp[wf.in4(wf.n4,k,j,i)];
           bet_to_right[index]=bet[wf.in3(k,j,i)];
          }
       }
      MPI_Send(&u_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag1_8,MPI_COMM_WORLD);
      MPI_Send(&bet_to_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag2_8,MPI_COMM_WORLD);
    }



//----------------------area 1 in up part
  if(myid != numprocs-1)
    {
      MPI_Recv(&gam_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag3_1,MPI_COMM_WORLD, &status);
      MPI_Recv(&u_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag4_1,MPI_COMM_WORLD, &status);
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	for( int k = 1 ; k <= wf.n3/2 ; k++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
           index=(k-1)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+i-1;
           wf_tmp[wf.in4(wf.n4,k,j,i)] -= gam_from_right[index]*u_from_right[index]; 
          }}
    }
  #pragma omp parallel for  num_threads(nthreads) collapse(3)
  	for( int k = 1 ; k <= wf.n3/2 ; k++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
           for( int l = wf.n4-1 ; l >= 1 ; l-- )
             {
              wf_tmp[ wf.in4(l,k,j,i)] -= gam[ wf.in4(l+1,k,j,i) ]*wf_tmp[wf.in4(l+1,k,j,i)];
     }}}

  if (myid != 0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int k = 1 ; k <= wf.n3/2 ; k++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
            index=(k-1)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+i-1;
            gam_to_left[index] = gam[wf.in4(1,k,j,i)];
      	    u_to_left[index] = wf_tmp[wf.in4(1,k,j,i)];
          }
       }
      MPI_Send(&gam_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag3_1,MPI_COMM_WORLD);
      MPI_Send(&u_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag4_1,MPI_COMM_WORLD);
    }

//----------------------area 2 in up part
  if(myid != numprocs-1)
    {
      MPI_Recv(&gam_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag3_2,MPI_COMM_WORLD, &status);
      MPI_Recv(&u_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag4_2,MPI_COMM_WORLD, &status);
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int k = 1 ; k <= wf.n3/2 ; k++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
           index=(k-1)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+i-wf.n1/2-1;
           wf_tmp[wf.in4(wf.n4,k,j,i)] -= gam_from_right[index]*u_from_right[index]; 
          }}
    }
  #pragma omp parallel for  num_threads(nthreads) collapse(3)
  	for( int k = 1 ; k <= wf.n3/2 ; k++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
           for( int l = wf.n4-1 ; l >= 1 ; l-- )
             {
              wf_tmp[ wf.in4(l,k,j,i)] -= gam[ wf.in4(l+1,k,j,i) ]*wf_tmp[wf.in4(l+1,k,j,i)];
     }}}

  if (myid != 0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int k = 1 ; k <= wf.n3/2 ; k++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
            index=(k-1)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+i-wf.n1/2-1;
            gam_to_left[index] = gam[wf.in4(1,k,j,i)];
      	    u_to_left[index] = wf_tmp[wf.in4(1,k,j,i)];
          }
       }
      MPI_Send(&gam_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag3_2,MPI_COMM_WORLD);
      MPI_Send(&u_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag4_2,MPI_COMM_WORLD);
    }

//----------------------area 3 in up part
  if(myid != numprocs-1)
    {
      MPI_Recv(&gam_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag3_3,MPI_COMM_WORLD, &status);
      MPI_Recv(&u_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag4_3,MPI_COMM_WORLD, &status);
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int k = 1 ; k <= wf.n3/2 ; k++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
           index=(k-1)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-1;
           wf_tmp[wf.in4(wf.n4,k,j,i)] -= gam_from_right[index]*u_from_right[index]; 
          }}
    }
  #pragma omp parallel for  num_threads(nthreads) collapse(3)
  	for( int k = 1 ; k <= wf.n3/2 ; k++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
           for( int l = wf.n4-1 ; l >= 1 ; l-- )
             {
              wf_tmp[ wf.in4(l,k,j,i)] -= gam[ wf.in4(l+1,k,j,i) ]*wf_tmp[wf.in4(l+1,k,j,i)];
     }}}

  if (myid != 0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int k = 1 ; k <= wf.n3/2 ; k++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
            index=(k-1)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-1;
            gam_to_left[index] = gam[wf.in4(1,k,j,i)];
      	    u_to_left[index] = wf_tmp[wf.in4(1,k,j,i)];
          }
       }
      MPI_Send(&gam_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag3_3,MPI_COMM_WORLD);
      MPI_Send(&u_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag4_3,MPI_COMM_WORLD);
    }

//----------------------area 4 in up part
  if(myid != numprocs-1)
    {
      MPI_Recv(&gam_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag3_4,MPI_COMM_WORLD, &status);
      MPI_Recv(&u_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag4_4,MPI_COMM_WORLD, &status);
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int k = 1 ; k <= wf.n3/2 ; k++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
           index=(k-1)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-wf.n1/2-1;
           wf_tmp[wf.in4(wf.n4,k,j,i)] -= gam_from_right[index]*u_from_right[index]; 
          }}
    }
  #pragma omp parallel for  num_threads(nthreads) collapse(3)
  	for( int k = 1 ; k <= wf.n3/2 ; k++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
           for( int l = wf.n4-1 ; l >= 1 ; l-- )
             {
              wf_tmp[ wf.in4(l,k,j,i)] -= gam[ wf.in4(l+1,k,j,i) ]*wf_tmp[wf.in4(l+1,k,j,i)];
     }}}

  if (myid != 0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int k = 1 ; k <= wf.n3/2 ; k++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
            index=(k-1)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-wf.n1/2-1;
            gam_to_left[index] = gam[wf.in4(1,k,j,i)];
      	    u_to_left[index] = wf_tmp[wf.in4(1,k,j,i)];
          }
       }
      MPI_Send(&gam_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag3_4,MPI_COMM_WORLD);
      MPI_Send(&u_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag4_4,MPI_COMM_WORLD);
    }


//----------------------area 5 in up part
  if(myid != numprocs-1)
    {
      MPI_Recv(&gam_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag3_5,MPI_COMM_WORLD, &status);
      MPI_Recv(&u_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag4_5,MPI_COMM_WORLD, &status);
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	for( int k = 1+wf.n3/2 ; k <= wf.n3 ; k++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
           index=(k-1-wf.n3/2)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+i-1;
           wf_tmp[wf.in4(wf.n4,k,j,i)] -= gam_from_right[index]*u_from_right[index]; 
          }}
    }
  #pragma omp parallel for  num_threads(nthreads) collapse(3)
  	for( int k = 1+wf.n3/2 ; k <= wf.n3 ; k++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
           for( int l = wf.n4-1 ; l >= 1 ; l-- )
             {
              wf_tmp[ wf.in4(l,k,j,i)] -= gam[ wf.in4(l+1,k,j,i) ]*wf_tmp[wf.in4(l+1,k,j,i)];
     }}}

  if (myid != 0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int k = 1+wf.n3/2 ; k <= wf.n3 ; k++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
            index=(k-1-wf.n3/2)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+i-1;
            gam_to_left[index] = gam[wf.in4(1,k,j,i)];
      	    u_to_left[index] = wf_tmp[wf.in4(1,k,j,i)];
          }
       }
      MPI_Send(&gam_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag3_5,MPI_COMM_WORLD);
      MPI_Send(&u_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag4_5,MPI_COMM_WORLD);
    }

//----------------------area 6 in up part
  if(myid != numprocs-1)
    {
      MPI_Recv(&gam_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag3_6,MPI_COMM_WORLD, &status);
      MPI_Recv(&u_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag4_6,MPI_COMM_WORLD, &status);
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int k = 1+wf.n3/2 ; k <= wf.n3 ; k++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
           index=(k-1-wf.n3/2)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+i-wf.n1/2-1;
           wf_tmp[wf.in4(wf.n4,k,j,i)] -= gam_from_right[index]*u_from_right[index]; 
          }}
    }
  #pragma omp parallel for  num_threads(nthreads) collapse(3)
  	for( int k = 1+wf.n3/2 ; k <= wf.n3 ; k++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
           for( int l = wf.n4-1 ; l >= 1 ; l-- )
             {
              wf_tmp[ wf.in4(l,k,j,i)] -= gam[ wf.in4(l+1,k,j,i) ]*wf_tmp[wf.in4(l+1,k,j,i)];
     }}}

  if (myid != 0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int k = 1+wf.n3/2 ; k <= wf.n3 ; k++ )
     for( int j = 1 ; j <= wf.n2/2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
            index=(k-1-wf.n3/2)*wf.n2*wf.n1/4+(j-1)*wf.n1/2+i-wf.n1/2-1;
            gam_to_left[index] = gam[wf.in4(1,k,j,i)];
      	    u_to_left[index] = wf_tmp[wf.in4(1,k,j,i)];
          }
       }
      MPI_Send(&gam_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag3_6,MPI_COMM_WORLD);
      MPI_Send(&u_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag4_6,MPI_COMM_WORLD);
    }

//----------------------area 7 in up part
  if(myid != numprocs-1)
    {
      MPI_Recv(&gam_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag3_7,MPI_COMM_WORLD, &status);
      MPI_Recv(&u_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag4_7,MPI_COMM_WORLD, &status);
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int k = 1+wf.n3/2 ; k <= wf.n3 ; k++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
           index=(k-1-wf.n3/2)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-1;
           wf_tmp[wf.in4(wf.n4,k,j,i)] -= gam_from_right[index]*u_from_right[index]; 
          }}
    }
  #pragma omp parallel for  num_threads(nthreads) collapse(3)
  	for( int k = 1+wf.n3/2 ; k <= wf.n3 ; k++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
           for( int l = wf.n4-1 ; l >= 1 ; l-- )
             {
              wf_tmp[ wf.in4(l,k,j,i)] -= gam[ wf.in4(l+1,k,j,i) ]*wf_tmp[wf.in4(l+1,k,j,i)];
     }}}

  if (myid != 0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int k = 1+wf.n3/2 ; k <= wf.n3 ; k++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1/2 ; i++ )
   	  {
            index=(k-1-wf.n3/2)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-1;
            gam_to_left[index] = gam[wf.in4(1,k,j,i)];
      	    u_to_left[index] = wf_tmp[wf.in4(1,k,j,i)];
          }
       }
      MPI_Send(&gam_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag3_7,MPI_COMM_WORLD);
      MPI_Send(&u_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag4_7,MPI_COMM_WORLD);
    }

//----------------------area 8 in up part
  if(myid != numprocs-1)
    {
      MPI_Recv(&gam_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag3_8,MPI_COMM_WORLD, &status);
      MPI_Recv(&u_from_right[0],num_comm,MPI_DOUBLE_COMPLEX,right,tag4_8,MPI_COMM_WORLD, &status);
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int k = 1+wf.n3/2 ; k <= wf.n3 ; k++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
           index=(k-1-wf.n3/2)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-wf.n1/2-1;
           wf_tmp[wf.in4(wf.n4,k,j,i)] -= gam_from_right[index]*u_from_right[index]; 
          }}
    }
  #pragma omp parallel for  num_threads(nthreads) collapse(3)
  	for( int k = 1+wf.n3/2 ; k <= wf.n3 ; k++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
           for( int l = wf.n4-1 ; l >= 1 ; l-- )
             {
              wf_tmp[ wf.in4(l,k,j,i)] -= gam[ wf.in4(l+1,k,j,i) ]*wf_tmp[wf.in4(l+1,k,j,i)];
     }}}

  if (myid != 0)
    {
#pragma omp parallel for  num_threads(nthreads) collapse(3)
	 for( int k = 1+wf.n3/2 ; k <= wf.n3 ; k++ )
     for( int j = 1+wf.n2/2 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1+wf.n1/2 ; i <= wf.n1 ; i++ )
   	  {
            index=(k-1-wf.n3/2)*wf.n2*wf.n1/4+(j-wf.n2/2-1)*wf.n1/2+i-wf.n1/2-1;
            gam_to_left[index] = gam[wf.in4(1,k,j,i)];
      	    u_to_left[index] = wf_tmp[wf.in4(1,k,j,i)];
          }
       }
      MPI_Send(&gam_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag3_8,MPI_COMM_WORLD);
      MPI_Send(&u_to_left[0],num_comm,MPI_DOUBLE_COMPLEX,left,tag4_8,MPI_COMM_WORLD);
    }



//-------------------final step
#pragma omp parallel for  num_threads(nthreads) collapse(4)
   for( int l = 1 ; l <= wf.n4 ; l++ )
   for( int k = 1 ; k <= wf.n3 ; k++ )
     {
     for( int j = 1 ; j <= wf.n2 ; j++ )
       {
   	for( int i = 1 ; i <= wf.n1 ; i++ )
   	  {
           wf.wave[ wf.in4(l, k , j , i ) ] = wf_tmp[wf.in4(l,k,j,i)];
   	  }}} 
}// end of tridag_MPI


inline void Tridag_FastAM( const complex &tridag_lower_diagonal , const vector<complex> &tridag_middle_diagonal , const complex &tridag_upper_diagonal, const vector<complex> &wf_1D_rightside ,
		  vector<complex> &wf_1D_tmp2 , vector<complex> &gam , int length)
{   


  //  if ( tridag_middle_diagonal[ 1 ] == 0. )
  //   nerror( "Error: tridag_middle_diagonal[1]==0 in tridag!" );
    




  complex bet=tridag_middle_diagonal[1];
  wf_1D_tmp2[ 1 ] = wf_1D_rightside[ 1 ]/bet;
  for( unsigned int index = 2 ; index <= length ; ++index )
    {
      gam[ index ] = tridag_upper_diagonal/bet;
      bet  = tridag_middle_diagonal[ index ];
      bet -= tridag_lower_diagonal*gam[ index ];
      //      if (bet==0.0)
      //	 nerror("Error: bet==0 in tridag!");
	

      wf_1D_tmp2[ index ] = ( wf_1D_rightside[ index ]-tridag_lower_diagonal*wf_1D_tmp2[ index-1 ] )/bet;
    }
  


  for( unsigned int index = length-1 ; index >= 1 ; --index )
       wf_1D_tmp2[ index ] -= gam[ index+1 ]*wf_1D_tmp2[ index+1 ];
 
  
  
}// end of tridag



/**
 *  The tridag routine to solve the threediagonal system of linear equation for ADAPATATIVE MESH
 */
inline void TridagAM( const vector<complex> &tridag_lower_diagonal , const vector<complex> &tridag_middle_diagonal , const vector<complex>&tridag_upper_diagonal, const vector<complex> &wf_1D_rightside ,
		      vector<complex> &wf_1D_tmp2 , vector<complex> &gam, int length )
{   


/**
 *  forth substitution
 */

  complex bet=tridag_middle_diagonal[ 1 ];
  wf_1D_tmp2[ 1 ] = wf_1D_rightside[ 1 ]/bet;
  for( unsigned int index = 2 ; index <= length; ++index )
    {
      gam[ index ] = tridag_upper_diagonal[ index-1 ]/bet;
      bet  = tridag_middle_diagonal[ index ];
      bet -= tridag_lower_diagonal[ index ]*gam[ index ];
      /*      if ( bet == 0.0 )
	{
	 cout<<"Error: bet==0 in tridag! at "<<index<<"\n";
         cout<<"tridag_middle_diagonal="<<tridag_middle_diagonal[ index ]<<"\n";
         cout<<"tridag_lower_diagonal="<<tridag_lower_diagonal[ index ]<<"\n";
         cout<<"gam="<<gam[ index ]<<"\n";
	 } */

      wf_1D_tmp2[ index ] = ( wf_1D_rightside[ index ]-tridag_lower_diagonal[ index ]*wf_1D_tmp2[ index-1 ] )/bet;
    }
  
/**
 *  back substitution
 */
  for( unsigned int index = length-1 ; index >= 1 ; --index )
       wf_1D_tmp2[ index ] -= gam[ index+1 ]*wf_1D_tmp2[ index+1 ];
 
}// end of tridag











/*
void tridag(complex a,vector<complex> b,complex c,vector<complex> r,vector<complex> &u, int n)
{
    int j;
    complex bet, *gam;
    gam=(complex *)malloc((n+2)*sizeof (complex));
    if(!gam) cout<<"\nerror de colocacion en la asignacion complex";
    if(b[1]==0.0) nrerror ("error 1 en tridag");
    u[1]=r[1]/(bet=b[1]);
    for(j=2;j<=n;j++)
    {
		gam[j]=c/bet;
		bet=b[j]-a*gam[j];
		if (bet==0.0) nrerror("error 2 en tridag");
		u[j]=(r[j]-a*u[j-1])/bet;
    }
    for(j=(n-1);j>=1;j--)
		u[j]-=gam[j+1]*u[j+1];
	
    //for(j=1;j<=n;j++)
    //  printf("%d %e\n",j,real(u[j]));
    free(gam);
    //    printf("I am here\n");
	
}
*/

inline void skip_comments(std::istream & is)
{
  bool comments=true;
  char c;
  std::string linebuffer;
  while (comments) {
    is.get(c);
    if ('#'==c)
      getline(is,linebuffer);
    else {
        is.unget();
        comments = false;
    }
  }
}

inline void set_precision(std::ostream &os, size_t p) {
  os.precision(p);
  os.setf(std::ios::scientific, std::ios::floatfield);
  os.setf(std::ios::showpoint);
}


/**
 *  convert from integer to string
 */

inline string Conversion_To_String( int i )
{
  ostringstream ost;
  ost << i;
  return ost.str();
} // to_string integer


/**
 *  convert from double to string
 */
 
inline string Conversion_To_String( double i )
{
  ostringstream ost;
  ost << i;
  return ost.str();
} // to_string double
 


#endif	/* UTILS_H */
