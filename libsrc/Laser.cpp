// $Id$
#include <iostream>
#include <string>
#include <vector>
#include <fstream>      // for Field::Field(filename)
#include "utils.h"      // for skip_comments(std::istream & is)
#include "constant.h"   // for value of pi
#include "Laser.h"

using std::cout;
using std::cerr;
using std::vector;
using std::string;
using std::endl;

double Laser::Ex(const double &t) const { return E(0, t); }; /**<setup alias*/
double Laser::Ey(const double &t) const { return E(1, t); }; /**<setup alias*/


/* --------------------- Helper function --------------------*/
void Integrate_E_trapezoid_rule(vector<double> &E, vector<double>& t, vector<double>& A) {
  size_t points=E.size();
  if( t.size() != points ) {
    cerr<<"Integrate_E_trapezoid_rule: length of E does not match length of t! ("<<points<<" != "<<t.size()<<")\n";
    exit(1);
  }
  // resize A vector to match length of E (both for x and y, respectively)
  A.resize(points);

  double vector_potential_temp = A[0];
  size_t index=1;
  while( index < points ) {
    // integrate with trapezoidal rule.
    double time_step=t[index]-t[index-1];
    vector_potential_temp += - (E[index]+E[index-1])/2.*time_step*lightC_au;
    A[index]=vector_potential_temp;
    ++index;
  } 
  if(index<2) { // clean up A if nothing has been done.
    A.clear();
  }
}

void Integrate_A_trapezoid_rule(vector<double> &A, vector<double>& t, vector<double>& Integral_A) {
  size_t points=A.size();
  if( t.size() != points ) {
    cerr<<"Integrate_A_trapezoid_rule: length of A does not match length of t! ("<<points<<" != "<<t.size()<<")\n";
    exit(1);
  }
  // resize Integral_A vector to match length of A (both for x and y, respectively)
  Integral_A.resize(points);

  Integral_A[points-1] = Integral_A[0]; // HACK: final value handed over via Integral_A[0] instead of Integral_A[points-1]
  double potential_integral_temp = Integral_A[points-1];
  size_t index=points-1;
  while( index > 0 ) {
    // integrate with trapezoidal rule.
    double time_step=t[index]-t[index-1];
    potential_integral_temp += (A[index]+A[index-1])/2.*time_step;
    Integral_A[index-1]=potential_integral_temp;
    --index;
  } 
  if(points<2) { // clean up Integral_A if nothing has been done.
    Integral_A.clear();
  }
}


void Integrate_A2_trapezoid_rule(vector<double> &A, vector<double>& t, vector<double>& Integral_A2) {
  size_t points=A.size();
  if( t.size() != points ) {
    cerr<<"Integrate_A_trapezoid_rule: length of A does not match length of t! ("<<points<<" != "<<t.size()<<")\n";
    exit(1);
  }
  // resize Integral_A vector to match length of A (both for x and y, respectively)
  Integral_A2.resize(points);

  Integral_A2[points-1] = Integral_A2[0]; // HACK: final value handed over via Integral_A[0] instead of Integral_A[points-1]
  double potential_integral_temp = Integral_A2[points-1];
  size_t index=points-1;
  while( index > 0 ) {
    // integrate with trapezoidal rule.
    double time_step=t[index]-t[index-1];
    potential_integral_temp += (A[index]*A[index]+A[index-1]* A[index-1] )/2.*time_step;
    Integral_A2[index-1]=potential_integral_temp;
    --index;
  }
  if(points<2) { // clean up Integral_A if nothing has been done.
    Integral_A2.clear();
  }
}


/* --------------------- Implementation details --------------------*/

Field::Field(Laser &l1, double initial_Ax = 0, double initial_Ay = 0)
// for shorthand notation Ex := E[0] :
: Ex(E[0]), Ey(E[1]), Ax(A[0]), Ay(A[1]),
  Integral_Ax(Integral_A[0]), Integral_Ay(Integral_A[1]),Integral_Ax2(Integral_A2[0]), Integral_Ay2(Integral_A2[1])
{
  linear_polarized = ( l1.ellipticity()==0 && l1.pol_angle() == 0);

  size_t steps;
  steps = ( size_t ) std::ceil( l1.laser_duration/l1.time_step );

  t.resize( steps , 0 );
  Ax.push_back(initial_Ax);
  Ay.push_back(initial_Ay);

  Integral_Ax.push_back(0);    // HACK: Integral_A(t_f) := 0; But how else??
  Integral_Ay.push_back(0);    // HACK: Integral_A(t_f) := 0;

  Integral_Ax2.push_back(0);    // HACK: Integral_A(t_f) := 0; But how else??
  Integral_Ay2.push_back(0);    // HACK: Integral_A(t_f) := 0;

  for ( int k=0; k<2; ++k ) {
    E[k].resize( steps , 0 );
    A[k].resize( steps , 0 );
  };

  size_t index=0;
  while( index < steps )
  {
    // the simulation time:
    t[index] = index*l1.time_step;
    
    for ( int k=0; k<2; ++k ) {
      E[k][index]=l1.E(k, t[index]);
    }

    ++index;
  }

  // Integrate E to get A
  for ( int k=0; k<2; ++k ) {
    Integrate_E_trapezoid_rule(E[k],t,A[k]);
  }

  // Integrate A to get Intregral_A
  for ( int k=0; k<2; ++k ) {
    Integrate_A_trapezoid_rule(A[k],t,Integral_A[k]);
  }

  for ( int k=0; k<2; ++k ) {
    Integrate_A2_trapezoid_rule(A[k],t,Integral_A2[k]);
  }


};

Field::Field(Laser &l1, Laser &l2, double delay, double initial_Ax = 0, double initial_Ay = 0)
// for shorthand notation Ex := E[0] :
: Ex(E[0]), Ey(E[1]), Ax(A[0]), Ay(A[1]),
  Integral_Ax(Integral_A[0]), Integral_Ay(Integral_A[1]),Integral_Ax2(Integral_A2[0]), Integral_Ay2(Integral_A2[1])
{
  linear_polarized = ( l1.ellipticity()==0 && l2.ellipticity()==0 &&
                       l1.pol_angle() ==0 && l2.pol_angle() ==0 );
  double min_step=std::min(l1.time_step, l2.time_step);
  double max_step=std::max(l1.time_step, l2.time_step);


  size_t index=0;         // the time index
  Ax.push_back(initial_Ax);
  Ay.push_back(initial_Ay);

  Integral_Ax.push_back(0);    // HACK: Integral_A(t_f) := 0; But how else??
  Integral_Ay.push_back(0);    // HACK: Integral_A(t_f) := 0;

  Integral_Ax2.push_back(0);    // HACK: Integral_A(t_f) := 0; But how else??
  Integral_Ay2.push_back(0);    // HACK: Integral_A(t_f) := 0;


  size_t steps;
  Laser *l;

  double center12_offset=(l1.laser_duration-l2.laser_duration)/2;
  double l1_start=0;
  double l2_start=center12_offset + delay;
  // phase I: one laser has startet, the other not
  // decide which one:
  if (l2_start > 0) { // l2 did not start yet
    l=&l1;
  } else {
    l=&l2;
    l2_start=0;
    l1_start=-center12_offset - delay;
  }
  double l1_end=l1_start + l1.laser_duration;
  double l2_end=l2_start + l2.laser_duration;

  double time_step=l->time_step;

  //double first_starts=std::min( l1_start, l2_start);
  double first_starts = 0;      // decision: first always starts at t=0;
  double time=first_starts;     // start when the first starts
  double second_starts=std::max( l1_start, l2_start);
  double first_ends=std::min( l1_end, l2_end);
  double second_ends=std::max( l1_end, l2_end);

  double endPhase1=std::min( second_starts, first_ends );
  //cout<<"l1_start: "<<l1_start<<endl;
  //cout<<"l2_start: "<<l2_start<<endl;

  if (second_starts > first_starts) { // one did not start yet
  steps = ( size_t ) std::ceil( (endPhase1-first_starts)/time_step );
  //cout<<"1: "<<endPhase1<<" - "<<first_starts<<endl;
  for(size_t i=0; i < steps ; ++i)
  {
    // the simulation time:
    t.push_back(time);
    
    for ( int k=0; k<2; ++k ) {
      double tmp=l->E(k, time); E[k].push_back(tmp);
    }

    time+=time_step; ++index;
  }
  }

  // phase II: both lasers: either on or off
  double endPhase2=std::max(second_starts, first_ends);
  bool both_on=( (first_ends-second_starts) > 0);

  time=endPhase1-first_starts;
  // take the smallest time_step of both lasers if both are on
  time_step=both_on ? min_step : max_step;   // min step if both_on
  steps = ( size_t ) std::ceil( (endPhase2 - endPhase1 )/time_step );

  //cout<<"2: "<<endPhase2<<" - "<<endPhase1<<endl;
  for(size_t i=0; i < steps ; ++i)
  {
    // the simulation time:
    t.push_back(time);
    
    for ( int k=0; k<2; ++k ) {
      if(both_on) {
        double tmp=l1.E(k, time-l1_start) + l2.E(k, time-l2_start);
        E[k].push_back(tmp);
      } else {
        E[k].push_back(0);                        // no laser is on
      }

    }   // loop over coordinates
    
    time+=time_step; ++index;
  }

  // phase III: one laser already off
  time=endPhase2-first_starts;
  // decide which one:
  double start_of_second_ender;
  if( l1_end > l2_end) {
    l=&l1;
    start_of_second_ender=l1_start;
  } else {
    l=&l2;
    start_of_second_ender=l2_start;
  }
  // get the time_step of the laser that is on now
  time_step=l->time_step;

  if (second_ends>endPhase2) {
    steps = ( size_t ) std::ceil( (second_ends - endPhase2 )/time_step );
  //cout<<"3: "<<second_ends<<" - "<<endPhase2<<endl;
    for(size_t i=0; i < steps ; ++i)
    {
      // the simulation time:
      t.push_back(time);
      
      for ( int k=0; k<2; ++k ) {
        double tmp=l->E(k, time-start_of_second_ender);
        E[k].push_back(tmp);
      }
      
      time+=time_step; ++index;
    }
  }
  // final time:
  //time=second_ends-first_starts;

  // Integrate E to get A
  for ( int k=0; k<2; ++k ) {
    Integrate_E_trapezoid_rule(E[k],t,A[k]);
  }

  // Integrate A to get Intregral_A
  for ( int k=0; k<2; ++k ) {
    Integrate_A_trapezoid_rule(A[k],t,Integral_A[k]);
  }

  for ( int k=0; k<2; ++k ) {
    Integrate_A2_trapezoid_rule(A[k],t,Integral_A2[k]);
  }


};

inline void build_derivs(vector<double> &t_data, vector<double> &F, vector<double> &dF) {

  const size_t points=F.size();
  if(t_data.size() != points ) {
    cerr<<"Field/build_derivs: time-vector size is different from Field-vector!"<<endl;
    cerr<<"t_data.size()="<<t_data.size()<<"\n";
    cerr<<"F.size()="<<F.size()<<"\n";
    exit(1);
  };

  if(dF.size() != points ) dF.resize(points );

  // build table of derivatives from F[i] and t[i]
  for (size_t i=0; i < points; ++i) {
    if(0==i) {
      dF[i]=(F[i+1]-F[i])/(t_data[i+1]-t_data[i]);
    } else if(points==i) {
      dF[i]=(F[i]-F[i-1])/(t_data[i]-t_data[i-1]);
    } else {
      double left =(F[i]-F[i-1])/(t_data[i]-t_data[i-1]);
      double right=(F[i+1]-F[i])/(t_data[i+1]-t_data[i]);
      dF[i]=(left+right)/2;
    }
  }

};

/// Spline code from Przemek
//
inline double hermite(const double t, const double y1, const double d1, const double y2, const double d2){
  // result=y1 for t=0;
  // result=y2 for t=1;
  double t2, t3, result;
  t2=t*t;
  t3=t2*t;
  result=(2*t3 - 3*t2 +1)*y1 + (t3 - 2*t2 +t)*d1 ;
  result+=  (-2*t3 + 3*t2)*y2 + (t3-t2)*d2;
  return result;
}

inline void build_spline(vector<double> &t_data, vector<double> &F, vector<double> &dF, double dt, vector<double> &R) {

  const size_t points=F.size();
  if(t_data.size() != points || dF.size() != points ) {
    cerr<<"Field/build_derivs: time-vector size is different from Field-vector!"<<endl;
    cerr<<"t_data.size()="<<t_data.size()<<"\n";
    cerr<<"F.size()="<<F.size()<<"\n";
    cerr<<"dF.size()="<<dF.size()<<"\n";
    exit(1);
  };

  const size_t steps=(size_t) round((t_data.back()-t_data[0])/dt) + 1;
  R.resize(steps );
  R[0]=F[0];

  double time=t_data[0]+dt; size_t jR=1;
  for (size_t i=1; i < steps; ++i) {

    while (time > t_data[jR]) { 
      // make sure we don't fall off the edge of data!
      if ( jR < points ) {
        ++jR;
      } else break;     // in case time is an epsilon over t[points-1]
    }

    double scale=(time - t_data[jR-1])/(t_data[jR] - t_data[jR-1]);
    R[i]=hermite(scale, F[jR-1], dF[jR-1], F[jR], dF[jR] );
    time+=dt;
  }

};

Field::Field(string filename, double dt, double scaling = 1, double initial_Ax = 0, double initial_Ay = 0)
// for shorthand notation Ex := E[0] :
: Ex(E[0]), Ey(E[1]), Ax(A[0]), Ay(A[1]),
  Integral_Ax(Integral_A[0]), Integral_Ay(Integral_A[1]), Integral_Ax2(Integral_A2[0]), Integral_Ay2(Integral_A2[1])
{
  /* Step (1): read in Data */
  vector<double> t_data, v[2];

  linear_polarized=false;             // assumption
  string tag("Field(\""+filename+"\"): ");
  std::ifstream ifs(filename.c_str());
  if(ifs.fail()) {
    cerr<<tag<<"Could not open file"<<"\n";
    exit (1);
  }
  skip_comments(ifs);     // skip_comments defined in utils.h

  size_t lines_read=0;

  string linebuffer;
  double _t, _vx, _vy;
  while(!ifs.fail()) {

    std::getline(ifs, linebuffer);
    if (ifs.eof()) {  // normal file end
      break;
    }
    if(linebuffer.empty()) {
      cerr<<tag<<"skipping empty line\n";
      continue;
    }
    std::stringstream parser(linebuffer);

    if (0==lines_read) {   // check if we have two or three columns
      parser>>_t>>_vx>>_vy;      // try to read three columns
      if (parser.fail()) {       // if this fails
        std::stringstream parser2(linebuffer);
        parser2>>_t>>_vx;        // try to read just two columns
        if (!parser2.fail()) {   // great, worked
          linear_polarized=true; // adjust assumption
          parser.clear();        // clear failure state of parser
        }
      }
    } else {

      parser>>_t>>_vx;
      if(!linear_polarized)
        parser>>_vy;
      else _vy=0;

    }
    if (parser.fail()){ // tedious error checking
      cerr<<tag<<"could not read data line "<<t.size()<<"\n";
      if(t_data.size()!=0) {
        cerr<<"     continuing with the data read so far\n";
        break;
      } else exit(1);
    }
    t_data.push_back(_t);
    v[0].push_back(_vx*scaling);
    v[1].push_back(_vy*scaling);
    lines_read++;
  }
  if (0==lines_read) {
    cerr<<tag<<"Could not read any data from file\n";
    exit(1);
  }

  if(!linear_polarized) { // print out a message
    cout<<tag<<"non-linear polarization"<<"\n";
  } else {
    cout<<tag<<"linear polarization"<<"\n";
  }

  /* Step (2): Interpolate to match demanded dt */
  /* Step (2.1): build up array of derivatives */

  // vectors for derivatives, will get resized by build_derivs:
  vector<double> dEdt[2];
  build_derivs(t_data, v[0], dEdt[0]);
  build_derivs(t_data, v[1], dEdt[1]);

  /* Step (2.0): build up new t-axis based on dt */
  const size_t steps=(size_t) round((t_data.back()-t_data[0])/dt) + 1;
  t.resize(steps );
  t[0]=t_data[0];
  for (size_t i=1; i < steps; ++i) {
    t[i]=t[i-1]+dt;
  }
  /* Step (2.1): build up new E based on spline */
  if (steps != 0) {
  build_spline(t_data, v[0], dEdt[0], dt, E[0]);
  Ax.push_back(initial_Ax);
  Integral_Ax.push_back(0);    // HACK: Integral_A(t_f) := 0; But how else??
  Integral_Ax2.push_back(0);    // HACK: Integral_A(t_f) := 0; But how else??

  build_spline(t_data, v[1], dEdt[1], dt, E[1]);
  Ay.push_back(initial_Ay);
  Integral_Ay.push_back(0);    // HACK: Integral_A(t_f) := 0;
  Integral_Ay2.push_back(0);    // HACK: Integral_A(t_f) := 0;
  }

  /* Step (3): Integrate E to get A */
  // Ex = E[0] and Ey = E[1]
  for ( int k=0; k<2; ++k ) {
     Integrate_E_trapezoid_rule(E[k],t,A[k]);
  };

  // Integrate A to get Intregral_A
  for ( int k=0; k<2; ++k ) {
    Integrate_A_trapezoid_rule(A[k],t,Integral_A[k]);
  }

  for ( int k=0; k<2; ++k ) {
    Integrate_A2_trapezoid_rule(A[k],t,Integral_A2[k]);
  }


};

