/*
 * fluid.cc
 * 
 * Fluid simulator using Lennard-Jones potential and Verlet algorithm. 
 * 
 * Copyright (c) 07/2004, 08/2005 by Wolfgang Wieser, 
 *      email: > wwieser -a- gmx -*- de <
 * 
 * Version 0.9b (2006-03-26)
 * 
 * The above email is for bug reports and improvements of the code, ONLY. 
 * It is not for general discussion or special questions about the code or 
 * about what is being done. 
 * Use plaintext emails; HTML mail may be considered as spam. 
 * 
 * This file may be distributed and/or modified under the terms of the 
 * GNU General Public License version 2 as published by the Free Software 
 * Foundation. 
 * 
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 * 
 */

//@@--------------------------------------------------------------------
// NOTE: In order to set simulation parameters, some lines in the 
//       source code need to be edited. 
//       The primary user-changeable parameters are set in lines 
//       which are enclosed by a comment beginning with "//@@----" just 
//       like these lines you are currently reading. Just search through 
//       the code, tune some parameters and then run "make run". 
//@@--------------------------------------------------------------------

#include <stdio.h>
#include "alloc.h"
#include <math.h>

// NOTE: If you do not have sys/time.h and/or gettimeofday(), remove this 
//       #include statement and eliminate the runtime measurement which is 
//       using gettimeofday() by commenting out the corresponding lines. 
//       It is not needed for the sim to run properly; just "eye candy". 
#include <sys/time.h>

#include <assert.h>

// Setting this to 1 will additionally use integer arithmetic for the 
// initial distance check. Expect little speed increase (maybe better for 
// larger box sizes and small cutoff radi). 
#define USE_INT_ARITH 0


// Calculation in reduced sizes: 
//  [ Lennard-Jones potential: V(r) = V0 * (r/sigma)^-12 - (r/sigma)^-6 ]
// --------------------------------------------
//  m     -> 1
//  x,y,z -> (x,y,z) / sigma
//  t     -> t * sqrt(eps/(m*sigma^2))
//  rho   -> rho * sigma^3
//  v     -> v * sqrt(m/epsilon)   (speed)
//  F     -> F * sigma/epsilon
//  p     -> p * sigma^3/epsilon   (pressure)
//  E     -> E / epsilon
//  T     -> T * k_B/epsilon
//  V0    -> 4 
// --------------------------------------------

template<typename T>inline T SQR(T x)
	{  return(x*x);  }
template<typename T>inline T POW3(T x)
	{  return(x*x*x);  }

template<typename T>inline T VLENGTH3(const T *x)
	{  return(sqrt(SQR(x[0])+SQR(x[1])+SQR(x[2])));  }

#undef MIN
template<typename T>inline T MIN(T a,T b)
	{  return(a<b ? a : b);  }
inline double MINabs(double a,double b)
	{  return(fabs(a)<fabs(b) ? a : b);  }
inline int MINabs(int a,int b)
	{  return(abs(a)<abs(b) ? a : b);  }

// Random number in range -1..1: 
inline double rand11()
{
	int r=rand();
	double dr=double(r)/RAND_MAX;
	return((dr+dr)-1.0);
}


#if USE_INT_ARITH
typedef int intpos_t;
#endif


class FluidSimulator
{
	public:
		// Representing a particle: 
		struct Particle
		{
			double r[3];   // location
			double v[3];   // speed
#if USE_INT_ARITH
			intpos_t ir[3];  // integer grid coordinates (approx for r)
#endif
		};
		
		// Represent acceleration: 
		// This is not stored with the particle as it is not part of the 
		// state of the particle. 
		struct Accel
		{
			double a[3];
		};
		
	private:
		// Particle array and number of particles. 
		Particle *P,*Pend;   // Pend=P+np
		int np;
		// Current accelerations: 
		Accel *A;
		
		// Simulation time step:
		double dt;
		
		// Simulation box size: 
		double L;
		
		// Force radius: compute PP force for all particles which are 
		// less distant then this: 
		double cutoff_r;
		
		// Time average of the PCF: 
		int npcf;
		double pcf_rmax;
		double *pcf_sum;
		double *pcf_sum2;
		int pcf_samples;   // number of samlples accumulated so far
		
#if USE_INT_ARITH
		// Integer coos scale factors: 
		double r2ir;
		// Integer corrs are in range 0..irmax-1
		intpos_t irmax;
#endif
		
	public:
		// Some simulation statistics: 
		int nskip0;
		int nskip1;
		int ncalc;
		int nborder;
		double r_min;
		
	private:
		// Calculate the potential energy. 
		// This is a Lennard-Jones potential: 
		// Argument: 1/r^2. 
		static inline double potential_1r2(double r2)
		{
			double r6=r2*r2*r2;
			return(4.0*r6*(r6-1.0));
		}
		// Argument: r
		// Note: The min of the potential is at r=1.1225. 
		static inline double potential(double r)
		{  return(potential_1r2(1.0/(r*r)));  }
		
		// Calc 1/|r| * \partial V(|r|) / \partial |r|: 
		// Argument: 1/r^2. 
		static inline double potentialD_1r2(double r2)
		{
			double r6=r2*r2*r2;
			return(48.0*r6*r2*(0.5-r6));
		}
		
		// Calc r+dv/dr = -24*(2*r^-12-r^-6). Argument r. 
		static inline double potential_w(double r)
		{
			double r6=1.0/POW3(SQR(r));  // r^-6
			return(24*r6*(1-2*r6));
		}
		
#if USE_INT_ARITH
		// Computer integer grid coos: 
		inline void _SetIR(Particle *p)
		{
			p->ir[0]=intpos_t(p->r[0]*r2ir);
			p->ir[1]=intpos_t(p->r[1]*r2ir);
			p->ir[2]=intpos_t(p->r[2]*r2ir);
		}
#endif
		
		// Calculate the current accelerations: 
		void _CalcAccel();
		
	public:
		FluidSimulator();
		~FluidSimulator();
		
		inline int NP() const
			{  return(np);  }
		
		inline int get_npcf() const
			{  return(npcf);  }
		
		inline double get_pcf_rmax() const
			{  return(pcf_rmax);  }
		
		inline double get_rho() const
			{  return(np/POW3(L));  }
		
		inline double get_L() const
			{  return(L);  }
		
		// Clear state, delete all particles. 
		void Clear();
		
		// Initialize fluid: 
		// Pass number of parameters, simulation time step, 
		// particle density, force cutoff
		void Initialize(int n,double _dt,double _rho,double _cutoff_r,
			double _initial_speed);
		// Initialize time-accumulation of PCF: 
		void InitializePCF(int _npcf,double _rmax);
		
		// Compute a simulation step: 
		void SimulationStep();
		
		// Compute total kinetic energy (returned) and impulse (stored in 
		// ptot[]). 
		double CalcEnergyAndImpulse(double *ptot);
		
		// Calculate the temperature from the kinetic energy. 
		// This makes use of E_kin = 3/2 * N * k_B * T. 
		double Ekin2Temp(double E_kin)
			{  return(E_kin/(1.5*np));  }
		
		// Compute the pair correlation function: 
		// Returns discrete PCF in pcf[] (if non-NULL) which has n entries 
		// for the distances 0..rmax. 
		// If sum_pcf, sum_pcf2 are non-NULL, _add_ the current PCF, 
		// to *sum_pcf and add the squares to *sum_pcf2. This is used by 
		// AccumulatePCF() to accumulate the PCF and the sigma estimates. 
		void ComputePCF(double *pcf,int n,double rmax,double *sum_pcf=NULL,
			double *sum_pcf2=NULL);
		
		// Accumulate PCF for time averaging. 
		// Internally calls ComputePCF(). 
		// Will also store the current PCF in curr_pcf if non-NULL. 
		void AccumulatePCF(double *curr_pcf=NULL);
		
		// Compute averaged PCF (accumulated data by AccumulatePCF()) 
		// and return PCF with sigma in pcf,pcf_s. Max sigma value is 
		// returned as retval. 
		double AveragedPCF(double *pcf,double *pcf_s);
		
		// Compute complete potential energy using PCF. 
		double ComputeEpot(const double *pcf,int n,double rmax);
		
		// Compute the pressure using PCF and current temperature. 
		double ComputePressure(const double *pcf,int n,double rmax,
			double T_curr);
		
		// Plot an X-Y dump of all particles; mainly useful for debugging 
		// reasons. 
		void PlotAllParticles();
};


void FluidSimulator::Clear()
{
	P=FREE(P);
	Pend=NULL;
	np=0;
	A=FREE(A);
	
	npcf=0;
	pcf_sum=FREE(pcf_sum);
	pcf_sum2=FREE(pcf_sum2);
	pcf_samples=0;
}


void FluidSimulator::InitializePCF(int _npcf,double _rmax)
{
	pcf_sum=FREE(pcf_sum);
	pcf_sum2=FREE(pcf_sum2);
	
	npcf=_npcf;
	pcf_rmax=_rmax;
	pcf_sum=ALLOC<double>(npcf);
	pcf_sum2=ALLOC<double>(npcf);
	pcf_samples=0;
	
	for(int i=0; i<npcf; i++)
	{
		pcf_sum[i]=0.0;
		pcf_sum2[i]=0.0;
	}
}

void FluidSimulator::Initialize(int n,double _dt,double _rho,double _cutoff_r,
	double _initial_speed)
{
	Clear();
	np=n;
	P=ALLOC<Particle>(np);
	Pend=P+np;
	A=ALLOC<Accel>(np);
	
	dt=_dt;
	L=pow(np/_rho,1.0/3.0);
	cutoff_r=_cutoff_r;
	
#if USE_INT_ARITH
	irmax=1<<13;  // make sure irmax*irmax < max value
	r2ir=irmax/L;
#endif
	
	// We have n particles in a box of size L x L x L. 
	// To distribute them regularly on a rectangular grid, this grid 
	// needs to have the size gsize: 
	int gsize=(int)(ceil(pow(double(n),1.0/3.0))+0.5);
	int gsize3=gsize*gsize*gsize;
	assert(n<=gsize3);
	
	// Particle displacelemnt as multiple of initial placement grid size. 
	// Never use values >0.5. 
	double displacement_fact=0.1;
	
	// Initial speed in x,y,z direction is in this range (for each 
	// coordinate): 
	double initial_speed=_initial_speed;
	
	// Precompute...
	double L_gsize=L/gsize;
	fprintf(stderr,"Initial grid size: %g (%d) (L=%g)\n",
		L_gsize,gsize,L);
	assert(L_gsize*sqrt(3)>1.225);  // -> See the min of the potential. 
	
	// Set initial locations and velocities: 
	// Also sum up total impulse. 
	double ptot[3]={0.0,0.0,0.0};
	for(Particle *p=P; p<Pend; p++)
	{
		// The particles are put on a grid: 
		// Make sure they are distributed regularly: 
		int grid_idx=int( double(gsize3)*(p-P)/n +0.5);
		int gx = grid_idx % gsize;  grid_idx/=gsize;
		int gy = grid_idx % gsize;  grid_idx/=gsize;
		int gz = grid_idx;  assert(grid_idx<gsize);
		p->r[0] = L_gsize*( gx+0.5 + rand11()*displacement_fact );
		p->r[1] = L_gsize*( gy+0.5 + rand11()*displacement_fact );
		p->r[2] = L_gsize*( gz+0.5 + rand11()*displacement_fact );
		
#if USE_INT_ARITH
		_SetIR(p);
#endif
		
		// Initial velocity: 
		p->v[0]= rand11()*initial_speed;
		p->v[1]= rand11()*initial_speed;
		p->v[2]= rand11()*initial_speed;
		
		// Sum up impulse: 
		ptot[0]+=p->v[0];
		ptot[1]+=p->v[1];
		ptot[2]+=p->v[2];
	}
	
	// Change speed for total impulse 0: 
	ptot[0]/=np;
	ptot[1]/=np;
	ptot[2]/=np;
	for(Particle *p=P; p<Pend; p++)
	{
		p->v[0]-=ptot[0];
		p->v[1]-=ptot[1];
		p->v[2]-=ptot[2];
	}
	
#if 0
	// Testing purposes only: 
	P[0].r[0]=L-1;
	P[0].r[1]=
	P[0].r[2]=0.5*L;
	P[1].r[0]=0+1;
	P[1].r[1]=
	P[1].r[2]=0.5*L;
	
	P[0].v[0]=
	P[0].v[1]=
	P[0].v[2]=
	P[1].v[0]=
	P[1].v[1]=
	P[1].v[2]=0.0;
#endif
	
	// Calculate kinetic energy and accumulated impulse (as a check): 
	double Ekin=CalcEnergyAndImpulse(ptot);
	fprintf(stderr,"Initially: E_kin=%g, |P|=%.2g, T=%g\n",
		Ekin,VLENGTH3(ptot),Ekin2Temp(Ekin));
	
	// Compute initial values for the acceleration: 
	_CalcAccel();
}


double FluidSimulator::CalcEnergyAndImpulse(double *ptot)
{
	ptot[0]=0.0;
	ptot[1]=0.0;
	ptot[2]=0.0;
	
	double Ekin=0.0;
	for(Particle *p=P; p<Pend; p++)
	{
		ptot[0]+=p->v[0];
		ptot[1]+=p->v[1];
		ptot[2]+=p->v[2];
		
		Ekin += SQR(p->v[0])+SQR(p->v[1])+SQR(p->v[2]);
	}
	Ekin*=0.5;
	
	return(Ekin);
}


void FluidSimulator::_CalcAccel()
{
	for(Accel *a=A,*Aend=A+np; a<Aend; a++)
	{
		a->a[0]=0.0;
		a->a[1]=0.0;
		a->a[2]=0.0;
	}
	
	double Lcutoff_r=L-cutoff_r;
	double cutoff_r2=cutoff_r*cutoff_r;
	double L2=0.5*L;
	
#if USE_INT_ARITH
	intpos_t iL=irmax;  //intpos_t(L*r2ir);
	intpos_t iL2=irmax/2;  //intpos_t(L2*r2ir);
	intpos_t icutoff_r=intpos_t(cutoff_r*r2ir)+3;
	intpos_t iLcutoff_r=iL-icutoff_r;
	intpos_t icutoff_r2=intpos_t(cutoff_r2*r2ir*r2ir)+16;
	//fprintf(stderr,">%d<\n",icutoff_r2);
#endif
	
	nskip0=0;
	nskip1=0;
	ncalc=0;
	nborder=0;
	r_min=1e30;
	
	int i=0;
	double wp[3];
#if USE_INT_ARITH
	intpos_t iwp[3];
#endif
	for(Particle *p=P; p<Pend; p++,i++)
	{
		// Calculate the current acceleration for particle *p: 
		
		// We have periodic boundary conditions here. 
		// In case the particle is enough far away from the box 
		// borders, there is no need for special treatment. 
#if USE_INT_ARITH
		bool on_border=0;
		if(p->ir[0]<icutoff_r || p->ir[1]<icutoff_r || p->ir[2]<icutoff_r || 
		   p->ir[0]>=iLcutoff_r || p->ir[1]>=iLcutoff_r || p->ir[2]>=iLcutoff_r)
		{
			on_border=
				!(p->r[0]>=cutoff_r && p->r[1]>=cutoff_r && p->r[2]>=cutoff_r && 
			      p->r[0]<Lcutoff_r && p->r[1]<Lcutoff_r && p->r[2]<Lcutoff_r );
		}
#else
		bool on_border=
			!(p->r[0]>=cutoff_r && p->r[1]>=cutoff_r && p->r[2]>=cutoff_r && 
		      p->r[0]<Lcutoff_r && p->r[1]<Lcutoff_r && p->r[2]<Lcutoff_r );
#endif
		if(on_border)
		{
			// Need special treatment for periodic BC: wrap around edge. 
			// Particle position wrapped around edge: 
			wp[0] = p->r[0]<L2 ? p->r[0]+L : p->r[0]-L;
			wp[1] = p->r[1]<L2 ? p->r[1]+L : p->r[1]-L;
			wp[2] = p->r[2]<L2 ? p->r[2]+L : p->r[2]-L;
#if USE_INT_ARITH
			iwp[0] = p->ir[0]<iL2 ? p->ir[0]+iL : p->ir[0]-iL;
			iwp[1] = p->ir[1]<iL2 ? p->ir[1]+iL : p->ir[1]-iL;
			iwp[2] = p->ir[2]<iL2 ? p->ir[2]+iL : p->ir[2]-iL;
#endif
			++nborder;
		}
		
		int j=i+1;
		for(Particle *k=p+1; k<Pend; k++,j++)
		{
#if USE_INT_ARITH
			intpos_t idr[3]={
				k->ir[0]-p->ir[0],
				k->ir[1]-p->ir[1],
				k->ir[2]-p->ir[2] };
			if(on_border)
			{
				idr[0] = MINabs( idr[0], k->ir[0]-iwp[0] );
				idr[1] = MINabs( idr[1], k->ir[1]-iwp[1] );
				idr[2] = MINabs( idr[2], k->ir[2]-iwp[2] );
			}
			
			//if(idr[0]>=icutoff_r || idr[1]>=icutoff_r || idr[2]>=icutoff_r) 
			//{  ++nskip0;  continue;  }
			int ir2=SQR(idr[0])+SQR(idr[1])+SQR(idr[2]);
			if(ir2>=icutoff_r2)
			{  ++nskip0;  continue;  }
#endif
			
			double dr[3]={
				k->r[0]-p->r[0],
				k->r[1]-p->r[1],
				k->r[2]-p->r[2] };
			if(on_border)
			{
				dr[0] = MINabs( dr[0], k->r[0]-wp[0] );
				dr[1] = MINabs( dr[1], k->r[1]-wp[1] );
				dr[2] = MINabs( dr[2], k->r[2]-wp[2] );
			}
			
			//if(dr[0]>=cutoff_r || dr[1]>=cutoff_r || dr[2]>=cutoff_r) 
			//{  ++nskip0;  continue;  }  // <-- Leave away for speed!
			
			double r2=SQR(dr[0])+SQR(dr[1])+SQR(dr[2]);
			if(r2>=cutoff_r2)
			{  ++nskip1;  continue;  }
			//if(r2<0.7*0.7)  fprintf(stderr,"COLLIDE: %g\n",r2);
			double dr_V_r=potentialD_1r2(1.0/r2);  // (\partial_r V) / r
			double r=sqrt(r2);
			
			double tmp;
			tmp=dr_V_r*dr[0];  A[i].a[0]+=tmp;  A[j].a[0]-=tmp;
			tmp=dr_V_r*dr[1];  A[i].a[1]+=tmp;  A[j].a[1]-=tmp;
			tmp=dr_V_r*dr[2];  A[i].a[2]+=tmp;  A[j].a[2]-=tmp;
			
			if(r_min>r)  r_min=r;
			++ncalc;
//fprintf(stderr,"%d <-> %d: %g\n",i,j,dr_V_r);
		}
	}
}


void FluidSimulator::ComputePCF(double *pcf,int n,double rmax,
	double *sum_pcf,double *sum_pcf2)
{
	int cnt[n+1];
	for(int i=0; i<=n; i++)
		cnt[i]=0;
	
	double Lcutoff_r=L-cutoff_r;
	double L2=0.5*L;
	
	double wp[3];
	double rmax2=SQR(rmax*(n+2)/n);
	for(Particle *p=P; p<Pend; p++)
	{
		bool on_border=
			!(p->r[0]>=cutoff_r && p->r[1]>=cutoff_r && p->r[2]>=cutoff_r && 
		      p->r[0]<Lcutoff_r && p->r[1]<Lcutoff_r && p->r[2]<Lcutoff_r );
		if(on_border)
		{
			// Need special treatment for periodic BC: wrap around edge. 
			// Particle position wrapped around edge: 
			wp[0] = p->r[0]<L2 ? p->r[0]+L : p->r[0]-L;
			wp[1] = p->r[1]<L2 ? p->r[1]+L : p->r[1]-L;
			wp[2] = p->r[2]<L2 ? p->r[2]+L : p->r[2]-L;
		}
		for(Particle *k=p+1; k<Pend; k++)
		{
			double dr[3]={
				k->r[0]-p->r[0],
				k->r[1]-p->r[1],
				k->r[2]-p->r[2] };
			if(on_border)
			{
				dr[0] = MINabs( dr[0], k->r[0]-wp[0] );
				dr[1] = MINabs( dr[1], k->r[1]-wp[1] );
				dr[2] = MINabs( dr[2], k->r[2]-wp[2] );
			}
			double r=SQR(dr[0])+SQR(dr[1])+SQR(dr[2]);
			if(r>rmax2)  continue;
			r=sqrt(r);
			int idx=int(n*r/rmax);  // truncate to int
			if(idx>n)  continue;
			++cnt[idx];
		}
	}
	
	double rho = get_rho();
	double N0_delta_fact = 2.0*M_PI/3.0*np*rho;
	for(int i=0; i<n; i++)
	{
		double r0=rmax*(i)/n;
		double r1=rmax*(i+1)/n;
		double PCF = double(cnt[i]) / ( N0_delta_fact*(POW3(r1)-POW3(r0)) );
		if(pcf)       pcf[i]=PCF;
		if(sum_pcf)   sum_pcf[i]+=PCF;
		if(sum_pcf2)  sum_pcf2[i]+=SQR(PCF);
	}
}


void FluidSimulator::AccumulatePCF(double *curr_pcf)
{
	ComputePCF(curr_pcf,npcf,pcf_rmax,pcf_sum,pcf_sum2);
	++pcf_samples;
}


double FluidSimulator::AveragedPCF(double *pcf,double *pcf_s)
{
	double max_s=0.0;
	for(int i=0; i<npcf; i++)
	{
		double p = pcf_sum[i]/pcf_samples;
		double p2 = pcf_sum2[i]/pcf_samples;
		pcf[i]=p;
		pcf_s[i]=sqrt((p2-SQR(p))/pcf_samples);
		if(max_s<pcf_s[i])
		{  max_s=pcf_s[i];  }
	}
	return(max_s);
}


double FluidSimulator::ComputeEpot(const double *pcf,int n,double rmax)
{
	// Basically, this does the numerical integration 
	// 0.5*N*rho * integral 0..rmax potential(r)*pcf(r)*4*pi*r^2 dr
	
	double sum=0.0;
	for(int i=0; i<n; i++)
	{
		double ri=rmax*(i+0.5)/n;
		sum += potential(ri)*pcf[i]*SQR(ri);
	}
	// Multiply with 0.5*N*rho * 4*pi * dr. 
	sum *= (2*np)*M_PI*get_rho()*rmax/n;
	
	// Add analytic correction due to finite range of PCF. 
	double rm_3=1.0/POW3(rmax);
	sum += 8.0/3.0*M_PI*np*get_rho()*rm_3*(SQR(rm_3)/3-1);
	
	return(sum);
}


double FluidSimulator::ComputePressure(const double *pcf,int n,double rmax,
	double T_curr)
{
	// First, do the numerical integration 
	// -N*rho/6 * integral 0..rmax potential_w(r)*pcf(r)*4*pi*r^2 dr
	// with potential_w = r * \partial potential(r) / \partial r
	
	double sum=0.0;
	for(int i=0; i<n; i++)
	{
		double ri=rmax*(i+0.5)/n;
		sum += potential_w(ri)*pcf[i]*SQR(ri);
	}
	// Multiply with N*rho/6 * 4*pi * dr. 
	sum *= -np*M_PI/1.5*get_rho()*rmax/n;
	
	// Add N * k_b * T: 
	sum += np*T_curr;
	
	// Add analytic correction due to finite range of PCF. 
	double rm_3=1.0/POW3(rmax);
	sum += 16.0/3.0*M_PI*np*SQR(get_rho())*rm_3*(SQR(rm_3)/1.5-1);
	
	// Now, sum=p*V. 
	return(sum/POW3(L));
}


void FluidSimulator::SimulationStep()
{
	// Calculate a simulation step. 
	double dt2=0.5*dt;
	Accel *a=A;

//@@--------------------------------------------------------------------
// User tunable parameter: 
// Set a cooling factor to damp velocities and hence take kinetic 
// energy out of the system. Use the empty define if you wish to 
// temperature to stay constant (normally recommended). 
// Usable cooling factors are in range 0.99 to 0.9999. 
// Do not use negative values or values larger than 1. 
//#define COOL_FACT *0.997
#define COOL_FACT
//@@--------------------------------------------------------------------

	for(Particle *p=P; p<Pend; p++,a++)
	{
		for(int i=0; i<3; i++)
		{
			double adt2 = a->a[i]*dt2;
			p->r[i] += p->v[i]*dt + adt2*dt;
assert(finite(p->r[i]));
if(p->r[i]<-L || p->r[i]>L+L)  fprintf(stderr,"OOPS: p[%d].r[%d]=%g\n",p-P,i,p->r[i]);
			     if(p->r[i]<0.0)   p->r[i]=L-fmod(-p->r[i],L);
			else if(p->r[i]>=L)    p->r[i]=fmod(p->r[i],L);
assert(p->r[i]>=0.0 && p->r[i]<=L);
			p->v[i] = (p->v[i] + adt2) COOL_FACT;
		}
#if USE_INT_ARITH
		_SetIR(p);
#endif
	}
	
	_CalcAccel();
	
	a=A;
	for(Particle *p=P; p<Pend; p++,a++)
	{	
		p->v[0] = (p->v[0] + a->a[0]*dt2) COOL_FACT;
		p->v[1] = (p->v[1] + a->a[1]*dt2) COOL_FACT;
		p->v[2] = (p->v[2] + a->a[2]*dt2) COOL_FACT;
	}
}


// Plot mainly useful for debugging reasons: 
void FluidSimulator::PlotAllParticles()
{
	printf("set title \"Particle dump xy\"\n");
	//printf("set xlabel \"MCS\"\n");
	//printf("set ylabel \"M\"\n");
	printf(
		"plot [0:1] [0:1] '-' title \"particles (x,y)\", "
		"'-' title \"|v|/const\" with linespoints, "
		"'-' title \"PCF (avg)\" with errorbars, "
		"'-' title \"PCF (avg)\" with lines, "
		"'-' title \"PCF (curr)\" with linespoints\n");
	
	double vmax=1e-10;
	for(Particle *p=P; p<Pend; p++)
	{
		printf("%g %g\n",p->r[0]/L,p->r[1]/L);
		double v=sqrt(SQR(p->v[0])+SQR(p->v[1])+SQR(p->v[2]));
		if(vmax<v)  vmax=v;
	}
	printf("e\n");
	
	for(Particle *p=P; p<Pend; p++)
	{
		double v=sqrt(SQR(p->v[0])+SQR(p->v[1])+SQR(p->v[2]));
		printf("%g %g\n",(p-P+0.5)/np,0.3*v/vmax);
	}
	printf("e\n");
	
	/*
	const int npcf=50;
	double pcf[npcf];
	ComputePCF(pcf,npcf,3.0);
	double pcfmax=0.0;
	for(int i=0; i<npcf; i++)  if(pcfmax<pcf[i])  pcfmax=pcf[i];
	for(int i=0; i<npcf; i++)
	{
		printf("%g %g\n",(i+0.5)/npcf,0.3*pcf[i]/pcfmax+0.5);
	}
	printf("e\n");
	*/
	
	double pcfmax=0.0;
	for(int i=0; i<npcf; i++)  if(pcfmax<pcf_sum[i])  pcfmax=pcf_sum[i];
	pcfmax/=pcf_samples;
	double scale=0.3/pcfmax;
	
	for(int i=0; i<npcf; i++)
	{
		double p = pcf_sum[i]/pcf_samples;
		double p2 = pcf_sum2[i]/pcf_samples;
		printf("%g %g %g\n",
			(i+0.5)/npcf,scale*p+0.4,
			sqrt((p2-SQR(p))/pcf_samples));
	}
	printf("e\n");
	
	for(int i=0; i<npcf; i++)
	{
		double p = pcf_sum[i]/pcf_samples;
		printf("%g %g\n",
			(i+0.5)/npcf,scale*p+0.4);
	}
	printf("e\n");
	
	double pcf[npcf];
	ComputePCF(pcf,npcf,pcf_rmax);
	pcfmax=0.0;
	for(int i=0; i<npcf; i++)  if(pcfmax<pcf[i])  pcfmax=pcf[i];
	scale=0.3/pcfmax;
	
	for(int i=0; i<npcf; i++)
	{
		printf("%g %g\n",(i+0.5)/npcf,scale*pcf[i]+0.6);
	}
	printf("e\n");
}


FluidSimulator::FluidSimulator()
{
	P=Pend=NULL;
	np=0;
	A=NULL;
	
	npcf=0;
	pcf_sum=NULL;
	pcf_sum2=NULL;
}

FluidSimulator::~FluidSimulator()
{
	Clear();
}


void SimFluid()
{
	timeval start_tv;
	gettimeofday(&start_tv,NULL);
	
	FluidSimulator sim;
	
//@@--------------------------------------------------------------------
// Set the fluid simulation paramters: 
//   sim.Initialize(
//     Number of particles to simulate (N),
//     ODE integration time step (dt),
//     particle density (rho),
//     potential cutoff radius (cutoff_r), 
//     initial particle speed (initial_speed) );
//     
// Note that the box size is calculated automatically from the particle 
// density and the number of particles. 
// 
// Make sure that the cutoff diamenter is smaller than the box size. 
// Do not use density values rho>1.3 (or 1.5) or you force particles 
// into the highly repulsive range of the potential. 
// 
// Note that you need to decrease the time step (dt) if you increase 
// the initial particle speed. If you get lots of dumps like 
// "OOPS: p[190].r[2]=-33101.9", then abort the sim and reduce dt 
// significantly. 
// 
// Set the pair correlation function parameters: 
//   sim.InitializePCF(
//     number of PCF points (npcf),
//     max radius (rmax) );
// 
// The pair correlation function is computed by dividing the 
// particle distances into npcf distance slots from 0 to rmax 
// (i.e. each slot has size rmax/npcf). 
// Do not use too large values for npcf as the PCF is determined 
// by summing up over all particle pairs whose distance is less 
// than rmax (i.e. max N/2*(N-1) many). 
// Note also that is does not make sense to use rmax>L/2 when L 
// is the length of the box (L = (N/rho)^(1/3)). 
// 
	//sim.Initialize(108,/*dt=*/5e-3,/*rho=*/0.5,/*cutoff_r=*/2.5,
	//	/*initial_speed=*/1.0 /**2.780155*/);  // Dr. Stintzing's parameters
	//sim.InitializePCF(200,2.5); // 3.5 for 500, 3.0 for 108 particles
	
	//sim.Initialize(500,/*dt=*/5e-3,/*rho=*/1.3,/*cutoff_r=*/3.0,
	//	/*initial_speed=*/1.0);
	sim.Initialize(500,/*dt=*/5e-3,/*rho=*/1.2,/*cutoff_r=*/3.0,
		/*initial_speed=*/3.0);
	sim.InitializePCF(200,3.5);
	
	//sim.Initialize(1000,/*dt=*/0.5e-3,/*rho=*/0.95,/*cutoff_r=*/3.5,
	//	/*initial_speed=*/0.8);  // 1.0
	//sim.InitializePCF(300,3.5);
	
	//sim.Initialize(1000,/*dt=*/1e-5,/*rho=*/0.95,/*cutoff_r=*/3.5,
	//	/*initial_speed=*/1000);
	//sim.InitializePCF(200,3.5);
	
	//sim.Initialize(1000,/*dt=*/1e-3,/*rho=*/1.05,/*cutoff_r=*/3.5,
	//	/*initial_speed=*/10);
	//sim.InitializePCF(200,3.5);
//@@--------------------------------------------------------------------

//@@--------------------------------------------------------------------
// Note the comments below. PCF = pair correlation function. 
// 
// Note: In case you (do not) want to follow the simulation every <dumpfreq> 
// cycles, uncomment (comment out) the two lines 
// with "//sim.PlotAllParticles();".
// 
// sim.PlotAllParticles() will plot a xy-dump of the particle positions and 
// two or three vertical curves: 
// On the top, the current pair correlation function, 
// in the center the averaged pair correlation funciton with estimated 
//   standard deviation and
// on the bottom a normalized velocity distribution (particle serial 
// number versus speed). 
// 
	int max_steps=100000;  // Max number of <dt> steps to calculate. 
	int pcf_samples=8;   // Compute PCF every this many dt cycles. 
	int pcf_skip=400;    // Skip these many steps before PCF accumulation
	int dumpfreq=100;    // Dump status every this many dt cycles. 
	int currfreq1=10;    // Record p,T,E/N every this many dt cycles after 
	                     // pcf_skip have been skipped. 
	int currfreq0=2;     // Record p,T,E/N every this many dt cycles for 
	                     // iterations 0..pcf_skip. 
	
	// Stop if PCF function accumulation if max. sigma is smaller than this. 
	double pcf_sigma=0.015;
	
	// Output file for current state dump; will be plotted in the end. 
	const char *curr_file="curr.dat";
//@@--------------------------------------------------------------------
	
	FILE *curr_fp=fopen(curr_file,"w");
	if(!curr_fp)
	{  fprintf(stderr,"Failed to open %s.\n",curr_file);  exit(1);  }
	
	sim.PlotAllParticles();
	
	double pcf[sim.get_npcf()];
	double pcf_s[sim.get_npcf()];
	double curr_pcf[sim.get_npcf()];  // <-- Current snapshot. 
	
	// For averaging. 
	double p_sum=0;
	double E_sum=0;
	
	int nacc=0,iter;
	double Ekin_sum=0.0;
	double max_s=-1.0;
	for(iter=0; iter<max_steps; iter++)
	{
		sim.SimulationStep();
		
		bool do_sample = (iter>pcf_skip && !(iter%pcf_samples));
		bool do_dump = (!(iter%dumpfreq) || iter+1==max_steps);
		bool do_write_curr = !(iter%(iter<=pcf_skip ? currfreq0 : currfreq1));
		
		if(do_sample)
		{  sim.AccumulatePCF(curr_pcf);  }
		else if(do_dump || do_write_curr)
		{  sim.ComputePCF(curr_pcf,sim.get_npcf(),sim.get_pcf_rmax());  }
		
		double ptot[3];
		double Epot=-1e10,Ekin=-1e10,E_N=-1e10,T_curr=-1e10,p_curr=-1e10;
		
		if(do_write_curr || do_sample || do_dump)
		{
			// Now as we have the current PCF, we can compute the current 
			// potential energy. 
			Epot=sim.ComputeEpot(curr_pcf,
				sim.get_npcf(),sim.get_pcf_rmax());
			// The kinetic energy is simple since we know all the particle 
			// impulses. 
			Ekin=sim.CalcEnergyAndImpulse(ptot);
			
			// Total energy per particle and the temperature: 
			E_N=(Ekin+Epot)/sim.NP();
			T_curr=sim.Ekin2Temp(Ekin);
			
			// Finally, compute current pressure: 
			p_curr=sim.ComputePressure(curr_pcf,
				sim.get_npcf(),sim.get_pcf_rmax(),T_curr);
		}
		
		if(do_write_curr)  // Write current state: 
		{  fprintf(curr_fp,"%d %g %g %g\n",iter,T_curr,E_N,p_curr);  }
		
		if(do_sample)
		{
			p_sum+=p_curr;
			E_sum+=E_N;
			
			Ekin_sum+=Ekin;
			++nacc;
			
			if(nacc>=5 && !(nacc%5))
			{
				// Compute averaged PCF: 
				max_s=sim.AveragedPCF(pcf,pcf_s);
				if(max_s<pcf_sigma)  break;
			}
		}
		
		if(do_dump)
		{
			fflush(curr_fp);
			// Dump some information to the user: 
			fprintf(stderr,
				"  nskip=%d+%d, ncalc=%d, border=%.1f%%; rmin=%g; max_s=%g\n",
				sim.nskip0,sim.nskip1,sim.ncalc,
				100.0*double(sim.nborder)/sim.NP(),
				sim.r_min,max_s);
			fprintf(stderr,"Iter[%d]: |P|=%.2g, E_kin=%g, T=%g, "
				"E_pot=%g, E/N=%g, p=%g\n",
				iter,VLENGTH3(ptot),Ekin,T_curr,
				Epot,E_N,p_curr);
			sim.PlotAllParticles();
		}
	}
	
	fclose(curr_fp);
	
	timeval end_tv;
	gettimeofday(&end_tv,NULL);
	int64_t start_time=start_tv.tv_sec*1000000LL+start_tv.tv_usec;
	int64_t end_time=end_tv.tv_sec*1000000LL+end_tv.tv_usec;
	long msec=(end_time-start_time+500)/1000;
	
	// Plot the PCF: 
	printf("set title \"PCF: N=%d, rho=%g, iter=%d, pcf_samps=%d  "
		"<E_kin>=%g, T=%g  (%.1f sec)\"\n",
		sim.NP(),sim.get_rho(),iter,nacc,
		Ekin_sum/nacc,sim.Ekin2Temp(Ekin_sum/nacc),
		double(msec/1000.0));
	printf("set xlabel \"r\"\n");
	//printf("set ylabel \"PCF\"\n");
	printf("plot 1 notitle, '-' notitle with errorbars, "
		"'-' title \"PCF\" with lines\n");
	for(int i=0; i<sim.get_npcf(); i++)
	{
		printf("%g %g %g\n",
			(i+0.5)/sim.get_npcf()*sim.get_pcf_rmax(),
			pcf[i],pcf_s[i]);
	}
	printf("e\n");
	for(int i=0; i<sim.get_npcf(); i++)
	{
		printf("%g %g\n",
			(i+0.5)/sim.get_npcf()*sim.get_pcf_rmax(),
			pcf[i]);
	}
	printf("e\n");
	
	// Plot current temperature, energy and pressure versus time. 
	// Several plots of the first iterations specified by percent[]. 
	int percent[]={1,5,10,20,50,100,-1};
	for(int i=0; percent[i]>0; i++)
	{
		int plot_end=percent[i]*(iter+99)/100;
		printf("set title \"Current temperature, energy and pressure "
			"(%s %d steps)\"\n",percent[i]==100 ? "all" : "first",plot_end);
		printf("set xlabel \"t\"\n");
		printf("plot [0:%d] "
			"\"%s\" using ($1):($2) title \"T (avg %.3g)\" with lines,"
			"\"%s\" using ($1):($3) title \"E/N (avg %.3g)\" with lines,"
			"\"%s\" using ($1):($4) title \"p (avg %.3g)\" with lines,"
			"%g notitle with dots lt 1,"
			"%g notitle with dots lt 2,"
			"%g notitle with dots lt 3"
			"\n",plot_end,
			curr_file,sim.Ekin2Temp(Ekin_sum/nacc),
			curr_file,E_sum/nacc,
			curr_file,p_sum/nacc,
			sim.Ekin2Temp(Ekin_sum/nacc),
			E_sum/nacc,
			p_sum/nacc);
	}
}


int main()
{
//@@--------------------------------------------------------------------
// These are the gnuplot commands which set the output: 
	printf("set term postsc color solid\n");
	printf("set outp \"particles.ps\"\n");
//@@--------------------------------------------------------------------
	
	SimFluid();
	
	return(0);
}
