#include <iostream>	// std::cout
#include <cmath>	// pow
#include <assert.h>	// Assertion errors
#include <ctime>	// Time program

// Empirical fit to secondary electron emission equation as in Stangeby
// Electron temperature is expected to be in units of electron volts!
double sec(double Te, char material)
{

	if(Te>0.0)
	{
		double result,x=log10(Te);
		if(material=='c' || material=='C' || material == 'g' || material == 'G' ) // Graphite or Carbon
		{

			result = - 0.0849*x*x*x + 0.1149*x*x + 0.7428*x - 1.341;
		}
		else if(material=='w' || material == 'W')
		{
			result = - 0.0765*x*x*x + 0.1521*x*x + 0.724*x - 1.4755;
		}
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		else if(material=='f' || material == 'F' )
		{
            result = - 0.0903*x*x*x + 0.1813*x*x + 0.6368*x - 1.2668;
			//result = 1.3*2.72*2.72*(2*Te/400)*exp(-2*sqrt(2*Te/400));
		}
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		else if(material=='b' || material == 'B')//Eq 3.5 Stangeby with E0=2*Te 
		{
            result = - 0.0906*x*x*x + 0.1018*x*x + 0.717*x - 1.4729;
			//result = 0.5*2.72*2.72*(2*Te/200)*exp(-2*sqrt(2*Te/200));
		    //result = - 0.0849*x*x*x + 0.1149*x*x + 0.7428*x - 1.341;
		}
        else if(material=='l' || material == 'L')
		{
            result = - 0.0648*x*x*x - 0.1358*x*x + 0.9887*x - 1.3253;
		}
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		else 
		{
			std::cout << "Error: Incorrect material information in sec" << std::endl;
		}
		if((material=='c')||(material=='C')||(material=='W')||(material=='w')){return pow(10.0,result);}
		else{return result;}
	}
	else return 0.0;
}

void WarnOnce(bool &T, std::string Message){
	if(T){
		std::cout << "\n\n*[W]* Warning! " << Message;
		T = false;
	}
}


double maxwellian(double E, double T)
{

	double PI = 3.14159265359;
	return 2*sqrt(E/(PI*T))*exp(-E/T)/T;
}

double ionback(double E, char isotope, char material, int flag)
{
	double M1,M2,Z1=1.0,Z2;
	double A1,A2,A3,A4,A5,A6;
// Hydrogen isotope
	if(isotope=='h')
	{
		M1 = 1.0;
	}
	else if(isotope=='d')
	{
		M1 = 2.0;
	}
	else if(isotope=='t')
	{
		M1 = 3.0;
	}
	else
	{
		std::cout << "Error: Incorrect isotope information in ionback" << std::endl;
		exit(1);
	}
	if(material=='c' || material == 'C')
	{
		M2 = 12.0;
		Z2 = 6.0;
// flag=0 returns energy
		if(flag==0)
		{
			A1 = 0.4484; A2 = 27.16; A3 = 15.66; A4 = 0.6598; A5 = 7.967; A6 = 1.822;
		}
// flag=1 returns fraction of particles
		else if(flag==1)
		{
            A1 = 0.6192; A2 = 20.01; A3 = 8.922; A4 = 0.6669; A5 = 1.864; A6 = 1.899;
		}
		else
		{
			std::cout << "Error: Incorrect flag in ionback" << std::endl;
			exit(1);
		}
	}
	else if(material=='w' || material == 'W')
	{
		M2 = 183.84;
		Z2 = 74.0;
// flag=0 returns energy
		if(flag==0)
		{
			A1 = 0.6831; A2 = 27.16; A3 = 15.66; A4 = 0.6598; A5 = 7.967; A6 = 1.822;
		}
// flag=1 returns fraction of particles
		else if(flag==1)
		{
            A1 = 0.8250; A2 = 21.41; A3 = 8.606; A4 = 0.6425; A5 = 1.907; A6 = 1.927;
		}
		else
		{
			std::cout << "Error: Incorrect flag in ionback" << std::endl;
			exit(1);
		}
	}
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	else if(material=='f' || material == 'F')//I just use the same values as for tungsten cause I could't find reference
	{
		M2 = 55.845;
		Z2 = 26.0;
// flag=0 returns energy
		if(flag==0)
		{
			A1 = 0.6831; A2 = 27.16; A3 = 15.66; A4 = 0.6598; A5 = 7.967; A6 = 1.822;
		}
// flag=1 returns fraction of particles
		else if(flag==1)
		{
            A1 = 0.8250; A2 = 21.41; A3 = 8.606; A4 = 0.6425; A5 = 1.907; A6 = 1.927;
		}
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    }		
			//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	else if(material=='b' || material == 'B')//I use the values for lithium and take the average with carbon
	{
		M2 = 6.941;
		Z2 = 3.0;
// flag=0 returns energy
		if(flag==0)
		{
			A1 = 0.4222; A2 = 3.092; A3 = 13.17; A4 = 0.5393; A5 = 4.464; A6 = 1.877;
		}
// flag=1 returns fraction of particles
		else if(flag==1)
		{
            A1 = 0.5173; A2 = 2.549; A3 = 5.325; A4 = 0.5719; A5 = 1.094; A6 = 1.933;
		}
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		else
		{
			std::cout << "Error: Incorrect flag in ionback" << std::endl;
			exit(1);
		}

		
		/*M2 = 12.0;
		Z2 = 6.0;
// flag=0 returns energy
		if(flag==0)
		{
			A1 = 0.4484; A2 = 27.16; A3 = 15.66; A4 = 0.6598; A5 = 7.967; A6 = 1.822;
		}
// flag=1 returns fraction of particles
		else if(flag==1)
		{
            A1 = 0.6192; A2 = 20.01; A3 = 8.922; A4 = 0.6669; A5 = 1.864; A6 = 1.899;
		}
		else
		{
			std::cout << "Error: Incorrect flag in ionback" << std::endl;
			exit(1);
		}*/
	}

	else if(material=='l' || material == 'L')
	{
		M2 = 6.941;
		Z2 = 3.0;
// flag=0 returns energy
		if(flag==0)
		{
			A1 = 0.4222; A2 = 3.092; A3 = 13.17; A4 = 0.5393; A5 = 4.464; A6 = 1.877;
		}
// flag=1 returns fraction of particles
		else if(flag==1)
		{
            A1 = 0.5173; A2 = 2.549; A3 = 5.325; A4 = 0.5719; A5 = 1.094; A6 = 1.933;
		}
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		else
		{
			std::cout << "Error: Incorrect flag in ionback" << std::endl;
			exit(1);
		}
	}
	else 
	{
		std::cout << "Error: Incorrect material information in ionback" << std::endl;
		exit(1);
	}
	double epsilon = 0.032534*E*M2/((M1+M2)*Z1*Z2*pow(pow(Z1,2.0/3.0)+pow(Z2,2.0/3.0),0.5));
	double R = A1*log(A2*epsilon+exp(1.0))/(1.0+A3*pow(epsilon,A4)+A5*pow(epsilon,A6));
	return R;
}

double backscatter(double Te, double Ti, double mi, double Vion, char material, 
				   double &RE, double &RN)
{
clock_t begin = clock();

double PI = 3.14159265359;
if(material!='b' && material !='B')
 {
    
	if(Ti>0.0)
	{
		int n,nmax=20000;
		double E,Emax=20*Te,h=Emax/double(nmax),energytot=0.0,particletot=0.0;
		for(n=1;n<=nmax;n++)
		{
			E = n*h;
			if((n==1)||(n==nmax))
			{
				energytot += maxwellian(E,Ti)*E*sqrt(2.0*E/mi)*(1.0+Vion*Te/E)
					*ionback(E+Vion*Te,'h',material,0);
				particletot += maxwellian(E,Ti)*sqrt(2.0*E/mi)*(1.0+Vion*Te/E)
					*ionback(E+Vion*Te,'h',material,1);
			}
			else if(n%2==0)
			{
				energytot += 2.0*maxwellian(E,Ti)*E*sqrt(2.0*E/mi)*(1.0+Vion*Te/E)
					*ionback(E+Vion*Te,'h',material,0);
				particletot += 2.0*maxwellian(E,Ti)*sqrt(2.0*E/mi)*(1.0+Vion*Te/E)
					*ionback(E+Vion*Te,'h',material,1);
			}
			else
			{
				energytot += 4.0*maxwellian(E,Ti)*E*sqrt(2.0*E/mi)*(1.0+Vion*Te/E)
					*ionback(E+Vion*Te,'h',material,0);
				particletot += 4.0*maxwellian(E,Ti)*sqrt(2.0*E/mi)*(1.0+Vion*Te/E)
					*ionback(E+Vion*Te,'h',material,1);
			}
		}
		energytot = h*energytot/3.0;
		particletot = h*particletot/3.0;
// Normalise
		RE = energytot/(2*Ti*sqrt(8*Ti/(PI*mi))*(1.0+Vion*Ti/Te));
		RN = particletot/(sqrt(8*Ti/(PI*mi))*(1.0+Vion*Ti/Te));
		if(RN>1.0) RN = 1.0;
	}
	else
	{
		RE = 0.0;
		RN = 0.0;
	}
 }
 if(material=='b' && material !='B')	
 {
     if(Ti>0.0)
	{
		int n,nmax=20000;
		double E,Emax=20*Te,h=Emax/double(nmax),energytot=0.0,particletot=0.0;
		for(n=1;n<=nmax;n++)
		{
			E = n*h;
			if((n==1)||(n==nmax))
			{
				energytot += maxwellian(E,Ti)*E*sqrt(2.0*E/mi)*(1.0+Vion*Te/E)
					*ionback(E+Vion*Te,'h',material,0);
				particletot += maxwellian(E,Ti)*sqrt(2.0*E/mi)*(1.0+Vion*Te/E)
					*ionback(E+Vion*Te,'h',material,1);
			}
			else if(n%2==0)
			{
				energytot += 2.0*maxwellian(E,Ti)*E*sqrt(2.0*E/mi)*(1.0+Vion*Te/E)
					*ionback(E+Vion*Te,'h',material,0);
				particletot += 2.0*maxwellian(E,Ti)*sqrt(2.0*E/mi)*(1.0+Vion*Te/E)
					*ionback(E+Vion*Te,'h',material,1);
			}
			else
			{
				energytot += 4.0*maxwellian(E,Ti)*E*sqrt(2.0*E/mi)*(1.0+Vion*Te/E)
					*ionback(E+Vion*Te,'h',material,0);
				particletot += 4.0*maxwellian(E,Ti)*sqrt(2.0*E/mi)*(1.0+Vion*Te/E)
					*ionback(E+Vion*Te,'h',material,1);
			}
		}
		energytot = h*energytot/3.0;
		particletot = h*particletot/3.0;
// Normalise
		RE = energytot/(2*Ti*sqrt(8*Ti/(PI*mi))*(1.0+Vion*Ti/Te));
		RN = particletot/(sqrt(8*Ti/(PI*mi))*(1.0+Vion*Ti/Te));
		if(RN>1.0) RN = 1.0;
		energytot=0.0;
        particletot=0.0;
        for(n=1;n<=nmax;n++)
		{
			E = n*h;
			if((n==1)||(n==nmax))
			{
				energytot += maxwellian(E,Ti)*E*sqrt(2.0*E/mi)*(1.0+Vion*Te/E)
					*ionback(E+Vion*Te,'h','c',0);
				particletot += maxwellian(E,Ti)*sqrt(2.0*E/mi)*(1.0+Vion*Te/E)
					*ionback(E+Vion*Te,'h','c',1);
			}
			else if(n%2==0)
			{
				energytot += 2.0*maxwellian(E,Ti)*E*sqrt(2.0*E/mi)*(1.0+Vion*Te/E)
					*ionback(E+Vion*Te,'h','c',0);
				particletot += 2.0*maxwellian(E,Ti)*sqrt(2.0*E/mi)*(1.0+Vion*Te/E)
					*ionback(E+Vion*Te,'h','c',1);
			}
			else
			{
				energytot += 4.0*maxwellian(E,Ti)*E*sqrt(2.0*E/mi)*(1.0+Vion*Te/E)
					*ionback(E+Vion*Te,'h','c',0);
				particletot += 4.0*maxwellian(E,Ti)*sqrt(2.0*E/mi)*(1.0+Vion*Te/E)
					*ionback(E+Vion*Te,'h','c',1);
			}
		}
		energytot = h*energytot/3.0;
		particletot = h*particletot/3.0;
// Normalise
		RE = (RE+energytot/(2*Ti*sqrt(8*Ti/(PI*mi))*(1.0+Vion*Ti/Te)))/2;
		RN = (RN+particletot/(sqrt(8*Ti/(PI*mi))*(1.0+Vion*Ti/Te)))/2;
		if(RN>1.0) RN = 1.0;
	}
	else
	{
		RE = 0.0;
		RN = 0.0;
	}
 }    
//	assert(RE < 1);
//	assert(RN < 1);
//	if(RE != RE){
//		std::cout << "\nError! In void backscatter, RE = " << RE;
//		throw std::exception();
//	}
//	if(RN != RN){
//		std::cout << "\nError! In void backscatter, RN = " << RN;
//		throw std::exception();
//	}
	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
	return elapsd_secs;
}

double round_to_digits(double value, int digits){
    if (value == 0.0) // otherwise it will return 'nan' due to the log10() of zero
        return 0.0;

    double factor = pow(10.0, digits - ceil(log10(fabs(value))));
    return round(value * factor) / factor;   
}

double LambertW(const double z) {
  int i; 
  const double eps=4.0e-16, em1=0.3678794411714423215955237701614608; 
  double p,e,t,w;
//  if (dbgW) fprintf(stderr,"LambertW: z=%g\n",z);
  if (z<-em1 || std::isinf(z) || std::isnan(z)) { 
	fprintf(stderr,"LambertW: bad argument %g, exiting.\n",z); // exit(1); 
	return 0.0;
  }
  if (0.0==z) return 0.0;
  if (z<-em1+1e-4) { // series near -em1 in sqrt(q)
    double q=z+em1,r=sqrt(q),q2=q*q,q3=q2*q;
    return 
     -1.0
     +2.331643981597124203363536062168*r
     -1.812187885639363490240191647568*q
     +1.936631114492359755363277457668*r*q
     -2.353551201881614516821543561516*q2
     +3.066858901050631912893148922704*r*q2
     -4.175335600258177138854984177460*q3
     +5.858023729874774148815053846119*r*q3
     -8.401032217523977370984161688514*q3*q;  // error approx 1e-16
  }
  /* initial approx for iteration... */
  if (z<1.0) { /* series near 0 */
    p=sqrt(2.0*(2.7182818284590452353602874713526625*z+1.0));
    w=-1.0+p*(1.0+p*(-0.333333333333333333333+p*0.152777777777777777777777)); 
  } else 
    w=log(z); /* asymptotic */
  if (z>3.0) w-=log(w); /* useful? */
  for (i=0; i<10; i++) { /* Halley iteration */
    e=exp(w); 
    t=w*e-z;
    p=w+1.0;
    t/=e*p-0.5*(p+1.0)*t/p; 
    w-=t;
    if (fabs(t)<eps*(1.0+fabs(w))) return w; /* rel-abs error */
  }
  /* should never get here */
  fprintf(stderr,"LambertW: No convergence at z=%g, exiting.\n",z); 
  exit(1);
}

