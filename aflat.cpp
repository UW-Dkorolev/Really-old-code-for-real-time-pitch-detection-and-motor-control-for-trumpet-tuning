#include<cstdio>
#include<cmath>
#include<cstdlib>
#include<sys/types.h>
#include<errno.h>
#include<sys/wait.h>
#include<signal.h>
#include <malloc.h>
#include <unistd.h>
//#include<unistd.h>
#include<cstring>
#include <dirent.h>
#include <fcntl.h>
#include <assert.h> 
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#define PATH "/home/pi/Desktop/cher.wav"
#define TWO_PI (6.2831853071795864769252867665590057683943L)  

#define BCM2708_PERI_BASE  0x20000000
#define GPIO_BASE		(BCM2708_PERI_BASE + 0x200000) /* GPIO controller */
#define PWM_BASE		(BCM2708_PERI_BASE + 0x20C000) /* PWM controller */
#define CLOCK_BASE		(BCM2708_PERI_BASE + 0x101000)

#define	PWM_CTL  0
#define	PWM_RNG1 4
#define	PWM_DAT1 5

#define	PWMCLK_CNTL 40
#define	PWMCLK_DIV  41
#define PAGE_SIZE (4*1024)
#define BLOCK_SIZE (4*1024)


// I/O access
volatile unsigned *gpio;
volatile unsigned *pwm;
volatile unsigned *clk;

// GPIO setup macros. Always use INP_GPIO(x) before using OUT_GPIO(x) or SET_GPIO_ALT(x,y)
#define INP_GPIO(g) *(gpio+((g)/10)) &= ~(7<<(((g)%10)*3))
#define OUT_GPIO(g) *(gpio+((g)/10)) |=  (1<<(((g)%10)*3))
#define SET_GPIO_ALT(g,a) *(gpio+(((g)/10))) |= (((a)<=3?(a)+4:(a)==4?3:2)<<(((g)%10)*3))

#define GPIO_SET *(gpio+7)  // sets   bits which are 1 ignores bits which are 0
#define GPIO_CLR *(gpio+10) // clears bits which are 1 ignores bits which are 0

// map 4k register memory for direct access from user space and return a user space pointer to it
volatile unsigned *mapRegisterMemory(int base)
{
	static int mem_fd = 0;
	char *mem, *map;
	
	/* open /dev/mem */
	if (!mem_fd) {
		if ((mem_fd = open("/dev/mem", O_RDWR|O_SYNC) ) < 0) {
			printf("can't open /dev/mem \n");
			exit (-1);
		}
	}
	
	/* mmap register */
	
	// Allocate MAP block
	if ((mem = malloc(BLOCK_SIZE + (PAGE_SIZE-1))) == NULL) {
		printf("allocation error \n");
		exit (-1);
	}
	
	// Make sure pointer is on 4K boundary
	if ((unsigned long)mem % PAGE_SIZE)
		mem += PAGE_SIZE - ((unsigned long)mem % PAGE_SIZE);
	
	// Now map it
	map = (char *)mmap(
		(caddr_t)mem,
		BLOCK_SIZE,
		PROT_READ|PROT_WRITE,
		MAP_SHARED|MAP_FIXED,
		mem_fd,
		base
	);
	
	if ((long)map < 0) {
		printf("mmap error %d\n", (int)map);
		exit (-1);
	}
	
	// Always use volatile pointer!
	return (volatile unsigned *)map;
}

// set up a memory regions to access GPIO, PWM and the clock manager
void setupRegisterMemoryMappings()
{
	gpio = mapRegisterMemory(GPIO_BASE);
	pwm = mapRegisterMemory(PWM_BASE);
	clk = mapRegisterMemory(CLOCK_BASE);
}

void setServo(int percent)
{
	int bitCount;
	unsigned int bits = 0;

	// 32 bits = 2 milliseconds
	bitCount = 16 + 16 * percent / 100;
	if (bitCount > 32) bitCount = 32;
	if (bitCount < 1) bitCount = 1;
	bits = 0;
	while (bitCount) {
		bits <<= 1;
		bits |= 1;
		bitCount--;
	}
	*(pwm + PWM_DAT1) = bits;
}

// init hardware
void initHardware()
{
	// mmap register space
	setupRegisterMemoryMappings();
	
	// set PWM alternate function for GPIO18
	SET_GPIO_ALT(18, 5);

	// stop clock and waiting for busy flag doesn't work, so kill clock
	*(clk + PWMCLK_CNTL) = 0x5A000000 | (1 << 5);
	usleep(10);  

	// set frequency
	// DIVI is the integer part of the divisor
	// the fractional part (DIVF) drops clock cycles to get the output frequency, bad for servo motors
	// 320 bits for one cycle of 20 milliseconds = 62.5 us per bit = 16 kHz
	int idiv = (int) (19200000.0f / 16000.0f);
	if (idiv < 1 || idiv > 0x1000) {
		printf("idiv out of range: %x\n", idiv);
		exit(-1);
	}
	*(clk + PWMCLK_DIV)  = 0x5A000000 | (idiv<<12);
	
	// source=osc and enable clock
	*(clk + PWMCLK_CNTL) = 0x5A000011;

	// disable PWM
	*(pwm + PWM_CTL) = 0;
	
	// needs some time until the PWM module gets disabled, without the delay the PWM module crashs
	usleep(10);  
	
	// filled with 0 for 20 milliseconds = 320 bits
	*(pwm + PWM_RNG1) = 320;
	
	// 32 bits = 2 milliseconds, init with 1 millisecond
	setServo(-90);
	
	// start PWM1 in serializer mode
	*(pwm + PWM_CTL) = 3;
}

/*----------------------------------------------------------------------------
   fft.c - fast Fourier transform and its inverse (both recursively)
   Copyright (C) 2004, Jerome R. Breitenbach.  All rights reserved.

   The author gives permission to anyone to freely copy, distribute, and use
   this file, under the following conditions:
      - No changes are made.
      - No direct commercial advantage is obtained.
      - No liability is attributed to the author for any damages incurred.
  ----------------------------------------------------------------------------*/

/******************************************************************************
 * This file defines a C function fft that, by calling another function       *
 * fft_rec (also defined), calculates an FFT recursively.  Usage:             *
 *   fft(N, x, X);                                                            *
 * Parameters:                                                                *
 *   N: number of points in FFT (must equal 2^n for some integer n >= 1)      *
 *   x: pointer to N time-domain samples given in rectangular form (Re x,     *
 *      Im x)                                                                 *
 *   X: pointer to N frequency-domain samples calculated in rectangular form  *
 *      (Re X, Im X)                                                          *
 * Similarly, a function ifft with the same parameters is defined that        *
 * calculates an inverse FFT (IFFT) recursively.  Usage:                      *
 *   ifft(N, x, X);                                                           *
 * Here, N and X are given, and x is calculated.                              *
 ******************************************************************************/


/* function prototypes */
void fft(int N, double (*x)[2], double (*X)[2]);
void fft_rec(int N, int offset, int delta,
             double (*x)[2], double (*X)[2], double (*XX)[2]);
void ifft(int N, double (*x)[2], double (*X)[2]);

/* FFT */
void fft(int N, double (*x)[2], double (*X)[2])
{
  /* Declare a pointer to scratch space. */
  double (*XX)[2] = (double (*)[2])malloc(2 * N * sizeof(double));

  /* Calculate FFT by a recursion. */
  fft_rec(N, 0, 1, x, X, XX);

  /* Free memory. */
  free(XX);
}

/* FFT recursion */
void fft_rec(int N, int offset, int delta,
             double (*x)[2], double (*X)[2], double (*XX)[2])
{
  int N2 = N/2;            /* half the number of points in FFT */
  int k;                   /* generic index */
  double cs, sn;           /* cosine and sine */
  int k00, k01, k10, k11;  /* indices for butterflies */
  double tmp0, tmp1;       /* temporary storage */

  if(N != 2)  /* Perform recursive step. */
    {
      /* Calculate two (N/2)-point DFT's. */
      fft_rec(N2, offset, 2*delta, x, XX, X);
      fft_rec(N2, offset+delta, 2*delta, x, XX, X);

      /* Combine the two (N/2)-point DFT's into one N-point DFT. */
      for(k=0; k<N2; k++)
        {
          k00 = offset + k*delta;    k01 = k00 + N2*delta;
          k10 = offset + 2*k*delta;  k11 = k10 + delta;
          cs = cos(TWO_PI*k/(double)N); sn = sin(TWO_PI*k/(double)N);
          tmp0 = cs * XX[k11][0] + sn * XX[k11][1];
          tmp1 = cs * XX[k11][1] - sn * XX[k11][0];
          X[k01][0] = XX[k10][0] - tmp0;
          X[k01][1] = XX[k10][1] - tmp1;
          X[k00][0] = XX[k10][0] + tmp0;
          X[k00][1] = XX[k10][1] + tmp1;
        }
    }
  else  /* Perform 2-point DFT. */
    {
      k00 = offset; k01 = k00 + delta;
      X[k01][0] = x[k00][0] - x[k01][0];
      X[k01][1] = x[k00][1] - x[k01][1];
      X[k00][0] = x[k00][0] + x[k01][0];
      X[k00][1] = x[k00][1] + x[k01][1];
    }
}

/* IFFT */
void ifft(int N, double (*x)[2], double (*X)[2])
{
  int N2 = N/2;       /* half the number of points in IFFT */
  int i;              /* generic index */
  double tmp0, tmp1;  /* temporary storage */

  /* Calculate IFFT via reciprocity property of DFT. */
  fft(N, X, x);
  x[0][0] = x[0][0]/N;    x[0][1] = x[0][1]/N;
  x[N2][0] = x[N2][0]/N;  x[N2][1] = x[N2][1]/N;
  for(i=1; i<N2; i++)
    {
      tmp0 = x[i][0]/N;       tmp1 = x[i][1]/N;
      x[i][0] = x[N-i][0]/N;  x[i][1] = x[N-i][1]/N;
      x[N-i][0] = tmp0;       x[N-i][1] = tmp1;
    }
}

// linear interpolate x in an array
// inline

float interp1( float x, float a[], int n )
{
    if( x <= 0 ){  
	return a[0];}
    if( x >= n - 1 ){  
	return a[n-1];}
	int j = (int)x;
    return a[j] + (x - j) * (a[j+1] - a[j]);
}

    // linear interpolate array a[] -> array b[]
void inter1parray( float a[], int n, float b[], int m )
{
    float step = (float)( n - 1 ) / (m - 1);
    int j = 0;
	for(j = 0; j < m; j ++ ){
        b[j] = interp1( j*step, a, n );
    }
}

//..............................................................................
    // parabola through 3 points, -1 < x < 1
float parabola( float x, float f_1, float f0, float f1 )
{
    if( x <= -1 )  return f_1; 
    if( x >= 1 )  return f1; 
    float l = f0 - x * (f_1 - f0);
    float r = f0 + x * (f1 - f0);
    return (l + r + x * (r - l)) / 2;
}

    // quadratic interpolate x in an array
float interp2( float x, float a[], int n )
{
    if( x <= .5  ||  x >= n - 1.5 )
        return interp1( x, a, n );
    int j = (int)( x + .5 );
    float t = 2 * (x - j);  // -1 .. 1
    return parabola( t, (a[j-1] + a[j]) / 2, a[j], (a[j] + a[j+1]) / 2 );
}

    // quadratic interpolate array a[] -> array b[]
void interp2array( float a[], int n, float b[], int m )
{
    float step = (float)( n - 1 ) / (m - 1);
    int j = 0;
	for(j = 0; j < m; j ++ ){
        b[j] = interp2( j*step, a, n );
    }
}

int main(){
char cmd[60];
double ADC_OUTPUT[1024][2];
int tau[384]={0},i = 0,aaa = 201,zz = 0,ppp=0;
int N = 1024; //window size
double SNAC[384] ={0},pitch_est = 0,m[384] = {0},sdf[384] ={0};
float SNACcher[21] = {0},SNACINTERP[201] = {0};
double notes[25] = {233.08,246.94,261.63,277.18,293.66,311.13,329.63,349.23,369.99,392.00,415.3,440.0,466.16,493.88,523.25,554.37,587.33,622.25,659.26,698.46,739.99,830.61,880.0,932.33};
double error[25] = {0},minsofar = 0;
static const int avlim = 3; //the number of iterations of avloop
int SignalPower[avlim] = {0};//should be the same size as avlim
int threshold = 700; // this decides whether or not the trumpet is being played and thus whether or not the system should take action
int loc = 0;//this will be used for a search through notes[] later
double pitcherr = 0;
int beep = 21;
int dd = 0;
float gg = 0;
int pos = 60; //position of the servo (-90 seems to correspond to all the way left)
int avloop = 0;
double maxnow = 0;
int zzmax = 0;
initHardware(); //get that servo ready to go
sleep(2);
for(ppp=0;ppp<50;ppp++)
{
	pitch_est = 0;
	printf("pos = %d\n",pos);
	setServo(pos);
	system("cd /");
	system("cd /home/pi/Desktop");
	system("arecord -f cd -t wav -d 1 -D plughw:1,0 cher.wav");


	// Buffers etc..
	char ChunkID[4], Format[4], Subchunk1ID[4],Subchunk2ID[4];
	int ChunkSize,Subchunk1Size, SampleRate, ByteRate,Subchunk2Size;
	short AudioFormat, NumChannels, BlockAlign, BitsPerSample;
	short *Data;
	//=(short *)malloc();
	
	// Read the wave file
	FILE *fhandle=fopen(PATH,"rb");
	fread(ChunkID,1,4,fhandle);
	fread(&ChunkSize,4,1,fhandle);
	fread(Format,1,4,fhandle);
	fread(Subchunk1ID,1,4,fhandle);
	fread(&Subchunk1Size,4,1,fhandle);
	fread(&AudioFormat,2,1,fhandle);
	fread(&NumChannels,2,1,fhandle);
	fread(&SampleRate,4,1,fhandle);
	fread(&ByteRate,4,1,fhandle);
	fread(&BlockAlign,2,1,fhandle);
	fread(&BitsPerSample,2,1,fhandle);
	fread(&Subchunk2ID,1,4,fhandle);
	fread(&Subchunk2Size,4,1,fhandle);
	Data=new short [Subchunk2Size/(BitsPerSample/8)]; // Create an element for every sample
	fread(Data,BitsPerSample/8,Subchunk2Size/(BitsPerSample/8),fhandle); // Reading raw audio data
		fclose(fhandle);
	for(avloop=0;avloop<avlim;avloop++)
	{
		memset(SNACcher,0,21);memset(SNACINTERP,0,201);
		memset(tau,0,384);
		memset(SNAC,0,384);memset(m,0,384);memset(sdf,0,384);

		for(i=12;i<=1034;i++)

		{
			ADC_OUTPUT[i-12][0] = (double)Data[i+1023*avloop];
			ADC_OUTPUT[i-12][1] = 0;
		}
		maxnow = 0;
		zzmax = 0;
		for(i=0;i<1024;i++)
		{
			if(ADC_OUTPUT[i][0]>maxnow)
			{
				maxnow = ADC_OUTPUT[i][0];
				zzmax = i;
			}
		}
		printf("ADCOUT_max = %f at sample %d\n",maxnow,zzmax);
		maxnow =0;
		zzmax = 0;
		i = 0;
		SignalPower[avloop] = 0;
		for(i=0;i<1023;i++)//loop through the ADC_OUTPUT array
			{
				SignalPower[avloop] += ADC_OUTPUT[i][0]*ADC_OUTPUT[i][0];
			}
		SignalPower[avloop] = SignalPower[avloop]/1023; //Get the average power of the signal over the frame by dividing by the number of terms
		printf("sigpower[%d] is %d\n",avloop,SignalPower[avloop]);
		if(SignalPower[avloop] <= threshold)
		{
			goto end; //skip to the end if the signal power is too low
		}
			
		
			for(zz = 0;zz< (384);zz++)
		{
			//printf("%d\n", zz);
			tau[zz] = zz; //time variable 
		}




	
		char file[FILENAME_MAX];  /* name of data file */
		/*N =  number of points in FFT */
		double (*x)[2];           /* pointer to time-domain samples */
		double (*X)[2];           /* pointer to frequency-domain samples */
		double dummy;             /* scratch variable */



		/* Allocate time- and frequency-domain memory. */
		x = (double (*)[2])malloc(2 * N * sizeof(double));
		X = (double (*)[2])malloc(2 * N * sizeof(double));


		/* Calculate FFT. */
		fft(N, ADC_OUTPUT, X);

		//making the m variable for Square Difference Function (SDF) 
		for (i = 0; i<=(511); i++) 
			{m[0] = m[0] + 2*(ADC_OUTPUT[i][0])*(ADC_OUTPUT[i][0]);}

		for (zz=1; zz<=tau[382]; zz++) //still making m
			{m[zz] = m[zz-1] - (ADC_OUTPUT[zz-1][0])*(ADC_OUTPUT[zz-1][0]) - (ADC_OUTPUT[512-zz][0])*(ADC_OUTPUT[512-zz][0]);}
		//m will be used later
	//free(ADC_OUTPUT);

	for(zz = 0; zz<=(1023);zz++) // This sets X equal to the PSD
	{
		X[zz][0] = X[zz][0]*X[zz][0] + X[zz][1]*X[zz][1];
		X[zz][1] = 0;
	}
	//for(i=0; i<100; i++) //printf("\n   k=%d: %12f %12f", i, X[i][0], X[i][1]);
	/* Clear time-domain samples and calculate IFFT. */
	for(i=0; i<N; i++) 
		{x[i][0] = x[i][1] = 0;}
	ifft(N, x, X); //Through Wiener-Khinchin, this sets x = autocorrelation.
	//for(i=0; i<100; i++) printf("\n   k=%d: %12f", i, x[i][0]);
	 //so far so good  


	for (zz=0;zz<=383;zz++) //calculating the SDF
	{
		sdf[zz] = (m[zz] - 2*x[zz][0]);
	}

	for (zz=0;zz<=382;zz++) // calculating the Specially Normalized Autocorrelation (SNAC) 
	// don't calculate at zz = 383 because m = 0 there
	{
		SNAC[zz] = 1-(sdf[zz]/m[zz]);
	}  
	maxnow = 0;
	zzmax = 0;
	for(zz=1;zz<383;zz++)
	{
		if(SNAC[zz]>maxnow)
		{
			maxnow = SNAC[zz];
			zzmax = zz;
		}
	}
	printf("SNAC max is %f at sample %d\n",SNAC[zzmax],zzmax); 
	for (zz=0;zz<=383;zz++) //Thresholding the SNAC to get rid of unnecessary peaks
	{
		if (SNAC[zz]<=.8)
		{
			SNAC[zz]=0;
		}
	}
	 //for(i=0; i<120; i++) //printf("\n   k=%d: snac = %12f", i, SNAC[i]);
	 //brute force find the first peak
	dd = 0; //dummy variable
	zz = 35; //starting at 36th element of SNAC because that corresponds to 1200Hz, the highest frequency we would ever expect from someone using this.
	while (dd == 0)
	{
		if (SNAC[zz] == 0)
		{
			if(zz>382)
			{
				printf("We broke because all of SNAC is zero\n");
				pitch_est = pitch_est + 0.0;
				goto end;
			}
			zz++;
			continue;
		}
		else if  (SNAC[zz]>=SNAC[zz-1] && SNAC[zz] >= SNAC[zz+1])
		{
			printf("We broke by finding a maximum\n");
			break;
		}
		else if (zz<383)
		{
			//printf("zz in the nonzero,nonmax = %d\n",zz);
			zz++;
		}
		else 
		{
			printf("We broke because of else\n");
			goto end;
		}
	}
	printf("zz = %d\n",zz); 
	for (i = 0;i<22;i++) //Making a variable that is part of the SNAC array
	{
		SNACcher[i] = (float)SNAC[zz-10+i];
	}
	

	interp2array(SNACcher,beep,SNACINTERP,aaa); //interpolate around that peak to get more accuracy
	//for(i=0;i<=10;i++)
	//{
		//printf("SNACINTERP[%d] = %f SNACcher[%d] = %f\n",i,SNACINTERP[i],i,SNACcher[i]);
	//}
	//for(i=0;i<=10;i++)
	//{
	//printf("sdf[%d] = %f m[%d] = %f\n",i,sdf[i],i,m[i]);
	//}
	i = 0;
	dd = 0;
	while (dd == 0) //BRUTE FORCE way to find the peak more accurately
	{ 
		if (SNACINTERP[i] == 0)
		{
			if(i>200)
			{
				break;
			}
			i++;
			continue;
		}
		else if (SNACINTERP[i]>=SNACINTERP[i-1] && SNACINTERP[i] >= SNACINTERP[i+1])
		{
			break;
		}
		else if (i<200)
		{
			i++;
		}
		else 
		{
			break;
		}
	 }
		printf("i = %d and is used to calculate gg\n",i);
		gg = (float)zz;
		printf("interpolated max is at %f\n", gg);
		//map this i value to the location of the peak by relating to zz, the less accurate peak around which interpolation was centered
		gg = gg +(0.1*(i-100)); //mapping back to time difference space
		//printf("final peak location is at = %f (SNAC without FFT is getting destroyed and turned to zero)\n",gg);
		printf("gg = %f and is used to caluclate pitch_est[i]\n",gg);
		pitch_est += 2*44100/gg;
		printf("pitch_est[%d] is %f\n",avloop,pitch_est);
		end:
		printf("end\n");
	}
	loc = 0; //going to use this variable to sum up how many frames are "voiced" (aka signal power is greater than a threshold)
	for(i=0;i<avlim;i++)
	{
		if(SignalPower[i] > threshold)
			SignalPower[i] = 1;
		else
			SignalPower[i] = 0;
	loc = loc+SignalPower[i]; //this will tell us how many frames are voiced
	}
	if (loc == 0)
	{
		continue;
	}
	for(i=0;i<avlim;i++)
	{
		printf("Power[%d] > thresh = %d\n",i,SignalPower[i]);
	}
	printf("loc = %d\n",loc);
	pitch_est = pitch_est/(double)loc;
	//find what note you are playing
	for (i = 0;i<25;i++)
	{
		error[i] = abs(notes[i]-pitch_est);
		//printf("error[i] = %f\n",error[i]);
	}
	minsofar = error[0];
	loc = 0;//now using loc to find closest note to what you are playing
	for (i = 1;i<25;i++)
	{
		if (error[i] < minsofar)
		{
			minsofar = error[i];
			loc = i;
		}
		//printf("min = %f loc = %d",minsofar,loc);
	}
	//printf("\n");
	//printf("loc = %d notes[loc] = %f\n",loc,notes[loc]);
	pitcherr = 3986.3137*log10(pitch_est/notes[loc]);
	printf("Pitch Estimate = %d Hz\n",(int)pitch_est);	
	printf("Correct pitch = %d Hz Pitch error = %d cents\n",(int)notes[loc],(int)pitcherr);
	if (abs(pitcherr) >(double)10.0)
	{
		printf("pitcherr was greater than 10\n");
		pos = pos + pitcherr/0.4324;
	}
	else
	{
		printf("pitcherr <= 10 good job!\n");
	}
	if (pos >=60)
	{
		if (pitcherr > (double)10.0)
		{
			pos = 0;
		}
		else
		{
			pos = 60;
		}	
	}	
	if (pos <=-90)
	{
		pos = -90;
	}	
	
	//sleep(500);
	

}
}
