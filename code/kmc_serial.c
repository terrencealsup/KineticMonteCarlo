/*
does KMC for specified initial profiles + parameters L, K, rate type (adatom or dH) and times
recording average height at 0 and specified times T_1, T_2,...,T_{num_times}
Input:
command line: seed number 
files: parameters.txt
     
Output: ./h.txt, a (num_times +1) x L array, storing avg. height profile at 0, T_1,T_2,...
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "mt64.h"

#ifndef uniform64
#define uniform64() genrand64_real3()
#endif

typedef struct {
	int height;
	double Lrate, Rrate;
	int Lnbr, Rnbr;
} crystal_site;

int my_itoa(int val, char* buf);
double getRate(crystal_site *h, int i, int whichNbr, double K, int is_dH);
void write2file(char *name, double *array, int dim1, int dim2);

int main(int argc, char *argv[])
{
	long m;
	int i, j, k, delta, i2, whichNbr, num_ts_ctr;
	double z1, z2, eta, timetonext, t, ratesum;
	int siteupdatelist[4];
	const double PI = 3.14159265358979323846264338327950288;
	int is_dH;
	int L;
	//int scale_factor = L;
	int nsmpls;
	double K;
	int num_times;

	
	
	// parse seed from command line
	long long unsigned int sd;
	sd = atoi(argv[1]);
	init_genrand64(sd);
	char seed[20];
	my_itoa(sd, seed);

	
	// read in parameter values
	char parameters[100] = "./parameters.txt";

	FILE *fid;
	fid = fopen(parameters, "r");
	fscanf(fid, "%d %d %lf %d %d", &is_dH, &L, &K, &nsmpls, &num_times);
	printf("is_dH = %d,  L=%d, K=%f, nsmpls = %d, num_times =%d\n", is_dH, L, K, nsmpls, num_times);
	int scale_factor = L;
	//double T[num_times];
	double *T;
	T = (double *)malloc(num_times*sizeof(double));
	for (i = 0; i < num_times; i++) {
		//fscanf(fid, "\n%lf", &T[i]);
		fscanf(fid, " %lf", &T[i]);
		printf("%e\n", T[i]);
		T[i] *= pow(L, 4);
	}
	fclose(fid);

	
   
	//double hini[L];
	//double avh[(num_times + 1)*L]; //record avh at each time T_i including 0
	double *hini, *avh;
	hini = (double *)malloc(L*sizeof(double));
	avh = (double *)malloc((num_times+1)*L*sizeof(double));
	for (i = 0; i < L; i++) {
		hini[i] = scale_factor*sin(2 * PI*((double)i) / ((double)L));
	}
	 
	for (j = 0; j<(num_times+1)*L; j++) {
			avh[j] = 0.0;
	}

	//crystal_site h[L];
	crystal_site *h;
	h = (crystal_site *)malloc(L*sizeof(crystal_site));
	h[0].Lnbr = L - 1;
	h[0].Rnbr = 1;
	h[L - 1].Lnbr = L - 2;
	h[L - 1].Rnbr = 0;
	for (i = 1; i<L - 1; i++) {
		h[i].Lnbr = i - 1;
		h[i].Rnbr = i + 1;
	}


	clock_t start, end;
	double cpuTime;

	start = clock();
	for (m = 0; m<nsmpls; m++) {
		printf("on %dth sample\n", m);
		for (i = 0; i < L; i++) {		
			h[i].height = (int) floor(hini[i]);
			if (uniform64() < hini[i] - h[i].height)
				h[i].height++;		
			avh[i] += (double)h[i].height/((double) scale_factor*nsmpls);
		}
		ratesum = 0;
		for (i = 0; i<L; i++) {
			h[i].Rrate = getRate(h, i, 0, K, is_dH);
			h[i].Lrate = getRate(h, i, 1, K, is_dH);
			ratesum += (h[i].Lrate + h[i].Rrate);
		}

		t = 0;
		for (k = 0; k < num_times; k++) { 
			int num_ts_ctr = 0;
			printf("before, t = %lf, rate average = %lf\n", t,ratesum/L);
			//this while loop performs KMC over each time interval, i.e. from T[k-1] (or 0) to T[k]
			while (t < T[k]) {
				timetonext = -log(uniform64()) / ratesum;
				if (t + timetonext > T[k])
					break;
				num_ts_ctr++;

				eta = ratesum*uniform64();
				z2 = 0;
				i = -1;
				//different algorithms for different rates (to save time in the case of adatom rates)
				if(is_dH){
					while (z2<eta) {
						i++;
						if (i == L&& eta - z2 < 1e-5) {
						printf("warning! slight overflow. z2 = %f,eta=%f\n", z2, eta);
						i = L - 1;
						break;
						}
						whichNbr = 1;
						z2 += h[i].Lrate;
						if (z2 < eta) {
							whichNbr = 0;
							z2 += h[i].Rrate;
						}
					}
				}
				else{
					while(z2<eta){ 
						i++; 
						z2+= 2*h[i].Rrate; 
                    } 
					whichNbr=0;
					if(uniform64() > 0.5)
						whichNbr=1;
				}
				siteupdatelist[0] = i;
				if (whichNbr == 1) {
					siteupdatelist[1] = h[siteupdatelist[0]].Lnbr;
					siteupdatelist[2] = h[siteupdatelist[1]].Lnbr;
					siteupdatelist[3] = h[siteupdatelist[0]].Rnbr;
				}
				else {
					siteupdatelist[1] = h[siteupdatelist[0]].Rnbr;
					siteupdatelist[2] = h[siteupdatelist[0]].Lnbr;
					siteupdatelist[3] = h[siteupdatelist[1]].Rnbr;
				}

				h[siteupdatelist[0]].height--;
				h[siteupdatelist[1]].height++;

				for (j = 0; j< 4; j++) {
					i2 = siteupdatelist[j];
					h[i2].Rrate = getRate(h, i2, 0, K, is_dH);
					if(is_dH == 1)
						h[i2].Lrate = getRate(h, i2, 1, K, is_dH);
					else
						h[i2].Lrate = h[i2].Rrate;
				}
				ratesum = 0;
				for (i = 0; i < L; i++)
					ratesum += (h[i].Lrate + h[i].Rrate);
				t += timetonext;
			}
			printf("after, t = %lf, rate average = %lf, num ts = %d\n", t,ratesum/L, num_ts_ctr);

			
			for (i = 0; i < L; i++){ 
				avh[(k+1)*L + i] += (double)h[i].height / ((double)scale_factor*nsmpls);
			}
			
		}
	}

	end = clock();
	cpuTime = (end - start) / (CLOCKS_PER_SEC);
	printf("\n");
    printf("KMC CPU time: %g minutes\n", cpuTime / 60.0);
	
	char str_avh[100] = "./h.txt";
	write2file(str_avh, avh, num_times+1, L);

	return 0;
}

int my_itoa(int val, char* buf)
{
	const unsigned int radix = 10;

	char* p;
	unsigned int a;        //every digit
	int len;
	char* b;            //start of the digit char
	char temp;
	unsigned int u;

	p = buf;

	if (val < 0)
	{
		*p++ = '-';
		val = 0 - val;
	}
	u = (unsigned int)val;

	b = p;

	do
	{
		a = u % radix;
		u /= radix;

		*p++ = a + '0';

	} while (u > 0);

	len = (int)(p - buf);

	*p-- = 0;

	//swap
	do
	{
		temp = *p;
		*p = *b;
		*b = temp;
		--p;
		++b;

	} while (b < p);

	return len;
}

double getRate(crystal_site *h, int i, int whichNbr, double K, int is_dH) {
	double rate;
	if (is_dH) {
		int h0, hL, hLL, hR, hRR, z1, z2, z3;
		double dH;
		h0 = h[i].height;
		hL = h[h[i].Lnbr].height;
		hLL = h[h[h[i].Lnbr].Lnbr].height;
		hR = h[h[i].Rnbr].height;
		hRR = h[h[h[i].Rnbr].Rnbr].height;

		if (whichNbr == 0) { //jump right
			z1 = h0 - hL;
			z2 = hR - h0;
			z3 = hRR - hR;
			dH = 6 - 2 * z1 + 4 * z2 - 2 * z3;
		}
		else { //jump left
			z1 = hL - hLL;
			z2 = h0 - hL;
			z3 = hR - h0;
			dH = 6 + 2 * z1 - 4 * z2 + 2 * z3;
		}
		rate = exp(-0.5*K*dH);
	}
	else {
		int hL, hR, h0, z1, z2;
		double cord;
		hL = h[h[i].Lnbr].height;
		hR = h[h[i].Rnbr].height;
		h0 = h[i].height;
		z1 = h0 - hL; z2 = hR - h0;
		cord = (double)(z2 - z1 + 1);
		rate = 0.5*exp(-2 * K*cord);
	}
	return rate;
}

void write2file(char *name, double *array, int dim1, int dim2) {
	FILE *fid;
	int i, j;
	fid = fopen(name, "w");
	for (i = 0; i<dim1; i++) {
		for (j = 0; j < dim2; j++)
			fprintf(fid, "%20.14lf ", array[i*dim2 + j]);
		fprintf(fid, "\n");
	}
	fclose(fid);
}

