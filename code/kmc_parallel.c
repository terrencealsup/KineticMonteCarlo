/*******************************************************************************
File: kmc_parallel.c

Authors: Anya Katsevich
         Terrence Alsup

High Performance Computing 2019, Final Project
Date: May 17, 2019

Run 1d Kinetic Monte Carlo in parallel using an MPI implementation.
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
#include "mt64.h"

// Random number generation.
#ifndef uniform64
#define uniform64() genrand64_real3()
#endif

const double PI = 3.14159265358979323846264338327950288;

// A crystal site.
typedef struct {
	int height;
	double rate;
} crystal_site;


int my_itoa(int val, char* buf);
void myprint(long N, char* name, crystal_site* h);
double getRate(crystal_site *h, int i, double K);
void initialize_lattice(crystal_site* h, int L_loc, int rank, int L, double K, int P, double* R_locs, MPI_Comm comm);
int KMC_section(int section_num, crystal_site *h, double t_stop, double* R_locs, double K, int L_loc, int rank, int itr_num);
int drawSite(crystal_site* h, double R_loc, int section_num, int L_loc);
void rateMax(crystal_site* h, int rank, int L_loc);




/*
Initialize the lattice for the block of size L_loc = L/P
*/
void initialize_lattice(crystal_site* h, int L_loc, int rank, int L, double K, int P, double* R_locs, MPI_Comm comm) {

	MPI_Status status1;

	//initialize heights
	for (int i = 1; i < L_loc + 1; i++) {
		int j = L_loc * rank + i;
		double x = ((double) j) / L;

		/*
		if (j > 0 && j < L/2) {
			//double heightVal = exp(8 - 1/x - 1/(0.5 - x));
			h[i].height = (int) floor(0.1 * L * heightVal);
		} else {
			h[i].height = 0;
		}
		*/

		double heightVal = L*sin(2 * PI* x);
		h[i].height = (int)floor(heightVal);
		if (uniform64() < heightVal - h[i].height)
			h[i].height++;
	}



	MPI_Send(&(h[1].height), 1, MPI_INT, (rank - 1 + P) % P, 0, comm); // send h[1] to the left
	MPI_Send(&(h[L_loc].height), 1, MPI_INT, (rank + 1) % P, 1, comm); // send h[L_loc] to the right

	MPI_Recv(&(h[L_loc + 1].height), 1, MPI_INT, (rank + 1) % P, 0, comm, &status1); //receive to L_loc + 1 from the right
	MPI_Recv(&(h[0].height), 1, MPI_INT, (rank - 1 + P) % P, 1, comm, &status1); //receive to 0 from the left



	for (int i = 1; i <= L_loc / 2; i++) {
		h[i].rate = getRate(h, i, K);
		R_locs[0] += 2 * h[i].rate;
	}
	for (int i = L_loc / 2+1; i <= L_loc; i++) {
		h[i].rate = getRate(h, i, K);
		R_locs[1] += 2 * h[i].rate;
	}


}


double getRate(crystal_site *h, int i, double K) {
	double rate;

	int hL, hR, h0, z1, z2;
	double cord;
	hL = h[i - 1].height;
	hR = h[i + 1].height;
	h0 = h[i].height;
	z1 = h0 - hL; z2 = hR - h0;
	cord = (double)(z2 - z1 + 1);
	rate = 0.5*exp(-2 * K*cord);

	return rate;
}



int KMC_section(int section_num, crystal_site *h, double t_stop, double* R_locs, double K, int L_loc, int rank, int itr_num) {

	//printf("Rank %d: at beginning of KMC section %d: R_loc%d = %f (passed in)\n", rank, section_num, section_num, R_locs[section_num]);


	double t = 0;
	double timetonext;
	int i, j, i2, whichNbr, siteupdatelist[4];
	int numEvents = 0;
	while (t < t_stop) {
		//numEvents = 0;
		timetonext = -log(uniform64()) / R_locs[section_num]; //draw from Exp(R_loc)
		if (t + timetonext > t_stop)
			break;

		//draw from distribution{ j w.p.r_j / R_loc, j = 0,...,L_loc - 1 }

		i = drawSite(h, R_locs[section_num], section_num, L_loc);

		if (i ==0 || i == L_loc + 1) printf("BAD!!\n");
		whichNbr = 0;
		if (uniform64() > 0.5)
			whichNbr = 1;

			numEvents++;
			siteupdatelist[0] = i;


			if (whichNbr == 1) {
				siteupdatelist[1] = i - 1;
				siteupdatelist[2] = i - 2;
				siteupdatelist[3] = i + 1;
			}
			else {
				siteupdatelist[1] = i + 1;
				siteupdatelist[2] = i - 1;
				siteupdatelist[3] = i + 2;
			}
			h[siteupdatelist[0]].height--;
			h[siteupdatelist[1]].height++;


			for (j = 0; j< 4; j++) {
				i2 = siteupdatelist[j];
				///printf("processor %d: i2 = %d\n", rank, i2);
				if (i2 > 0 && i2 < L_loc + 1) {
					double oldRate = 2 * h[i2].rate;
					h[i2].rate = getRate(h, i2, K);
          if(i2 <= L_loc/2){
                 R_locs[0] -= oldRate;
                 R_locs[0] += 2*h[i2].rate;
          }
          else{
                R_locs[1] -= oldRate;
                R_locs[1] += 2*h[i2].rate;
          }

          }

       }

		t += timetonext;

	}

	return numEvents;
  // printf("Rank %d: At end of KMC section %d: R_loc%d = %f\n\n", rank, section_num, section_num,R_locs[section_num]);


}

void rateMax(crystal_site* h, int rank, int L_loc) {

	int max_idx = 1;
	double max_rate = 0;
	for (int k = 1; k <= L_loc; k++) {
		if (h[k].rate >= max_rate) {
			max_rate = h[k].rate;
			max_idx = k;
		}
	}
	printf("Rank %d: new max rate = %f at index %d, heights = %d %d %d\n", rank, max_rate, max_idx, h[max_idx - 1].height, h[max_idx].height, h[max_idx + 1].height);


}

int drawSite(crystal_site* h, double R_loc, int section_num, int L_loc) {


	double eta = R_loc*uniform64();

	double z2 = 0;

	int i;
	if (section_num == 0)
		i = 0;
	else
		i = L_loc / 2 ;

	while (z2<eta) {
		i++;
		z2 += 2 * h[i].rate;
	}

	return i;
}

int main(int argc, char * argv[]) {

	int rank; // Which processor is executing currently.
	int P;		// The number of blocks/processors.

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &P);

	// parse seed from command line
	long long unsigned int sd;
	sd = atoi(argv[1]);
	init_genrand64(sd + rank);
	char seed[20];
	my_itoa(sd, seed);


	int L, L_loc, numEvents;
	double R_max0, R_max1, c, n, t_stop, K, Tfinal, cpuTime;
	clock_t start, end;
	char parameters[100] = "./parameters.txt";

	FILE *fid;
	fid = fopen(parameters, "r");
	fscanf(fid, "%d %lf %lf %lf", &L, &K, &Tfinal, &c);

	fclose(fid);


	Tfinal *= pow(L, 4);
	if (rank == 0)
		printf("L=%d, K=%f, T=%e, c = %f\n", L, K, Tfinal, c);
	n = c * L;

	// Initialize all of the blocks.
	L_loc = L / P; // Get the size of the block.

	crystal_site* h = (crystal_site *)malloc((L_loc + 2) * sizeof(crystal_site));
	double* R_locs = (double *)calloc(sizeof(double), 2);
	initialize_lattice(h, L_loc, rank, L, K, P, R_locs, MPI_COMM_WORLD);

	printf("Rank %d, initial R_loc0 = %f, initial R_loc1=%f\n\n", rank, R_locs[0], R_locs[1]);
	char str_hInit[100];
	snprintf(str_hInit, 100, "hInit%02d.txt", rank);

	myprint(L_loc, str_hInit, h);

	start = clock();
	double t = 0;  // Keep track of the time.
	int num_itr = 0;

	int total_events = 0;

	while (t < Tfinal) {
	//while (num_itr < 40) {


		MPI_Allreduce(&(R_locs[0]), &R_max0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		t_stop = n / R_max0;
			if (t + t_stop > Tfinal)
			  break;


   total_events += KMC_section(0, h, t_stop, R_locs, K, L_loc, rank, num_itr);


		t += t_stop;



		R_locs[1] -= 2 * h[L_loc].rate;
		R_locs[1] -= 2 * h[L_loc - 1].rate;
 		//COMMUNICATE BETWEEN PROCESSORS
		MPI_Comm comm = MPI_COMM_WORLD;
		MPI_Status status1;

		MPI_Send(&(h[0].height), 1, MPI_INT, (rank - 1 + P) % P, 123, comm);
		MPI_Send(&(h[1].height), 1, MPI_INT, (rank - 1 + P) % P, 124, comm);

		MPI_Recv(&(h[L_loc].height), 1, MPI_INT, (rank + 1) % P, 123, comm, &status1);
		MPI_Recv(&(h[L_loc + 1].height), 1, MPI_INT, (rank + 1) % P, 124, comm, &status1);


/*
			int sum = 0;
			for (int i = 1; i < L_loc + 1; i++) {
				sum += h[i].height;
			}
			printf("sum of heights = %d\n", sum);
*/


		h[L_loc].rate = getRate(h, L_loc, K);
		h[L_loc - 1].rate = getRate(h, L_loc - 1, K);



		R_locs[1] += 2 * h[L_loc].rate;
		R_locs[1] += 2 * h[L_loc - 1].rate;


		MPI_Allreduce(&(R_locs[1]), &R_max1, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);

		t_stop = n / R_max1;


		total_events += KMC_section(1, h, t_stop, R_locs, K, L_loc, rank, num_itr);



		t += t_stop;


		R_locs[0] -= 2 * h[1].rate;
		R_locs[0] -= 2 * h[2].rate;



		//MPI_Barrier(MPI_COMM_WORLD);
		//COMMUNICATE BETWEEN PROCESSORS

		MPI_Send(&(h[L_loc].height), 1, MPI_INT, (rank + 1) % P, 123, comm);
		MPI_Send(&(h[L_loc + 1].height), 1, MPI_INT, (rank + 1) % P, 124, comm);

		MPI_Recv(&(h[0].height), 1, MPI_INT, (rank - 1 + P) % P, 123, comm, &status1);
		MPI_Recv(&(h[1].height), 1, MPI_INT, (rank - 1 + P) % P, 124, comm, &status1);


/*
			sum = 0;
			for (int i = 1; i < L_loc + 1; i++) {
				sum += h[i].height;
			}
			printf("sum of heights = %d\n", sum);
			*/

		h[1].rate = getRate(h, 1, K);
		h[2].rate = getRate(h, 2, K);


		R_locs[0] += 2 * h[1].rate;
		R_locs[0] += 2 * h[2].rate;


		if (rank == 0 && num_itr%1000 == 0)
			printf("num_itr = %d, R_max1 = %f, new t = %f\n",  num_itr, R_max1, t);


		num_itr++;
	}

	end = clock();
	cpuTime = (end - start) / (CLOCKS_PER_SEC);
	if (rank == 0) {
		printf("\n");
		printf("KMC CPU time: %g minutes\n", cpuTime / 60.0);
	}

	int all_processor_events;
	MPI_Reduce(&total_events, &all_processor_events, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if (rank == 0) {
		double avg = ((double) all_processor_events) / P;
		printf("Average number of events per processor per time = %f\n", avg/Tfinal);
	}



	char str_hFinal[100];
	snprintf(str_hFinal, 100, "hFinal%02d.txt", rank);
	myprint(L_loc, str_hFinal, h);

	MPI_Finalize();

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

void myprint(long N, char* name, crystal_site* h) {
	FILE *fid;

  fid = fopen(name, "w+");

  for (long k = 0; k < N; k++)
		fprintf(fid, "%d ", h[k].height);
  fclose(fid);

}
