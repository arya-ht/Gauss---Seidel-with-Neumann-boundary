#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// created by Arya HajiTaheri
//#pragma warning(disable:4996)
#define tMax 10
#define dx (M_PI/10)
#define xMAX M_PI
#define x0 0
#define D 1.0
#define F 0.0
#define n 1.0

#define lambda .4
#define dt (lambda*dx*dx/D)
//function declaration
double U0(double);
double exact(double, double);
void Gauss_S_N();

int main(void)
{
	Gauss_S_N();
	return 0;
}

void Gauss_S_N()
{   // init
	int k;
	printf("Enter k:\n");
	scanf("%d", &k);
	int N = (int)((xMAX - x0) / dx), M = (int)(tMax / dt);
	double u[M + 1][N + 1], error = 0.0;


	int i, j, t;
	FILE *error_f, *estimate_f, *exact_f;
	estimate_f = fopen("Estimate.csv", "w");
	for (t = 0; t <= M; t++)
	{
		for (i = 1; i<N; i++)
		{
			u[t][i] = U0(i*dx);
		}
	}
	for (i = 0; i <= M; i++)
	{
		//u[i][0] = 0; // we don't know the left bound ... we use ghost nodes instead.
		u[i][N] = 0;
	}

	for (i = 0; i <= M; i++)
	{
		for (j = 0; j <= N; j++)
		{
			printf("%.8lf ", u[i][j]);
			fprintf(estimate_f, "%lf, ", u[i][j]);
		}
		printf("\n");
		fprintf(estimate_f, "\n");
	}
	printf("\n\n");
	exact_f = fopen("Exact.csv", "w");

	for (i = 0; i <= M; i++)
	{
		for (j = 0; j <= N; j++)
		{
			printf("%.8lf ", exact(j*dx, i*dt));
			fprintf(exact_f, "%lf, ", exact(j*dx, i*dt));
		}
		printf("\n");
		fprintf(exact_f, "\n");
	}
	error_f = fopen("Error.csv", "w");
	for (i = 0; i <= M; i++)
	{
		for (j = 0; j <= N; j++)
		{
			error = fabs((u[i][j] - exact(j*dx, i*dt)) / u[i][j]) * 100;
			if (u[i][j] == 0) // because error would be infinity (0/0)
			{
				error = 0;
			}

			fprintf(error_f, "%lf, ", error);
		}

		fprintf(error_f, "\n");
	}
	fclose(exact_f);
	fclose(estimate_f);
	fclose(error_f);
}
double U0(double x)
{
	return cos((n - .5)*x);
}
double exact(double x, double t)
{
	return exp(-1 * (n - .5)*(n - .5)*t)*cos((n - .5)*x);
}

