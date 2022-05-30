#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/time.h>
#define VALORPRECISIONP 0.01
#define VALORPRECISIONN -0.01

//Para calcular tiempo
double dwalltime()
{
	double sec;
	struct timeval tv;
	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}

//Funcion para obtener un valor random
double randFP(double min, double max) 
{ 
	double range = (max - min); 
	double div = RAND_MAX / range; 
	return min + (rand() / div); 
}

//Variables
float *M, *Maux;
float *swap;
int * convergioV;
int convergio;
int N,T;
float divCuatro, divSeis, divNueve;
int nroIteraciones;
pthread_barrier_t barrera;


void * funcion(void * arg)
{
	int tid= *(int *)arg;
	int inicio;
	int final;
	float comparacion;
	float primerValor;

	inicio= tid * N / T;
	final= inicio + N /T;
	if (tid==0)
	{
		inicio++;
	}
	if (tid==T-1)
	{
		final--;
	}
	while(!convergio)
	{
		//Primera y segunda esquina
		if (tid == 0)
		{
			Maux[0]= (M[0] + M[1] + M[N] + M[N+1]) * divCuatro;
			for (int i = 1; i < N-1; i++)
			{
				Maux[i] =   (M[i-1] +        M[i] +      M[i+1] +
                        M[N + i-1] +    M[N + i] +  M[N + i+1]) * divSeis;
			}
			Maux[N-1]= (M[N-2] + M[N-1] + M[2*N - 2] + M[2*N -1]) * divCuatro;
		}
		//Tercera y cuarta esquina
		if (tid == T-1)
		{
			Maux[N*(N-1)]= (M[N*(N-2)] + M[(N*(N-2))+1] + M[N*(N-1)] + M[(N*(N-1))+1]) * divCuatro;
			for (int i = 1; i < N-1; i++)
			{
				Maux[(N-1)*N + i] = ( M[(N-1-1)*N + i-1] +    M[(N-1-1)*N + i] + M[(N-1-1)*N + i+1]
                            +   M[(N-1)*N + i-1] +      M[(N-1)*N + i] + M[(N-1)*N + i+1]) * divSeis;
			}
			//Maux[(N*N)-1]= (M[((N-1)*N - 2)] + M[((N-1)*N -1)] + M[(N*N)-2] + M[(N*N)-1]) * divCuatro;
			Maux[(N-1)*N+ N-1] = ( M[(N-1)*N + N-1-1] + M[(N-1)*N + N-1] + M[(N-1-1)*N + N-1] + M[(N-1-1)*N + N-1-1] ) * divCuatro;
		}
		for (int i = inicio; i < final; i++)
		{
			Maux[i*N] = ( M[(i-1)*N] +    M[(i-1)*N + 1] +
                        M[i*N] +        M[i*N + 1] +
                        M[(i+1)*N] +    M[(i+1)*N + 1] ) * divSeis;

            Maux[i*N + (N-1)] = ( M[i*N - 2] +    M[i*N-1] +
                                M[ i*N + (N-2)] +        M[i*N + (N-1)] +
                                M[(i+1)*N + (N-2)] +    M[(i+1)*N + (N-1)] ) * divSeis;

			for (int j = 1; j < N-1; j++)
			{
				Maux[i*N+j] = ( M[(i-1)*N+ (j-1)] +     M[(i-1)*N+(j)] +        M[(i-1)*N+ (j+1)] 
                        +       M[(i)*N+ (j-1)] +       M[(i)*N+ j] +           M[(i)*N+  (j+1)]
                        +       M[(i+1)*N+ (j-1)] +     M[(i+1)*N+ (j)] +       M[(i+1)*N+ (j+1)]) * divNueve;
			}	
		}
		pthread_barrier_wait(&barrera);
		primerValor= Maux[0];
		convergioV[tid]= 1;
		for (int i = tid * N * N / T; i < tid * N * N / T + N * N / T; i++)
		{
			comparacion= primerValor - Maux[i];
			if ((comparacion > VALORPRECISIONP) || (comparacion < VALORPRECISIONN))
			{
				convergioV[tid]= 0;
				break;
			}
		}
		pthread_barrier_wait(&barrera);
		if (tid == 0)
		{
			swap= M;
			M= Maux;
			Maux= swap;
			convergio= 1;
			for (int i = 1; i < T; i++)
			{
				if (!convergioV[i])
				{
					convergio= 0;
					break;
				}
			}
			nroIteraciones++;
		}
		pthread_barrier_wait(&barrera);
	}


	pthread_exit(NULL);
}


int main(int argc, char const *argv[])
{
	T=atoi(argv[1]);
	N=atoi(argv[2]);
	pthread_t misPthread[T];
	int threads_id[T];
	double timetick;

	divCuatro= 1.0/4.0;
	divSeis= 1.0/6.0;
	divNueve= 1.0/9.0;
	convergio= 0;
	nroIteraciones= 0;

	M=(float *)malloc(sizeof(float)*N*N);
	Maux=(float *)malloc(sizeof(float)*N*N);
	convergioV=(int *)malloc(sizeof(int)*T);

	pthread_barrier_init(&barrera,NULL,T);

	for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            M[i*N+j] = randFP(0.0,1.0);
        }
    }

	timetick= dwalltime();
	for (int id = 0; id < T; id++)
	{
		threads_id[id]= id;
		pthread_create(&misPthread[id],NULL,&funcion,(void *)&threads_id[id]);
	}

	for (int i = 0; i < T; i++)
	{
		pthread_join(misPthread[i],NULL);
	}

	for (int i = 0; i < N*N; i++)
    {
        printf("%f, ",M[i]);
    }
	printf("\nTiempo en segundos %f segundos y numero de iteraciones %d\n",dwalltime() - timetick,nroIteraciones);

	free(M);
	free(Maux);
	free(convergioV);

	return 0;
}