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

int N,T;
float divDos,divTres;
float * V;
float * Vauxiliar;
float * swap;
int convergio;
int * convergioV;
pthread_barrier_t barrera;

void * funcion(void * arg)
{
	int tid=*(int *) arg;
	int inicio;
	int final;
	float comparacion;
	float primerValor;
	while(!convergio)
	{
		//Defino el inicio y fin de cada tarea segun su id
		if ((tid == 0) || (tid == T - 1))
		{
			if (tid == 0)
			{
				inicio= (tid * N / T) + 1;
			}else
			{
				inicio= (tid * N / T);
			}
			if (tid == T - 1)
			{
				final= (tid * N / T + N / T) - 1;
			}else
			{
				final= (tid * N / T + N / T);
			}
		}else
		{
			inicio= tid * N / T;
			final= tid * N / T + N / T;
		}
		//Valor 0 del vector
		if (tid == 0)
		{
			Vauxiliar[0]= (V[0] + V[1]) * divDos;
		}
		//Procesamiento general
		for (int i = inicio; i < final; i++)
		{
			Vauxiliar[i]= (V[i-1] + V[i] + V[i+1]) * divTres;
		}
		//Valor N - 1 del vector
		if (tid == T - 1)
		{
			Vauxiliar[N-1]= (V[N-1] + V[N-2]) * divDos;
		}
		pthread_barrier_wait(&barrera);
		primerValor= Vauxiliar[0];
		convergioV[tid]= 1;
		for (int i = tid * N / T; i < tid * N / T + N / T; i++)
		{
			comparacion= primerValor - Vauxiliar[i];
			if ((comparacion > VALORPRECISIONP) || (comparacion < VALORPRECISIONN))
			{
				convergioV[tid]= 0;
				break;				
			}
		}
		pthread_barrier_wait(&barrera);
		if (tid == 0)
		{
			swap= V;
			V= Vauxiliar;
			Vauxiliar= swap;
			convergio= 1;
			for (int i = 1; i < T; i++)
			{
				if (!convergioV[i])
				{
					convergio= 0;
					break;
				}
			}
			//printf("una iteracion");
		}
		pthread_barrier_wait(&barrera);
		/*while(tid % x == 0)
		{
			cantidad[tid]+= cantidad[tid+j];
			x *= 2;
			if(x < T)
			{
				pthread_barrier_wait(&barreras[tid/x]);
				j *=2;
			}else
			{
				break;
			}
		}*/
	}
	

	pthread_exit(NULL);
}

int main(int argc, char const *argv[])
{
	T= atoi(argv[1]);
	N= atoi(argv[2]);
	pthread_t misPthread[T];
	int threads_id[T];
	double timetick;

	divDos= 1.0/2.0;
	divTres= 1.0/3.0;
	convergio= 0;

	V=(float *)malloc(sizeof(float)*N);
	Vauxiliar=(float *)malloc(sizeof(float)*N);
	convergioV=(int *)malloc(sizeof(int)*T);

	for (int i = 0; i < N; i++)
	{
		V[i]= randFP(0.0,1.0);
	}

	
	pthread_barrier_init(&barrera,NULL,T);

	timetick= dwalltime();
	for (int id = 0; id < T; id++)
	{
		threads_id[id]= id;
		pthread_create(&misPthread[id],NULL,&funcion,(void *)&threads_id[id]);
	}

	for (int i = 0; i < T; ++i)
	{
		pthread_join(misPthread[i],NULL);
	}

	printf("Tiempo en segundos %f\n",dwalltime() - timetick);
	/*printf("Vector resultante:\n");
	for (int i = 0; i < N; i++)
	{
		printf("%f, ",V[i]);
	}*/

	free(V);

	return 0;
}