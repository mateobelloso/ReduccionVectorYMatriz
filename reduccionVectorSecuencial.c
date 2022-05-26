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

int main(int argc, char const *argv[])
{
	int N= atoi(argv[1]);
	double timetick;
	float divDos= 1.0/2.0;
	float divTres= 1.0/3.0;
	float * Vauxiliar;
	float * Vsecuencial;
	float * swap;
	float primerValor;
	float comparacion;
	int nroIteraciones= 0;

	Vauxiliar=(float *)malloc(sizeof(float)*N);
	Vsecuencial=(float *)malloc(sizeof(float)*N);

	for (int i = 0; i < N; i++)
	{
		Vsecuencial[i]= randFP(0.0,1.0);
	}

	printf("Calculando...\n");

	//SECUENCIAL
	int convergioSecuencial= 0;
	timetick= dwalltime();
	while(!convergioSecuencial)
	{
		//Procesamiento
		Vauxiliar[0]= (Vsecuencial[0] + Vsecuencial[1]) >> 1;
		for (int i = 1; i < N - 1; i++)
		{
			Vauxiliar[i]= (Vsecuencial[i-1] + Vsecuencial[i] + Vsecuencial[i+1]) * divTres;
		}
		Vauxiliar[N-1]= (Vsecuencial[N-1] + Vsecuencial[N-2]) >> 1;
		//Swapea los vectores
		swap= Vsecuencial;
		Vsecuencial= Vauxiliar;
		Vauxiliar= swap;
		//Evalua si convergio el vector
		convergioSecuencial= 1;
		primerValor= Vsecuencial[0];
		for (int i = 0; i < N; i++)
		{
			comparacion= primerValor - Vsecuencial[i];
			if ((comparacion > VALORPRECISIONP) || (comparacion < VALORPRECISIONN))
			{
				convergioSecuencial= 0;
				break;				
			}
		}
		nroIteraciones++;
	}

	printf("REDUCCION DE VECTOR SECUENCIAL: Tiempo en segundos %f y numero de iteraciones %d\n",dwalltime() - timetick,nroIteraciones);

	/*printf("Vector resultante:\n");
	for (int i = 0; i < N; i++)
	{
		printf("%f, ",Vsecuencial[i]);
	}*/

	free(Vsecuencial);

	return 0;
}