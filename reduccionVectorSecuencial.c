#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/time.h>
#define VALORPRECISIONP 0.01
#define VALORPRECISIONN -0.01

/***************************************
 FUNCION PARA CALCULAR TIEMPO
 ***************************************/
//Para calcular tiempo
double dwalltime()
{
	double sec;
	struct timeval tv;
	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}


/***************************************
 	FUNCION QUE RETORNA UN VALOR RANDOM
 ***************************************/
double randFP(double min, double max) 
{ 
	double range = (max - min); 
	double div = RAND_MAX / range; 
	return min + (rand() / div); 
}


/***************************************
			FUNCION MAIN
 ***************************************/
int main(int argc, char const *argv[])
{
	/* DECLARACION DE VARIABLES */
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


	/* ALOCACION DE MEMORIA DE LOS VECTORES */
	Vauxiliar=(float *)malloc(sizeof(float)*N);
	Vsecuencial=(float *)malloc(sizeof(float)*N);

	/* INICIALIZACION DEL VECTOR */
	for (int i = 0; i < N; i++)
	{
		Vsecuencial[i]= randFP(0.0,1.0);
	}

	printf("Calculando...\n");

	/* COMIENZA LA REDUCCION */
	int convergioSecuencial= 0;
	timetick= dwalltime();	//ARRANCA A CONTAR EL TIEMPO
	while(!convergioSecuencial)
	{
		//PROMEDIO DEL PRIMER VALOR
		Vauxiliar[0]= (Vsecuencial[0] + Vsecuencial[1]) *divDos;
		//PROMEDIO VALORES INTERMEDIOS
		for (int i = 1; i < N - 1; i++)
		{
			Vauxiliar[i]= (Vsecuencial[i-1] + Vsecuencial[i] + Vsecuencial[i+1]) * divTres;
		}
		//PROMEDIO ULTIMO VALOR
		Vauxiliar[N-1]= (Vsecuencial[N-1] + Vsecuencial[N-2]) *divDos;

		//SWAPEO DE VECTORES
		swap= Vsecuencial;
		Vsecuencial= Vauxiliar;
		Vauxiliar= swap;

		//CHEQUEO DE CONVERGENCIA DEL VECTOR
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

		//INCREMENTA NUMERO DE ITERACIONES
		nroIteraciones++;
	}

	//IMPRIME RESULTADO EN TIEMPO Y NUMERO DE ITERACIONES
	printf("REDUCCION DE VECTOR SECUENCIAL: Tiempo en segundos %f y numero de iteraciones %d\n",dwalltime() - timetick,nroIteraciones);


	//DESCOMENTAR SI SE QUIERE IMPRIMIR EL VECTOR AL QUE CONVERGIO
	/*printf("Vector resultante:\n");
	for (int i = 0; i < N; i++)
	{
		printf("%f, ",Vsecuencial[i]);
	}*/

	//LIBERACION DE MEMORIA
	free(Vsecuencial);
	free(Vauxiliar);

	return 0;
}