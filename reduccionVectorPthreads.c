#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/time.h>
#define VALORPRECISIONP 0.01 	//VALOR DE PRECISION POSITIVO
#define VALORPRECISIONN -0.01 	//VALOR DE PRECISION NEGATIVO

/***************************************
 FUNCION PARA CALCULAR TIEMPO
 ***************************************/
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

//DECLARACION DE VARIABLES
int N,T;
float divDos,divTres;
float * V;
float * Vauxiliar;
float * swap;
int convergio;
int * convergioV;
pthread_barrier_t barrera;
int nroIteraciones;

//FUNCION QUE VA A REALIZAR CADA TAREA
void * funcion(void * arg)
{
	int tid=*(int *) arg;	//RECUPERA SU ID
	int inicio;
	int final;
	int inicioC;
	int finalC;
	float comparacion;
	float primerValor;

	//INICIO Y FINAL PARA EL FOR QUE CHEQUEA CONVERGENCIA
	inicioC= tid * N / T;
	finalC= inicioC + N / T;

	//DEFINO EL INICIO Y FIN DE CADA TAREA SEGUN SU ID
	if ((tid == 0) || (tid == T - 1))
	{
		//LA TAREA 0 VA A COMENZAR EN EL VALOR 1 DEL VECTOR YA QUE EL PRIMER ELEMENTO LO PROCESA APARTE
		if (tid == 0)
		{
			inicio= (tid * N / T) + 1;
			final= (inicio - 1 + N / T);
			//LA TAREA 0 PROCESA EL PRIMER VALOR
			Vauxiliar[0]= (V[0] + V[1]) * divDos;
		}else	//LA ULTIMA TAREA VA A PROCESAR HASTA EL ANTEULTIMO VALOR YA QUE EL ULTIMO VALOR LO PROCESA APARTE
		{
			inicio= (tid * N / T);
			final= (inicio + N / T) - 1;
		}
	}else	//LAS DEMAS TAREAS DEFINEN SU INICIO Y FINAL SEGUN SU ID
	{
		inicio= inicioC;
		final= finalC;
	}

	//TODOS LOS PROCESOS TIENE QUE ESPERAR A LA PRIMER TAREA QUE TERMINE DE PROCESAR LA PRIMERA POSICION
	pthread_barrier_wait(&barrera);

	//MIENTRAS NO CONVERGA
	while(!convergio)
	{
		//PROCESAMIENTO GENERAL
		for (int i = inicio; i < final; i++)
		{
			Vauxiliar[i]= (V[i-1] + V[i] + V[i+1]) * divTres;
		}
		//LA ULTIMA TAREA PROCESA EL ULTIMO VALOR
		if (tid == T - 1)
		{
			Vauxiliar[N-1]= (V[N-1] + V[N-2]) * divDos;
		}

		//CADA TAREA SE GUARDA EL PRIMER VALOR CON EL CUAL VAN A CHEQUEAR SI SU PARTE DEL VECTOR CONVERGE
		primerValor= Vauxiliar[0];

		convergioV[tid]= 1;
		for (int i = inicioC; i < finalC; i++)
		{
			comparacion= primerValor - Vauxiliar[i];
			//SI LA COMPARACION DA MAYOR QUE EL VALOR DE PRECISION POSITIVO O DA MENOR QUE EL VALOR DE PRECISION NEGATIVO SIGNIFICA QUE EL VECTOR NO CONVERGIO
			if ((comparacion > VALORPRECISIONP) || (comparacion < VALORPRECISIONN))
			{
				convergioV[tid]= 0;
				break;				
			}
		}
		//TODAS LAS TAREAS ESPERAN A LAS DEMAS QUE TERMINEN DE CHEQUEAR LA CONVERGENCIA DE SU PARTE
		pthread_barrier_wait(&barrera);
		//LA TAREA 0 REALIZA EL SWAPEO DE LOS VECTORES Y VERIFICA SI TODAS LAS TAREAS CONVERGIERON
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
			nroIteraciones++;
			if(!convergio)
			{
				//LA TAREA 0 PROCESA EL PRIMER VALOR
				Vauxiliar[0]= (V[0] + V[1]) * divDos;
			}
		}
		//LAS DEMAS TAREAS ESPERAN QUE LA TAREA 0 HAGA EL SWAPEO Y CHEQUEE LA CONVERGENCIA TOTAL
		pthread_barrier_wait(&barrera);
	}
	

	pthread_exit(NULL);
}


/***************************************
			FUNCION MAIN
 ***************************************/
int main(int argc, char const *argv[])
{
	//DECLARACION Y INICIALIZACION DE VARIABLES
	T= atoi(argv[1]);
	N= atoi(argv[2]);
	pthread_t misPthread[T];
	int threads_id[T];
	double timetick;
	double tiempoEnSeg;

	//VARIABLES PARA REALIZAR LA DIVISION
	divDos= 1.0/2.0;
	divTres= 1.0/3.0;
	//VARIABLE COMPARTIDA QUE CONTROLA SI EL VECTOR CONVERGIO
	convergio= 0;
	//VARIABLE QUE CUENTA LA CANTIDAD DE ITERACIONES
	nroIteraciones= 0;

	//ALOCACION DE MEMORIA PARA LOS VECTORES
	V=(float *)malloc(sizeof(float)*N);
	Vauxiliar=(float *)malloc(sizeof(float)*N);
	convergioV=(int *)malloc(sizeof(int)*T);

	//INICIALIZACION DEL VECTOR
	for (int i = 0; i < N; i++)
	{
		V[i]= randFP(0.0,1.0);
	}

	//INICIALIZACION DE LA BARRERA
	pthread_barrier_init(&barrera,NULL,T);

	//ARRANCA A CONTAR EL TIEMPO Y CREO LAS TAREAS PASANDOLES SU ID
	timetick= dwalltime();
	for (int id = 0; id < T; id++)
	{
		threads_id[id]= id;
		pthread_create(&misPthread[id],NULL,&funcion,(void *)&threads_id[id]);
	}

	//REALIZA EL JOIN DE CADA TAREA CREADA
	for (int i = 0; i < T; ++i)
	{
		pthread_join(misPthread[i],NULL);
	}

	tiempoEnSeg= dwalltime() - timetick;
	//IMPRIME EL TIEMPO EN SEGUNDOS Y LA CANTIDAD DE ITERACIONES
	printf("Tiempo en segundos %f y numero de iteraciones %d\n",tiempoEnSeg,nroIteraciones);

	//DESCOMENTAR SI SE QUIERE IMPRIMIR EL VECTOR RESULTADO
	/*printf("Vector resultante:\n");
	for (int i = 0; i < N; i++)
	{
		printf("%f, ",V[i]);
	}*/

	//LIBERA LAS VARIABLES ALOCADAS
	free(V);
	free(Vauxiliar);
	free(convergioV);

	return 0;
}