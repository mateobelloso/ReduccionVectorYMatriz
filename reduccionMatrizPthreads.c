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
float *M, *Maux;
float *swap;
int * convergioV;
int convergio;
int N,T;
float divCuatro, divSeis, divNueve;
int nroIteraciones;
pthread_barrier_t barrera;

//FUNCION QUE EJECUTA CADA TAREA
void * funcion(void * arg)
{
	int tid= *(int *)arg;
	int inicio;
	int final;
	float comparacion;
	float primerValor;

	//INICIO Y FINAL DE FILA A PROCESAR DE CADA TAREA
	inicio= tid * N / T;
	final= inicio + N /T;
	if (tid==0)	//LA TAREA 0 VA A COMENZAR A PROCESAR 1 FILA DESPUES YA QUE LA PRIMERA LA HACE APARTE
	{
		inicio++;
	}
	if (tid==T-1)	//LA TAREA ULTIMA VA A TERMINAR DE PROCESAR 1 FILA ANTES A QUE LA ULTIMA FILA LA HACE APARTE
	{
		final--;
	}
	//MIENTRA NO CONVERGA
	while(!convergio)
	{
		//PRIMERA Y SEGUNDA ESQUINA Y BANDA LATERAL SUPERIOR
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

		//PROCESAMIENTO DE VALORES INTERMEDIOS
		for (int i = inicio; i < final; i++)
		{
			//BANDA LATERAL IZQUIERDA
			Maux[i*N] = ( M[(i-1)*N] +    M[(i-1)*N + 1] +
                        M[i*N] +        M[i*N + 1] +
                        M[(i+1)*N] +    M[(i+1)*N + 1] ) * divSeis;

            //VALORES INTERMEDIOS
			for (int j = 1; j < N-1; j++)
			{
				Maux[i*N+j] = ( M[(i-1)*N+ (j-1)] +     M[(i-1)*N+(j)] +        M[(i-1)*N+ (j+1)] 
                        +       M[(i)*N+ (j-1)] +       M[(i)*N+ j] +           M[(i)*N+  (j+1)]
                        +       M[(i+1)*N+ (j-1)] +     M[(i+1)*N+ (j)] +       M[(i+1)*N+ (j+1)]) * divNueve;
			}	

			//BANDA LATERAL DERECHA
            Maux[i*N + (N-1)] = ( M[i*N - 2] +    M[i*N-1] +
                                M[ i*N + (N-2)] +        M[i*N + (N-1)] +
                                M[(i+1)*N + (N-2)] +    M[(i+1)*N + (N-1)] ) * divSeis;
		}
		//TERCERA Y CUARTA ESQUINA Y BANDA LATERAL INFERIOR
		if (tid == T-1)
		{
			Maux[N*(N-1)]= (M[N*(N-2)] + M[(N*(N-2))+1] + M[N*(N-1)] + M[(N*(N-1))+1]) * divCuatro;
			for (int i = 1; i < N-1; i++)
			{
				Maux[(N-1)*N + i] = ( M[(N-1-1)*N + i-1] +    M[(N-1-1)*N + i] + M[(N-1-1)*N + i+1]
                            +   M[(N-1)*N + i-1] +      M[(N-1)*N + i] + M[(N-1)*N + i+1]) * divSeis;
			}
			Maux[(N-1)*N+ N-1] = ( M[(N-1)*N + N-1-1] + M[(N-1)*N + N-1] + M[(N-1-1)*N + N-1] + M[(N-1-1)*N + N-1-1] ) * divCuatro;
		}
		//ESPERO QUE TODAS LAS TAREAS TERMINEN SU PROCESAMIENTO YA QUE NECESITO EL RESULTADO DE LA PRIMERA POSICION PARA CHEQUEAR CONVERGENCIA
		pthread_barrier_wait(&barrera);

		//ME GUARDO EL RESULTADO DEL PRIMER VALOR
		primerValor= Maux[0];

		//CHEQUEO SI MI PARTE CONVERGIO
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
		//ESPERO QUE TODAS LAS TAREAS TERMINEN DE CHEQUEAR LA CONVERGENCIA DE SU PARTE
		pthread_barrier_wait(&barrera);

		//LA TAREA 0 REALIZA EL SWAPEO Y CHEQUEA SI TODAS LAS TAREAS CONVERGIERON
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
		//LAS DEMAS TAREAS ESPERAN A LA TAREA 0 QUE SWAPE LAS MATRICES Y EL CHEQUEO DE CONVERGENCIA TOTAL
		pthread_barrier_wait(&barrera);
	}


	pthread_exit(NULL);
}


/***************************************
			FUNCION MAIN
 ***************************************/
int main(int argc, char const *argv[])
{
	T=atoi(argv[1]);
	N=atoi(argv[2]);
	pthread_t misPthread[T];
	int threads_id[T];
	double timetick;
	double tiempoEnSeg;

	divCuatro= 1.0/4.0;
	divSeis= 1.0/6.0;
	divNueve= 1.0/9.0;
	convergio= 0;
	nroIteraciones= 0;

	//ALOCACION DE MEMORIA PARA LAS MATRICES Y EL VECTOR DE CONVERGENCIA
	M=(float *)malloc(sizeof(float)*N*N);
	Maux=(float *)malloc(sizeof(float)*N*N);
	convergioV=(int *)malloc(sizeof(int)*T);

	//INICIALIZACION DE LA BARRERA
	pthread_barrier_init(&barrera,NULL,T);

	//INICIALIZACION DE LA MATRIZ
	for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            M[i*N+j] = randFP(0.0,1.0);
        }
    }

    //ARRANCO A CONTAR EL TIEMPO, CREO LAS TAREAS Y LES PASO SU ID Y ARRANCAN A PROCESAR
	timetick= dwalltime();
	for (int id = 0; id < T; id++)
	{
		threads_id[id]= id;
		pthread_create(&misPthread[id],NULL,&funcion,(void *)&threads_id[id]);
	}

	//JOIN DE CADA TAREA
	for (int i = 0; i < T; i++)
	{
		pthread_join(misPthread[i],NULL);
	}

	//PARO DE CONTAR EL TIEMPO SE TERMINO EL PROCESAMIENTO
	tiempoEnSeg= dwalltime() - timetick;
	//DESCOMENTAR SI SE QUIERE IMPRIMIR EL RESULTADO

	/*
	printf("Resultado:\n");
	for (int i = 0; i < N; i++)
    {
        for(int j=0; j<N; j++)
        {
            printf("%f, ",M[i*N+j]);
        }
        printf("\n");
    }*/
	printf("\nTiempo en segundos %f segundos y numero de iteraciones %d\n",tiempoEnSeg,nroIteraciones);

	free(M);
	free(Maux);
	free(convergioV);

	return 0;
}