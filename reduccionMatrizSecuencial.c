#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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

//  Funcion del profe para valores randoms
double randFP(double min, double max) 
{   
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div); 
}

int main(int argc, char *argv[])
{
    //Variables
    float *M, *Maux;
    float *swap;
    int converge;
    int N;
    double timetick;
    double diferencia;
    float divCuatro, divSeis, divNueve;
    int nroIteraciones= 0;

    //  Inicializacion de variables y vector
    N = atoi(argv[1]);
    M = (float *) malloc (sizeof(float)*N*N);
    Maux = (float *) malloc (sizeof(float)*N*N);
    converge = 0;    
    divCuatro = 1.0/4.0;
    divSeis = 1.0/6.0;
    divNueve = 1.0/9.0;


    //Inicializar matriz
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            M[i*N+j] = randFP(0.0,1.0);
        }
    }

    //  Comenzar
    printf("Calculando... \n");
    timetick = dwalltime();
    while(!converge)
    {
        //  Calculo de todos los promedios   
        // CALCULO DE VERTICES

        // Primer vertice primera fila
        Maux[0] = ( M[0] + M[1] + M[N] + M[N+1] ) * divCuatro;
        
        // Ultimo vertice primera fila
        Maux[N-1] = ( M[N-2] + M[N-1] + M[N + N - 2] + M[N + N - 1] ) * divCuatro;
        
        // Primer vertice ultima fila
        Maux[(N-1)*N] = ( M[(N-2)*N] + M[(N-2)*N + 1] + M[(N-1)*N] + M[(N-1)*N + 1] ) * divCuatro; 
        
        // Ultimo vertice ultima fila
        Maux[(N-1)*N+ N-1] = ( M[(N-1)*N + N-1-1] + M[(N-1)*N + N-1] + M[(N-1-1)*N + N-1] + M[(N-1-1)*N + N-1-1] ) * divCuatro; 
        

        //  CALCULO DE MATRIZ 'INTERIOR' Y BANDAS LATERALES
        for(int i=1; i<N-1; i++)
        {
            //Primera y ultima fila
            Maux[i] =   (M[i-1] +        M[i] +      M[i+1] +
                        M[1*N + i-1] +    M[1*N + i] +  M[1*N + i+1]) * divSeis;
            
            Maux[(N-1)*N + i] = ( M[(N-1-1)*N + i-1] +    M[(N-1-1)*N + i] + M[(N-1-1)*N + i+1]
                            +   M[(N-1)*N + i-1] +      M[(N-1)*N + i] + M[(N-1)*N + i+1]) * divSeis;

            //Primera y ultima columna
            Maux[i*N] = ( M[(i-1)*N] +    M[(i-1)*N + 1] +
                        M[i*N] +        M[i*N + 1] +
                        M[(i+1)*N] +    M[(i+1)*N + 1] ) * divSeis;

            Maux[i*N + (N-1)] = ( M[i*N - 2] +    M[i*N-1] +
                                M[ i*N + (N-2)] +        M[i*N + (N-1)] +
                                M[(i+1)*N + (N-2)] +    M[(i+1)*N + (N-1)] ) * divSeis;

            for(int j=1; j<N-1; j++)
            {
                Maux[i*N+j] = ( M[(i-1)*N+ (j-1)] +     M[(i-1)*N+(j)] +        M[(i-1)*N+ (j+1)] 
                        +       M[(i)*N+ (j-1)] +       M[(i)*N+ j] +           M[(i)*N+  (j+1)]
                        +       M[(i+1)*N+ (j-1)] +     M[(i+1)*N+ (j)] +       M[(i+1)*N+ (j+1)]) * divNueve;
            } 
            
        }
        
        //  CALCULO DE BANDAS LATERALES
        // Primera y ultima fila
        /*for(int i=1; i<N-1; i++ )
        {
            //Primera y ultima fila
            Maux[i] =   (M[i-1] +        M[i] +      M[i+1] +
                        M[1*N + i-1] +    M[1*N + i] +  M[1*N + i+1]) * divSeis;
            
            Maux[(N-1)*N + i] = ( M[(N-1-1)*N + i-1] +    M[(N-1-1)*N + i] + M[(N-1-1)*N + i+1]
                            +   M[(N-1)*N + i-1] +      M[(N-1)*N + i] + M[(N-1)*N + i+1]) * divSeis;

            //Primera y ultima columna
            Maux[i*N] = ( M[(i-1)*N] +    M[(i-1)*N + 1] +
                        M[i*N] +        M[i*N + 1] +
                        M[(i+1)*N] +    M[(i+1)*N + 1] ) * divSeis;

            Maux[i*N + (N-1)] = ( M[(i-1)*N - 1] +    M[(i-1)*N] +
                                M[i*N - 1] +        M[i*N] +
                                M[(i+1)*N - 1] +    M[(i+1)*N] ) * divSeis;
        }*/
        
        //  Swapeo de matrices
        swap = Maux;
        Maux = M;
        M = swap;

        //  Chequear convergencia
        converge = 1;
        //printf("VECTOR V[0]= %.2f   ",V[0]);
        for (int i = 1; i < N*N; i++)
        {
            diferencia= M[0] - M[i];
            if ((diferencia > VALORPRECISIONP) || (diferencia < VALORPRECISIONN))
            {
                converge= 0;
                break;
            }
        }
        nroIteraciones++;
        /*for(int i=0; i<N ; i++)
        {
            for(int j=0; j<N; j++)
            {
                diferencia = M[0] - M[i*N + j];
                if( (diferencia > valorPresicion) || (diferencia < -valorPresicion))
                {
                    converge = 0;
                    i= N;
                    break;
                }
            }            
        }*/
    }

    for (int i = 0; i < N*N; i++)
    {
        printf("%f, ",M[i]);
    }
    printf("\nTiempo en segundos: %f y numero de iteraciones: %d\n", dwalltime() - timetick, nroIteraciones);
    
}