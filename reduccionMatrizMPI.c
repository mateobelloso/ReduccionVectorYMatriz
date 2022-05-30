#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include<mpi.h>
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


int miID;
int nrProcesos;
int N;
int convergio;
int convergioaux;
float *MRec;
float *M;
float *Maux;
float *filaVecAbajo, *filaVecArriba;
float primerValor;
float divCuatro= 1.0/4.0;
float divSeis= 1.0/6.0;
float divNueve= 1.0/9.0;
long int cantidad;

void fProcesoRoot()
{
    float vecinoDer;
    float comparacion;
    double timetick;
    int nroIteraciones= 0;

    int inicio= 1;
    int final= (N / nrProcesos) - 1; 
    
    M=(float *)malloc(sizeof(float)*N*N);
    //Inicializar matriz
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            M[i*N+j] = randFP(0.0,1.0);
        }
    }
    
    //primerValor= Vaux[0];
    /*for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            printf("%f, ",M[i*N+j]);
        }
        printf("\n");
    }*/

    timetick= dwalltime();
    while(!convergio)
    {
        MPI_Scatter(M,N*N/nrProcesos,MPI_FLOAT,MRec,N*N/nrProcesos,MPI_FLOAT,0,MPI_COMM_WORLD);
        MPI_Send(&MRec[N*N/nrProcesos-N],N,MPI_FLOAT,miID+1,99,MPI_COMM_WORLD);
        MPI_Recv(filaVecAbajo,N,MPI_FLOAT,miID+1,99,MPI_COMM_WORLD,MPI_STATUS_IGNORE);



        /*printf("Soy el proceso %d y me llego esta fila:\n",miID);
        for (int i = 0; i < N; i++)
        {
            printf("%f, ",filaVecAbajo[i]);
        }*/

        //printf("llegue vecinoDer: %f",vecinoDer);

        Maux[0]= (MRec[0] + MRec[1] + MRec[N] + MRec[N+1]) * divCuatro;
        for (int i = 1; i < N-1; i++)
        {
            Maux[i] =   (MRec[i-1] +        MRec[i] +      MRec[i+1] +
                    MRec[N + i-1] +    MRec[N + i] +  MRec[N + i+1]) * divSeis;
        }
        Maux[N-1]= (MRec[N-2] + MRec[N-1] + MRec[2*N - 2] + MRec[2*N -1]) * divCuatro;

        for (int i = inicio; i < final; i++)
        {
            Maux[i*N] = ( MRec[(i-1)*N] +    MRec[(i-1)*N + 1] +
                        MRec[i*N] +        MRec[i*N + 1] +
                        MRec[(i+1)*N] +    MRec[(i+1)*N + 1] ) * divSeis;

            Maux[i*N + (N-1)] = ( MRec[i*N - 2] +    MRec[i*N-1] +
                                MRec[ i*N + (N-2)] +        MRec[i*N + (N-1)] +
                                MRec[(i+1)*N + (N-2)] +    MRec[(i+1)*N + (N-1)] ) * divSeis;

            for (int j = 1; j < N-1; j++)
            {
                Maux[i*N+j] = ( MRec[(i-1)*N+ (j-1)] +     MRec[(i-1)*N+(j)] +        MRec[(i-1)*N+ (j+1)] 
                        +       MRec[(i)*N+ (j-1)] +       MRec[(i)*N+ j] +           MRec[(i)*N+  (j+1)]
                        +       MRec[(i+1)*N+ (j-1)] +     MRec[(i+1)*N+ (j)] +       MRec[(i+1)*N+ (j+1)]) * divNueve;
            }   
        }

        //Ultima fila de la tarea 0 SACAR EL FOR
        for (int i = final; i < final+1; i++)
        {
            Maux[i*N] = ( MRec[(i-1)*N] +    MRec[(i-1)*N + 1] +
                        MRec[i*N] +        MRec[i*N + 1] +
                        filaVecAbajo[0] +    filaVecAbajo[1] ) * divSeis;

            Maux[i*N + (N-1)] = ( MRec[i*N - 2] +    MRec[i*N-1] +
                                MRec[ i*N + (N-2)] +        MRec[i*N + (N-1)] +
                                filaVecAbajo[N-2] +    filaVecAbajo[N-1] ) * divSeis;
            for (int j = 1; j < N-1; j++)
            {
                Maux[i*N+j] = ( MRec[(i-1)*N+ (j-1)] +     MRec[(i-1)*N+(j)] +        MRec[(i-1)*N+ (j+1)] 
                        +       MRec[(i)*N+ (j-1)] +       MRec[(i)*N+ j] +           MRec[(i)*N+  (j+1)]
                        +       filaVecAbajo[j-1] +     filaVecAbajo[j] +       filaVecAbajo[j+1]) * divNueve;
            }
        }

        MPI_Bcast(&Maux[0],1,MPI_FLOAT,0,MPI_COMM_WORLD);

        convergioaux=1;
        for (int i = 0; i < N * N / nrProcesos; i++)
        {
            comparacion= Maux[0] - Maux[i];
            if ((comparacion > VALORPRECISIONP) || (comparacion < VALORPRECISIONN))
            {
                convergioaux= 0;
                break;
            }
        }

        MPI_Allreduce(&convergioaux,&convergio,1,MPI_INT,MPI_LAND,MPI_COMM_WORLD);

        MPI_Gather(Maux,N*N/nrProcesos,MPI_FLOAT,M,N*N/nrProcesos,MPI_FLOAT,0,MPI_COMM_WORLD);

        nroIteraciones++;
    }

    printf("Resultado:\n");
    for (int i = 0; i < N; i++)
    {
        for(int j=0; j<N; j++)
        {
            printf("%f, ",M[i*N+j]);
        }
        printf("\n");
    }
    printf("\nEl tiempo en segundos es %f y se realizaron %d iteraciones\n",dwalltime()-timetick,nroIteraciones);

}

void fProcesoDelMedio()
{
    float vecinoIzq, vecinoDer;
    float comparacion;

    int inicio= 1;
    int final= N / nrProcesos - 1;
    while(!convergio)
    {
        MPI_Scatter(M,N*N/nrProcesos,MPI_FLOAT,MRec,N*N/nrProcesos,MPI_FLOAT,0,MPI_COMM_WORLD);
        MPI_Recv(filaVecArriba,N,MPI_FLOAT,miID-1,99,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Send(&MRec[0],N,MPI_FLOAT,miID-1,99,MPI_COMM_WORLD);
        MPI_Send(&MRec[N*N/nrProcesos-N],N,MPI_FLOAT,miID+1,99,MPI_COMM_WORLD);
        MPI_Recv(filaVecAbajo,N,MPI_FLOAT,miID+1,99,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        Maux[0]= (  filaVecArriba[0] + filaVecArriba[1] +
                    MRec[0] +     MRec[1] +
                    MRec[N] +   MRec[N+1] ) * divSeis;
        for (int i = 1; i < N-1; i++)
        {
            Maux[i] =   (filaVecArriba[i-1] + filaVecArriba[i] + filaVecArriba[i+1] +
                        MRec[i-1] +        MRec[i] +      MRec[i+1] +
                        MRec[N + i-1] +    MRec[N + i] +  MRec[N + i+1]) * divNueve;
        }
        Maux[N-1]= (    filaVecArriba[N-2] + filaVecArriba[N-1] +
                        MRec[N-2] +     MRec[N-1] +
                        MRec[2*N - 2] + MRec[2*N -1]) * divSeis;

        for (int i = inicio; i < final; i++)
        {
            Maux[i*N] = ( MRec[(i-1)*N] +    MRec[(i-1)*N + 1] +
                        MRec[i*N] +        MRec[i*N + 1] +
                        MRec[(i+1)*N] +    MRec[(i+1)*N + 1] ) * divSeis;

            Maux[i*N + (N-1)] = ( MRec[i*N - 2] +    MRec[i*N-1] +
                                MRec[ i*N + (N-2)] +        MRec[i*N + (N-1)] +
                                MRec[(i+1)*N + (N-2)] +    MRec[(i+1)*N + (N-1)] ) * divSeis;

            for (int j = 1; j < N-1; j++)
            {
                Maux[i*N+j] = ( MRec[(i-1)*N+ (j-1)] +     MRec[(i-1)*N+(j)] +        MRec[(i-1)*N+ (j+1)] 
                        +       MRec[(i)*N+ (j-1)] +       MRec[(i)*N+ j] +           MRec[(i)*N+  (j+1)]
                        +       MRec[(i+1)*N+ (j-1)] +     MRec[(i+1)*N+ (j)] +       MRec[(i+1)*N+ (j+1)]) * divNueve;
            }   
        }

        //Ultima fila de las tareas en el medio SACAR EL FOR
        for (int i = final; i < final+1; i++)
        {
            Maux[i*N] = ( MRec[(i-1)*N] +    MRec[(i-1)*N + 1] +
                        MRec[i*N] +        MRec[i*N + 1] +
                        filaVecAbajo[0] +    filaVecAbajo[1] ) * divSeis;

            Maux[i*N + (N-1)] = ( MRec[i*N - 2] +    MRec[i*N-1] +
                                MRec[ i*N + (N-2)] +        MRec[i*N + (N-1)] +
                                filaVecAbajo[N-2] +    filaVecAbajo[N-1] ) * divSeis;
            for (int j = 1; j < N-1; j++)
            {
                Maux[i*N+j] = ( MRec[(i-1)*N+ (j-1)] +     MRec[(i-1)*N+(j)] +        MRec[(i-1)*N+ (j+1)] 
                        +       MRec[(i)*N+ (j-1)] +       MRec[(i)*N+ j] +           MRec[(i)*N+  (j+1)]
                        +       filaVecAbajo[j-1] +     filaVecAbajo[j] +       filaVecAbajo[j+1]) * divNueve;
            }
        }

        MPI_Bcast(&Maux[0],1,MPI_FLOAT,0,MPI_COMM_WORLD);

        convergioaux=1;
        for (int i = 0; i < N * N / nrProcesos; i++)
        {
            comparacion= Maux[0] - Maux[i];
            if ((comparacion > VALORPRECISIONP) || (comparacion < VALORPRECISIONN))
            {
                convergioaux= 0;
                break;
            }
        }

        MPI_Allreduce(&convergioaux,&convergio,1,MPI_INT,MPI_LAND,MPI_COMM_WORLD);

        MPI_Gather(Maux,N*N/nrProcesos,MPI_FLOAT,M,N*N/nrProcesos,MPI_FLOAT,0,MPI_COMM_WORLD);
    }
}

void fProcesoUltimo()
{
    float vecinoIzq, vecinoDer;
    float comparacion;

    int inicio= 1;
    int final= N / nrProcesos - 1;
    while(!convergio)
    {
        MPI_Scatter(M,N*N/nrProcesos,MPI_FLOAT,MRec,N*N/nrProcesos,MPI_FLOAT,0,MPI_COMM_WORLD);
        MPI_Recv(filaVecArriba,N,MPI_FLOAT,miID-1,99,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Send(&MRec[0],N,MPI_FLOAT,miID-1,99,MPI_COMM_WORLD);
        //MPI_Send(&VRec[N/nrProcesos-1],1,MPI_FLOAT,miID+1,99,MPI_COMM_WORLD);
        //MPI_Recv(&vecinoDer,1,MPI_FLOAT,miID+1,99,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        Maux[0]= (  filaVecArriba[0] + filaVecArriba[1] +
                    MRec[0] +     MRec[1] +
                    MRec[N] +   MRec[N+1] ) * divSeis;
        for (int i = 1; i < N-1; i++)
        {
            Maux[i] =   (filaVecArriba[i-1] + filaVecArriba[i] + filaVecArriba[i+1] +
                        MRec[i-1] +        MRec[i] +      MRec[i+1] +
                        MRec[N + i-1] +    MRec[N + i] +  MRec[N + i+1]) * divNueve;
        }
        Maux[N-1]= (    filaVecArriba[N-2] + filaVecArriba[N-1] +
                        MRec[N-2] +     MRec[N-1] +
                        MRec[2*N - 2] + MRec[2*N -1]) * divSeis;

        for (int i = inicio; i < final; i++)
        {
            Maux[i*N] = ( MRec[(i-1)*N] +    MRec[(i-1)*N + 1] +
                        MRec[i*N] +        MRec[i*N + 1] +
                        MRec[(i+1)*N] +    MRec[(i+1)*N + 1] ) * divSeis;

            Maux[i*N + (N-1)] = ( MRec[i*N - 2] +    MRec[i*N-1] +
                                MRec[ i*N + (N-2)] +        MRec[i*N + (N-1)] +
                                MRec[(i+1)*N + (N-2)] +    MRec[(i+1)*N + (N-1)] ) * divSeis;

            for (int j = 1; j < N-1; j++)
            {
                Maux[i*N+j] = ( MRec[(i-1)*N+ (j-1)] +     MRec[(i-1)*N+(j)] +        MRec[(i-1)*N+ (j+1)] 
                        +       MRec[(i)*N+ (j-1)] +       MRec[(i)*N+ j] +           MRec[(i)*N+  (j+1)]
                        +       MRec[(i+1)*N+ (j-1)] +     MRec[(i+1)*N+ (j)] +       MRec[(i+1)*N+ (j+1)]) * divNueve;
            }   
        }

        //Ultima fila de la tarea ultima SACAR EL FOR
        //Tercera esquina
        Maux[final*N]= (    MRec[(final-1)*N] + MRec[(final-1)*N+1] +
                            MRec[final*N] +     MRec[(final*N)+1]) * divCuatro;
        for (int i = 1; i < N-1; i++)
        {
            Maux[(final*N)+i] = (   MRec[(final-1)*N-1+i] +   MRec[(final-1)*N+i] +    MRec[(final-1)*N+1+i] +
                                MRec[final*N-1+i] +           MRec[(final*N)+i] +      MRec[(final*N)+1+i]) * divSeis;
        }

        //Cuarta esquina
        Maux[(final*N)+N-1]= (  MRec[(final*N)-2] +     MRec[(final*N)-1] +
                                MRec[(final*N)+N-2] +   MRec[(final*N)+N-1]) * divCuatro;

        MPI_Bcast(&Maux[0],1,MPI_FLOAT,0,MPI_COMM_WORLD);

        convergioaux=1;
        for (int i = 0; i < N * N / nrProcesos; i++)
        {
            comparacion= Maux[0] - Maux[i];
            if ((comparacion > VALORPRECISIONP) || (comparacion < VALORPRECISIONN))
            {
                convergioaux= 0;
                break;
            }
        }

        MPI_Allreduce(&convergioaux,&convergio,1,MPI_INT,MPI_LAND,MPI_COMM_WORLD);

        MPI_Gather(Maux,N*N/nrProcesos,MPI_FLOAT,M,N*N/nrProcesos,MPI_FLOAT,0,MPI_COMM_WORLD);
    }
}

int main(int argc, char* argv[]){
    N= atoi(argv[1]);
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&miID);
    MPI_Comm_size(MPI_COMM_WORLD,&nrProcesos);

    MRec=(float *)malloc(sizeof(float)*N*N/nrProcesos);
    Maux=(float *)malloc(sizeof(float)*N*N/nrProcesos);
    filaVecAbajo=(float *)malloc(sizeof(float)*N);
    filaVecArriba=(float *)malloc(sizeof(float)*N);
    convergio=0;

    
    if((miID != 0) && (miID != nrProcesos-1))
    {
        fProcesoDelMedio();
    }
    else
    {
        if (miID == 0)
        {
            fProcesoRoot();
            free(M);
        }else
        {
            fProcesoUltimo();
        }
    }

    free(MRec);
    free(Maux);
    
    MPI_Finalize();
    return 0;
}