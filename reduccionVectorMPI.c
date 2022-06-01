#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include<mpi.h>
#define VALORPRECISIONP 0.01    //VALOR DE PRECISION POSITIVO
#define VALORPRECISIONN -0.01   //VALOR DE PRECISION NEGATIVO

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


int miID;
int nrProcesos;
int N;
int convergio;
int convergioaux;
float *VRec;
float *V;
float *Vaux;
float primerValor;
float divDos= 1.0/2.0;
float divTres= 1.0/3.0;

void fProcesoRoot()
{
    float vecinoDer;
    float comparacion;
    double timetick;
    double tiempoEnSeg;
    int nroIteraciones= 0;
    
    V=(float *)malloc(sizeof(float)*N);
    for (int i = 0; i < N; i++)
    {
        V[i]= randFP(0.0,1.0);
    }

    timetick= dwalltime();
    while(!convergio)
    {
        MPI_Scatter(V,N/nrProcesos,MPI_FLOAT,VRec,N/nrProcesos,MPI_FLOAT,0,MPI_COMM_WORLD);
        MPI_Send(&VRec[N/nrProcesos-1],1,MPI_FLOAT,miID+1,99,MPI_COMM_WORLD);
        MPI_Recv(&vecinoDer,1,MPI_FLOAT,miID+1,99,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        Vaux[0]= (VRec[0] + VRec[1]) * divDos;
        
        for (int i = 1; i < ((N/nrProcesos)-1); i++)
        {
            Vaux[i]= (VRec[i-1] + VRec[i] + VRec[i+1]) * divTres;
        }
        Vaux[N/nrProcesos-1]= (VRec[N/nrProcesos-2] + VRec[N/nrProcesos-1] + vecinoDer) * divTres;

        MPI_Bcast(&Vaux[0],1,MPI_FLOAT,0,MPI_COMM_WORLD);

        
        convergioaux=1;
        for(int i = 1; i < N/nrProcesos; i++)
        {
            comparacion= Vaux[0] - Vaux[i];
            if ((comparacion > VALORPRECISIONP) || (comparacion < VALORPRECISIONN))
            {
                convergioaux= 0;
                break;              
            }
        }

        MPI_Allreduce(&convergioaux,&convergio,1,MPI_INT,MPI_LAND,MPI_COMM_WORLD);

        MPI_Gather(Vaux,N/nrProcesos,MPI_FLOAT,V,N/nrProcesos,MPI_FLOAT,0,MPI_COMM_WORLD);
        nroIteraciones++;
    }

    tiempoEnSeg= dwalltime() - timetick;
    //DESCOMENTAR SI SE QUIERE IMPRIMIR EL VECTOR RESULTADO
    /*printf("Resultado:\n");
    for (int i = 0; i < N; i++)
    {
        printf("%f, ",V[i]);
    }*/
    printf("\nEl tiempo en segundos es %f y se realizaron %d iteraciones\n",tiempoEnSeg,nroIteraciones);
}

void fProcesoDelMedio()
{
    float vecinoIzq, vecinoDer;
    float comparacion;
    while(!convergio)
    {
        MPI_Scatter(V,N/nrProcesos,MPI_FLOAT,VRec,N/nrProcesos,MPI_FLOAT,0,MPI_COMM_WORLD);
        MPI_Send(&VRec[0],1,MPI_FLOAT,miID-1,99,MPI_COMM_WORLD);
        MPI_Send(&VRec[N/nrProcesos-1],1,MPI_FLOAT,miID+1,99,MPI_COMM_WORLD);
        MPI_Recv(&vecinoIzq,1,MPI_FLOAT,miID-1,99,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Recv(&vecinoDer,1,MPI_FLOAT,miID+1,99,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        
        Vaux[0]= (vecinoIzq + VRec[0] + VRec[1]) * divTres;
        for (int i = 1; i < ((N/nrProcesos)-1); i++)
        {
            Vaux[i]= (VRec[i-1] + VRec[i] + VRec[i+1]) * divTres;
        }
        Vaux[N/nrProcesos-1]= (VRec[N/nrProcesos-2] + VRec[N/nrProcesos-1] + vecinoDer) * divTres;

        MPI_Bcast(&primerValor,1,MPI_FLOAT,0,MPI_COMM_WORLD);
        
        convergioaux=1;
        for(int i = 1; i < N/nrProcesos; i++)
        {
            comparacion= primerValor - Vaux[i];
            if ((comparacion > VALORPRECISIONP) || (comparacion < VALORPRECISIONN))
            {
                convergioaux= 0;
                break;              
            }
        }

        MPI_Allreduce(&convergioaux,&convergio,1,MPI_INT,MPI_LAND,MPI_COMM_WORLD);

        MPI_Gather(Vaux,N/nrProcesos,MPI_FLOAT,V,N/nrProcesos,MPI_FLOAT,0,MPI_COMM_WORLD);

        
    }
}

void fProcesoUltimo()
{
    float vecinoIzq, vecinoDer;
    float comparacion;
    while(!convergio)
    {
        MPI_Scatter(V,N/nrProcesos,MPI_FLOAT,VRec,N/nrProcesos,MPI_FLOAT,0,MPI_COMM_WORLD);
        MPI_Send(&VRec[0],1,MPI_FLOAT,miID-1,99,MPI_COMM_WORLD);
        MPI_Recv(&vecinoIzq,1,MPI_FLOAT,miID-1,99,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        
        Vaux[0]= (vecinoIzq + VRec[0] + VRec[1]) * divTres;
        for (int i = 1; i < ((N/nrProcesos)-1); i++)
        {
            Vaux[i]= (VRec[i-1] + VRec[i] + VRec[i+1]) * divTres;
        }
        Vaux[N/nrProcesos-1]= (VRec[N/nrProcesos-2] + VRec[N/nrProcesos-1]) * divDos;

        MPI_Bcast(&primerValor,1,MPI_FLOAT,0,MPI_COMM_WORLD);
        
        convergioaux=1;
        for(int i = 1; i < N/nrProcesos; i++)
        {
            comparacion= primerValor - Vaux[i];
            if ((comparacion > VALORPRECISIONP) || (comparacion < VALORPRECISIONN))
            {
                convergioaux= 0;
                break;              
            }
        }

        MPI_Allreduce(&convergioaux,&convergio,1,MPI_INT,MPI_LAND,MPI_COMM_WORLD);

        MPI_Gather(Vaux,N/nrProcesos,MPI_FLOAT,V,N/nrProcesos,MPI_FLOAT,0,MPI_COMM_WORLD);

        
    }
}


/***************************************
            FUNCION MAIN
 ***************************************/
int main(int argc, char* argv[]){
    N= atoi(argv[1]);
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&miID);
    MPI_Comm_size(MPI_COMM_WORLD,&nrProcesos);

    VRec=(float *)malloc(sizeof(float)*N/nrProcesos);
    Vaux=(float *)malloc(sizeof(float)*N/nrProcesos);
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
        }else
        {
            fProcesoUltimo();
        }
    }
    
    MPI_Finalize();
    return 0;
}