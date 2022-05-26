#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<mpi.h>

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
char message[20];

void fProcesoTipoA()
{
    strcpy(message,"Hola MPI");
    for (int i = 1; i < nrProcesos; ++i)
    {
        MPI_Send(message,strlen(message)+1,MPI_CHAR,i,99,MPI_COMM_WORLD);
    }
    printf("%s soy el proceso %d\n",message,miID);
}
void fProcesoTipoB()
{
    MPI_Recv(message,20,MPI_CHAR,0,99,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    printf("%s soy el proceso %d\n",message,miID);
}

int main(int argc, char* argv[]){
    N= atoi(argv[1]);
    float *VRec;
    float *V;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&miID);
    MPI_Comm_size(MPI_COMM_WORLD,&nrProcesos);

    VRec=(float *)malloc(sizeof(float)*N/nrProcesos);
    

    if (miID==0)
    {
        V=(float *)malloc(sizeof(float)*N);        
        for (int i = 0; i < N; i++)
        {
            V[i]= randFP(0.0,1.0);
        }
        for (int i = 0; i < N; i++)
        {
            printf("%f\n",V[i]);
        }
    }

    MPI_Scatter(V,N/nrProcesos,MPI_FLOAT,VRec,N/nrProcesos,MPI_FLOAT,0,MPI_COMM_WORLD);

    printf("\nSoy el proceso %d y recibi este vector:\n",miID);
    for (int i = 0; i < N/nrProcesos; i++)
    {
        printf("%f\n",VRec[i]);
    }

    /*
    if(miID== 0)
    {
        fProcesoTipoA();
    }
    else
    {
        fProcesoTipoB();
    }
    */
    MPI_Finalize();
    return 0;
}