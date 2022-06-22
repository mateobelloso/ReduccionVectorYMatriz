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

//FUNCION PROCESO ROOT
void fProcesoRoot(int N,int nrProcesos,int miID)
{
    int convergio;
    int convergioaux;
    float *VRec;
    float *V;
    float *Vaux;
    float *swap;
    float primerValor;
    float divDos= 1.0/2.0;
    float divTres= 1.0/3.0;
    float vecinoDer;
    float comparacion;
    double timetick;
    double tiempoEnSeg;
    int nroIteraciones= 0;
    
    //ALOCACION DE MEMORIA PARA EL VECTOR DE TAMAÑO N
    V=(float *)malloc(sizeof(float)*N);

    //ALOCACION DE MEMORIA DE LOS VECTORES QUE RECIBE CADA PROCESO Y SU VECTOR AUXILIAR
    VRec=(float *)malloc(sizeof(float)*N/nrProcesos);
    Vaux=(float *)malloc(sizeof(float)*N/nrProcesos);

    //VARIABLE PARA INDICAR SI LA MATRIZ CONVERGIO
    convergio=0;

    //INICIALIZACION DEL VECTOR
    for (int i = 0; i < N; i++)
    {
        V[i]= randFP(0.0,1.0);
    }

    //COMIENZO A CONTAR TIEMPO ARRANCA EL PROCESAMIENTO
    timetick= dwalltime();
    //DIVIDO PARA CADA PROCESO LA CANTIDAD DE VALORES DEL VECTOR QUE LE TOCAN LO VAN A RECIBIR EN LA VARIABLE VRec
    MPI_Scatter(V,N/nrProcesos,MPI_FLOAT,VRec,N/nrProcesos,MPI_FLOAT,0,MPI_COMM_WORLD);
    while(!convergio)
    {
        //LE ENVIO A MI PROCESO SIGUIENTE MI ULTIMO VALOR
        MPI_Send(&VRec[N/nrProcesos-1],1,MPI_FLOAT,miID+1,99,MPI_COMM_WORLD);
        //RECIBO DE MI PROCESO SIGUIENTE SU PRIME VALOR
        MPI_Recv(&vecinoDer,1,MPI_FLOAT,miID+1,99,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        //PROCESAMIENTO PRIMERA POSICION
        Vaux[0]= (VRec[0] + VRec[1]) * divDos;
        
        //PROCESAMIENTO VALORES INTERMEDIOS
        for (int i = 1; i < ((N/nrProcesos)-1); i++)
        {
            Vaux[i]= (VRec[i-1] + VRec[i] + VRec[i+1]) * divTres;
        }

        //PROCESAMIENTO DE MI ULTIMO VALOR
        Vaux[N/nrProcesos-1]= (VRec[N/nrProcesos-2] + VRec[N/nrProcesos-1] + vecinoDer) * divTres;

        //LE ENVIO A TODOS LOS PROCESOS EL RESULTADO DE MI PRIMERA POSICION
        MPI_Bcast(&Vaux[0],1,MPI_FLOAT,0,MPI_COMM_WORLD);

        //CHEQUEO SI CONVIRGIO MI PARTE
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

        swap= VRec;
        VRec= Vaux;
        Vaux= swap;

        //A TRAVES DEL ALLREDUCE SE CHEQUEA LA CONVERGENCIA DE TODOS LOS PROCESOS Y EL RESULTADO LO OBTIENEN TODOS
        MPI_Allreduce(&convergioaux,&convergio,1,MPI_INT,MPI_LAND,MPI_COMM_WORLD);
        nroIteraciones++;
    }

    //REALIZO EL GATTER PARA OBTENER EL RESULTADO DE CADA PROCESO Y LO RECIBO EN LA VARIABLE DEL VECTOR
    MPI_Gather(VRec,N/nrProcesos,MPI_FLOAT,V,N/nrProcesos,MPI_FLOAT,0,MPI_COMM_WORLD);

    //TERMINO EL PROCESAMIENTO PARO DE CONTAR
    tiempoEnSeg= dwalltime() - timetick;
    //DESCOMENTAR SI SE QUIERE IMPRIMIR EL VECTOR RESULTADO
    /*printf("Resultado:\n");
    for (int i = 0; i < N; i++)
    {
        printf("%f, ",V[i]);
    }*/
    printf("\nEl tiempo en segundos es %f y se realizaron %d iteraciones\n",tiempoEnSeg,nroIteraciones);
    free(VRec);
    free(Vaux);
    free(V);
}

//FUNCION PROCESOS DEL MEDIO
void fProcesoDelMedio(int N,int nrProcesos,int miID)
{
    int convergio;
    int convergioaux;
    float *VRec;
    float *V;
    float *Vaux;
    float *swap;
    float primerValor;
    float divDos= 1.0/2.0;
    float divTres= 1.0/3.0;
    float vecinoIzq, vecinoDer;
    float comparacion;

    //ALOCACION DE MEMORIA DE LOS VECTORES QUE RECIBE CADA PROCESO Y SU VECTOR AUXILIAR
    VRec=(float *)malloc(sizeof(float)*N/nrProcesos);
    Vaux=(float *)malloc(sizeof(float)*N/nrProcesos);

    //VARIABLE PARA INDICAR SI LA MATRIZ CONVERGIO
    convergio=0;

    //RECIBO EN M VARIABLE VRec LOS VALORES QUE TENGO QUE PROCESAR
    MPI_Scatter(V,N/nrProcesos,MPI_FLOAT,VRec,N/nrProcesos,MPI_FLOAT,0,MPI_COMM_WORLD);
    while(!convergio)
    {
        //LE ENVIO A MI PROCESO ANTERIOR MI PRIMERA POSICION
        MPI_Send(&VRec[0],1,MPI_FLOAT,miID-1,99,MPI_COMM_WORLD);
        //LE ENVIO A MI PROCESO SIGUIENTE MI ULTIMA POSICION
        MPI_Send(&VRec[N/nrProcesos-1],1,MPI_FLOAT,miID+1,99,MPI_COMM_WORLD);
        //RECIBO DE MI PROCESO ANTERIOR SU ULTIMA POSICION
        MPI_Recv(&vecinoIzq,1,MPI_FLOAT,miID-1,99,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        //RECIBO DE MI PROCESO SIGUIENTE SU PRIMERA POSICION
        MPI_Recv(&vecinoDer,1,MPI_FLOAT,miID+1,99,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        //PROCESO MI PRIMER VALOR
        Vaux[0]= (vecinoIzq + VRec[0] + VRec[1]) * divTres;

        //PROCESAMIENTO DE LOS VALORES INTERMEDIOS
        for (int i = 1; i < ((N/nrProcesos)-1); i++)
        {
            Vaux[i]= (VRec[i-1] + VRec[i] + VRec[i+1]) * divTres;
        }

        //PROCESAMIENTO DE MI ULTIMO VALOR
        Vaux[N/nrProcesos-1]= (VRec[N/nrProcesos-2] + VRec[N/nrProcesos-1] + vecinoDer) * divTres;

        //RECIBO EL RESULTADO DEL PRIMER VALOR DE LA MATRIZ
        MPI_Bcast(&primerValor,1,MPI_FLOAT,0,MPI_COMM_WORLD);
        
        //CHEQUEO SI CONVERGIO MI PARTE
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

        swap= VRec;
        VRec= Vaux;
        Vaux= swap;

        //SE CHEQUEA LA CONVERGENCIA TOTAL
        MPI_Allreduce(&convergioaux,&convergio,1,MPI_INT,MPI_LAND,MPI_COMM_WORLD);     
    }

    //SE ENVIA LA PARTE PROCESADA AL PROCESO 0 PARA QUE OBTENGA TODO EL VECTOR PROCESADO
    MPI_Gather(VRec,N/nrProcesos,MPI_FLOAT,V,N/nrProcesos,MPI_FLOAT,0,MPI_COMM_WORLD);
    free(VRec);
    free(Vaux); 
}

void fProcesoUltimo(int N,int nrProcesos,int miID)
{
    int convergio;
    int convergioaux;
    float *VRec;
    float *V;
    float *Vaux;
    float *swap;
    float primerValor;
    float divDos= 1.0/2.0;
    float divTres= 1.0/3.0;
    float vecinoIzq, vecinoDer;
    float comparacion;

    //ALOCACION DE MEMORIA DE LOS VECTORES QUE RECIBE CADA PROCESO Y SU VECTOR AUXILIAR
    VRec=(float *)malloc(sizeof(float)*N/nrProcesos);
    Vaux=(float *)malloc(sizeof(float)*N/nrProcesos);

    //VARIABLE PARA INDICAR SI LA MATRIZ CONVERGIO
    convergio=0;

    //RECIBO EN MI VRec LOS VALORES QUE ME TOCAN PROCESAR
    MPI_Scatter(V,N/nrProcesos,MPI_FLOAT,VRec,N/nrProcesos,MPI_FLOAT,0,MPI_COMM_WORLD);
    while(!convergio)
    {
        //LE ENVIO A MI PROCESO ANTERIOR MI PRIMERA POSICION
        MPI_Send(&VRec[0],1,MPI_FLOAT,miID-1,99,MPI_COMM_WORLD);

        //RECIBO DE MI PROCESO ANTERIOR SU ULTIMA POSICION
        MPI_Recv(&vecinoIzq,1,MPI_FLOAT,miID-1,99,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        //PROCESAMIENTO DE MI PRIMER VALOR
        Vaux[0]= (vecinoIzq + VRec[0] + VRec[1]) * divTres;

        //PROCESAMIENTO DE VALORES INTERMEDIOS
        for (int i = 1; i < ((N/nrProcesos)-1); i++)
        {
            Vaux[i]= (VRec[i-1] + VRec[i] + VRec[i+1]) * divTres;
        }

        //PROCESAMIENTO DE EL ULTIMO VALOR DEL VECTOR
        Vaux[N/nrProcesos-1]= (VRec[N/nrProcesos-2] + VRec[N/nrProcesos-1]) * divDos;

        //RECIBO EL RESULTADO DE LA PRIMERA POSICION
        MPI_Bcast(&primerValor,1,MPI_FLOAT,0,MPI_COMM_WORLD);
        
        //CHEQUEO SI MI PARTE CONVERGIO
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

        swap= VRec;
        VRec= Vaux;
        Vaux= swap;

        //SE OBTIENE EL RESULTADO DE LA CONVERGENCIA TOTAL
        MPI_Allreduce(&convergioaux,&convergio,1,MPI_INT,MPI_LAND,MPI_COMM_WORLD);        
    }

    //SE LE MANDA AL PROCESO 0 EL RESULTADO DEL VECTOR PARA QUE ESTE OBTENGA EL RESULTADO TOTAL
    MPI_Gather(VRec,N/nrProcesos,MPI_FLOAT,V,N/nrProcesos,MPI_FLOAT,0,MPI_COMM_WORLD);
    free(VRec);
    free(Vaux);
}


/***************************************
            FUNCION MAIN
 ***************************************/
int main(int argc, char* argv[]){
    int miID;
    int nrProcesos;
    int N;
    N= atoi(argv[1]);
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&miID);
    MPI_Comm_size(MPI_COMM_WORLD,&nrProcesos);

    //SI NO SOY EL PROCESO 0 NI EL ULTIMO SOY UN PROCESO DEL MEDIO
    if((miID != 0) && (miID != nrProcesos-1))
    {
        fProcesoDelMedio(N,nrProcesos,miID);
    }
    else
    {
        //SOY EL PROCESO 0
        if (miID == 0)
        {
            fProcesoRoot(N,nrProcesos,miID);
        }else   //SOY EL ULTIMO PROCESO
        {
            fProcesoUltimo(N,nrProcesos,miID);
        }
    }
    
    MPI_Finalize();
    return 0;
}