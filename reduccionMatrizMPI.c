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

void fProcesoRoot(int N,int nrProcesos,int miID)
{
    int convergio;
    int convergioaux;
    float *MRec;
    float *M;
    float *Maux;
    float *filaVecAbajo;
    float *swap;
    float primerValor;
    float divCuatro= 1.0/4.0;
    float divSeis= 1.0/6.0;
    float divNueve= 1.0/9.0;
    float comparacion;
    double timetick;
    double tiempoEnSeg;
    int nroIteraciones= 0;

    //INICIO Y FINAL DE FILAS A PROCESAR
    int inicio= 1;
    int final= (N / nrProcesos) - 1; 
    
    //ALOCACION EN MEMORIA DE LA MATRIZ DE N x N ELEMENTOS
    M=(float *)malloc(sizeof(float)*N*N);
    //ALOCACION DE MEMORIA PARA LAS MATRICES QUE RECIBE CADA PROCESO Y LA MATRIZ AUXILIAR
    Maux=(float *)malloc(sizeof(float)*N*N/nrProcesos);
    MRec=(float *)malloc(sizeof(float)*N*N/nrProcesos);

    //ALOCACION DEL VECTOR PARA LAS FILAS RECIBIDAS DE LOS PROCESOS VECINOS
    filaVecAbajo=(float *)malloc(sizeof(float)*N);
    convergio=0;

    //INICIALIZACION DE MATRIZ
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            M[i*N+j] = randFP(0.0,1.0);
        }
    }

    //COMIENZO A CONTAR EL TIEMPO ARRANCA EL PROCESAMIENTO
    timetick= dwalltime();
    //REALIZO LA DIVISION DE CANTIDAD DE ELEMENTOS PARA CADA PROCESO
    MPI_Scatter(M,N*N/nrProcesos,MPI_FLOAT,MRec,N*N/nrProcesos,MPI_FLOAT,0,MPI_COMM_WORLD);
    while(!convergio)
    {
        //LE ENVIO A MI PROCESO SIGUIENTE MI ULTIMA FILA
        MPI_Send(&MRec[N*N/nrProcesos-N],N,MPI_FLOAT,miID+1,99,MPI_COMM_WORLD);
        //RECIBO DE MI PROCESO SIGUIENTE SU PRIMERA FILA
        MPI_Recv(filaVecAbajo,N,MPI_FLOAT,miID+1,99,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        //PROCESAMIENTO PRIMERA ESQUINA
        Maux[0]= (MRec[0] + MRec[1] + MRec[N] + MRec[N+1]) * divCuatro;

        //PROCESAMIENTO BANDA LATERAL SUPERIOR
        for (int i = 1; i < N-1; i++)
        {
            Maux[i] =   (MRec[i-1] +        MRec[i] +      MRec[i+1] +
                    MRec[N + i-1] +    MRec[N + i] +  MRec[N + i+1]) * divSeis;
        }

        //PROCESAMIENTO SEGUNDA ESQUINA
        Maux[N-1]= (MRec[N-2] + MRec[N-1] + MRec[2*N - 2] + MRec[2*N -1]) * divCuatro;

        //PROCESAMIENTO VALORES INTERMEDIOS
        for (int i = inicio; i < final; i++)
        {
            //PROCESAMIENTO VALORES BANDA LATERAL IZQUIERDA
            Maux[i*N] = ( MRec[(i-1)*N] +    MRec[(i-1)*N + 1] +
                        MRec[i*N] +        MRec[i*N + 1] +
                        MRec[(i+1)*N] +    MRec[(i+1)*N + 1] ) * divSeis;

            //PROCESAMIENTO VALORES BANDA LATERAL DERECHA
            Maux[i*N + (N-1)] = ( MRec[i*N - 2] +    MRec[i*N-1] +
                                MRec[ i*N + (N-2)] +        MRec[i*N + (N-1)] +
                                MRec[(i+1)*N + (N-2)] +    MRec[(i+1)*N + (N-1)] ) * divSeis;

            //PROCESAMIENTO VALORES INTERMEDIOS
            for (int j = 1; j < N-1; j++)
            {
                Maux[i*N+j] = ( MRec[(i-1)*N+ (j-1)] +     MRec[(i-1)*N+(j)] +        MRec[(i-1)*N+ (j+1)] 
                        +       MRec[(i)*N+ (j-1)] +       MRec[(i)*N+ j] +           MRec[(i)*N+  (j+1)]
                        +       MRec[(i+1)*N+ (j-1)] +     MRec[(i+1)*N+ (j)] +       MRec[(i+1)*N+ (j+1)]) * divNueve;
            }   
        }

        //PROCESAMIENTO ULTIMA FILA, TERCERA Y CUARTA ESQUINA Y BANDA LATERAL INFERIOR -- EL FOR SE PUEDE SACAR SOLO ITERA UNA VEZ
        Maux[final*N] = ( MRec[(final-1)*N] +    MRec[(final-1)*N + 1] +
                    MRec[final*N] +        MRec[final*N + 1] +
                    filaVecAbajo[0] +    filaVecAbajo[1] ) * divSeis;

        for (int j = 1; j < N-1; j++)
        {
            Maux[final*N+j] = ( MRec[(final-1)*N+ (j-1)] +     MRec[(final-1)*N+(j)] +        MRec[(final-1)*N+ (j+1)] 
                    +       MRec[(final)*N+ (j-1)] +       MRec[(final)*N+ j] +           MRec[(final)*N+  (j+1)]
                    +       filaVecAbajo[j-1] +     filaVecAbajo[j] +       filaVecAbajo[j+1]) * divNueve;
        }

        Maux[final*N + (N-1)] = ( MRec[final*N - 2] +    MRec[final*N-1] +
                            MRec[ final*N + (N-2)] +        MRec[final*N + (N-1)] +
                            filaVecAbajo[N-2] +    filaVecAbajo[N-1] ) * divSeis;

        //LE ENVIO A TODOS LOS PROCESOS EL RESULTADO DE LA PRIMERA POSICION
        MPI_Bcast(&Maux[0],1,MPI_FLOAT,0,MPI_COMM_WORLD);

        //CHEQUEO SI CONVERGIO MI PARTE
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

        swap= MRec;
        MRec= Maux;
        Maux= swap;

        //SE REALIZA LA OPERACION LAND EN LA VARIABLE DE CONVERGENCIA PARA CHEQUEAR SI TODOS LOS PROCESOS CONVERGIERON Y EL RESULTADO LO OBTIENEN TODOS LOS PROCESOS
        MPI_Allreduce(&convergioaux,&convergio,1,MPI_INT,MPI_LAND,MPI_COMM_WORLD);

        nroIteraciones++;
    }

    //REALIZO EL GATHER DONDE RECIBO EL RESULTADO DE TODOS LOS PROCESOS Y LO ALMACENO EN LA MATRIZ ORIGINAL
    MPI_Gather(MRec,N*N/nrProcesos,MPI_FLOAT,M,N*N/nrProcesos,MPI_FLOAT,0,MPI_COMM_WORLD);

    //PARO DE CONTAR TERMINO EL PROCESAMIENTO
    tiempoEnSeg= dwalltime() - timetick;
    //DESCOMENTAR SI SE QUIERE IMPRIMIR EL RESULTADO
    /*printf("Resultado:\n");
    for (int i = 0; i < N; i++)
    {
        for(int j=0; j<N; j++)
        {
            printf("%f, ",M[i*N+j]);
        }
        printf("\n");
    }*/
    printf("\nEl tiempo en segundos es %f y se realizaron %d iteraciones\n",tiempoEnSeg,nroIteraciones);

    free(M);
    free(Maux);
    free(MRec);

}

void fProcesoDelMedio(int N,int nrProcesos,int miID)
{
    int convergio;
    int convergioaux;
    float *MRec;
    float *M;
    float *Maux;
    float *filaVecAbajo, *filaVecArriba;
    float *swap;
    float primerValor;
    float divCuatro= 1.0/4.0;
    float divSeis= 1.0/6.0;
    float divNueve= 1.0/9.0;
    float comparacion;

    //ALOCACION DE MEMORIA PARA LAS MATRICES QUE RECIBE CADA PROCESO Y LA MATRIZ AUXILIAR
    MRec=(float *)malloc(sizeof(float)*N*N/nrProcesos);
    Maux=(float *)malloc(sizeof(float)*N*N/nrProcesos);

    //ALOCACION DEL VECTOR PARA LAS FILAS RECIBIDAS DE LOS PROCESOS VECINOS
    filaVecAbajo=(float *)malloc(sizeof(float)*N);
    filaVecArriba=(float *)malloc(sizeof(float)*N);
    convergio=0;

    int inicio= 1;
    int final= N / nrProcesos - 1;
    //RECIBO LA CANTIDAD DE VALORES QUE TENGO QUE PROCESAR
    MPI_Scatter(M,N*N/nrProcesos,MPI_FLOAT,MRec,N*N/nrProcesos,MPI_FLOAT,0,MPI_COMM_WORLD);
    while(!convergio)
    {
        
        //RECIBO DE MI PROCESO ANTERIOR SU ULTIMA FILA
        MPI_Recv(filaVecArriba,N,MPI_FLOAT,miID-1,99,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        //LE ENVIO A MI PROCESO ANTERIOR MI PRIMERA FILA
        MPI_Send(&MRec[0],N,MPI_FLOAT,miID-1,99,MPI_COMM_WORLD);
        //LE ENVIO A MI PROCESO SIGUIENTE MI ULTIMA FILA
        MPI_Send(&MRec[N*N/nrProcesos-N],N,MPI_FLOAT,miID+1,99,MPI_COMM_WORLD);
        //RECIBO DE MI PROCESO SIGUIENTE SU PRIMERA FILA
        MPI_Recv(filaVecAbajo,N,MPI_FLOAT,miID+1,99,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        //PROCESO MI PRIMER FILA
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

        //PROCESAMIENTO DE VALORES INTERMEDIOS
        for (int i = inicio; i < final; i++)
        {
            //PROCESAMIENTO VALORES BANDA LATERAL IZQUIERDA
            Maux[i*N] = ( MRec[(i-1)*N] +    MRec[(i-1)*N + 1] +
                        MRec[i*N] +        MRec[i*N + 1] +
                        MRec[(i+1)*N] +    MRec[(i+1)*N + 1] ) * divSeis;

            //PROCESAMIENTO VALORES INTERMEDIOS
            for (int j = 1; j < N-1; j++)
            {
                Maux[i*N+j] = ( MRec[(i-1)*N+ (j-1)] +     MRec[(i-1)*N+(j)] +        MRec[(i-1)*N+ (j+1)] 
                        +       MRec[(i)*N+ (j-1)] +       MRec[(i)*N+ j] +           MRec[(i)*N+  (j+1)]
                        +       MRec[(i+1)*N+ (j-1)] +     MRec[(i+1)*N+ (j)] +       MRec[(i+1)*N+ (j+1)]) * divNueve;
            } 

            //PROCESAMIENTO VALORES BANDA LATERAL DERECHA
            Maux[i*N + (N-1)] = ( MRec[i*N - 2] +    MRec[i*N-1] +
                                MRec[ i*N + (N-2)] +        MRec[i*N + (N-1)] +
                                MRec[(i+1)*N + (N-2)] +    MRec[(i+1)*N + (N-1)] ) * divSeis;  
        }

        //ULTIMA FILA DEL PROCESO -- EL FOR NO ES NECESARIO SE PUEDE SACAR
        Maux[final*N] = ( MRec[(final-1)*N] +    MRec[(final-1)*N + 1] +
                    MRec[final*N] +        MRec[final*N + 1] +
                    filaVecAbajo[0] +    filaVecAbajo[1] ) * divSeis;

        for (int j = 1; j < N-1; j++)
        {
            Maux[final*N+j] = ( MRec[(final-1)*N+ (j-1)] +     MRec[(final-1)*N+(j)] +        MRec[(final-1)*N+ (j+1)] 
                    +       MRec[(final)*N+ (j-1)] +       MRec[(final)*N+ j] +           MRec[(final)*N+  (j+1)]
                    +       filaVecAbajo[j-1] +     filaVecAbajo[j] +       filaVecAbajo[j+1]) * divNueve;
        }

        Maux[final*N + (N-1)] = ( MRec[final*N - 2] +    MRec[final*N-1] +
                            MRec[ final*N + (N-2)] +        MRec[final*N + (N-1)] +
                            filaVecAbajo[N-2] +    filaVecAbajo[N-1] ) * divSeis;

        //RECIBO EL RESULTADO DEL PRIMER VALOR PARA CHEQUEAR CONVERGENCIA
        MPI_Bcast(&primerValor,1,MPI_FLOAT,0,MPI_COMM_WORLD);

        //CHEQUEO SI MI PARTE CONVERGIO
        convergioaux=1;
        for (int i = 0; i < N * N / nrProcesos; i++)
        {
            comparacion= primerValor - Maux[i];
            if ((comparacion > VALORPRECISIONP) || (comparacion < VALORPRECISIONN))
            {
                convergioaux= 0;
                break;
            }
        }

        swap= MRec;
        MRec= Maux;
        Maux= swap;

        //SE CHEQUEA LA CONVERGENCIA TOTAL DE LA MATRIZ
        MPI_Allreduce(&convergioaux,&convergio,1,MPI_INT,MPI_LAND,MPI_COMM_WORLD);
    }
    //SE LE ENVIA EL RESULTADO AL PROCESO 0 PARA QUE ESTE TENGA EL RESULTADO TOTAL DE LA MATRIZ
    MPI_Gather(MRec,N*N/nrProcesos,MPI_FLOAT,M,N*N/nrProcesos,MPI_FLOAT,0,MPI_COMM_WORLD);

    free(MRec);
    free(Maux);
}

void fProcesoUltimo(int N,int nrProcesos,int miID)
{
    int convergio;
    int convergioaux;
    float *MRec;
    float *M;
    float *Maux;
    float *filaVecArriba;
    float *swap;
    float primerValor;
    float divCuatro= 1.0/4.0;
    float divSeis= 1.0/6.0;
    float divNueve= 1.0/9.0;
    float comparacion;


    //ALOCACION DE MEMORIA PARA LAS MATRICES QUE RECIBE CADA PROCESO Y LA MATRIZ AUXILIAR
    MRec=(float *)malloc(sizeof(float)*N*N/nrProcesos);
    Maux=(float *)malloc(sizeof(float)*N*N/nrProcesos);

    //ALOCACION DEL VECTOR PARA LAS FILAS RECIBIDAS DE LOS PROCESOS VECINOS
    filaVecArriba=(float *)malloc(sizeof(float)*N);
    convergio=0;

    int inicio= 1;
    int final= N / nrProcesos - 1;
    //RECIBO LOS ELEMENTOS QUE TENGO QUE PROCESAR
    MPI_Scatter(M,N*N/nrProcesos,MPI_FLOAT,MRec,N*N/nrProcesos,MPI_FLOAT,0,MPI_COMM_WORLD);
    while(!convergio)
    {
        //RECIBO DE MI PROCESO ANTERIOR SU ULTIMA FILA
        MPI_Recv(filaVecArriba,N,MPI_FLOAT,miID-1,99,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        //LE ENVIO A MI PROCESO ANTERIOR MI PRMERA FILA
        MPI_Send(&MRec[0],N,MPI_FLOAT,miID-1,99,MPI_COMM_WORLD);

        //PROCESAMIENTO DE MI PRIMERA FILA
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

        //PROCESAMIENTO VALORES INTERMEDIOS
        for (int i = inicio; i < final; i++)
        {
            //PROCESAMIENTO BANDA LATERAL IZQUIERDA
            Maux[i*N] = ( MRec[(i-1)*N] +    MRec[(i-1)*N + 1] +
                        MRec[i*N] +        MRec[i*N + 1] +
                        MRec[(i+1)*N] +    MRec[(i+1)*N + 1] ) * divSeis;

            //PROCESAMIENTO BANDA LATERAL DERECHA
            Maux[i*N + (N-1)] = ( MRec[i*N - 2] +    MRec[i*N-1] +
                                MRec[ i*N + (N-2)] +        MRec[i*N + (N-1)] +
                                MRec[(i+1)*N + (N-2)] +    MRec[(i+1)*N + (N-1)] ) * divSeis;

            //PROCESAMIENTO VALORES INTERMEDIOS
            for (int j = 1; j < N-1; j++)
            {
                Maux[i*N+j] = ( MRec[(i-1)*N+ (j-1)] +     MRec[(i-1)*N+(j)] +        MRec[(i-1)*N+ (j+1)] 
                        +       MRec[(i)*N+ (j-1)] +       MRec[(i)*N+ j] +           MRec[(i)*N+  (j+1)]
                        +       MRec[(i+1)*N+ (j-1)] +     MRec[(i+1)*N+ (j)] +       MRec[(i+1)*N+ (j+1)]) * divNueve;
            }   
        }

        //ULTIMA FILA 
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

        //RECIBO EL RESULTADO DE LA PRIMERA POSICION
        MPI_Bcast(&primerValor,1,MPI_FLOAT,0,MPI_COMM_WORLD);

        //CHEQUEO SI CONVERGE MI PARTE
        convergioaux=1;
        for (int i = 0; i < N * N / nrProcesos; i++)
        {
            comparacion= primerValor - Maux[i];
            if ((comparacion > VALORPRECISIONP) || (comparacion < VALORPRECISIONN))
            {
                convergioaux= 0;
                break;
            }
        }

        swap= MRec;
        MRec= Maux;
        Maux= swap;

        //SE CHEQUEA SI LA MATRIZ CONVERGIO TOTALMENTE
        MPI_Allreduce(&convergioaux,&convergio,1,MPI_INT,MPI_LAND,MPI_COMM_WORLD);
    }
    //SE LE ENVIA AL PROCESO 0 EL RESULTADO DE LO PROCESADO PARA QUE TENGA LA MATRIZ TOTAL
    MPI_Gather(MRec,N*N/nrProcesos,MPI_FLOAT,M,N*N/nrProcesos,MPI_FLOAT,0,MPI_COMM_WORLD);

    free(MRec);
    free(Maux);
}


/***************************************
            FUNCION MAIN
 ***************************************/
int main(int argc, char* argv[]){
    int N,nrProcesos,miID;
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
        //SI SOY EL PROCESO 0
        if (miID == 0)
        {
            fProcesoRoot(N,nrProcesos,miID);
        }else
        {
            fProcesoUltimo(N,nrProcesos,miID);
        }
    }
    
    MPI_Finalize();
    return 0;
}