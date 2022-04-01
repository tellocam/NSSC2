#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <chrono>
//#include <omp.h>
#include "mpi.h"

void DirichletBC(std::vector<double> &u, const int N) {
    //u(x,1)=sin(2*pi*x)*sinh(2*pi)
    //u(x,0)=u(0,y)=u(1,y)=0 already to initialization of u
    double h=1/double(N-1);
    const double sh=sinh(2*M_PI);
    int offs=N*(N-1);
    for(int i=0;i<N;i++) {        
        u[offs+i]=sin(2*M_PI*h*i)*sh;
    }
}

double fun(const double x, const double y) {    
    return sin(2*M_PI*x)*sinh(2*M_PI*y);
}
std::vector<double> forcing(const int N) {
    std::vector<double> b(N*N,0);
    double x, y, h=1/double(N-1);
    for(int i=0;i<N*N;i++) {
        x=(i%N)*h;
        y=floor(i/N)*h;
        b[i]=(4*M_PI*M_PI)*fun(x,y);
    }
    return b;
}
std::vector<double> exactSolution(const int N) {
    std::vector<double> u(N*N,0);
    double x, y, h=1/double(N-1);
    for(int i=0;i<N*N;i++) {
        x=(i%N)*h;
        y=floor(i/N)*h;
        u[i]=fun(x,y);
    }
    return u;
}

bool isActive(const int i, const int N) {
    return (i%N!=0 && (i+1)%N!=0 && int(floor(i/N))>0 && int(floor(i/N))<(N-1));
}

std::vector<double> JacobiMethod(std::vector<double> u, std::vector<double> &b, const int N, const int iter, double &time, int myid, int numprocs) {
    double h=1/double(N-1);
    
    double aiikonst=1/(4*M_PI*M_PI*(h*h)+4);
    int i, j, k, l, sendLength, recvPos=0;
    std::vector<double> u2(N*N);
    
    std::vector<int> borderRows;
    std::vector<int> recvCounts, recvOrder;
    MPI_Status status;
    MPI_Request request, request2;
    
    
    for (i=0;i<=numprocs;i++) {
        borderRows.push_back(i*(N-2)/numprocs);
    }
    if (myid==0) {
        for (i=0;i<numprocs;i++) {
            sendLength=N*(borderRows[i+1]-borderRows[i]);
            recvCounts.push_back(sendLength);
            
        }
        recvOrder.push_back(0);
        for (i=1;i<numprocs;i++) {
            recvPos=recvPos+recvCounts[i-1];
            recvOrder.push_back(recvPos);            
        }
    }
       
    
    u2=u;
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    for(j=0;j<iter;j++) {

        for(l=borderRows[myid]+1;l<=borderRows[myid+1];l++) {
            for(k=1;k<N-1;k++){
                i=l*N+k;
                u2[i]=(b[i]*(h*h)+u[i-N]+ u[i-1]+ u[i+1]+ u[i+N])*aiikonst;                    
            }        
        }
        if(numprocs>1) {
            if(myid==0) {
                
                MPI_Isend(&u2[N*borderRows[myid+1]], N, MPI_DOUBLE, myid+1, 0, MPI_COMM_WORLD, &request);
                
                MPI_Recv(&u2[N*(borderRows[myid+1]+1)], N, MPI_DOUBLE, myid+1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                
                MPI_Wait(&request, &status);
            } else if(myid==numprocs-1) {
                
                MPI_Isend(&u2[N*(borderRows[myid]+1)], N, MPI_DOUBLE, myid-1, 0, MPI_COMM_WORLD, &request);
                
                MPI_Recv(&u2[N*(borderRows[myid])], N, MPI_DOUBLE, myid-1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
          
                MPI_Wait(&request, &status);    
            } else {
                MPI_Isend(&u2[N*borderRows[myid+1]], N, MPI_DOUBLE, myid+1, 0, MPI_COMM_WORLD, &request2);
                MPI_Isend(&u2[N*(borderRows[myid]+1)], N, MPI_DOUBLE, myid-1, 0, MPI_COMM_WORLD, &request);
                
                MPI_Recv(&u2[N*(borderRows[myid+1]+1)], N, MPI_DOUBLE, myid+1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                MPI_Recv(&u2[N*(borderRows[myid])], N, MPI_DOUBLE, myid-1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                
                MPI_Wait(&request2, &status);
                MPI_Wait(&request, &status);
            }
        }
        

        u=u2;
    }

    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    std::chrono::duration<double> timespan= std::chrono::duration_cast<::std::chrono::duration<double>>(t2-t1);
    time=timespan.count()/double(iter);
    
    if(numprocs>1) {
        if (myid == 0) {   
            MPI_Gatherv(&u[N*(borderRows[myid]+1)], N*(borderRows[myid+1]-borderRows[myid]), MPI_DOUBLE, &u[N], &recvCounts[0], &recvOrder[0]  , MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
            
        } else {
            MPI_Gatherv(&u[N*(borderRows[myid]+1)], N*(borderRows[myid+1]-borderRows[myid]), MPI_DOUBLE, NULL , NULL , NULL, NULL, 0, MPI_COMM_WORLD);
        
        }
    }    
    
    return u;
}

std::vector<double> Errors(const std::vector<double> &u, const std::vector<double> &uP, const std::vector<double> &b, const int N) {
    //u being the numeric Solution, uP being the exact Solution, b the right hand side of the PDE
    double h=1/double(N-1);
    double Au; //Au = left side of PDE
    double residual, absolut; 
    double resSQRsum=0, resMaxError=0, absSQRsum=0, absMaxError=0;
    std::vector<double> Errors(4,0);
    
    for(int i=N+1;i<(N*(N-1)-1);i++) {
        if(isActive(i, N)) {
            Au=(4*M_PI*M_PI+4/(h*h))*u[i]-(u[i-N]+ u[i-1]+ u[i+1]+ u[i+N])/(h*h);
            residual=pow(Au-b[i],2);            
            absolut=pow(u[i]-uP[i],2);
            
            resSQRsum+=residual;
            absSQRsum+=absolut;
            if(residual>resMaxError) {
                resMaxError=residual;
            }
            if(absolut>absMaxError) {
                absMaxError=absolut;
            }
        }
    }
    Errors[0]=sqrt(resSQRsum);
    Errors[1]=sqrt(resMaxError);
    Errors[2]=sqrt(absSQRsum);
    Errors[3]=sqrt(absMaxError);

    return Errors;
}

int main(int argc, char* argv[]) {

    if(argc!=3) {
        std::cout << "Usage: " << argv[0] << " resolution(int) iteration(int)";
        return 0;
    }
    int N=atoi(argv[1]);
    int iter=atoi(argv[2]);
    
    if(N<1 || iter<1) {
        return 0;
    }
    
    double time;
    std::vector<double> u(N*N, 0), b, uS, uP, errs(4);  
    DirichletBC(u, N);
    b=forcing(N);
    std::cout << std::setprecision(15) << std::scientific;
    int myid, numprocs;

    
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs); 

    uS=JacobiMethod(u, b, N, iter, time, myid, numprocs);

    MPI_Finalize();
    
    uP=exactSolution(N);
    
    if (myid==0) {
        errs=Errors(uS, uP, b, N);
    std::cout << " resErrorEucl: " << errs[0] << " resErrorMax: " << errs[1] << " totalErrorEucl: " << errs[2] << " totalErrorMax: " << errs[3];
    std::cout << " time: " << time << " \n"; 
    }
    
    return 0;
}