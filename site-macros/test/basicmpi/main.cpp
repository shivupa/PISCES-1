#include <iostream>
#include <mpi.h>

int main () {
   MPI_Init(0,0);
   
   int rank, size; 
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   using namespace std;
   cout << rank << " of " << size << endl;

   MPI_Finalize();
   return 0;
}
