#include "mpi.h"
#include "lbm.hpp"
#include "System.hpp"

int main(int argc, char *argv[])
{
    int nx = 256;
    int ny = 64;
    int nz = 64;
    double start, end;
    MPI_Init(&argc, &argv);

    System s1(nx, ny, nz);
    s1.Initialize();
    s1.Monitor();
    MPI_Barrier(MPI_COMM_WORLD);
    LBM l1(MPI_COMM_WORLD, s1.sx, s1.sy, s1.sz, s1.tau, s1.rho0, s1.u0);

    l1.Initialize();

    l1.MPIoutput(0);

    for (int it = 1; it <= 10; it++)
    {
        // l1.Collision();
        l1.exchange();

        l1.Streaming();
        l1.Update();
        end = MPI_Wtime();
        if (it % 1 == 0)
        {
            l1.MPIoutput(it / 1);
        }
    }

    MPI_Finalize();

    return 0;
}