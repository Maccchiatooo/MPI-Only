#ifndef _LBM_H_
#define _LBM_H_

#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define q 27
#define dim 3

struct CommHelper
{

    MPI_Comm comm;
    int rx, ry, rz;
    int me;
    int px, py, pz;
    int up, down, left, right, front, back, frontup, frontdown, frontleft, frontright, frontleftup, frontleftdown, frontrightup, frontrightdown, backup, backdown, backleft, backright, backleftup, backrightup, backleftdown, backrightdown, leftup, leftdown, rightup, rightdown;

    CommHelper(MPI_Comm comm_)
    {
        comm = comm_;
        int nranks;
        MPI_Comm_size(comm, &nranks);
        MPI_Comm_rank(comm, &me);

        rx = std::pow(1.0 * nranks, 1.0 / 3.0);
        while (nranks % rx != 0)
            rx++;

        rz = std::sqrt(1.0 * (nranks / rx));
        while ((nranks / rx) % rz != 0)
            rz++;

        ry = nranks / rx / rz;

        // printf("rx=%i,ry=%i,rz=%i\n", rx, ry, rz);
        px = me % rx;
        pz = (me / rx) % rz;
        py = (me / rx / rz);

        left = px == 0 ? -1 : me - 1;
        leftup = (px == 0 || pz == rz - 1) ? -1 : me - 1 + rx;
        rightup = (px == rx - 1 || pz == rz - 1) ? -1 : me + 1 + rx;
        leftdown = (px == 0 || pz == 0) ? -1 : me - 1 - rx;
        rightdown = (px == rx - 1 || pz == 0) ? -1 : me + 1 - rx;
        right = px == rx - 1 ? -1 : me + 1;
        down = pz == 0 ? -1 : me - rx;
        up = pz == rz - 1 ? -1 : me + rx;

        front = py == 0 ? -1 : me - rx * rz;
        frontup = (py == 0 || pz == rz - 1) ? -1 : me - rx * rz + rx;
        frontdown = (py == 0 || pz == 0) ? -1 : me - rx * rz - rx;
        frontleft = (py == 0 || px == 0) ? -1 : me - rx * rz - 1;
        frontright = (py == 0 || px == rx - 1) ? -1 : me - rx * rz + 1;
        frontleftdown = (py == 0 || px == 0 || pz == 0) ? -1 : me - rx * rz - rx - 1;
        frontrightdown = (py == 0 || px == rx - 1 || pz == 0) ? -1 : me - rx * rz - rx + 1;
        frontrightup = (py == 0 || px == rx - 1 || pz == rz - 1) ? -1 : me - rx * rz + rx + 1;
        frontleftup = (py == 0 || px == 0 || pz == rz - 1) ? -1 : me - rx * rz + rx - 1;

        back = py == ry - 1 ? -1 : me + rx * rz;
        backup = (py == ry - 1 || pz == rz - 1) ? -1 : me + rx * rz + rx;
        backdown = (py == ry - 1 || pz == 0) ? -1 : me + rx * rz - rx;
        backleft = (py == ry - 1 || px == 0) ? -1 : me + rx * rz - 1;
        backright = (py == ry - 1 || px == rx - 1) ? -1 : me + rx * rz + 1;
        backleftdown = (py == ry - 1 || px == 0 || pz == 0) ? -1 : me + rx * rz - rx - 1;
        backrightdown = (py == ry - 1 || px == rx - 1 || pz == 0) ? -1 : me + rx * rz - rx + 1;
        backrightup = (py == ry - 1 || px == rx - 1 || pz == rz - 1) ? -1 : me + rx * rz + rx + 1;
        backleftup = (py == ry - 1 || px == 0 || pz == rz - 1) ? -1 : me + rx * rz + rx - 1;

        MPI_Barrier(MPI_COMM_WORLD);
    }
    template <class ViewType>
    void isend_irecv(int partner, ViewType send_buffer, ViewType recv_buffer, MPI_Request *request_send, MPI_Request *request_recv)
    {
        MPI_Irecv(recv_buffer.data(), recv_buffer.size(), MPI_DOUBLE, partner, 1, comm, request_recv);
        MPI_Isend(send_buffer.data(), send_buffer.size(), MPI_DOUBLE, partner, 1, comm, request_send);
    }
};
struct LBM
{

    CommHelper comm;

    int mpi_active_requests;

    int glx;
    int gly;
    int glz;
    int ghost = 3;
    // exact nodes
    // 256/1,64/2,64/1
    int ex = glx / comm.rx;
    int ey = gly / comm.ry;
    int ez = glz / comm.rz;
    // include ghost nodes
    // 256+2*3,32+2*3,64+2*3
    int lx = ex + 2 * ghost;
    int ly = ey + 2 * ghost;
    int lz = ez + 2 * ghost;
    // 3,3,3  262,35,67
    // 256,32,64
    int l_s[3] = {ghost, ghost, ghost};
    int l_e[3] = {lx - ghost, ly - ghost, lz - ghost};

    int x_lo, x_hi, y_lo, y_hi, z_lo, z_hi;
    double rho0 = 1.0;
    double mu;
    double cs2 = 1.0 / 3.0;
    double tau0 = 0.12;
    double u0 = 0.1;

    double *f, *ft, *fb, *ua, *va, *wa, *rho, *p, *t;
    int *e, *usr, *ran, *bb;

    LBM(MPI_Comm comm_, int sx, int sy, int sz, double &tau, double &rho0, double &u0) : comm(comm_), glx(sx), gly(sy), glz(sz), tau0(tau), rho0(rho0), u0(u0)
    {
        ghost = 3;

        f = (double *)malloc(sizeof(double) * q * lx * ly * lz);
        ft = (double *)malloc(sizeof(double) * q * lx * ly * lz);
        fb = (double *)malloc(sizeof(double) * q * lx * ly * lz);

        ua = (double *)malloc(sizeof(double) * lx * ly * lz);
        va = (double *)malloc(sizeof(double) * lx * ly * lz);
        wa = (double *)malloc(sizeof(double) * lx * ly * lz);

        rho = (double *)malloc(sizeof(double) * lx * ly * lz);
        p = (double *)malloc(sizeof(double) * lx * ly * lz);

        e = (int *)malloc(sizeof(int) * q * dim);
        t = (double *)malloc(sizeof(double) * q);
        usr = (int *)malloc(sizeof(int) * lx * ly * lz);
        ran = (int *)malloc(sizeof(int) * lx * ly * lz);
        bb = (int *)malloc(sizeof(int) * q);

        x_lo = ex * comm.px;
        x_hi = ex * (comm.px + 1) - 1;
        y_lo = ey * comm.py;
        y_hi = ey * (comm.py + 1) - 1;
        z_lo = ez * comm.pz;
        z_hi = ez * (comm.pz + 1) - 1;
    };

    void Initialize();
    void Collision();
    void setup_subdomain();
    void pack();
    void exchange();
    void unpack();
    void Streaming();
    void Update();
    void MPIoutput(int n);
    void Output(int n);
};
#endif