#include <mpi_dispatcher/mpi_dispatcher.hpp>
#include <thread>
#include <random>
#include <iostream>

using namespace pMPI;

int dumb_task_counter;

void dumb_task(double seconds, int jobid, int rank) {
    std::cout << "[" << rank << "] running job " << jobid << " " << seconds << " seconds..." << std::flush;
    std::this_thread::sleep_for(std::chrono::milliseconds(int(seconds * 1000)));
    ++dumb_task_counter;
    std::cout << "done." << std::endl;
};

int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);

    std::random_device rd;
    boost::mpi::communicator world;
    std::mt19937 gen(100000);
    std::uniform_real_distribution<double> dist(0, 0.001);
    size_t ROOT = 0;
    int rank = world.rank();

    try {

        MPIWorker worker(world, ROOT);
        int ntasks = 7;
        dumb_task_counter = 0;

        std::unique_ptr<MPIMaster> disp;

        if (world.rank() == ROOT) {
            disp.reset(new MPIMaster(world, ntasks, true));
            disp->order();
            std::cout << "ordered" << std::endl;
        };
        world.barrier();

        for (; !worker.is_finished();) {
            if (rank == ROOT) disp->order();
            worker.receive_order();
            if (worker.is_working()) {
                dumb_task(dist(gen), worker.current_job(), rank);
                worker.report_job_done();
            };
            if (rank == ROOT) std::cout << "--> stack size = " << disp->JobStack.size() << " --> worker stack size =" << disp->WorkerStack.size() << std::endl;
            if (rank == ROOT) disp->check_workers();
        };
        if (rank == ROOT) disp.release();

        world.barrier();
        MPI_Allreduce(MPI_IN_PLACE, &dumb_task_counter, 1, MPI_INT, MPI_SUM, world);
        if(dumb_task_counter != ntasks) {
            std::cout << "ntasks = " << ntasks
                      << ", dumb_task_counter = "
                      << dumb_task_counter << std::endl;
            return EXIT_FAILURE;
        }
    } // end try
    catch (std::exception &e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    };

    world.barrier();
/*
    if (world.size() > 1) {
        int t = 10;
        //DEBUG(world.rank() << " " << t);
        if (world.rank() == 0) {
            MPI_Request req;
            t = 45;
            //world.isend(1, 1, t);
            MPI_Isend(&t, 1, MPI_INT, 1, 1, world, &req);
        }
        else if (world.rank() == 1) {
            MPI_Request req;
            MPI_Irecv(&t, 1, MPI_INT, 0, 1, world, &req); //  = world.irecv(0, 1, t);
            MPI_Wait(&req, MPI_STATUS_IGNORE);
        };
        world.barrier();
        std::cout << world.rank() << " " << t << std::endl;
    };

// do it again
    try {

        MPIWorker worker(world, ROOT);
        int ntasks = 9;
        dumb_task_counter = 0;

        std::unique_ptr<MPIMaster> disp;

        if (world.rank() == ROOT) {
            disp.reset(new MPIMaster(world, ntasks, true));
            disp->order();
        };
        world.barrier();

        for (; !worker.is_finished();) {
            if (rank == ROOT) disp->order();
            worker.receive_order();
            if (worker.is_working()) {
                dumb_task(dist(gen), 0);
                worker.report_job_done();
            };
            if (rank == ROOT) disp->check_workers();
        };
        std::cout << "Huy! " << std::endl;
        if (rank == ROOT) disp.release();


    } // end try
    catch (std::exception &e) { return EXIT_FAILURE; };
*/

    MPI_Finalize();
    return EXIT_SUCCESS;
}