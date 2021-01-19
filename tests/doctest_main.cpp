
#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest/doctest.h>
#include <mpi.h>

int main(int argc, char** argv) {
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);
  doctest::Context context;

    context.setOption("order-by", "name");            // sort the test cases by their name
    context.applyCommandLine(argc, argv);
    context.setOption("no-breaks", true);             // don't break in the debugger when assertions fail
    int res = context.run(); // run
    if(context.shouldExit()) // important - query flags (and --exit) rely on the user doing this
        return res;          // propagate the result of the tests
    return res ; 
}
