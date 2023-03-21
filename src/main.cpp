#include <iostream>
#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>
#include "../include/ShallowWater.h"

namespace po = boost::program_options;

int main(int argc, char** argv) {
    boost::timer::auto_cpu_timer timer("%t sec CPU, %w sec real\n");
    double dt = 0.1; double t = 80.0;
    int nx = 100; int ny = 100;
    int ic = 3; bool loopBlas = 0;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("dt", po::value<double>(&dt)->default_value(0.1), "Time-step to use.")
        ("T", po::value<double>(&t)->default_value(100.0), "Total integration time.")
        ("Nx", po::value<int>(&nx)->default_value(100), "Number of grid points in x.")
        ("Ny", po::value<int>(&ny)->default_value(100), "Number of grid points in y.")
        ("ic", po::value<int>(&ic)->default_value(3), "Index of the initial condition to use (1-4).");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 1;
    }

    ShallowWater sw(nx, ny, 1.0, 1.0);
    sw.setInitialConditions(ic, 10.0);
    sw.timeIntegrate(loopBlas, dt, t);
    // sw.test();
    sw.exportGrid("output.txt");

    return 0;
}
