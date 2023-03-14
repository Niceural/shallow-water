#include <iostream>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main(int argc, char** argv) {
    double dt = 0.0; double t = 0.0;
    int nx = 0; int ny = 0;
    int ic = 1;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("dt", po::value<double>(&dt)->default_value(0.1), "Time-step to use.")
        ("T", po::value<double>(&t)->default_value(100.0), "Total integration time.")
        ("Nx", po::value<int>(&nx)->default_value(100), "Number of grid points in x.")
        ("Ny", po::value<int>(&ny)->default_value(100), "Number of grid points in y.")
        ("ic", po::value<int>(&ic)->default_value(1), "Index of the initial condition to use (1-4).");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 1;
    }



    return 0;
}
