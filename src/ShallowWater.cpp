#include "../include/ShallowWater.h"

ShallowWater::ShallowWater(
    const int nx, const int ny,
    const double dx, const double dy
):
    _dx(dx), _dy(dy),
    _grid(MultiQuantityMatrix(nx, ny, 9)),
    _fd(FiniteDifference(nx, ny, dx, dy))
{
    if (
        nx < 6 || ny < 6 ||
        dx < 0.0 || dy < 0.0
    ) throw std::invalid_argument("Invalid argument.");
}

void ShallowWater::setInitialConditions(const int ic, const double meanH) {
    // set U and V to nil
    for (int id = 0; id < _grid.size(); id++) {
        _grid.set(id, 0, 0.0);
        _grid.set(id, 3, 0.0);
    }

    // set H to the initial surface height for each test case
    switch (ic) {
        case 1:
            // plane waves propagating in x
            for (int i = 0; i < _grid.m(); i++) {
                double x = i * _dx;
                x -= 50.0;
                x *= x;
                double h = meanH + std::exp(-x / 25.0);
                for (int j = 0; j < _grid.n(); j++) {
                    _grid.set(i, j, 6, h);
                }
            }
            break;

        case 2:
            // plane waves propagating in y
            for (int j = 0; j < _grid.n(); j++) {
                double y = j * _dy;
                y -= 50.0;
                y *= y;
                double h = meanH + std::exp(-y / 25.0);
                for (int i = 0; i < _grid.m(); i++) {
                    _grid.set(i, j, 6, h);
                }
            }
            break;

        case 3:
            // single droplet
            for (int j = 0; j < _grid.n(); j++) {
                for (int i = 0; i < _grid.m(); i++) {
                    double x = _dx * i; double y = _dy * j;
                    x -= 50.0; y -= 50.0;
                    x *= x; y *= y;
                    double h = meanH + std::exp(-(x + y) / 25.0);
                    _grid.set(i, j, 6, h);
                }
            }
            break;

        case 4:
            // double droplet
            for (int j = 0; j < _grid.n(); j++) {
                for (int i = 0; i < _grid.m(); i++) {
                    double x = _dx * i; double y = _dy * j;
                    double h = meanH + 
                        std::exp(-((x-25.0)*(x-25.0) + (y-25.0)*(y-25.0)) / 25.0) +
                        std::exp(-((x-75.0)*(x-75.0) + (y-75.0)*(y-75.0)) / 25.0);
                    _grid.set(i, j, 6, h);
                }
            }
            break;

        default:
            throw std::invalid_argument("Invalid initial condition");
            break;
    }
}

void ShallowWater::timeIntegrate(const bool loopBlas, const double dt, const double t) {
    const int numPoints = _grid.size();

    MultiQuantityMatrix tgrid(_grid.m(), _grid.n(), 9);

    MultiQuantityMatrix k1(_grid.m(), _grid.n(), 3);
    MultiQuantityMatrix k2(_grid.m(), _grid.n(), 3);
    MultiQuantityMatrix k3(_grid.m(), _grid.n(), 3);
    MultiQuantityMatrix k4(_grid.m(), _grid.n(), 3);

    double ct = 0.0; // current time
    while (ct < t) {
        //--------- k1
        // finite difference
        _fd.centralDifference(loopBlas, _grid);
        std::cout << _grid.get(40, 60, 7) << std::endl;
        // k1
        for (int id = 0; id < numPoints; id++) {
            k1.set(id, 0, -_grid.get(id,0)*_grid.get(id,1) -_grid.get(id,3)*_grid.get(id,2) -CONST_G*_grid.get(id,7));
            k1.set(id, 1, -_grid.get(id,0)*_grid.get(id,4) -_grid.get(id,3)*_grid.get(id,5) -CONST_G*_grid.get(id,8));
            k1.set(id, 2, -_grid.get(id,0)*_grid.get(id,7) -_grid.get(id,6)*_grid.get(id,1) -_grid.get(id,3)*_grid.get(id,8) -_grid.get(id,6)*_grid.get(id,5));
        }

        //--------- k2
        // update U, V, and H
        for (int id = 0; id < numPoints; id++) {
            tgrid.set(id, 0, _grid.get(id,0) +0.5*dt*k1.get(id,0));
            tgrid.set(id, 3, _grid.get(id,3) +0.5*dt*k1.get(id,1));
            tgrid.set(id, 6, _grid.get(id,6) +0.5*dt*k1.get(id,2));
        }
        // central difference
        _fd.centralDifference(loopBlas, tgrid);
        // k2
        for (int id = 0; id < numPoints; id++) {
            k2.set(id, 0, -tgrid.get(id,0)*tgrid.get(id,1) -tgrid.get(id,3)*tgrid.get(id,2) -CONST_G*tgrid.get(id,7));
            k2.set(id, 1, -tgrid.get(id,0)*tgrid.get(id,4) -tgrid.get(id,3)*tgrid.get(id,5) -CONST_G*tgrid.get(id,8));
            k2.set(id, 2, -tgrid.get(id,0)*tgrid.get(id,7) -tgrid.get(id,6)*tgrid.get(id,1) -tgrid.get(id,3)*tgrid.get(id,8) -tgrid.get(id,6)*tgrid.get(id,5));
        }

        //--------- k3
        // update U, V, and H
        for (int id = 0; id < numPoints; id++) {
            tgrid.set(id, 0, _grid.get(id,0) +0.5*dt*k2.get(id,0));
            tgrid.set(id, 3, _grid.get(id,3) +0.5*dt*k2.get(id,1));
            tgrid.set(id, 6, _grid.get(id,6) +0.5*dt*k2.get(id,2));
        }
        // finite difference
        _fd.centralDifference(loopBlas, tgrid);
        // k3
        for (int id = 0; id < numPoints; id++) {
            k3.set(id, 0, -tgrid.get(id,0)*tgrid.get(id,1) -tgrid.get(id,3)*tgrid.get(id,2) -CONST_G*tgrid.get(id,7));
            k3.set(id, 1, -tgrid.get(id,0)*tgrid.get(id,4) -tgrid.get(id,3)*tgrid.get(id,5) -CONST_G*tgrid.get(id,8));
            k3.set(id, 2, -tgrid.get(id,0)*tgrid.get(id,7) -tgrid.get(id,6)*tgrid.get(id,1) -tgrid.get(id,3)*tgrid.get(id,8) -tgrid.get(id,6)*tgrid.get(id,5));
        }

        //--------- k4
        // update U, V, and H
        for (int id = 0; id < numPoints; id++) {
            tgrid.set(id, 0, _grid.get(id,0) +dt*k3.get(id,0));
            tgrid.set(id, 3, _grid.get(id,3) +dt*k3.get(id,1));
            tgrid.set(id, 6, _grid.get(id,6) +dt*k3.get(id,2));
        }
        // finite difference
        _fd.centralDifference(loopBlas, tgrid);
        // k4
        for (int id = 0; id < numPoints; id++) {
            k4.set(id, 0, -tgrid.get(id,0)*tgrid.get(id,1) -tgrid.get(id,3)*tgrid.get(id,2) -CONST_G*tgrid.get(id,7));
            k4.set(id, 1, -tgrid.get(id,0)*tgrid.get(id,4) -tgrid.get(id,3)*tgrid.get(id,5) -CONST_G*tgrid.get(id,8));
            k4.set(id, 2, -tgrid.get(id,0)*tgrid.get(id,7) -tgrid.get(id,6)*tgrid.get(id,1) -tgrid.get(id,3)*tgrid.get(id,8) -tgrid.get(id,6)*tgrid.get(id,5));
        }

        // update U, V, and H
        for (int id = 0; id < numPoints; id++) {
            _grid.add(id, 0, (k1.get(id,0) +2.0*k2.get(id,0) +2.0*k3.get(id,0) +k4.get(id,0)) *dt/6.0);
            _grid.add(id, 3, (k1.get(id,1) +2.0*k2.get(id,1) +2.0*k3.get(id,1) +k4.get(id,1)) *dt/6.0);
            _grid.add(id, 6, (k1.get(id,2) +2.0*k2.get(id,2) +2.0*k3.get(id,2) +k4.get(id,2)) *dt/6.0);
        }

        ct += dt;
    }
}

void ShallowWater::exportGrid(const std::string& fname) {
    std::ofstream file;
    file.open(fname);

    for (int j = 0; j < _grid.n(); j++) {
        double y = j * _dy;
        for (int i = 0; i < _grid.m(); i++) {
            double x = i * _dx;
            file << x << " " << y << " ";
            file << _grid.get(i, j, 0) << " ";
            file << _grid.get(i, j, 3) << " ";
            file << _grid.get(i, j, 6) << std::endl;
        }
        file << std::endl;
    }

    file.close();
}

void ShallowWater::test() {
    std::string fname = "output.txt";
    std::ofstream file;
    file.open(fname);

    for (int i = 0; i < _grid.m(); i++) {
        for(int j = 0; j < _grid.n(); j++) {
            _grid.set(i, j, 6, j*2.0);
        }
    }

    _fd.centralDifference(0, _grid);

    for (int j = 0; j < _grid.n(); j++) {
        double y = j * _dy;
        for (int i = 0; i < _grid.m(); i++) {
            double x = i * _dx;
            file << x << " " << y << " ";
            file << _grid.get(i, j, 0) << " ";
            file << _grid.get(i, j, 3) << " ";
            file << _grid.get(i, j, 8) << std::endl;
        }
        file << std::endl;
    }

    file.close();
}
