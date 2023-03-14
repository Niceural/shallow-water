#define BOOST_TEST_MODULE InitialCondition
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include "../include/ShallowWater.h"
#include "../include/matrices/GeneralMatrix.h"
namespace utf = boost::unit_test;
#define TOLERANCE 0.000001

BOOST_AUTO_TEST_CASE(TestCase1, *utf::tolerance(TOLERANCE))
{
    const double dt = 0.84234; const double t = 232.2739;
    const int nx = 546; const int ny = 957;
    const int ic = 1;
    ShallowWater sw = ShallowWater(dt, t, nx, ny, ic);
    sw.setInitialConditions();

    GeneralMatrix h = sw.getH();

    // values generated with MATLAB
    const double expected[] = {};

    for (int i = 0; i < h.size(); i++) {
        BOOST_TEST(h[i] == expected[i]);
    }
}

BOOST_AUTO_TEST_CASE(TestCase2, *utf::tolerance(TOLERANCE))
{
    const double dt = 0.2079437; const double t = 272.274;
    const int nx = 964; const int ny = 157;
    const int ic = 2;
    ShallowWater sw = ShallowWater(dt, t, nx, ny, ic);
    sw.setInitialConditions();

    GeneralMatrix h = sw.getH();

    // values generated with MATLAB
    const double expected[] = {};

    for (int i = 0; i < h.size(); i++) {
        BOOST_TEST(h[i] == expected[i]);
    }
}

BOOST_AUTO_TEST_CASE(TestCase3, *utf::tolerance(TOLERANCE))
{
    const double dt = 0.2479271; const double t = 975.9274;
    const int nx = 485; const int ny = 800;
    const int ic = 3;
    ShallowWater sw = ShallowWater(dt, t, nx, ny, ic);
    sw.setInitialConditions();

    GeneralMatrix h = sw.getH();

    // values generated with MATLAB
    const double expected[] = {};

    for (int i = 0; i < h.size(); i++) {
        BOOST_TEST(h[i] == expected[i]);
    }
}

BOOST_AUTO_TEST_CASE(TestCase4, *utf::tolerance(TOLERANCE))
{
    const double dt = 2.9743; const double t = 648.2957;
    const int nx = 141; const int ny = 421;
    const int ic = 4;
    ShallowWater sw = ShallowWater(dt, t, nx, ny, ic);
    sw.setInitialConditions();

    GeneralMatrix h = sw.getH();

    // values generated with MATLAB
    const double expected[] = {};

    for (int i = 0; i < h.size(); i++) {
        BOOST_TEST(h[i] == expected[i]);
    }
}
