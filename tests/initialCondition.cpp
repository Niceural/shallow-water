#define BOOST_TEST_MODULE InitialCondition
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include "../include/ShallowWater.h"
namespace utf = boost::unit_test;
#define TOLERANCE 0.000001

BOOST_AUTO_TEST_CASE(TestCase1, *utf::tolerance(TOLERANCE))
{
    ShallowWater sw = ShallowWater(0.1, 13.5, 47, 108, 1);
    double exp, act;

    // check u(x, y, 0) = 0
    exp = 0.0;
    for (int i = 0; i < sw.get_nx(); i++) {
        for (int j = 0; j < sw.get_ny(); j++) {
            act = sw.get_u(i, j);
            BOOST_TEST(act == exp);
        }
    }

    // check v(x, y, 0) = 0
    exp = 0.0;
    for (int i = 0; i < sw.get_nx(); i++) {
        for (int j = 0; j < sw.get_ny(); j++) {
            act = sw.get_v(i, j);
            BOOST_TEST(act == exp);
        }
    }

    // check h(x, y, 0) = exp(-(x-50)^2 / 25)
    double dx = sw.get_dx();
    for (int i = 0; i < sw.get_nx(); i++) {
        double x = dx * i;
        exp = 10.0 + std::exp(-(x-50.0)*(x-50.0) / 25.0);
        for (int j = 0; j < sw.get_ny(); j++) {
            act = sw.get_h(i, j);
            BOOST_TEST(act == exp);
        }
    }
}

BOOST_AUTO_TEST_CASE(TestCase2, *utf::tolerance(TOLERANCE))
{
    ShallowWater sw = ShallowWater(0.1, 13.5, 47, 108, 2);
    double exp, act;

    // check u(x, y, 0) = 0
    exp = 0.0;
    for (int i = 0; i < sw.get_nx(); i++) {
        for (int j = 0; j < sw.get_ny(); j++) {
            act = sw.get_u(i, j);
            BOOST_TEST(act == exp);
        }
    }

    // check v(x, y, 0) = 0
    exp = 0.0;
    for (int i = 0; i < sw.get_nx(); i++) {
        for (int j = 0; j < sw.get_ny(); j++) {
            act = sw.get_v(i, j);
            BOOST_TEST(act == exp);
        }
    }

    // check h(x, y, 0) = 10 + exp(-(x-50)^2 / 25)
    double dy = sw.get_dy();
    for (int j = 0; j < sw.get_ny(); j++) {
        double y = dy * j;
        exp = 10.0 + std::exp(-(y-50.0)*(y-50.0) / 25.0);
        for (int i = 0; i < sw.get_nx(); i++) {
            act = sw.get_h(i, j);
            BOOST_TEST(act == exp);
        }
    }
}