#ifndef LEAPFROG_TEST
#define LEAPFROG_TEST
#include "../include/leapfrog.hpp"
#include "./googletest/googletest/include/gtest/gtest.h"
#include <cmath>

TEST( leapfrog, half_step_velocity )
{
    // Initial velocity
    double v[2] = {0.12, 0.27};

    // energy gradient
    // U = 4x ** 2 + x * sin(y)
    // dux = 8x + sin(y)
    // duy = x * cos(y)
    // let x=2, y=5
    double grad[2] = {8*2 + std::sin(5), 2*std::cos(5)};
    double eps=0.2;

    // compute half step velocity
    double vnew[2];
    leapfrog::half_step_velocity( vnew, v, grad, 2, eps);

    EXPECT_DOUBLE_EQ( -1.3841075725336864, vnew[0] );
    EXPECT_DOUBLE_EQ( 0.21326756290735477, vnew[1] );
}

TEST( leapfrog, step_state )
{
    // Initial state
    double x[2] = {1, 2};

    // Velocity and step size
    double vel[2] = {0.2, -0.7};
    double eps = 0.5;

    // Compute the next state
    double xnew[2];
    leapfrog::step_state( xnew, x, vel, 2, eps );

    EXPECT_DOUBLE_EQ( 1.1, xnew[0] );
    EXPECT_DOUBLE_EQ( 1.65, xnew[1] );
}

TEST( leapfrog, lfs )
{
    // Initial state and velocity
    double x[2] = {2, 5};
    double v[2] = {0.12, 0.27};

    // step size
    double eps=0.2;

    // energy gradient
    // U = 4x ** 2 + x * sin(y)
    // dux = 8x + sin(y)
    // duy = x * cos(y)
    std::function<void(double*,const double*,const size_t)> f_grad =
        [](double * out, const double * in, const size_t)
        {
            out[0] = 8*in[0] + std::sin( in[1] );
            out[1] = in[0] * std::cos( in[1] );
        };

    // compute the next velocity and state
    double xnew[2], vnew[2], work[2];
    leapfrog::lfs( xnew, vnew, work, x, v, f_grad, 2, eps );

    ASSERT_DOUBLE_EQ( 1.7231784854932628, xnew[0] );
    ASSERT_DOUBLE_EQ( 5.0426535125814711, xnew[1] );
    ASSERT_DOUBLE_EQ( -2.668054701866847, vnew[0] );
    ASSERT_DOUBLE_EQ( 0.15738604333738995, vnew[1] );
}

#endif
