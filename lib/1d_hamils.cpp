#include "../include/1d_hamils.hpp"
#include "../include/constants.hpp"

void hmc::exchange_grad_1d(std::valarray<double>& grad_out,
                      const std::valarray<double>& data)
{
    int halfsize = data.size() / 2;
    std::slice tslice(0, halfsize, 1);
    std::slice pslice(halfsize, halfsize, 1);

    std::valarray<double> cos_theta = cos(data[tslice]);
    std::valarray<double> sin_theta = sin(data[tslice]);
    std::valarray<double> cos_phi = cos(data[pslice]);
    std::valarray<double> sin_phi = sin(data[pslice]);

    std::valarray<double> cos_theta_left = cos_theta.cshift(1);
    std::valarray<double> sin_theta_left = sin_theta.cshift(1);
    std::valarray<double> cos_phi_left = cos_phi.cshift(1);
    std::valarray<double> sin_phi_left = sin_phi.cshift(1);

    std::valarray<double> cos_theta_right = cos_theta.cshift(-1);
    std::valarray<double> sin_theta_right = sin_theta.cshift(-1);
    std::valarray<double> cos_phi_right = cos_phi.cshift(-1);
    std::valarray<double> sin_phi_right = sin_phi.cshift(-1);

    std::valarray<double> t1, t2, t3, t4, t5, t6;
    t1 = -cos_theta_left * sin_theta * sin_phi * sin_phi_left;
    t2 = cos_theta * sin_phi * sin_theta_left * sin_phi_left;
    t3 = -cos_theta_right * sin_theta * sin_phi * sin_phi_right;
    t4 = cos_theta * sin_phi * sin_theta_right * sin_phi_right;
    grad_out[tslice] += t1 + t2 + t3 + t4;

    t1 = -cos_phi_left * sin_phi;
    t2 = cos_theta * cos_phi * cos_theta_left * sin_phi_left;
    t3 = cos_phi * sin_theta * sin_theta_left * sin_phi_left;
    t4 = -cos_phi_right * sin_phi;
    t5 = cos_theta * cos_phi * cos_theta_right * sin_phi_right;
    t6 = cos_phi * sin_theta * sin_theta_right * sin_phi_right;
    grad_out[pslice] += t1 + t2 + t3 + t4 + t5 + t6;
}

std::function<void(std::valarray<double>&, const std::valarray<double>&)>
hmc::gen_total_grad_1d( const HamiltonianOptions options )
{
    std::function<void(std::valarray<double>&, const std::valarray<double>&)> init_f =
        [](std::valarray<double>& grad_out, const std::valarray<double>& data)
        {grad_out = 0;};

    std::function<void(std::valarray<double>&, const std::valarray<double>&)> new_f;
    if ( options.J != 0 )
    {
        new_f =
            [init_f](std::valarray<double>& grad_out, const std::valarray<double>& data)
            {
                init_f(grad_out, data);
                exchange_grad_1d(grad_out, data);
            };
        init_f = new_f;
    }

    return init_f;
}

double hmc::exchange_energy_1d(const std::valarray<double>& data)
{
    int halfsize = data.size() / 2;
    std::slice tslice(0, halfsize, 1);
    std::slice pslice(halfsize, halfsize, 1);

    std::valarray<double> cos_theta = cos(data[tslice]);
    std::valarray<double> sin_theta = sin(data[tslice]);
    std::valarray<double> cos_phi = cos(data[pslice]);
    std::valarray<double> sin_phi = sin(data[pslice]);

    std::valarray<double> cos_theta_left = cos_theta.cshift(1);
    std::valarray<double> sin_theta_left = sin_theta.cshift(1);
    std::valarray<double> cos_phi_left = cos_phi.cshift(1);
    std::valarray<double> sin_phi_left = sin_phi.cshift(1);

    double t1 = (cos_phi * cos_phi_left).sum();
    double t2 = (cos_theta * cos_theta_left * sin_phi * sin_phi_left).sum();
    double t3 = (sin_theta * sin_phi * sin_theta_left * sin_phi_left).sum();

    return -(t1 + t2 + t3);
}

std::function<double(const std::valarray<double>&)>
hmc::gen_total_energy_1d( const HamiltonianOptions options, const double beta )
{
    std::function<double(const std::valarray<double>&)> init_f =
        [](const std::valarray<double>& data){return 0;};

    std::function<double(const std::valarray<double>&)> new_f;
    if ( options.J != 0 )
    {
        new_f =
            [init_f, beta](const std::valarray<double>& data)
            {return init_f(data) + exchange_energy_1d( data ) * beta; };
        init_f = new_f;
    }

    return init_f;
}
