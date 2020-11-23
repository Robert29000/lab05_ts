//
//  main.cpp
//  lab05
//
//  Created by Роберт Артур Меликян on 21/11/2020.
//  Copyright © 2020 Роберт Артур Меликян. All rights reserved.
//

#include <iostream>
#include <vector>
#include <random>


const int N = 20;
std::vector<std::pair<double, double>> values(20); // (x_i, t_i)

double func(double c, double d){
    double res = 0;
    for (auto i = 0; i < N; i++){
        auto x_i = values[i].first;
        auto t_i = values[i].second;
        double y_i = c * x_i + d;
        res += (y_i - t_i) * (y_i - t_i);
    }
    return res;
}

double dihotomia_method(double a_k, double b_k, double c){
    double d_l, d_r, f_l, f_r;
    const double eps = 0.1;
    const double delta = 0.01;
    do {
        d_l = 0.5 * (b_k + a_k) - delta;
        d_r = 0.5 * (b_k + a_k) + delta;
        f_l = func(c, d_l);
        f_r = func(c, d_r);
        if (f_r > f_l) {
            b_k = d_r;
        }
        else {
            a_k = d_l;
        }
    } while ((b_k - a_k) > eps);
    return ((b_k + a_k) / 2);
}

double golden_ratio_method(double a_k, double b_k){
    double l = std::abs(b_k - a_k);
    std::swap(a_k, b_k);
    a_k = std::fabs(a_k);
    b_k = std::fabs(b_k);
    const double e = 0.1;
    const double t = (std::sqrt(5) + 1) / 2;
    double c_k1 = a_k + (1 - 1/t)*b_k;
    double d_k1 = dihotomia_method(-5, 5, -c_k1);
    double c_k2 = a_k + b_k / t;
    double d_k2 = dihotomia_method(-5, 5, -c_k2);
    double f_k1 = func(-c_k1, d_k1);
    double f_k2 = func(-c_k2, d_k2);
    while (l > e){
        if (f_k1 < f_k2){
            b_k = c_k2;
            c_k2 = a_k + b_k - c_k1;
            d_k2 = dihotomia_method(-5, 5, -c_k2);
            f_k2 = func(-c_k2, d_k2);
        } else {
            a_k = c_k1;
            c_k1 = a_k + b_k - c_k2;
            d_k1 = dihotomia_method(-5, 5, -c_k1);
            f_k1 = func(-c_k1, d_k1);
        }
        if (c_k1 > c_k2){
            std::swap(c_k1, c_k2);
            std::swap(d_k1, d_k2);
            std::swap(f_k1, f_k2);
        }
        l = std::abs(b_k - a_k);
    }
    return -((b_k + a_k) / 2);
}


int main(int argc, const char * argv[]) {
    const double C = -10;
    const double D = 0;
    const double A = 10;
    const double a = -5;
    const double b = 0;
    auto step = (b - a)/N;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(-0.5, 0.5);
    int count = 0;
    for (double x = a; x < b; x += step){
        auto t_i = C * x + D + A * dis(gen);  // with errors
        //auto t_i = C * x + D;  // without erros
        values[count] = {x, t_i};
        count++;
    }
    

    auto c_r = golden_ratio_method(-20, -1);
    auto d_r = dihotomia_method(-5, 5, c_r);
    std::cout << c_r << " " << d_r << std::endl;
    for (auto i = 0; i < N; i++){
        std::cout << "(" << values[i].first << ";" << values[i].second << ") ";
    }
    std::cout << std::endl;
    return 0;
}
