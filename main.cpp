#include <iostream>
#include <vector>
#include <cmath>

#include <fstream>
#include <functional>

//template <typename T>
//class GridI {
//
//};
//
//template <typename T>
//class Grid {
//
//};
//
//template <typename T>
//class Interpolation {
//
//};
//
//template <typename T>
//class GridFunction {
//
//    Interpolation<T>* interpolation;
//    GridI<T>* grid
//};


// Сначала обусловимся на равномерную сетку
const int M = 1000;
const float L = 1;
const float h = L/M;

const int N = 5000;
const float T = 1;
const float t = T/N;

const std::string PATH = "../data";

float calc_f(float u){
    return u*u/2;
}

// Пока - схема Куранта-Изаксона-Рис
// Краевые условия - нулевая производная
void progress_KIR(const std::vector<float>& u0, std::vector<float>& u1){
    for(int i=0; i<M; i++){
        // approximate derivative
        auto dfdu = u0[i];
//        if (fabs(dfdu)*t/h > 1) std::cout << "Breaching stability on x=" << float(i)*h << '\n';
        if (dfdu > 0) {
            if(i != 0) u1[i] = u0[i] - t/h*(calc_f(u0[i]) - calc_f(u0[i-1]));
            else u1[i] = u0[i];
        } else {
            if(i != M) u1[i] = u0[i] - t/h*(calc_f(u0[i+1]) - calc_f(u0[i]));
            else u1[i] = u0[i];
        }
    }
}

// Лакса-Вендроффа
void progress_LW(const std::vector<float>& u0, std::vector<float>& u1){
    for(int i=0; i<M; i++){
        //predictor
        auto u0c = u0[i];
        auto u0r = i != M ? u0[i+1] : u0c;
        auto u0l = i != 0 ? u0[i-1] : u0c;

        auto ur = (u0c + u0r)/2 - t/2/h * (calc_f(u0r) - calc_f(u0c));
        auto ul = (u0c + u0l)/2 - t/2/h * (calc_f(u0c) - calc_f(u0l));
        // corrector
        u1[i] = u0c - t/h * (calc_f(ur) - calc_f(ul));
    }
}

void write(std::ofstream& fout, const std::vector<float>& u, float time){
    fout << time;
    for(int i=0; i<M; i++){
        fout << ',' << u[i];
    }
    fout << '\n';
}

void write2file(std::vector<float> u0, // we need copy there
                const std::function<void(const std::vector<float>&, std::vector<float>&)>& solver,
                const std::string& path)
{
    std::ofstream fout;
    fout.open(PATH + path);

    std::vector<float> u1; // buffer
    u1.resize(M);

//    std::cout << "-- Filling x-es top row\n";
    fout << 0;
    for (int i = 0; i < M; i++)
        fout << ',' << i * h;
    fout << '\n';

    for(int j=0; (2*j)<N; j++){
        if (j % (N*5/200) == 0) std::cout << (j / (N/200)) << "%\n";
        solver(u0, u1);
        write(fout, u1, float(2*j)*t);
        solver(u1, u0);
        write(fout, u0, float(2*j+1)*t);
    }

    fout.close();
}

int main() {
    std::vector<float> u0; // input

    u0.resize(M);

    std::cout << "Line-up calculation\n";

    for (int i = 0; i < M; i++)
        u0[i] = 1 + i * h;

    write2file(u0, progress_KIR, "/line_up/dataKIR.csv");
    write2file(u0, progress_LW, "/line_up/dataLW.csv");

    std::cout << "Sin calculation\n";

    for (int i = 0; i < M; i++)
        u0[i] = 1 + 2.f/3 * sinf(M_PI * 2 * i / M);

    write2file(u0, progress_KIR, "/sin/dataKIR.csv");
    write2file(u0, progress_LW, "/sin/dataLW.csv");


    std::cout << "Step-up calculation\n";

    for (int i = 0; i < M; i++)
        if (2*i > M) u0[i] = 2;
        else u0[i] = 1;

    write2file(u0, progress_KIR, "/stepup/dataKIR.csv");
    write2file(u0, progress_LW, "/stepup/dataLW.csv");


    std::cout << "Step-down calculation\n";

    for (int i = 0; i < M; i++)
        if (2*i < M) u0[i] = 2;
        else u0[i] = 1;

    write2file(u0, progress_KIR, "/stepdown/dataKIR.csv");
    write2file(u0, progress_LW, "/stepdown/dataLW.csv");

    return 0;
}
