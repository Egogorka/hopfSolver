#include <iostream>
#include <vector>
#include <cmath>

#include <fstream>
#include <functional>


struct Grid {
    double L; // Length
    double T; // Time
    int M; // Length discretion (amount of intervals, i.e. M+1 - amount of points)
    int N; // Time discretion

    [[nodiscard]] double h() const { return L/double(M); }
    [[nodiscard]] double t() const { return T/double(N); }

    void write(std::ostream& out) const {
        out << L << ',' << M << ',' << T << ',' << N << '\n';
    }
};

struct FuncGrid {
    Grid grid{};
    std::vector<std::vector<double>> func;

    explicit FuncGrid(const Grid& grid): grid(grid), func(){
        func.resize(grid.N+1, {});
        for(auto& subvec: func){
            subvec.resize(grid.M+1, {});
        }
    }


    void write(std::ostream& out) const {
        grid.write(out);
        for(int i=0; i<=grid.N; i++){
            for(int j=0; j<=grid.M; j++){
                if (j != 0) out << ',';
                out << func[i][j];
            }
            out << '\n';
        }
    }
};

//class Scheme {
//public:
//    /**
//     * Atm only 0 derivative border conditions are implemented
//     * TODO: implement other border conditions
//     */
//
//    /**
//     * @param fg - FunctionGrid to be solved. fg[0][i] must be filled with initial condition
//     */
//    void solve(FuncGrid& fg) = 0;
//};

double calc_f(double u){
    return u*u/2;
}

class Scheme_KIR {
public:
    static void solve(FuncGrid& fg, bool border_present = false) {
        std::cout << "KIR scheme start\n";
        for(int j=1; j<=fg.grid.N; j++){
            if (fg.grid.N > 10 && j % (fg.grid.N/10) == 0) std::cout << "KIR : Finished " << (j*10) / fg.grid.N * 10 << "%\n";
            solve_layer(fg, j, border_present);
        }
        std::cout << "KIR scheme end\n";
    }

protected:
    /**
     * @param fg
     * @param layer - layer to be filled
     */
    static void solve_layer(FuncGrid& fg, int layer, bool border_present = false) {
        auto t = fg.grid.t();
        auto h = fg.grid.h();
        auto M = fg.grid.M;

        std::vector<double>& u0 = fg.func[layer-1];
        std::vector<double>& u1 = fg.func[layer];

        if (!border_present){
            u1[0] = u0[0];
            u1[M] = u0[M];
        }

        for(int i=1; i<M; i++){
            // approximate derivative
            auto dfdu = u0[i];

            // if (fabs(dfdu)*t/h > 1) std::cout << "Breaching stability on x=" << double(i)*h << '\n';
            if (dfdu > 0) {
                u1[i] = dfdu - t/h*(calc_f(u0[i]) - calc_f(u0[i-1]));
            } else {
                u1[i] = dfdu - t/h*(calc_f(u0[i+1]) - calc_f(u0[i]));
            }
        }
    }
};

class Scheme_LW {
protected:
    static void solve_layer(FuncGrid& fg, int layer, bool border_present = false) {
        auto t = fg.grid.t();
        auto h = fg.grid.h();
        auto M = fg.grid.M;

        std::vector<double>& u0 = fg.func[layer-1];
        std::vector<double>& u1 = fg.func[layer];

        if (!border_present){
            u1[0] = u0[0];
            u1[M] = u0[M];
        }

        for(int i=1; i<M; i++){
            //predictor
            auto u0c = u0[i];
            auto u0r = u0[i+1];
            auto u0l = u0[i-1];

            auto ur = (u0c + u0r)/2 - t/2/h * (calc_f(u0r) - calc_f(u0c));
            auto ul = (u0c + u0l)/2 - t/2/h * (calc_f(u0c) - calc_f(u0l));
            // corrector
            u1[i] = u0c - t/h * (calc_f(ur) - calc_f(ul));
        }
    }

public:
    static void solve(FuncGrid& fg, bool border_present = false) {
        std::cout << "LW scheme start\n";
        for(int j=1; j<=fg.grid.N; j++){
            if (fg.grid.N > 10 && j % (fg.grid.N/10) == 0) std::cout << "LW : Finished " << (j*10) / fg.grid.N * 10 << "%\n";
            solve_layer(fg, j, border_present);
        }
        std::cout << "LW scheme end\n";
    }
};


// Сначала обусловимся на равномерную сетку
const int M = 100;
const double L = 1;

const int N = 500;
const double T = 1;

const std::string PATH = "../data";

double step(double a, double b, double x0, double x, double t){
    if (a < b) {
        if (x < x0 + a * t) return a;
        if (x >= x0 + b * t) return b;
        return (x - x0) / t;
    } else {
        if (x < x0 + (a+b)/2 * t) return a;
        return b;
    }
}

/**
 * @param a = u(x=0, t=0)
 * @param b = u(x=L, t=0)
 * @param x
 * @param t
 * @return u(x,t)
 */
double line(double a, double b, double L, double x, double t){
    return (L*a + (b-a)*x)/(L+(b-a)*t);
}

double soliton(double C, double x, double t){
    auto w = x/sqrt(t);
    return (w + sqrt(w*w - 4*C))/2/sqrt(t);
}


void write2file(const std::string& path, const std::function<void(std::ofstream&)>& f){
    std::ofstream fout;
    fout.open(PATH + path);
    f(fout);
    fout.close();
}

void calc_func(const Grid& grid, const std::string& path, const std::function<double(double,double)>& f, bool border = false){
    std::cout << "Step-up calculation for " << path << "\n";
    std::cout << "Grid parameters:";
    grid.write(std::cout);
    FuncGrid fg0(grid), fg_KIR(grid), fg_LW(grid);

    // Filling the analitycal solution
    for(int j=0; j<=grid.N; j++)
        for(int i=0; i<=grid.M; i++)
//            fg0.func[j][i] = step(1, 2, 0.5, double(i)*grid.h(), double(j)*grid.t());
            fg0.func[j][i] = f(double(i)*grid.h(), double(j)*grid.t());

    // Filling initial condition
    for(int i=0; i<=grid.M; i++)
        fg_LW.func[0][i] = fg_KIR.func[0][i] = f(double(i)*grid.h(), 0);

    if (border) {
        for (int i = 0; i <= grid.N; i++) {
            fg_LW.func[i][0] = fg_KIR.func[i][0] = f(0, double(i)*grid.t());
            fg_LW.func[i][grid.M] = fg_KIR.func[i][grid.M] = f(grid.L, double(i)*grid.t());
        }
    }

    Scheme_KIR::solve(fg_KIR, border);
    Scheme_LW::solve(fg_LW, border);

    write2file(path+"_dataKIR.csv", [&](std::ofstream& fout){fg_KIR.write(fout);});
    write2file(path+"_dataLW.csv", [&](std::ofstream& fout){fg_LW.write(fout);});
    write2file(path+"_dataAnalytic.csv", [&](std::ofstream& fout){fg0.write(fout);});
}

int main() {

    auto step_up = [&](double x, double t){
        return step(1, 2, 0.5, x, t);
    };

    auto line_up = [&](double x, double t){
        return step(1, 2, - 1*1.f/(2-1), x, t + 1.f/(2-1));
    };

    auto soliton_up = [&](double x, double t) {
        return soliton(-0.1, x, t+1.1f);
    };

//    step_up_calc({.L = 1, .T = 0.5, .M = 10, .N = 50}, "/stepup/1e1.0");
//    step_up_calc({.L = 1, .T = 0.5, .M = 20, .N = 100}, "/stepup/1e1.3");
//    step_up_calc({.L = 1, .T = 0.5, .M = 50, .N = 250}, "/stepup/1e1.6");
//    step_up_calc({.L = 1, .T = 0.5, .M = 100, .N = 500}, "/stepup/1e2.0");
//    step_up_calc({.L = 1, .T = 0.5, .M = 200, .N = 1000}, "/stepup/1e2.3");
//    step_up_calc({.L = 1, .T = 0.5, .M = 500, .N = 2500}, "/stepup/1e2.6");
//    step_up_calc({.L = 1, .T = 0.5, .M = 1000, .N = 5000}, "/stepup/1e3.0");
//    step_up_calc({.L = 1, .T = 0.5, .M = 2000, .N = 10000}, "/stepup/1e3.3");
//    step_up_calc({.L = 1, .T = 0.5, .M = 5000, .N = 25000}, "/stepup/1e3.6");

//    calc_func({.L = 1, .T = 0.5, .M = 10, .N = 50}, "/line_up/1e1.0", line_up);
//    calc_func({.L = 1, .T = 0.5, .M = 20, .N = 100}, "/line_up/1e1.3", line_up);
//    calc_func({.L = 1, .T = 0.5, .M = 50, .N = 250}, "/line_up/1e1.6", line_up);
//    calc_func({.L = 1, .T = 0.5, .M = 100, .N = 500}, "/line_up/1e2.0", line_up);
//    calc_func({.L = 1, .T = 0.5, .M = 200, .N = 1000}, "/line_up/1e2.3", line_up);
//    calc_func({.L = 1, .T = 0.5, .M = 500, .N = 2500}, "/line_up/1e2.6", line_up);
//    calc_func({.L = 1, .T = 0.5, .M = 1000, .N = 5000}, "/line_up/1e3.0", line_up);

    calc_func({.L = 1, .T = 0.5, .M = 2, .N = 3}, "/soliton_up/1e0.3", soliton_up, true);
    calc_func({.L = 1, .T = 0.5, .M = 5, .N = 8}, "/soliton_up/1e0.6", soliton_up, true);
    calc_func({.L = 1, .T = 0.5, .M = 10, .N = 15}, "/soliton_up/1e1.0", soliton_up, true);
    calc_func({.L = 1, .T = 0.5, .M = 20, .N = 30}, "/soliton_up/1e1.3", soliton_up, true);
    calc_func({.L = 1, .T = 0.5, .M = 50, .N = 75}, "/soliton_up/1e1.6", soliton_up, true);
    calc_func({.L = 1, .T = 0.5, .M = 100, .N = 150}, "/soliton_up/1e2.0", soliton_up, true);
    calc_func({.L = 1, .T = 0.5, .M = 200, .N = 300}, "/soliton_up/1e2.3", soliton_up, true);
    calc_func({.L = 1, .T = 0.5, .M = 500, .N = 750}, "/soliton_up/1e2.6", soliton_up, true);
    calc_func({.L = 1, .T = 0.5, .M = 1000, .N = 1500}, "/soliton_up/1e3.0", soliton_up, true);

//    calc_func({.L = 1, .T = 0.5, .M = 2000, .N = 10000}, "/line_up/1e3.3", line_up);
//    calc_func({.L = 1, .T = 0.5, .M = 5000, .N = 25000}, "/line_up/1e3.6", line_up);




//    for (int i = 0; i < M; i++)
//        u0[i] = step(1, 2, 0.5f, h*double(i), 0);
//
//    write2file(u0, progress_KIR, "/stepup/10_3dataKIR.csv");
//    write2file(u0, progress_LW, "/stepup/10_3dataLW.csv");
//    write2file()
//
//    std::cout << "Step-down calculation\n";
//
//    for (int i = 0; i < M; i++)
//        u0[i] = step(2, 1, 0.5f, h*double(i), 0);
//
//    write2file(u0, progress_KIR, "/stepdown/10_3dataKIR.csv");
//    write2file(u0, progress_LW, "/stepdown/10_3dataLW.csv");

    return 0;
}
