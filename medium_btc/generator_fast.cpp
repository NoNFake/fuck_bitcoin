// @author yurii
// Компиляция: g++ -O3 generator_fast.cpp -o generator_fast -lgmpxx -lgmp
#include <iostream>
#include <chrono>
#include <gmpxx.h>

using namespace std;

struct Point {
    mpz_class x, y;
    bool found;
};


bool is_prime(const mpz_class& n) {
    return mpz_probab_prime_p(n.get_mpz_t(), 25) > 0;
}


Point find_generator(mpz_class P) {
    mpz_class p_minus_1_div_2 = (P - 1) / 2;
    mpz_class p_plus_1_div_4 = (P + 1) / 4;
    mpz_class v, res;

    // Начинаем перебор с P - 1 вниз, чтобы x тоже был гигантским!
    for (mpz_class x = P - 1; x > 0; --x) {
        v = (x * x * x + 7) % P;
        
        mpz_powm(res.get_mpz_t(), v.get_mpz_t(), p_minus_1_div_2.get_mpz_t(), P.get_mpz_t());
        
        if (res == 1) {
            mpz_class y;
            mpz_powm(y.get_mpz_t(), v.get_mpz_t(), p_plus_1_div_4.get_mpz_t(), P.get_mpz_t());
            return {x, y, true};
        }
    }
    return {0, 0, false};
}


int main() {
    cout << "[*] looking for curves with parameters ranging from 1_000_000...\n\n";

    mpz_class START("1000000000000000");

    for (mpz_class Pt = START; Pt < START + 5000; ++Pt) {
        if (Pt % 12 != 11 || !is_prime(Pt)) continue;

        cout << "[*] checking  P = " << Pt << "...\n";
        
        // start timer
        auto start_time = chrono::high_resolution_clock::now();

        mpz_class N = Pt + 1;
        Point G = find_generator(Pt);
        
        //  end timer
        auto end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> calc_time = end_time - start_time;
        
      
        if (G.found) {
            cout <<" [+] luck! found by " << calc_time.count() << " s\n";
            cout << "    P = " << Pt << "\n";
            cout << "    N = " << N << "\n";
            cout << "    G = (x: " << G.x << ", y: " << G.y << ")\n\n";
            //break; // Если нужно найти больше вариантов - просто закомментируй этот break
        }
        
        
    }

    return 0;
}