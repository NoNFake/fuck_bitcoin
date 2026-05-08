// @author yurii etc: NoNFake
// g++ -O3 generator.cpp -o generator -lgmpxx -lgmp

#include <iostream>
#include <chrono>
#include <cstdint>
#include <gmpxx.h>

using namespace std;


struct Point {
    mpz_class x, y;
    bool found;
};


// faster right to  Left binary method
mpz_class pow_mod(mpz_class base, mpz_class exp, mpz_class mod) {
    mpz_class res = 1;

    base  = base % mod;

    while (exp > 0) {
        if (exp % 2 == 1 ) {
            res = (res * base) % mod;
        }

        base = (base * base) % mod;

        exp /= 2;
    }

    return res;
}



bool is_prime(const mpz_class& n) {
    return mpz_probab_prime_p(n.get_mpz_t(), 25) > 0;
}


mpz_class get_curve_order(mpz_class P) {
    mpz_class N = 1; // Учитываем точку в бесконечности
    mpz_class p_minus_1_div_2 = (P - 1) / 2;
    mpz_class v, res;

    // ВНИМАНИЕ: Этот цикл делает P шагов. 
    // Для P ~ 100_000_000 отработает быстро, для 256-битных чисел зависнет!
    for (mpz_class x = 0; x < P; ++x) {
        v = (x * x * x + 7) % P;
        
        if (v == 0) {
            N += 1;
        } else {
            mpz_powm(res.get_mpz_t(), v.get_mpz_t(), p_minus_1_div_2.get_mpz_t(), P.get_mpz_t());
            if (res == 1) {
                N += 2;
            }
        }
    }
    return N;
}


Point find_generator(mpz_class P) {
    mpz_class p_minus_1_div_2 = (P - 1) / 2;
    mpz_class p_plus_1_div_4 = (P + 1) / 4;
    mpz_class v, res;

    for (mpz_class x = P / 3; x < P; ++x) {
        v = (x * x * x + 7) % P;
        
        // res = v^((P-1)/2) % P
        mpz_powm(res.get_mpz_t(), v.get_mpz_t(), p_minus_1_div_2.get_mpz_t(), P.get_mpz_t());
        
        if (res == 1) {
            if (P % 4 == 3) {
                mpz_class y;
                // y = v^((P+1)/4) % P
                mpz_powm(y.get_mpz_t(), v.get_mpz_t(), p_plus_1_div_4.get_mpz_t(), P.get_mpz_t());
                return {x, y, true};
            }
        }
    }
    return {0, 0, false};
}



int main() {
    cout << "[*] looking for curves with parameters ranging from 1_000_000...\n\n";

    mpz_class START("1000000000000000");

    for (mpz_class Pt = START; Pt < START + 5000; ++Pt) {
        // Нам нужны только простые , дающие остаток 3 при делении на 4
        if (Pt % 4 != 3 || !is_prime(Pt)) continue;

        cout << "[*] checking  P = " << Pt << "...\n";
        
        // start timer
        auto start_time = chrono::high_resolution_clock::now();

        mpz_class N = get_curve_order(Pt);
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