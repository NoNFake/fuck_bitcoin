// @author yurii etc: NoNFake
// g++ -O3 generator.cpp -o generator

#include <iostream>
#include <chrono>
#include <cstdint>

using namespace std;


struct Point {
    uint64_t x, y;
    bool found;
};


// faster right to  Left binary method
uint64_t pow_mod(uint64_t base, uint64_t exp, uint64_t mod) {
    uint64_t res = 1;

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



bool is_prime(uint64_t n) {
    if (n < 2) return false;
    if (n == 2 || n == 3) return true;
    if (n % 2 == 0 || n % 3 == 0) return false;
    for (uint64_t i = 5; i * i <= n; i += 6) {
        if (n % i == 0 || n % (i + 2) == 0) return false;
    }
    return true;
}


uint64_t get_curve_order(uint64_t P) {
    uint64_t N = 1;
    for (uint64_t x = 0; x < P; ++x) {
        uint64_t x2 = (x * x) % P;
        uint64_t v = ((x2 * x) % P + 7) % P;
        
        if (v == 0) {
            N += 1;
        } else {
            if (pow_mod(v, (P - 1) / 2, P) == 1) {
                N += 2;
            }
        }
    }
    return N;
}


Point find_generator(uint64_t P) {
    for (uint64_t x = 1; x < P; ++x) {
        uint64_t x2 = (x * x) % P;
        uint64_t v = ((x2 * x) % P + 7) % P;
        
        if (pow_mod(v, (P - 1) / 2, P) == 1) {
            if (P % 4 == 3) {
                uint64_t y = pow_mod(v, (P + 1) / 4, P);
                return {x, y, true};
            }
        }
    }
    return {0, 0, false};
}



int main() {
    cout << "[*] looking for curves with parameters ranging from 1_000_000...\n\n";

    uint64_t START = 1000000;

    for (uint64_t Pt = START; Pt < START + 5000; ++Pt) {
        // Нам нужны только простые , дающие остаток 3 при делении на 4
        if (Pt % 4 != 3 || !is_prime(Pt)) continue;

        cout << "[*] checking  P = " << Pt << "...\n";
        
        // start timer
        auto start_time = chrono::high_resolution_clock::now();

        uint64_t N = get_curve_order(Pt);
        Point G = find_generator(Pt);
        
        //  end timer
        auto end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> calc_time = end_time - start_time;
        
      
        if (G.found) {
            cout <<" [+] Успех! Найдено за " << calc_time.count() << " сек\n";
            cout << "    P = " << Pt << "\n";
            cout << "    N = " << N << "\n";
            cout << "    G = (x: " << G.x << ", y: " << G.y << ")\n\n";
            //break; // Если нужно найти больше вариантов - просто закомментируй этот break
        }
        
        
    }

    return 0;
}