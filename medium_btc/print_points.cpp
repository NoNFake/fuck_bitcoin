// @author yurii etc: NoNFake
// g++ -O3 print_points_gmp.cpp -o print_points_gmp -lgmpxx -lgmp

#include <iostream>
#include <cstdint>
#include <iomanip>
#include <stdexcept>
#include <gmpxx.h>

using namespace std;

struct Point {
    mpz_class x, y;
    bool is_infinity;
};


mpz_class sub_mod(mpz_class a, mpz_class b, mpz_class p) {
    return (a + p - (b % p)) % p;
}


mpz_class safe_mod(mpz_class a, mpz_class p) {
    mpz_class res = a % p;
    if (res < 0) res += p;
    return res;
}


mpz_class pow_mod(mpz_class base, mpz_class exp, mpz_class mod) {
    mpz_class res;
    mpz_powm(res.get_mpz_t(), base.get_mpz_t(), exp.get_mpz_t(), mod.get_mpz_t());
    return res;
}


mpz_class mod_inverse(mpz_class n, mpz_class p) {
    mpz_class res;
    if (mpz_invert(res.get_mpz_t(), n.get_mpz_t(), p.get_mpz_t()) == 0) {
        throw runtime_error("There is no inverse element (division by zero)!");
    }
    return res;
}

Point add_points(Point P1, Point P2, mpz_class p) {
    if (P1.is_infinity) return P2;
    if (P2.is_infinity) return P1;

    mpz_class s = 0; 

    if (P1.x == P2.x && P1.y == P2.y) {
        if (P1.y == 0) return {0, 0, true}; 
        
        mpz_class num = safe_mod(3 * P1.x * P1.x, p); // a = 0 для Биткоина, так что +a не пишем
        mpz_class den = mod_inverse(safe_mod(2 * P1.y, p), p);
        s = safe_mod(num * den, p);
    } else {
        if (P1.x == P2.x) return {0, 0, true}; 
        
        mpz_class num = safe_mod(P2.y - P1.y, p);
        mpz_class den = mod_inverse(safe_mod(P2.x - P1.x, p), p);
        s = safe_mod(num * den, p);
    }

    mpz_class x3 = safe_mod(s * s - P1.x - P2.x, p);
    mpz_class y3 = safe_mod(s * (P1.x - x3) - P1.y, p);

    return {x3, y3, false};
}

int main() {
    // Берем "боевую" кривую из твоего лога (не суперсингулярную!)
    // mpz_class P("115792089237316195423570985008687907853269984665640564039457584007908834671663");
    // N 1000932
    // mpz_class Gx("55066263022277343669578718895168534326250603453777594175500187360389116729240");
    // mpz_class Gy("32670510020758816978083085130507043184471273380659243275938904335757337482424");

    mpz_class P("1000000000004783");
    mpz_class Gx("1000000000004782");
    mpz_class Gy("25166234176929");
    
    
    Point G = {Gx, Gy, false};


    int max_i = 100200;
    Point current = G; // Начинаем с 1*G

    cout << "-------------------------------------------\n";
    cout << " Bitcoin Curve secp256k1 y^2 = x^3 + 7 \n";
    cout << " module P = \n" << P << "\n";
    cout << "-------------------------------------------\n\n";

    for (int i = 100000; i <= max_i; ++i) {
        if (current.is_infinity) {
            cout << "i:" << i << " :  point in inf (O)\n";
        } else {
            cout << "i:" << i << " :  x:" << current.x << "  y:" << current.y << "\n";
        }
        
        // Магия оптимизации: для вычисления (i+1)*G мы просто прибавляем G к текущей точке
        current = add_points(current, G, P);
    }

    return 0;
}