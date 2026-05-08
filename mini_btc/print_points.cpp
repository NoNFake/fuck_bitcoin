// @author yurii etc: NoNFake
// g++ -O3 print_points.cpp -o print_points

#include <iostream>
#include <cstdint>
#include <iomanip>


using namespace std;

struct Point {
    uint64_t x, y;
    bool is_infinity;
};


uint64_t sub_mod(uint64_t a, uint64_t b, uint64_t p) {
    return (a + p - (b % p)) % p;
}


uint64_t pow_mod(uint64_t base, uint64_t exp, uint64_t mod) {
    uint64_t res = 1;
    base = base % mod;
    while (exp > 0) {
        if (exp % 2 == 1) res = (res * base) % mod;
        base = (base * base) % mod;
        exp /= 2;
    }
    return res;
}


uint64_t mod_inverse(uint64_t n, uint64_t p) {
    return pow_mod(n, p - 2, p);
}


Point add_points(Point P1, Point P2, uint64_t p) {
    // Если одна из точек - бесконечность, возвращаем другую
    if (P1.is_infinity) return P2;
    if (P2.is_infinity) return P1;

    uint64_t s = 0; // Наклон (slope)

    if (P1.x == P2.x && P1.y == P2.y) {
        // Случай 1: Удвоение точки (P1 + P1)
        if (P1.y == 0) return {0, 0, true}; // Касательная вертикальна
        
        // s = (3 * x^2) / (2 * y) mod p
        uint64_t num = (3 * P1.x % p * P1.x) % p;
        uint64_t den = mod_inverse((2 * P1.y) % p, p);
        s = (num * den) % p;
    } else {
        // Случай 2: Сложение разных точек
        if (P1.x == P2.x) return {0, 0, true}; // Точки лежат на одной вертикали
        
        // s = (y2 - y1) / (x2 - x1) mod p
        uint64_t num = sub_mod(P2.y, P1.y, p);
        uint64_t den = mod_inverse(sub_mod(P2.x, P1.x, p), p);
        s = (num * den) % p;
    }

    // Считаем новые координаты
    // x3 = s^2 - x1 - x2
    uint64_t x3 = sub_mod(sub_mod((s * s) % p, P1.x, p), P2.x, p);
    // y3 = s * (x1 - x3) - y1
    uint64_t y3 = sub_mod((s * sub_mod(P1.x, x3, p)) % p, P1.y, p);

    return {x3, y3, false};
}

int main() {
    // Берем "боевую" кривую из твоего лога (не суперсингулярную!)
    uint64_t P = 1000931; 
    // N 1000932
    Point G = {2, 727727, false}; // Это 1*G

    int max_i = 100200;
    Point current = G; // Начинаем с 1*G

    cout << "-------------------------------------------\n";
    cout << " curve y^2 = x^3 + 7 (mod " << P << ")\n";
    cout << " base point G = (" << G.x << ", " << G.y << ")\n";
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