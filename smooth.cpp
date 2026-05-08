// @author yurii
// Компиляция: g++ -O3 smooth.cpp -o smooth -lgmpxx -lgmp
#include <iostream>
#include <vector>
#include <gmpxx.h>
#include <chrono>
#include <stdexcept>

using namespace std;

// ── Вспомогательная математика ──
mpz_class safe_mod(mpz_class a, mpz_class m) {
    mpz_class res = a % m;
    if (res < 0) res += m;
    return res;
}

mpz_class mod_inverse(mpz_class a, mpz_class m) {
    mpz_class res;
    if (mpz_invert(res.get_mpz_t(), a.get_mpz_t(), m.get_mpz_t()) == 0) {
        throw runtime_error("[-] Нет обратного элемента");
    }
    return res;
}

struct Point {
    mpz_class x, y;
    bool is_infinity; 
};

// ── Арифметика кривой ──
class EllipticCurve {
public:
    mpz_class a, b, p;

    EllipticCurve(mpz_class a_val, mpz_class b_val, mpz_class p_val) : a(a_val), b(b_val), p(p_val) {}

    Point add(Point P, Point Q) {
        if (P.is_infinity) return Q;
        if (Q.is_infinity) return P;

        mpz_class m = 0;
        if (P.x == Q.x) {
            if (P.y != Q.y || P.y == 0) return {0, 0, true};
            
            mpz_class num = safe_mod(3 * safe_mod(P.x * P.x, p) + a, p);
            mpz_class den = mod_inverse(safe_mod(2 * P.y, p), p);
            m = safe_mod(num * den, p);
        } else {
            mpz_class num = safe_mod(Q.y - P.y, p);
            mpz_class den = mod_inverse(safe_mod(Q.x - P.x, p), p);
            m = safe_mod(num * den, p);
        }

        mpz_class rx = safe_mod(safe_mod(m * m, p) - P.x - Q.x, p);
        mpz_class ry = safe_mod(m * safe_mod(P.x - rx, p) - P.y, p);
        return {rx, ry, false};
    }

    Point mul(mpz_class k, Point P) {
        Point result = {0, 0, true};
        Point addend = P;
        while (k > 0) {
            if ((k % 2) != 0) result = add(result, addend);
            addend = add(addend, addend);
            k /= 2; 
        }
        return result;
    }
};

bool is_prime(const mpz_class& n) {
    return mpz_probab_prime_p(n.get_mpz_t(), 25) > 0;
}

// ── НОВАЯ ФУНКЦИЯ: Поиск ИДЕАЛЬНОЙ базовой точки ──
Point find_perfect_generator(mpz_class P, mpz_class N) {
    EllipticCurve curve(0, 7, P);
    vector<int> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29};
    
    mpz_class p_minus_1_div_2 = (P - 1) / 2;
    mpz_class p_plus_1_div_4 = (P + 1) / 4;
    mpz_class v, res;

    // Идем с конца, проверяем точки
    for (mpz_class x = P - 1; x > 0; --x) {
        v = (x * x * x + 7) % P;
        mpz_powm(res.get_mpz_t(), v.get_mpz_t(), p_minus_1_div_2.get_mpz_t(), P.get_mpz_t());
        
        if (res == 1) { // Если y существует
            mpz_class y;
            mpz_powm(y.get_mpz_t(), v.get_mpz_t(), p_plus_1_div_4.get_mpz_t(), P.get_mpz_t());
            Point G = {x, y, false};
            
            // ПРОВЕРКА КАЧЕСТВА: Точка должна генерировать все подгруппы!
            bool perfect = true;
            for (int q : primes) {
                if (curve.mul(N / q, G).is_infinity) {
                    perfect = false;
                    break; 
                }
            }
            
            // Если точка идеальна - возвращаем её
            if (perfect) return G;
        }
    }
    return {0, 0, true};
}

int main() {
    cout << "=================================================\n";
    cout << "  ГЕНЕРАТОР ГИГАНТСКИХ ГЛАДКИХ КРИВОЙ (~256 бит)\n";
    cout << "=================================================\n\n";

    mpz_class N_base = 1;
    mpz_class temp;

    mpz_pow_ui(temp.get_mpz_t(), mpz_class(2).get_mpz_t(), 60); N_base *= temp;
    mpz_pow_ui(temp.get_mpz_t(), mpz_class(3).get_mpz_t(), 40); N_base *= temp;
    mpz_pow_ui(temp.get_mpz_t(), mpz_class(5).get_mpz_t(), 20); N_base *= temp;
    mpz_pow_ui(temp.get_mpz_t(), mpz_class(7).get_mpz_t(), 15); N_base *= temp;
    mpz_pow_ui(temp.get_mpz_t(), mpz_class(11).get_mpz_t(), 10); N_base *= temp;
    mpz_pow_ui(temp.get_mpz_t(), mpz_class(13).get_mpz_t(), 10); N_base *= temp;
    mpz_pow_ui(temp.get_mpz_t(), mpz_class(17).get_mpz_t(), 5); N_base *= temp;
    mpz_pow_ui(temp.get_mpz_t(), mpz_class(19).get_mpz_t(), 5); N_base *= temp;
    mpz_pow_ui(temp.get_mpz_t(), mpz_class(23).get_mpz_t(), 5); N_base *= temp;
    mpz_pow_ui(temp.get_mpz_t(), mpz_class(29).get_mpz_t(), 5); N_base *= temp;

    cout << "[*] Базовое супер-гладкое число N_base сформировано.\n";
    cout << "[*] Начинаем перебор множителя k для поиска простого P...\n\n";

    auto start_time = chrono::high_resolution_clock::now();

    for (unsigned long k = 1; k < 1000000; ++k) {
        mpz_class N = N_base * k;
        mpz_class P = N - 1;

        if ((P % 12 == 11) && is_prime(P)) {
            // Используем новую функцию поиска идеальной точки
            Point G = find_perfect_generator(P, N);

            if (!G.is_infinity) {
                auto end_time = chrono::high_resolution_clock::now();
                chrono::duration<double> calc_time = end_time - start_time;
                
                EllipticCurve curve(0, 7, P);
                mpz_class secret_k("313370000000000000000000000000000000000000000000000000001337");
                Point Q = curve.mul(secret_k, G);

                cout << "[+] ИДЕАЛЬНАЯ МИШЕНЬ НАЙДЕНА ЗА " << calc_time.count() << " СЕК!\n";
                cout << "--------------------------------------------------------\n";
                cout << "mpz_class p(\"" << P << "\");\n";
                cout << "mpz_class group_order(\"" << N << "\");\n";
                cout << "Point G = {mpz_class(\"" << G.x << "\"), mpz_class(\"" << G.y << "\"), false};\n";
                cout << "Point Q = {mpz_class(\"" << Q.x << "\"), mpz_class(\"" << Q.y << "\"), false};\n";
                cout << "--------------------------------------------------------\n";
                break;
            }
        }
    }

    return 0;
}