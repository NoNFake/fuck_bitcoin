// @author yurii
// g++ -O3 pohlig_hellman.cpp -o pohlig_hellman -lgmpxx -lgmp
#include <iostream>
#include <iostream>
#include <vector>
#include <map>
#include <stdexcept>
#include <iomanip>
#include <gmpxx.h>

using namespace std;



mpz_class safe_mod(mpz_class a, mpz_class m) {
    mpz_class res = a % m;
    if (res < 0) res += m;
    return res;
}




mpz_class mod_inverse(mpz_class a, mpz_class m) {
    mpz_class res;
    if (mpz_invert(res.get_mpz_t(), a.get_mpz_t(), m.get_mpz_t()) == 0) {
        throw runtime_error("[-] There is no reverse element");
    }
    return res;
}

mpz_class pow_int(mpz_class base, unsigned long exp) {
    mpz_class res;
    mpz_pow_ui(res.get_mpz_t(), base.get_mpz_t(), exp);
    return res;
}



// ── Арифметика на эллиптической кривой ──
struct Point {
    mpz_class x, y;
    bool is_infinity;

    bool operator==(const Point& other) const {
        if (is_infinity && other.is_infinity) return true;
        if (is_infinity || other.is_infinity) return false;
        return x == other.x && y == other.y;
    }
};


class EllipticCurve {
public:
    mpz_class a, b, p;

    EllipticCurve(mpz_class a_val, mpz_class b_val, mpz_class p_val) : a(a_val), b(b_val), p(p_val) {}

    Point add(Point P, Point Q) {
        if (P.is_infinity) return Q;
        if (Q.is_infinity) return P;

        mpz_class m = 0;
        if (P.x == Q.x) {
            if (P.y != Q.y || P.y == 0) return {0, 0, true}; // P + (-P) = O
            
            // Удвоение точки
            mpz_class num = safe_mod(3 * safe_mod(P.x * P.x, p) + a, p);
            mpz_class den = mod_inverse(safe_mod(2 * P.y, p), p);
            m = safe_mod(num * den, p);
        } else {
            // Сложение разных точек
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
            // Замена битовых операций на совместимые с GMP
            if ((k % 2) != 0) result = add(result, addend);
            addend = add(addend, addend);
            k /= 2; 
        }
        return result;
    }

    Point neg(Point P) {
        if (P.is_infinity) return P;
        return {P.x, safe_mod(-P.y, p), false};
    }

    mpz_class order_of(Point P) {
        Point cur = P;
        mpz_class n = 1;
        while (!cur.is_infinity) {
            cur = add(cur, P);
            n++;
        }
        return n;
    }
};



map<mpz_class, unsigned long> factorize(mpz_class n) {
    map<mpz_class, unsigned long> factors;
    mpz_class d = 2;
    while (d * d <= n) {
        while (n % d == 0) {
            factors[d]++;
            n /= d;
        }
        d++;
    }
    if (n > 1) factors[n]++;
    return factors;
}



mpz_class crt(const vector<mpz_class>& residues, const vector<mpz_class>& moduli) {
    mpz_class M = 1;
    for (const mpz_class& m : moduli) M *= m;
    mpz_class x = 0;
    
    for (size_t i = 0; i < residues.size(); ++i) {
        mpz_class Mi = M / moduli[i];
        mpz_class yi = mod_inverse(Mi, moduli[i]);
        x = safe_mod(x + residues[i] * safe_mod(Mi * yi, M), M);
    }
    return x;
}



// Pohlig-Hempz_classman ─
mpz_class pohlig_hellman(EllipticCurve& curve, Point G, Point Q, mpz_class group_order) {
    auto factors = factorize(group_order);
    
    cout << "\n  Order factorization " << group_order << ": ";
    bool first = true;
    for (auto const& [q, e] : factors) {
        if (!first) cout << " x ";
        if (e > 1) cout << q << "^" << e;
        else cout << q;
        first = false;
    }
    cout << "\n";

    vector<mpz_class> residues;
    vector<mpz_class> moduli;

    for (auto const& [q, e] : factors) {
        mpz_class qe = pow_int(q, e);
        cout << "\n  subgroup q=" << q << ", e=" << e << ", q^e=" << qe << "\n";

        mpz_class group_order_div_q = group_order / q;
        Point gamma = curve.mul(group_order_div_q, G);
        mpz_class k_i = 0;

        for (unsigned long j = 0; j < e; ++j) {
            mpz_class q_pow_j_plus_1 = pow_int(q, j + 1);
            mpz_class coeff = group_order / q_pow_j_plus_1;
            
            Point term = curve.mul(k_i, G);
            Point diff = curve.add(Q, curve.neg(term));
            Point delta = curve.mul(coeff, diff);

            mpz_class d_j = -1;
            Point test = {0, 0, true};
            
            for (mpz_class d = 0; d < q; ++d) {
                if (test == delta) {
                    d_j = d;
                    break;
                }
                test = curve.add(test, gamma);
            }

            if (d_j == -1) {
                throw runtime_error("d_j not found! Perhaps point Q does not lie on this subgroup.");
            }

            mpz_class q_pow_j = pow_int(q, j);
            k_i = k_i + d_j * q_pow_j;
            cout << "    j=" << j << ": d_j=" << d_j << " -> k mod " << q_pow_j_plus_1 << " = " << k_i << "\n";
        }

        residues.push_back(k_i % qe);
        moduli.push_back(qe);
        cout << "  -> k ≡ " << (k_i % qe) << " (mod " << qe << ")\n";
    }

    return crt(residues, moduli);
}


int main() {
    cout << "============================================================\n";
    cout << "  EC Discrete Log — Pohlig-Hellman + CRT (GMP EDITION)\n";
    cout << "============================================================\n";

    // Инициализация через строки или обычные числа, mpz_class скушает всё!
    mpz_class p("1000000000004783");
    mpz_class a("0");
    mpz_class b("7");
    EllipticCurve curve(a, b, p);

    Point G = {mpz_class("1000000000004782"), mpz_class("25166234176929"), false};
    Point Q = {mpz_class("725828768659486"), mpz_class("364722906798421"), false}; 
    mpz_class group_order("1000000000001784");

    cout << "\ncurve : y^2 = x^3 + " << a << "x + " << b << "  (mod " << p << ")\n";
    cout << "G      : (" << G.x << ", " << G.y << ")\n";
    cout << "Q = k*G: (" << Q.x << ", " << Q.y << ")\n";

    cout << "\n[1] calculate the group order #E(F_" << p << ")...\n";
    cout << "    #E(F_" << p << ") = " << group_order << "\n";

    cout << "\n[2] Pohlig-Hellman...\n";
    mpz_class k = 0;
    try {
        k = pohlig_hellman(curve, G, Q, group_order);
    } catch(const exception& e) {
        cerr << "\n[!] error: " << e.what() << "\n";
        return 1; 
    }

    cout << "\n[3] checking...\n";
    Point check = curve.mul(k, G);
    bool ok = (check == Q);
    cout << "    k = " << k << "\n";
    if (check.is_infinity) cout << "    k*G = (point in inf)\n";
    else cout << "    k*G = (" << check.x << ", " << check.y << ")\n";
    cout << "    Matches с Q: " << (ok ? "✓ yes" : "✗ no") << "\n";

    cout << "\n============================================================\n";
    cout << "  answer: k = " << k << "\n";
    cout << "============================================================\n";

    cout << "\ni*G table (showing relevant points):\n";
    Point pt = {0, 0, true};
    
    // Ограничиваем вывод, чтобы консоль не умерла от огромных N
    mpz_class limit;
    if (group_order > 200000) limit = 200000;
    else limit = group_order;
    
    for (mpz_class i = 1; i <= limit; ++i) {
        pt = pt.is_infinity ? G : curve.add(pt, G);
        
        string marker = " <- k (searching)";
        if (pt == Q) marker += " <- Q";
        
        if (pt.is_infinity) {
            // Выводим только если это бесконечность
            cout << "  i=" << setw(3) << i.get_str() << ": (Infinity)          " << marker << "\n";
        } else {
            // Выводим только ту самую искомую точку, чтобы не спамить
            if (i == k) {
                cout << "  i=" << setw(3) << i.get_str() << ": (" << pt.x << ", " << pt.y << ")  " << marker << "\n";
            }
        }
    }
    
    if (group_order > limit) cout << "  ... and " << (group_order - limit) << " points ...\n";

    return 0;
}