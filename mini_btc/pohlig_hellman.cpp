// @author yurii
// g++ -O3 pohlig_hellman.cpp -o pohlig_hellman
#include <iostream>
#include <vector>
#include <map>
#include <stdexcept>
#include <iomanip>


using namespace std;

typedef int64_t ll;


ll safe_mod(ll a, ll m) {
    return (a % m + m) % m;
}



ll ext_gcd(ll a, ll b, ll &x, ll &y) {
    if (b == 0) {
        x = 1; y = 0;
        return a;
    }
    ll x1, y1;
    ll d = ext_gcd(b, safe_mod(a, b), x1, y1);
    x = y1;
    y = x1 - y1 * (a / b);
    return d;
}



ll mod_inverse(ll a, ll m) {
    ll x, y;
    ll g = ext_gcd(a, m, x, y);
    if (g != 1) throw runtime_error("[-] There is no reverse element");
    return safe_mod(x, m);
}


ll pow_int(ll base, ll exp) {
    ll res = 1;
    while (exp > 0) {
        if (exp & 1) res *= base;
        base *= base;
        exp >>= 1;
    }
    return res;
}



// ── Арифметика на эллиптической кривой ──
struct Point {
    ll x, y;
    bool is_infinity;

    bool operator==(const Point& other) const {
        if (is_infinity && other.is_infinity) return true;
        if (is_infinity || other.is_infinity) return false;
        return x == other.x && y == other.y;
    }
};


class EllipticCurve {
public:
    ll a, b, p;

    EllipticCurve(ll a_val, ll b_val, ll p_val) : a(a_val), b(b_val), p(p_val) {}

    Point add(Point P, Point Q) {
        if (P.is_infinity) return Q;
        if (Q.is_infinity) return P;

        ll m = 0;
        if (P.x == Q.x) {
            if (P.y != Q.y || P.y == 0) return {0, 0, true}; // P + (-P) = O
            
            // Удвоение точки
            ll num = safe_mod(3 * safe_mod(P.x * P.x, p) + a, p);
            ll den = mod_inverse(safe_mod(2 * P.y, p), p);
            m = safe_mod(num * den, p);
        } else {
            // Сложение разных точек
            ll num = safe_mod(Q.y - P.y, p);
            ll den = mod_inverse(safe_mod(Q.x - P.x, p), p);
            m = safe_mod(num * den, p);
        }

        ll rx = safe_mod(safe_mod(m * m, p) - P.x - Q.x, p);
        ll ry = safe_mod(m * safe_mod(P.x - rx, p) - P.y, p);
        return {rx, ry, false};
    }

    Point mul(ll k, Point P) {
        Point result = {0, 0, true};
        Point addend = P;
        while (k > 0) {
            if (k & 1) result = add(result, addend);
            addend = add(addend, addend);
            k >>= 1;
        }
        return result;
    }

    Point neg(Point P) {
        if (P.is_infinity) return P;
        return {P.x, safe_mod(-P.y, p), false};
    }

    ll order_of(Point P) {
        Point cur = P;
        ll n = 1;
        while (!cur.is_infinity) {
            cur = add(cur, P);
            n++;
        }
        return n;
    }
};



map<ll, ll> factorize(ll n) {
    map<ll, ll> factors;
    ll d = 2;
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



ll crt(const vector<ll>& residues, const vector<ll>& moduli) {
    ll M = 1;
    for (ll m : moduli) M *= m;
    ll x = 0;
    
    for (size_t i = 0; i < residues.size(); ++i) {
        ll Mi = M / moduli[i];
        ll yi = mod_inverse(Mi, moduli[i]);
        x = safe_mod(x + residues[i] * safe_mod(Mi * yi, M), M);
    }
    return x;
}

// Pohlig-Hellman ─
ll pohlig_hellman(EllipticCurve& curve, Point G, Point Q, ll group_order) {
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

    vector<ll> residues;
    vector<ll> moduli;

    for (auto const& [q, e] : factors) {
        ll qe = pow_int(q, e);
        cout << "\n  subgroup q=" << q << ", e=" << e << ", q^e=" << qe << "\n";

        Point gamma = curve.mul(group_order / q, G);
        ll k_i = 0;

        for (ll j = 0; j < e; ++j) {
            ll coeff = group_order / pow_int(q, j + 1);
            Point term = curve.mul(k_i, G);
            Point diff = curve.add(Q, curve.neg(term));
            Point delta = curve.mul(coeff, diff);

            ll d_j = -1;
            Point test = {0, 0, true};
            
            for (ll d = 0; d < q; ++d) {
                if (test == delta) {
                    d_j = d;
                    break;
                }
                test = curve.add(test, gamma);
            }

            if (d_j == -1) {
                throw runtime_error("d_j not found! Perhaps point Q does not lie on this subgroup.");
            }

            k_i = k_i + d_j * pow_int(q, j);
            cout << "    j=" << j << ": d_j=" << d_j << " -> k mod " << pow_int(q, j + 1) << " = " << k_i << "\n";
        }

        residues.push_back(k_i % qe);
        moduli.push_back(qe);
        cout << "  -> k ≡ " << (k_i % qe) << " (mod " << qe << ")\n";
    }

    return crt(residues, moduli);
}


int main() {
    cout << "============================================================\n";
    cout << "  EC Discrete Log — Pohlig-Hellman + CRT\n";
    cout << "============================================================\n";

    ll p = 1000931;
    ll a = 0;
    ll b = 7;
    EllipticCurve curve(a, b, p);

    Point G = {2, 727727, false};
    Point Q = {822105, 694873, false}; 
    // Внимание: если для Q=(544, 709) не существует такого k, 
    // что Q=k*G на этой кривой, алгоритм выбросит ошибку "Не найдено d_j".





    cout << "\ncurve : y^2 = x^3 + " << a << "x + " << b << "  (mod " << p << ")\n";
    cout << "G      : (" << G.x << ", " << G.y << ")\n";
    cout << "Q = k*G: (" << Q.x << ", " << Q.y << ")\n";

    cout << "\n[1] calculate the group order #E(F_" << p << ")...\n";
    ll group_order = curve.order_of(G);
    cout << "    #E(F_" << p << ") = " << group_order << "\n";

    cout << "\n[2] Pohlig-Hellman...\n";
    ll k = 0;
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

    cout << "\ni*G table (the first 100 points are shown to save space):\n";
    Point pt = {0, 0, true};
    
    // Ограничиваем вывод 100 строками, чтобы не зависло
    ll limit = min(group_order, (ll)200000); 
    
    for (ll i = 1; i <= limit; ++i) {
        pt = pt.is_infinity ? G : curve.add(pt, G);
        
        // string marker = (i == k) ? " <- k (searching)" : "";
        string marker =  " <- k (searching)" ;
        if (pt == Q) marker += " <- Q";
        
        if (pt.is_infinity) {
            cout << "  i=" << setw(3) << i << ": (Infinity)          " << marker << "\n";
        } else {
            if ( i == k ) {
                cout << "  i=" << setw(3) << i << ": (" << pt.x << ", " << pt.y << ")  " << marker << "\n";
\
            }
        }
    }
    if (group_order > 100) cout << "  ... and " << group_order - 100 << " points ...\n";

    return 0;
}