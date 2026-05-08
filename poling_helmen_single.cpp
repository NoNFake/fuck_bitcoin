//  g++ -O3 -march=native test.cpp -o pohlig_hellman -lgmpxx -lgmp
#include <iostream>
#include <stdexcept>
#include <gmp.h>

using namespace std;
mpz_t curve_a, curve_b, curve_p;

void init_curve(const char* a, const char* b, const char* p) {
    mpz_init_set_str(curve_a, a, 10);
    mpz_init_set_str(curve_b, b, 10);
    mpz_init_set_str(curve_p, p, 10);
}



// ── Арифметика на эллиптической кривой ──
struct Point {
    mpz_t x, y;
    bool is_infinity;


    Point() {
        mpz_init(x);
        mpz_init(y);

        is_infinity = true;
    }

    ~Point() {
        mpz_clear(x);
        mpz_clear(y);
    }


    void set(const Point& other) {
        mpz_set(x, other.x);
        mpz_set(y, other.y);
        is_infinity = other.is_infinity;
    }



    void set_str(const char* x_str, const char* y_str) {
        mpz_set_str(x, x_str, 10);
        mpz_set_str(y, y_str, 10);
        is_infinity = false;
    }


    bool equals(const Point& other) const {
        if (is_infinity && other.is_infinity) return true;
        if (is_infinity || other.is_infinity) return false;
        return mpz_cmp(x, other.x) == 0 && mpz_cmp(y, other.y) == 0;
    }

};






mpz_t tmp_num, tmp_den, tmp_m, tmp_rx, tmp_ry, tmp_diff;


void init_math(){
    mpz_inits(tmp_num, tmp_den, tmp_m, tmp_rx, tmp_ry, tmp_diff, NULL);
}



void add_points(Point& res, const Point& P, const Point& Q) {
    if (P.is_infinity) { res.set(Q); return; }
    if (Q.is_infinity) { res.set(P); return; }


    if (mpz_cmp(P.x, Q.x) == 0) {
        if (mpz_cmp(P.y, Q.y) != 0 || mpz_cmp_ui(P.y, 0) == 0) {
            res.is_infinity = true;
            return;
        }


        mpz_mul(tmp_num, P.x, P.x);
        mpz_mod(tmp_num, tmp_num, curve_p);
        mpz_mul_ui(tmp_num, tmp_num, 3);
        mpz_add(tmp_num, tmp_num, curve_a);
        mpz_mod(tmp_num, tmp_num, curve_p);

        mpz_mul_ui(tmp_den, P.y, 2);
        mpz_mod(tmp_den, tmp_den, curve_p);

    } else {
        mpz_sub(tmp_num, Q.y, P.y);
        mpz_mod(tmp_num, tmp_num, curve_p);
        
        mpz_sub(tmp_den, Q.x, P.x);
        mpz_mod(tmp_den, tmp_den, curve_p);
    }




    if (mpz_invert(tmp_den, tmp_den, curve_p) == 0) {
        throw runtime_error("no invrese");
    }

    mpz_mul(tmp_m, tmp_num, tmp_den);
    mpz_mod(tmp_m, tmp_m, curve_p);


    // rx = m^2 - P.x - Q.x
    mpz_mul(tmp_rx, tmp_m, tmp_m);
    mpz_sub(tmp_rx, tmp_rx, P.x);
    mpz_sub(tmp_rx, tmp_rx, Q.x);
    mpz_mod(tmp_rx, tmp_rx, curve_p);
    if (mpz_sgn(tmp_rx) < 0) mpz_add(tmp_rx, tmp_rx, curve_p);


    // ry = m*(P.x - rx) - P.y
    mpz_sub(tmp_diff, P.x, tmp_rx);
    mpz_mul(tmp_ry, tmp_m, tmp_diff);
    mpz_sub(tmp_ry, tmp_ry, P.y);
    mpz_mod(tmp_ry, tmp_ry, curve_p);
    if (mpz_sgn(tmp_ry) < 0) mpz_add(tmp_ry, tmp_ry, curve_p);

    mpz_set(res.x, tmp_rx);
    mpz_set(res.y, tmp_ry);
    res.is_infinity = false;

}



// optimize pollard`s rho
void pollards_rho_fast(mpz_t result, const Point& G, const Point& Q, const mpz_t N) {
    if (Q.is_infinity) {
        mpz_set_ui(result, 0);
        return;
    }

    Point T, H;
    T.set(G);
    H.set(G);

    mpz_t aT, bT, aH, bH, subset;
    mpz_inits(aT, bT, aH, bH, subset, NULL);
    mpz_set_ui(aT, 1); mpz_set_ui(bT, 0);
    mpz_set_ui(aH, 1); mpz_set_ui(bH, 0);


    auto next_step = [&](Point& P, mpz_t& a, mpz_t& b) {
        if (P.is_infinity) return;
        mpz_fdiv_r_ui(subset, P.x, 3); 
        unsigned long s = mpz_get_ui(subset);

        if (s == 0) {
            mpz_add_ui(a, a, 1);
            mpz_mod(a, a, N);
            add_points(P, P, G);
        } else if (s == 1) {
            mpz_mul_ui(a, a, 2); mpz_mod(a, a, N);
            mpz_mul_ui(b, b, 2); mpz_mod(b, b, N);
            add_points(P, P, P);
        } else {
            mpz_add_ui(b, b, 1);
            mpz_mod(b, b, N);
            add_points(P, P, Q);
        }
    };


    cout << "    [!] started optimize pollards rho... \n";


    long long iterations = 0;
    while (true) {
        // Черепаха делает 1 шаг
        next_step(T, aT, bT);
        
        // Заяц делает 2 шага
        next_step(H, aH, bH);
        next_step(H, aH, bH);
        
        iterations++;
        if (iterations % 500000 == 0) { 
            cout << "      ... done " << iterations << " interations...\r" << flush;
        }

        if (T.equals(H)) {
            if (T.is_infinity) throw runtime_error("[-] in inf.");
            break;
        }
    }


    cout << "\n      [+] colision found by " << iterations << " interations!\n";

    mpz_t r, r_inv, num;
    mpz_inits(r, r_inv, num, NULL);

    mpz_sub(r, bT, bH);
    mpz_mod(r, r, N);
    if (mpz_sgn(r) < 0) mpz_add(r, r, N);

    if (mpz_cmp_ui(r, 0) == 0) throw runtime_error("[-] Division by zero at the end.");

    mpz_invert(r_inv, r, N);

    // num = (aH - aT) mod N
    mpz_sub(num, aH, aT);
    mpz_mod(num, num, N);
    if (mpz_sgn(num) < 0) mpz_add(num, num, N);

    // result = (num * r_inv) mod N
    mpz_mul(result, num, r_inv);
    mpz_mod(result, result, N);

    mpz_clears(aT, bT, aH, bH, subset, r, r_inv, num, NULL);

}



int main() {
    cout << "============================================================\n";
    cout << "  optimize hacjk (clear GMP C)\n";
    cout << "============================================================\n\n";


    init_curve(
        "0", 
        "7", 
        "115792089237316195423570985008687907853269984665640564039457584007908834671663"
    );

    init_math();

    Point G, Q;
    G.set_str(
        "55066263022277343669578718895168534326250603453777594175500187360389116729240",
        "32670510020758816978083085130507043184471273380659243275938904335757337482424"
    );



    Q.set_str(
        "71715483259960251104886616086749964585141029199321453995874288114660996885376",
        "40560706241088654877709886733221375330368142163939632832049960868661138127382"
    );


    mpz_t group_order;
    mpz_init_set_str(group_order, "115792089237316195423570985008687907852837564279074904382605163141518161494337", 10);

    mpz_t secret_key;
    mpz_init(secret_key);



    try {
        pollards_rho_fast(secret_key, G, Q, group_order);
        gmp_printf("\n[+] secret key found:\n    k = %Zd\n", secret_key);
    } catch(const exception& e) {
        cerr << "\n[!] error: " << e.what() << "\n";
    }


    mpz_clear(group_order);
    mpz_clear(secret_key);
    return 0;

}   