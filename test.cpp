//  g++ -O3 -march=native -funroll-loops -ffast-math \
    test.cpp -o pohlig_hellman \
    -lgmpxx -lgmp -lpthread -ltbb


#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <mutex>
#include <shared_mutex>  
#include <unordered_map>
#include <random>
#include <atomic>
#include <chrono>
#include <gmp.h>



using namespace std;
mpz_t curve_a, curve_b, curve_p, group_order;

void init_curve(const char* a, const char* b, const char* p) {
    mpz_init_set_str(curve_a, a, 10);
    mpz_init_set_str(curve_b, b, 10);
    mpz_init_set_str(curve_p, p, 10);
}




struct Point {
    mpz_t x, y;
    bool is_infinity;

    Point()                  { mpz_inits(x, y, NULL); is_infinity = true; }
    ~Point()                 { mpz_clears(x, y, NULL); }

    void set(const Point& o) {
        mpz_set(x, o.x); mpz_set(y, o.y);
        is_infinity = o.is_infinity;
    }
    void set_str(const char* xs, const char* ys) {
        mpz_set_str(x, xs, 10);
        mpz_set_str(y, ys, 10);
        is_infinity = false;
    }
    bool equals(const Point& o) const {
        if (is_infinity && o.is_infinity) return true;
        if (is_infinity || o.is_infinity)  return false;
        return mpz_cmp(x, o.x) == 0 && mpz_cmp(y, o.y) == 0;
    }
};

Point G, Q;

struct Ctx {
    mpz_t num, den, m, rx, ry, diff, sub;
    Ctx()  { mpz_inits(num, den, m, rx, ry, diff, sub, NULL); }
    ~Ctx() { mpz_clears(num, den, m, rx, ry, diff, sub, NULL); }
};




inline void point_add(Point& R, const Point& P, const Point& Pq, Ctx& c) {
    if (P.is_infinity)  { R.set(Pq); return; }
    if (Pq.is_infinity) { R.set(P);  return; }

    if (mpz_cmp(P.x, Pq.x) == 0) {
        if (mpz_cmp(P.y, Pq.y) != 0 || mpz_sgn(P.y) == 0)
            { R.is_infinity = true; return; }
        mpz_mul(c.num, P.x, P.x); mpz_mod(c.num, c.num, curve_p);
        mpz_mul_ui(c.num, c.num, 3);
        mpz_add(c.num, c.num, curve_a);
        mpz_mod(c.num, c.num, curve_p);
        mpz_mul_ui(c.den, P.y, 2);
        mpz_mod(c.den, c.den, curve_p);
    } else {
        mpz_sub(c.num, Pq.y, P.y); mpz_mod(c.num, c.num, curve_p);
        mpz_sub(c.den, Pq.x, P.x); mpz_mod(c.den, c.den, curve_p);
    }

    mpz_invert(c.den, c.den, curve_p);
    mpz_mul(c.m, c.num, c.den); mpz_mod(c.m, c.m, curve_p);

    mpz_mul(c.rx, c.m, c.m);
    mpz_sub(c.rx, c.rx, P.x); mpz_sub(c.rx, c.rx, Pq.x);
    mpz_mod(c.rx, c.rx, curve_p);
    if (mpz_sgn(c.rx) < 0) mpz_add(c.rx, c.rx, curve_p);

    mpz_sub(c.diff, P.x, c.rx);
    mpz_mul(c.ry, c.m, c.diff);
    mpz_sub(c.ry, c.ry, P.y);
    mpz_mod(c.ry, c.ry, curve_p);
    if (mpz_sgn(c.ry) < 0) mpz_add(c.ry, c.ry, curve_p);

    mpz_set(R.x, c.rx); mpz_set(R.y, c.ry);
    R.is_infinity = false;
}

void point_mul(Point& R, const Point& P, const mpz_t k, Ctx& c) {
    R.is_infinity = true;
    Point tmp; tmp.set(P);
    size_t bits = mpz_sizeinbase(k, 2);
    for (int i = (int)bits - 1; i >= 0; --i) {
        if (!R.is_infinity) point_add(R, R, R, c);
        if (mpz_tstbit(k, i))  point_add(R, R, tmp, c);
    }
}






static constexpr int WALK_R = 32;
static constexpr unsigned long DIST_MASK = (1UL << 20) - 1;


Point precomp[WALK_R];
mpz_t  pre_a[WALK_R], pre_b[WALK_R]; // coefficients


void init_precomp() {
    Ctx c;
    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, 0xDEADBEEF);

    for (int i = 0; i < WALK_R; i++) {
        mpz_init(pre_a[i]); mpz_init(pre_b[i]);
        mpz_urandomm(pre_a[i], rng, group_order);
        mpz_urandomm(pre_b[i], rng, group_order);

        Point aG, bQ;
        point_mul(aG, G, pre_a[i], c);
        point_mul(bQ, Q, pre_b[i], c);
        point_add(precomp[i], aG, bQ, c);
    }
    gmp_randclear(rng);
}


inline void walk_step(Point& P, mpz_t& a, mpz_t& b, Ctx& c) {
    unsigned idx = (unsigned)(mpz_get_ui(P.x) & (WALK_R - 1));
    mpz_add(a, a, pre_a[idx]); mpz_mod(a, a, group_order);
    mpz_add(b, b, pre_b[idx]); mpz_mod(b, b, group_order);
    point_add(P, P, precomp[idx], c);
}



// Distinguished point

inline bool is_dp(const Point& P) {
    return !P.is_infinity && (mpz_get_ui(P.x) & DIST_MASK) == 0;
}


shared_mutex                                   g_rwmutex;
unordered_map<string, pair<string,string>>     g_table;
ofstream                                       g_file;
atomic<bool>                                   g_found{false};
mpz_t                                          g_key;
atomic<uint64_t>                               g_iters{0};

atomic<unsigned long long> total_iterations(0);




void worker(int tid) {
    Ctx c;
    Point P, aG, bQ;
    mpz_t a, b, ac, bc, r, ri, num;
    mpz_inits(a, b, ac, bc, r, ri, num, NULL);

    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, (unsigned long)time(nullptr) ^ tid * 6364136223846793005ULL);

    uint64_t local = 0;

    while (!g_found) {
        // random start: P = a·G + b·Q
        mpz_urandomm(a, rng, group_order);
        mpz_urandomm(b, rng, group_order);
        if (!mpz_sgn(a)) mpz_add_ui(a, a, 1);
        if (!mpz_sgn(b)) mpz_add_ui(b, b, 1);

        point_mul(aG, G, a, c);
        point_mul(bQ, Q, b, c);
        point_add(P, aG, bQ, c);

        // Random walk -> distinguished point
        while (!g_found) {
            walk_step(P, a, b, c);
            ++local;

            if (local % (1u << 16) == 0) {
                g_iters += (1u << 16); local = 0;
            }

            if (!is_dp(P)) continue;

            // ──  DP ─────────────────────────────────────────────
            g_iters += local; local = 0;
            char* xs = mpz_get_str(nullptr, 16, P.x);
            string key(xs); free(xs);

            // read lock
            {
                shared_lock<shared_mutex> rl(g_rwmutex);
                auto it = g_table.find(key);
                if (it == g_table.end()) {
                    rl.unlock(); 
                    char* as = mpz_get_str(nullptr, 16, a);
                    char* bs = mpz_get_str(nullptr, 16, b);
                    {
                        unique_lock<shared_mutex> wl(g_rwmutex);
                        g_table[key] = {as, bs};
                        g_file << key << ' ' << as << ' ' << bs << '\n';
                        if (g_table.size() % 500 == 0) g_file.flush();
                    }
                    free(as); free(bs);
                    break; 
                }
                mpz_set_str(ac, it->second.first.c_str(),  16);
                mpz_set_str(bc, it->second.second.c_str(), 16);
            }
            if (g_found) break;


            // k = (a - ac) * (bc - b)^{-1} mod n
            mpz_sub(r, bc, b);  mpz_mod(r, r, group_order);
            if (mpz_sgn(r) < 0) mpz_add(r, r, group_order);

            if (mpz_sgn(r)) {
                mpz_invert(ri, r, group_order);
                mpz_sub(num, a, ac); mpz_mod(num, num, group_order);
                if (mpz_sgn(num) < 0) mpz_add(num, num, group_order);

                unique_lock<shared_mutex> wl(g_rwmutex);
                if (!g_found) {
                    mpz_mul(g_key, num, ri);
                    mpz_mod(g_key, g_key, group_order);
                    g_found = true;
                }
            }
            break;
        }
    }

    g_iters += local;
    mpz_clears(a, b, ac, bc, r, ri, num, NULL);
    gmp_randclear(rng);
}




void monitor() {
    auto t0 = chrono::steady_clock::now();
    while (!g_found) {
        this_thread::sleep_for(chrono::seconds(5));
        if (g_found) break;
        double sec = chrono::duration<double>(
            chrono::steady_clock::now() - t0).count();
        size_t pts;
        {
            shared_lock<shared_mutex> rl(g_rwmutex);
            pts = g_table.size();
        }
        printf("\r  iters: %llu | DPs: %zu | speed: %.2f M/s   ",
               (unsigned long long)g_iters.load(), pts,
               g_iters.load() / sec / 1e6);
        fflush(stdout);
    }
}



int main() {
    puts("============================================================");
    puts("  Optimized Pollard Rho — van Oorschot–Wiener, secp256k1");
    puts("============================================================\n");

    mpz_inits(curve_a, curve_b, curve_p, group_order, g_key, NULL);

    // secp256k1
    mpz_set_str(curve_a, "0", 10);
    mpz_set_str(curve_b, "7", 10);
    mpz_set_str(curve_p,
        "115792089237316195423570985008687907853269984665640564039"
        "457584007908834671663", 10);
    mpz_set_str(group_order,
        "115792089237316195423570985008687907852837564279074904382"
        "605163141518161494337", 10);

    G.set_str(
        "55066263022277343669578718895168534326250603453777594175500187360389116729240",
        "32670510020758816978083085130507043184471273380659243275938904335757337482424");
    Q.set_str(
        "71715483259960251104886616086749964585141029199321453995874288114660996885376",
        "40560706241088654877709886733221375330368142163939632832049960868661138127382");

    printf("[*] Initializing %d precomputed walk points...\n", WALK_R);
    init_precomp();

    {
        ifstream f("progress.txt");
        if (f) {
            string x, a, b;
            while (f >> x >> a >> b) g_table[x] = {a, b};
            printf("[*] Loaded %zu distinguished points from disk\n",
                   g_table.size());
        }
    }
    g_file.open("progress.txt", ios::app);

    unsigned n = max(1u, thread::hardware_concurrency());
    printf("[*] Launching %u threads\n\n", n);

    thread mon(monitor);
    vector<thread> pool;
    pool.reserve(n);
    for (unsigned i = 0; i < n; ++i) pool.emplace_back(worker, i);
    
    for (auto& t : pool) t.join();
    mon.join();

    if (g_found) {
        printf("\n=========================================\n");
        gmp_printf("[+] Secret key found:\n    k = %Zd\n", g_key);
        printf("=========================================\n");
    }

    g_file.flush();
    g_file.close();

    for (int i = 0; i < WALK_R; i++) {
        mpz_clears(pre_a[i], pre_b[i], NULL);
    }
    mpz_clears(curve_a, curve_b, curve_p, group_order, g_key, NULL);

    return 0;
}