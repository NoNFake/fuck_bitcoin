//  g++ -O3 -march=native test.cpp -o pohlig_hellman -lgmpxx -lgmp -pthread
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <thread>
#include <mutex>
#include <unordered_map>
#include <string>
#include <random>
#include <gmp.h>
#include <atomic>


using namespace std;
mpz_t curve_a, curve_b, curve_p, group_order;

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


Point G, Q;


struct ThreadMathContext {
    mpz_t tmp_num, tmp_den, tmp_m, tmp_rx, tmp_ry, tmp_diff, subset;

    ThreadMathContext() {
        mpz_inits(tmp_num, tmp_den, tmp_m, tmp_rx, tmp_ry, tmp_diff, subset, NULL);
    }
    
    ~ThreadMathContext() {
        mpz_clears(tmp_num, tmp_den, tmp_m, tmp_rx, tmp_ry, tmp_diff, subset, NULL);
    }
};




// mpz_t tmp_num, tmp_den, tmp_m, tmp_rx, tmp_ry, tmp_diff;


// void init_math(){
//     mpz_inits(tmp_num, tmp_den, tmp_m, tmp_rx, tmp_ry, tmp_diff, NULL);
// }



void add_points(Point& res, const Point& P1, const Point& P2, ThreadMathContext& ctx) {
    if (P1.is_infinity) { res.set(P2); return; }
    if (P2.is_infinity) { res.set(P1); return; }


    if (mpz_cmp(P1.x, P2.x) == 0) {
        if (mpz_cmp(P1.y, P2.y) != 0 || mpz_cmp_ui(P1.y, 0) == 0) {
            res.is_infinity = true;
            return;
        }

        mpz_mul(ctx.tmp_num, P1.x, P1.x);
        mpz_mod(ctx.tmp_num, ctx.tmp_num, curve_p);
        mpz_mul_ui(ctx.tmp_num, ctx.tmp_num, 3);
        mpz_add(ctx.tmp_num, ctx.tmp_num, curve_a);
        mpz_mod(ctx.tmp_num, ctx.tmp_num, curve_p);

        mpz_mul_ui(ctx.tmp_den, P1.y, 2);
        mpz_mod(ctx.tmp_den, ctx.tmp_den, curve_p);

    } else {
        mpz_sub(ctx.tmp_num, P2.y, P1.y);
        mpz_mod(ctx.tmp_num, ctx.tmp_num, curve_p);
        
        mpz_sub(ctx.tmp_den, P2.x, P1.x);
        mpz_mod(ctx.tmp_den, ctx.tmp_den, curve_p);
    }

    if (mpz_invert(ctx.tmp_den, ctx.tmp_den, curve_p) == 0) {
        throw runtime_error("No inverse found in point addition");
    }




    mpz_mul(ctx.tmp_m, ctx.tmp_num, ctx.tmp_den);
    mpz_mod(ctx.tmp_m, ctx.tmp_m, curve_p);

    // rx = m^2 - P.x - Q.x
    mpz_mul(ctx.tmp_rx, ctx.tmp_m, ctx.tmp_m);
    mpz_sub(ctx.tmp_rx, ctx.tmp_rx, P1.x);
    mpz_sub(ctx.tmp_rx, ctx.tmp_rx, P2.x);
    mpz_mod(ctx.tmp_rx, ctx.tmp_rx, curve_p);
    if (mpz_sgn(ctx.tmp_rx) < 0) mpz_add(ctx.tmp_rx, ctx.tmp_rx, curve_p);

    // ry = m*(P.x - rx) - P.y
    mpz_sub(ctx.tmp_diff, P1.x, ctx.tmp_rx);
    mpz_mul(ctx.tmp_ry, ctx.tmp_m, ctx.tmp_diff);
    mpz_sub(ctx.tmp_ry, ctx.tmp_ry, P1.y);
    mpz_mod(ctx.tmp_ry, ctx.tmp_ry, curve_p);
    if (mpz_sgn(ctx.tmp_ry) < 0) mpz_add(ctx.tmp_ry, ctx.tmp_ry, curve_p);

    mpz_set(res.x, ctx.tmp_rx);
    mpz_set(res.y, ctx.tmp_ry);
    res.is_infinity = false;

}

void multiply_point(Point& res, const Point& P, const mpz_t k, ThreadMathContext& ctx) {
    Point R, tempP;
    tempP.set(P);
    
    size_t bits = mpz_sizeinbase(k, 2);
    for (int i = bits - 1; i >= 0; --i) {
        add_points(R, R, R, ctx);
        if (mpz_tstbit(k, i)) {
            add_points(R, R, tempP, ctx);
        }
    }
    res.set(R);
}


void next_step_thread_safe(Point& P, mpz_t& a, mpz_t& b, ThreadMathContext& ctx) {
    if (P.is_infinity) return;
    
    mpz_fdiv_r_ui(ctx.subset, P.x, 3); 
    unsigned long s = mpz_get_ui(ctx.subset);

    if (s == 0) {
        mpz_add_ui(a, a, 1);
        mpz_mod(a, a, group_order);
        add_points(P, P, G, ctx);
    } else if (s == 1) {
        mpz_mul_ui(a, a, 2); mpz_mod(a, a, group_order);
        mpz_mul_ui(b, b, 2); mpz_mod(b, b, group_order);
        add_points(P, P, P, ctx);
    } else {
        mpz_add_ui(b, b, 1);
        mpz_mod(b, b, group_order);
        add_points(P, P, Q, ctx);
    }
}


bool is_distinguished(const Point& P) {
    if (P.is_infinity) return false;
    unsigned long mask = 0xFFFFF; // 20 bit
    return (mpz_get_ui(P.x) & mask) == 0;
}



mutex db_mutex;
unordered_map<string, pair<string, string>> dp_table;
ofstream save_file;
bool key_found = false;
mpz_t final_secret_key;

atomic<unsigned long long> total_iterations(0);

void monitor_thread() {
    auto start_time = chrono::steady_clock::now();
    
    while (!key_found) {
        this_thread::sleep_for(chrono::milliseconds(500));
        if (key_found) break;

        auto now = chrono::steady_clock::now();
        chrono::duration<double> elapsed = now - start_time;
        
        double speed = (elapsed.count() > 0) ? (total_iterations.load() / elapsed.count() / 1000000.0) : 0;
        
        size_t current_table_size = 0;
        {
            lock_guard<mutex> lock(db_mutex);
            current_table_size = dp_table.size();
        }

        cout << "\r      ... completed: " << total_iterations.load() 
             << " steps | points: " << current_table_size 
             << " | speed: " << speed << " mill/s ... " << flush;
    }
}


void worker_thread(int thread_id) {
    ThreadMathContext ctx;
    Point P, aG, bQ;
    mpz_t a, b, a_col, b_col, r, r_inv, num;
    mpz_inits(a, b, a_col, b_col, r, r_inv, num, NULL);

    // Уникальный сид для потока
    random_device rd;
    unsigned long seed = rd() ^ (thread_id * 1000) ^ time(NULL);
    gmp_randstate_t rand_state;
    gmp_randinit_default(rand_state);
    gmp_randseed_ui(rand_state, seed);


    unsigned long long local_steps = 0;


    while (!key_found) {
        // 1. Прыгаем в случайную точку: P = a*G + b*Q
        mpz_urandomm(a, rand_state, group_order);
        mpz_urandomm(b, rand_state, group_order);
        
        if (mpz_cmp_ui(a, 0) == 0) mpz_add_ui(a, a, 1);
        if (mpz_cmp_ui(b, 0) == 0) mpz_add_ui(b, b, 1);

        multiply_point(aG, G, a, ctx);
        multiply_point(bQ, Q, b, ctx);
        add_points(P, aG, bQ, ctx);

        // 2. Блуждаем, пока не найдем отличимую точку
        while (!key_found) {
            next_step_thread_safe(P, a, b, ctx);
            

            local_steps++;
            if (local_steps >= 50000) {
                total_iterations += local_steps;
                local_steps = 0;
            }


            if (is_distinguished(P)) {
                total_iterations += local_steps;
                local_steps = 0;


                lock_guard<mutex> lock(db_mutex);
                if (key_found) break; // Кто-то другой уже нашел ключ

                char* x_str = mpz_get_str(NULL, 16, P.x);
                string hash_key(x_str);
                
                if (dp_table.count(hash_key)) {
                    // КОЛЛИЗИЯ НАЙДЕНА!
                    mpz_set_str(a_col, dp_table[hash_key].first.c_str(), 16);
                    mpz_set_str(b_col, dp_table[hash_key].second.c_str(), 16);

                    // r = (b_col - b) mod N
                    mpz_sub(r, b_col, b);
                    mpz_mod(r, r, group_order);
                    if (mpz_sgn(r) < 0) mpz_add(r, r, group_order);

                    if (mpz_cmp_ui(r, 0) != 0) {
                        mpz_invert(r_inv, r, group_order);

                        // num = (a - a_col) mod N
                        mpz_sub(num, a, a_col);
                        mpz_mod(num, num, group_order);
                        if (mpz_sgn(num) < 0) mpz_add(num, num, group_order);

                        // result = (num * r_inv) mod N
                        mpz_mul(final_secret_key, num, r_inv);
                        mpz_mod(final_secret_key, final_secret_key, group_order);
                        
                        key_found = true;
                    }
                } else {
                    // Новая точка — сохраняем
                    char* a_str = mpz_get_str(NULL, 16, a);
                    char* b_str = mpz_get_str(NULL, 16, b);
                    
                    dp_table[hash_key] = {a_str, b_str};
                    save_file << hash_key << " " << a_str << " " << b_str << "\n";
                    save_file.flush(); // Важно! Сразу пишем на диск
                    
                    free(a_str); free(b_str);
                }
                free(x_str);
                
                // Начинаем новый случайный путь
                break; 
            }
        }
    }

    total_iterations += local_steps;
    mpz_clears(a, b, a_col, b_col, r, r_inv, num, NULL);
    gmp_randclear(rand_state);
}



int main() {
    cout << "============================================================\n";
    cout << "     Multi-Threaded Pollard's Rho (Van Oorschot-Wiener)\n";
    cout << "============================================================\n\n";

    mpz_inits(curve_a, curve_b, curve_p, group_order, final_secret_key, NULL);
    mpz_set_str(curve_a, "0", 10);
    mpz_set_str(curve_b, "7", 10);
    mpz_set_str(curve_p, "115792089237316195423570985008687907853269984665640564039457584007908834671663", 10);
    mpz_set_str(group_order, "115792089237316195423570985008687907852837564279074904382605163141518161494337", 10);

    G.set_str(
        "55066263022277343669578718895168534326250603453777594175500187360389116729240",
        "32670510020758816978083085130507043184471273380659243275938904335757337482424"
    );

    Q.set_str(
        "71715483259960251104886616086749964585141029199321453995874288114660996885376",
        "40560706241088654877709886733221375330368142163939632832049960868661138127382"
    );


    ifstream in_file("progress.txt");
    if (in_file.is_open()) {
        string x, a, b;
        while (in_file >> x >> a >> b) {
            dp_table[x] = {a, b};
        }
        in_file.close();
        cout << "[*] loaded " << dp_table.size() << " points fro progress.txt\n";
    }
    save_file.open("progress.txt", ios::app);

    unsigned int num_cores = thread::hardware_concurrency();
    if (num_cores == 0) num_cores = 4;
    cout << "[*] started on " << num_cores << " cores\n";


    thread monitor(monitor_thread);


   vector<thread> threads;
    for (unsigned int i = 0; i < num_cores; ++i) {
        threads.push_back(thread(worker_thread, i));
    }

    for (auto& th : threads) {
        th.join();
    }

    monitor.join();

    if (key_found) {
        cout << "\n=========================================\n";
        gmp_printf("[+] secret key found :\n    k = %Zd\n", final_secret_key);
        cout << "=========================================\n";
    }

    save_file.close();
    mpz_clears(curve_a, curve_b, curve_p, group_order, final_secret_key, NULL);
    return 0;

}   