#include <iostream>
#include <fstream>
#include <ctime>
#include <vector>

using namespace std;

typedef unsigned int uint;

int col; //Число столбцов
int str; //Число строк
vector<int> k; //Мощности множеств P_i
vector<bool> chain; //True, если цепь. False, если антицепь.
int num_elem = 0;
int max_anti_chain = 0;
int curr_anti_chain = 0;

const int bit_size = sizeof(uint) * 8;

ofstream out("output.txt");
ifstream in("input.txt");

uint ones;

class Bit_vector{
    int len;
public:
    uint *data;
    
    Bit_vector();
    
    Bit_vector(int len);
    
    Bit_vector(Bit_vector &B);
    
    void set_len(int len);
    
    void del_str_from_Rb();
    
    void del_nei_col(int i);
    
    void set_one(int i);
    
    void set_zero(int i);
    
    int get_len();
    
    bool get(int i);
    
    bool is_nil();
    
    void do_support_str(int j);
    
    void del_str(int j);
    
    //~Bit_vector();
};

void build_subtree();

int create_node(int j);


uint *mask;

void init_mask(){
    mask = new uint[bit_size];
    mask[0] = 1;
    for (int i = 1; i < bit_size; ++i){
        mask[i] = mask[i-1] << 1;
    }
}

Bit_vector::Bit_vector(){}

Bit_vector::Bit_vector(int len){
    this->len = len;
    if (len % bit_size == 0){
        data = new uint[len / bit_size];
        for (int i = 0; i < len / bit_size; ++i){
            data[i] = 0;
        }
    }
    else{
        data = new uint[len / bit_size + 1];
        for (int i = 0; i < len / bit_size + 1; ++i){
            data[i] = 0;
        }
    }
}

int Bit_vector::get_len(){
    return len;
}

Bit_vector::Bit_vector(Bit_vector &B){
    len = B.len;
    if (len % bit_size == 0){
        data = new uint[len / bit_size];
        for (int i = 0; i < len / bit_size; ++i){
            data[i] = B.data[i];
        }
    }
    else{
        data = new uint[len / bit_size + 1];
        for (int i = 0; i < len / bit_size + 1; ++i){
            data[i] = B.data[i];
        }
    }
}

void Bit_vector::set_len(int len){
    this->len = len;
    if (len % bit_size == 0){
        data = new uint[len / bit_size];
        for (int i = 0; i < len / bit_size; ++i){
            data[i] = 0;
        }
    }
    else{
        data = new uint[len / bit_size + 1];
        for (int i = 0; i < len / bit_size + 1; ++i){
            data[i] = 0;
        }
    }
}

bool Bit_vector::get(int i){
    int n = i / bit_size;
    return (data[n] & mask[i % bit_size]) > 0;
}



Bit_vector H(col), R(str), D(col), *M, *S, *Mst;




void Bit_vector::del_str_from_Rb(){
    for (int i = 0; i <  (col - 1) / bit_size + 1; ++i){
        R.data[i] &= ~data[i];
    }
}

void Bit_vector::set_one(int i){
    int n = i / bit_size;
    data[n] |= mask[i % bit_size];
}

void Bit_vector::set_zero(int i){
    int n = i / bit_size;
    data[n] &= ~mask[i % bit_size];
}

/*void Bit_vector::del_nei_col(int i){
 int start = i/k;
 for (int j = start*k; j < (start + 1) * k; ++j){
 this->set_zero(j);
 }
 }*/
void Bit_vector::del_nei_col(int i){
    
    int pi_prev = 0;
    
    for (int pi : k){
        if (i < pi_prev + pi){
            for (int j = pi_prev; j < pi + pi_prev; ++j){
                this->set_zero(j);
            }
            return;
        }
        else{
            pi_prev += pi;
        }
    }
}


bool Bit_vector::is_nil(){
    for (int i = 0; i < (len-1) / bit_size + 1; ++i){
        if (data[i]) return false;
    }
    return true;
}

void Bit_vector::do_support_str(int j){
    if (chain[j]){
        for (int i = 0; i < (str - 1)/bit_size + 1; ++i){
            data[i] = R.data[i] & (~Mst[j+1].data[i]) & Mst[j].data[i];
        }
    }
    else{
        for (int i = 0; i < (str - 1)/bit_size + 1; ++i){
            data[i] = R.data[i] & Mst[j].data[i];
        }
    }
}

void Bit_vector::del_str(int j){
    for (int i = 0; i < (str - 1)/bit_size + 1; ++i){
        data[i] = ~Mst[j].data[i] & data[i];
    }
}

bool is_forg_col(int j1, int j2){
    for (int i = 0; i < (str - 1)/bit_size + 1; ++i){
        if (S[j2].data[i] != (Mst[j1].data[i] & S[j2].data[i])) return false;
    }
    return true;
}

void clean_prev_S(int j_prev, int j_curr){
    for (int i = 0; i < (str - 1)/bit_size + 1; ++i){
        S[j_prev].data[i] = S[j_prev].data[i] & (~Mst[j_curr].data[i]);
    }
}



/*Bit_vector::~Bit_vector(){
 delete data;
 }*/

int find_min_weight(){
    int h, t;
    int w1 = col;
    int st = -1;
    int w, k;
    for (int s = 0; s < str; ++s){
        w = 0;
        if (!R.get(s)) continue;
        for (int i = 0; i < (col-1) / bit_size + 1; ++i){
            t = 1 << (bit_size-2);
            h = M[s].data[i];
            k = bit_size - 2;
            while (t > 0){
                if (t <= h){
                    h &= ~mask[k];
                    ++w;
                }
                t >>= 1;
                --k;
            }
        }
        if (w < w1){
            w1 = w;
            st = s;
        }
    }
    return st;
}


/*void print_H(){
 for (int i = 0; i < col/k; ++i){
 bool f = true;
 for (int j = 0; j < k; ++j){
 if (H.get(i*k + j)){
 f = false;
 out << j << " ";
 }
 }
 if (f) out << k-1 << " ";
 }
 out << "\n";
 }*/

void print_H(){
    int c = 0;
    for (int pi : k){
        bool f = true;
        for (int j = 0; j < pi; ++j){
            if (H.get(c + j)){
                f = false;
                out << j << " ";
            }
        }
        if (chain[c]){
            if (f) out << pi-1 << " ";
        }
        else{
            if (f) out << -1 << " ";
        }
        c = c + pi;
    }
    out << "\n";
}

void print(Bit_vector B){
    for (int i = 0; i < B.get_len(); ++i){
        if (B.get(i)) cout << i << " ";
    }
    cout << "\n";
}

void build_subtree(){
    int i;
    if ((i = find_min_weight()) < 0){
        return;
    }
    for (int j = 0; j < col; ++j){
        if (!D.get(j)) continue;
        if (!M[i].get(j)) continue;
        D.set_zero(j);
        Bit_vector H1(H), D1(D), R1(R), *S1;
        S1 = new Bit_vector[col];
        for (int p = 0; p < col; ++p){
            int len = S[p].get_len();
            S1[p].set_len(len);
            if (len % bit_size == 0){
                S1[p].data = new uint[len / bit_size];
                for (int t = 0; t < len / bit_size; ++t){
                    S1[p].data[t] = S[p].data[t];
                }
            }
            else{
                S1[p].data = new uint[len / bit_size + 1];
                for (int t = 0; t < len / bit_size + 1; ++t){
                    S1[p].data[t] = S[p].data[t];
                }
            }
        }
        D.del_nei_col(j);
        if (!chain[j]){
            ++curr_anti_chain;
        }
        
        if (max_anti_chain != -1){
            if (curr_anti_chain >= max_anti_chain){
                for (int p = 0; p < col; ++p){
                    if (!chain[p]){
                        D.set_zero(p);
                    }
                }
            }
        }
        if (create_node(j)){
            if (R.is_nil()){
                print_H();
                ++num_elem;
            }
            else{
                build_subtree();
            }
            D = D1;
            H = H1;
            R = R1;
            for (int p = 0; p < col; ++p){
                int len = S[p].get_len();
                if (len % bit_size == 0){
                    for (int t = 0; t < len / bit_size; ++t){
                        S[p].data[t] = S1[p].data[t];
                    }
                }
                else{
                    for (int t = 0; t < len / bit_size + 1; ++t){
                        S[p].data[t] = S1[p].data[t];
                    }
                }
            }
            delete S1;
        }
        else{
            D = D1;
            H = H1;
            R = R1;
            for (int p = 0; p < col; ++p){
                int len = S[p].get_len();
                if (len % bit_size == 0){
                    for (int t = 0; t < len / bit_size; ++t){
                        S[p].data[t] = S1[p].data[t];
                    }
                }
                else{
                    for (int t = 0; t < len / bit_size + 1; ++t){
                        S[p].data[t] = S1[p].data[t];
                    }
                }
            }
            delete S1;
        }
        
        if (!chain[j]){
            --curr_anti_chain;
        }
    }
}

int create_node(int j){
    S[j].do_support_str(j);
    if (S[j].is_nil()) return 0;
    for (int i=0; i < col; ++i){
        if (i == j){
            continue;
        }
        if (H.get(i)){
            clean_prev_S(i, j);
            if (S[i].is_nil()) return 0;
        }
    }
    R.del_str(j);
    for (int i = 0; i < col; ++i){
        if (D.get(i)){
            if (is_forg_col(i, j)){
                D.set_zero(i);
            }
        }
    }
    H.set_one(j);
    return 1;
}


int main(){
    int n, m, **L, P_i;
    in >> max_anti_chain >> m >> n;
    str = m;
    //col = k * n;
    col = 0;
    for (int i = 0; i < n; ++i){
        in >> P_i;
        k.push_back(P_i);
        col += P_i;
    }
    int ch;
    for (int pi : k){
        in >> ch;
        if (ch == 1){
            for (int i = 0; i < pi; ++i){
                chain.push_back(true);
            }
        }
        else{
            for (int i = 0; i < pi; ++i){
                chain.push_back(false);
            }
        }
    }
    init_mask();
    H = Bit_vector(col);
    R = Bit_vector(str);
    D = Bit_vector(col);
    
    for (int i = 0; i < col; ++i){
        D.set_one(i);
    }
    
    for (int i = 0; i < str; ++i){
        R.set_one(i);
    }
    
    M = new Bit_vector[str];
    S = new Bit_vector[col];
    Mst = new Bit_vector[col];
    
    L = new int*[str];
    
    for (int i = 0; i < str; ++i){
        L[i] = new int[n];
        for (int j = 0; j < n; ++j){
            in >> L[i][j];
        }
    }
    
    for (int i = 0; i < str; ++i){
        M[i].set_len(col);
    }
    
    for (int i = 0; i < str; ++i){
        int c = 0;
        for (int j = 0; j < n; ++j){
            if (chain[c]){
                for (int s = 0; s < k[j]; ++s){
                    if (L[i][j] > 0) M[i].set_one(c + s);
                    L[i][j] -= 1;
                }
            }
            else{
                for (int s = 0; s < k[j]; ++s){
                    M[i].set_one(c + s);
                }
                M[i].set_zero(c + L[i][j]);
            }
            c = c + k[j];
        }
    }
    
    
    for (int i = 0; i < col; ++i){
        Mst[i].set_len(str);
        for (int j = 0; j < str; ++j){
            if (M[j].get(i)) Mst[i].set_one(j);
        }
    }
    
    for (int i = 0; i < col; ++i){
        S[i].set_len(str);
    }
    
    in.close();
    
    int start = clock();
    
    if (max_anti_chain == 0){
        for (int j = 0; j < col; ++j){
            if (!chain[j]){
                D.set_zero(j);
            }
        }
    }
    
    build_subtree();
    
    int finish = clock();
    
    //cout << finish - start << "\n";
    
    //cout << num_elem <<"\n";
    out.close();
    
    
}


















