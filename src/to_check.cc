#include <bits/stdc++.h>
using namespace std;

using ll = long long;

ll modpow(ll a, ll e, ll mod) {
    ll r = 1;
    while (e) {
        if (e & 1) r = (r * a) % mod;
        a = (a * a) % mod;
        e >>= 1;
    }
    return r;
}

// ---------- Polynômes sur F_p ----------
struct Poly {
    vector<int> c; // coefficients croissants

    Poly(int n = 0) : c(n) {}

    int deg() const {
        for (int i = (int)c.size() - 1; i >= 0; i--)
            if (c[i]) return i;
        return 0;
    }
};

int mod(int x, int p) {
    x %= p;
    if (x < 0) x += p;
    return x;
}

Poly trim(Poly a) {
    while (!a.c.empty() && a.c.back() == 0) a.c.pop_back();
    return a;
}

// addition
Poly add(const Poly &a, const Poly &b, int p) {
    Poly r(max(a.c.size(), b.c.size()));
    for (size_t i = 0; i < r.c.size(); i++) {
        int va = i < a.c.size() ? a.c[i] : 0;
        int vb = i < b.c.size() ? b.c[i] : 0;
        r.c[i] = mod(va + vb, p);
    }
    return trim(r);
}

// multiplication modulo f
Poly mulmod(Poly a, Poly b, const Poly &modf, int p) {
    Poly r(a.c.size() + b.c.size() + 5);

    for (size_t i = 0; i < a.c.size(); i++)
        for (size_t j = 0; j < b.c.size(); j++)
            r.c[i + j] = mod(r.c[i + j] + a.c[i] * b.c[j], p);

    // réduction modulo modf (division polynomiale simplifiée)
    int df = modf.deg();

    for (int i = (int)r.c.size() - 1; i >= df; i--) {
        if (r.c[i] == 0) continue;
        int coef = r.c[i];
        for (int j = 0; j <= df; j++) {
            r.c[i - df + j] = mod(r.c[i - df + j] - coef * modf.c[j], p);
        }
        r.c[i] = 0;
    }

    return trim(r);
}

// exponentiation x^(p^k) mod f
Poly frob(Poly x, ll exp, const Poly &f, int p) {
    Poly base = x;
    Poly res;
    res.c = {1};

    while (exp) {
        if (exp & 1) res = mulmod(res, base, f, p);
        base = mulmod(base, base, f, p);
        exp >>= 1;
    }
    return res;
}

// gcd polynômes
Poly poly_gcd(Poly a, Poly b, int p) {
    auto modpoly = [&](Poly a, Poly b) {
        while (a.deg() >= b.deg()) {
            int d = a.deg() - b.deg();
            int coef = a.c.back() * modpow(b.c.back(), p - 2, p) % p;

            for (int i = 0; i <= b.deg(); i++) {
                a.c[i + d] = mod(a.c[i + d] - coef * b.c[i], p);
            }
            a = trim(a);
        }
        return a;
    };

    while (!b.c.empty()) {
        Poly r = modpoly(a, b);
        a = b;
        b = r;
    }
    return a;
}

// test irréductibilité Rabin
bool irreducible(Poly f, int p, int n) {
    Poly x;
    x.c = {0, 1};

    for (int r = 2; r <= n; r++) {
        if (n % r == 0) {
            Poly t = frob(x, modpow(p, n / r, LLONG_MAX), f, p);
            Poly diff = add(t, x, p);

            if (poly_gcd(f, diff, p).deg() > 0)
                return false;
        }
    }
    return true;
}

// recherche brute d'un polynôme irréductible
Poly find_irred(int p, int n) {
    int limit = 1;
    for (int i = 0; i < n; i++) limit *= p;

    for (int mask = 0; mask < limit; mask++) {
        Poly f(n + 1);
        f.c[n] = 1;

        int t = mask;
        for (int i = 0; i < n; i++) {
            f.c[i] = t % p;
            t /= p;
        }

        if (irreducible(f, p, n))
            return f;
    }

    return {};
}

// factorisation naïve q = p^n
bool decompose(ll q, ll &p, int &n) {
    for (ll i = 2; i * i <= q; i++) {
        int cnt = 0;
        ll t = q;
        while (t % i == 0) {
            t /= i;
            cnt++;
        }
        if (t == 1) {
            p = i;
            n = cnt;
            return true;
        }
    }
    return false;
}

int main() {
    ll q;
    cin >> q;

    ll p;
    int n;

    if (!decompose(q, p, n)) {
        cout << "q n'est pas une puissance de nombre premier\n";
        return 0;
    }

    cout << "p = " << p << ", n = " << n << "\n";

    Poly f = find_irred(p, n);

    cout << "Polynome irreductible : ";
    for (int i = f.deg(); i >= 0; i--) {
        cout << f.c[i] << "x^" << i;
        if (i) cout << " + ";
    }
    cout << "\n";
}
