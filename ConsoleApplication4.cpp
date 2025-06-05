#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

using namespace std;

class Wielomian {
private:
    vector<double> wsp;

public:
    Wielomian(const vector<double>& wspolczynniki) : wsp(wspolczynniki) {
        if (wsp.empty())
            throw invalid_argument("Wielomian nie może być pusty.");
        while (wsp.size() > 1 && abs(wsp.back()) < 1e-12)
            wsp.pop_back();
    }

    int stopien() const {
        return wsp.size() - 1;
    }

    string toString() const {
        ostringstream oss;
        oss << "W(x) = ";
        bool pierwsze = true;

        for (int i = stopien(); i >= 0; --i) {
            double c = wsp[i];
            if (abs(c) < 1e-12) continue;

            if (!pierwsze) oss << (c >= 0 ? " + " : " - ");
            else if (c < 0) oss << "-";

            if (abs(c) != 1 || i == 0) oss << abs(c);
            if (i > 0) oss << "x" << (i > 1 ? "^" + to_string(i) : "");

            pierwsze = false;
        }

        return pierwsze ? "W(x) = 0" : oss.str();
    }

    double operator()(double x) const {
        double wynik = 0;
        for (int i = wsp.size() - 1; i >= 0; --i)
            wynik = wynik * x + wsp[i];
        return wynik;
    }

    Wielomian operator+(const Wielomian& o) const {
        size_t n = max(wsp.size(), o.wsp.size());
        vector<double> wynik(n, 0);
        for (size_t i = 0; i < wsp.size(); ++i) wynik[i] += wsp[i];
        for (size_t i = 0; i < o.wsp.size(); ++i) wynik[i] += o.wsp[i];
        return Wielomian(wynik);
    }

    Wielomian operator-(const Wielomian& o) const {
        size_t n = max(wsp.size(), o.wsp.size());
        vector<double> wynik(n, 0);
        for (size_t i = 0; i < wsp.size(); ++i) wynik[i] += wsp[i];
        for (size_t i = 0; i < o.wsp.size(); ++i) wynik[i] -= o.wsp[i];
        return Wielomian(wynik);
    }

    Wielomian operator*(const Wielomian& o) const {
        vector<double> wynik(wsp.size() + o.wsp.size() - 1, 0);
        for (size_t i = 0; i < wsp.size(); ++i)
            for (size_t j = 0; j < o.wsp.size(); ++j)
                wynik[i + j] += wsp[i] * o.wsp[j];
        return Wielomian(wynik);
    }

    Wielomian& operator+=(const Wielomian& o) { return *this = *this + o; }
    Wielomian& operator-=(const Wielomian& o) { return *this = *this - o; }
    Wielomian& operator*=(const Wielomian& o) { return *this = *this * o; }
};

int main() {
    Wielomian w1({ 1, 2, 3 });     // 3x^2 + 2x + 1
    Wielomian w2({ -1, 0, 1 });    // x^2 - 1

    cout << w1.toString() << endl;
    cout << w2.toString() << endl;
    cout << "Suma:      " << (w1 + w2).toString() << endl;
    cout << "Różnica:   " << (w1 - w2).toString() << endl;
    cout << "Iloczyn:   " << (w1 * w2).toString() << endl;
    cout << "Wartość w1(2) = " << w1(2.0) << endl;

    w1 += w2;
    cout << "w1 += w2:  " << w1.toString() << endl;

    return 0;
}
