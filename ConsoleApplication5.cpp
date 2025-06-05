#include <iostream>
#include <string>
#include <stdexcept>
#include <algorithm>

using namespace std;

/** Klasa reprezentująca sekwencję RNA */
class RNASequence;

/** Klasa reprezentująca sekwencję białkową */
class ProteinSequence;

/** Klasa reprezentująca sekwencję DNA */
class DNASequence {
private:
    string identifier;     ///< Identyfikator sekwencji
    string data;           ///< Sekwencja zasad DNA
    static const string VALID_CHARS; ///< Dozwolone znaki zasad

public:
    /**
     * Konstruktor sprawdzający poprawność sekwencji.
     * @param id Identyfikator sekwencji
     * @param seq Sekwencja zasad
     */
    DNASequence(string id, string seq) : identifier(id), data(seq) {
        for (char c : data) {
            if (VALID_CHARS.find(c) == string::npos) // npos - stała reprezentująca brak pozycji (ang. no position)
                throw invalid_argument("Nieprawidłowy znak w sekwencji DNA.");
        }
    }

    /**
     * Zwraca długość sekwencji.
     * @return Liczba znaków w sekwencji
     */
    int length() const {
        return data.length();
    }

    /**
     * Zwraca reprezentację FASTA.
     * @return Sekwencja w formacie FASTA
     */
    string toString() const {
        return ">" + identifier + "\n" + data;
    }

    /**
     * Mutuje zasadę w podanej pozycji.
     * @param position Pozycja do zmiany
     * @param value Nowa zasada
     */
    void mutate(int position, char value) {
        if (position < 0 || position >= data.length())
            throw out_of_range("Pozycja poza zakresem.");
        if (VALID_CHARS.find(value) == string::npos)
            throw invalid_argument("Nieprawidłowy znak mutacji.");
        data[position] = value;
    }

    /**
     * Szuka motywu w sekwencji.
     * @param motif Motyw do znalezienia
     * @return Indeks pierwszego wystąpienia lub -1
     */
    int findMotif(const string& motif) const {
        size_t pos = data.find(motif);
        return pos != string::npos ? static_cast<int>(pos) : -1;
    }

    /**
     * Zwraca nić komplementarną.
     * @return Komplementarna sekwencja DNA
     */
    string complement() const {
        string result = data;
        for (char& c : result) {
            switch (c) {
            case 'A': c = 'T'; break;
            case 'T': c = 'A'; break;
            case 'C': c = 'G'; break;
            case 'G': c = 'C'; break;
            }
        }
        return result;
    }

    /**
     * Transkrybuje nić matrycową DNA do RNA.
     * @return Obiekt RNASequence
     */
    RNASequence transcribe() const;
};

const string DNASequence::VALID_CHARS = "ATCG";

/** Klasa reprezentująca sekwencję RNA */
class RNASequence {
private:
    string identifier;     ///< Identyfikator sekwencji
    string data;           ///< Sekwencja zasad RNA
    static const string VALID_CHARS; ///< Dozwolone znaki zasad

public:
    /**
     * Konstruktor sprawdzający poprawność sekwencji.
     * @param id Identyfikator
     * @param seq Sekwencja RNA
     */
    RNASequence(string id, string seq) : identifier(id), data(seq) {
        for (char c : data) {
            if (VALID_CHARS.find(c) == string::npos)
                throw invalid_argument("Nieprawidłowy znak w sekwencji RNA.");
        }
    }

    /**
     * Zwraca długość sekwencji.
     */
    int length() const {
        return data.length();
    }

    /**
     * Zwraca format FASTA.
     */
    string toString() const {
        return ">" + identifier + "\n" + data;
    }

    /**
     * Mutuje zasadę RNA.
     */
    void mutate(int position, char value) {
        if (position < 0 || position >= data.length())
            throw out_of_range("Pozycja poza zakresem.");
        if (VALID_CHARS.find(value) == string::npos)
            throw invalid_argument("Nieprawidłowy znak mutacji.");
        data[position] = value;
    }

    /**
     * Znajduje motyw w RNA.
     */
    int findMotif(const string& motif) const {
        size_t pos = data.find(motif);
        return pos != string::npos ? static_cast<int>(pos) : -1;
    }

    /**
     * Zwraca komplementarną nić RNA.
     */
    string complement() const {
        string result = data;
        for (char& c : result) {
            switch (c) {
            case 'A': c = 'U'; break;
            case 'U': c = 'A'; break;
            case 'C': c = 'G'; break;
            case 'G': c = 'C'; break;
            }
        }
        return result;
    }

    /**
     * Transluje RNA na uproszczoną sekwencję białkową.
     */
    ProteinSequence transcribe() const;
};

const string RNASequence::VALID_CHARS = "AUCG";

/** Klasa reprezentująca sekwencję białkową */
class ProteinSequence {
private:
    string identifier;     ///< Identyfikator białka
    string data;           ///< Sekwencja aminokwasów
    static const string VALID_CHARS; ///< Dozwolone aminokwasy

public:
    /**
     * Konstruktor białka.
     */
    ProteinSequence(string id, string seq) : identifier(id), data(seq) {
        for (char c : data) {
            if (VALID_CHARS.find(c) == string::npos)
                throw invalid_argument("Nieprawidłowy znak w sekwencji białka.");
        }
    }

    /**
     * Zwraca długość białka.
     */
    int length() const {
        return data.length();
    }

    /**
     * Format FASTA.
     */
    string toString() const {
        return ">" + identifier + "\n" + data;
    }

    /**
     * Mutuje aminokwas.
     */
    void mutate(int position, char value) {
        if (position < 0 || position >= data.length())
            throw out_of_range("Pozycja poza zakresem.");
        if (VALID_CHARS.find(value) == string::npos)
            throw invalid_argument("Nieprawidłowy znak mutacji.");
        data[position] = value;
    }

    /**
     * Znajduje motyw aminokwasowy.
     */
    int findMotif(const string& motif) const {
        size_t pos = data.find(motif);
        return pos != string::npos ? static_cast<int>(pos) : -1;
    }
};

const string ProteinSequence::VALID_CHARS = "ACDEFGHIKLMNPQRSTVWY"; // 20 standardowych aminokwasów

// Transkrypcja nici matrycowej DNA → RNA
RNASequence DNASequence::transcribe() const {
    string rna;
    for (char c : data) {
        switch (c) {
        case 'A': rna += 'U'; break;
        case 'T': rna += 'A'; break;
        case 'C': rna += 'G'; break;
        case 'G': rna += 'C'; break;
        }
    }
    return RNASequence(identifier + "_RNA", rna);
}

// Uproszczona translacja: każde 3 nukleotydy → 'X'
ProteinSequence RNASequence::transcribe() const {
    string protein;
    for (size_t i = 0; i + 2 < data.length(); i += 3) {
        protein += 'X'; // uproszczenie — normalnie kodon → aminokwas
    }
    return ProteinSequence(identifier + "_protein", protein);
}

/** Funkcja główna — demonstruje użycie klas */
int main() {
    try {
        DNASequence dna("seq1", "TACGGA");
        cout << dna.toString() << endl;
        cout << "Komplementarna nić: " << dna.complement() << endl;

        RNASequence rna = dna.transcribe();
        cout << rna.toString() << endl;

        ProteinSequence protein = rna.transcribe();
        cout << protein.toString() << endl;

        dna.mutate(0, 'A');
        cout << "Po mutacji DNA: " << dna.toString() << endl;

        cout << "Motyw 'CGG' w DNA na pozycji: " << dna.findMotif("CGG") << endl;

    }
    catch (const exception& e) {
        cerr << "Błąd: " << e.what() << endl;
    }

    return 0;
}
