#include <iostream>
#include <string>
#include <stdexcept>
#include <unordered_map>

using namespace std;

/** Klasa reprezentująca sekwencję RNA */
class RNASequence;

/** Klasa reprezentująca sekwencję białkową */
class ProteinSequence;

/** Klasa reprezentująca sekwencję DNA */
class DNASequence {
private:
    string identifier;
    string data;
    static const string VALID_CHARS;

public:
    /**
     * Konstruktor sprawdzający poprawność sekwencji.
     * @param id Identyfikator sekwencji
     * @param seq Sekwencja zasad
     */
    DNASequence(string id, string seq) : identifier(id), data(seq) {
        for (char c : data) {
            if (VALID_CHARS.find(c) == string::npos)
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
    string identifier;
    string data;
    static const string VALID_CHARS;

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

    ProteinSequence transcribe() const;
};

const string RNASequence::VALID_CHARS = "AUCG";

/** Klasa reprezentująca sekwencję białkową */
class ProteinSequence {
private:
    string identifier;
    string data;
    static const string VALID_CHARS;

public:
    ProteinSequence(string id, string seq) : identifier(id), data(seq) {
        for (char c : data) {
            if (VALID_CHARS.find(c) == string::npos)
                throw invalid_argument("Nieprawidłowy znak w sekwencji białka.");
        }
    }

    /**
     * Zwraca długość białka.
     */
    int length() const { return data.length(); }

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

    int findMotif(const string& motif) const {
        size_t pos = data.find(motif);
        return pos != string::npos ? static_cast<int>(pos) : -1;
    }
};

const string ProteinSequence::VALID_CHARS = "ACDEFGHIKLMNPQRSTVWY";

// Transkrypcja DNA na RNA
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

// Translacja RNA na białko
ProteinSequence RNASequence::transcribe() const {
    static const unordered_map<string, char> codonTable = {  // mapa zaproponowana przez ChatGPT
        {"UUU", 'F'}, {"UUC", 'F'}, {"UUA", 'L'}, {"UUG", 'L'},
        {"CUU", 'L'}, {"CUC", 'L'}, {"CUA", 'L'}, {"CUG", 'L'},
        {"AUU", 'I'}, {"AUC", 'I'}, {"AUA", 'I'}, {"AUG", 'M'},
        {"GUU", 'V'}, {"GUC", 'V'}, {"GUA", 'V'}, {"GUG", 'V'},
        {"UCU", 'S'}, {"UCC", 'S'}, {"UCA", 'S'}, {"UCG", 'S'},
        {"CCU", 'P'}, {"CCC", 'P'}, {"CCA", 'P'}, {"CCG", 'P'},
        {"ACU", 'T'}, {"ACC", 'T'}, {"ACA", 'T'}, {"ACG", 'T'},
        {"GCU", 'A'}, {"GCC", 'A'}, {"GCA", 'A'}, {"GCG", 'A'},
        {"UAU", 'Y'}, {"UAC", 'Y'}, {"UAA", '*'}, {"UAG", '*'},
        {"CAU", 'H'}, {"CAC", 'H'}, {"CAA", 'Q'}, {"CAG", 'Q'},
        {"AAU", 'N'}, {"AAC", 'N'}, {"AAA", 'K'}, {"AAG", 'K'},
        {"GAU", 'D'}, {"GAC", 'D'}, {"GAA", 'E'}, {"GAG", 'E'},
        {"UGU", 'C'}, {"UGC", 'C'}, {"UGA", '*'}, {"UGG", 'W'},
        {"CGU", 'R'}, {"CGC", 'R'}, {"CGA", 'R'}, {"CGG", 'R'},
        {"AGU", 'S'}, {"AGC", 'S'}, {"AGA", 'R'}, {"AGG", 'R'},
        {"GGU", 'G'}, {"GGC", 'G'}, {"GGA", 'G'}, {"GGG", 'G'}
    };

    string protein;
    for (size_t i = 0; i + 2 < data.length(); i += 3) {
        string codon = data.substr(i, 3);
        auto it = codonTable.find(codon);
        if (it != codonTable.end()) {
            if (it->second == '*') break;
            protein += it->second;
        }
        else {
            protein += 'X'; // Nieznany kodon
        }
    }
    return ProteinSequence(identifier + "_protein", protein);
}

// Testy
int main() {
    try {
        DNASequence dna("seq1", "TACGGCATTGAA");
        cout << dna.toString() << endl;
        cout << "Komplementarna nić: " << dna.complement() << endl;

        RNASequence rna = dna.transcribe();
        cout << rna.toString() << endl;

        ProteinSequence protein = rna.transcribe();
        cout << protein.toString() << endl;

        dna.mutate(0, 'A');
        cout << "Po mutacji DNA: " << dna.toString() << endl;

        cout << "Motyw 'GGC' w DNA na pozycji: " << dna.findMotif("GGC") << endl;

    }
    catch (const exception& e) {
        cerr << "Błąd: " << e.what() << endl;
    }

    return 0;
}
