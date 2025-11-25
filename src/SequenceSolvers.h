#pragma once
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <map>
#include <gmpxx.h>
#include <iostream>
#include <memory>
#include <iomanip>

// Use GMP for arbitrary precision arithmetic
using BigInt = mpz_class;
using BigRat = mpq_class;

// ============================================================================
// HELPERS
// ============================================================================

inline BigRat absRat(const BigRat& r) {
    return (r < 0) ? -r : r;
}

inline std::vector<BigInt> getFactSequence(size_t n) {
    std::vector<BigInt> facts;
    BigInt f = 1;
    for(size_t i=0; i<n; ++i) {
        if (i > 0) f *= static_cast<unsigned long>(i);
        facts.push_back(f);
    }
    return facts;
}

inline bool isPrime(const BigInt& n) {
    if (n <= 1) return false;
    return mpz_probab_prime_p(n.get_mpz_t(), 25) > 0;
}

inline BigInt getNthPrime(long n) {
    if (n < 1) return 2;
    BigInt p = 2;
    for (long i = 1; i < n; ++i) mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
    return p;
}

inline long getPrimeIndex(const BigInt& n) {
    if (!isPrime(n)) return -1;
    BigInt p = 2;
    long idx = 1;
    while (p < n) {
        mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
        idx++;
    }
    return (p == n) ? idx : -1;
}

// Solves Ax = b for x
inline std::vector<BigRat> solveSystem(std::vector<std::vector<BigRat>> mat, std::vector<BigRat> vec) {
    size_t n = mat.size();
    if (n == 0 || mat[0].size() != n || vec.size() != n) return {};

    // augment
    for (size_t i = 0; i < n; ++i) mat[i].push_back(vec[i]);

    // gaussian elimination
    for (size_t i = 0; i < n; ++i) {
        // pivot
        size_t p = i;
        while (p < n && mat[p][i] == 0) p++;
        if (p == n) return {}; // singular
        std::swap(mat[i], mat[p]);

        BigRat div = mat[i][i];
        for (size_t j = i; j <= n; ++j) mat[i][j] /= div;

        for (size_t k = 0; k < n; ++k) {
            if (k != i) {
                BigRat mul = mat[k][i];
                for (size_t j = i; j <= n; ++j) mat[k][j] -= mul * mat[i][j];
            }
        }
    }
    std::vector<BigRat> res;
    for (size_t i = 0; i < n; ++i) res.push_back(mat[i][n]);
    return res;
}

inline std::vector<BigRat> solveHomogeneousSystem(std::vector<std::vector<BigRat>>& mat) {
    if (mat.empty()) return {};
    size_t rows = mat.size();
    size_t cols = mat[0].size();
    size_t pivot_row = 0;
    std::vector<size_t> pivot_cols;
    
    for (size_t col = 0; col < cols && pivot_row < rows; ++col) {
        size_t sel = pivot_row;
        while (sel < rows && mat[sel][col] == 0) sel++;
        if (sel == rows) continue;
        std::swap(mat[sel], mat[pivot_row]);
        pivot_cols.push_back(col);
        BigRat pivot_val = mat[pivot_row][col];
        for (size_t j = col; j < cols; ++j) mat[pivot_row][j] /= pivot_val;
        for (size_t i = 0; i < rows; ++i) {
            if (i != pivot_row && mat[i][col] != 0) {
                BigRat factor = mat[i][col];
                for (size_t j = col; j < cols; ++j) mat[i][j] -= factor * mat[pivot_row][j];
            }
        }
        pivot_row++;
    }

    if (pivot_cols.size() < cols) {
        std::vector<BigRat> sol(cols, 0);
        int free_col = -1;
        for(int c = (int)cols-1; c >= 0; --c) {
             bool is_pivot = false;
             for(size_t p : pivot_cols) if((int)p == c) is_pivot = true;
             if(!is_pivot) { free_col = c; break; }
        }
        if(free_col == -1) return {};
        sol[free_col] = 1;
        for (int i = (int)pivot_row - 1; i >= 0; --i) {
            int p_col = -1;
            for(size_t j=0; j<cols; ++j) {
                if(mat[i][j] != 0) { p_col = j; break; }
            }
            if (p_col == -1 || p_col > free_col) continue;
            BigRat sum = 0;
            for (size_t j = p_col + 1; j < cols; ++j) sum += mat[i][j] * sol[j];
            sol[p_col] = -sum;
        }
        return sol;
    }
    return {};
}

struct SolverResult {
    BigRat nextValue;
    std::string formula;
    std::string formattedOutput;
    double confidence; 
    std::string unsimplifiedNum = "";
    std::string unsimplifiedDen = "";
    std::string explanation = "";
};

// ============================================================================
// RANKING LOGIC (TIER SYSTEM)
// ============================================================================

inline int getSolverTier(const std::string& name) {
    // Tier 1: Absolute Fundamentals
    if (name.find("Constant") != std::string::npos) return 1;
    if (name.find("Look-and-Say") != std::string::npos) return 1;
    if (name.find("n^n") != std::string::npos) return 1;
    if (name.find("Standard Math") != std::string::npos) return 1;
    if (name.find("Perfect Power") != std::string::npos) return 1;
    if (name.find("Geometric") != std::string::npos) return 1;
    
    // Only trust "Linear" if it's NOT a recurrence (i.e. Arithmetic Progression)
    if (name.find("Linear") != std::string::npos && name.find("Recurrence") == std::string::npos) return 1;
    if (name.find("Visual") != std::string::npos) return 1;
    if (name.find("Alternating") != std::string::npos) return 1;
    if (name.find("Indexed Primes") != std::string::npos) return 1;

    // Tier 2: Strong Mathematical Models
    if (name.find("Polynomial") != std::string::npos && name.find("Weighted") == std::string::npos) return 2;
    if (name.find("Coupled Step") != std::string::npos) return 2;
    if (name.find("x_n =") != std::string::npos) return 2; // Explicit formulas
    if (name.find("Simple Linear Recurrence") != std::string::npos) return 2; // Tribonacci/Fibonacci
    if (name.find("Variable Linear Recurrence") != std::string::npos) return 2;

    // Weighted patterns inheriting strength
    if (name.find("Weighted") != std::string::npos) {
        if (name.find("Perfect Power") != std::string::npos) return 2;
        if (name.find("Geometric") != std::string::npos) return 2;
        if (name.find("Simple Linear Recurrence") != std::string::npos) return 2;
        // Note: Generic "Weighted: Linear Recurrence" falls to Tier 6 to avoid overfitting
        return 6;
    }

    // Tier 3: Structural/Pattern
    if (name.find("Interleaved") != std::string::npos && name.find("Differences") == std::string::npos) return 3;
    if (name.find("Cyclic") != std::string::npos) return 3;

    // Tier 4: Recursive/Derivatives
    if (name.find("Recurrence") != std::string::npos) return 4;
    if (name.find("Holonomic") != std::string::npos) return 4;
    if (name.find("Prime Differences") != std::string::npos) return 4;
    if (name.find("Interleaved Differences") != std::string::npos) return 4; // Specifically demoted below Cyclic

    // Tier 6: Complex/Composite
    if (name.find("Fraction Pattern") != std::string::npos) return 6;
    
    return 7; // Unknown
}

inline bool compareResults(const SolverResult& a, const SolverResult& b) {
    // 1. Significant confidence difference overrides everything
    if (std::abs(a.confidence - b.confidence) > 0.05) {
        return a.confidence > b.confidence;
    }

    // 2. Tier System
    int tierA = getSolverTier(a.formula);
    int tierB = getSolverTier(b.formula);
    if (tierA != tierB) {
        return tierA < tierB;
    }

    // 3. Simplicity (Shorter formula / description)
    if (a.formula.length() != b.formula.length()) {
        return a.formula.length() < b.formula.length();
    }

    // 4. Fallback to raw confidence
    return a.confidence > b.confidence;
}

class SequenceSolver {
public:
    virtual ~SequenceSolver() = default;
    virtual std::string getName() const = 0;
    virtual std::vector<SolverResult> solve(const std::vector<BigRat>& series) = 0;
};

// ============================================================================
// SOLVERS
// ============================================================================

class PerfectPowerSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Perfect Powers"; }
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        if (series.size() < 3) return {};
        
        for (int p = 2; p <= 5; ++p) {
            std::vector<BigRat> roots;
            bool allPowers = true;
            for (const auto& val : series) {
                if (val < 0) { allPowers = false; break; } 
                BigInt num = val.get_num();
                BigInt den = val.get_den();
                BigInt rNum, rDen;
                if (mpz_root(rNum.get_mpz_t(), num.get_mpz_t(), p) == 0) { allPowers = false; break; }
                if (mpz_root(rDen.get_mpz_t(), den.get_mpz_t(), p) == 0) { allPowers = false; break; }
                roots.push_back(BigRat(rNum, rDen));
            }
            
            if (allPowers) {
                BigRat diff = roots[1] - roots[0];
                bool linear = true;
                for(size_t i=1; i<roots.size(); ++i) {
                    if (roots[i] - roots[i-1] != diff) { linear = false; break; }
                }
                if (linear) {
                    BigRat nextRoot = roots.back() + diff;
                    BigInt nextValNum, nextValDen;
                    mpz_pow_ui(nextValNum.get_mpz_t(), nextRoot.get_num().get_mpz_t(), p);
                    mpz_pow_ui(nextValDen.get_mpz_t(), nextRoot.get_den().get_mpz_t(), p);
                    return {{BigRat(nextValNum, nextValDen), "Perfect Power (n^" + std::to_string(p) + ")", "", 1.0, "", "", "Sequence is a perfect power of a linear progression."}};
                }
            }
        }
        return {};
    }
};

class PolynomialSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Polynomial"; }
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        if (series.size() < 3) return {};
        bool allSame = true;
        for(size_t i=1; i<series.size(); ++i) if(series[i] != series[0]) allSame = false;
        if(allSame) return {{series[0], "Constant", "", 1.0, "", "", "All terms are equal."}};

        std::vector<std::vector<BigRat>> table;
        table.push_back(series);
        int depth = 0;
        size_t max_depth = series.size() - 1;

        while (table.back().size() > 1 && depth < (int)max_depth) {
            const auto& row = table.back();
            bool constantRow = true;
            for (size_t i = 1; i < row.size(); ++i) {
                if (row[i] != row[0]) { constantRow = false; break; }
            }
            if (constantRow) {
                BigRat nextVal = row[0]; 
                for (int i = depth - 1; i >= 0; --i) nextVal = nextVal + table[i].back();
                double conf = 0.97;
                if (depth > 2 && (size_t)depth >= series.size() - 2) conf = 0.85;
                std::string name = (depth == 1) ? "Linear" : "Polynomial (Degree " + std::to_string(depth) + ")";
                return {{nextVal, name, "", conf, "", "", "Polynomial sequence."}};
            }
            std::vector<BigRat> nextRow;
            for(size_t i = 0; i < row.size() - 1; ++i) nextRow.push_back(row[i+1] - row[i]);
            table.push_back(nextRow);
            depth++;
        }
        return {};
    }
};

class VisualSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Visual/String Pattern"; }
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        size_t N = series.size();
        if (N < 3) return {};
        if (N >= 4) {
            BigRat d1 = series[1] - series[0];
            BigRat d2 = series[2] - series[1];
            if ((d1 == 0 && d2 != 0) || (d1 != 0 && d2 == 0)) {
                bool altZero = true;
                BigRat step = (d1 == 0) ? d2 : d1;
                for(size_t i=1; i<N-1; ++i) {
                    BigRat d = series[i+1] - series[i];
                    if (d == 0 && series[i] - series[i-1] == 0) altZero = false;
                    if (d != 0 && d != step) altZero = false;
                }
                if (altZero) {
                    BigRat lastDiff = series.back() - series[N-2];
                    BigRat nextDiff = (lastDiff == 0) ? step : 0;
                    return {{series.back() + nextDiff, "Visual: Step-Repeat", "", 1.0, "", "", "Values repeat once then step."}};
                }
            }
        }
        BigRat diff = series[1] - series[0];
        bool isLin = true;
        for(size_t i=1; i<N; ++i) if (series[i] - series[i-1] != diff) isLin = false;
        if (isLin) return {{series.back() + diff, "Linear", "", 1.0, "", "", "Simple linear progression."}};
        
        if (N >= 3) {
            bool alt = true;
            for(size_t i=2; i<N; ++i) if (series[i] != series[i-2]) alt = false;
            if (alt) {
                return {{series[N-2], "Alternating (Period 2)", "", 0.99, "", "", "Values alternate A, B, A..."}};
            }
        }
        return {};
    }
};

class PowerSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Geometric"; }
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        if (series.size() < 3) return {};
        
        for(size_t skip=0; skip <= 2 && skip + 3 <= series.size(); ++skip) {
            BigRat s0 = series[skip];
            BigRat s1 = series[skip+1];
            if (s0 == 0) continue;
            
            BigRat ratio = s1 / s0;
            if (ratio == 1) continue; 
            
            bool match = true;
            for(size_t i=skip; i<series.size()-1; ++i) {
                if (series[i] == 0 || series[i+1] / series[i] != ratio) { match = false; break; }
            }
            
            if (match) {
                double conf = (skip == 0) ? 1.0 : 0.99; 
                std::string name = (skip == 0) ? "Geometric" : "Geometric (with offset)";
                return {{series.back() * ratio, name, "", conf, "", "", "Geometric progression."}};
            }
        }
        return {};
    }
};

class StandardMathSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Standard Math"; }
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        if (series.size() < 3) return {};
        int offset = (series[0] == 1 && series[1] == 4) ? 1 : 0;
        bool isPower = true;
        for(size_t i=0; i<series.size(); ++i) {
            unsigned long idx = i + offset;
            BigInt p; mpz_ui_pow_ui(p.get_mpz_t(), idx, idx);
            if (series[i].get_num() != p) { isPower = false; break; }
        }
        if (isPower) {
            unsigned long nextIdx = series.size() + offset;
            BigInt res; mpz_ui_pow_ui(res.get_mpz_t(), nextIdx, nextIdx);
            return {{BigRat(res), "n^n Sequence", "", 1.0, "", "", "Each term is n raised to the power of n."}};
        }
        return {};
    }
};

class DigitSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Digit Patterns"; }
    std::string lookAndSayNext(std::string s) {
        if (s.empty()) return "";
        std::stringstream ss;
        for (size_t i = 0; i < s.length(); i++) {
            int count = 1;
            while (i + 1 < s.length() && s[i] == s[i+1]) { i++; count++; }
            ss << count << s[i];
        }
        return ss.str();
    }
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        bool isLAS = true;
        for(size_t i=0; i<series.size()-1; ++i) {
             if (series[i].get_den() != 1) { isLAS = false; break; }
             if (lookAndSayNext(series[i].get_num().get_str()) != series[i+1].get_num().get_str()) { isLAS = false; break; }
        }
        if (isLAS && series.size() > 2) {
            return {{BigRat(BigInt(lookAndSayNext(series.back().get_num().get_str()))), "Look-and-Say Sequence", "", 1.0, "", "", "Each term describes digits of previous term."}};
        }
        return {};
    }
};

class BerlekampMasseySolver : public SequenceSolver {
public:
    std::string getName() const override { return "Linear Recurrence"; }
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        if (series.size() < 4) return {};
        int N = series.size();
        
        bool inputsAreIntegers = true;
        for(const auto& s : series) if(s.get_den() != 1) inputsAreIntegers = false;

        std::vector<BigRat> C = {1}, B = {1};
        int L = 0, m = -1;
        BigRat b(1);

        for (int n = 0; n < N; ++n) {
            BigRat d(0);
            for (size_t i = 0; i < C.size(); ++i) 
                if (n >= (int)i) d += C[i] * series[n - i];
            if (d == 0) continue;
            std::vector<BigRat> T = C;
            BigRat coef = d / b;
            int shift = n - m;
            std::vector<BigRat> scaledB(shift, 0);
            for(const auto& val : B) scaledB.push_back(val * coef);
            if (scaledB.size() > T.size()) T.resize(scaledB.size(), 0);
            for(size_t i=0; i<scaledB.size(); ++i) T[i] -= scaledB[i];
            if (2 * L <= n) { L = n + 1 - L; m = n; B = C; b = d; }
            C = T;
        }
        if (L == 0 || 2 * L > N + 1) return {}; 
        BigRat nextVal(0);
        BigRat c0 = C[0];
        
        bool allSmallInts = true;
        for (size_t i = 1; i < C.size(); ++i) {
            if (series.size() >= i) nextVal += ((C[i] * -1) / c0) * series[series.size() - i];
            BigRat coef = (C[i] * -1) / c0;
            if (coef.get_den() != 1 || abs(coef.get_num().get_si()) > 5) allSmallInts = false;
        }

        double ratio = (double)N / (double)L;
        double conf = 0.5;

        std::string name = "Linear Recurrence";
        if (allSmallInts) name = "Simple Linear Recurrence";

        if (ratio >= 3.0) conf = 0.96; 
        else if (ratio >= 2.0) {
            if (allSmallInts) conf = 0.99; 
            else conf = 0.85; 
        } else {
            conf = 0.4; 
        }

        if (L <= 2 && conf < 0.8) conf += 0.2;

        if (inputsAreIntegers && nextVal.get_den() > 1) {
            conf = 0.1; 
        }

        return {{nextVal, name, "", conf, "", "", "Linear recurrence found via Berlekamp-Massey."}}; 
    }
};

class RecursiveSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Recursive Non-Linear"; }
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        std::vector<SolverResult> results;
        size_t N = series.size();
        if (N < 3) return {};
        
        bool sq = true, sylv = true;
        for(size_t i=1; i<N; ++i) {
             if(series[i] != series[i-1]*series[i-1]) sq=false;
             if(series[i] != (series[i-1]*series[i-1]) + series[i-1]) sylv=false;
        }
        if (sq) results.push_back({series.back()*series.back(), "x_n = x_{n-1}^2", "", 1.0, "", "", ""});
        if (sylv) results.push_back({(series.back()*series.back())+series.back(), "x_n = x_{n-1}^2 + x_{n-1}", "", 1.0, "", "", ""});
        
        if (N >= 4) {
            bool lag1 = true, lag2 = true;
            for(size_t i=2; i<N; ++i) {
                if (series[i] != (series[i-1] * series[i-2]) - series[i-2]) lag1 = false;
                if (series[i] != (series[i-1] * series[i-2]) - series[i-1]) lag2 = false;
            }
            if (lag1) results.push_back({(series.back() * series[N-2]) - series[N-2], "x_n = x_{n-1}x_{n-2} - x_{n-2}", "", 0.99, "", "", ""});
            if (lag2) results.push_back({(series.back() * series[N-2]) - series.back(), "x_n = x_{n-1}x_{n-2} - x_{n-1}", "", 0.99, "", "", ""});
        }
        if (N >= 3) {
             BigRat num = series[2] - series[1];
             BigRat den = series[1] - series[0];
             if (den != 0) {
                 BigRat A = num / den;
                 BigRat B = series[1] - A * series[0];
                 bool match = true;
                 for(size_t i=1; i<N; ++i) {
                     if (series[i] != A * series[i-1] + B) { match = false; break; }
                 }
                 if (match) {
                     double conf = 0.98;
                     if (N == 3) {
                         if (A.get_den() == 1 && B.get_den() == 1 && absRat(A) < 20 && absRat(B) < 20) conf = 0.90;
                         else conf = 0.5;
                     }
                     if (conf > 0.6)
                        results.push_back({A * series.back() + B, "Affine Recurrence", "", conf, "", "", "x_n = A*x_{n-1} + B"});
                 }
             }
        }
        if (N >= 5) {
             std::vector<std::vector<BigRat>> mat;
             std::vector<BigRat> vec;
             for(size_t i=1; i<=4; ++i) { 
                 if (i >= N) break;
                 std::vector<BigRat> row;
                 BigRat idx = static_cast<unsigned long>(i); 
                 row.push_back(idx * series[i-1]);
                 row.push_back(series[i-1]);
                 row.push_back(idx);
                 row.push_back(1);
                 mat.push_back(row);
                 vec.push_back(series[i]);
             }
             if (mat.size() == 4) {
                 std::vector<BigRat> sol = solveSystem(mat, vec);
                 if (!sol.empty()) {
                     BigRat A = sol[0], B = sol[1], C = sol[2], D = sol[3];
                     bool match = true;
                     for(size_t i=1; i<N; ++i) {
                         BigRat idx = static_cast<unsigned long>(i);
                         BigRat pred = series[i-1] * (A*idx + B) + (C*idx + D);
                         if (pred != series[i]) { match = false; break; }
                     }
                     if (match) {
                         BigRat nextIdx = static_cast<unsigned long>(N);
                         BigRat nextVal = series.back() * (A*nextIdx + B) + (C*nextIdx + D);
                         double conf = (N >= 7) ? 0.995 : 0.95; 
                         results.push_back({nextVal, "Variable Linear Recurrence", "", conf, "", "", "x_n = x_{n-1}(An+B) + (Cn+D)"});
                     }
                 }
             }
        }
        if (N >= 3) {
            bool isAP = true;
            BigRat diff = series[1] - series[0];
            for(size_t i=1; i<N; ++i) if(series[i] - series[i-1] != diff) isAP = false;
            if (!isAP) {
                bool fibSum = true, fibSub1 = true, fibSub2 = true;
                for(size_t i=2; i<N; ++i) {
                    if (series[i] != series[i-1] + series[i-2]) fibSum = false;
                    if (series[i] != series[i-1] - series[i-2]) fibSub1 = false;
                    if (series[i] != series[i-2] - series[i-1]) fibSub2 = false;
                }
                if (fibSum) results.push_back({series.back() + series[N-2], "x_n = x_{n-1} + x_{n-2}", "", 0.99, "", "", ""});
                if (fibSub1) results.push_back({series.back() - series[N-2], "x_n = x_{n-1} - x_{n-2}", "", 0.99, "", "", ""});
                if (fibSub2) results.push_back({series[N-2] - series.back(), "x_n = x_{n-2} - x_{n-1}", "", 0.99, "", "", ""});
            }
        }
        if (N >= 3) {
            bool varInc = true;
            BigRat k = (series[1] - series[0]) - 1;
            for(size_t i=1; i<N; ++i) {
                if (series[i] != series[i-1] + BigRat(static_cast<unsigned long>(i)) + k) { varInc = false; break; }
            }
            if (varInc) results.push_back({series.back() + BigRat(static_cast<unsigned long>(N)) + k, "x_n = x_{n-1} + n + K", "", 1.0, "", "", ""});
        }
        if (N >= 4) {
            bool coupled = true;
            for(size_t i=2; i<N; ++i) {
                size_t offset = 2 + (i % 2); 
                if (i < offset) continue; 
                if (series[i] != series[i-1] + series[i-offset]) { coupled = false; break; }
            }
            if (coupled) {
                size_t offset = 2 + (N % 2);
                double conf = (N >= 6) ? 0.99 : 0.96; 
                results.push_back({series.back() + series[N-offset], "Coupled Step Recurrence", "", conf, "", "", ""});
            }
        }
        return results;
    }
};

class HolonomicSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Holonomic"; }
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        int N = series.size();
        if (N < 5) return {}; 
        std::vector<std::pair<int, int>> configs = {
            {1, 1}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {3, 0}, {3, 1}
        };
        for (auto [R, D] : configs) {
            int num_coeffs = (R + 1) * (D + 1);
            if (N - R < num_coeffs) continue;
            std::vector<std::vector<BigRat>> matrix;
            for (int k = 0; k < N - R; ++k) {
                std::vector<BigRat> row;
                BigInt n_val = k; 
                for (int r = 0; r <= R; ++r) {
                    BigRat a_term = series[k + r];
                    BigInt n_pow = 1;
                    for (int d = 0; d <= D; ++d) { row.push_back(a_term * BigRat(n_pow)); n_pow *= n_val; }
                }
                matrix.push_back(row);
            }
            auto coeffs = solveHomogeneousSystem(matrix);
            if (!coeffs.empty()) {
                BigInt n_target = N - R;
                BigRat rhs = 0, lhs_coeff = 0; 
                int coeff_idx = 0;
                for (int r = 0; r <= R; ++r) {
                    BigRat p_r_val = 0; BigInt n_pow = 1;
                    for (int d = 0; d <= D; ++d) { p_r_val += coeffs[coeff_idx++] * BigRat(n_pow); n_pow *= n_target; }
                    if (r == R) lhs_coeff = p_r_val;
                    else rhs -= p_r_val * series[n_target.get_si() + r];
                }
                if (lhs_coeff != 0) return {{rhs / lhs_coeff, "Holonomic (Order " + std::to_string(R) + ")", "", 0.90, "", "", ""}};
            }
        }
        return {};
    }
};

class PrimeSolver : public SequenceSolver {
    std::unique_ptr<PolynomialSolver> poly;
public:
    PrimeSolver() { poly = std::make_unique<PolynomialSolver>(); }
    std::string getName() const override { return "Primes"; }
    
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        if (series.size() < 3) return {};
        for(auto& s : series) if(s.get_den() != 1) return {};

        std::vector<BigRat> indices;
        bool allPrime = true;
        for(const auto& s : series) {
            BigInt val = s.get_num();
            if (val <= 0) { allPrime = false; break; }
            long idx = getPrimeIndex(val);
            if (idx == -1) { allPrime = false; break; }
            indices.push_back(BigRat(idx));
        }

        if (allPrime) {
            auto results = poly->solve(indices);
            if (!results.empty() && results[0].confidence > 0.5) {
                BigInt nextIdx = results[0].nextValue.get_num();
                if (nextIdx > 0) return {{BigRat(getNthPrime(nextIdx.get_si())), "Indexed Primes", "", 0.995, "", "", "Sequence of indexed primes."}};
            }
        }
        
        std::vector<BigRat> diffs;
        bool diffsArePrime = true;
        bool decreasing = true; 
        for(size_t i=0; i<series.size()-1; ++i) {
            BigRat d = series[i+1] - series[i];
            if (d > 0) decreasing = false;
            BigInt num = d.get_num();
            if (num < 0) num = -num; 
            if (!isPrime(num)) { diffsArePrime = false; break; }
            diffs.push_back(d);
        }
        
        if (diffsArePrime && !diffs.empty()) {
            std::vector<BigRat> absDiffs;
            for(auto& d : diffs) absDiffs.push_back(absRat(d));
            
            std::vector<BigRat> diffIndices;
            for(auto& d : absDiffs) diffIndices.push_back(BigRat(getPrimeIndex(d.get_num())));
            
            auto results = poly->solve(diffIndices);
            if(!results.empty() && results[0].confidence > 0.5) {
                BigInt nextIdx = results[0].nextValue.get_num();
                BigInt nextPrimeVal = getNthPrime(nextIdx.get_si());
                BigRat nextDiff = BigRat(nextPrimeVal);
                if (decreasing && diffs.back() < 0) nextDiff = -nextDiff;
                else if (!decreasing && diffs.back() > 0) nextDiff = nextDiff;
                else if (diffs.back() < 0) nextDiff = -nextDiff;
                
                double conf = (series.size() <= 4) ? 0.85 : 0.99;
                return {{series.back() + nextDiff, "Prime Differences", "", conf, "", "", "Differences are primes."}}; 
            }
        }
        return {};
    }
};

class CyclicOpSolver : public SequenceSolver {
    struct Op { char type; BigRat val; }; 
    
    double getComplexity(const BigRat& r) {
        BigInt n = abs(r.get_num());
        BigInt d = r.get_den();
        if (n > 1000000 || d > 1000000) return 1e9;
        return n.get_d() + d.get_d();
    }

public:
    std::string getName() const override { return "Cyclic Operations"; }
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        size_t N = series.size();
        if (N < 3) return {};

        std::vector<SolverResult> results;

        for (size_t L = 1; L <= N/2 + 1; ++L) {
            std::vector<Op> cycleOps(L);
            bool cycleValid = true;
            bool allSimple = true; 

            for (size_t step = 0; step < L; ++step) {
                std::vector<std::pair<BigRat, BigRat>> pairs;
                for (size_t i = step; i < N - 1; i += L) {
                    pairs.push_back({series[i], series[i+1]});
                }
                if (pairs.empty()) continue; 

                bool foundOp = false;
                
                // 1. Check Square
                bool isSq = true;
                for (auto& p : pairs) {
                    if (p.first == 0 || p.second != p.first * p.first) isSq = false;
                }
                if (isSq) {
                    cycleOps[step] = {'^', 2};
                    foundOp = true;
                } 
                else {
                    // 2. Check Multiply
                    bool isMul = true;
                    BigRat mulVal;
                    if (pairs[0].first == 0) {
                        if (pairs[0].second == 0) mulVal = 0;
                        else isMul = false;
                    } else {
                        mulVal = pairs[0].second / pairs[0].first;
                    }

                    if (isMul) {
                        for(auto& p : pairs) {
                            if (p.first == 0) {
                                if (p.second != 0) isMul = false;
                            } else if (p.second / p.first != mulVal) {
                                isMul = false;
                            }
                        }
                    }

                    // 3. Check Add
                    bool isAdd = true;
                    BigRat addVal = pairs[0].second - pairs[0].first;
                    for (auto& p : pairs) {
                        if (p.second - p.first != addVal) isAdd = false;
                    }

                    if (isMul && isAdd) {
                        if (getComplexity(mulVal) <= getComplexity(addVal)) {
                            cycleOps[step] = {'*', mulVal};
                            if (mulVal.get_den() > 10) allSimple = false;
                        } else {
                            cycleOps[step] = {'+', addVal};
                            if (addVal.get_den() > 10) allSimple = false;
                        }
                        foundOp = true;
                    } else if (isMul) {
                        cycleOps[step] = {'*', mulVal};
                        // Allow simple fractions like 0.5 (1/2), 2.5 (5/2)
                        if (mulVal.get_den() > 10) allSimple = false;
                        foundOp = true;
                    } else if (isAdd) {
                        cycleOps[step] = {'+', addVal};
                        if (addVal.get_den() > 10) allSimple = false;
                        foundOp = true;
                    }
                }

                if (!foundOp) {
                    cycleValid = false; 
                    break;
                }
            }

            if (cycleValid) {
                size_t stepIdx = (N - 1) % L;
                Op op = cycleOps[stepIdx];
                BigRat nextVal;
                if (op.type == '+') nextVal = series.back() + op.val;
                else if (op.type == '*') nextVal = series.back() * op.val;
                else if (op.type == '^') nextVal = series.back() * series.back();

                double confidence = 0.5;
                if (N >= 3 * L) confidence = 0.99;
                else if (N >= 2 * L) confidence = 0.97;
                else {
                    // Relaxed Partial Cycle: e.g. Length 5, Cycle 2 -> 2.5 cycles
                    if (allSimple && N >= L + 2) confidence = 0.96;
                    else confidence = 0.6; 
                }

                results.push_back({nextVal, "Cyclic Operations (Length " + std::to_string(L) + ")", "", confidence, "", "", "Operations repeat every L steps."});
            }
        }
        return results;
    }
};

class InterleavedSolver : public SequenceSolver {
    std::unique_ptr<PolynomialSolver> poly;
    std::unique_ptr<PowerSolver> power;
    std::unique_ptr<VisualSolver> visual;
    std::unique_ptr<RecursiveSolver> recursive; 
public:
    InterleavedSolver() { 
        poly = std::make_unique<PolynomialSolver>(); 
        power = std::make_unique<PowerSolver>();
        visual = std::make_unique<VisualSolver>();
        recursive = std::make_unique<RecursiveSolver>();
    }
    std::string getName() const override { return "Interleaved Sequence"; }
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        if (series.size() < 4) return {};
        std::vector<BigRat> even, odd;
        for(size_t i=0; i<series.size(); ++i) {
            if (i % 2 == 0) even.push_back(series[i]);
            else odd.push_back(series[i]);
        }
        auto solveSub = [&](const std::vector<BigRat>& s) -> SolverResult {
            auto r1 = poly->solve(s);
            if (!r1.empty() && r1[0].confidence > 0.5) return r1[0];
            auto r2 = power->solve(s);
            if (!r2.empty() && r2[0].confidence > 0.5) return r2[0];
            auto r3 = visual->solve(s);
            if (!r3.empty() && r3[0].confidence > 0.5) return r3[0];
            auto r4 = recursive->solve(s); 
            if (!r4.empty() && r4[0].confidence > 0.5) return r4[0];
            return {0, "", "", 0.0};
        };
        auto bestEven = solveSub(even);
        auto bestOdd = solveSub(odd);
        if (bestEven.confidence >= 0.5 && bestOdd.confidence >= 0.5) {
            double avgConf = 0.96; 
            
            BigRat predicted;
            if (series.size() % 2 == 0) predicted = bestEven.nextValue; 
            else predicted = bestOdd.nextValue;
            return {{predicted, "Interleaved Sequence", "", avgConf, "", "", "Two interleaved sequences found."}};
        }
        return {};
    }
};

class DifferenceSolver : public SequenceSolver {
    std::unique_ptr<InterleavedSolver> inter;
    std::unique_ptr<PolynomialSolver> poly;
public:
    DifferenceSolver() {
        inter = std::make_unique<InterleavedSolver>();
        poly = std::make_unique<PolynomialSolver>();
    }
    std::string getName() const override { return "Interleaved Differences"; }
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        if (series.size() < 5) return {};
        std::vector<BigRat> diffs;
        for(size_t i=0; i<series.size()-1; ++i) diffs.push_back(series[i+1] - series[i]);
        
        auto res = inter->solve(diffs);
        if (!res.empty() && res[0].confidence > 0.8) {
            return {{series.back() + res[0].nextValue, "Interleaved Differences", "", res[0].confidence * 0.95, "", "", "Differences form an interleaved sequence."}};
        }
        auto pRes = poly->solve(diffs);
        if (!pRes.empty() && pRes[0].confidence > 0.8) {
             return {{series.back() + pRes[0].nextValue, "Polynomial Differences", "", pRes[0].confidence * 0.95, "", "", ""}};
        }
        return {};
    }
};

class SimpleRobustSolver : public SequenceSolver {
    std::vector<std::unique_ptr<SequenceSolver>> solvers;
public:
    SimpleRobustSolver() {
        solvers.push_back(std::make_unique<PolynomialSolver>());
        solvers.push_back(std::make_unique<PerfectPowerSolver>()); 
        solvers.push_back(std::make_unique<PowerSolver>());
        solvers.push_back(std::make_unique<VisualSolver>());
        solvers.push_back(std::make_unique<StandardMathSolver>());
        solvers.push_back(std::make_unique<DigitSolver>());
        solvers.push_back(std::make_unique<BerlekampMasseySolver>());
        solvers.push_back(std::make_unique<RecursiveSolver>());
        solvers.push_back(std::make_unique<HolonomicSolver>());
        solvers.push_back(std::make_unique<PrimeSolver>());
    }
    std::string getName() const override { return "Simple Robust Solver"; }
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        std::vector<SolverResult> all;
        for(auto& s : solvers) {
            auto res = s->solve(series);
            if(!res.empty()) all.insert(all.end(), res.begin(), res.end());
        }
        std::sort(all.begin(), all.end(), compareResults);
        return all;
    }
};

class WeightedSolver : public SequenceSolver {
    std::unique_ptr<SequenceSolver> sub;
public:
    WeightedSolver() {
        sub = std::make_unique<SimpleRobustSolver>(); 
    }
    std::string getName() const override { return "Weighted Pattern"; }
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        if (series.size() < 3) return {};
        std::vector<SolverResult> allResults;

        // Helper
        auto process = [&](std::vector<BigRat>& wSeries, std::string name, int p) {
             auto res = sub->solve(wSeries);
             if (!res.empty() && res[0].confidence > 0.8) {
                 double finalConf = res[0].confidence;
                 // Penalize weighted fits on short sequences
                 if (series.size() < 6) finalConf = std::min(finalConf, 0.80);
                 else finalConf = std::min(0.90, finalConf);

                 if (p == 1 && series.size() >= 4) finalConf += 0.02;

                 // PREPEND "Weighted: " so Tier logic can see it
                 std::string finalName = "Weighted: " + res[0].formula;
                 return std::make_pair(true, SolverResult{res[0].nextValue, finalName, "", finalConf, "", "", ""});
             }
             return std::make_pair(false, SolverResult{0,"","",0.0});
        };

        // 1. Standard n^p weights
        for (int p = 1; p <= 3; ++p) {
             std::vector<BigRat> wSeries;
             bool possible = true;
             for(size_t i=0; i<series.size(); ++i) {
                 BigInt n = static_cast<unsigned long>(i + 1);
                 BigInt weight; mpz_ui_pow_ui(weight.get_mpz_t(), n.get_ui(), p);
                 if (weight == 0) { possible = false; break; }
                 wSeries.push_back(series[i] * BigRat(weight));
             }
             if (!possible) continue;
             
             auto r = process(wSeries, "n^" + std::to_string(p), p);
             if(r.first) {
                 BigInt nextN = static_cast<unsigned long>(series.size() + 1);
                 BigInt w; mpz_ui_pow_ui(w.get_mpz_t(), nextN.get_ui(), p);
                 r.second.nextValue = r.second.nextValue / BigRat(w);
                 allResults.push_back(r.second);
             }
        }

        // 2. Fibonacci Weights
        for(int offset = -2; offset <= 2; ++offset) {
             std::vector<BigRat> wSeries;
             bool possible = true;
             BigInt f1 = 1, f2 = 1; 
             std::vector<BigInt> fibs;
             BigInt a = 1, b = 1;
             for(size_t i=0; i<=series.size()+5; ++i) {
                 fibs.push_back(a);
                 BigInt next = a + b;
                 a = b; b = next;
             }
             for(size_t i=0; i<series.size(); ++i) {
                 long idx = (long)i + offset;
                 if (idx < 0 || idx >= (long)fibs.size()) { possible = false; break; }
                 if (fibs[idx] == 0) { possible = false; break; }
                 wSeries.push_back(series[i] * BigRat(fibs[idx]));
             }
             if (!possible) continue;

             auto res = sub->solve(wSeries);
             if (!res.empty() && res[0].confidence > 0.8) {
                 long nextIdx = (long)series.size() + offset;
                 BigInt nextFib = fibs[nextIdx];
                 double conf = (series.size() < 6) ? 0.80 : std::min(0.90, res[0].confidence);
                 std::string finalName = "Weighted: " + res[0].formula;
                 allResults.push_back({res[0].nextValue / BigRat(nextFib), finalName, "", conf, "", "", ""});
             }
        }

        // 3. Factorial Weights
        {
            std::vector<BigRat> wSeries;
            auto facts = getFactSequence(series.size() + 2);
            bool possible = true;
            for(size_t i=0; i<series.size(); ++i) {
                if (facts[i+1] == 0) { possible = false; break; }
                wSeries.push_back(series[i] * BigRat(facts[i+1]));
            }
            if (possible) {
                auto res = sub->solve(wSeries);
                if (!res.empty() && res[0].confidence > 0.8) {
                    double conf = (series.size() < 6) ? 0.80 : std::min(0.90, res[0].confidence);
                    std::string finalName = "Weighted: " + res[0].formula;
                    allResults.push_back({res[0].nextValue / BigRat(facts[series.size()+1]), finalName, "", conf, "", "", ""});
                }
            }
        }

        if (allResults.empty()) return {};
        std::sort(allResults.begin(), allResults.end(), compareResults);
        return {allResults[0]};
    }
};

class RobustSolverWrapper : public SequenceSolver {
    std::vector<std::unique_ptr<SequenceSolver>> solvers;
public:
    RobustSolverWrapper() {
        solvers.push_back(std::make_unique<PolynomialSolver>());
        solvers.push_back(std::make_unique<PerfectPowerSolver>());
        solvers.push_back(std::make_unique<PowerSolver>());
        solvers.push_back(std::make_unique<VisualSolver>());
        solvers.push_back(std::make_unique<StandardMathSolver>());
        solvers.push_back(std::make_unique<DigitSolver>());
        solvers.push_back(std::make_unique<BerlekampMasseySolver>());
        solvers.push_back(std::make_unique<RecursiveSolver>());
        solvers.push_back(std::make_unique<HolonomicSolver>());
        solvers.push_back(std::make_unique<PrimeSolver>());
        solvers.push_back(std::make_unique<CyclicOpSolver>());
        solvers.push_back(std::make_unique<DifferenceSolver>());
    }
    std::string getName() const override { return "Robust Wrapper"; }
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        std::vector<SolverResult> all;
        for(auto& s : solvers) {
            auto res = s->solve(series);
            if(!res.empty()) all.insert(all.end(), res.begin(), res.end());
        }
        std::sort(all.begin(), all.end(), compareResults);
        return all;
    }
};

class FractionSolver : public SequenceSolver {
    std::unique_ptr<RobustSolverWrapper> sub;
    std::unique_ptr<VisualSolver> visual;
    std::unique_ptr<InterleavedSolver> inter;
    std::unique_ptr<WeightedSolver> weight;
public:
    FractionSolver() {
        sub = std::make_unique<RobustSolverWrapper>();
        visual = std::make_unique<VisualSolver>();
        inter = std::make_unique<InterleavedSolver>();
        weight = std::make_unique<WeightedSolver>();
    }
    std::string getName() const override { return "Fraction Pattern"; }
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        bool hasFrac = false;
        std::vector<BigRat> nums, dens;
        for(const auto& s : series) { 
            if(s.get_den() != 1) hasFrac = true;
            nums.push_back(BigRat(s.get_num())); 
            dens.push_back(BigRat(s.get_den())); 
        }
        if (!hasFrac) return {};

        std::vector<SolverResult> results;
        auto vNum = visual->solve(nums);
        auto vDen = visual->solve(dens);
        if (!vNum.empty() && !vDen.empty()) {
            BigRat nextVal = vNum[0].nextValue / vDen[0].nextValue;
            double conf = std::min(vNum[0].confidence, vDen[0].confidence);
            SolverResult res = {nextVal, "Fraction Pattern", "", conf, "", "", "Numerator and Denominator follow visual patterns."};
            res.unsimplifiedNum = vNum[0].nextValue.get_str();
            res.unsimplifiedDen = vDen[0].nextValue.get_str();
            results.push_back(res);
        }
        auto solvePart = [&](const std::vector<BigRat>& p) {
            auto r = sub->solve(p);
            if (!r.empty() && r[0].confidence > 0.5) return r;
            return inter->solve(p); 
        };
        auto rNum = solvePart(nums);
        auto rDen = solvePart(dens);
        if (!rNum.empty() && !rDen.empty()) {
            BigRat nextVal = rNum[0].nextValue / rDen[0].nextValue;
            double conf = (rNum[0].confidence + rDen[0].confidence) / 2.0;
            
            // Check for overfitting (if sub-patterns are complex)
            bool simple = true;
            auto check = [](const std::string& n) {
                return n.find("Linear")!=std::string::npos || n.find("Constant")!=std::string::npos || n.find("Geometric")!=std::string::npos;
            };
            if (!check(rNum[0].formula) || !check(rDen[0].formula)) simple = false;

            if (simple) conf = 0.99;
            else conf = 0.90; // Strictly capped confidence for complex fractional fits

            SolverResult finalRes = {nextVal, "Fraction Pattern", "", conf, "", "", "Numerator and Denominator solved separately."};
            finalRes.unsimplifiedNum = rNum[0].nextValue.get_str();
            finalRes.unsimplifiedDen = rDen[0].nextValue.get_str();
            results.push_back(finalRes);
        }
        auto wRes = weight->solve(series);
        if (!wRes.empty()) {
            if (wRes[0].confidence > 0.9) {
                results.insert(results.begin(), wRes[0]); 
            } else if (results.empty() || wRes[0].confidence > results[0].confidence + 0.1) {
                results.push_back(wRes[0]);
            } else {
                results.push_back(wRes[0]); 
            }
        }
        std::sort(results.begin(), results.end(), compareResults);
        return results;
    }
};

class GrandUnifiedSolver : public SequenceSolver {
    std::vector<std::unique_ptr<SequenceSolver>> solvers;
public:
    GrandUnifiedSolver() {
        solvers.push_back(std::make_unique<PolynomialSolver>());
        solvers.push_back(std::make_unique<PerfectPowerSolver>());
        solvers.push_back(std::make_unique<VisualSolver>());
        solvers.push_back(std::make_unique<PowerSolver>());
        solvers.push_back(std::make_unique<StandardMathSolver>());
        solvers.push_back(std::make_unique<DigitSolver>());
        solvers.push_back(std::make_unique<BerlekampMasseySolver>());
        solvers.push_back(std::make_unique<RecursiveSolver>());
        solvers.push_back(std::make_unique<HolonomicSolver>());
        solvers.push_back(std::make_unique<PrimeSolver>());
        solvers.push_back(std::make_unique<CyclicOpSolver>());
        solvers.push_back(std::make_unique<InterleavedSolver>());
        solvers.push_back(std::make_unique<DifferenceSolver>());
        solvers.push_back(std::make_unique<FractionSolver>());
        solvers.push_back(std::make_unique<WeightedSolver>());
    }
    std::string getName() const override { return "Grand Unified Solver"; }
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        std::vector<SolverResult> allResults;
        for(auto& s : solvers) {
            auto res = s->solve(series);
            allResults.insert(allResults.end(), res.begin(), res.end());
        }
        std::sort(allResults.begin(), allResults.end(), compareResults);
        return allResults;
    }
};