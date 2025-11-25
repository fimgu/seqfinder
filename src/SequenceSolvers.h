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

// Helper to get absolute value of BigRat
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

inline bool compareResults(const SolverResult& a, const SolverResult& b) {
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
                if (depth > 1 && (size_t)depth >= series.size() - 1) conf = 0.5;
                return {{nextVal, "Polynomial (Degree " + std::to_string(depth) + ")", "", conf, "", "", "Polynomial sequence."}};
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

        // 1. Step Repeat (1, 3, 3, 5, 5...)
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
                    return {{series.back() + nextDiff, "Visual: Step-Repeat", "", 0.99, "", "", "Values repeat once then step."}};
                }
            }
        }

        // 2. Simple Linear
        BigRat diff = series[1] - series[0];
        bool isLin = true;
        for(size_t i=1; i<N; ++i) if (series[i] - series[i-1] != diff) isLin = false;
        if (isLin) return {{series.back() + diff, "Linear", "", 0.95, "", "", "Simple linear progression."}};
        
        // 3. Alternating (A, B, A...)
        if (N >= 3) {
            bool alt = true;
            for(size_t i=2; i<N; ++i) if (series[i] != series[i-2]) alt = false;
            if (alt) {
                return {{series[N-2], "Alternating (Period 2)", "", 0.95, "", "", "Values alternate A, B, A..."}};
            }
        }
        return {};
    }
};

class PowerSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Geometric"; }
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        if (series.size() < 3 || series[0] == 0) return {};
        BigRat ratio = series[1] / series[0];
        if (ratio == 1) return {}; 
        for(size_t i=1; i<series.size()-1; ++i) {
            if (series[i] == 0 || series[i+1] / series[i] != ratio) return {};
        }
        return {{series.back() * ratio, "Geometric", "", 1.0, "", "", "Geometric progression."}};
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
        for (size_t i = 1; i < C.size(); ++i) {
            if (series.size() >= i) nextVal += ((C[i] * -1) / c0) * series[series.size() - i];
        }
        return {{nextVal, "Linear Recurrence", "", 0.85, "", "", "Linear recurrence found via Berlekamp-Massey."}}; 
    }
};

class RecursiveSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Recursive Non-Linear"; }
    
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        std::vector<SolverResult> results;
        size_t N = series.size();
        if (N < 3) return {};
        
        // 1. Basic Relations
        bool sq = true, sylv = true;
        for(size_t i=1; i<N; ++i) {
             if(series[i] != series[i-1]*series[i-1]) sq=false;
             if(series[i] != (series[i-1]*series[i-1]) + series[i-1]) sylv=false;
        }
        if (sq) results.push_back({series.back()*series.back(), "x_n = x_{n-1}^2", "", 1.0, "", "", ""});
        if (sylv) results.push_back({(series.back()*series.back())+series.back(), "x_n = x_{n-1}^2 + x_{n-1}", "", 1.0, "", "", ""});
        
        // 2. Lag 2
        if (N >= 3) {
            bool lag1 = true, lag2 = true;
            for(size_t i=2; i<N; ++i) {
                if (series[i] != (series[i-1] * series[i-2]) - series[i-2]) lag1 = false;
                if (series[i] != (series[i-1] * series[i-2]) - series[i-1]) lag2 = false;
            }
            if (lag1) results.push_back({(series.back() * series[N-2]) - series[N-2], "x_n = x_{n-1}x_{n-2} - x_{n-2}", "", 1.0, "", "", ""});
            if (lag2) results.push_back({(series.back() * series[N-2]) - series.back(), "x_n = x_{n-1}x_{n-2} - x_{n-1}", "", 1.0, "", "", ""});
        }

        // 3. Affine Recurrence: x_n = A * x_{n-1} + B
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
                         // Penalize short sequences unless coefficients are small integers
                         if (A.get_den() == 1 && B.get_den() == 1 && 
                             absRat(A) < 20 && absRat(B) < 20) conf = 0.90;
                         else conf = 0.5;
                     }
                     if (conf > 0.6)
                        results.push_back({A * series.back() + B, "Affine Recurrence", "", conf, "", "", "x_n = A*x_{n-1} + B"});
                 }
             }
        }

        // 4. Variable Linear Recurrence: x_n = x_{n-1}(Ai + B) + (Ci + D)
        if (N >= 5) { // Need at least 4 equations
             std::vector<std::vector<BigRat>> mat;
             std::vector<BigRat> vec;
             for(size_t i=1; i<=4; ++i) { 
                 if (i >= N) break;
                 std::vector<BigRat> row;
                 // FIX: Explicit cast to unsigned long for BigRat conversion
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
                         // FIX: Explicit cast
                         BigRat idx = static_cast<unsigned long>(i);
                         BigRat pred = series[i-1] * (A*idx + B) + (C*idx + D);
                         if (pred != series[i]) { match = false; break; }
                     }
                     if (match) {
                         // FIX: Explicit cast
                         BigRat nextIdx = static_cast<unsigned long>(N);
                         BigRat nextVal = series.back() * (A*nextIdx + B) + (C*nextIdx + D);
                         results.push_back({nextVal, "Variable Linear Recurrence", "", 0.99, "", "", "x_n = x_{n-1}(An+B) + (Cn+D)"});
                     }
                 }
             }
        }

        // 5. Fibonacci/Subtraction
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

        // 6. Variable Increment
        if (N >= 3) {
            bool varInc = true;
            BigRat k = (series[1] - series[0]) - 1;
            for(size_t i=1; i<N; ++i) {
                if (series[i] != series[i-1] + BigRat(static_cast<unsigned long>(i)) + k) { varInc = false; break; }
            }
            if (varInc) results.push_back({series.back() + BigRat(static_cast<unsigned long>(N)) + k, "x_n = x_{n-1} + n + K", "", 1.0, "", "", ""});
        }

        // 7. Coupled Step
        if (N >= 4) {
            bool coupled = true;
            for(size_t i=2; i<N; ++i) {
                size_t offset = 2 + (i % 2); 
                if (i < offset) continue; 
                if (series[i] != series[i-1] + series[i-offset]) { coupled = false; break; }
            }
            if (coupled) {
                size_t offset = 2 + (N % 2);
                results.push_back({series.back() + series[N-offset], "Coupled Step Recurrence", "", 1.0, "", "", ""});
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
                if (nextIdx > 0) return {{BigRat(getNthPrime(nextIdx.get_si())), "Indexed Primes", "", 0.99, "", "", "Sequence of indexed primes."}};
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
                return {{series.back() + nextDiff, "Prime Differences", "", 0.94, "", "", "Differences are primes."}}; 
            }
        }
        return {};
    }
};

class CyclicOpSolver : public SequenceSolver {
    struct Op { char type; BigRat val; }; 
    void generateCombinations(size_t L, size_t currentStep, std::vector<Op>& currentOps, 
        const std::vector<std::vector<Op>>& possibleOpsPerStep, std::vector<std::vector<Op>>& allCombinations) {
        if (currentStep == L) {
            allCombinations.push_back(currentOps);
            return;
        }
        for (const auto& op : possibleOpsPerStep[currentStep]) {
            currentOps[currentStep] = op;
            generateCombinations(L, currentStep + 1, currentOps, possibleOpsPerStep, allCombinations);
        }
    }
public:
    std::string getName() const override { return "Cyclic Operations"; }
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        size_t N = series.size();
        if (N < 4) return {}; 

        for (size_t L = 1; L <= N/2; ++L) {
            std::vector<std::vector<Op>> possibleOpsPerStep(L);
            bool possible = true;
            for (size_t i = 0; i < L; ++i) {
                std::vector<Op> ops;
                // Check +
                BigRat diff = series[i+1] - series[i];
                bool constDiff = true;
                for (size_t k = i; k < N-1; k += L) {
                    if (series[k+1] - series[k] != diff) constDiff = false;
                }
                if (constDiff) ops.push_back({'+', diff});

                // Check *
                if (series[i] != 0) {
                    BigRat ratio = series[i+1] / series[i];
                    bool constRatio = true;
                    for (size_t k = i; k < N-1; k += L) {
                        if (series[k] == 0 || series[k+1] / series[k] != ratio) constRatio = false;
                    }
                    if (constRatio) ops.push_back({'*', ratio});
                }
                
                // Check ^ (Power)
                if (series[i] > 0 && series[i+1] > 0) {
                     // Simple check for square/cube
                     if (series[i+1] == series[i]*series[i]) {
                         bool constPow = true;
                         for(size_t k=i; k<N-1; k+=L) if(series[k+1] != series[k]*series[k]) constPow = false;
                         if(constPow) ops.push_back({'^', 2});
                     }
                }

                if (ops.empty()) { possible = false; break; }
                possibleOpsPerStep[i] = ops;
            }

            if (possible) {
                std::vector<Op> currentOps(L);
                std::vector<std::vector<Op>> allCombinations;
                generateCombinations(L, 0, currentOps, possibleOpsPerStep, allCombinations);
                
                if (!allCombinations.empty()) {
                    auto& ops = allCombinations[0];
                    BigRat nextVal = series.back();
                    size_t stepIdx = (N - 1) % L; 
                    Op op = ops[stepIdx];
                    if (op.type == '+') nextVal = nextVal + op.val;
                    else if (op.type == '*') nextVal = nextVal * op.val;
                    else if (op.type == '^') nextVal = nextVal * nextVal; // Assuming ^2
                    
                    return {{nextVal, "Cyclic Operations (Length " + std::to_string(L) + ")", "", 0.96, "", "", "Operations repeat every L steps."}};
                }
            }
        }
        return {};
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
            auto r4 = recursive->solve(s); // Check recursive
            if (!r4.empty() && r4[0].confidence > 0.5) return r4[0];
            return {0, "", "", 0.0};
        };

        auto bestEven = solveSub(even);
        auto bestOdd = solveSub(odd);
        if (bestEven.confidence >= 0.5 && bestOdd.confidence >= 0.5) {
            double avgConf = (bestEven.confidence + bestOdd.confidence) / 2.0;
            if (bestEven.confidence > 0.85 && bestOdd.confidence > 0.85) avgConf = 0.999; 
            else if (avgConf > 0.9) avgConf = 0.98; 
            
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
             
             auto res = sub->solve(wSeries);
             if (!res.empty() && res[0].confidence > 0.8) {
                  BigInt nextN = static_cast<unsigned long>(series.size() + 1);
                  BigInt w; mpz_ui_pow_ui(w.get_mpz_t(), nextN.get_ui(), p);
                  double finalConf = res[0].confidence;
                  if (p == 1) finalConf += 0.05; // Boost simple n*x_n
                  allResults.push_back({res[0].nextValue / BigRat(w), "Weighted Pattern (n^" + std::to_string(p) + ")", "", finalConf, "", "", ""});
             }
        }

        // 2. Fibonacci Weights (Cube/Fib)
        for(int offset = -2; offset <= 2; ++offset) {
             std::vector<BigRat> wSeries;
             bool possible = true;
             BigInt f1 = 1, f2 = 1; 
             // Generate fib sequence
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
                  allResults.push_back({res[0].nextValue / BigRat(nextFib), "Weighted Pattern (Fibonacci)", "", res[0].confidence, "", "", ""});
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
                    allResults.push_back({res[0].nextValue / BigRat(facts[series.size()+1]), "Weighted Pattern (Factorial)", "", res[0].confidence, "", "", ""});
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

        // 1. Visual/String check on nums/dens
        auto vNum = visual->solve(nums);
        auto vDen = visual->solve(dens);
        if (!vNum.empty() && !vDen.empty()) {
            BigRat nextVal = vNum[0].nextValue / vDen[0].nextValue;
            double conf = (vNum[0].confidence + vDen[0].confidence) / 2.0;
            SolverResult res = {nextVal, "Fraction Pattern", "", conf, "", "", "Numerator and Denominator follow visual patterns."};
            res.unsimplifiedNum = vNum[0].nextValue.get_str();
            res.unsimplifiedDen = vDen[0].nextValue.get_str();
            results.push_back(res);
        }

        // 2. Robust check on nums/dens
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
            
            if (conf > 0.9) conf = 0.99; 
            else conf = std::min(0.99, conf + 0.05);

            SolverResult finalRes = {nextVal, "Fraction Pattern", "", conf, "", "", "Numerator and Denominator solved separately."};
            finalRes.unsimplifiedNum = rNum[0].nextValue.get_str();
            finalRes.unsimplifiedDen = rDen[0].nextValue.get_str();
            results.push_back(finalRes);
        }
        
        // 3. Weighted Solver check
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