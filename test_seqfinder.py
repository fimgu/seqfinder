import subprocess
import json
import math
import sys
import os

# ============================================================================
# CONFIGURATION
# ============================================================================

EXECUTABLE_PATH = os.path.join(os.getcwd(), "seqfinder_cli.exe")
TOLERANCE = 0.01

# ============================================================================
# ANSI COLORS
# ============================================================================

class Colors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

# ============================================================================
# TEST CASES
# ============================================================================

TEST_CASES = [
    # Prime Solver Tests
    ([2, 3, 5, 7, 11], 13, "Primes: Simple"),
    ([3, 7, 13, 19], 29, "Primes: With Gaps"),
    ([13, 15, 18, 23, 30, 41], 54, "Primes: Differences are primes"),

    # Polynomial Solver Tests
    ([2, 4, 6, 8], 10, "Polynomial: Linear"),
    ([1, 4, 9, 16], 25, "Polynomial: Quadratic (Squares)"),
    ([1, 8, 27, 64, 125], 216, "Polynomial: Cubic (Cubes)"),
    ([0, 1, 3, 6, 10], 15, "Polynomial: Triangular Numbers"),
    ([1, 5, 14, 30, 55], 91, "Polynomial: Square Pyramidal Numbers"),
    ([1, 2, 5, 10, 17], 26, "Polynomial: n^2 + 1"),
    ([0, 6, 24, 60, 120], 210, "Polynomial: n^3 - n"),

    # Fraction Solver Tests
    (["1/2", "1/3", "1/4"], 0.2, "Fractions: Simple (1/n)"),
    (["5/7", "7/7", "7/9", "9/9", "9/11"], 1.0, "Fractions: Complex Unsimplified"), # this one
    ([0.1, 0.2, 0.3, 0.4], 0.5, "Fractions: Numerator Pattern (Decimals)"),
    (["1/1", "1/2", "1/4", "1/8"], 0.0625, "Fractions: Geometric (1/2 ratio)"),
    (["1/2", "2/3", "3/4", "4/5"], 0.8333333333333334, "Fractions: n/(n+1)"),

    # Power Solver Tests
    ([2, 4, 8, 16], 32, "Power: Powers of 2"),
    ([3, 9, 27, 81], 243, "Power: Powers of 3"),
    ([1, 1, 2, 4, 8, 16], 32, "Power: Powers of 2 starting with 1"),
    ([1, 4, 27, 256], 3125, "Power: n^n"),

    # Multiplier Solver Tests
    ([1, 3, 9, 27], 81, "Multiplier: Constant (*3)"),
    ([5, 10, 20, 40], 80, "Multiplier: Double (*2)"),
    ([100, 50, 25], 12.5, "Multiplier: Fractional (*0.5)"),
    ([-2, 4, -8, 16], -32, "Multiplier: Negative common ratio"),

    # Interleaved Solver Tests
    ([1, 10, 2, 20, 3, 30], 4, "Interleaved: Two Sequences"),
    ([2, 5, 4, 10, 6, 15], 8, "Interleaved: Alternating Pattern"),
    ([1, 1, 2, 1, 3, 1, 4], 1, "Interleaved: Constant and increasing"),
    ([1, 2, 1, 4, 1, 8], 1, "Interleaved: Constant and geometric"),

    # Operation Cycle Solver Tests
    ([10, 15, 13, 18, 16], 21, "Cycle: Add/Subtract (+5, -2)"),
    ([2, 4, 2, 4, 2], 4, "Cycle: Multiply/Divide (*2, /2)"),
    ([1, 2, 3, 6, 7, 14], 15, "Cycle: Add 1, Multiply 2"),
    ([10, 5, 15, 7.5, 17.5], 8.75, "Cycle: Divide 2, Add 10"),

    # Linear Recurrence (Fibonacci) Tests
    ([1, 1, 2, 3, 5, 8], 13, "Recurrence: Fibonacci"),
    ([1, 1, 2, 4, 7, 13], 24, "Recurrence: Tribonacci"),
    ([14, 5, -9, -14, -5, 9], 14, "Recurrence: Recursive Subtraction"),
    ([2, 1, 3, 4, 7, 11], 18, "Recurrence: Lucas Numbers"),
    ([1, 2, 3, 5, 8, 13], 21, "Recurrence: Fibonacci starting 1,2"),

    # Difference Pattern Tests
    ([5, 8, 11, 14], 17, "Difference: Constant (3)"),
    ([1, 2, 4, 7, 11], 16, "Difference: Increasing (+1, +2, +3...)"),
    ([1, 3, 6, 10, 15], 21, "Difference: Arithmetic Progression (Triangular)"),
    ([1, 2, 6, 24, 120], 720, "Difference: Factorials (0!, 1!, 2!, 3!, 4!, 5!)"),
    ([1, 2, 4, 8, 15, 26], 42, "Difference: Quadratic differences"),

    # Digit Solver Tests
    ([11, 22, 33], 44, "Digits: Palindromes"),
    ([111, 222, 333], 444, "Digits: Repeating"),
    ([1, 10, 100, 1000], 10000, "Digits: Powers of 10"),

    # Edge Cases
    ([-5, -3, -1, 1], 3, "Edge: Negative Numbers"),
    ([1.5, 3.0, 4.5, 6.0], 7.5, "Edge: Decimal Numbers"),
    ([1000000, 2000000, 3000000], 4000000, "Edge: Large Numbers"),
    ([0, 0, 0, 0], 0, "Edge: All zeros"),
    ([1, 1, 1, 1], 1, "Edge: Constant sequence"),
    ([1, 0, 1, 0], 1, "Edge: Alternating zeros and ones"),

    # Complex Patterns
    ([1, 5, 3, 7, 5, 9, 7], 11, "Complex: Alternating Arithmetic"),
    ([1, 2, 4, 7, 11, 16], 22, "Complex: Cumulative Difference"),
    ([1, 4, 9, 16, 25], 36, "Complex: Position Dependent"),
    ([3, 3, 9, 18, 81, 108], 6561, "Complex: Odd Squared, Even Multiplied"),
    (["3/8", "5/8", "5/10", "7/10", "7/12"], 0.75, "Complex: Alternating Fractions"),
    ([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], 21, "Complex: Long sequence"),
    ([1, 2, 6, 42, 1806], 3263442, "Complex: x_n = x_{n-1} * (x_{n-1} + 1)"),
    ([1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181, 6765, 10946], 17711, "Complex: Very long Fibonacci"),
    ([1, 2, 4, 5, 7, 8, 10, 11], 13, "Complex: Skip one number"),
    ([1, 4, 27, 256, 3125], 46656, "Complex: n^n"),
    ([1, 11, 21, 1211, 111221], 312211, "Complex: Look-and-say sequence"),

    # User Requested Test Cases (New)
    ([4, 3, 7, 9, 10, 27], 13, "Cycle: Alternating +3, *3"),
    ([7, 6, 8, 24, 6, 5, 7], 21, "Cycle: -1, +2, *3, /4"),
    ([11, 9, 6, 1, -6, -17], -30, "Difference: Subtracting Primes"),
    ([102, 34, 17, 18, 6, 3], 4, "Cycle: /3, /2, +1"),
    (["1/4", "3/4", "3/6", "5/6", "5/8"], 7/8, "Fractions: Complex Pattern"),
    ([2, 1, 2, 5, 13, "73/2"], 114.5, "Complex: Variable Multiplier (*0.5, *1+1...)"),
    ([1, 2, 5, 11, 26, 59], 137, "Complex: Derived Sequence"),
    ([1, "1/4", "1/9", "5/64", "8/125", "13/216"], 21/343, "Fractions: Fib / Cube"),
    ([1, 4, 9, "64/5", "125/8", "216/13"], 343/21, "Fractions: Cube / Fib"),
    ([3, 2, 4, 4, 12, 36], 396, "Complex: (Prev * Prev2) - Prev2"),
    ([2, "3/2", 5, "23/2", 28, "135/2"], 163, "Complex: Mixed Fractions Sequence"),
    (["2/3", "1/6", "5/6", "1/12", "53/60", "1/20"], 127/140, "Interleaved: Harmonic Logic"),
    ([9, 10, 6, 15, 4, "45/2"], 8/3, "Interleaved: Mixed Ops (8/3)"),
    ([3, 2, 8, 0.5, 23, 2], 68, "Interleaved: +5,+15... & +/-1.5"),
    ([18, 6, 36, 24, 8, 64], 42 + 2/3, "Cycle: /3, ^2, *2/3"),
    ([3, 7, 13, 19, 29, 37], 43, "Primes: Gap Pattern"),
    ([4, 5, 9, 13, 22, 31], 53, "Complex: Variable Increment Sum"),
    ([1, "9/10", "5/6", "11/14", "3/4", "13/18"], 7/10, "Fractions: Complex Descending"),
    ([17, 23, 22, 16, 18, 24, 21], 15, "Difference: Alt +6 / Count"),
]

# ============================================================================
# TEST RUNNER
# ============================================================================

def run_test(sequence, expected, name):
    """Run a single test case"""
    print(f"Running test: {name}...", flush=True)
    # Convert inputs to strings for the CLI
    args = [str(x) for x in sequence]
    
    try:
        result = subprocess.run(
            [EXECUTABLE_PATH] + args, 
            capture_output=True, 
            text=True, 
            check=False
        )
        
        if result.returncode != 0:
            print(f"{Colors.FAIL}[FAIL] {name}{Colors.ENDC}")
            print(f"  Error: Executable returned code {result.returncode}")
            return False

        try:
            output_json = json.loads(result.stdout)
        except json.JSONDecodeError:
            print(f"{Colors.FAIL}[FAIL] {name}{Colors.ENDC}")
            print(f"  Error: Could not parse JSON output.")
            print(f"  Raw output: {result.stdout.strip()}")
            return False

        if "error" in output_json:
            print(f"{Colors.FAIL}[FAIL] {name}{Colors.ENDC}")
            print(f"  Solver returned error: {output_json['error']}")
            return False

        actual = output_json.get("next_number")
        pattern_found = output_json.get("pattern", "Unknown")
        formula = output_json.get("formula", "")
        
        if actual is None:
            print(f"{Colors.FAIL}[FAIL] {name}{Colors.ENDC}")
            print(f"  No next_number provided in output.")
            return False

        # Helper to parse potential fraction string to float
        def parse_val(val):
            if isinstance(val, str) and '/' in val:
                try:
                    n, d = val.split('/')
                    return float(n) / float(d)
                except ValueError:
                    return float(val) # Fallback
            return float(val)

        actual_float = parse_val(actual)

        if math.isclose(actual_float, expected, abs_tol=TOLERANCE):
            print(f"{Colors.OKGREEN}[PASS] {name}{Colors.ENDC} (Pattern: {pattern_found})")
            return True
        else:
            # Check candidates
            candidates = output_json.get("candidates", [])
            for cand in candidates:
                c_next = cand.get("next")
                if c_next is not None:
                    c_next_float = parse_val(c_next)
                    if math.isclose(c_next_float, expected, abs_tol=TOLERANCE):
                        print(f"{Colors.OKGREEN}[PASS] {name}{Colors.ENDC} (Found in candidates: {cand.get('pattern')})")
                        return True

            print(f"{Colors.FAIL}[FAIL] {name}{Colors.ENDC}")
            print(f"  Expected: {expected}")
            print(f"  Got:      {actual} (Parsed: {actual_float}) (Pattern: {pattern_found})")
            print(f"  Formula:  {formula}")
            print(f"  Candidates checked: {len(candidates)}")
            return False

    except Exception as e:
        print(f"{Colors.FAIL}[FAIL] {name} - Exception: {e}{Colors.ENDC}")
        return False

# ============================================================================
# MAIN
# ============================================================================

def main():
    print(f"{Colors.HEADER}Running Sequence Solver Tests...{Colors.ENDC}")
    print(f"Executable: {EXECUTABLE_PATH}\n")
    
    passed = 0
    total = len(TEST_CASES)
    
    for seq, expected, name in TEST_CASES:
        if run_test(seq, expected, name):
            passed += 1
            
    print(f"\n{Colors.HEADER}Summary:{Colors.ENDC}")
    
    if passed == total:
        print(f"{Colors.OKGREEN}All {total} tests passed!{Colors.ENDC}")
        sys.exit(0)
    else:
        print(f"{Colors.FAIL}{total - passed} tests failed.{Colors.ENDC}")
        print(f"{Colors.OKGREEN}{passed} tests passed.{Colors.ENDC}")
        sys.exit(1)

if __name__ == "__main__":
    main()