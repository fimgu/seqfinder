import requests
import json
import math
import sys
import time
from test_seqfinder import TEST_CASES

# Configuration
URL = "http://127.0.0.1:5000/predict"
TOLERANCE = 0.01

# Colors
class Colors:
    HEADER = '\033[95m'
    OKGREEN = '\033[92m'
    FAIL = '\033[91m'
    WARNING = '\033[93m'
    ENDC = '\033[0m'

def parse_val(val):
    if isinstance(val, str) and '/' in val:
        try:
            n, d = val.split('/')
            return float(n) / float(d)
        except ValueError:
            return float(val)
    return float(val)

def verify_all():
    print(f"{Colors.HEADER}Verifying {len(TEST_CASES)} sequences against GUI Endpoint...{Colors.ENDC}")
    print(f"Target: {URL}\n")
    
    failures = []
    passed_count = 0

    # Ensure server is reachable
    try:
        requests.get("http://127.0.0.1:5000/", timeout=2)
    except requests.exceptions.ConnectionError:
        print(f"{Colors.FAIL}Error: Could not connect to Flask server.{Colors.ENDC}")
        print("Please run 'python app.py' in a separate terminal first.")
        sys.exit(1)

    for sequence, expected, name in TEST_CASES:
        seq_str = ", ".join(str(x) for x in sequence)
        
        try:
            response = requests.post(URL, json={'sequence': seq_str}, timeout=5)
            
            if response.status_code != 200:
                print(f"{Colors.FAIL}[FAIL] {name}{Colors.ENDC}: HTTP {response.status_code}")
                failures.append((name, f"HTTP {response.status_code}"))
                continue
                
            data = response.json()
            
            if 'error' in data:
                print(f"{Colors.FAIL}[FAIL] {name}{Colors.ENDC}: API Error - {data['error']}")
                failures.append((name, data['error']))
                continue
                
            actual_raw = data.get('next_number')
            pattern = data.get('pattern', 'Unknown')
            candidates = data.get('candidates', [])
            
            if actual_raw is None:
                print(f"{Colors.FAIL}[FAIL] {name}{Colors.ENDC}: No result returned")
                failures.append((name, "No result"))
                continue
                
            actual_float = parse_val(actual_raw)
            
            # STRICT CHECK: Top result must match
            if math.isclose(actual_float, expected, abs_tol=TOLERANCE):
                print(f"{Colors.OKGREEN}[PASS] {name}{Colors.ENDC} -> {actual_raw} ({pattern})")
                passed_count += 1
            else:
                # It failed strict check. Did we find it in candidates?
                found_in_candidates = False
                cand_pattern = ""
                for cand in candidates:
                    c_val = parse_val(cand.get('next'))
                    if math.isclose(c_val, expected, abs_tol=TOLERANCE):
                        found_in_candidates = True
                        cand_pattern = cand.get('pattern')
                        break

                print(f"{Colors.FAIL}[FAIL] {name}{Colors.ENDC}")
                print(f"  Input:    {seq_str}")
                print(f"  Expected: {expected}")
                print(f"  Got Top:  {actual_raw} ({pattern})")
                
                if found_in_candidates:
                    print(f"  {Colors.WARNING}Note: Correct answer found in candidates ({cand_pattern}), but was ranked lower.{Colors.ENDC}")
                    failures.append((name, f"Ranking Error: Got {actual_raw}, Expected {expected} (Found in candidates)"))
                else:
                    failures.append((name, f"Calculation Error: Got {actual_raw}, Expected {expected}"))
                
        except Exception as e:
            print(f"{Colors.FAIL}[FAIL] {name}{Colors.ENDC}: Exception - {e}")
            failures.append((name, str(e)))

    print("\n" + "="*50)
    if failures:
        print(f"{Colors.FAIL}Summary: {passed_count}/{len(TEST_CASES)} passed.{Colors.ENDC}")
        print(f"{len(failures)} failures detected.")
        sys.exit(1)
    else:
        print(f"{Colors.OKGREEN}Summary: {passed_count}/{len(TEST_CASES)} passed.{Colors.ENDC}")
        print("All strict tests passed! The GUI is showing the correct top results.")
        sys.exit(0)

if __name__ == "__main__":
    verify_all()