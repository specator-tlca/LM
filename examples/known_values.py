"""
Known values of Ramanujan tau function and Hecke eigenvalues
for testing and verification
"""

# Ramanujan tau function τ(p) for first primes
RAMANUJAN_TAU = {
    2: -24,
    3: 252,
    5: -4830,
    7: -16744,
    11: 534612,
    13: -577738,
    17: -6905934,
    19: 10661420,
    23: -94196280,
    29: 915195642,
    31: -1127792134,
    37: -12946668786,
    41: 59717358024,
    43: -58263944292,
    47: 495868873020,
    53: -4107941704602,
    59: 12523346350488,
    61: 18664062532404,
    67: -111853766027598,
    71: 216857579623160,
    73: -430299993094738,
    79: 2268778246060424,
    83: -5268635602263100,
    89: 17290748347478024,
    97: -88706983025139602
}

# Convert to lambda(p) = tau(p) / p^(11/2)
def tau_to_lambda(p, tau):
    from math import pow
    return tau / pow(p, 11.0/2.0)

# First cusp form of level 11, weight 2
# a_p values (not normalized)
LEVEL11_WEIGHT2_AP = {
    2: -2,
    3: -1,
    5: 1,
    7: -2,
    11: 1,  # This is special (bad prime)
    13: 4,
    17: -2,
    19: 0,
    23: -1,
    29: 0,
    31: 7,
    37: 3,
    41: -8,
    43: -6,
    47: 8,
    53: -6,
    59: 5,
    61: 12,
    67: -7,
    71: -3,
    73: 4,
    79: -4,
    83: -1,
    89: 7,
    97: 2
}

# Convert to normalized eigenvalues
def ap_to_lambda(p, ap, weight=2):
    from math import pow
    return ap / pow(p, (weight-1)/2.0)


def create_test_files():
    """Create CSV files with known values"""
    import csv
    
    # Ramanujan lambda
    with open('ramanujan_known.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['p', 'tau', 'lambda', 'type'])
        for p, tau in RAMANUJAN_TAU.items():
            lam = tau_to_lambda(p, tau)
            writer.writerow([p, tau, lam, 'tau'])
    print("Created ramanujan_known.csv")
    
    # Level 11 form
    with open('level11_known.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['p', 'ap', 'lambda', 'type'])
        for p, ap in LEVEL11_WEIGHT2_AP.items():
            lam = ap_to_lambda(p, ap, 2)
            writer.writerow([p, ap, lam, 'lambda'])
    print("Created level11_known.csv")

if __name__ == "__main__":
    create_test_files()
    
    # Show some values
    print("\nRamanujan λ(p) = τ(p)/p^(11/2):")
    for p in [2, 3, 5, 7, 11]:
        tau = RAMANUJAN_TAU[p]
        lam = tau_to_lambda(p, tau)
        print(f"p={p}: τ(p)={tau}, λ(p)={lam:.8f}")
    
    print("\nLevel 11, weight 2 normalized eigenvalues:")
    for p in [2, 3, 5, 7, 11]:
        ap = LEVEL11_WEIGHT2_AP[p]
        lam = ap_to_lambda(p, ap)
        print(f"p={p}: a_p={ap}, λ(p)={lam:.8f}")
