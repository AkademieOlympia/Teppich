#!/usr/bin/env python3
"""
Test-Skript zum Prüfen, ob alle Module korrekt importiert werden können.
"""

print("Teste SageMath-Import...")
try:
    from sage.all import *
    print("✓ SageMath erfolgreich importiert")
except ImportError as e:
    print(f"✗ SageMath-Import fehlgeschlagen: {e}")
    exit(1)

print("\nTeste numpy...")
try:
    import numpy as np
    print("✓ numpy erfolgreich importiert")
except ImportError as e:
    print(f"✗ numpy-Import fehlgeschlagen: {e}")
    exit(1)

print("\nTeste pandas...")
try:
    import pandas as pd
    print("✓ pandas erfolgreich importiert")
except ImportError as e:
    print(f"✗ pandas-Import fehlgeschlagen: {e}")
    exit(1)

print("\nTeste matplotlib...")
try:
    import matplotlib.pyplot as plt
    print("✓ matplotlib erfolgreich importiert")
except ImportError as e:
    print(f"✗ matplotlib-Import fehlgeschlagen: {e}")
    exit(1)

print("\nTeste SageMath-Funktionen...")
try:
    # Test einiger SageMath-Funktionen, die im Skript verwendet werden
    CC = ComplexField(100)
    test_matrix = zero_matrix(CC, 2, 2)
    print("✓ SageMath-Funktionen (CC, zero_matrix) funktionieren")
except Exception as e:
    print(f"✗ SageMath-Funktionen fehlgeschlagen: {e}")
    exit(1)

print("\nTeste nth_prime...")
try:
    p1 = nth_prime(1)
    p2 = nth_prime(2)
    print(f"✓ nth_prime funktioniert (p1={p1}, p2={p2})")
except Exception as e:
    print(f"✗ nth_prime fehlgeschlagen: {e}")
    exit(1)

print("\n" + "="*50)
print("✓ ALLE TESTS ERFOLGREICH!")
print("="*50)
print("\nDu kannst jetzt das Hauptskript ausführen mit:")
print("  sage -python 'Hoffbauer Antiprimes.py'")
print("oder")
print("  ./run_with_sage.sh")
