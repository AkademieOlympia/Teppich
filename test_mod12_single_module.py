#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test-Skript für mod12_single_module
Führt den Code aus dem Notebook aus, um Probleme zu identifizieren.
"""

# SageMath Imports - Nutze die schnellsten verfügbaren Bibliotheken
from sage.all import *
import csv
import sys

# SageMath's prime_range() ist deutlich schneller als manueller Sieb
# SageMath's Integer-Typ ist optimiert für große Zahlen
# factor() ist optimiert für große Zahlen
print("✓ SageMath erfolgreich geladen!")
print(f"SageMath Version: {version()}")
print(f"Python Version: {sys.version}")
sys.stdout.flush()  # Stelle sicher, dass die Ausgabe sofort angezeigt wird

# Konfiguration
GROUP_ORDER = [1, 5, 7, 11]

# Beispiel: Obergrenze n
n = 1000
out_prefix = "mod12_single_module"

def factorization_string(factor_dict):
    """
    Erzeugt eine schöne Faktorisierungsdarstellung wie 5^3*17*29^2.
    Nutzt SageMath's factor() Ergebnis.
    """
    parts = []
    # SageMath's factor() gibt ein Factorization-Objekt zurück
    # Wir konvertieren es zu sortierten Paaren (p, e)
    for p, e in sorted(factor_dict.items()):
        if e == 1:
            parts.append(str(p))
        else:
            parts.append(f"{p}^{e}")
    return "*".join(parts)

def factorization_string_from_sage_factor(fac):
    """
    Konvertiert SageMath's Factorization-Objekt zu String.
    """
    if isinstance(fac, Factorization):
        parts = []
        for p, e in sorted(fac):
            if e == 1:
                parts.append(str(p))
            else:
                parts.append(f"{p}^{e}")
        return "*".join(parts)
    else:
        # Fallback für dict
        return factorization_string(fac)

def generate_group_numbers(n, primes, residue):
    """
    Generiere alle Zahlen <= n, deren Primfaktoren nur aus 'primes' stammen (eine Restklasse).
    Vermeidet Duplikate durch Erzwingen nicht-abnehmender Primzahl-Indizes in der Rekursion.
    
    Args:
        n: Obergrenze
        primes: Liste von Primzahlen (eine Restklasse mod 12)
        residue: Restklasse mod 12 (1, 5, 7 oder 11)
    
    Returns:
        Liste von Dictionaries mit Informationen über jede generierte Zahl
    """
    records = []
    
    # Precompute Indizes als 1-basierte Positionen innerhalb dieser Restklasse
    # primes[i] hat Index i+1
    def dfs(start_i, current_val, idx_list, fac_counter):
        # Speichere aktuelles Produkt wenn > 1 (1 ausschließen)
        if current_val > 1:
            # Nutze SageMath's factor() für saubere Faktorisierung
            # Aber wir haben bereits fac_counter, also nutzen wir das
            fac_str = factorization_string(fac_counter)
            
            records.append({
                "group_mod12": residue,
                "value": current_val,
                "indices": tuple(idx_list),  # wiederholte Indizes für Multiplizität
                "min_index": idx_list[0] if idx_list else None,
                "max_index": idx_list[-1] if idx_list else None,
                "len": len(idx_list),
                "factorization": fac_str,
            })
        
        # Versuche Multiplikation mit Primzahlen ab start_i (Wiederholungen erlauben)
        for i in range(start_i, len(primes)):
            p = primes[i]
            # Nutze SageMath Integer für präzise Arithmetik
            if current_val > n // p:
                break
            
            # Multipliziere mit p einmal
            idx_list.append(i + 1)  # 1-basierter Index in dieser Gruppe
            fac_counter[p] = fac_counter.get(p, 0) + 1
            dfs(i, current_val * p, idx_list, fac_counter)  # i (nicht i+1) => Wiederholungen erlauben
            # Backtrack
            fac_counter[p] -= 1
            if fac_counter[p] == 0:
                del fac_counter[p]
            idx_list.pop()
    
    dfs(0, 1, [], {})
    return records

def write_csv(path, rows, fieldnames):
    """Schreibe CSV-Datei mit UTF-8 Kodierung."""
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)

# Hauptausführung
print("\n" + "="*60)
print("HAUPTPROGRAMM STARTET")
print("="*60 + "\n")

# Prüfe Eingabe
if n < 1:
    raise ValueError("n muss >= 1 sein.")

print(f"Starte Berechnung für n = {n}...")
print(f"Nutze SageMath's optimierte Primzahl-Funktionen.\n")
sys.stdout.flush()

# Primzahlen bis n - Nutze SageMath's prime_range() (viel schneller!)
all_primes = prime_range(n + 1)
print(f"Gefunden: {len(all_primes)} Primzahlen <= {n}")
sys.stdout.flush()

# Gruppiere Primzahlen nach Restklasse mod 12 (nur 1, 5, 7, 11)
groups = {r: [] for r in GROUP_ORDER}
for p in all_primes:
    r = p % 12
    if r in groups:
        groups[r].append(p)

print("\nVerteilung nach Restklassen mod 12:")
for r in GROUP_ORDER:
    print(f"  Klasse {r}: {len(groups[r])} Primzahlen")
sys.stdout.flush()

# Generiere Datensätze pro Gruppe
all_records = []
for r in GROUP_ORDER:
    ps = groups[r]
    if not ps:
        continue
    print(f"\nGeneriere Zahlen für Restklasse {r}...")
    sys.stdout.flush()
    records = generate_group_numbers(n, ps, r)
    all_records.extend(records)
    print(f"  {len(records)} Zahlen generiert")
    sys.stdout.flush()

print(f"\nGesamt: {len(all_records)} Zahlen <= {n} erzeugt")
sys.stdout.flush()

# Füge zusätzliche Spalten für CSV-Lesbarkeit hinzu
for rec in all_records:
    # Indizes-Tupel -> String für CSV-Lesbarkeit
    rec["indices_str"] = " ".join(map(str, rec["indices"]))
    del rec["indices"]

# Ausgabe-Schema
fields = [
    "group_mod12",
    "value",
    "factorization",
    "indices_str",
    "len",
    "min_index",
    "max_index",
]

# Sortierung 1: Nach Modul-Gruppenordnung, dann nach Index-Muster (lexikographisch), dann nach Wert
def parse_indices(s):
    """Parse Indizes-String zurück zu Tupel."""
    return tuple(int(x) for x in s.split()) if s else tuple()

group_rank = {r: i for i, r in enumerate(GROUP_ORDER)}

by_index = sorted(
    all_records,
    key=lambda rec: (
        group_rank.get(rec["group_mod12"], 999),
        parse_indices(rec["indices_str"]),
        rec["value"],
    )
)

# Sortierung 2: Nach numerischem Wert, dann Gruppe, dann Faktorisierung für Stabilität
by_value = sorted(
    all_records,
    key=lambda rec: (
        rec["value"],
        group_rank.get(rec["group_mod12"], 999),
        rec["factorization"],
    )
)

print("Sortierungen abgeschlossen.")
sys.stdout.flush()

# Schreibe CSV-Dateien
out1 = f"{out_prefix}_by_index.csv"
out2 = f"{out_prefix}_by_value.csv"

write_csv(out1, by_index, fields)
write_csv(out2, by_value, fields)

print(f"\n✓ Fertig. {len(all_records)} Zahlen <= {n} erzeugt.")
print(f"CSV 1 (Index-Ordnung): {out1}")
print(f"CSV 2 (Wert-Ordnung):  {out2}")
sys.stdout.flush()

# Zeige Beispiel-Ausgabe
print("\n" + "="*60)
print("BEISPIEL-AUSGABE")
print("="*60 + "\n")

try:
    import pandas as pd
    try:
        df = pd.read_csv(out2)
        print("Erste 20 Einträge (sortiert nach Wert):")
        print(df.head(20).to_string(index=False))
    except FileNotFoundError:
        print(f"Warnung: Datei {out2} nicht gefunden. Zeige manuell aus by_value:")
        if len(by_value) > 0:
            for i, rec in enumerate(by_value[:20]):
                print(f"{i+1:3d}. Gruppe {rec['group_mod12']:2d} | Wert: {rec['value']:6d} | {rec['factorization']:20s} | Indizes: {rec['indices_str']}")
        else:
            print("Keine Daten verfügbar.")
except ImportError:
    print("pandas nicht verfügbar. Zeige manuell:")
    if len(by_value) > 0:
        for i, rec in enumerate(by_value[:20]):
            print(f"{i+1:3d}. Gruppe {rec['group_mod12']:2d} | Wert: {rec['value']:6d} | {rec['factorization']:20s} | Indizes: {rec['indices_str']}")
    else:
        print("Keine Daten verfügbar.")

print("\n" + "="*60)
print("PROGRAMM ABGESCHLOSSEN")
print("="*60)
sys.stdout.flush()
