# QZG Test-Prozedur: Dokumentation

## Übersicht

Dieses Skript implementiert die formale Testprozedur zur quantitativen Prüfung der QZG-Behauptung (Quaternionische Zustandsgleichung des Primzahl-Gases). Es basiert auf:

1. **Lemma**: Exakte Fixierung eines Eigenwerts durch Nullraum-Bedingung
2. **Definition**: Quaternion-Formfaktor Q(k) als prüfbare, endliche Größe
3. **Minimale Testprozedur**: Falsifizierbar und robust

## Installation

Das Skript benötigt:
- SageMath (für numerische Berechnungen)
- NumPy
- Pandas
- SciPy

Installation:
```bash
# SageMath installieren (falls nicht vorhanden)
# Dann:
sage -pip install numpy pandas scipy
```

## Verwendung

### Basis-Ausführung

```bash
sage -python QZG_Test_Prozedur.py
```

### Parameter anpassen

Bearbeiten Sie die `run_full_test_procedure()` Funktion am Ende der Datei:

```python
results = run_full_test_procedure(
    K=1000,              # Startwert für k (Intervall [K, 2K])
    m=5,                 # Fenstergröße (Fensterradius)
    epsilon=1e-3,        # Kopplungskonstante
    beta=0.01,           # Temperatur-Parameter
    taper_length=2.0,    # Taper-Länge für Gewichte
    M_permutations=1000, # Anzahl Permutationen für Nullmodell
    run_robustness=True  # Robustheitstest durchführen
)
```

## Komponenten

### 1. Lemma-Implementierung

Die Funktion `construct_dirac_block()` konstruiert den Dirac-Operator D mit der Nullraum-Bedingung:
- `R * e_0 = 0` und `R* * e_0 = 0`
- Dies garantiert, dass `t_k` ein exakter Eigenwert bleibt

Die Funktion `verify_lemma()` verifiziert die Aussage des Lemmas numerisch.

### 2. Quaternion-Formfaktor Q(k)

Die Funktion `compute_formfactor_Q()` berechnet:

```
Q(k) = (1/Z_k) * sum_{i≠j in W_k} w_ij * sgn(g(i), g(j); g(k))
```

wobei:
- `g(n)` die Primzahl-Klasse (E, A, B, C) basierend auf `p_n mod 12` ist
- `sgn(X, Y; Z)` die Quaternion-Orientierung kodiert (+1, -1, oder 0)
- `w_ij` Gewichte sind (entweder konstant oder Taper-Gewichte)

### 3. Testprozedur

Die Testprozedur besteht aus 4 Schritten:

1. **Daten erzeugen**: Für k in [K, 2K] werden berechnet:
   - Q(k) und Q_0(k) (Formfaktoren)
   - G(k) (lokaler Gap = p_{k+1} - p_k)
   - S_k(β) (spektrale Entropie, optional)

2. **Nullmodell**: Permutation der Klassenlabels innerhalb jedes Fensters (M Wiederholungen)

3. **Teststatistiken**:
   - Korrelation `corr(Q(k), G(k))`
   - Regression: `G(k) = a + b*Q(k) + c*log(p_k) + Fehler`
   - Vergleich mit Permutationsverteilung → p-Wert

4. **Robustheitstest** (optional): Wiederholung für mehrere Parameter-Kombinationen

## Ausgabe

Das Skript erzeugt folgende Dateien:

- `QZG_test_data.csv`: Hauptdaten (k, p_k, Q_k, Q_0_k, G_k, S_k)
- `QZG_permutation_data.csv`: Permutationsdaten für Nullmodell
- `QZG_robustness_results.csv`: Robustheitstest-Ergebnisse (falls aktiviert)

## Interpretation der Ergebnisse

### Statistische Signifikanz

- **Signifikant**: `|corr(Q(k), G(k))| > 0.1` und `p < 0.05`
  → Tao-kompatibler Claim: "Q(k) zeigt eine statistisch signifikante, robuste Assoziation mit G(k)"

- **Nicht signifikant**: `p >= 0.05` oder `|corr| <= 0.1`
  → Die Metapher sollte entschärft werden

### Robustheit

Ein echter Effekt muss **stabil** bleiben über verschiedene Parameter-Kombinationen:
- Verschiedene Fenstergrößen (m)
- Verschiedene Kopplungskonstanten (ε)
- Verschiedene Temperaturen (β)
- Verschiedene Gewichtsschemata (taper_length)

## Hinweise

1. **Riemann-Nullstellen**: Das Skript verwendet aktuell eine Näherung für t_k. Für echte Tests sollten echte Nullstellen aus `zeros1.gz` oder ähnlichen Quellen geladen werden.

2. **Performance**: Die Berechnung der spektralen Entropie kann für große Fenster langsam sein. Setzen Sie `compute_entropy=False` für schnellere Tests.

3. **Permutationen**: Für robuste p-Werte sollten mindestens M=1000 Permutationen durchgeführt werden.

## Erweiterungen

Mögliche Erweiterungen:
- Laden echter Riemann-Nullstellen aus Dateien
- Parallelisierung der Permutationen
- Visualisierung der Ergebnisse (Plots)
- Integration in das LaTeX-Paper
