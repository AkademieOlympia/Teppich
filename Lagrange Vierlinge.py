import csv
from sage.all import is_prime, is_squarefree

def generate_energiedoku_csv(limit, filename="energiedoku_viertupel.csv"):
    """
    Ermittelt Zahlen bis 'limit', die als Norm von quaternionischen 
    Primfaktoren (Viertupel e, a, b, c) fungieren.
    """
    results = []
    
    # Wir iterieren durch die natürlichen Zahlen und prüfen auf Quadratfreiheit
    # Im Kontext der #Energiedoku suchen wir oft nach Primzahlen, 
    # die sich als e^2 + a^2 + b^2 + c^2 darstellen lassen.
    
    for n in range(1, limit + 1):
        if is_squarefree(n):
            # In der quaternionischen Zahlentheorie (Lagrange) lässt sich jede 
            # natürliche Zahl als Summe von 4 Quadraten darstellen.
            # Wir extrahieren hier beispielhaft eine Zerlegung für das Viertupel.
            # Für eine vollständige Liste aller Zerlegungen wäre der Rechenaufwand hoch.
            
            # Beispielhafte Zerlegung (e, a, b, c) finden
            dict_decomp = four_number_sum(n)
            if dict_decomp:
                e, a, b, c = dict_decomp
                results.append([n, e, a, b, c])
        
        # Zwischenspeichern alle 100.000 Schritte für Effizienz
        if n % 100000 == 0:
            print(f"Fortschritt: {n} Zahlen geprüft...")

    # Speichern in CSV
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['n', 'e_family', 'a_family', 'b_family', 'c_family'])
        writer.writerows(results)
    
    print(f"Datei {filename} erfolgreich erstellt.")

def four_number_sum(n):
    """Hilfsfunktion zur Zerlegung in 4 Quadrate (e^2 + a^2 + b^2 + c^2 = n)"""
    try:
        from sage.all import Integer
        
        # Konvertiere zu Sage Integer falls nötig
        if not isinstance(n, Integer):
            n = Integer(n)
        
        # Versuche verschiedene Sage-Methoden
        # Methode 1: Direkte Methode (falls verfügbar)
        if hasattr(n, 'is_sum_of_four_squares'):
            res = n.is_sum_of_four_squares(ext=True)
            return list(res)
        
        # Methode 2: Brute-Force für kleine Zahlen
        n_val = int(n)
        if n_val < 10000:
            sqrt_n = int(n_val ** 0.5)
            for e in range(sqrt_n + 1):
                for a in range(sqrt_n + 1):
                    for b in range(sqrt_n + 1):
                        c_sq = n_val - e**2 - a**2 - b**2
                        if c_sq >= 0:
                            c = int(c_sq ** 0.5)
                            if c**2 == c_sq:
                                return [e, a, b, c]
        
        # Methode 3: Klassischer Algorithmus
        if n_val == 0:
            return [0, 0, 0, 0]
        
        sqrt_n = int(n_val ** 0.5)
        for e in range(sqrt_n, -1, -1):
            remainder = n_val - e**2
            if remainder == 0:
                return [e, 0, 0, 0]
            
            sqrt_r = int(remainder ** 0.5)
            for a in range(sqrt_r, -1, -1):
                remainder2 = remainder - a**2
                if remainder2 == 0:
                    return [e, a, 0, 0]
                
                sqrt_r2 = int(remainder2 ** 0.5)
                for b in range(sqrt_r2, -1, -1):
                    remainder3 = remainder2 - b**2
                    if remainder3 >= 0:
                        c = int(remainder3 ** 0.5)
                        if c**2 == remainder3:
                            return [e, a, b, c]
        
        return None
        
    except Exception as e:
        print(f"Fehler bei Zerlegung von {n}: {e}")
        return None

# Interaktive Ausführung mit Dialogeingabe
if __name__ == "__main__":
    print("=== Energiedoku Viertupel Generator ===")
    
    # Eingabe des Limits
    while True:
        try:
            limit_input = input("Bitte geben Sie das Limit ein (z.B. 1000000): ")
            limit = int(limit_input)
            if limit > 0:
                break
            else:
                print("Bitte geben Sie eine positive Zahl ein.")
        except ValueError:
            print("Ungültige Eingabe. Bitte geben Sie eine ganze Zahl ein.")
    
    # Optionale Eingabe des Dateinamens
    filename_input = input("Dateiname (Enter für Standard 'energiedoku_viertupel.csv'): ").strip()
    filename = filename_input if filename_input else "energiedoku_viertupel.csv"
    
    print(f"\nStarte Berechnung bis {limit}...")
    print("(Hinweis: Dies kann einige Zeit in Anspruch nehmen)\n")
    
    generate_energiedoku_csv(limit, filename)