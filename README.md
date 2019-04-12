# Numerische Methoden für technische Physiker

In diesem Repository befinden sich sämtliche Übungsbeispiele die ich im Rahmen meines Physikstudiums (TU Wien) im Kurs 'Numerische Methoden für Technische Physiker 2' erstellt habe.

In jedem Übungsordner befindet sich die zu lösende Aufgabenstellung (als pdf), ein FORTRAN95 Skript das die Aufgabenstellung löst, sowie weitere Hilffskripte zum Plotten mit GnuPlot, notwendige Daten, usw.
Sämtlicher Code der hier zugägnlich ist, wurde von mir selbst geschrieben. (In UE7 befinden sich Daten von NMR Bildern die mit Erlaubnis von Prof. Karsten Held (TU Wien) bereitgestellt wurden, ich möchte ihm dafür nochmals persöhnlich danken.)

Die numerischen Aufgabenstellungen der einzelnen Übungsbeispiele stellen klassische Probleme der Physik dar:
- UE1:
  - Rekursive Berechnung von Von-Neumann Funktionen und Binominalkoeffizienten
  - Lösung der zeitabhängigen Schrödingergleichung (mithilfe des Thohmas Algorithmus für tridiagonale Gleichungssysteme)
- U2: Simulation des Wannier-Stark Effekts, d.h. Lokalisierung von Elektronen im Festkörper bei angelegtem, äußeren elektrischen Feld (Lösung des Eigenwertproblem)
- U3:
  - Lösung der Laplacegleichung für ein elektrisches Potential (mit Gauss Algorithmus)
  - Lösung der zeitabhängigen Schrödingergleichung mit Split-Operator Methode (mit FFT)
- UE4: Lösung des Keplerproblems (Lösung der DGL mit Newton Verfahren und RK4 Verfahren)
- UE5: Berechnung der Radialfunktionen des Heliumatoms (Lösung der SGL mit Numerov Methode)
- UE6: Simulation des 2d Ising Modells und Abschätzung der kritischen Temperatur für den Phasenübergang einfacher Ferromagneten (mit Markov Ketten Monte Carlo Simulation)
- UE7: Segmentierung einer Kernspinresonanzaufnahme (NMR) des Gehirns (graue Masse, weiße Masse, Schädel, etc.) (mittels Simulated Annealing)

Einige meiner scripts verwenden subroutinen aus allgemein zugänglichen Bibliotheken (z.B. LAPACK, CERNLIB), darauf weise ich falls nötig hin.

Viel Spaß beim Basteln mit meinen Programmierbeispielen :-) .
