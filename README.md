Par_Comp
========

Parallel Computing - Matrix multiplication speedup


# Abgabe
bis einschließlich 17. November 2013

# Ziel
der Aufgabe ist es, eine Routine für eine optimierte Matrix-Matrix Multiplikation auf SuperMUC zu erstellen (nur für den Spezialfall von **quadratischen Matrizen**). 

# Ausgangspunkt 
ist die nicht optimierte Version mit 3 geschachtelten Schleifen: mmult.c

Versuchen Sie für N=1000x1000 Matrizen eine möglichst gute Leistung zu erreichen, experimentieren sie mit verschiedenen Optimierungsmöglichkeiten, wie etwa:
- Compiler: Intel compiler _icc_ vs _gcc_
- Compiler flags _(-O2/-O3/...)_ und Annotationen; siehe auch
http://www.lrz.de/services/compute/supermuc/programming/#TOC4
- Einführung von Blocking / Tiling
- Änderung der Anordnung der Matrizen im Speicher
- Loop Unrolling
- SSE/AVX assembler / intrinsics

