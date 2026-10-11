# resortes-mott

Herramienta de cálculo general de resortes helicoidales de compresión
(Python 3, sin dependencias externas) para **analizar** resortes existentes o
**diseñar** nuevos por tanteo (método 1, sección 19-6).

Fórmulas tomadas de MOTT, Robert L. *Diseño de Elementos de Máquinas*, 4ta ed.,
Pearson, 2006 — Capítulo 19 "Resortes" (desde la pág. 735): ecuaciones (19-1) a
(19-8), tablas 19-2 (calibres), 19-3 (materiales), 19-4 (módulos G y E) y
figuras 19-8 a 19-15 (esfuerzos admisibles y pandeo).

## Uso

```bash
python3 resortes_mott.py               # modo interactivo (analizar o diseñar)
python3 resortes_mott.py -h            # referencia de fórmulas y API
```

Desde código:

```python
from resortes_mott import analizar, disenar

analizar(material="A228", Dw=0.055, DE=0.561, N=10, Lf=1.75, Fo=14.0,
         extremos="escuadrados_rectificados")
disenar(material="A227", F1=50, L1=70, F2=120, L2=55, DI_min=18, unidades="mm")
```

## Nota

Los esfuerzos admisibles son aproximaciones lineales de las figuras 19-8..19-13
del libro; para un diseño definitivo verificar las gráficas del texto.
