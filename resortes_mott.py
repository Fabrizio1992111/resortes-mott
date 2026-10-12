#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=============================================================================
 RESORTES HELICOIDALES DE COMPRESIÓN — Herramienta general de cálculo
 Fórmulas de: MOTT, Robert L. "Diseño de Elementos de Máquinas", 4ta ed.
              Pearson, 2006 — Capítulo 19 "Resortes" (desde la pág. 735)
=============================================================================
 Fórmulas implementadas (ecuaciones del libro):
   (19-1)   k   = ΔF/ΔL;  (19-1a) k = (Fo-Fi)/(Li-Lo);  (19-1b) k = Fo/(Lf-Lo)
          C     = Dm / Dw                          índice del resorte
          Na    = N-2 (escuadrados y rectificados), N (planos),
                  N-1 (planos y rectificados)       espiras activas
          Ls    = Dw · N                           longitud comprimida
   (19-2)   α   = arctan(p / (π · Dm))             ángulo de paso
   (19-3)   DEs = √(Dm² + (p² - Dw²)/4) + Dw       DE en estado comprimido
          cc    = (Lo - Ls) / Na                   holgura de espira (> Dw/10)
   (19-4)   τ   = 8·K·F·Dm / (π·Dw³) = 8·K·F·C / (π·Dw²)    esfuerzo cortante
   (19-5)   K   = (4C - 1)/(4C - 4) + 0.615/C      factor Wahl
   (19-6)   f   = 8·F·Dm³·Na / (G·Dw⁴) = 8·F·C³·Na / (G·Dw)  deflexión
   (19-7)   Dw  = ∛(3.06·Fo·Dm / τd)               diámetro tentativo (diseño)
   (19-8)   Na  = G·Dw / (8·k·C³)                  nº de espiras activas
 Tablas 19-2 (calibres), 19-3 (materiales), 19-4 (módulos G y E).
 Pandeo: Fig. 19-15 (curvas A, B, C según extremos).
-----------------------------------------------------------------------------
 USO:
   python3 resortes_mott.py              # modo interactivo (analizar o diseñar)
   desde código:
      from resortes_mott import analizar, disenar
      analizar(material="A228", Dw=0.055, DE=0.561, N=10, Lf=1.75,
               Fo=14.0, extremos="escuadrados_rectificados")
      disenar(material="A227", F1=50, L1=70, F2=120, L2=55,
              DI_min=18, unidades="mm")
 Nota: los esfuerzos admisibles de las figuras 19-8..19-13 son gráficos del
 libro; aquí se usan aproximaciones lineales (ver TAU_ANCORS). Para diseño
 definitivo, verificar las gráficas del texto.
=============================================================================
"""

import math
import sys

PI = math.pi

# --------------------------------------------------------------------------
# Tabla 19-4: módulos de elasticidad (psi y GPa)
# --------------------------------------------------------------------------
G_PSI = {
    "A227": 11.5e6, "A228": 11.85e6, "A229": 11.2e6,
    "A231": 11.2e6, "A401": 11.2e6, "A313": 11.2e6,
    "B134": 10.0e6, "B159": 5.0e6, "B197": 6.0e6,
}
G_GPA = {
    "A227": 79.3, "A228": 81.7, "A229": 77.2,
    "A231": 77.2, "A401": 77.2, "A313": 77.2,
}

# --------------------------------------------------------------------------
# Tabla 19-2: calibres y diámetros de alambre (pulg)
# --------------------------------------------------------------------------
# Alambre de instrumentos musicales (ASTM A228)
GAG_MUSIC = [
    ("7/0", 0.003), ("6/0", 0.004), ("5/0", 0.005), ("4/0", 0.006),
    ("3/0", 0.007), ("2/0", 0.008), ("0", 0.009), ("1", 0.010), ("2", 0.011),
    ("3", 0.012), ("4", 0.013), ("5", 0.014), ("6", 0.016), ("7", 0.018),
    ("8", 0.020), ("9", 0.022), ("10", 0.024), ("11", 0.026), ("12", 0.029),
    ("13", 0.031), ("14", 0.033), ("15", 0.035), ("16", 0.037), ("17", 0.039),
    ("18", 0.041), ("19", 0.043), ("20", 0.045), ("21", 0.047), ("22", 0.049),
    ("23", 0.051), ("24", 0.055), ("25", 0.059), ("26", 0.063), ("27", 0.067),
    ("28", 0.071), ("29", 0.075), ("30", 0.080), ("31", 0.085), ("32", 0.090),
    ("33", 0.095), ("34", 0.100), ("35", 0.106), ("36", 0.112), ("37", 0.118),
    ("38", 0.124), ("39", 0.130), ("40", 0.138),
]
# Calibre U.S. para alambre de acero (W&M)
GAG_US = [
    ("7/0", 0.4900), ("6/0", 0.4615), ("5/0", 0.4305), ("4/0", 0.3938),
    ("3/0", 0.3625), ("2/0", 0.3310), ("0", 0.3065), ("1", 0.2830),
    ("2", 0.2625), ("3", 0.2437), ("4", 0.2253), ("5", 0.2070), ("6", 0.1920),
    ("7", 0.1770), ("8", 0.1620), ("9", 0.1483), ("10", 0.1350), ("11", 0.1205),
    ("12", 0.1055), ("13", 0.0915), ("14", 0.0800), ("15", 0.0720), ("16", 0.0625),
    ("17", 0.0540), ("18", 0.0475), ("19", 0.0410), ("20", 0.0348), ("21", 0.0317),
    ("22", 0.0286), ("23", 0.0258), ("24", 0.0230), ("25", 0.0204), ("26", 0.0181),
    ("27", 0.0173), ("28", 0.0162), ("29", 0.0150), ("30", 0.0140), ("31", 0.0132),
    ("32", 0.0128), ("33", 0.0118), ("34", 0.0104), ("35", 0.0095), ("36", 0.0090),
    ("37", 0.0085), ("38", 0.0080), ("39", 0.0075), ("40", 0.0070),
]
# Diámetros métricos preferidos (mm)
WIRE_MM = [13.0, 12.0, 11.0, 10.0, 9.0, 8.5, 8.0, 7.0, 6.5, 6.0, 5.5,
           5.0, 4.8, 4.5, 4.0, 3.8, 3.5, 3.0, 2.8, 2.5, 2.0, 1.8,
           1.6, 1.4, 1.2, 1.0, 0.90, 0.80, 0.70, 0.65, 0.60, 0.55,
           0.50, 0.45, 0.40, 0.35, 0.30, 0.28, 0.25, 0.22, 0.20, 0.18]

# --------------------------------------------------------------------------
# Esfuerzos admisibles aproximados (figs. 19-8..19-13, fuente Carlson)
# puntos ancla: (diámetro pulg, τ promedio psi)
# --------------------------------------------------------------------------
TAU_ANCORS = {
    "A227": [(0.031, 58e3), (0.125, 46e3), (0.250, 40e3)],
    "A228": [(0.010, 150e3), (0.055, 135e3), (0.108, 125e3), (0.148, 120e3)],
    "A229": [(0.062, 95e3), (0.148, 85e3), (0.250, 78e3)],
    "A231": [(0.055, 130e3), (0.0625, 145e3), (0.148, 130e3), (0.250, 120e3)],
    "A401": [(0.055, 140e3), (0.0625, 155e3), (0.148, 140e3), (0.250, 130e3)],
    "A313": [(0.031, 100e3), (0.107, 90e3), (0.250, 80e3)],
}
# Factor τ_ligero/τ_promedio (libro: A228 150/135 pág. 748; A231 170/145 pág. 750)
TAU_LIG = {"A227": 1.12, "A228": 1.111, "A229": 1.12,
           "A231": 1.172, "A401": 1.15, "A313": 1.12}

EXTREMOS_OPC = {
    1: ("escuadrados_rectificados", "escuadrados y rectificados (Na=N-2)"),
    2: ("escuadrados", "escuadrados sin rectificar (Na=N-1)"),
    3: ("planos_rectificados", "planos y rectificados (Na=N-1)"),
    4: ("planos", "planos sin rectificar (Na=N)"),
}


# --------------------------------------------------------------------------
# Conversión de unidades
# --------------------------------------------------------------------------
def mm_a_pulg(mm):
    return mm / 25.4


def pulg_a_mm(inch):
    return inch * 25.4


def psi_a_MPa(psi):
    return psi / 145.0377


def MPa_a_psi(MPa):
    return MPa * 145.0377


# --------------------------------------------------------------------------
# Fórmulas del capítulo 19
# --------------------------------------------------------------------------
def wahl(C):
    """(19-5) Factor Wahl: K = (4C-1)/(4C-4) + 0.615/C."""
    return (4 * C - 1) / (4 * C - 4) + 0.615 / C


def esfuerzo(F, Dm, Dw, K=None):
    """(19-4) τ = 8·K·F·Dm/(π·Dw³) = 8·K·F·C/(π·Dw²)."""
    C = Dm / Dw
    if K is None:
        K = wahl(C)
    return 8 * K * F * Dm / (PI * Dw ** 3)


def deflexion(F, Dm, Dw, Na, G):
    """(19-6) f = 8·F·Dm³·Na/(G·Dw⁴)."""
    return 8 * F * Dm ** 3 * Na / (G * Dw ** 4)


def constante(F1, L1, F2, L2):
    """(19-1a) k = (F1-F2)/(L2-L1)  entre dos estados."""
    return (F1 - F2) / (L2 - L1)


def Na_desde_k(k, Dm, Dw, G):
    """(19-8) Na = G·Dw/(8·k·C³)  (diseño: se conoce k)."""
    C = Dm / Dw
    return G * Dw / (8 * k * C ** 3)


def Dw_tentativo(Fo, Dm, tau_d, K=1.2):
    """(19-7) Dw = ∛(3.06·Fo·Dm/τd)  (K nominal 1.2, pág. 750)."""
    return (3.06 * Fo * Dm / tau_d) ** (1 / 3)


def espiras_activas(N, extremos):
    """Na según tipo de extremos (pág. 739)."""
    if extremos == "escuadrados_rectificados":
        return N - 2
    if extremos in ("escuadrados", "planos_rectificados"):
        return N - 1
    return N  # planos


def Ls_de(Dw, N, extremos):
    """Longitud comprimida. Escuadrados: Ls = Dw·N (pág. 747)."""
    if extremos == "planos":
        return Dw * (N + 1)
    return Dw * N


def angulo_paso(p, Dm):
    """(19-2) Ángulo de paso α = arctan(p/(π·Dm)), en grados."""
    return math.degrees(math.atan(p / (PI * Dm)))


def pandeo(Lf, Dm, fo_op, extremos):
    """Chequeo de pandeo, Fig. 19-15.
    Curva A: empotrados (escuadrados y rectificados); B: fijo-articulado;
    C: ambos articulados. A: sin pandeo si Lf/Dm < 5.2 (pág. 748);
    crítico ≈ 0.20·(8/(Lf/Dm))² (anclado al ejemplo de la pág. 746);
    B y C se escalan con la relación de columnas de Euler (4:2:1)."""
    ratio = Lf / Dm
    curva = {"escuadrados_rectificados": "A", "escuadrados": "A",
             "planos_rectificados": "B", "planos": "C"}.get(extremos, "A")
    if curva in ("A", "B") and ratio < 5.2:
        return ratio, curva, None, True
    crit = 0.20 * (8.0 / ratio) ** 2
    if curva == "B":
        crit /= 2.0
    elif curva == "C":
        crit /= 4.0
    return ratio, curva, crit, (fo_op / Lf) <= crit


def holgura(Lo, Ls, Na, Dw):
    """Holgura de espira cc = (Lo - Ls)/Na; debe ser > Dw/10 (pág. 740)."""
    cc = (Lo - Ls) / Na
    return cc, cc > Dw / 10


def tau_adm(material, Dw_pulg, servicio="promedio"):
    """Esfuerzo cortante admisible (psi) aproximado (figs. 19-8..19-13)."""
    pts = TAU_ANCORS[material]
    if Dw_pulg <= pts[0][0]:
        tau = pts[0][1]
    elif Dw_pulg >= pts[-1][0]:
        tau = pts[-1][1]
    else:
        tau = pts[-1][1]
        for (d1, t1), (d2, t2) in zip(pts, pts[1:]):
            if d1 <= Dw_pulg <= d2:
                tau = t1 + (t2 - t1) * (Dw_pulg - d1) / (d2 - d1)
                break
    if servicio == "ligero":
        return tau * TAU_LIG.get(material, 1.12)
    if servicio == "severo":
        return tau * 0.85
    return tau


def calibre_cercano(tabla, Dw):
    """(calibre, diámetro) del tamaño estándar igual o inmediatamente mayor."""
    for g, d in tabla:
        if d >= Dw - 1e-9:
            return g, d
    return tabla[-1]


def _fmt(titulo):
    print("\n" + "=" * 66)
    print(titulo)
    print("=" * 66)


# --------------------------------------------------------------------------
# ANÁLISIS GENERAL: dado un resorte (geometría + material + carga),
# calcula todas sus características (cap. 19).
# --------------------------------------------------------------------------
def analizar(material=None, G=None, Dw=None, Dm=None, DE=None, DI=None,
             N=None, Na=None, extremos="escuadrados_rectificados",
             Lf=None, Fo=None, Lo=None, Fi=None, Li=None,
             tau_adm_psi=None, servicio="promedio", unidades="pulg",
             imprimir=True):
    """Análisis general de un resorte helicoidal de compresión.

    Parámetros (al menos Dw y un diámetro Dm/DE/DI, y Na o N):
      material   clave ASTM (A227, A228, A229, A231, A401, A313)
      G          módulo en cortante propio (en la unidad de presión usada);
                 si se da material, G se toma de la tabla 19-4
      Dw         diámetro del alambre
      Dm|DE|DI   diámetro medio / exterior / interior (cualquiera de los tres)
      N          número total de espiras (con extremos, define Na)
      Na         número de espiras activas (alternativa a N)
      extremos   "escuadrados_rectificados" | "escuadrados" |
                 "planos_rectificados" | "planos"
      Lf         longitud libre (opcional; permite Lo, pandeo, holgura)
      Fo, Lo     fuerza y longitud de operación (Lo=0 -> se calcula con Lf)
      Fi, Li     segundo estado (fuerza/longitud instalada), opcional
      tau_adm_psi esfuerzo admisible propio (psi), opcional
      servicio   "ligero" | "promedio" | "severo"
      unidades   "pulg" (lb, pulg, psi) o "mm" (N, mm, MPa)
    Devuelve un dict con todos los resultados. Las unidades de salida
    coinciden con las de entrada (τ en psi si unidades="pulg", MPa si "mm").
    """
    si = unidades == "mm"
    # --- módulo G ---------------------------------------------------------
    if G is None:
        if material is None:
            raise ValueError("Indica material (p. ej. A228) o G propio.")
        G = (G_GPA[material] * 1000.0) if si else G_PSI[material]
    # --- diámetros --------------------------------------------------------
    if Dw is None:
        raise ValueError("Falta Dw.")
    if Dm is None:
        if DE is not None:
            Dm = DE - Dw
        elif DI is not None:
            Dm = DI + Dw
        else:
            raise ValueError("Falta Dm, DE o DI.")
    DE_c, DI_c = Dm + Dw, Dm - Dw
    C = Dm / Dw
    K = wahl(C)
    # --- espiras ----------------------------------------------------------
    if Na is None:
        if N is None:
            raise ValueError("Falta N o Na.")
        Na = espiras_activas(N, extremos)
    # --- constante y longitudes -------------------------------------------
    k = G * Dw ** 4 / (8 * Dm ** 3 * Na)          # desde geometría (19-6)
    # medidas directas del usuario tienen prioridad sobre la geometría
    if Fo is not None and Fi is not None and Li is not None and Lo is not None:
        k = constante(Fo, Lo, Fi, Li)              # (19-1a)
    elif Fo is not None and Lo is not None and Lf is not None:
        k = Fo / (Lf - Lo)                         # (19-1b)
    if Fo is not None and Lo is None and Lf is not None:
        Lo = Lf - Fo / k
    N_c = N if N is not None else Na + 2           # sólo informativo
    Ls = Ls_de(Dw, N_c if N is not None else Na + 2, extremos)
    # --- cargas -----------------------------------------------------------
    res = {"Dw": Dw, "Dm": Dm, "DE": DE_c, "DI": DI_c, "C": C, "K": K,
           "Na": Na, "N": N_c, "G": G, "k": k, "Ls": Ls, "Lf": Lf,
           "unidades": unidades, "material": material, "extremos": extremos}
    if Fo is not None:
        f_o = Fo / k
        tau_o = esfuerzo(Fo, Dm, Dw, K)
        res.update(Fo=Fo, fo=f_o, tau_o=tau_o, Lo=(Lf - f_o) if Lf else Lo)
        if Lf:
            res["Lo"] = Lf - f_o
    if Fi is not None:
        f_i = Fi / k
        tau_i = esfuerzo(Fi, Dm, Dw, K)
        res.update(Fi=Fi, fi=f_i, tau_i=tau_i,
                   Li=(Lf - f_i) if Lf else Li)
    if Lf:
        Fs = k * (Lf - Ls)
        tau_s = esfuerzo(Fs, Dm, Dw, K)
        res.update(Fs=Fs, tau_s=tau_s, fo_solid=Lf - Ls)
        ratio, curva, crit, ok = pandeo(Lf, Dm, res.get("fo", Fs / k),
                                        extremos)
        res.update(pandeo_ratio=ratio, pandeo_curva=curva,
                   pandeo_crit=crit, pandeo_ok=ok)
        if res.get("Lo") is not None:
            cc, cc_ok = holgura(res["Lo"], Ls, Na, Dw)
            res.update(cc=cc, cc_ok=cc_ok)
        p = (Lf - 2 * Dw) / Na if extremos == "escuadrados_rectificados" \
            else None
        if p:
            res.update(paso=p, alfa=angulo_paso(p, Dm))
    res["D_hueco_min"] = DE_c + Dw / 10
    # --- esfuerzo admisible ----------------------------------------------
    if tau_adm_psi is None and material:
        tau_adm_psi = tau_adm(material, mm_a_pulg(Dw) if si else Dw, servicio)
    res["tau_adm_psi"] = tau_adm_psi
    if tau_adm_psi is not None:
        res["tau_adm"] = psi_a_MPa(tau_adm_psi) if si else tau_adm_psi
    if imprimir:
        _imprimir_analisis(res)
    return res


def _imprimir_analisis(r):
    si = r["unidades"] == "mm"
    uL = "mm" if si else "pulg"
    uF = "N" if si else "lb"
    uT = "MPa" if si else "psi"

    def p3(x):
        return f"{x:.3f}"

    _fmt("ANÁLISIS DE RESORTE HELICOIDAL DE COMPRESIÓN")
    print(f"Geometría: Dw={r['Dw']:.4g} {uL}, Dm={r['Dm']:.4g} {uL}, "
          f"DE={r['DE']:.4g} {uL}, DI={r['DI']:.4g} {uL}")
    print(f"           Na={r['Na']:.2f}, Ls={r['Ls']:.3f} {uL}, "
          f"extremos: {r['extremos'].replace('_', ' ')}")
    if r.get("Lf"):
        print(f"           Lf={r['Lf']:.3f} {uL}")
    print(f"\nC = Dm/Dw = {r['C']:.2f};   K (19-5) = {r['K']:.3f}")
    print(f"k (desde geometría) = G·Dw⁴/(8·Dm³·Na) = {r['k']:.3f} {uF}/{uL}")
    if r.get("fo") is not None:
        print(f"\nA la carga Fo = {r['Fo']:.3f} {uF}:")
        print(f"  fo (19-6) = 8·Fo·Dm³·Na/(G·Dw⁴) = {r['fo']:.3f} {uL}")
        if r.get("Lf"):
            print(f"  Lo = Lf - fo = {r['Lo']:.3f} {uL}")
        print(f"  τo (19-4) = 8·K·Fo·Dm/(π·Dw³) = {r['tau_o']:,.0f} {uT}")
    if r.get("Fi") is not None:
        print(f"\nA la fuerza instalada Fi = {r['Fi']:.3f} {uF}:")
        print(f"  fi = {r['fi']:.3f} {uL};  Li = {r['Li']:.3f} {uL};  "
              f"τi = {r['tau_i']:,.0f} {uT}")
    if r.get("Fs") is not None:
        print(f"\nEn longitud comprimida Ls = {r['Ls']:.3f} {uL}:")
        print(f"  Fs = k(Lf - Ls) = {r['Fs']:.2f} {uF};  "
              f"τs = {r['tau_s']:,.0f} {uT}")
    if r.get("paso"):
        print(f"\nPaso p = {r['paso']:.3f} {uL};  α (19-2) = "
              f"{r['alfa']:.1f}° (< 12° recomendado)")
    if r.get("tau_adm") is not None:
        print(f"\nEsfuerzo admisible aprox.: {r['tau_adm']:,.0f} {uT}")
        if r.get("tau_o") is not None:
            ok = r["tau_o"] <= r["tau_adm"]
            print(f"  τo vs τd  -> {'ACEPTABLE' if ok else 'NO ACEPTABLE'}")
        if r.get("tau_s") is not None:
            ok = r["tau_s"] <= r["tau_adm"] * TAU_LIG.get(
                r["material"] or "A227", 1.12)
            print(f"  τs vs τmáx (ligero) -> "
                  f"{'ACEPTABLE' if ok else 'NO ACEPTABLE'}")
    if r.get("pandeo_ok") is not None:
        crit = ("sin pandeo (Lf/Dm < 5.2)" if r["pandeo_crit"] is None else
                f"crítico fo/Lf ≈ {r['pandeo_crit']:.3f}")
        print(f"\nPandeo (fig. 19-15, curva {r['pandeo_curva']}): "
              f"Lf/Dm = {r['pandeo_ratio']:.2f}, {crit} -> "
              f"{'OK' if r['pandeo_ok'] else 'RIESGO DE PANDEO'}")
    if r.get("cc") is not None:
        print(f"Holgura de espira cc = (Lo-Ls)/Na = {r['cc']:.4f} {uL} > "
              f"Dw/10 = {r['Dw'] / 10:.4f} {uL} -> "
              f"{'OK' if r['cc_ok'] else 'NO OK'}")
    print(f"Diámetro mínimo del orificio: DE + Dw/10 = "
          f"{r['D_hueco_min']:.3f} {uL}")


# --------------------------------------------------------------------------
# DISEÑO GENERAL: por tanteo (método 1, sección 19-6 del libro)
# --------------------------------------------------------------------------
def disenar(material, F1, L1, F2, L2, extremos="escuadrados_rectificados",
           servicio="promedio", unidades="pulg", hueco=None, DI_min=None,
           Dm_estimado=None, G=None, imprimir=True):
    """Diseño general de un resorte helicoidal de compresión.

    Parámetros:
      F1, L1   fuerza y longitud del primer estado (p. ej. instalada)
      F2, L2   fuerza y longitud del segundo estado (p. ej. operación;
               debe ser la fuerza máxima Fo)
      material clave ASTM;  servicio: ligero | promedio | severo
      hueco    diámetro del orificio disponible (opcional)
      DI_min   diámetro interior mínimo exigido (opcional)
      Dm_estimado  Dm tentativo del usuario (opcional; si no, se deduce
               del hueco/DI_min o se escanea)
    Devuelve dict con el diseño ( Dw, Dm, DE, DI, Na, N, Ls, k, Lf,
    esfuerzos y verificaciones ).
    """
    si = unidades == "mm"
    G = (G_GPA[material] * 1000.0) if si else G_PSI[material]
    # estado de operación = fuerza máxima
    if F2 < F1:
        F1, F2, L1, L2 = F2, F1, L2, L1
    Fo, Lo, Fi, Li = F2, L2, F1, L1
    k = constante(Fo, Lo, Fi, Li)                    # (19-1a)
    Lf = Li + Fi / k                                 # pág. 749
    fo = Fo / k
    if fo >= Lf:
        raise ValueError("Datos inconsistentes: fo ≥ Lf.")
    # candidatos a Dw (de fino a grueso: busca el más esbelto factible)
    if si:
        candidatos = sorted(d for d in WIRE_MM if d >= 0.2)
    else:
        tabla = GAG_MUSIC if material == "A228" else GAG_US
        candidatos = sorted(d for _, d in tabla)
    Dw_t = None
    resultados = []

    def _seleccionar_N(Na_t, Dw, Dm):
        """Mejor N entero: k_i ≈ k con Ls < Lo y holgura suficiente."""
        mejor = None
        for Nn in range(math.ceil(Na_t) + 1, math.ceil(Na_t) + 4):
            Na_i = espiras_activas(Nn, extremos)
            k_i = G * Dw ** 4 / (8 * Dm ** 3 * Na_i)
            Ls_i = Ls_de(Dw, Nn, extremos)
            if Ls_i >= Lo:
                continue
            cc_i = (Lo - Ls_i) / Na_i
            if cc_i <= Dw / 10:
                continue
            err = abs(k_i - k) / k
            if mejor is None or err < mejor[0]:
                mejor = (err, Nn, Na_i, k_i, Ls_i, cc_i)
        return mejor

    for Dw in candidatos:
        tau_d_psi = tau_adm(material, mm_a_pulg(Dw) if si else Dw, servicio)
        tau_d = psi_a_MPa(tau_d_psi) if si else tau_d_psi
        # límites de Dm: C∈[5,12]; hueco acota por arriba; DI_min por abajo
        Dm_lo = 5.0 * Dw * 1.001
        if DI_min is not None:
            Dm_lo = max(Dm_lo, DI_min + Dw + (0.5 if si else 0.02))
        Dm_hi = 12.0 * Dw
        if hueco is not None:
            Dm_hi = min(Dm_hi, hueco - Dw - Dw / 10)
        if Dm_estimado is not None:
            dms = [Dm_estimado]
        else:
            if Dm_hi < Dm_lo:
                continue
            paso = 0.1 if si else 0.01
            n = int((Dm_hi - Dm_lo) / paso) + 1
            dms = [Dm_lo + i * paso for i in range(max(n, 1))]
        for Dm in dms:
            if Dm <= Dw:
                continue
            C = Dm / Dw
            if C < 5 or C > 12:
                continue
            K = wahl(C)
            tau = esfuerzo(Fo, Dm, Dw, K)
            if tau > tau_d:
                if Dw_t is None:
                    Dw_t = Dw_tentativo(Fo, Dm, tau_d)
                continue
            Na = G * Dw ** 4 / (8 * k * Dm ** 3)     # (19-8) equivalente
            sel = _seleccionar_N(Na, Dw, Dm)
            if sel is None:
                continue
            err, N, Na_i, k_i, Ls_i, cc_i = sel
            if Dw_t is None:
                Dw_t = Dw_tentativo(Fo, Dm, tau_d)
            # τs con los valores finales (N entero) y chequeo τs ≤ τmáx
            tau_s_c = esfuerzo(k_i * (Lf - Ls_i), Dm, Dw, K)
            solid_ok = tau_s_c <= tau_d * TAU_LIG.get(material, 1.12)
            resultados.append(dict(Dw=Dw, Dm=Dm, C=C, K=K, tau=tau,
                                   tau_d=tau_d, Na=Na, N=N, Na_i=Na_i,
                                   Ls=Ls_i, cc=cc_i, k_real=k_i, err_k=err,
                                   tau_s=tau_s_c, solid_ok=solid_ok))
            break                                   # mejor Dm para este Dw
    if not resultados:
        raise ValueError("Sin solución factible: revisar espacio, cargas "
                         "o material.")
    # preferir el primer Dw que además cumpla τs ≤ τmáx a sólido;
    # si ninguno, usar el primero factible (el informe avisará)
    con_solid = [t for t in resultados if t["solid_ok"]]
    elegido = con_solid[0] if con_solid else resultados[0]
    Dw, Dm = elegido["Dw"], elegido["Dm"]
    C, K = elegido["C"], elegido["K"]
    Na_t, N, Na_i = elegido["Na"], elegido["N"], elegido["Na_i"]
    k_i, Ls_i, cc_i, err = (elegido["k_real"], elegido["Ls"],
                            elegido["cc"], elegido["err_k"])
    tau_o, tau_s = elegido["tau"], elegido["tau_s"]
    DE, DI = Dm + Dw, Dm - Dw
    Fs = k_i * (Lf - Ls_i)
    # valores con Na fraccional (19-8), como los problemas del libro
    Ls_t = Ls_de(Dw, Na_t + 2, extremos)
    Fs_t = k * (Lf - Ls_t)
    tau_s_t = esfuerzo(Fs_t, Dm, Dw, K)
    ratio, curva, crit, ok = pandeo(Lf, Dm, fo, extremos)
    res = dict(material=material, servicio=servicio, extremos=extremos,
               unidades=unidades, k=k, Lf=Lf, Fi=Fi, Li=Li, Fo=Fo, Lo=Lo,
               Dw=Dw, Dm=Dm, DE=DE, DI=DI, C=C, K=K, Na=Na_i, N=N,
               Na_teo=Na_t, Ls=Ls_i, k_real=k_i, err_k=err, cc=cc_i,
               tau_o=tau_o, tau_s=tau_s, Fs=Fs, tau_d=elegido["tau_d"], G=G,
               Ls_teo=Ls_t, Fs_teo=Fs_t, tau_s_teo=tau_s_t,
               solid_ok=elegido["solid_ok"], Dw_t=Dw_t, pandeo_ratio=ratio,
               pandeo_curva=curva, pandeo_ok=ok, hueco=hueco, DI_min=DI_min,
               tabla_tanteo=resultados)
    if imprimir:
        _imprimir_diseno(res)
    return res


def _imprimir_diseno(r):
    si = r["unidades"] == "mm"
    uL = "mm" if si else "pulg"
    uF = "N" if si else "lb"
    uT = "MPa" if si else "psi"

    _fmt("DISEÑO DE RESORTE HELICOIDAL DE COMPRESIÓN (método 1, secc. 19-6)")
    print(f"Material: ASTM {r['material']};  servicio: {r['servicio']};  "
          f"G = {r['G']:.4g} {uT}")
    print(f"Requisitos: Fo = {r['Fo']:.3f} {uF} @ Lo = {r['Lo']:.3f} {uL};  "
          f"Fi = {r['Fi']:.3f} {uF} @ Li = {r['Li']:.3f} {uL}")
    print(f"k (19-1a) = {r['k']:.3f} {uF}/{uL};  Lf = {r['Lf']:.3f} {uL};  "
          f"fo = {r['Lf'] - r['Lo']:.3f} {uL}")
    print(f"Dw tentativo (19-7) ≈ {r['Dw_t']:.4f} {uL}")
    print("\nTanteo (Dw estándar; el primero factible se selecciona):")
    for t in r["tabla_tanteo"][:6]:
        print(f"  Dw={t['Dw']:.3g} -> Dm={t['Dm']:.3g}, C={t['C']:.1f}, "
              f"τ={t['tau']:,.0f} {uT} ≤ τd={t['tau_d']:,.0f} {uT}, "
              f"Na={t['Na']:.1f}, Ls={t['Ls']:.2f} {uL}, "
              f"τs={t['tau_s']:,.0f} {uT}: "
              f"{'OK' if t['solid_ok'] else 'no (τs > τmáx)'}")
    print(f"\nSelección: Dw = {r['Dw']:.4g} {uL}, Dm = {r['Dm']:.4g} {uL}, "
          f"C = {r['C']:.2f}, K (19-5) = {r['K']:.3f}")
    print(f"Na teórico (19-8) = {r['Na_teo']:.2f} -> N = {r['N']:.0f} "
          f"(Na = {r['Na']:.2f}), extremos {r['extremos'].replace('_', ' y ')}, "
          f"k real = {r['k_real']:.3f} {uF}/{uL} "
          f"(error {r['err_k'] * 100:.1f}%)")
    if not r.get("solid_ok"):
        print("  Nota: ningún candidato cumple τs ≤ τmáx a sólido; "
              "verificar figuras 19-8..19-13 del libro.")
    print(f"DE = {r['DE']:.3f} {uL};  DI = {r['DI']:.3f} {uL};  "
          f"Ls = {r['Ls']:.3f} {uL}")
    if abs(r["Na_teo"] - r["Na"]) > 0.01:
        print(f"Con Na teórico ({r['Na_teo']:.2f}), como en los problemas "
              f"del libro: Ls = {r['Ls_teo']:.3f} {uL}, "
              f"Fs = {r['Fs_teo']:.2f} {uF}, τs = {r['tau_s_teo']:,.0f} {uT}")
    print(f"\nVerificaciones:")
    print(f"  τo (19-4) = {r['tau_o']:,.0f} {uT}  vs  τd ≈ "
          f"{r['tau_d']:,.0f} {uT} -> "
          f"{'SEGURO' if r['tau_o'] <= r['tau_d'] else 'REVISAR'}")
    tmax = r['tau_d'] * TAU_LIG.get(r["material"], 1.12)
    print(f"  τs (en Ls) = {r['tau_s']:,.0f} {uT}  vs  τmáx ≈ "
          f"{tmax:,.0f} {uT} -> "
          f"{'SEGURO (sin fluencia)' if r['tau_s'] <= tmax else 'REVISAR'}")
    print(f"  Pandeo (fig. 19-15): Lf/Dm = {r['pandeo_ratio']:.2f} -> "
          f"{'OK' if r['pandeo_ok'] else 'RIESGO'}")
    print(f"  Holgura cc = {r['cc']:.4f} {uL} > Dw/10 = "
          f"{r['Dw'] / 10:.4f} {uL} -> OK")
    if r.get("hueco"):
        ok = r["DE"] + r["Dw"] / 10 <= r["hueco"]
        print(f"  Orificio {r['hueco']:.3f} {uL}: DE + Dw/10 = "
              f"{r['DE'] + r['Dw'] / 10:.3f} {uL} -> "
              f"{'CABE' if ok else 'NO CUADRA'}")
    if r.get("DI_min"):
        print(f"  DI = {r['DI']:.3f} {uL} vs DI mín. {r['DI_min']:.3f} {uL} "
              f"-> {'OK' if r['DI'] >= r['DI_min'] else 'NO'}")


# --------------------------------------------------------------------------
# Interfaz interactiva
# --------------------------------------------------------------------------
def _preguntar(texto, default=None, cast=float):
    s = input(texto).strip()
    if not s and default is not None:
        return default
    if cast is str:
        return s
    return cast(s)


def _elegir(texto, opciones):
    print(texto)
    for n, desc in opciones.items():
        if isinstance(desc, tuple):
            desc = f"{desc[1]}  ({desc[0]})"
        print(f"  [{n}] {desc}")
    while True:
        c = int(input("Opción: ").strip())
        if c in opciones:
            return c


def interactivo():
    print("=" * 66)
    print(" RESORTES HELICOIDALES DE COMPRESIÓN — Mott, cap. 19 (pág. 735+)")
    print("=" * 66)
    op = _elegir("¿Qué desea hacer?", {
        1: "Analizar un resorte existente (geometría -> características)",
        2: "Diseñar un resorte nuevo (cargas y espacio -> dimensiones)",
    })
    u = _elegir("Unidades", {1: "SI (mm, N, MPa)",
                             2: "Inglés (pulg, lb, psi)"})
    unidades, uL, uF = ("mm", "mm", "N") if u == 1 else ("pulg", "pulg", "lb")
    if op == 1:
        mat = input("Material ASTM (A227, A228, A229, A231, A401, A313) "
                    "o vacío para G propio: ").strip().upper()
        G = None
        if not mat:
            G = _preguntar(f"Ingrese G en {uL} (mismo sistema de presión): ")
            mat = None
        Dw = _preguntar(f"Diámetro del alambre Dw ({uL}): ")
        d_op = _elegir("Diámetro conocido", {1: "Dm (medio)", 2: "DE (exterior)",
                                             3: "DI (interior)"})
        dv = _preguntar("Valor: ")
        kw = {("Dm" if d_op == 1 else "DE" if d_op == 2 else "DI"): dv}
        e_op = _elegir("Espiras", {1: "N total + tipo de extremos",
                                   2: "Na (espiras activas) directo"})
        N = Na = None
        extremos = "escuadrados_rectificados"
        if e_op == 1:
            N = _preguntar("Número total de espiras N: ")
            x = _elegir("Extremos", EXTREMOS_OPC)
            extremos = EXTREMOS_OPC[x][0]
        else:
            Na = _preguntar("Na: ")
        Lf = _preguntar(f"Longitud libre Lf ({uL}) [0 si no se conoce]: ",
                        default=0.0) or None
        Fo = _preguntar(f"Fuerza de operación Fo ({uF}) [0 si no aplica]: ",
                        default=0.0) or None
        Fo = Fo if Fo else None
        Fo_arg = Fo
        x = _elegir("Servicio (para comparar esfuerzos)", {
            1: "ligero (estático / ≤10 000 ciclos)",
            2: "promedio (típico, ≤10⁶ ciclos)",
            3: "severo (>10⁶ ciclos, impacto)"})
        servicio = {1: "ligero", 2: "promedio", 3: "severo"}[x]
        analizar(material=mat or None, G=G, Dw=Dw, N=N, Na=Na,
                 extremos=extremos, Lf=Lf or None, Fo=Fo_arg,
                 servicio=servicio, unidades=unidades, **kw)
    else:
        mat = input("Material ASTM (A227, A228, A229, A231, A401, A313): "
                    ).strip().upper()
        F1 = _preguntar(f"Fuerza 1 ({uF}): ")
        L1 = _preguntar(f"Longitud 1 ({uL}): ")
        F2 = _preguntar(f"Fuerza 2 ({uF}): ")
        L2 = _preguntar(f"Longitud 2 ({uL}): ")
        x = _elegir("Servicio", {1: "ligero", 2: "promedio", 3: "severo"})
        servicio = {1: "ligero", 2: "promedio", 3: "severo"}[x]
        x = _elegir("Extremos", EXTREMOS_OPC)
        extremos = EXTREMOS_OPC[x][0]
        x = _elegir("Restricción de espacio", {
            1: "Diámetro del orificio (hueco)",
            2: "Diámetro interior mínimo (DI)",
            3: "Diámetro medio estimado (Dm)",
            4: "Sin restricción"})
        hueco = DI_min = Dm_est = None
        if x == 1:
            hueco = _preguntar(f"Diámetro del hueco ({uL}): ")
        elif x == 2:
            DI_min = _preguntar(f"DI mínimo ({uL}): ")
        elif x == 3:
            Dm_est = _preguntar(f"Dm estimado ({uL}): ")
        disenar(material=mat, F1=F1, L1=L1, F2=F2, L2=L2,
                extremos=extremos, servicio=servicio, unidades=unidades,
                hueco=hueco, DI_min=DI_min, Dm_estimado=Dm_est)


def main():
    if len(sys.argv) > 1 and sys.argv[1] in ("-h", "--help"):
        print(__doc__)
        return
    if sys.stdin.isatty() or "--interactivo" in sys.argv:
        interactivo()
    else:
        # sin TTY (p. ej. tubería): muestra ayuda
        print(__doc__)


if __name__ == "__main__":
    main()
