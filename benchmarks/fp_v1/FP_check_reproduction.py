"""
Comprobador de reproducibilidad del benchmark FP (regresion numerica, NO validacion
fisica). Regenera la tabla principal y la compara con la referencia CONGELADA por
TOLERANCIA mixta (abs+rel), robusto a plataforma. NO modifica la referencia (solo lectura).

Comprueba:
  - mismas filas (mismos w): sin ausentes, sin adicionales, sin duplicadas;
  - todos los valores finitos (rechaza NaN/inf) en ambos CSV;
  - columnas de RESULTADO dentro de tolerancia: |b-a| <= ATOL + RTOL*|a|
    (RTOL=1e-4, ATOL=1e-12; el termino ATOL hace la comparacion segura si |a|~0);
  - columnas de DIAGNOSTICO (error de extrapolacion, pseudo-u) solo finitas + informe,
    sin fallar (pueden variar mas entre BLAS/LAPACK; no son el resultado cientifico).
Uso:   python3 FP_check_reproduction.py        (exit 0 = PASS, 1 = FAIL)
"""
import csv, subprocess, sys
RTOL, ATOL = 1e-4, 1e-12
REF = "reference/fp_benchmark_reference.csv"
OUT = "out/fp_benchmark.csv"
RESULT_COLS = ['kappa_local','cH','kloc_over_2pi','kbar_over_2pi','T_upstream','T_total',
               'Tup_over_kloc2pi','Tup_over_kbar2pi','Ttot_over_kloc2pi']
DIAG_COLS   = ['pseudo_u_max','fit_rms_resid','order_dispersion','T_up_linear','T_up_cubic']

def load(path):
    with open(path, 'r') as f:           # solo lectura; nunca se escribe la referencia
        rows = list(csv.DictReader(f))
    return rows

def finite(x):
    try: v = float(x)
    except (TypeError, ValueError): return False, None
    return (v == v and abs(v) != float('inf')), v  # v==v descarta NaN

fails = []
print("Regenerando out/fp_benchmark.csv ...")
subprocess.run([sys.executable, __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "FP_benchmark_table.py")], check=True,
               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

ref_rows, cur_rows = load(REF), load(OUT)

# --- integridad de filas: duplicados, ausentes, adicionales ---
def keys(rows, name):
    ks = []
    for r in rows:
        ok, v = finite(r.get('w'))
        if not ok: fails.append(f"[{name}] fila con w no numerico: {r.get('w')!r}"); continue
        ks.append(round(v, 9))
    dups = {k for k in ks if ks.count(k) > 1}
    if dups: fails.append(f"[{name}] w duplicados: {sorted(dups)}")
    return set(ks), ks
ref_set, _ = keys(ref_rows, 'ref'); cur_set, _ = keys(cur_rows, 'out')
missing = ref_set - cur_set; extra = cur_set - ref_set
if missing: fails.append(f"filas AUSENTES respecto a la referencia: w={sorted(missing)}")
if extra:   fails.append(f"filas ADICIONALES no en la referencia: w={sorted(extra)}")

refd = {round(float(r['w']),9): r for r in ref_rows}
curd = {round(float(r['w']),9): r for r in cur_rows}

# --- comparacion por columna ---
maxdev = 0.0; worst = None
for w in sorted(ref_set & cur_set):
    for c in RESULT_COLS + DIAG_COLS:
        okr, a = finite(refd[w].get(c)); okc, b = finite(curd[w].get(c))
        if not okr: fails.append(f"w={w:.0f} col={c}: valor NO finito en referencia ({refd[w].get(c)})")
        if not okc: fails.append(f"w={w:.0f} col={c}: valor NO finito en salida ({curd[w].get(c)})")
        if not (okr and okc): continue
        if c in RESULT_COLS:
            dev = abs(b-a)
            if dev > ATOL + RTOL*abs(a):
                fails.append(f"w={w:.0f} col={c}: |{b:.10g}-{a:.10g}|={dev:.2e} > tol")
            rel = dev/abs(a) if a != 0 else dev
            if rel > maxdev: maxdev, worst = rel, (w, c, a, b)

print(f"\nColumnas de resultado: max desviacion relativa vs referencia = {maxdev:.2e}  "
      f"(RTOL={RTOL:.0e}, ATOL={ATOL:.0e})")
if worst: print(f"  peor caso: w={worst[0]:.0f} col={worst[1]} ref={worst[2]:.10g} cur={worst[3]:.10g}")
print("Diagnosticos (pseudo_u, error extrapolacion): solo se exige finitud; no fallan la prueba.")

print("\nSenal cientifica (de la referencia): T_upstream/(kbar/2pi) por w:")
for w in sorted(ref_set):
    val = float(refd[w]['Tup_over_kbar2pi'])
    print(f"  w={w:.0f}: {val:.4f}  (desv. de 1 = {100*(val-1):+.2f}%)")

if fails:
    print(f"\nFAIL ({len(fails)} problema/s):")
    for m in fails: print(f"  - {m}")
    sys.exit(1)
print(f"\nPASS: reproduccion dentro de tolerancia; filas integras; sin NaN/inf. "
      f"Ruido numerico ({maxdev:.1e}) << senal fisica (0.5-0.75%).")
sys.exit(0)
