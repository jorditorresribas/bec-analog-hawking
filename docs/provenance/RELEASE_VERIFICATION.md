# Verificación de la versión pública (v1.4, 2026-11-17)

La versión pública es una copia reorganizada de la carpeta de trabajo original. Cambios de
código: solo rutas (cabecera de 1 línea para importar `solver/`; rutas absolutas → relativas
al repo en 4 gates de `residual_anatomy/gates/`). Ninguna línea de cálculo se modificó.

Re-ejecución en la estructura nueva (macOS 12, Python 3.9.6, numpy 2.0.2, scipy 1.13.1),
comparando con las salidas originales:

| Prueba | Comando | Resultado | Salidas vs originales |
|---|---|---|---|
| Benchmark FP v1.0 | `python3 benchmarks/fp_v1/FP_check_reproduction.py` | PASS (desviación máx. 0.0e+00) | `out/fp_benchmark.csv` idéntico byte a byte |
| E1 (Larré et al. 2012) | `python3 residual_anatomy/gates/gate_E1_lrcp.py` | INTEGRIDAD PASS · CLASIFICACIÓN PASS (46 s) | 3 CSV idénticos byte a byte |
| T2 (fondo consistente) | `python3 residual_anatomy/diag_T2_qp_consistent.py` | integridad PASS (62 s) | 2 CSV idénticos byte a byte |
| Test de Goldstone | `python3 residual_anatomy/diag_goldstone.py` | exit 0 (414 s) | CSV idéntico byte a byte |

Los manifiestos JSON de esa re-ejecución difieren de los originales solo en fecha, duración y
rutas de dependencias; se guardan en `release_rerun/` junto con el log. En `out/` se conservan
los manifiestos originales, que son el registro de las ejecuciones citadas en los documentos.
