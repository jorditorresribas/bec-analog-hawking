"""
DIAGNOSTICO (no gate, sin commit, 2026-11-18) — complemento de diag_recati_norm.py (A): simetria del corchete
INTEGRADO en omega sobre la ventana de Hawking que usa density_correlation (lam=1, tanh w=2, L=40, N=4800).
Compara la lectura literal \"+ c.c.\" (2 Re B) con la simetrizacion del codigo (B + B^dag).
"""
import sys,io,contextlib,numpy as np
sys.path.insert(0,__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.dirname(__import__('os').path.abspath(__file__))),'solver'))
with contextlib.redirect_stdout(io.StringIO()): import correlations as C
L,N,NC=40.0,4800,200; idx=np.linspace(0,N-1,NC).astype(int)
om_all=np.linspace(0.004,0.6,150); dw=om_all[1]-om_all[0]
lit=np.zeros((NC,NC)); cod=np.zeros((NC,NC)); used=[]
for om in om_all:
    md={}
    for ch in ('u','d1','d2'):
        r=C.scatter_field(L,N,om,ch,2.0,1.0,'horizon')
        if r['Rin_neg'] is None or r['nchan']!=2: md=None; break
        u=r['u'][idx]; vv=r['v_mode'][idx]; md[ch]=dict(u=u,v=vv,r=u+vv)
    if md is None: continue
    used.append(om)
    B=sum(np.outer(np.conj(md[I]['v']),md[I]['r']) for I in ('u','d1'))+np.outer(np.conj(md['d2']['u']),md['d2']['r'])
    lit+=2*B.real*dw; cod+=np.real(B+np.conj(B.T))*dw
x=np.linspace(-L,L,N)[idx]; m=np.abs(x)<=25; M=np.ix_(m,m)
print(f"omegas usadas {len(used)}: [{used[0]:.3f}, {used[-1]:.3f}]  (resto fuera de la ventana de Hawking: omitidas por el codigo)")
print(f"literal (2 Re B integrado): asimetria max|L-L^T|/max|L| = {np.abs(lit-lit.T)[M].max()/np.abs(lit[M]).max():.2e}")
print(f"codigo  (B+B^dag integrado): asimetria = {np.abs(cod-cod.T)[M].max()/np.abs(cod[M]).max():.1e}")
print(f"max|codigo - literal|/max|literal| (|x|<=25) = {np.abs(cod-lit)[M].max()/np.abs(lit[M]).max():.2e}")
sym=0.5*(lit+lit.T); print(f"max|codigo - sim(literal)|/max = {np.abs(cod-sym)[M].max()/np.abs(lit[M]).max():.1e}")
