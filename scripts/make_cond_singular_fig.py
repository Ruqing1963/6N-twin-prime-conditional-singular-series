#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build the obs-vs-prediction showdown figure (4 panels: 6dN = 30,42,60,210) from
css_plot_S10.csv (produced by cond_singular_series.py with MAXK=10).
The 210 panel shades the systematic high-omega second-order residual.
"""
import csv, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
rows=list(csv.DictReader(open('../data/css_plot_S10.csv')))
def series(sixd):
    om=[];obs=[];pred=[]
    for r in rows:
        if int(r['6dN'])==sixd:
            om.append(int(r['omega'])); obs.append(float(r['obs_r'])); pred.append(float(r['pred_r']))
    return np.array(om),np.array(obs),np.array(pred)

fig,axes=plt.subplots(1,4,figsize=(19,4.6),sharey=True)
panels=[(30,'#2ca25f','lockdown primes: 29,31'),
        (42,'#185FA5','flanked by primes 41,43'),
        (60,'#777777','no small-prime lock'),
        (210,'#c0392b','flanked by 209=11x19')]
for ax,(sixd,c,note) in zip(axes,panels):
    om,obs,pred=series(sixd)
    ax.plot(om,obs,'o-',color=c,lw=2.2,ms=7,label='observed $r$')
    ax.plot(om,pred,'s--',color=c,lw=1.6,ms=6,alpha=.65,label='first-order $\\mathfrak{S}_{2,\\omega}$')
    if sixd==210:
        ax.fill_between(om,pred,obs,where=(obs<pred),color='#c0392b',alpha=.18,label='2nd-order residual')
    ax.axhline(1,color='gray',ls=':',lw=1)
    ax.set_title(f'$6\\Delta N={sixd}$\n{note}',fontsize=10.5)
    ax.set_xlabel(r'$\omega_{>3}(N)$',fontsize=11)
    ax.grid(alpha=.25); ax.set_xticks([1,2,3,4,5,6])
    ax.legend(fontsize=8,loc='best')
axes[0].set_ylabel(r'relative preference $r(d\mid\omega)$',fontsize=11)
plt.suptitle('First-order conditional singular series $\\mathfrak{S}_{2,\\omega}(d)$ vs observation in $S_{10}$ '
             '(23,988,173 twin centres)',fontsize=13,y=1.04)
plt.figtext(0.5,-0.04,'Theory $=C_2(d)\\times$ lockdown-allowance over the stratum\'s real centre factors. '
            'Matches 30/42/60 and the direction of 210; the systematic high-$\\omega$ residual on 210 '
            '(shaded, stable across $S_9$/$S_{10}$) is the second-order open problem.',
            ha='center',fontsize=9,style='italic')
plt.tight_layout()
plt.savefig('fig_paper3_cond_singular.pdf',bbox_inches='tight')
plt.savefig('fig_paper3_cond_singular.png',dpi=160,bbox_inches='tight')
print("figure saved")
