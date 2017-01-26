#!/usr/bin/env python

import linmod as lm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import LoadData
# from LoadData import annot
from LoadData import bestTrx
from LoadData import bestTrxToSymbol
from LoadData import group
from LoadData import transcriptAnnot
from LoadData import x
from LoadData import xg


## -----------------------------------------------------------------
theKeys = bestTrx + " " + bestTrx.index.to_series()
tkey = bestTrx
gkey = bestTrxToSymbol

fscores = lm.linmod(
    y = xg.transpose(),
    x = group,
    ref = "VIP"
)["F"]

def regroup(groupings, val):
    preout = (groupings == val)
    out = preout.copy()
    out.ix[preout] = val
    out.ix[~preout] = "other"
    return out

tscores = {
    g : lm.linmod(
        y = xg.transpose(),
        x = regroup(group, g)
    )["t"].ix[:, 0]
    for g in group.unique()
}

useKeys = theKeys.ix[(set(fscores.dropna().index) &
                      set(theKeys.index))]
for g in tscores:
    useKeys = useKeys.ix[(set(tscores[g].dropna().index) &
                          set(useKeys.index))]
useTrx = bestTrx[useKeys.index]

fscores.index = bestTrx[fscores.index]
for g in tscores:
    tscores[g].index = bestTrx[tscores[g].index]


## -----------------------------------------------------------------
out = pd.DataFrame({
    "F" : fscores.ix[useTrx],
    "t_exc" : tscores["excitatory"].ix[useTrx],
    "t_pv" : tscores["PV"].ix[useTrx],
    "t_vip" : tscores["VIP"].ix[useTrx]
})[["F", "t_exc", "t_pv", "t_vip"]]

out.to_csv("grcm38_hyp_tests.tsv", sep="\t", index=True, header=True)
