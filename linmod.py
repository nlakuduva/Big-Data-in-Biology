#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp
from scipy import stats


def dichotomize(x, interact=False, ref=None):
    if type(x) == type(pd.Series()):
        x = pd.DataFrame(x)
    levels = {cn : x[cn].unique() for cn in x.columns}
    if ref is None:
        ref = {}
    if "keys" not in dir(ref):
        ref = {x.columns[0] : ref}
    indicators = []
    for cn in x.columns:
        if cn not in ref:
            ref[cn] = levels[cn][0]
        refpos = (levels[cn] == ref[cn]).nonzero()[0][0]
        newinds = levels[cn][list(range(refpos)) +
                             list(range(refpos+1, len(levels[cn])))]
        if len(x.columns) > 1:
            newinds = [cn+"_"+ind for ind in newinds]
        indicators.extend(newinds)
    out = pd.DataFrame(np.zeros([x.shape[0], len(indicators)]),
                       index = x.index,
                       columns = indicators)
    for cn in x.columns:
        for level in levels[cn]:
            if level != ref[cn]:
                if len(x.columns) > 1:
                    out[cn+"_"+level] = (x[cn] == level) + 0.0
                else:
                    out[level] = (x[cn] == level) + 0.0
    if interact:
        cnout = out.columns.copy()
        for i in range(len(cnout)-1):
            cni = out.columns[i]
            for j in range(i+1, len(cnout)):
                cnj = out.columns[j]
                if np.sum(out[cni] * out[cnj]) > 0:
                    out[cni + ":" + cnj] = out[cni] * out[cnj]
    return out


def linmod(y, x, ridge=0, p=True,
           interact=False, ref=None, **kwargs):
    if not np.all(x.dtypes == "float64"):
        x = dichotomize(x, interact, ref)
    n = 1.0 * y.shape[0]
    ## y should have features in columns, sampling units in rows
    yorig = y
    xorig = x
    ## first remove means from y and x
    y = yorig.subtract(yorig.mean(axis=0))
    x = xorig.subtract(xorig.mean(axis=0))
    ## make sure two DataFrames are properly aligned
    y = y.ix[x.index]
    ## sum of squares for y
    ynorm = np.sqrt((y**2).sum(axis=0))
    ## ridge parameter matrix
    tau = np.diag(ridge * np.ones(x.shape[1]))
    ## hat matrix
    xt = x.transpose()
    hatMatrix = np.dot(np.linalg.inv(np.dot(xt, x) + (tau**2)), xt)
    yt = y.transpose()
    ## regression coefficients beta
    beta = pd.DataFrame(np.dot(hatMatrix, y),
                        columns=y.columns, index=x.columns)
    ## project y onto x to get predicted values yhat
    yhat = np.dot(x, beta)
    yhatnorm = np.sqrt((yhat**2).sum(axis=0))
    ## pearson correlation rho between y and yhat (magnitude)
    rho = (y * yhat).sum(axis=0) / (ynorm * yhatnorm)
    ## ...and sign...
    if beta.shape[0] == 1:
        rho = np.sign(beta.ix[0]) * rho
    df = n - 1.0 - x.shape[1]
    ## standard error
    se = np.sqrt((1.0 - rho**2) / df) * ynorm
    out = {"coef" : beta.transpose(),
           "rho" : rho,
           "se" : se,
           "resid" : y - yhat}
    ## calculate F and t statistics plus p values
    if p and ridge <= 1e-10:
        xnorm = np.sqrt((x**2).sum(axis=0))
        featse = {k : np.sqrt((1-rho**2) / df) * (ynorm / xnorm[k])
                  for k in xnorm.index}
        if x.shape[1] > 1:
            featrho = {x.columns[i] : linmod(x.ix[:, [i]],
                                             x.ix[:, list(range(i)) +
                                                  list(range(i+1, x.shape[1]))],
                                             ridge = ridge,
                                             p = False)["rho"][0]
                       for i in range(x.shape[1])}
            featse = {k : featse[k] / np.sqrt(1-featrho[k]**2)
                      for k in xnorm.index}
        featse = pd.DataFrame(featse)
        featse = featse[beta.index.values]
        nGroups = x.shape[1] + 1
        out["F"] = (1./x.shape[1]) * ((nGroups - n) +
                                      ((n - 1) *
                                       ((yorig.var(axis=0, ddof=1) /
                                         out["resid"].var(axis=0,
                                                          ddof=nGroups)))))
        out["p_F"] = pd.Series(stats.f.sf(out["F"], x.shape[1], df),
                               index = out["F"].index)
        out["t"] = beta.transpose() / featse
        out["p_t"] = pd.DataFrame(2 * stats.t.cdf(-abs(out["t"]), df),
                                  index=y.columns, columns=x.columns)
    return out


def regress(y, x, ridge=0, p=True, **kwargs):
    results = linmod(y, pd.DataFrame(x), ridge, p, **kwargs)
    out = pd.DataFrame({
        "coefficient" : results["coef"].ix[:, 0],
        "rho" : results["rho"],
        "rho^2" : results["rho"]**2,
        "se" : results["coef"].ix[:, 0] / results["t"].ix[:, 0],
        "p" : results["p_F"]
    })
    out["intercept"] = y.mean(axis=0) - out["coefficient"] * x.mean()
    out = out[["coefficient", "intercept", "rho", "rho^2", "se", "p"]]
    return out
