#!/usr/bin/env python3
"""Pure-Python transcription of the MACS3 3.0.5 Poisson paths.

`MACS3.Signal.Prob` is Cython, so importing it needs a compiled MACS3. This
module is a line-for-line Python copy of the same functions, transcribed from
upstream commit c544319, so we can generate parity vectors for the Rust port
before (and independently of) a working MACS3 build.

Transcribed functions:
  * poisson_cdf / __poisson_cdf / __poisson_cdf_large_lambda
  * __poisson_cdf_Q / __poisson_cdf_Q_large_lambda
  * log10_poisson_cdf_P_large_lambda / log10_poisson_cdf_Q_large_lambda
  * poisson_cdf_inv / poisson_pdf
  * binomial_pdf / binomial_cdf / binomial_cdf_inv / pduplication

Every arithmetic step, comparison and rounding is preserved. Do not "improve"
anything here: this file is the specification, not a better implementation.
"""

from math import exp, log, log1p, fabs

LSTEP = 200
EXPTHRES = exp(LSTEP)
EXPSTEP = exp(-1 * LSTEP)
bigx = 20


def __poisson_cdf(k, a):
    nextcdf = exp(-1 * a)
    cdf = nextcdf
    for i in range(1, k + 1):
        lastcdf = nextcdf
        nextcdf = lastcdf * a / i
        cdf = cdf + nextcdf
    if cdf > 1.0:
        return 1.0
    return cdf


def __poisson_cdf_large_lambda(k, a):
    num_parts = int(a / LSTEP)
    lastexp = exp(-1 * (a % LSTEP))
    nextcdf = EXPSTEP
    cdf = 0
    num_parts -= 1
    for i in range(1, k + 1):
        lastcdf = nextcdf
        nextcdf = lastcdf * a / i
        cdf += nextcdf
        if nextcdf > EXPTHRES or cdf > EXPTHRES:
            if num_parts >= 1:
                cdf *= EXPSTEP
                nextcdf *= EXPSTEP
                num_parts -= 1
            else:
                cdf *= lastexp
                lastexp = 1
    for i in range(num_parts):
        cdf *= EXPSTEP
    cdf *= lastexp
    return cdf


def __poisson_cdf_Q(k, a):
    cdf = 0.0
    nextcdf = exp(-1 * a)
    for i in range(1, k + 1):
        lastcdf = nextcdf
        nextcdf = lastcdf * a / i
    i = k + 1
    while nextcdf > 0.0:
        lastcdf = nextcdf
        nextcdf = lastcdf * a / i
        cdf += nextcdf
        i += 1
    return cdf


def __poisson_cdf_Q_large_lambda(k, a):
    num_parts = int(a / LSTEP)
    lastexp = exp(-1 * (a % LSTEP))
    nextcdf = EXPSTEP
    cdf = 0.0
    num_parts -= 1
    for i in range(1, k + 1):
        lastcdf = nextcdf
        nextcdf = lastcdf * a / i
        if nextcdf > EXPTHRES:
            if num_parts >= 1:
                nextcdf *= EXPSTEP
                num_parts -= 1
            else:
                raise Exception("Unexpected error")
    i = k + 1
    while nextcdf > 0.0:
        lastcdf = nextcdf
        nextcdf = lastcdf * a / i
        cdf += nextcdf
        i += 1
        if nextcdf > EXPTHRES or cdf > EXPTHRES:
            if num_parts >= 1:
                cdf *= EXPSTEP
                nextcdf *= EXPSTEP
                num_parts -= 1
            else:
                cdf *= lastexp
                lastexp = 1
    for i in range(num_parts):
        cdf *= EXPSTEP
    cdf *= lastexp
    return cdf


def logspace_add(logx, logy):
    if logx > logy:
        return logx + log1p(exp(logy - logx))
    return logy + log1p(exp(logx - logy))


def log10_poisson_cdf_P_large_lambda(k, lbd):
    residue = 0
    ln_lbd = log(lbd)
    m = k
    sum_ln_m = 0
    for i in range(1, m + 1):
        sum_ln_m += log(i)
    logx = m * ln_lbd - sum_ln_m
    residue = logx
    while m > 1:
        m -= 1
        logy = logx - ln_lbd + log(m)
        pre_residue = residue
        residue = logspace_add(pre_residue, logy)
        if fabs(pre_residue - residue) < 1e-10:
            break
        logx = logy
    return round((residue - lbd) / 2.302585092994046, 5)


def log10_poisson_cdf_Q_large_lambda(k, lbd):
    residue = 0
    logx = 0
    ln_lbd = log(lbd)
    m = k + 1
    sum_ln_m = 0
    for i in range(1, m + 1):
        sum_ln_m += log(i)
    logx = m * ln_lbd - sum_ln_m
    residue = logx
    while True:
        m += 1
        logy = logx + ln_lbd - log(m)
        pre_residue = residue
        residue = logspace_add(pre_residue, logy)
        if fabs(pre_residue - residue) < 1e-5:
            break
        logx = logy
    return round((residue - lbd) / log(10), 5)


def poisson_cdf(n, lam, lower=False, log10=False):
    assert lam > 0.0
    if log10:
        if lower:
            return log10_poisson_cdf_P_large_lambda(n, lam)
        return log10_poisson_cdf_Q_large_lambda(n, lam)
    if lower:
        if lam > 700:
            return __poisson_cdf_large_lambda(n, lam)
        return __poisson_cdf(n, lam)
    if lam > 700:
        return __poisson_cdf_Q_large_lambda(n, lam)
    return __poisson_cdf_Q(n, lam)


def poisson_cdf_inv(cdf, lam, maximum=1000):
    assert lam < 740
    if cdf < 0 or cdf > 1:
        raise Exception("CDF must >= 0 and <= 1")
    elif cdf == 0:
        return 0
    newval = exp(-1 * lam)
    sum2 = newval
    for i in range(1, maximum + 1):
        sumold = sum2
        newval = newval * lam / i
        sum2 = sum2 + newval
        if sumold <= cdf and cdf <= sum2:
            return i
    return maximum


def binomial_pdf(x, a, b):
    if a < 1:
        return 0.0
    elif x < 0 or a < x:
        return 0.0
    elif b == 0:
        if x == 0:
            return 1.0
        return 0.0
    elif b == 1:
        if x == a:
            return 1.0
        return 0.0
    if x > a - x:
        p = 1 - b
        mn = a - x
        mx = x
    else:
        p = b
        mn = x
        mx = a - x
    pdf = 1
    t = 0
    for q in range(1, mn + 1):
        pdf *= (a - q + 1) * p / (mn - q + 1)
        if pdf < 1e-100:
            while pdf < 1e-3:
                pdf /= 1 - p
                t -= 1
        if pdf > 1e+100:
            while pdf > 1e3 and t < mx:
                pdf *= 1 - p
                t += 1
    for i in range(mx - t):
        pdf *= 1 - p
    return pdf


def _binomial_cdf_r(x, a, b):
    argmax = int(a * b)
    if x < 0:
        return 1
    elif a < x:
        return 0
    elif b == 0:
        return 0
    elif b == 1:
        return 1
    if x < argmax:
        seedpdf = binomial_pdf(argmax, a, b)
        pdf = seedpdf
        cdf = pdf
        for i in range(argmax - 1, x, -1):
            pdf /= (a - i) * b / (1 - b) / (i + 1)
            if pdf == 0.0:
                break
            cdf += pdf
        pdf = seedpdf
        i = argmax
        while True:
            pdf *= (a - i) * b / (1 - b) / (i + 1)
            if pdf == 0.0:
                break
            cdf += pdf
            i += 1
        cdf = min(1, cdf)
        return cdf
    else:
        pdf = binomial_pdf(x + 1, a, b)
        cdf = pdf
        i = x + 1
        while True:
            pdf *= (a - i) * b / (1 - b) / (i + 1)
            if pdf == 0.0:
                break
            cdf += pdf
            i += 1
        cdf = min(1, cdf)
        return cdf


def _binomial_cdf_f(x, a, b):
    argmax = int(a * b)
    if x < 0:
        return 0
    elif a < x:
        return 1
    elif b == 0:
        return 1
    elif b == 1:
        return 0
    if x > argmax:
        seedpdf = binomial_pdf(argmax, a, b)
        pdf = seedpdf
        cdf = pdf
        for i in range(argmax - 1, -1, -1):
            pdf /= (a - i) * b / (1 - b) / (i + 1)
            if pdf == 0.0:
                break
            cdf += pdf
        pdf = seedpdf
        for i in range(argmax, x):
            pdf *= (a - i) * b / (1 - b) / (i + 1)
            if pdf == 0.0:
                break
            cdf += pdf
        cdf = min(1, cdf)
        return cdf
    else:
        pdf = binomial_pdf(x, a, b)
        cdf = pdf
        for i in range(x - 1, -1, -1):
            pdf /= (a - i) * b / (1 - b) / (i + 1)
            if pdf == 0.0:
                break
            cdf += pdf
        cdf = min(1, cdf)
        return cdf


def binomial_cdf(x, a, b, lower=True):
    if lower:
        return _binomial_cdf_f(x, a, b)
    return _binomial_cdf_r(x, a, b)


def binomial_sf(x, a, b, lower=True):
    if lower:
        return 1.0 - _binomial_cdf_f(x, a, b)
    return 1.0 - _binomial_cdf_r(x, a, b)


def binomial_cdf_inv(cdf, a, b):
    if cdf < 0 or cdf > 1:
        raise Exception("CDF must >= 0 or <= 1")
    cdf2 = 0.0
    for x in range(0, a + 1):
        pdf = binomial_pdf(x, a, b)
        cdf2 = cdf2 + pdf
        if cdf < cdf2:
            return x
    return a


def pduplication(pmf, N_obs):
    import numpy as np
    n = pmf.shape[0]
    sf = np.float32(0.0)
    for p in pmf:
        sf = np.float32(sf + np.float32(binomial_sf(2, N_obs, float(p))))
    return np.float32(sf / np.float32(n))
