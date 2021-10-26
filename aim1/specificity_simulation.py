import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from statsmodels import api
from scipy import stats
from scipy.optimize import minimize

# if using frequentist
y = np.abs(np.random.normal(0, 0.7, 1200))

def MLE_func(parameters):
    sigma = parameters
    ll = np.sum(stats.halfnorm.logpdf(y, 0, sigma))
    neg_ll = -1 * ll
    return neg_ll

mle_model = minimize(MLE_func, np.array([0.2]), method='L-BFGS-B')
mle_model.x

# if using bayesian
y = np.abs(np.random.normal(0, 0.7, 1200))
x = np.random.poisson(y.mean(),50)

import pymc3 as pm
import theano
with pm.Model() as m:
    sigma = pm.Uniform('sigma',lower=0,upper=1)
    nc = pm.HalfNormal('nc',sigma=sigma,observed=y)
    rate = pm.math.sum(nc)/len(y)
    c = pm.Poisson('c',mu=rate,observed=x)

with m:
    step = pm.NUTS()
    trace = pm.sample(500,step=step,return_inferencedata=False,cores=1)

import arviz as az
with m:
    az.plot_trace(trace)
    df = az.summary(trace,round_to=2)