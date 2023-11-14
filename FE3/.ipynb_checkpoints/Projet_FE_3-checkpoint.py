#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 19:57:44 2023

@author: williamarnaud
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 17:37:06 2023
 
@author: S084339
"""
 
"""
https://www.kaggle.com/code/eneszvo/time-series-forecasting-es-arima-var
Vector Autoregression Model (VAR)
Similar to the AR model but using multiple time series. As opposite to previous models which are utilizing univariate time series, where the signal has only a single time-dependent variable, the VAR model using multivariate time series. Basically, it means that each signal depends not only on its past but also on the previous values of some other signals. The simplest model VAR(1) with two time series  
yt and xt can be written as
 
yt = a1 + b11 yt-1 + b12 xt-1 + eyt
xt = a2 + b21 yt-1 + b22 xt-1 + ext
 
where a and b are coefficients, and e is error.
 
VAR model requires stationarity of signals. However, when signals xt and yt are I(1) (non-stationary with the order of integration 1)
and if there is a o such that (yt - o xt) is stationary, xt and yt are cointegrated and a VEC model can be used.
 
Cointegration can be tested using augmented Engle-Granger cointegration test
where the null hypothesis assumes that there is no cointegration.
It means that if the p-value is less than 0.05, then we can reject a null hypothesis and say that signals are cointegrated.
 
"""

# CE SITE !!
# https://www.machinelearningplus.com/time-series/vector-autoregression-examples-python/
 
 
 
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.tsa.stattools import acf
from statsmodels.graphics.tsaplots import plot_acf
from statsmodels.tsa.arima.model import ARIMA
from statsmodels.graphics.tsaplots import plot_pacf
from statsmodels.graphics.tsaplots import plot_predict
import numpy as np
import os
from statsmodels.tsa.api import VAR, SVAR
from statsmodels.tsa.stattools import adfuller
from statsmodels.tsa.seasonal import seasonal_decompose
from statsmodels.tools.eval_measures import mse,rmse
from fredapi import Fred
from statsmodels.tsa.stattools import grangercausalitytests
from statsmodels.tsa.vector_ar.vecm import coint_johansen
import yfinance as yf
 
 
def get_FRED_series(ticker):
    FRED_API_KEY = "cc8c5ce8b5f37ae0255d927b30d50982"
    fred = Fred(api_key=FRED_API_KEY)
    data = fred.get_series(ticker)
    data = data.dropna()
    data = pd.DataFrame(data)
    data.index = pd.to_datetime(data.index, utc=True)
    return data
 
def get_yahoo_series(ticker):
    data = yf.Ticker(ticker).history(period='max')["Close"]
    data = data.dropna()
    data = pd.DataFrame(data)
    data.index = pd.to_datetime(data.index, utc=True)
    return data

x = get_yahoo_series('GOOG')

def ADF_reg_detail(series):
    """
   Parameters
    ----------
    series : TYPE
        DESCRIPTION.
 
    Returns
    -------
    None.
   
    Rappel
    -------
    Dans les regs : x1 c'est le Rho, const le coef, x2 le coef trend!
    p-value : ProbabilitÃ© de H0, si p-value <0.05 rejeter H0, sinon H0 ok
 
    """
    reg = adfuller(series, regression='ct', regresults=True, store=True)
    rho_coef = reg[3].resols.params[0]
    rho_pv = reg[3].resols.pvalues[0]
    cst_coef = reg[3].resols.params[1]
    cst_pv = reg[3].resols.pvalues[1]
    trend_coef = reg[3].resols.params[2]
    trend_pv = reg[3].resols.pvalues[2]
    print('====================')
    print('Regression CT : ')
    print('====================')
    print('Rho coef :', rho_coef)
    print('Rho p-value :', rho_pv)
    print('Constant coef :', cst_coef)
    print('Constant p-value :', cst_pv)
    print('Trend coef :', cst_coef)
    print('Trend p-value :', cst_pv)
    print("")
    print("--------------------")
    if trend_pv < 0.05 :
        print('Strong evidence against the null hypothesis (Trend coefficient =0 )')
        print('Reject the null hypothesis')
        print('The series has a trend')
        if rho_pv <0.05 :
            print("rho diff 0, the series has no UR, but a deterministic trend")
        else :
            print("rho=0, the series has a UR, and has a deterministic trend")
        return
    if trend_pv > 0.05 :
        print('Weak evidence against the null hypothesis (Trend coefficient =0 )')
        print('Fail to reject the null hypothesis')
        print('The series has no trend')
    print("")
    print("")
    reg = adfuller(series, regression='c', regresults=True, store=True)
    rho_coef = reg[3].resols.params[0]
    rho_pv = reg[3].resols.pvalues[0]
    cst_coef = reg[3].resols.params[1]
    cst_pv = reg[3].resols.pvalues[1]
    print('====================')
    print('Regression C : ')
    print('====================')
    print('Rho coef :', rho_coef)
    print('Rho p-value :', rho_pv)
    print('Constant coef :', cst_coef)
    print('Constant p-value :', cst_pv)
    print("")
    print("--------------------")
    if cst_pv < 0.05 :
        print('Strong evidence against the null hypothesis (constant coefficient = 0)')
        print('Reject the null hypothesis')
        print('The series has a constant')
        if rho_pv < 0.05 :
            print("rho diff 0, the series has no UR, but a constant and no trend.")
        else :
            print("rho=0, the series has a UR, and has a constant, but no trend")
        return
    if cst_pv > 0.05 :
        print('Weak evidence against the null hypothesis (constant coefficient =0)')
        print('Fail to reject the null hypothesis')
        print('The series has no constant')
    print("")
    print("")
    reg = adfuller(series, regression='n', regresults=True, store=True)
    rho_coef = reg[3].resols.params[0]
    rho_pv = reg[3].resols.pvalues[0]
    print('====================')
    print('Regression NC : ')
    print('====================')
    print('Rho coef :', rho_coef)
    print('Rho p-value :', rho_pv)
    print("")
    print("--------------------")
    if rho_pv < 0.05 :
        print('Strong evidence against the null hypothesis (UR coefficient = 0)')
        print('Reject the null hypothesis')
        print('The series is stationary, no UR, without a constant, nor trend.')
    if rho_pv > 0.05 :
        print('Weak evidence against the null hypothesis (UR coefficient =0)')
        print('Fail to reject the null hypothesis')
        print('Data has a unit root and is non stationary')
    print("")
    print("")
 
 
 
def johansen_coint_test(df):
    """
    Eigenvalues: In the Johansen cointegration test, the eigenvalues represent the characteristic roots of a matrix formed during the test.
    These eigenvalues provide valuable information about the cointegration structure of the time series.
   
    Sorting Eigenvalues: The eigenvalues are typically sorted in descending order, which means the largest eigenvalue is listed first.
   
    Comparing Eigenvalues to Critical Values: To determine the number of cointegrating relationships, you need to compare the eigenvalues to critical
    values from statistical tables or generated using Monte Carlo simulations.
    These critical values depend on the sample size, significance level, and the number of time series.
   
    Interpretation based on Eigenvalues and Critical Values:        H0 : r=0, H1: r=1
 
        * Case 1: If all eigenvalues are less than the critical values, this suggests that there are no cointegrating relationships.
        In other words, the time series are not cointegrated.
        * Case 2: If at least one eigenvalue exceeds the critical values, it suggests the presence of cointegration.
            - If only one eigenvalue exceeds the critical values, it indicates the existence of a single cointegrating relationship.
            - If multiple eigenvalues exceed the critical values, you'll need to examine the pattern of eigenvalues.
            The number of eigenvalues that exceed the critical values determines the rank of cointegration.
                - If all but one eigenvalue exceeds the critical values, there is one cointegrating relationship.
                - If all but two eigenvalues exceed the critical values, there are two cointegrating relationships.
                - And so on...
        * Case 3: If all eigenvalues exceed the critical values, it suggests that all variables are cointegrated, but the test may not be able to identify the exact number of cointegrating relationships.
   
    """
    result = coint_johansen(df, det_order=0, k_ar_diff=1)
    print("===================================================")
    print("Trace Statistic:")
    print("--------------------------------------------------")
    print(result.trace_stat)  # Trace statistic
    print("\nCritical Values at 90%-95%-99% Confidence Level:")
    print("--------------------------------------------------")
    print(result.trace_stat_crit_vals)  # Critical values
    print("===================================================")
    print("\nEigenvalues:")
    print("--------------------------------------------------")
    print(result.eig)  # Eigenvalue
    print("\n Max Eigenvalues Statistic:")
    print("--------------------------------------------------")
    print(result.max_eig_stat)  # Eigenvalue
    print("\nCritical Values at 90%-95%-99% Confidence Level:")
    print("--------------------------------------------------")
    print(result.max_eig_stat_crit_vals)  # Critical values
    print("===================================================")
 
 
 
def select_VAR_order(df, nb_orders_tested=50):
    l = []
    for p in range(nb_orders_tested):
        model = VAR(df)
        results = model.fit(p)
        l.append((p, results.aic, results.bic))
    l = pd.DataFrame(l)
    l.columns = ["p", "AIC", "BIC"]
    l = l.sort_values(by="AIC",ascending=True)
    print(l[:5])
    print("==================================")
    l = l.sort_values(by="BIC",ascending=True)
    print(l[:5])
 
 
 
def predict_VAR(df, nb_order=10, prop_train=0.9):
        train = df[:round(prop_train*len(df))]
        test = df[round((prop_train)*len(df)):]
        model = VAR(df)
        results = model.fit(nb_order)
        print(results.summary())
        lag_order = results.k_ar
        forecast = results.forecast(y=train.values[-lag_order:],steps = round((1-prop_train)*len(df)))
        df_forecast = pd.DataFrame(forecast, index=test.index, columns=['SPY_forecast', 'IR10_forecast', 'IR3_forecast'])
        alls = df.merge(df_forecast, how='left', left_index=True, right_index=True)
       
        plt.figure(1)
        alls[["SPY","SPY_forecast"]].plot()
        plt.figure(1)
        alls[["IR10", "IR10_forecast"]].plot()
        plt.figure(3)
        alls[["IR3", "IR3_forecast"]].plot()
       
        #model evaluation :
        RMSE1 = rmse(alls.dropna()["SPY"], alls.dropna()["SPY_forecast"] )
        print("RMSE for SPY is :", RMSE1)
        RMSE2 = rmse(alls.dropna()["IR10"], alls.dropna()["IR10_forecast"] )
        print("RMSE for IR10 is :", RMSE2)
        RMSE3 = rmse(alls.dropna()["IR3"], alls.dropna()["IR3_forecast"] )
        print("RMSE for IR3 is :", RMSE3)
 
 
def grangers_causation_matrix(data, test='ssr_chi2test', verbose=False, maxlag=30):   
    
    """
    H0 : the coefficient of the corresponding causation is 0; That is the X does not cause Y.
    If P-values in the table < significance level (0.05) --> H0 can be rejected .
    """
    variables = data.columns
    df = pd.DataFrame(np.zeros((len(variables), len(variables))), columns=variables, index=variables)
    for c in df.columns:
        for r in df.index:
            test_result = grangercausalitytests(data[[r, c]], maxlag=maxlag, verbose=False)
            p_values = [round(test_result[i+1][0][test][1],4) for i in range(maxlag)]
            if verbose: print(f'Y = {r}, X = {c}, P Values = {p_values}')
            min_p_value = np.min(p_values)
            df.loc[r, c] = min_p_value
    df.columns = [var + '_x' for var in variables]
    df.index = [var + '_y' for var in variables]
    return df 
 


#%%
 
 
# SP500 (daily close)
series_SPY =  np.log(get_FRED_series('SP500')["2013":"2023"])
 
#Market Yield on U.S. Treasury Securities at 10-Year Constant Maturity, Quoted on an Investment Basis (DGS10)
series_ir10 = np.log(get_FRED_series('DGS10')["2013":"2023"])
 
#3-Month Treasury Bill Secondary Market Rate
series_ir3 = np.log(get_FRED_series('DTB3')["2013":"2023"].resample('D').interpolate(method='linear'))

series3 = np.log(get_FRED_series('GEPUCURRENT')["2013":"2023"].resample('D').interpolate(method='linear'))

series4 = np.log(get_FRED_series('T5YIE')["2013":"2023"].resample('D').interpolate(method='linear'))

series5 = np.log(get_FRED_series('OVXCLS')["2013":"2023"].resample('D').interpolate(method='linear'))

series6 = np.log(get_FRED_series('MEDCPIM158SFRBCLE')["2013":"2023"].resample('D').interpolate(method='linear'))

series7 = np.log(get_yahoo_series('META')["2013":"2023"])

series8 = np.log(get_yahoo_series('GOOG')["2013":"2023"])

series9 = np.log(get_yahoo_series('INFL.SW')["2013":"2023"])

series10 = np.log(get_yahoo_series('TSLA')["2013":"2023"])




merged = pd.concat([series7, series8, series10], axis=1)
merged = merged["2020":]
merged = merged.replace([np.inf, -np.inf], np.nan)
merged = merged.apply(lambda x: x.fillna(x.mean()), axis=1)
merged.columns = ['1', '2', "4"]
merged.plot()

johansen_coint_test(merged)


for i in range(len(merged.columns)):
    print(ADF_reg_detail(merged.iloc[:,i]))


johansen_coint_test(merged)

grangers_causation_matrix(merged)





