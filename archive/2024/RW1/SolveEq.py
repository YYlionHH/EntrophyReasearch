# This file is used to solve the differential equation and fit the data

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit

# # Define the differential equation
# def threecoe(t, S, a1, a2, a3, Sf):
#     return a1 * (S - Sf) + a2 * (S - Sf)**2 + a3 * (S - Sf)**3

# Generate a differential equation model of any order
def generate_diff_eq(t, S, a, Sf):
    # Initialize the result of the equation to 0
    equation_result = 0

    # Loop through each coefficient in the list 'a'
    for i in range(len(a)):
        # Calculate the term for the current coefficient and add it to the result
        term = a[i] * (S - Sf) ** (i + 1)
        equation_result += term
        
    return equation_result

# Function to solve the differential equation, now can handle a coefficient list of any length
def solve_diff_eq(t, a, S0, Sf):
    # Use a lambda function to adapt generate_diff_eq to the requirements of solve_ivp
    sol = solve_ivp(lambda t, S: generate_diff_eq(t, S, a, Sf), [t[0], t[-1]], [S0], t_eval=t)
    return sol.y[0]

# Fitting function, adapts to a coefficient list of any length
def fit_data(t, S, S0, Sf, order):
    # Initialize the guess for coefficients
    initial_guess = [1] * order
    
    # Use a lambda function to adapt solve_diff_eq to the requirements of curve_fit
    popt, pcov = curve_fit(lambda t, *a: solve_diff_eq(t, a, S0, Sf), t, S, p0=initial_guess)
    
    # 打印拟合得到的系数
    print(f'Coefficients: {popt}')
    
    # 计算参数的标准误差
    perr = np.sqrt(np.diag(pcov))
    
    # 打印参数的标准误差
    print(f'Standard errors: {perr}')
    
    # 使用拟合得到的系数预测y值
    ypred = solve_diff_eq(t, popt, S0, Sf)
    
    return ypred


#def main():
    
        
# if __name__ == '__main__':
#     main()