# This file is used to solve the differential equation and fit the data

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit, fsolve

from Static import *

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
    initial_guess = [-1] * order
    
    # Use a lambda function to adapt solve_diff_eq to the requirements of curve_fit
    popt, pcov = curve_fit(lambda t, *a: solve_diff_eq(t, a, S0, Sf), t, S, p0=initial_guess, maxfev=100000)
    
    # 打印拟合得到的系数
    print(f'Coefficients: {popt}')
    
    # 计算参数的标准误差
    perr = np.sqrt(np.diag(pcov))
    
    # 打印参数的标准误差
    print(f'Standard errors: {perr}')
    
    #find the roots of fitted equation
    roots = find_roots(popt, Sf, S0)
    print('roots:' + str(roots))
    
    # 使用拟合得到的系数预测y值
    ypred = solve_diff_eq(t, popt, S0, Sf)
    
    return ypred

def find_roots(a, Sf, S_guess, search_range=(0, 10), guess_density=100):
    unique_roots = []
    tolerance = 1e-5  # 用于检查根是否唯一的容差
    guesses = np.linspace(search_range[0], search_range[1], guess_density)

    for guess in guesses:
        try:
            root = fsolve(lambda S: generate_diff_eq(0, S, a, Sf), guess)[0]
            if not any(np.isclose(root, uroot, atol=tolerance) for uroot in unique_roots):
                unique_roots.append(root)
        except (ValueError, RuntimeError):
            continue  # 如果找不到根或者遇到错误，就跳过

    return unique_roots

def main():
    # 读取Excel文件
    df0 = pd.read_excel('RWdata300.xlsx')

    # 提取entropy_ave列和time列
    y0 = np.array(df0['entropy_ave'].tolist())
    yerr0 = np.array(df0['entropy_std'].tolist())
    x = df0['time'].tolist()

    docurvefit(y0, T0)
    
    print(f'yerr0 errors: {np.average(yerr0)}')
    
    
    plt.fill_between(x, y0 - yerr0, y0 + yerr0, color=(1, 0, 0, 0.5))  # 红色误差阴影
    #plt.fill_between(x, y1 - yerr1, y1 + yerr1, color=(0, 1, 0, 0.5))  # 绿色误差阴影
    # plt.fill_between(x[::100], y2[::100] - yerr2[::100], y2[::100] + yerr2[::100], color=(0, 0, 1, 0.5))  # 蓝色误差阴影
    plt.plot(x, y0, 'b-', label='T0')
    #plt.plot(x, y1, 'y-', label='T1')
    # plt.plot(time_values, y2, 'g-', label='T2')
    plt.legend()
    plt.title('Error Bars, ' + str(order) + 'th order')
    plt.show()
        
if __name__ == '__main__':
    main()