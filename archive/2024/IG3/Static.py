# this program will do some static analysis on the Gibbs'energy
# this program needs the GibbsEnergy111.py
# corresponding random walk simulation

from GibbsEnergyIG3 import *

#define some new constants
NumberofSamples = 50

# the temperature will change to T0
T0 = 300
T1 = 400
T2 = 600

# define a function to do the static analysis
# do the simlutate for 50 times in GibbsEnergy11.py, the calculate the average value of every points
# create the error and standard deviation of the Gibbs'energy
def StaticAnalysis(tempreture):
    # to store the every results by Kinematics()
    results = []
    
    # do the Kinematics() many times
    for i in range(NumberofSamples):
        
        print(i + 1 ," of ", NumberofSamples )
        vx,vy,vz,max_vx,max_vy,max_vz = initialV(tempreture)
        result = Kinematics(tempreture,vx,vy,vz,max_vx,max_vy,max_vz)
        results.append(result)
        if i == 0:
            data = {"time": time_values, "entropy": result}
        else:
            data["entropy_" + str(i)] = result
            
    # calculate the average value of every points
    AveS = np.mean(results, axis=0)
    StDS = np.std(results, axis=0)
    
    data["entropy_ave"] = AveS
    data["entropy_std"] = StDS
    
    # 将字典转换为 DataFrame
    df = pd.DataFrame(data)

    # 将 DataFrame 保存到 Excel 文件中
    dataname = 'RWdata' + str(tempreture) + '.xlsx'
    df.to_excel(dataname, index=False)
    return AveS, StDS

# this function is used to plot the static analysis
def docurvefit(AveS, x,tempreture):
    
    # due to the curve_fit function in python must use two array with same length
    # then we have to use the average value to do the curve fit

    # plot all the results
    plt.scatter(x, AveS, s=2, zorder=2)
    
    S0 = AveS[0]
    Sfestimate  = AveS[-1]

    # calculate the predicted values with the optimal parameters
    y_pred = fit_data(x, AveS, S0, Sfestimate, order)
    
    plt.title('Curfit of Gibbs Energy with particles ' + str(N))

    # plot the fitted curve
    plt.plot(x, y_pred, 'r-', zorder=3)
    
    plt.show()
    
    
    
# this function is used to plot the uncertainty of the change of entropy
def plotuncertainty(dS, yerr0, yerr1, tempreture0, tempreture1):
    
    # calculate the uncertainty of the change of entropy
    dSerror = np.sqrt(yerr0**2 + yerr1**2)
    
    #   do curve fit of flat line
    # define the model function
    def flatline(time, a, b):
        return a * time + b

    # use curve_fit to find the optimal parameters
    popt, pcov = curve_fit(flatline, time_values, dSerror)
    
    # calculate the standard error of the parameters
    perr = np.sqrt(np.diag(pcov))

    # print the standard error of the parameters
    print(f'Standard errors: {perr}')

    # print the optimal parameters
    print(f'Optimal parameters: {popt}')
    
    # plot the curve fit
    plt.plot(time_values, flatline(time_values, *popt), 'r-', label='Curve Fit')
    
    # plot the uncertainty
    plt.plot(time_values, dS, 'b-', label='T')
    plt.fill_between(time_values[::100], dS[::100] - dSerror[::100], dS[::100] + dSerror[::100], color=(1, 0, 0, 0.5))
    plt.title('Uncertainty of change of entropy from temperature' + str(tempreture0) + ' to ' + str(tempreture1))
    plt.xlabel('Time (s)')
    plt.ylabel('delta S (J/K)')
    plt.legend()
    plt.show()
        

def main():
    T = T0
    
    # draw the error bar plot
    x = time_values
    y0, yerr0 = StaticAnalysis(T)
    
    print(f'yerr0 errors: {np.average(yerr0)}')
    
    # do and draw curvefit 
    docurvefit(y0,time_values, T)
    
    plt.xlabel('Time (s)')
    plt.ylabel('S: Entropy(J/K)')
    plt.plot(time_values, y0, 'b-', label='data')
    plt.fill_between(x, y0 - yerr0, y0 + yerr0, color=(1, 0, 0, 0.5))  # 红色误差阴影，覆盖整个数据范围
    plt.legend()
    plt.title('order ' + str(order) + ' Random work Error Bars' + 'N= ' + str(N) + ' in' + str(NumberofSamples)  + ' times')
    plt.show()
    
if __name__ == "__main__":
    main()