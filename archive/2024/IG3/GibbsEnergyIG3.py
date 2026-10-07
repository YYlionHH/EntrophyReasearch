# This program use the Gibbs free energy equation to calculate the change in Gibbs free energy of a reaction.
# this program assume particles expand to larger box at time t0
# there is a barrier in the middle of the box, with a small square hole in it
# all particles will be be on one side of the wall, and the opening will allow them to slowly migrate to both sides
# 

# Import modules
import numpy as np 
from scipy.optimize import curve_fit
from scipy import special, optimize 

import matplotlib.pyplot as plt
from decimal import Decimal
import math as math
import pandas as pd
import csv
from tqdm import tqdm

from SolveEq import *

# Define constants
L = 1.0 # length of the box in meters
T = 300 # Temperature of the system in Kelvin
k = 1.38e-23 # Boltzmann constant in J/K
m = 1.67e-27 # Mass of a hyrdrogen particle in kg

N = 50000 # Number of particles in the system
tmax = 1 # Maximum time in ss
t0 = 0 #tmax * 0.2  # the time when the box disappear, it should less than tmax
t1 = tmax * 0.6
dt = 1e-4 # Time step in s
# the length of the square hole
L_S = L/20

# larger box length, less wave range
# Initial length of the tiny box in meters
L0 = L/20

# the order of formular when do the curvefit
order = 3

# create a time list to store the time values
time_values = np.arange(0, tmax, dt)
    
def initialV(Temperature):
     
    # Calculate the standard deviation for the Maxwell-Boltzmann distribution
    sigma = np.sqrt(k * Temperature / m)

    # Generate random speeds
    vx = np.random.normal(loc=0, scale=sigma, size=N)
    vy = np.random.normal(loc=0, scale=sigma, size=N)
    vz = np.random.normal(loc=0, scale=sigma, size=N)

    # find maximum velocity
    max_vx = np.max(np.abs(vx))
    max_vy = np.max(np.abs(vy))
    max_vz = np.max(np.abs(vz))
    
    return vx,vy,vz,max_vx,max_vy,max_vz

def kinematics(x, y, z, vx, vy, vz, t):
    # Returns the updated position and velocity after one time step
    # x, y, z are the positions in m, vx, vy, vz are the velocities in m/s, dt is the time step in s
    # at time t0, the tiny box disappear, and the particles are in the larger box
    
    if t < t0: # If the time is before t0, the particles are in the larger box
        # LLower = L/2 - L0/2
        # LUpper = L/2 + L0/2
        LLower = 0
        LUpper = L0
        
    # elif t >= t0 and t < t1: # If the time is between t0 and t1, the particles are in the tiny box
    #     LLower = L/2 - L0
    #     LUpper = L/2 + L0
    else:
        LLower = 0 # set a varible to classify whether the box disappeared, this is lower bound
        LUpper = L # set a varible to classify whether the box disappeared, this is upper bound
        
    x += vx * dt            # Update x position by kinematics
    if x <= LLower:         # Check if the particle hits the left wall
        x = 2*LLower - x    # Make the position positive
        vx = -vx            # Change the sign of the x velocity
    elif x >= LUpper:       # Check if the particle hits the right wall
        x = 2*LUpper - x    # Reflect the position from the wall
        vx = -vx            # Change the sign of the velocity

    y += vy * dt            # Update y position by kinematics
    if y <= LLower:         # Check if the particle hits the bottom wall
        y = 2*LLower - y    # Make the position positive
        vy = -vy            # Change the sign of the velocity
    elif y >= LUpper:       # Check if the particle hits the top wall
        y = 2 * LUpper - y  # Reflect the position from the wall
        vy = -vy            # Change the sign of the velocity
        

    z += vz * dt            # Update z position by kinematics
    if z <= LLower:         # Check if the particle hits the back wall
        z = 2*LLower - z    # Make the position positive
        vz = -vz            # Change the sign of the velocity
    elif z >= LUpper:       # Check if the particle hits the front wall
        z = 2 * LUpper - z  # Reflect the position from the wall
        vz = -vz            # Change the sign of the velocity
        
    # if the particle hit the middle of the box
    if (x - vx * dt >= L/2 and x < L/2) or (x - vx * dt <= L/2 and x > L/2) :
        
        def line(x0 ,veo ,pos):
            a = veo/vx
            b = pos - a * x
            return a*x0 + b
        
        if line(L/2 , vy ,y) > L/2 + L_S or line(L/2 , vy ,y) < L/2 - L_S or line(L/2 , vz ,z) > L/2 + L_S or line(L/2 , vz ,z) < L/2 - L_S:
            x = L - x
            vx = -vx
        
    return x, y, z, vx, vy, vz


# Define a function to calculate the entropy by Gibbs formula
def Gibbs(n):
    # Returns the Gibbs entropy of a system with N particles and n cells
    # k is the Boltzmann constant in J/K
    
    # use a if statement to avoid the error of log(0)
    if n == 0 or n == N:
        return 0
    
    else:
        # Calculate the probability of a particle being in a cell
        Pro = n/N

        # Calculate the entropy by Gibbs formula
        # the reason of times 1e22 is easy to do the curvefiting
        S = -k * Pro * np.log(Pro) * 1e22

    return S

# Define a function to calculate the entropy by Gibbs formula
def entropy(x):
    # Initialize the sum to zero
    Left, Right = 0,0

    for j in range(N):
        # Find the cell index for the particle
        
        # chech if the postion of all dimension is in the range of the box
        # or on the boundary
        if x[j] > 0 and x[j] < L/2:
            Left += 1
        else:
            Right += 1
            
    # Calculate the entropy by Gibbs formula for each cell
    S = Gibbs(Left) + Gibbs(Right)
             
    return S

    
def Kinematics(T,vx, vy, vz, max_vx, max_vy, max_vz):
    
    # Define a list to store the entropy values at each time step
    entropy_values = []
    
    # Initialize the position arrays
    #generate random posion for x,y,z, and velocity for x,y,z
    x = np.random.uniform(0, L0, size=N)
    y = np.random.uniform(0, L , size=N)
    z = np.random.uniform(0, L , size=N)
    
    
    print(f"max_vx: {max_vx}, max_vy: {max_vy}, max_vz: {max_vz} with T = {T} K")
    
    # Calculate total iterations
    total_iterations = int(len(time_values))

    # Create a progress bar
    progress_bar = tqdm(total=total_iterations, desc="Processing", bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt}')

    # Kinematics loop
    t = 0 # Initial time in s
    for j in range(total_iterations): # Loop until maximum time is reached
        
        # Update the position and velocity of each particle by kinematics
        for i in range(N): 
            x[i], y[i], z[i], vx[i], vy[i], vz[i] = kinematics(x[i], y[i], z[i], vx[i], vy[i], vz[i], t)
            
        # Call the s_script function and append the result to the s_script_values list
        entropy_values.append(entropy(x))

        t += dt # Update time by time step
        
        # Update the progress bar
        progress_bar.update(1)

    # Close the progress bar
    progress_bar.close()
    return entropy_values


def main():
    vx, vy, vz, max_vx, max_vy, max_vz = initialV(T)
    entropy_values = Kinematics(T,vx, vy, vz, max_vx, max_vy, max_vz)

    #print(entropy_values[:len(entropy_values)//100])
    print(len(entropy_values))
    print(len(time_values))
    
    # 创建一个字典，其中键是列名，值是数据
    data = {"time": time_values, "entropy": entropy_values}

    # 将字典转换为 DataFrame
    df = pd.DataFrame(data)

    # 将 DataFrame 保存到 Excel 文件中
    df.to_excel('scripts2.xlsx', index=False)

    # do curvefit with entropy_values and time_values
    # define the data 
    # 步骤 1: 找到大于 t0 的第一个值的索引
    index = next(i for i, v in enumerate(time_values) if v > t0)

    # 步骤 2: 生成 time_values1
    time_values1 = time_values[index:]
    entropy_values1 = entropy_values[index:]
    
    # use the data to estimate the parameters
    S0 = entropy_values1[1]
    Sf = entropy_values1[-1]
    
    # fit_data(t, S, S0, Sf, order)
    y_pred = fit_data(time_values1, entropy_values1, S0, Sf, order)
    
    # 绘制拟合后的模型曲线和原始数据
    plt.xlabel('Time (s)')
    plt.ylabel('S: Entropy(J/K)')
    plt.plot(time_values1, y_pred, 'r-', label='fit')
    plt.plot(time_values, entropy_values, 'b-', label='data')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()