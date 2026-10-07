#This program use the Gibbs free energy equation to calculate the change in Gibbs free energy of a reaction.
# this program assume particles expand to larger box at time t0
# the difference of GibbsEnergy is This one with three coefficent
# this program will simulate the random walk of particles in a box and calculate the entropy of the system
# the main idea of random walk is that every partle wiil move same distance but random direction
# generate the veolicity in sphere coordinate, and then convert it to the x,y,z coordinate

# Import modules
import numpy as np
import matplotlib.pyplot as plt 
from decimal import Decimal 
import math as math 
import pandas as pd
import csv
from tqdm import tqdm

from SolveEq import *

# Define constants
L = 1.0 # length of the box in meters
N = 1000 # Number of particles in the system
T = 300 # Temperature of the system in Kelvin
k = 1.38e-23 # Boltzmann constant in J/K
m = 1.67e-27 # Mass of a hyrdrogen particle in kg
dt = 1e-1 # Time step in s defalut: 1e-5
tmax = 1e2 # Maximum time in ss defalut: 1e-2
t0 = 0#tmax * 0.4  # the time when the box disappear, it should less than tmax
t1 = tmax * 0.6

# larger box length, less wave range
# Initial length of the tiny box in meters, it should be odd number 
# defalut: 1/10 of the box length
L0 = L/10 

# step length of random walk
# defalut: 1/50
l = L/200

#This part refer to SolveEq.py
order = 3

# # Initialize the position arrays
# #generate random posion for x,y,z, and velocity for x,y,z
# # x = np.random.uniform(0, L0, size=N)
# # y = np.random.uniform(0, L0, size=N)
# # z = np.random.uniform(0, L0, size=N)
# x = np.random.uniform(L/2-L0/2, L/2+L0/2, size=N)
# y = np.random.uniform(L/2-L0/2, L/2+L0/2, size=N)
# z = np.random.uniform(L/2-L0/2, L/2+L0/2, size=N)   


def initialV(Temperature):
     
    # # Calculate the standard deviation for the Maxwell-Boltzmann distribution
    # sigma = np.sqrt(k * Temperature / m)

    # # Generate random speeds
    # vtotal = np.random.normal(loc=0, scale=sigma, size=N)
    
    vtotal = l/dt
    
    # find maximum velocity
    max_vtotal = np.max(np.abs(vtotal))
    
    return vtotal,max_vtotal

# Define cells for the entropy calculation
# number of cell in x axis
# one side
# cell_x  = int(L/(1.3001*L0))

# define the number of cells in x axis
cell_x  = int(L/(0.5*L0))

# create a time list to store the time values
time_values = np.arange(0, tmax, dt)

def kinematics(x, y, z, vtotal, dt, t):
    # Returns the updated position and velocity after one time step
    # x, y, z are the positions in m, vx, vy, vz are the velocities in m/s, dt is the time step in s
    # at time t0, the tiny box disappear, and the particles are in the larger box
    
    if t < t0: # If the time is before t0, the particles are in the larger box
        LLower = L/2 - L0/2
        LUpper = L/2 + L0/2
        # LLower = 0
        # LUpper = L/3
        
    # elif t >= t0 and t < t1: # If the time is between t0 and t1, the particles are in the tiny box
    #     LLower = L/2 - L0
    #     LUpper = L/2 + L0
    else:
        LLower = 0 # set a varible to classify whether the box disappeared, this is lower bound
        LUpper = L # set a varible to classify whether the box disappeared, this is upper bound
        
    # generate random direction of the velocity
    
    # generate random direction of the velocity
    theta = 2 * np.pi * np.random.rand()
    phi = np.arccos(2 * np.random.rand() - 1)

    # calculate the velocity in x, y, z direction
    vx = vtotal * np.sin(phi) * np.cos(theta)
    vy = vtotal * np.sin(phi) * np.sin(theta)
    vz = vtotal * np.cos(phi)
    
    x += vx * dt            # Update x position by kinematics
    if x <= LLower:         # Check if the particle hits the left wall
        x = 2*LLower - x    # Make the position positive
    elif x >= LUpper:       # Check if the particle hits the right wall
        x = 2*LUpper - x  # Reflect the position from the wall

    y += vy * dt            # Update y position by kinematics
    if y <= LLower:         # Check if the particle hits the bottom wall
        y = 2*LLower - y    # Make the position positive
    elif y >= LUpper:       # Check if the particle hits the top wall
        y = 2 * LUpper - y  # Reflect the position from the wall
        
    z += vz * dt            # Update z position by kinematics
    if z <= LLower:         # Check if the particle hits the back wall
        z = 2*LLower - z    # Make the position positive
    elif z >= LUpper:       # Check if the particle hits the front wall
        z = 2 * LUpper - z  # Reflect the position from the wall
    
    return x, y, z


    
    
# Define a function to calculate the entropy by Gibbs formula
def Gibbs(n):
    # Returns the Gibbs entropy of a system with N particles and n cells
    # k is the Boltzmann constant in J/K
    
    # use a if statement to avoid the error of log(0)
    if (n == 0).all():
        return 0
    
    else:
        # Find the indices of the zeros in n
        zero_indices = n == 0

        # Calculate the probability of a particle being in a cell
        Pro = n/N
        Pro[zero_indices] = 1e-10  # replace zeros with a small number

        # Calculate the entropy by Gibbs formula
        # the reason of times 1e22 is easy to do the curvefiting
        S = -k * Pro * np.log(Pro) *1e22

        # Set the entropy of the zero cells to 0
        S[zero_indices] = 0
    return S

# Define a function to calculate the entropy by Gibbs formula
def entropy(x,y,z,cell_x):
    # Initialize the sum to zero
    S = 0
    
    # Initialize the cells to zero
    cells = np.zeros((cell_x,cell_x,cell_x))
    
    for j in range(N):
        # confirm the index of the cell
        x_index = int(x[j]//L0)
        y_index = int(y[j]//L0)
        z_index = int(z[j]//L0)
        
        # check if the index is out of the range
        if x_index >= cell_x:
            x_index = cell_x - 1
        if y_index >= cell_x:
            y_index = cell_x - 1
        if z_index >= cell_x:
            z_index = cell_x - 1# Find the cell index for the particle
        cells[x_index, y_index, z_index] += 1
            
    # Calculate the entropy by Gibbs formula for each cell
    S += np.sum(Gibbs(cells))
             
    return S

    
def Kinematics(T):
    
    # Define a list to store the entropy values at each time step
    entropy_values = []
    
    # Initialize the position arrays
    #generate random posion for x,y,z, and velocity for x,y,z
    # x = np.random.uniform(0, L0, size=N)
    # y = np.random.uniform(0, L0, size=N)
    # z = np.random.uniform(0, L0, size=N)
    x = np.random.uniform(L/2-L0/2, L/2+L0/2, size=N)
    y = np.random.uniform(L/2-L0/2, L/2+L0/2, size=N)
    z = np.random.uniform(L/2-L0/2, L/2+L0/2, size=N)
    
    vtotal, vtotal_max = initialV(T)
    print(f"vtotal_max: {vtotal_max}, with T = {T} K") 

    # Calculate total iterations
    total_iterations = len(time_values)

    # Create a progress bar
    progress_bar = tqdm(total=total_iterations, desc="Processing", bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt}')

    # Kinematics loop
    t = 0 # Initial time in s
    for j in range(total_iterations): # Loop until maximum time is reached
        
        # Update the position and velocity of each particle by kinematics
        for i in range(N): 
            x[i], y[i], z[i] = kinematics(x[i], y[i], z[i], vtotal, dt, t)
            
        # Call the s_script function and append the result to the s_script_values list
        entropy_values.append(entropy(x,y,z,cell_x))

        t += dt # Update time by time step
        
        # Update the progress bar
        progress_bar.update(1)

    # Close the progress bar
    progress_bar.close()
    return entropy_values



def main():
    
    # Check if the time step is optimally small
    # Chech if the change of position larger than tiny box
    if L0  < l :
        print("The parameter look like not reasonable, let l < L0")
    
    entropy_values = Kinematics(T)

    #print(entropy_values[:len(entropy_values)//100])
    #print(entropy_values)
    
    # 创建一个字典，其中键是列名，值是数据
    data = {"time": time_values, "entropy": entropy_values}
    #data = {"time": time_values[::10], "entropy": entropy_values[::10]}
    print(f'time_values {len(time_values)}')
    print(f'entropy_values {len(entropy_values)}')
    
    # 将字典转换为 DataFrame
    df = pd.DataFrame(data)

    # 将 DataFrame 保存到 Excel 文件中
    df.to_excel('randomwalk.xlsx', index=False)

    # do curvefit with entropy_values and time_values
    # define the data 
    entropy_values1 = entropy_values[int(t0//dt):]
    time_values1 = np.linspace(t0 , tmax, len(entropy_values1))

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