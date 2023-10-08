import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation

data_d  = pd.read_csv("data/densidad.dat", delimiter='\t', skip_blank_lines=False)

data_d.columns = ["t", "x", "u"]
# print(data_d.head())
# print(data_2.head())
# print(data_d.shape[0])
print(data_d['u'].max())

vacios_d = data_d[pd.isna(data_d['t'])]
bounds_d = np.array([-1] + list(vacios_d.index))
# print(bounds_d)
desfases_d = np.zeros(len(bounds_d), dtype=int)
desfases_d[::2] = 1
bounds_d = np.add(bounds_d, desfases_d)
bounds_d = bounds_d[:-1]
# for i in range(len(bounds_d)//2)[::2]:
#     print(bounds_d[i], bounds_d[i+1])

subsets_d = [data_d.iloc[bounds_d[i]:bounds_d[i+1], :] for i in range(len(bounds_d))[::2]]

# Create a figure and axis for the animation
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

# Create a line object for the plot
line, = ax1.plot([], [], '-', markersize=5)

def init():
    line.set_data([], [])
    return (line,)

# Define the animation function
def animate(i):
    # print(f"Animating frame {i}")
    # Get the x and y data_d for the current timestep
    x = subsets_d[i]['x'].values
    u = subsets_d[i]['u'].values
    
    # Update the data_d for the line
    line.set_data(x, u)
    
    # Set the x and y limits of the axis
    min_u = data_d['u'].min()
    max_u = data_d['u'].max()
    min_x = 0
    max_x = 100 
    ax1.set_xlim([min_x, max_x])
    ax1.set_ylim([min_u, max_u])
    # Get time
    t = subsets_d[i]['t'].values[0]
    # Set the title of the axis to the current timestamp
    ax1.set_title("Densidad, t={} s".format(t))
    return line,
    

num_frames = len(subsets_d)
dt = (subsets_d[1]['t'].values[0] - subsets_d[0]['t'].values[0])*1000

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=(range(num_frames)), interval= 100)

###############################

data  = pd.read_csv("data/presion.dat", delimiter='\t', skip_blank_lines=False)

data.columns = ["t", "x", "u"]
# print(data.head())
# print(data_2.head())
# print(data.shape[0])

vacios = data[pd.isna(data["t"])]
bounds = np.array([-1] + list(vacios.index))
# print(bounds)
desfases = np.zeros(len(bounds), dtype=int)
desfases[::2] = 1
bounds = np.add(bounds, desfases)
bounds = bounds[:-1]
# for i in range(len(bounds)//2)[::2]:
#     print(bounds[i], bounds[i+1])

subsets = [data.iloc[bounds[i]:bounds[i+1], :] for i in range(len(bounds))[::2]]



# Create a line object for the plot
line2, = ax2.plot([], [], '-', markersize=5)

def init2():
    line2.set_data([], [])
    return (line2,)


# Define the animation function
def animate_2(i):
    # print(f"Animating frame {i}")
    # Get the x and y data for the current timestep
    x = subsets[i]['x'].values
    u = subsets[i]['u'].values
    
    # Update the data for the line
    line2.set_data(x, u)
    
    # Set the x and y limits of the axis
    min_u = data['u'].max()
    max_u = data['u'].min()
    min_x = 0
    max_x = 100 
    ax2.set_xlim([min_x, max_x])
    ax2.set_ylim([min_u, max_u])
    # Get time
    t = subsets[i]['t'].values[0]
    # Set the title of the axis to the current timestamp
    ax2.set_title("Presión, t={} s".format(t))
    return line2,
    

num_frames = len(subsets)
dt = (subsets[1]['t'].values[0] - subsets[0]['t'].values[0])*1000

anim2 = animation.FuncAnimation(fig, animate_2, init_func=init2, frames=(range(num_frames)), interval= 100)

plt.show()