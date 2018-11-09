import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

#Tracklength sensitivity plot
fig, ax = plt.subplots()
data = np.loadtxt('./sensitivity/cutoff_length_sensitivity_data')
plt.scatter(data[:,0],data[:,1])
plt.title('Variation due to Tracklength, 100 Rays, No Deadzone')
plt.ylabel('$k$')
plt.xlabel('Tracklength (cm)')
plt.show()

#Deadzone sensitivity plot
fig, ax = plt.subplots()
data = np.loadtxt('./sensitivity/dead_zone_sensitivity_data')
plt.scatter(data[:,0],data[:,1])
plt.title('Variation due to Deadzone, 100 Rays, Tracklength = 300 cm')
plt.ylabel('$k$')
plt.xlabel('Deadzone (cm)')
plt.show()

#nrays sensitivity plot
fig, ax = plt.subplots()
data = np.loadtxt('./sensitivity/nrays_sensitivity_data_2')
plt.scatter(data[:,0],data[:,1])
plt.title('Variation due to Number of Rays, 100 Rays, Tracklength = 300 cm, Deadzone = 50 cm')
plt.ylabel('$k$')
plt.xlabel('Number of rays')
plt.show()
