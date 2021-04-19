import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

std = pd.read_csv('data/new/2.csv', index_col = None, header = None)[0:5000];
std[0]=(1.0*(std[0] - std[0][0]) /1000000000).astype('float');

time = std[0].values

## Acceleration

std_acc = std[std[1]==' ACC_UN']

acc_x = std_acc[2].values
acc_y = std_acc[3].values
acc_z = std_acc[4].values

acc_t = std_acc[0].values

plt.figure(1)
plt.plot(acc_t, acc_x)
plt.plot(acc_t, acc_y)
plt.plot(acc_t, acc_z)
plt.legend(['acc x','acc y','acc z'], loc = 'upper right')
plt.xlabel('time (s)')
plt.ylabel('acceleration (m/s^2)')

## Gyroscope

std_gyro = std[std[1]==' GYRO_UN']

gyro_x = std_gyro[2].values
gyro_y = std_gyro[3].values
gyro_z = std_gyro[4].values

gyro_t = std_gyro[0].values

plt.plot(gyro_t, gyro_x)
plt.plot(gyro_t, gyro_y)
plt.plot(gyro_t, gyro_z)
plt.legend(['gyro x','gyro y','gyro z'], loc = 'upper right')
plt.xlabel('time (s)')
plt.ylabel('angular velocity (rad/s)')


## Stability Detection

g = np.mean(acc_z)

sig_acc = g
sig_gyro = 0 + 0.3

w_min = 0.3
w_max = 1.1

std_acc  = np.zeros(len(acc_t))
std_gyro = np.zeros(len(acc_t))

std_static  = np.zeros(len(acc_t))

i = 0
indexR = 66

while(acc_t[i] + w_max < acc_t[-1]):

    std_part = std[std[0] < (acc_t[i] + w_min)]
    std_part = std_part[std_part[0] > acc_t[i]]
    
    # get acceleration in the window range
    
    acc_part = std_part[std_part[1]==' ACC_UN']
    acc_part_x = acc_part[2].values
    acc_part_y = acc_part[3].values
    acc_part_z = acc_part[4].values
    
    
    std_acc_x = np.std(acc_part_x) # std of accx
    std_acc_y = np.std(acc_part_y) # std of accy
    std_acc_z = np.std(acc_part_z) # std of accz

    std_acc[i] = np.sqrt(std_acc_x**2 + std_acc_y**2 + std_acc_z**2)
    
    # get gyroscope in the window range
    
    gyro_part = std_part[std_part[1]==' GYRO_UN']
    gyro_part_x = gyro_part[2].values
    gyro_part_y = gyro_part[3].values
    gyro_part_z = gyro_part[4].values
    
    std_gyro_x = np.std(gyro_part_x) # std of gyrox
    std_gyro_y = np.std(gyro_part_y) # std of gyrox
    std_gyro_z = np.std(gyro_part_z) # std of gyrox

    std_gyro[i] = np.sqrt(std_gyro_x**2 + std_gyro_y**2 + std_gyro_z**2)
    
    
    if std_acc[i] < 2*sig_acc and std_gyro[i] < 3*sig_gyro:
        
        # increase to max window size
        
        std_part = std[std[0] < (time[i] + w_max)]
        std_part = std_part[std_part[0] > time[i]]

        # accleration in the max window rnage
        acc_part = std_part[std_part[1]==' ACC_UN']
        acc_part_x = acc_part[2].values
        acc_part_y = acc_part[3].values
        acc_part_z = acc_part[4].values
        
        std_acc_x = np.std(acc_part_x) # std of accx
        std_acc_y = np.std(acc_part_y) # std of accy
        std_acc_z = np.std(acc_part_z) # std of accz

        std_acc_ = np.sqrt(std_acc_x**2 + std_acc_y**2 + std_acc_z**2)
        
        # gyroscope in the max window range
        
        gyro_part = std_part[std_part[1]==' GYRO_UN']
        gyro_part_x = gyro_part[2].values
        gyro_part_y = gyro_part[3].values
        gyro_part_z = gyro_part[4].values
        
        std_gyro_x = np.std(gyro_part_x) # std of gyrox
        std_gyro_y = np.std(gyro_part_y) # std of gyrox
        std_gyro_z = np.std(gyro_part_z) # std of gyrox

        std_gyro_ = np.sqrt(std_gyro_x**2 + std_gyro_y**2 + std_gyro_z**2)
    
        
        
        ## update parameters
        if std_acc_ < 2*sig_acc and std_gyro_ < 3*sig_gyro:
            
            if std_acc_ < sig_acc:
                sig_acc = std_acc_
                
            if std_gyro_ < sig_gyro:
                sig_gyro = std_gyro_
                
            if acc_t[i] > indexR:
                std_static[i] = 1
        
    
    ##############

    i = i + 1

plt.figure(2)
plt.plot(acc_t,std_acc)

plt.figure(3)
plt.plot(acc_t,std_gyro)

plt.figure(4)
plt.plot(acc_t, acc_x)
plt.plot(acc_t, acc_y)
plt.plot(acc_t, acc_z)
plt.plot(acc_t,std_static*np.max(acc_z))
plt.legend(['acc x','acc y','acc z','detector'], loc = 'upper right')
plt.xlabel('time (s)')
plt.ylabel('acceleration (m/s^2)')
plt.title('Detection on Acceleration')

plt.figure(5)
plt.plot(gyro_t, gyro_x)
plt.plot(gyro_t, gyro_y)
plt.plot(gyro_t, gyro_z)
plt.plot(acc_t,std_static*np.max(gyro_x))
plt.legend(['gyro x','gyro y','gyro z','detector'], loc = 'upper right')
plt.xlabel('time (s)')
plt.ylabel('angular velocity (rad/s)')

windows = []

for i in range(len(std_static)-1):
    if (std_static[i+1] - std_static[i]) == 1:
        start = i + 1
       
    if (std_static[i+1] - std_static[i]) == -1:
        endd = i
    
        if (endd - start) > 15:
                windows.append(start)
                windows.append(endd)
                
res = np.zeros((int(len(windows)/2),3))

for i in range(len(res)):
    print('window1', [windows[2*i], windows[2*i+1]])
    res[i,0] = np.mean(acc_x[windows[2*i]:windows[2*i+1]])
    res[i,1] = np.mean(acc_y[windows[2*i]:windows[2*i+1]])
    res[i,2] = np.mean(acc_z[windows[2*i]:windows[2*i+1]])
    
    
import scipy.io as sio
sio.savemat('res.mat', {'res':res})