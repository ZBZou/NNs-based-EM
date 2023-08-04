import torch
import torch.nn as nn
import random
import math
import numpy as np
import cmath
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import torch.nn.utils as utils


def set_seed(seed):
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)

set_seed(0)

SampleSize = 1024 * 5

#Create 4x5x4x5 kernal
NumberOfUsers = 2
TimeLength = 3
taps_max = 3
taps_min = 1

Taps = torch.randint(taps_min, taps_max, (NumberOfUsers, TimeLength))
kernal = torch.zeros(SampleSize, NumberOfUsers, TimeLength, NumberOfUsers, TimeLength, dtype=torch.cfloat)
kernal_split = torch.zeros(SampleSize, NumberOfUsers, TimeLength, NumberOfUsers, TimeLength, 2)
        
for i in range(SampleSize):
    for u1 in range(NumberOfUsers):
        for t1 in range(TimeLength):
            for u2 in range(NumberOfUsers):
                for t2 in range(TimeLength):
                    taps = Taps[u2,t1];
                    if (t2 > t1) or (t2 < t1 - taps) :
                        continue
                    else :
                        val = complex(torch.randn(1), torch.randn(1))
                        kernal[i, u1,t1,u2,t2] = val
                        kernal_split[i, u1,t1,u2,t2,0] = val.real
                        kernal_split[i, u1,t1,u2,t2,1] = val.imag

numpy_array = kernal_split.numpy()
np.save('kernal_2x3x2x3.npy', numpy_array)

#Create 4x5x4x5 kernal
NumberOfUsers = 4
TimeLength = 5
taps_max = 5
taps_min = 3

Taps = torch.randint(taps_min, taps_max, (NumberOfUsers, TimeLength))
kernal = torch.zeros(SampleSize, NumberOfUsers, TimeLength, NumberOfUsers, TimeLength, dtype=torch.cfloat)
kernal_split = torch.zeros(SampleSize, NumberOfUsers, TimeLength, NumberOfUsers, TimeLength, 2)
        
for i in range(SampleSize):
    for u1 in range(NumberOfUsers):
        for t1 in range(TimeLength):
            for u2 in range(NumberOfUsers):
                for t2 in range(TimeLength):
                    taps = Taps[u2,t1];
                    if (t2 > t1) or (t2 < t1 - taps) :
                        continue
                    else :
                        val = complex(torch.randn(1), torch.randn(1))
                        kernal[i, u1,t1,u2,t2] = val
                        kernal_split[i, u1,t1,u2,t2,0] = val.real

numpy_array = kernal_split.numpy()
np.save('kernal_4x5x4x5.npy', numpy_array)

#Create 4x10x4x10 kernal
NumberOfUsers = 4
TimeLength = 10
taps_max = 6
taps_min = 2

Taps = torch.randint(taps_min, taps_max, (NumberOfUsers, TimeLength))
kernal = torch.zeros(SampleSize, NumberOfUsers, TimeLength, NumberOfUsers, TimeLength, dtype=torch.cfloat)
kernal_split = torch.zeros(SampleSize, NumberOfUsers, TimeLength, NumberOfUsers, TimeLength, 2)
        
for i in range(SampleSize):
    for u1 in range(NumberOfUsers):
        for t1 in range(TimeLength):
            for u2 in range(NumberOfUsers):
                for t2 in range(TimeLength):
                    taps = Taps[u2,t1];
                    if (t2 > t1) or (t2 < t1 - taps) :
                        continue
                    else :
                        val = complex(torch.randn(1), torch.randn(1))
                        kernal[i, u1,t1,u2,t2] = val
                        kernal_split[i, u1,t1,u2,t2,0] = val.real
                        kernal_split[i, u1,t1,u2,t2,1] = val.imag

numpy_array = kernal_split.numpy()
np.save('kernal_4x10x4x10.npy', numpy_array)


#Create 2x10x2x10 kernal
NumberOfUsers = 2
TimeLength = 10
taps_max = 6
taps_min = 2

Taps = torch.randint(taps_min, taps_max, (NumberOfUsers, TimeLength))
kernal = torch.zeros(SampleSize, NumberOfUsers, TimeLength, NumberOfUsers, TimeLength, dtype=torch.cfloat)
kernal_split = torch.zeros(SampleSize, NumberOfUsers, TimeLength, NumberOfUsers, TimeLength, 2)
        
for i in range(SampleSize):
    for u1 in range(NumberOfUsers):
        for t1 in range(TimeLength):
            for u2 in range(NumberOfUsers):
                for t2 in range(TimeLength):
                    taps = Taps[u2,t1];
                    if (t2 > t1) or (t2 < t1 - taps) :
                        continue
                    else :
                        val = complex(torch.randn(1), torch.randn(1))
                        kernal[i, u1,t1,u2,t2] = val
                        kernal_split[i, u1,t1,u2,t2,0] = val.real
                        kernal_split[i, u1,t1,u2,t2,1] = val.imag

numpy_array = kernal_split.numpy()
np.save('kernal_2x10x2x10.npy', numpy_array)


#Create 6x10x6x10 kernal
NumberOfUsers = 6
TimeLength = 10
taps_max = 6
taps_min = 2

Taps = torch.randint(taps_min, taps_max, (NumberOfUsers, TimeLength))
kernal = torch.zeros(SampleSize, NumberOfUsers, TimeLength, NumberOfUsers, TimeLength, dtype=torch.cfloat)
kernal_split = torch.zeros(SampleSize, NumberOfUsers, TimeLength, NumberOfUsers, TimeLength, 2)
        
for i in range(SampleSize):
    for u1 in range(NumberOfUsers):
        for t1 in range(TimeLength):
            for u2 in range(NumberOfUsers):
                for t2 in range(TimeLength):
                    taps = Taps[u2,t1];
                    if (t2 > t1) or (t2 < t1 - taps) :
                        continue
                    else :
                        val = complex(torch.randn(1), torch.randn(1))
                        kernal[i, u1,t1,u2,t2] = val
                        kernal_split[i, u1,t1,u2,t2,0] = val.real
                        kernal_split[i, u1,t1,u2,t2,1] = val.imag

numpy_array = kernal_split.numpy()
np.save('kernal_6x10x6x10.npy', numpy_array)


#Create 8x10x8x10 kernal
NumberOfUsers = 8
TimeLength = 10
taps_max = 6
taps_min = 2

Taps = torch.randint(taps_min, taps_max, (NumberOfUsers, TimeLength))
kernal = torch.zeros(SampleSize, NumberOfUsers, TimeLength, NumberOfUsers, TimeLength, dtype=torch.cfloat)
kernal_split = torch.zeros(SampleSize, NumberOfUsers, TimeLength, NumberOfUsers, TimeLength, 2)
        
for i in range(SampleSize):
    for u1 in range(NumberOfUsers):
        for t1 in range(TimeLength):
            for u2 in range(NumberOfUsers):
                for t2 in range(TimeLength):
                    taps = Taps[u2,t1];
                    if (t2 > t1) or (t2 < t1 - taps) :
                        continue
                    else :
                        val = complex(torch.randn(1), torch.randn(1))
                        kernal[i, u1,t1,u2,t2] = val
                        kernal_split[i, u1,t1,u2,t2,0] = val.real
                        kernal_split[i, u1,t1,u2,t2,1] = val.imag

numpy_array = kernal_split.numpy()
np.save('kernal_8x10x8x10.npy', numpy_array)

    
    
