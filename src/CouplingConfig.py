import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib.ticker import MaxNLocator

#Function to generate Christandl PST configuration
def Christandl_PST_coupling(J_max,N):
    PST_coupling = np.array([])
    if (N % 2 == 0):
        J_min = (2*J_max)/N
    else:
        J_min = (2*J_max)/(N*np.sqrt(1-(1/(N**2))))
    for i in range(1,N):
        PST_coupling = np.append(PST_coupling, J_min*np.sqrt(i*(N-i)))
    return PST_coupling

#ARGUMENTS PASSED
N = int(sys.argv[1])
custom = sys.argv[2]

# Set up the figure
plt.figure(figsize=(10,7))

# Create the axes
ax1=plt.subplot()

#Get length of argument list
coupling_length = len(sys.argv)-3
couplings = [sys.argv[i] for i in range(3,len(sys.argv))]

#Define coupling array
couplings = np.array(couplings,dtype=float)
couplings = np.reshape(couplings,(-1,N))
if (custom):
    #USING JS2D
    coupling_range = np.max(np.nonzero(couplings[0]))

    #Currently only works for 1st NN
    #But no clear way to show coupling config for second nearest neighbour anyway
    
    #---------------------------#
    #This needs to be made so each row decreases in size for smaller number of couplings for NN terms
    coupling_arr = np.zeros((coupling_range,N-1))
    #---------------------------#

    #Assuming symmetry in coupling matrix
    for j in range(1,coupling_range+1):
        for i in range(N-j):
            coupling_arr[j-1][i] = couplings[i][i+j]
else:
    #USING Js
    N = len(couplings)

#Get Christandl et al config for given N and J_max
ideal_coupling = Christandl_PST_coupling(np.max(couplings),N)

#SET SIZE OF THE AXIS TICKS
ax1.tick_params(axis='x', labelsize=20)
ax1.tick_params(axis='y', labelsize=20)

#SET X TICKS TO BE INTEGERS
ax1.xaxis.set_major_locator(MaxNLocator(integer=True))

#SET SIZE OF AXIS LABELS
ax1.set_ylabel("Normalised energy",fontsize=25)
ax1.set_xlabel('site number, $i$',fontsize=25)

#CREATE X RANGE 
sites = np.arange(1,N+1,1)

# DEFINE THE RANGE
ax1.set_xlim(1,N+1)

#PLOT THE COUPLING CONFIGURATIONS
plt.plot(sites[:-1]+0.5,coupling_arr[0],'--o',color='black',label='Optimised Coupling')
plt.plot(sites[:-1]+0.5,ideal_coupling, marker='o',color='orange',linestyle='dashdot',label='Christandl et al Coupling')
plt.tight_layout()
plt.legend()

#SAVE AS PNG PICTURE
plt.savefig('couplings.pdf',transparent=False)

plt.cla()
#PLOT AS FUNCTION OF COUPLING

for i in range(coupling_range):
    if i == 0:
        ax1.set_ylabel("Normalised energy",fontsize=25)
        ax1.set_xlabel('ith 1NN coupling',fontsize=25)
        plt.plot(sites[:-1]+0.5,coupling_arr[0],color='black',label='1st NN')
        plt.tight_layout()
        plt.legend()
        plt.savefig('1NNcouplings.pdf',transparent=False)
        plt.cla()
    elif i == 1:
        ax1.set_ylabel("Normalised energy",fontsize=25)
        ax1.set_xlabel('ith 2NN coupling',fontsize=25)
        plt.plot(sites[:-2]+0.5,coupling_arr[1][:-1],color='black',label='2nd NN')
        plt.tight_layout()
        plt.legend()
        plt.savefig('2NNcouplings.pdf',transparent=False)
        plt.cla()
    elif i == 2:
        plt.plot(sites[:-3]+0.5,coupling_arr[2][:-2],color='blue',label='3rd NN')

    else:
        plt.plot(sites[:-i]+(0.5),coupling_arr[i][:-i],label=str(i)+'th NN')


