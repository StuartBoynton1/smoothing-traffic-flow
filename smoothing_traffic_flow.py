from lqr_sdp import *
import array as arr
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


#run this command before running mosek
#python <MSKHome>/mosek/9.2/tools/platform/osx64x86/bin/install.py


# Mix or not
mix = 1
# 1.Optimal Control  2.FollowerStopper  3.PI with Saturation
controllerType = 1

brakeID = 6 - 1
# In the simulation, the ID of the AV is 20, and thus the brakeID needs to minus one    

#Chose number of Vehicles
N = 20

#Choose circumference of track (0) or equilibrium spacing of vehicles s_star (1)
#If circumference is chosen then you may ignore value of spacing and same for spacing
circumference_or_spacing = 0
circumference = 400
s_star = 20



if mix == 1 and controllerType == 1:
    gammaType = 2
    if gammaType == 1:
        gamma_s = 0.03
        gamma_v = 0.15
        gamma_u = 1
    elif gammaType == 2:
        gamma_s = 3
        gamma_v = 15
        gamma_u = 1
    
    #K = np.array( [[3.05214216883394, 2.04047486745414, 4.23441918486731, 0.816748753025211, 4.27938838608577, -0.0841721905231217, 3.58825708391328, -0.622013145554968, 2.55870322570934, -0.831571241306728, 1.51118202880832, -0.799391067556516, 0.648096599208884, -0.632100576735584, 0.0499895764298799, -0.426406463994733, -0.298173114761084, -0.249029740041272, -0.465166658625890, -0.130429139967252, -0.532302845404502, -0.0708478026818128, -0.564848465218233, -0.0532302254917162, -0.600421176828263, -0.0564142781996425, -0.651376872199535, -0.0638967241144484, -0.714263700496758, -0.0670317867827030, -0.779984326687099, -0.0648390708889171, -0.842126285544850, -0.0638888616892629, -0.905090752213027, -0.0796764272045653, -0.993965118029564, -0.131210521554255, -1.31725756289723, 5.23568807061450]])
    K = lqr_sdp(N,s_star,gamma_s,gamma_v,gamma_u,1)
    print(K)
    
alpha_k = 0.6
    
v_max = 30
acel_max = 5
dcel_max = -5
#Driver Model: OVM
alpha = 0.6
beta = 0.9
s_st = 5
s_go = 35


TotalTime = 100
Tstep = 0.01
NumStep = int(TotalTime/Tstep)


if mix:
    ActuationTime = 0  
else:
    ActuationTime = 9999



if circumference_or_spacing:
    circumference = s_star*N
elif circumference_or_spacing == 0:
    s_star = circumference/N

v_star  = (v_max/2) * (1-math.cos(math.pi * (s_star - s_st)/(s_go - s_st)))
s_ctr  = s_star
v_ctr  = (v_max/2) * (1-math.cos(math.pi * (s_ctr - s_st)/(s_go - s_st)))              #What is s_ctr and v_ctr


sd = 0 # Collision avoidance safe distance

S = np.zeros((NumStep,N,3))

#Initial state for each vehicle
dev_s = 0
dev_v = 0
co_v = 1.0
v_ini = co_v * v_star                                                                 # What is co_v ?????

var1 = np.linspace(circumference, s_star, N)
var2 = np.random.rand(N)*2*dev_s-dev_s

S[0, :, 0] = var1 + var2

#print(S[0,:,0])

var1 = v_ini*np.ones([N])
var2 = (np.random.rand(N)*2*dev_v-dev_v)
S[0, :, 1] =  var1 + var2

ID = np.zeros([N])

if mix == 1:
    ID[-1] = 1
    X = np.zeros([2*N, NumStep])


#Velocity difference
V_diff = np.zeros([NumStep,N])

#Following Distance
D_diff = np.zeros([NumStep,N])
temp = np.zeros([N])


#Average Speed
V_avg = np.zeros([NumStep,1])
v_cmd = np.zeros([NumStep,1])   #For controller's 2 and 3



##Theoretical Analysis

#...



##Simulation

sd_actuate = 0

for k in range(0,NumStep-2):
    
    #Update Acceleration
    temp[1:] = S[k,:-1,1]
    temp[0] = S[k,-1,1]

    
    V_diff[k,:] = temp-S[k,:,1]
    temp[0]=S[k,-1,0]+circumference
    temp[1:] = S[k,:-1,0]
    #print(temp)
    D_diff[k,:] = temp-S[k,:,0]
    #print(D_diff)
    cal_D = D_diff[k,:]
    cal_D[cal_D>s_go] = s_go
    cal_D[cal_D<s_st] = s_st
    
    
    #OVM Model
    
    acel2 = math.pi*(cal_D-s_st)/(s_go-s_st)
    acel1 = (1-np.cos(acel2))
    acel = alpha*(v_max/2*acel1-S[k,:,1])+beta*V_diff[k,:]
    acel[acel>acel_max] = acel_max
    acel[acel<dcel_max] = dcel_max
    
    #SD as ADAS to prevent crash
    temp[1:] = S[k,:-1,1]
    temp[0] = S[k,-1,1]
    acel_sd = (S[k,:,1]**2-temp**2)/2/D_diff[k,:]
    acel[acel_sd>abs(dcel_max)] = dcel_max
    
    S[k,:,2] = acel;
    
    
    
    if (k*Tstep>20) and (k*Tstep<22):
        S[k,brakeID,2]=-5
        
    
    if k>=ActuationTime/Tstep:
        
        if controllerType==1:
            X[np.arange(0,2*N,2),k] = D_diff[k,:]-s_ctr
            X[np.arange(1,2*N,2),k] = S[k,:,1]-v_ctr
            u = np.matmul(-K,X[:,k])
            
        elif controllerType==2:
            dx10 = 9.5
            dx20 = 10.75
            dx30 = 11
            
            dv_temp = min(S[k,N-2,1]-S[k,N-1,1],0)
            
            d1 = 1.5
            d2 = 1.0
            d3 = 0.5
            
            dx1 = dx10+dv_temp**2/2/d1
            dx2 = dx20+dv_temp**2/2/d2
            dx3 = dx30+dv_temp**2/2/d3
            
            dx = D_diff[k,N-1]
            v_temp = min(S[k,N-2,1],12)
            
            
            if dx<=1:
                v_cmd = 0
            elif dx<=d2:
                v_cmd = v_temp*(dx-dx1)/(dx2-dx1)
            elif dx<=d3:
                v_cmd = v_temp+(v_ctr-v_temp)*(dx-dx2)/(dx3-dx2)
            else:
                v_cmd = v_ctr
            
            u = alpha_k*(v_cmd-S[k,N-1,1])
            
        elif controllerType==3:
            gl = 7
            gu = 30
            v_catch = 1
            gamma_temp = 2
            
            if k-26/Tstep<=0:
                v_hisAvg = np.mean(S[:k+1,-1,1])
            else:
                v_hisAvg = np.mean(S[int(k-26/Tstep):k+1,-1,1])
            
            v_target = v_hisAvg + v_catch*min(max((D_diff[k,-1]-gl)/(gu-gl),0),1)
            alpha_temp = min(max((D_diff[k,-1]-max(2*V_diff[k,-1],4))/gamma_temp, 0),1)
            beta_temp = 1-0.5*alpha_temp
            v_cmd[k+1] = beta_temp*(alpha_temp*v_target+(1-alpha_temp)*S[k,-2,1])+(1-beta_temp)*v_cmd[k]
            u = alpha_k*(v_cmd[k+1]-S[k,-1,1])
            
        
        if u>acel_max:
            u=acel_max
        elif u<dcel_max:
            u=dcel_max
        
        
        if (S[k,N-1,1]**2-S[k,N-2,1]**2)/2/(S[k,N-2,0]-S[k,N-1,0]-sd)>abs(dcel_max):
            u=dcel_max
            
        S[k,-1,2] = u
 
        
    S[k+1,:,1] = S[k,:,1] + Tstep*S[k,:,2]
    S[k+1,:,0] = S[k,:,0] + Tstep*S[k,:,1]
    

x = np.zeros(NumStep)

for k in range(NumStep):
    V_avg[k] = np.mean(S[k,:,1])
    x[k] = k
    #print(V_avg[k])

       
       
#Display data

fig = plt.figure()

# syntax for 3-D projection
ax = plt.axes(projection ='3d')

#y = np.linspace(0,20,20)

for i in range(20):
    z = np.ones(NumStep-1)*i
    
    if i==19 and mix==1:
        ax.plot3D(z, x[:-1], S[:-1,19,1], 'red', linewidth=0.5)
        continue
    elif i==19 and mix==0:
        ax.plot3D(z, x[:-1], S[:-1,i,1], 'blue', linewidth=0.5)
        continue
    
    ax.plot3D(z, x[:-1], S[:-1,i,1], 'blue', linewidth=0.5)
         
#plt.plot(V_avg[:-2])
#plt.xlabel("Time")
#plt.ylabel("Average velocity of cars")
plt.show()
