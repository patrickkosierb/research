from torch import nn as nn
import torch
import torch.autograd
import numpy as np
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import pandas as pd
import os
import time
'''
Establish Neural network class for outlining the architecure employed, used as constructor for neural network object:
input size = 1: input is select collocation points of time data
hidden layers = 5, hidden nodes = 10: Arbitrarily decided, hyperparameters needs to be checked via validation (grid search etc)
output layer size = 2: estimation of dependent variables x1 and x2
'''
class UNN(nn.Module):
    def __init__(self, input_size = 1, output_size = 1, hidden_layers = 6, hidden_nodes = 10):
        super().__init__()
        self.input_layer  = nn.Linear(input_size, hidden_nodes) # first layer: input to hidden size: ie 1 to 10
        self.hidden_layers_list = nn.ModuleList([nn.Linear(hidden_nodes, hidden_nodes)] * hidden_layers) # hidden layers, hidden to hidden layer: ie 10 to 10
        self.output_layer = nn.Linear(hidden_nodes, output_size) # hidden to output layer: 10 to 2 (prediction)
        self.activation_function = nn.Tanh() # Construct activation fxn

    def forward(self, t):
        t = self.activation_function(self.input_layer(t)) #Input to first hidden layer
        for hidden_layer in self.hidden_layers_list:  #Passes the previous layer input through a hidden layer, followed by an activation, repeating for all hidden layers
            t = self.activation_function(hidden_layer(t))
        t = self.output_layer(t) # last hidden layer to output, no activations applied as this is the prediction
        return t

def rmse(y_true, y_pred):
    return np.sqrt(np.mean((y_true - y_pred) ** 2))


dataPath = "/mnt/c/Users/patri/OneDrive/Documents/papers/parameterID/data/"

x_True = pd.read_csv(dataPath+'v250000_a50_121mv2_clean.csv')
t = x_True['Time'].to_numpy()
x1_True = x_True['Position'].to_numpy()
x2_True = x_True['Velocity'].to_numpy()

# Savitzky-Golay filter on xdot
window_size = 31
poly_order = 8    # Typically 2 to 4
x2_True = savgol_filter(x2_True, window_length=window_size, polyorder=poly_order)

f_in = x_True['Force'].to_numpy()
   
tsample = 5000
t = t.reshape(-1,1)
t_data = t[0:tsample:5] # use 40% of total range as training data, take 10 points evenly from range
t_train = torch.tensor(t_data, dtype=torch.float32, device="cpu") #convert to torch tensor

x2_True = x2_True.reshape(-1,1)
x2_data = x2_True[0:tsample:5] # use 40% of total range as training data
x2_train = torch.tensor(x2_data, dtype=torch.float32, device="cpu") #convert to torch tensor
x2_tensor = torch.tensor(x2_True.reshape(-1,1), dtype=torch.float32, device="cpu") #convert to torch tensor

u = f_in[0:tsample:5]
u_data = u.reshape(-1,1)
u_train = torch.tensor(u_data, dtype=torch.float32, device="cpu") #convert to torch tensor
u_tensor = torch.tensor(f_in.reshape(-1,1), dtype=torch.float32, device="cpu") #convert to torch tensor


# Parameters

# initial 
cd_init = 1.0
fc_init = 1.0
f0_init=  1.0

# 200000 epoch checkpoint 
# cd_init = 11.389698028564453
# fc_init = 16.670991897583008
# f0_init=  -5.08684778213501


cd_est = nn.Parameter(torch.tensor([cd_init], requires_grad=True).float())
fc_est = nn.Parameter(torch.tensor([fc_init ], requires_grad=True).float())
f0_est = nn.Parameter(torch.tensor([f0_init], requires_grad=True).float())

pinn = UNN() # using the constructor, construct neural network object
pinn.register_parameter('cd', cd_est) #register parameter to be optimized in the neural network
pinn.register_parameter('fc', fc_est) #register parameter to be optimized in the neural network
pinn.register_parameter('f0', f0_est) #register parameter to be optimized in the neural network

optim = torch.optim.Adam(pinn.parameters(), lr=0.0001) # using the Adam algorithm, for first-order gradient-based optimization of stochastic objective functions
mse_loss = nn.MSELoss() # mean squared error object construction

x2_physics = torch.linspace(0, 1, 1000, requires_grad=True, device="cpu", dtype=torch.float32).reshape(-1,1) # range for prediction for ODE loss eval

epochs = 150000 # train 50000 iterations
loss_over_training = []
start_time = time.time()
for epoch in range(epochs):

    optim.zero_grad()
    u_pred1 = pinn(x2_physics)
    u_pred1 = u_pred1.t()
    #use NN to predict output over entire domain of interest
    u_pred = pinn(x2_train) 
    u_pred = u_pred.t()
     
    # make sure components are balanced and same units
    loss_d1 = u_train-u_pred1
    loss_d = u_train - u_pred 
    loss_phys = u_train - (100*cd_est*x2_train + fc_est*torch.sign(x2_train) + f0_est) # Bingham
    
    physicsLoss = torch.mean(loss_phys**2)/len(loss_phys)
    dataLoss = torch.mean(loss_d**2)/len(loss_d)
    dataLoss1 = torch.mean(loss_d1**2)/len(loss_d1)

    totalLoss = physicsLoss + dataLoss + dataLoss1

    totalLoss.backward()
    optim.step() # optimize using the specified algorithm and learning rate

    if epoch % 1000 == 0: # Print out the error(loss) and estimated parameters every 1000 optimization cycles
        print(f"Epochs = {epoch} of {epochs}, Loss = {float(totalLoss):.4f}, params = {cd_est.item(), fc_est.item(), f0_est.item()}")
        loss_over_training.append(totalLoss)

end_time = time.time()
cd_lls= 12.19649101898022
fc_lls= 14.816658010248124
f0_lls= -5.029104318005954

f_lls = (100*cd_lls*x2_data + fc_lls*np.sign(x2_data) + f0_lls)
# pred = [12.046374320983887, 15.031420707702637, -5.047383785247803] 
# f_pred = (100*pred[0]*x2_data + pred[1]*np.sign(x2_data) + pred[2])
f_pred = (100*cd_est.item()*x2_data + fc_est.item()*np.sign(x2_data) + f0_est.item())

# f = open("checkpoint_bingham.txt", "a")
# f.write("param=",cd_est.item(), fc_est.item(), f0_est.item())

print("RMSE LLS: ",rmse(f_lls,u))
print("RMSE after training: ",rmse(f_pred,u))
print("Elapsed: ", end_time-start_time)
plt.plot(t_data,u, label="Actual Force")
plt.plot(t_data,f_pred,label="PINN Force")
plt.plot(t_data,f_lls,label="LLS Force", alpha=0.5)
plt.xlabel("Time(s)")
plt.ylabel("Force(N)")
plt.legend(loc="upper right")
cwd = os.getcwd()
plt.savefig(cwd+"/forcePINN.jpg")