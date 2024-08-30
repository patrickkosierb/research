from torch import nn as nn
import torch
import torch.autograd
import numpy as np
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import pandas as pd
import os
'''
Establish Neural network class for outlining the architecure employed, used as constructor for neural network object:
input size = 1: input is select collocation points of time data
hidden layers = 5, hidden nodes = 10: Arbitrarily decided, hyperparameters needs to be checked via validation (grid search etc)
output layer size = 2: estimation of dependent variables x1 and x2
'''
class UNN(nn.Module):
    def __init__(self, input_size = 3, output_size = 2, hidden_layers = 6, hidden_nodes = 10):
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

data = pd.read_csv(dataPath+'v250000_a50_121mv2_clean.csv')
t = data['Time'].to_numpy()
x1_True = data['Position'].to_numpy()
x2_True = data['Velocity'].to_numpy()

# Savitzky-Golay filter on xdot
window_size = 31
poly_order = 8    # Typically 2 to 4
x2_True = savgol_filter(x2_True, window_length=window_size, polyorder=poly_order)

f_in = data['Force'].to_numpy()

# lls
c_lls = 8.86846727390849
k_lls = 0.597632140788303
alpha_lls = 30.837990991911244
beta_lls = -954.4497438208226
gamma_lls = 972.7733104878047
A_lls = 7.622000649825172
x0_lls = 5.388211694664788

# init
c_init=1.0
k_init=1.0
alpha_init=1.0
beta_init=1.0
gamma_init=1.0
A_init=1.0
x0_init=1.0
n_po = 2

tsample = 1000
t_temp = t

t = t.reshape(-1,1)
t_data = t[0:tsample:10] # use 40% of total range as training data, take 10 points evenly from range
t_train = torch.tensor(t_data, dtype=torch.float32, device="cpu") #convert to torch tensor

x1_True =100*x1_True.reshape(-1,1)
x1_data = x1_True[0:tsample:10] # use 40% of total range as training data
x1_train = torch.tensor(x1_data, dtype=torch.float32, device="cpu") #convert to torch tensor

x2_True =100*x2_True.reshape(-1,1)
x2_data = x2_True[0:tsample:10] # use 40% of total range as training data
x2_train = torch.tensor(x2_data, dtype=torch.float32, device="cpu") #convert to torch tensor

x_train = torch.cat((x1_train,x2_train), 1) 

u = f_in[0:tsample:10]
u_data = u.reshape(-1,1)
u_train = torch.tensor(u_data, dtype=torch.float32, device="cpu") #convert to torch tensor

c_est = nn.Parameter(torch.tensor([1.0], requires_grad=True).float())
k_est = nn.Parameter(torch.tensor([1.0 ], requires_grad=True).float())
alpha_est = nn.Parameter(torch.tensor([1.0], requires_grad=True).float())
beta_est = nn.Parameter(torch.tensor([1.0], requires_grad=True).float())
gamma_est = nn.Parameter(torch.tensor([1.0], requires_grad=True).float())
A_est = nn.Parameter(torch.tensor([1.0], requires_grad=True).float())
x0_est = nn.Parameter(torch.tensor([1.0], requires_grad=True).float())
n = 2

pinn = UNN() # using the constructor, construct neural network object
pinn.register_parameter('c', c_est) #register parameter to be optimized in the neural network
pinn.register_parameter('k', k_est) #register parameter to be optimized in the neural network
pinn.register_parameter('alpha', alpha_est) #register parameter to be optimized in the neural network
pinn.register_parameter('beta', beta_est) #register parameter to be optimized in the neural network
pinn.register_parameter('gamma', gamma_est) #register parameter to be optimized in the neural network
pinn.register_parameter('A', A_est) #register parameter to be optimized in the neural network
pinn.register_parameter('x0', x0_est) #register parameter to be optimized in the neural network

optim = torch.optim.Adam(pinn.parameters(), lr=0.0001) # using the Adam algorithm, for first-order gradient-based optimization of stochastic objective functions
mse_loss = nn.MSELoss() # mean squared error object construction

t_physics = torch.linspace(0, 2, 100, requires_grad=True, device="cpu", dtype=torch.float32).reshape(-1,1) # range for prediction for ODE loss eval

torch.autograd.set_detect_anomaly(True)

checkpoint = torch.load('checkpoint.pth')
pinn.load_state_dict(checkpoint['model_state_dict'])
optim.load_state_dict(checkpoint['optimizer_state_dict'])
epoch = checkpoint['epoch']

nn_train = torch.cat((x1_train,x2_train,t_physics), 1) 

epochs = 50000 # train 50000 iterations
loss_over_training = []
for epoch in range(epochs):

    optim.zero_grad()
    
    # use NN to predict output over entire domain of interest
    nn_out = pinn(nn_train)
    u_nn, z_nn = torch.split(nn_out.t(),1)

    loss_d = u_train - u_nn
    loss_fmr = u_train - (c_est*x2_train+k_est*(x1_train-x0_est)+alpha_est*z_nn) 
    
    dzdt = torch.autograd.grad(z_nn.sum(), t_physics, create_graph=True)[0]
    loss_dz = dzdt-(-gamma_est*abs(x2_train)*z_nn*pow(abs(z_nn),n-1)-beta_est*x2_train*pow(abs(z_nn),n)+A_est*x2_train)

    fmrLoss = torch.mean(loss_fmr**2)/len(loss_fmr)
    dzLoss = torch.mean(loss_dz**2)/len(loss_dz)
    nnLoss = torch.mean(loss_d**2)/len(loss_fmr)

    totalLoss = fmrLoss + dzLoss + nnLoss 

    totalLoss.backward()
    optim.step() # optimize using the specified algorithm and learning rate

    if epoch % 1000 == 0: # Print out the error(loss) and estimated parameters every 1000 optimization cycles
        print(f"Epochs = {epoch} of {epochs}, Loss = {float(totalLoss):.7f}, params = {c_est.item(), k_est.item(),alpha_est.item() ,beta_est.item(),gamma_est.item(),A_est.item(),x0_est.item()}")
        loss_over_training.append(totalLoss)
        # Save checkpoint
        checkpoint = {
            'epoch': epoch,
            'model_state_dict': pinn.state_dict(),
            'optimizer_state_dict': optim.state_dict(),
            # You can add more information such as loss, metrics, etc.
        }
        torch.save(checkpoint, 'checkpoint.pth')

torch.save(pinn.state_dict(), 'final_model.pth')

# z_pred = 0.1
# z_lls = 0.1
# f_pred = []
# f_lls = []

# for i in range(len(t_data)):
        # dt = t_temp[i]
        # dzdt = -gamma_est*abs(x2_train[i])*z_pred*pow(abs(z_pred),n-1)-beta_est*x2_train[i]*pow(abs(z_pred),n)+A_est*x2_train[i]
        # dzdt_lls = -gamma_lls*abs(x2_train[i])*z_lls*pow(abs(z_lls),n-1)-beta_lls*x2_train[i]*pow(abs(z_lls),n)+A_lls*x2_train[i]
        
        # z_pred += dzdt*dt
        # z_lls += dzdt_lls*dt

        # f_pred.append(c_est*x2_train+k_est*(x1_train[i]-x0_est)+alpha_est*z_nn)
        # f_lls.append(c_lls*x2_train+k_lls*(x1_train[i]-x0_lls)+alpha_lls*z_lls)
# print(f_lls)
# f_lls = np.array(f_lls)
# print("RMSE LLS: ",rmse(f_lls,u))
# print("RMSE PINN: ",rmse(f_pred,u))

# plt.plot(t_data,f_lls,label="LLS Force")
# # plt.plot(t_data,f_pred,label="Predicted Force")
# plt.plot(t_data,u, label="Actual Force")
# plt.legend()
# cwd = os.getcwd()
# plt.savefig(cwd+"/forcePINN.jpg")
    

""" METHOD 2 """ 
    # pred = pinn(x_train).t()
    # loss_d = u_train - pred 
    # z_pred=torch.zeros(len(t_data)).reshape(-1,1)
    # dzdt=0
    # n=2
    # z_pred[0] = 0.1
    
    # for i in range(len(t_data)):
    #     dt = t_temp[i]
    #     dzdt = -gamma_est*abs(x2_train[i])*z_pred[i].clone()*pow(abs(z_pred[i].clone()),n-1)-beta_est*x2_train[i]*pow(abs(z_pred[i].clone()),n)+A_est*x2_train[i]
    #     z_pred[i] += dzdt*dt

"""""" 