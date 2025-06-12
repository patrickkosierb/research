import time
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import savgol_filter
from scipy.optimize import leastsq
from geneticalgorithm import geneticalgorithm as ga

''' Helpers '''
def importData(path):
    MR_data = pd.read_csv(path)
    t, x, dxdt, volt, f = MR_data['Time'].to_numpy(), MR_data['Position'].to_numpy(), MR_data['Velocity'].to_numpy(), MR_data['Voltage'].to_numpy(), MR_data['Force'].to_numpy()
    # Savitzky-Golay filter on xdot
    window_size = 31
    poly_order = 8    # Typically 2 to 4
    xdot = savgol_filter(dxdt, window_length=window_size, polyorder=poly_order)
    return t, x, xdot, volt, f

def sign(z):
    return np.where(z < 0, -1, np.where(z == 0, 0, 1))

def rmse(y_true, y_pred):
    return np.sqrt(np.mean((y_true - y_pred) ** 2))

'''' Bingham '''
class Bingham:
    def __init__(self,measured,position,velocity,time):
        self.f = measured
        self.x = position
        self.xdot = velocity
        self.t = time
        
    def calcBingham(self,cd, fc, f0):
        force = []

        for i in range(len(self.t) - 1):
            force.append(cd*self.xdot[i]+fc*sign(self.xdot[i])+f0) if not (np.isinf(cd) or np.isinf(fc) or np.isinf(f0)) else force.append(0) 

        return np.array(force)

    def residualBingham(self,params):
        cd, fc, f0 = params
        return self.f[:-1] - self.calcBingham(cd,fc, f0)

    def llsBingham(self):
        initial = [0,0,0]
        return leastsq(self.residualBingham, initial)

    def plotBingham(self,param):
        cd, fc, f0 = param
        f_r = self.calcBingham(cd, fc, f0)
        print("RMSE Bingham: "+ str(rmse(self.f[0:len(f)-1],f_r)))
        plt.plot(f_r)
        plt.plot(self.f)
        plt.legend(["theoretical", "measured"])
        plt.savefig("lls_results/Bingham-lls.png")
    
''''Bouc-Wen '''
class BoucWen:
    def __init__(self,measured,position,velocity,time):
        self.f = measured
        self.x = position
        self.xdot = velocity
        self.t = time

    def calcBoucWenSimple(self,c, k, beta, gamma):
        force = []
        z=0.1
        dzdt=0
        n=2

        for i in range(len(self.t) - 1):
            dt = (self.t[i+1] - self.t[i])
            if np.abs(z) != 0:  # Check for division by zero or invalid values
                term1 = beta*abs(self.xdot[i])*z*pow(abs(z),n-1)
                term2 = gamma*self.xdot[i]*pow(abs(z),n)
                if not (np.isinf(term1) or np.isinf(term2)):
                    dzdt = self.xdot[i]-term1-term2
                else:
                    dzdt = 0  # Handle overflow
            else:
                dzdt = 0

            z += dzdt*dt
            force.append(c*self.xdot[i]+k*z)

        return np.array(force)
    

    def residualBoucWenSimple(self, params):
        c,k,beta,gamma = params
        return self.f[:-1] - self.calcBoucWenSimple(c,k, beta, gamma)

    def llsBoucWenSimple(self):
        initial = [1,1,1,1]
        return leastsq(self.residualBoucWenSimple, initial)

    def plotBoucWenSimple(self,params):
        c, k, beta, gamma = params
        f_r = self.calcBoucWenSimple(c,k,beta,gamma)
        print("RMSE Bouc-Wen Simple: "+ str(rmse(self.f[0:len(f)-1],f_r)))
        plt.plot(f_r)
        plt.plot(self.f)
        plt.legend(["theoretical", "measured"])
        plt.savefig("lls_results/BoucWenSimple-lls.png")

    def calcBoucWen(self,c, k, alpha, beta, gamma, A, x0): #making a standard BW one aswell with more param
        force = []
        z=0.1
        dzdt=0
        n=2

        for i in range(len(self.t) - 1):
            dt = (self.t[i+1] - self.t[i])
            if np.abs(z) != 0:  # Check for division by zero or invalid values
                term1 = gamma*abs(self.xdot[i])*z*pow(abs(z),n-1)
                term2 = beta*self.xdot[i]*pow(abs(z),n)
                term3 = A*self.xdot[i]
                if not (np.isinf(term1) or np.isinf(term2) or np.isinf(term3)):
                    dzdt = -term1-term2+term3
                else:
                    dzdt = 0  # Handle overflow
            else:
                dzdt = 0

            z += dzdt*dt
            force.append(c*self.xdot[i]+k*(self.x[i]-x0)+alpha*z)

        return np.array(force)
    
    def residualBoucWen(self, params):
        c, k, alpha, beta, gamma, A, x0 = params
        return self.f[:-1] - self.calcBoucWen(c, k, alpha, beta, gamma, A, x0)

    def llsBoucWen(self):
        initial = [1,1,1,1,1,1,1]
        return leastsq(self.residualBoucWen, initial)
    
    def fitnessFunction(self, params):
        c, k, alpha, beta, gamma, A, x0 = params
        loss = rmse(self.calcBoucWen(c, k, alpha, beta, gamma, A, x0),(self.f[:-1]))
        return loss

    def gaBoucWen(self):

        # param
        c_range = [-1000, 1000]
        k_range = [-1000, 1000]
        alpha_range = [-100, 100]
        beta_range = [-100, 100]
        gamma_range = [-100, 100]
        A_range =[-1000, 1000]
        x0_range =[-1000, 1000]

        # ga bounds
        varbound=np.array([c_range, k_range,alpha_range,beta_range,gamma_range,A_range,x0_range])
        vartype=np.array([['real'],['real'],['real'],['real'],['real'],['real'],['real']])
        num_hp=len(vartype)

        algorithm_param = {'max_num_iteration': 200,\
                   'population_size':50,\
                   'mutation_probability':0.4,\
                   'elit_ratio': 0.01,
                   'crossover_probability': 0.3,\
                   'parents_portion': 0.3,\
                   'crossover_type':'uniform',\
                   'max_iteration_without_improv':None}

        ga_model=ga(function=self.fitnessFunction, dimension=num_hp,variable_type_mixed=vartype,variable_boundaries=varbound, algorithm_parameters=algorithm_param, function_timeout=600)
        ga_model.run()
        return ga_model

    def plotBoucWenLLS(self,params):
        c, k, alpha, beta, gamma, A, x0 = params
        f_r = self.calcBoucWen(c, k, alpha, beta, gamma, A, x0)
        print("RMSE Bouc-Wen: "+ str(rmse(self.f[0:len(f)-1],f_r)))

        plt.plot(self.t[0:len(f)-1],f_r)
        plt.plot(self.t[0:len(f)-1],self.f[0:len(f)-1])
        plt.xlabel("Time(ms)")
        plt.ylabel("Force(N)")
        plt.legend(["LLS" ,"Actual Force"], loc="upper right")
        plt.savefig("results/BoucWen-LLS.png")
        plt.show()

    def plotBoucWenGA(self,params):
        c, k, alpha, beta, gamma, A, x0 = params
        f_r = self.calcBoucWen(c, k, alpha, beta, gamma, A, x0)
        print("RMSE Bouc-Wen: "+ str(rmse(self.f[0:len(f)-1],f_r)))

        plt.plot(self.t[0:len(f)-1],f_r)
        plt.plot(self.t[0:len(f)-1],self.f[0:len(f)-1])
        plt.xlabel("Time(ms)")
        plt.ylabel("Force(N)")
        plt.legend(["GA" ,"Actual Force"], loc="upper right")
        plt.savefig("results/BoucWen-GA.png")
        plt.show()

if __name__ == "__main__":
    
    ''''Import Data'''
    # t, x, xdot, volt, f = importData("../data/v250000_a50_121mv2_clean.csv")   # v250000_a50_121mv2_clean sample:5885 # data_voltagevarying
    t, x, xdot, volt, f = importData("data/v250000_a50_121mv2_clean.csv")   # v250000_a50_121mv2_clean sample:5885 # data_voltagevarying
    # index = 39418 # for data_voltagevarying.csv when looking at 0.5V
    index = 100

    t = t[index:len(t)-1]
    x = x[index:len(x)-1]*100
    xdot = xdot[index:len(xdot)-1]*100
    volt = volt[index:len(volt)-1]
    f = f[index:len(f)-1]
    
    # Boucwen Standard 
    bw = BoucWen(f,x,xdot,t)

    # run nlls
    # start_time = time.time()
    # optimalBoucWenLLS = bw.llsBoucWen()[0]
    # print("c: "+str(optimalBoucWenLLS[0])+"\nk: "+str(optimalBoucWenLLS[1])+"\nAlpha: "+str(optimalBoucWenLLS[2])+"\nBeta: "+str(optimalBoucWenLLS[3])+"\nGamma: "+str(optimalBoucWenLLS[4])+"\nA: "+str(optimalBoucWenLLS[5])+"\nx0: "+str(optimalBoucWenLLS[6]))
    # end_time = time.time()
    # execution_time = end_time - start_time
    # print("Execution time Bouc-Wen NLLS:", execution_time, "seconds")
    # bw.plotBoucWenLLS(optimalBoucWenLLS)
    
    # run ga
    start_time = time.time()
    optimalBoucWenGA = bw.gaBoucWen()
    gaParam = optimalBoucWenGA.output_dict["variable"]
    print("c: "+str(gaParam[0])+"\nk: "+str(gaParam[1])+"\nAlpha: "+str(gaParam[2])+"\nBeta: "+str(gaParam[3])+"\nGamma: "+str(gaParam[4])+"\nA: "+str(gaParam[5])+"\nx0: "+str(gaParam[6]))
    end_time = time.time()
    execution_time = end_time - start_time
    print("Execution time Bouc-Wen GA:", execution_time, "seconds")
    bw.plotBoucWenGA(gaParam)


    # # plot data
    # fig, axs = plt.subplots(4, 1, figsize=(10, 10))

    # axs[0].plot(t, x)
    # axs[0].set_title('Position (x)')
    # axs[0].set_xlabel('Time (t)')
    # axs[0].set_ylabel('x')

    # axs[1].plot(t, xdot)
    # axs[1].set_title('Velocity (xdot)')
    # axs[1].set_xlabel('Time (t)')
    # axs[1].set_ylabel('xdot')

    # axs[2].plot(t, volt)
    # axs[2].set_title('Voltage (volt)')
    # axs[2].set_xlabel('Time (t)')
    # axs[2].set_ylabel('volt')

    # axs[3].plot(t, f)
    # axs[3].set_title('Force (f)')
    # axs[3].set_xlabel('Time (t)')
    # axs[3].set_ylabel('f')

    # plt.tight_layout()
    # plt.show()
    # plt.clf() 


