import time
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import savgol_filter
from scipy.optimize import leastsq

''' Helpers '''
def importData(path):
    MR_data = pd.read_csv(path)
    t, x, dxdt, volt, f = MR_data['Time'].to_numpy(), MR_data['Position'].to_numpy(), MR_data['Velocity'].to_numpy(), MR_data['Voltage'].to_numpy(), MR_data['Force'].to_numpy()
    # Savitzky-Golay filter on xdot
    window_size = 31
    poly_order = 8    # Typically 2 to 4
    xdot = savgol_filter(dxdt, window_length=window_size, polyorder=poly_order)
    x=x
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

    def plotBoucWen(self,params):
        c, k, alpha, beta, gamma, A, x0 = params
        f_r = self.calcBoucWen(c, k, alpha, beta, gamma, A, x0)
        print("RMSE Bouc-Wen: "+ str(rmse(self.f[0:len(f)-1],f_r)))

        # c, k, alpha, beta, gamma, A, x0 = [22.43910026550293, 1.090404748916626, 0.8857629299163818, 0.34663310647010803, 1.1787058156187413e-06, 1.4282690286636353, 0.7886528968811035]
        # f_ml = self.calcBoucWen(c, k, alpha, beta, gamma, A, x0) #making a standard BW one aswell with more param

        # plt.plot(self.t[0:5885],f_ml)
        plt.plot(self.t[0:len(f)-1],f_r)
        plt.plot(self.t[0:len(f)-1],self.f[0:len(f)-1])
        plt.xlabel("Time(ms)")
        plt.ylabel("Force(N)")
        plt.legend(["LLS" ,"Actual Force"], loc="upper right")
        plt.savefig("lls_results/BoucWen-LLS.png")
        plt.show()

if __name__ == "__main__":
    
    ''''Import Data'''
    t, x, xdot, volt, f = importData("../data/data_voltagevarying.csv")   # v250000_a50_121mv2_clean sample:5885 # data_voltagevarying

    # plot data
    fig, axs = plt.subplots(4, 1, figsize=(10, 10))

    axs[0].plot(t, x)
    axs[0].set_title('Position (x)')
    axs[0].set_xlabel('Time (t)')
    axs[0].set_ylabel('x')

    axs[1].plot(t, xdot)
    axs[1].set_title('Velocity (xdot)')
    axs[1].set_xlabel('Time (t)')
    axs[1].set_ylabel('xdot')

    axs[2].plot(t, volt)
    axs[2].set_title('Voltage (volt)')
    axs[2].set_xlabel('Time (t)')
    axs[2].set_ylabel('volt')

    axs[3].plot(t, f)
    axs[3].set_title('Force (f)')
    axs[3].set_xlabel('Time (t)')
    axs[3].set_ylabel('f')

    plt.tight_layout()
    plt.show()
    plt.clf() 

    # Boucwen Standard 
    bw = BoucWen(f,x,xdot,t)

    start_time = time.time()

    optimalBoucWen = bw.llsBoucWen()[0]
    print("c: "+str(optimalBoucWen[0])+"\nk: "+str(optimalBoucWen[1])+"\nAlpha: "+str(optimalBoucWen[2])+"\nBeta: "+str(optimalBoucWen[3])+"\nGamma: "+str(optimalBoucWen[4])+"\nA: "+str(optimalBoucWen[5])+"\nx0: "+str(optimalBoucWen[6]))
    end_time = time.time()
    
    execution_time = end_time - start_time

    print("Execution time Bouc-Wen:", execution_time, "seconds")
    bw.plotBoucWen(optimalBoucWen)
