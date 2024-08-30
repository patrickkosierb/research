

def cost_function(params, f1, f_r):
    c_0a, c_0b, c_1a, c_1b, alpha_a, alpha_b, k_0, k_1, A, beta_gamma = params
    squared_errors = (f1[:-1] - f_r(c_0a, c_0b, c_1a, c_1b, alpha_a, alpha_b, k_0, k_1, A, beta_gamma))**2
    return np.sum(squared_errors)

def f_r(c_0a, c_0b, c_1a, c_1b, alpha_a, alpha_b, k_0, k_1, A, beta_gamma):
    beta = beta_gamma
    gamma = beta_gamma
    n = 2

    z = 0
    y = 0

    f_r = []

    for i in range(len(t) - 1):
        dt_i = (t[i + 1] - t[i])
        x_i = x[i]
        dxdt_i = dxdts[i]
        volt_i = volt[i]
        alpha = alpha_a + alpha_b * volt_i
        c_0 = c_0a + c_0b * volt_i
        c_1 = c_1a + c_1b * volt_i

        if c_0 + c_1 != 0:
            term1 = alpha * z
            #term2 = k_0 * (x_i - y)
            if not (np.isinf(k_0) or np.isinf(x_i) or np.isinf(y)):
                term2 = k_0 * (x_i - y)
            else:
                term2 = 0  # Handle overflow
            term3 = c_0 * dxdt_i

            # Check for overflow
            if not (np.isinf(term1) or np.isinf(term2) or np.isinf(term3)):
                dydt = (term1 + term2 + term3) / (c_0 + c_1)
            else:
                dydt = 0  # Handle overflow
        else:
            dydt = 0  # Handle division by zero
        y += dydt * dt_i


        if np.abs(z) != 0:  # Check for division by zero or invalid values
            term1 = -gamma * np.abs(dxdt_i - dydt) * z
            term2 = -beta * (dxdt_i - dydt)
            term3 = A * (dxdt_i - dydt)

            # Check for overflow
            if not (np.isinf(term1) or np.isinf(term2) or np.isinf(term3)):
                dzdt = term1 * np.power(np.abs(z), n - 1) - term2 * np.power(np.abs(z), n) + term3
            else:
                dzdt = 0  # Handle overflow
        else:
            dzdt = 0  # Handle division by zero or invalid values

        z += dzdt * dt_i

        #f_r.append(c_1 * dydt + k_1 * (x_i))
        if not (np.isinf(c_1) or np.isinf(dydt) or np.isinf(k_1)):
            f_r.append(c_1 * dydt + k_1 * (x_i))
        else:
            f_r.append(0)  # Handle overflow

    return f_r

from scipy.optimize import leastsq
alpha_a_initial = 1921.141 #N/m
alpha_b_initial = 5882.51 #N/Vm
c_0a_initial = 651.4718 #Ns/m
c_0b_initial = 1043.7559 #Ns/Vm
c_1a_initial = 2089.263 #Ns/m
c_1b_initial = 14384.918 #Ns/Vm
k_0_initial = 1940.405 #N/m
k_1_initial = 1.751268 #N/m
A_initial = 155.32 # /m
beta_gamma_initial = 36332.07 # /m^2
n = 2 #

# Define a residual function that computes the residuals (differences) between f1 and f_r
def residual(params, f1, f_r):
    c_0a, c_0b, c_1a, c_1b, alpha_a, alpha_b, k_0, k_1, A, beta_gamma = params
    return f[:-1] - f_r(c_0a, c_0b, c_1a, c_1b, alpha_a, alpha_b, k_0, k_1, A, beta_gamma)

# Initial guesses for the parameters
initial_guess = [c_0a_initial, c_0b_initial, c_1a_initial, c_1b_initial, alpha_a_initial, alpha_b_initial, k_0_initial, k_1_initial, A_initial, beta_gamma_initial]

# Perform the least squares optimization to find the best parameters
result = leastsq(residual, initial_guess, args=(f, f_r))

# Extract the optimal parameters
optimal_params = result[0]
print(optimal_params)