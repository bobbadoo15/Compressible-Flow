import numpy as np

def main(gamma,imax,M2,tol,vM2):
    # Call Newton-raphson function
    e,iter,M2 = NewtonRaphson(gamma,imax,M2,tol,vM2)
    # Display the results
    print()
    print("============================================")
    print(f"    Prandtl-Meyer Expansion in Region 2")
    print("============================================")
    print(f"Mach Speed   = {M2:10.3f}")
    print(f"Iterations   = {iter:6d}")
    print(f"Rel. % Error = {e:14.3e}")
    print("============================================")
    print()

def error(M2_new,M2):
    """ This calculates the error for convergence. """
    # Approximate percent relative error
    e = np.abs((M2_new-M2)/M2_new)*100
    return e

def f(gamma,M2,vM2):
    """ 
    This calculates the Prandtl-Meyer equation. The actual value is subtracted 
    so the f(x) = 0 because v(M2) = 56.3798. So, f(M2) = 0 = v(M2_new) - v(M2). 
    """
    v = np.sqrt((gamma+1)/(gamma-1))*np.arctan(np.sqrt(((gamma-1)/(gamma+1))*(M2**2-1))) \
        -np.arctan(np.sqrt(M2**2-1))
    return v - vM2

def df(gamma,M2):
    """ This finds the first derivative. """
    return (np.sqrt(M2**2-1)/(M2*(1+((gamma-1)/2)*M2**2)))

def NewtonRaphson(gamma,imax,M2,tol,vM2):
    """ Runs the Newton-Raphson iteration. """
    iter = 0
    while iter<=imax:
        # Call f(x)
        fx = f(gamma,M2,vM2)
        # Call f'(x)
        dfx = df(gamma,M2)
        # Check derivative
        if abs(dfx) < 1e-12:
            print("Error: Derivative is approximately zero.")
            break
        # Newton-Raphson equation
        M2_new = M2 - fx/dfx
        # Find error amount
        e = error(M2_new,M2)
        # Redefine M2
        M2 = M2_new
        # Check convergence
        if e < tol:
            break
        # Count the iteration
        iter += 1

    return e,iter,M2

if __name__ == "__main__":
    # Initial guess
    M2 = 3
    # Tolerence
    tol = 1e-6
    # Max interations
    imax = 1000
    # Variables/constants
    gamma = 1.4
    vM2 = np.deg2rad(56.3798)
    # Run Program
    main(gamma,imax,M2,tol,vM2)