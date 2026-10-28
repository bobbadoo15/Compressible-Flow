""" 
This finds beta in the beta-theta-M Relattion (if theta and M are known). 
    
GREEK UNICODE REFERENCE SHEET
Usage: print("\u03B1") will print α

--- LOWERCASE LETTERS ---
\u03B1 : α (alpha)      \u03B7 : η (eta)        \u03BD : ν (nu)         \u03C4 : τ (tau)
\u03B2 : β (beta)       \u03B8 : θ (theta)      \u03BE : ξ (xi)         \u03C5 : υ (upsilon)
\u03B3 : γ (gamma)      \u03B9 : ι (iota)       \u03BF : ο (omicron)    \u03C6 : φ (phi)
\u03B4 : δ (delta)      \u03BA : κ (kappa)      \u03C0 : π (pi)         \u03C7 : χ (chi)
\u03B5 : ε (epsilon)    \u03BB : λ (lambda)     \u03C1 : ρ (rho)        \u03C8 : ψ (psi)
\u03B6 : ζ (zeta)       \u03BC : μ (mu)         \u03C3 : σ (sigma)      \u03C9 : ω (omega)

--- UPPERCASE LETTERS ---
\u0391 : Α (Alpha)      \u0397 : Η (Eta)        \u039D : Ν (Nu)         \u03A4 : Τ (Tau)
\u0392 : Β (Beta)       \u0398 : Θ (Theta)      \u039E : Ξ (Xi)         \u03A5 : Υ (Upsilon)
\u0393 : Γ (Gamma)      \u0399 : Ι (Iota)       \u039F : Ο (Omicron)    \u03A6 : Φ (Phi)
\u0394 : Δ (Delta)      \u039A : Κ (Kappa)      \u03A0 : Π (Pi)         \u03A7 : Χ (Chi)
\u0395 : Ε (Epsilon)    \u039B : Λ (Lambda)     \u03A1 : Ρ (Rho)        \u03A8 : Ψ (Psi)
\u0396 : Ζ (Zeta)       \u039C : Μ (Mu)         \u03A3 : Σ (Sigma)      \u03A9 : Ω (Omega)
"""

import numpy as np
import math as m
import sys

def main():
    """ This runs the program, importing constants. """
    # Read greek info
    greek = Greek()
    # Find lambda
    l = lamb(gamma,M,theta)
    # Find chi
    c = chi(gamma,l,M,theta)
    # Find beta
    b = beta(c,gamma,l,M,theta)
    # Display the results
    results(b,greek)

def beta(c,gamma,l,M,theta):
    """ This finds beta using chi and lambda; also depends on weak/stron shock. """
    delta = float(input("Enter the desired delta value (0 = strong shock, 1 = weak shock): "))
    # Check for quit command first
    if not isinstance(delta, (int,float)):
        if delta in ["q", "quit"]:
            print("\nExiting program. Goodbye!\n")
            sys.exit()
        else:
            print(f"\nInvalid input. Please enter a number or quit program.")
    return np.rad2deg(m.atan((M**2 - 1 + (2*l*m.cos(((4*pi*delta)+m.acos(c))/3))) / (3*(1 + 0.5*(gamma-1)*M**2)*m.tan(theta))))

def chi(gamma,l,M,theta):
    """ This finds the chi value. """
    return l**(-3) * ((M**2 - 1)**3 - 9*(1 + 0.5*(gamma-1)*M**2)*(1 + 0.5*(gamma-1)*M**2 + 0.25*(gamma+1)*M**4)*(m.tan(theta))**2)

def Greek():
    """ This stores greek unicode values. """
    greek = {
    'alpha'  : {'lower': '\u03B1', 'upper': '\u0391'},
    'beta'   : {'lower': '\u03B2', 'upper': '\u0392'},
    'gamma'  : {'lower': '\u03B3', 'upper': '\u0393'},
    'delta'  : {'lower': '\u03B4', 'upper': '\u0394'},
    'epsilon': {'lower': '\u03B5', 'upper': '\u0395'},
    'zeta'   : {'lower': '\u03B6', 'upper': '\u0396'},
    'eta'    : {'lower': '\u03B7', 'upper': '\u0397'},
    'theta'  : {'lower': '\u03B8', 'upper': '\u0398'},
    'iota'   : {'lower': '\u03B9', 'upper': '\u0399'},
    'kappa'  : {'lower': '\u03BA', 'upper': '\u039A'},
    'lambda' : {'lower': '\u03BB', 'upper': '\u039B'},
    'mu'     : {'lower': '\u03BC', 'upper': '\u039C'},
    'nu'     : {'lower': '\u03BD', 'upper': '\u039D'},
    'xi'     : {'lower': '\u03BE', 'upper': '\u039E'},
    'omicron': {'lower': '\u03BF', 'upper': '\u039F'},
    'pi'     : {'lower': '\u03C0', 'upper': '\u03A0'},
    'rho'    : {'lower': '\u03C1', 'upper': '\u03A1'},
    'sigma'  : {'lower': '\u03C3', 'upper': '\u03A3'},  # Word-final lowercase variant is '\u03C2' (ς)
    'tau'    : {'lower': '\u03C4', 'upper': '\u03A4'},
    'upsilon': {'lower': '\u03C5', 'upper': '\u03A5'},
    'phi'    : {'lower': '\u03C6', 'upper': '\u03A6'},
    'chi'    : {'lower': '\u03C7', 'upper': '\u03A7'},
    'psi'    : {'lower': '\u03C8', 'upper': '\u03A8'},
    'omega'  : {'lower': '\u03C9', 'upper': '\u03A9'}
}
    return greek

def lamb(gamma,M,theta):
    """ This finds the lambda value. """
    return m.sqrt((M**2 - 1)**2 - 3*(1 + 0.5*(gamma-1)*M**2)*(1 + 0.5*(gamma+1)*M**2)*(m.tan(theta))**2)

def results(b,greek):
    """ This displays the results. """
    print()
    print("============================================")
    print(f"   {greek["beta"]["lower"]} in the {greek["beta"]["lower"]}-{greek["theta"]["lower"]}-M Relationship in Degrees")
    print("============================================")
    print(f"{greek["beta"]["lower"]} = {b:10.3f}")
    print("============================================")
    print()

if __name__ == "__main__":
    # Constants
    gamma = 1.4
    M = 5.443
    pi = np.pi
    theta = np.deg2rad(15)
    # Run Program
    main()