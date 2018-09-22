import numpy as np

# Helper Functions
def MatchArrayLengths(*Args):
    MaxLength = 1
    for Arg in Args:
        try:
            if (MaxLength > 1) and MaxLength != len(Arg):
                raise ValueError()
            elif MaxLength < len(Arg):
                MaxLength = len(Arg)
        except:
            pass
    ReturnValues = []
    for Arg in Args:
        ReturnValues.append(Arg*np.ones(MaxLength))
    return ReturnValues

# One to One binding Functions
def coth(z):
    return np.reciprocal(np.tanh(z))

def SS_AB(Kd, A0, B0, X_Y_Axes = False, Normalized = False):
    """
    One to One : A-B
    Kd, A0, B0
    X_Y_Axes : tuple of three elements from Kd, A0, B0 as str
        returns 2D array with first element on the "x-axis", second on the "y-axis" and third constant
    
    """
    Func = lambda Kd, A0, B0: (A0 + B0 + Kd - np.sqrt(np.power(A0 + B0 + Kd, 2.0) - 4.0 * A0 * B0))/2.0
    
    if not X_Y_Axes:
        Kd, A0, B0 = MatchArrayLengths(Kd, A0, B0)
        if Normalized:
            return Func(Kd, A0, B0)/np.min([A0, B0])
        else:
            return Func(Kd, A0, B0)
    else:
        x, y = np.meshgrid(vars()[X_Y_Axes[0]],vars()[X_Y_Axes[1]])
        Result = Func(Kd, x, y)
        if Normalized:
            Norm = x
            Norm[x>y] = y[x>y]
            Result = Result/Norm
        return Result
        
def SS_AB_2(Kd, A0, B0):
    """
    One to One : A-B
    Kd, A0, B0
    """
    Kd, A0, B0 = MatchArrayLengths(Kd, A0, B0)
    z = np.sqrt(np.power(A0-B0,2.0)+np.power(Kd,2)+2.0*(A0+B0)*Kd)
    return (2.0*A0*B0)/(A0+B0+Kd+z)

def SS_A_AB(Kd, A0, B0, X_Y_Axes = False, Normalized = False):
    """
    One to One : A-B
    Kd, A0, B0
    """
    # Kd, A0, B0 = MatchArrayLengths(Kd, A0, B0)
    Func = lambda Kd, A0, B0: A0 - SS_AB(Kd, A0, B0, X_Y_Axes = False, Normalized = False)

    if not X_Y_Axes:
        Kd, A0, B0 = MatchArrayLengths(Kd, A0, B0)
        if Normalized:
            return Func(Kd, A0, B0)/A0
        else:
            return Func(Kd, A0, B0)
    else:
        x, y = np.meshgrid(vars()[X_Y_Axes[0]],vars()[X_Y_Axes[1]])
        Result = Func(Kd, x, y)
        if Normalized:
            if X_Y_Axes[0] == 'A0':
                Norm = x
            else:
                Norm = y
            Result = Result/Norm
        return Result

def SS_B_AB(Kd, A0, B0, X_Y_Axes = False, Normalized = False):
    """
    One to One : A-B
    Kd, A0, B0
    """
    # Kd, A0, B0 = MatchArrayLengths(Kd, A0, B0)
    Func = lambda Kd, A0, B0: B0 - SS_AB(Kd, A0, B0, X_Y_Axes = False, Normalized = False)

    if not X_Y_Axes:
        Kd, A0, B0 = MatchArrayLengths(Kd, A0, B0)
        if Normalized:
            return Func(Kd, A0, B0)/B0
        else:
            return Func(Kd, A0, B0)
    else:
        x, y = np.meshgrid(vars()[X_Y_Axes[0]],vars()[X_Y_Axes[1]])
        Result = Func(Kd, x, y)
        if Normalized:
            if X_Y_Axes[0] == 'B0':
                Norm = x
            else:
                Norm = y
            Result = Result/Norm
        return Result

#def SS_B_AB(Kd, A0, B0, X_Y_Axes = False, Normalized = False):
#    """
#    One to One : A-B
#    Kd, A0, B0
#    """
#    Kd, A0, B0 = MatchArrayLengths(Kd, A0, B0)
#    return B0 - SS_AB(Kd, A0, B0)

def AB(Kd, A0, B0, kon, t):
    """
    One to One : A-B
    Kd, A0, B0, kon, t
    """
    t = np.array(t)
    if np.sum(t==0) > 0:
        t_contains_zeros = True
        t_zeroes = t==0
        t[t==0] = 1
    else:
        t_contains_zeros = False
    Kd, A0, B0, kon, t = MatchArrayLengths(Kd, A0, B0, kon, t)
#     z = np.tanh(0.5*koff*np.sqrt(np.power(A0-B0,2.0)+np.power(Kd,2)+2.0*(A0+B0)*Kd)*t)
#     return SS_AB(Kd, A0, B0)*z
    z = np.sqrt(np.power(A0-B0,2.0)+np.power(Kd,2)+2.0*(A0+B0)*Kd)
    Result = (2.0*A0*B0)/(A0+B0+Kd+z*coth(0.5*kon*z*t))
    if t_contains_zeros:
        Result[t_zeroes] = 0
    return Result

def A_AB(Kd, A0, B0, kon, t):
    """
    One to One : A-B
    Kd, A0, B0, kon, t
    """
    t = np.array(t)
#     if np.sum(t==0) > 0:
#         t_contains_zeros = True
#         t_zeroes = t==0
#         t[t==0] = 1
#     else:
#         t_contains_zeros = False
    Kd, A0, B0, kon, t = MatchArrayLengths(Kd, A0, B0, kon, t)
#     z = np.tanh(0.5*koff*np.sqrt(np.power(A0-B0,2.0)+np.power(Kd,2)+2.0*(A0+B0)*Kd)*t)
#     return SS_AB(Kd, A0, B0)*z
    Result = A0-AB(Kd, A0, B0, kon, t)
#     if t_contains_zeros:
#         Result[t_zeroes] = 0
    return Result


def B(Kd, A0, B0, kon, t):
    """
    One to One : A-B
    Kd, A0, B0, kon, t
    """
    t = np.array(t)
#     if np.sum(t==0) > 0:
#         t_contains_zeros = True
#         t_zeroes = t==0
#         t[t==0] = 1
#     else:
#         t_contains_zeros = False
    Kd, A0, B0, kon, t = MatchArrayLengths(Kd, A0, B0, kon, t)
#     z = np.tanh(0.5*koff*np.sqrt(np.power(A0-B0,2.0)+np.power(Kd,2)+2.0*(A0+B0)*Kd)*t)
#     return SS_AB(Kd, A0, B0)*z
    Result = B0-AB(Kd, A0, B0, kon, t)
#     if t_contains_zeros:
#         Result[t_zeroes] = 0
    return Result

def SS_AA(Kd, A0):
    """
    One to One : A-A
    Kd, A0
    """
    Kd, A0 = MatchArrayLengths(Kd, A0)
    return (4.0*A0 + Kd - np.sqrt(np.power(Kd, 2.0) + 8.0 * A0 * Kd))/8.0

def SS_A(Kd, A0):
    """
    One to One : A-A
    Kd, A0
    """
    Kd, A0 = MatchArrayLengths(Kd, A0)
    return A0-2.0*SS_AA(Kd, A0)

# def A(Kd, A0, kon, t):
#     t = np.array(t)
#     if np.sum(t==0) > 0:
#         t_contains_zeros = True
#         t_zeroes = t==0
#         t[t==0] = 1
#     else:
#         t_contains_zeros = False
#     koff = Kd*kon
#     Kd, A0, kon, t, koff = MatchArrayLengths(Kd, A0, kon, t, koff)
#     z = np.sqrt(koff*(8.0*A0*kon+koff))
#     Result = np.real((-koff+z*np.tanh(0.5*(t*z+(6.0+32.0*A0*kon/koff)*np.arctanh((4.0*A0*kon+koff)/z+0j))))/(4.0*kon))
#     if t_contains_zeros:
#         Result[t_zeroes] = A0[0]
#     return Result

# def AA(Kd, A0, kon, t):
#     t = np.array(t)
#     if np.sum(t==0) > 0:
#         t_contains_zeros = True
#         t_zeroes = t==0
#         t[t==0] = 1
#     else:
#         t_contains_zeros = False
#     Result = (A0-A(Kd, A0, kon, t))/2.0
#     if t_contains_zeros:
#         Result[t_zeroes] = 0
#     return Result

'''
def AA(Kd, A0, kon, t):
    """
    One to One : A-A
    Kd, A0, kon,t
    """
#     t = np.array(t)
#     if np.sum(t==0) > 0:
#         t_contains_zeros = True
#         t_zeroes = t==0
#         t[t==0] = 1
#     else:
#         t_contains_zeros = False
    koff = Kd*kon
    Kd, A0, kon, t, koff = MatchArrayLengths(Kd, A0, kon, t, koff)
    z = np.sqrt(koff*(koff+8.0*A0*kon))
    Result = np.real((koff+4*A0*kon-z*np.tanh(0.5*(z*t+2*np.arctanh((koff+0j+4.0*A0*kon)/z))))/(8.0*kon))
#     if t_contains_zeros:
#         Result[t_zeroes] = A0[0]
    return Result
'''
    
def AA(Kd, A0, kon, t):
    """
    One to One : A-A
    Kd, A0, kon,t
    """
#     t = np.array(t)
#     if np.sum(t==0) > 0:
#         t_contains_zeros = True
#         t_zeroes = t==0
#         t[t==0] = 1
#     else:
#         t_contains_zeros = False
    koff = Kd*kon
    Kd, A0, kon, t, koff = MatchArrayLengths(Kd, A0, kon, t, koff)
    z = np.sqrt(Kd*(Kd+8.0*A0))
    Result = np.real((Kd+4*A0-z*np.tanh(0.5*z*t*kon+np.arctanh((Kd+0j+4.0*A0)/z)))/(8.0))
#     if t_contains_zeros:
#         Result[t_zeroes] = A0[0]
    return Result

def A(Kd, A0, kon, t):
    """
    One to One : A-A
    Kd, A0, kon, t
    """
#     t = np.array(t)
#     if np.sum(t==0) > 0:
#         t_contains_zeros = True
#         t_zeroes = t==0
#         t[t==0] = 1
#     else:
#         t_contains_zeros = False
    Result = (A0-2.0*AA(Kd, A0, kon, t))
#     if t_contains_zeros:
#         Result[t_zeroes] = 0
    return Result