#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Friday 22 20:14:33 2023

@author: flaviaferrusmarimon
"""

from ctypes import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

## Setting LaTeX parameters for the plots
os.environ['PATH'] = os.environ['PATH'] + ':/Library/TeX/texbin'
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Computer Modern Serif"
})

##########################################
## Loading the function from the c code ##
##########################################

def loading_translators(fun_path = './rkf45_v2.so'):
    # Load the shared library containing the rgk45 function
    rgk45_lib = CDLL(fun_path)

    # Define the argument and return types of the rgk45 function
    rgk45_lib.rkf45.restype = c_double
    rgk45_lib.rkf45.argtypes = [
        POINTER(c_double),  # at: initial time
        POINTER(c_double),  # *x: initial value of x, y
        c_int,      # n: dimension of ode
        POINTER(c_double),  # h: initial step size
        c_int,      # sc: stepsize control
        c_double,  # tol: tolerance
        POINTER(c_double),  # atf: final time
        POINTER(c_double),  # aer: indicator of the integration stepsize evolution
        c_void_p,  # (*ode): pointer to the ode function
    ]
    print(type(rgk45_lib))
    return rgk45_lib 


#############################
## Defining the system ode ##
#############################

def my_ode(t, y, n, f):
    # Calculate the derivative of y
    # y is the input 6-dimensional vector, given by y = (x, y, v11, v12, v21, v22)
    omega = np.sqrt(2)
    epsilon = 0.01
    f[0] = y[1]
    f[1] = -omega**2 * np.sin(y[0]) + epsilon*np.sin(t)
    f[2] = y[4]
    f[3] = y[5]
    f[4] = - omega**2 * np.cos(y[0]) * y[2] 
    f[5] = - omega**2 * np.cos(y[0]) * y[3]
    return f
# Convert the Python function into a C function pointer
ODE_FUNC_TYPE = CFUNCTYPE(None, c_double, POINTER(c_double), c_int, POINTER(c_double))
ode_ptr = ODE_FUNC_TYPE(my_ode)
print(type(ode_ptr))  

def my_ode_2d(t, y, n, f):
    # Calculate the derivative of y
    # y is the input 2-dimensional vector, given by y = (x, y)
    omega = np.sqrt(2)
    epsilon = 0.01
    f[0] = y[1]
    f[1] = -omega**2 * np.sin(y[0]) + epsilon*np.sin(t)
    return f
# Convert the Python function into a C function pointer
#ODE_FUNC_TYPE = CFUNCTYPE(None, c_double, POINTER(c_double), c_int, POINTER(c_double))
ode_ptr_2d = ODE_FUNC_TYPE(my_ode_2d)
#print(type(ode_ptr))   


#######################################################
## Function that applies RKF given initial condition ##
#######################################################

def apply_RK(rgk45_lib, x_0, ode_ptr = ode_ptr, t_0 = 0.0, h_p=0.01, atf_p=0.1, tol_p=1e-6, aer_p =0.0):
    x = np.array(x_0, dtype=np.float64)
    print('Initial conditions:', x)
    n = len(x)
    print('ODE dimension:', n)
    h = c_double(h_p)
    #print('Initial stepsize:',h)
    sc = c_int(0)
    tol = c_double(tol_p)
    aer = c_double(aer_p)

    if h_p > 0:
        Niters = abs(atf_p-t_0)/h_p
        at = c_double(t_0)
        atf = c_double(atf_p)
    else: 
        Niters = -abs(atf_p-t_0)/h_p
        at = c_double(atf_p)
        atf = c_double(t_0) 

    y = np.empty(0)
    temps = np.empty(0)

    for i in range(int(Niters)):
        x_aux = rgk45_lib.rkf45(byref(at), x.ctypes.data_as(POINTER(c_double)), n, byref(h), sc, tol, byref(atf), byref(aer), ode_ptr)     
        ## Here we store the 6 dimensional value of x
        y = np.concatenate([y, x])
        ## Here we store the current time
        temps_aux  = [np.float32(at)]
        temps = np.concatenate([temps, temps_aux])

    print('Point at time:',x, at)
    return temps, y


def Newton_method(rgk45_lib, ode_ptr = ode_ptr, x_0=[0,0.01], t_0 = 0, v_0 = [1, 0, 0, 1], Nmax= 100, tol= 1e-4):
    
    niter = 0
    
    """
    xx = np.concatenate((x_0, v_0), axis = None)
    t_1, xx_1 = apply_RK(rgk45_lib, ode_ptr=ode_ptr, x_0= xx, t_0 = t_0, h_p=0.01, atf_p= t_0 + 0.1, tol_p=1e-6, aer_p =0.0)   
    x_1 = xx_1[::6]
    y_1 = xx_1[1::6]
    v_11 = xx_1[2::6]
    v_12 = xx_1[3::6]
    v_21 = xx_1[4::6]
    v_22 = xx_1[5::6] 
    tol_x = np.linalg.norm(x_0 - np.array([x_1[-1], y_1[-1] ]))
    x_0 = np.array([x_1[-1], y_1[-1]] )
    """
    tol_x = 1
       
    while niter<Nmax and tol_x>tol:
        print('----')
        print('Iteration:', niter)
        print('Input x_0:', x_0)        
        x_aux, t_aux = Newt_step(rgk45_lib=rgk45_lib, ode_ptr = ode_ptr, x_0 = x_0, t_0 = t_0, v_0= v_0)
        tol_x = np.linalg.norm(x_aux - x_0)
        x_0, t_0 = x_aux, t_aux
        print('Tolerance:', tol_x)
        niter += 1
        
    if niter < Nmax:
        print('Solution found!')
        #print(t_0, x_0, v_0)
        
    return x_0, t_0
     
    
        
def Newt_step(rgk45_lib, ode_ptr= ode_ptr, x_0=[0,0.01], t_0 = 0, v_0 = [1, 0, 0, 1]):   
    xx = np.concatenate((x_0, v_0), axis = None)   
    t_1, xx_1 = apply_RK(rgk45_lib, ode_ptr=ode_ptr, x_0= xx, t_0 = t_0, h_p=0.01, atf_p= t_0 + 0.1, tol_p=1e-6, aer_p =0.0)       
    x_1 = xx_1[::6]
    y_1 = xx_1[1::6]
    v_11 = xx_1[2::6]
    v_12 = xx_1[3::6]
    v_21 = xx_1[4::6]
    v_22 = xx_1[5::6]          
    v_matrix = np.array([[v_11[-1], v_12[-1]], [v_21[-1], v_22[-1]]])
    I = np.identity(2)    
    x_k1 = np.array([x_1[-1], y_1[-1]]) - np.matmul( (np.array([x_1[-1], y_1[-1]]) - x_0) , np.linalg.inv(v_matrix - I))
    return x_k1, t_1[-1]
    
    
def generate_orbits(rgk45_lib, x_0, t_0 = 0.0, atf_p=2*np.pi, NP = 20, thresh = 1e-4):
    """This function applies the map of 2 dimensions to the obtained 
    initial condition

    Args:
        rgk45_lib (_type_): _description_
        x_0 (_type_): _description_
        t_0 (float, optional): _description_. Defaults to 0.0.
        atf_p (_type_, optional): _description_. Defaults to 2*np.pi.
        NP (int, optional): _description_. Defaults to 20.
        thresh (_type_, optional): _description_. Defaults to 1e-4.

    Returns:
        _type_: _description_
    """
    xaux, yaux = np.empty(0), np.empty(0)   
    x_, y_ = np.empty(0), np.empty(0) 
    iter = 0
    
    while iter < NP: 
        xaux = apply_RK(rgk45_lib, ode_ptr = ode_ptr_2d, x_0= x_0, t_0 = t_0, h_p=0.01, atf_p=atf_p, tol_p=1e-6, aer_p =0.0)
        print(xaux)
        x_ = np.concatenate([x_, xaux[0]])
        y_ = np.concatenate([y_, xaux[1]])
    return x_, y_

def plot_orbit(x, y, picture_path = './Ex3_py_1.png'):
    plt.figure(figsize=(10,8))
    plt.scatter(x, y, s = 0.4)
    plt.xlabel('$x$', fontsize=15)
    plt.ylabel('$y$', fontsize=15)
    plt.title('Periodic orbit', fontsize=17.5)
    
    plt.xlim([min(x) - 0.001, max(x) + 0.001])
    plt.ylim([min(y) - 0.001, max(y) + 0.001])
    
    plt.savefig(picture_path)
    plt.show() 
    return 0

def plot_orbit_2(t, x, picture_path = './Ex3_py_1.png'):
    plt.figure(figsize=(10,8))
    plt.scatter(t, x, s = 0.4)
    plt.xlabel('$t$', fontsize=15)
    plt.ylabel('$x$', fontsize=15)
    plt.title('Periodic orbit', fontsize=17.5)
    
    plt.xlim([min(t) - 0.01, max(t) + 0.01])
    plt.ylim([min(x) - 0.01, max(x) + 0.01])
    
    plt.savefig(picture_path)
    plt.show() 
    return 0  

### MAIN
rgk45_lib = loading_translators(fun_path = './rkf45_v2.so')
x_0 = [0.01,0.0]
t_0 = 0
v_0 = [1,0,0,1]
NP= 25

print('---------')
print('Trying one step:')
# ode_ptr= ode_ptr(w = np.sqrt(2)),
x_1, t_1 = Newt_step(rgk45_lib=rgk45_lib, x_0 = x_0, t_0 = t_0, v_0= v_0)
print(x_1, t_1)
print('Tolerance first step:', np.linalg.norm( x_1[0] - x_0[0]))
print('---------')
print('Trying 15 steps:')
x_2, t_2 = Newton_method(rgk45_lib, x_0=x_0, t_0 = t_0, v_0 = v_0, Nmax= NP, tol= 1e-4)
print(x_2, t_2)
xx_1 = x_2[::6]
yy_1 = x_2[1::6]
print(len(xx_1))

#xx_, yy_ = generate_orbits(rgk45_lib, np.array([xx_1, yy_1]), t_0 = 0.0, atf_p=2*np.pi, NP = 20, thresh = 1e-4)
#plot_orbit(xx_, yy_, picture_path = './Ex1.png')


