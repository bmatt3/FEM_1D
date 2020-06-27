#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

def f(t):
    return ((t*(1-t))/2)

def main():
   t    = np.linspace(0.00, 1.00, 100)
   x, y = np.loadtxt('./data/Solution.dat', delimiter=',', unpack=True)
   plt.figure()
   plt.plot(x, y,    'b--', label='Approximate Solution')
   plt.plot(t, f(t), 'r--', label='True Solution')
   plt.xlabel('x')
   plt.ylabel('y')
   plt.show()
         
main()
