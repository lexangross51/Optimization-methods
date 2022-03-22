from math import *
import matplotlib.pyplot as plt
from scipy import interpolate
import numpy as np

dichotomy = [18, 26, 32, 38, 46, 52, 58]
golden_ratio = [14, 18, 23, 28, 33, 38, 42]
fibonacci = [14, 18, 23, 28, 33, 37, 42]
eps = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]
log_eps = [log10(x) for x in eps]

xnew = np.linspace(log_eps[0], log_eps[-1], 30)
f1 = interpolate.interp1d(log_eps, dichotomy, kind = 'cubic')
f2 = interpolate.interp1d(log_eps, golden_ratio, kind = 'cubic')
f3 = interpolate.interp1d(log_eps, fibonacci, kind = 'cubic')
plt.plot(log_eps, dichotomy, 'o', color='red')
plt.plot(xnew, f1(xnew), color='red', label = 'dichotomy')
plt.plot(log_eps, golden_ratio, 'o', color='green')
plt.plot(xnew, f2(xnew), color='green', label = 'golden ratio')
plt.plot(log_eps, fibonacci, 'o', color='blue')
plt.plot(xnew, f3(xnew), color='blue', label = 'fibonacci')
plt.legend()
plt.grid()
plt.show()