import matplotlib.pyplot as plt
import numpy as np

amp = 1
width = 60
w2 = width**2
dmax = 3.0*width
dmax2 = dmax**2

sols = []
print dmax2
for d2 in np.linspace(40000, 0, 1000):
    if d2 < dmax2:
        sols.append([d2, amp*np.exp(-d2/w2)])
    else:
        sols.append([d2, 0])

sols_arr = np.array(sols)
plt.plot(sols_arr[:, 0], sols_arr[:, 1])
plt.show()