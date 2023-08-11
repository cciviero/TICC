import matplotlib.pyplot as plt
import numpy as np

"""
Histogram:
#events - depth
all the required inputs are in "INPUT"!
"""

# Change the style of the plots
plt.style.use('ggplot')

# --------- INPUT
user_bins = np.arange(0, 600 + 1, 1)
Title = '1887 events, 1km bins'
sel_ev_fi = "./event_info/pdata_events.txt"
# --------- END INPUT

evs = np.loadtxt(sel_ev_fi, dtype='object', delimiter=',')

# ================== INVERTED Depths
ev_depths_inv = evs[:, 3].astype(float)
print "\n\n%s events were found!" % len(ev_depths_inv)

hist_inv, bins = np.histogram(ev_depths_inv, bins=user_bins)
width = 0.7 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2

plt.bar(center-width/2, hist_inv, align='center', width=width/2,
        label='inverted', edgecolor='none')

# ================== Catalog Depths
ev_depths_cat = evs[:, 11].astype(float)

hist_cat, bins = np.histogram(ev_depths_cat, bins=user_bins)
width = 0.7 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2

plt.bar(center, hist_cat, color='r', align='center', width=width/2,
        label='catalog', edgecolor='none')

# ================== PLOTTING
plt.xlabel('Depth', size=32, weight='bold')
plt.ylabel('#Events', size=32, weight='bold')
plt.xticks(size=24, weight='bold')
plt.yticks(size=24, weight='bold')
plt.legend(prop={'size': 24, 'weight': 'bold'})

plt.xlim(xmin=-5)

# Horizontal line, 10 events!
plt.axhline(10, 0, 1, c='k', ls='--', lw=2)
plt.annotate('10 events', fontsize=24, xy=(300, 10), xytext=(350, 100),
             arrowprops=dict(facecolor='black', shrink=0.05))

plt.title(Title, size=32, weight='bold')

print "=================="
for bi in range(len(hist_inv)):
    print "%s -- %s -- %s" % (bins[bi], hist_inv[bi], hist_cat[bi])

plt.show()
