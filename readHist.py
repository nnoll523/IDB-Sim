import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm

loci = 10000
runs = 200
data = np.zeros((runs,loci))

i = 0
for line in open('out.txt', 'rb'):
    hist = line.split()
    data[i,int(hist[0])-1] = int(hist[2])
    if (int(hist[0]) == loci):
        i+=1

for i in range(15):
    print sum(data[i,:])
    plt.plot(range(loci),data[i,:], 'o', color=cm.hsv(i/15.,1), label='Generation %i'%(5*i))
plt.xlabel("Block Width")
plt.ylabel("Count")
plt.legend()
plt.show()
