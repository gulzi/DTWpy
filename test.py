from pandas.io.data import DataReader
from datetime import datetime
from distance_measure import distancefunc
import dtwpy
import plotDTWalignment

f = DataReader("F", "yahoo", datetime(2000, 1, 1), datetime(2012, 1, 1))
f_2008 = f[f.index.year == 2008]
f_2009 = f[f.index.year == 2009]
x = f_2008.Volume.values
y = f_2009.Volume.values    

dist = distancefunc(name="manhattan")

window,distance1, path = dtwpy.dtw(x, y, dist=dist,windowtype="paliwal", windowsize=50,  pattern="symmetricP1", normalized=False,dist_only=False,cost=False)
#print(distance1)
#window,d,path = dtwpy.fastdtw(x, y, radius=10, dist=dist, pattern="symmetricP1", normalized=False)
#window,distance1, _ = dtwpy.dtw(x, y, dist=dist,windowtype="itakura", windowsize=50,  pattern="symmetricP1", normalized=False,dist_only=False)
plotDTWalignment.plotalignment(window,path, title="")
'''
lenght_of_time_series = [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000]
lenght_of_time_series_for_dtw = [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000,5000,6000,8500,10000]
time_taken_by_r = [0.4, .68, 1.28, 1.88, 2.8, 3.68, 5.81, 6.62, 8.12, 10.28, 12.56, 14.78, 17.3, 19.44, 22.58, 25.64, 29.2, 32.86, 36.76, 45.16]
time_taken_by_r_for_dtw = [0.4, .68, 1.28, 1.88, 2.8, 3.68, 5.81, 6.62,10.28,19.44,29.2,45.16]
time_taken_by_dtw = [.897, 3.5, 7.9, 14.0, 22.3, 34.76, 47.265, 61.79,88.778,127.743,182.067,495.006,720.476]
time_taken_by_fastdtw_radius_6 = [0.283, 0.586, 0.877, 1.194, 1.494, 1.812, 2.125, 2.463, 2.754, 3.073, 3.386, 3.725, 4.046, 4.404, 4.701, 5.052, 5.382, 5.693, 6.12, 6.393]
time_taken_by_fastdtw_radius_10 = [0.47, 0.98, 1.55, 2.11, 2.66, 3.21, 3.787, 4.348, 4.923, 5.454, 6.040, 6.626, 7.206, 7.755, 8.349, 8.905, 9.619, 10.061, 10.657, 11.282]
#distance1, path = dtwpy.fastdtw(x, y, dist=dist, radius=6, pattern="symmetricP0", normalized=False)

plt.figure(num=None, figsize=(5, 5), dpi=80, facecolor='w', edgecolor='k')
temp1x = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
temp1y = [0,1,2,3,4,5,5,7,8,9,10,12,12,13,14,15,16]
temp2x = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
temp2y = [0,1,1,3,4,5,6,7,7,8,9,9,10,11,12,13,15]
plt.plot(temp1x,temp1y,linewidth=2.0,label="correct",linestyle='--',color="black")
plt.plot(temp2x,temp2y,linewidth=2.0,label="incorrect",color="black")
plt.title("DTW boundary constraint")
plt.legend(loc=2)
plt.show()

plt.plot(lenght_of_time_series,time_taken_by_r,linewidth=2.0,label="R package's DTW")
plt.plot(lenght_of_time_series,time_taken_by_fastdtw_radius_6,linewidth=2.0,label="DTWpy's FastDTW with radius = 6")
plt.plot(lenght_of_time_series,time_taken_by_fastdtw_radius_10,linewidth=2.0,label="DTWpy's FastDTW with radius = 10")
plt.axis([500, 10000, 0, 50])
plt.grid(True)
plt.xlabel('lenght of time series')
plt.ylabel('Time taken')
plt.legend()
plt.title("DTWPy vs R Package for DTW")
plt.show()

plt.plot(lenght_of_time_series_for_dtw,time_taken_by_r_for_dtw,linewidth=2.0,label="R package's DTW")
plt.plot(lenght_of_time_series_for_dtw,time_taken_by_dtw,linewidth=2.0,label="DTWPy's DTW")
plt.axis([500,4000,0,65])
plt.grid(True)
plt.xlabel('lenght of time series')
plt.ylabel('Time taken')
plt.legend()
plt.title("DTWPy vs R Package for DTW")
plt.show()



# DTWPy code for testing in terminal

from distance_measure import distancefunc
import dtwpy
import cProfile, pstats, io
import numpy as np 
dist = distancefunc(name="manhattan")

idx=np.linspace(0,6.28,10000)
query= np.sin(idx)+np.random.uniform((10000,))/10
template = np.cos(idx)
pr = cProfile.Profile()
pr.enable()
distance1, path = dtwpy.dtw(query, template, dist=dist,windowtype="", windowsize="",  pattern="symmetric2", normalized=False,dist_only=False,cost=False)
pr.disable()
s = io.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
print(s.getvalue())



#mlpy testing
import cProfile, pstats, StringIO
import mlpy
import numpy as np
idx=np.linspace(0,6.28,5000)
query= np.sin(idx)+np.random.uniform((2000,))/10
template = np.cos(idx)
pr = cProfile.Profile()
pr.enable()
dist, cost, path = mlpy.dtw_std(query, template, dist_only=False)
pr.disable()
s = StringIO.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
print s.getvalue()



#R code for testing in terminal
idx<-seq(0,6.28,len=7000);
query<-sin(idx)+runif(7000)/10;

## A cosine is for template; sin and cos are offset by 25 samples
template<-cos(idx)

## Find the best match with the canonical recursion formula
library(dtw);
Rprof("out.out")
alignment<-dtw(query,template,keep=TRUE,step=symmetric1,distance.only=FALSE);
Rprof(NULL)
summaryRprof("out.out")
'''
