import matplotlib.pyplot as plt
import numpy as np

def plotwithQT(x, y, path,title):
    xaxis = [x1[0] for x1 in path]
    yaxis = [y1[1] for y1 in path]
    
    plt.figure(0)
    # plt.axes([0, x.size, 0, y.size])
    ax1 = plt.subplot2grid((5, 5), (0, 0), rowspan=4)
    x1 = np.arange(0, y.size)    
    ax1.set_ylabel('Template axis')            
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.plot(y, x1, color='green')
    
    ax1.set_ylim([0, y.size])
    
    ax3 = plt.subplot2grid((5, 5), (4, 1), colspan=4)
    ax3.plot(x)
    plt.setp(ax3.get_yticklabels(), visible=False)
    ax3.set_xlabel('Query axis')
    ax3.set_xlim([0, x.size])
    
    ax2 = plt.subplot2grid((5, 5), (0, 1), colspan=4, rowspan=4)
    ax2.plot(xaxis,yaxis, label=r"symmetric1",color='red')
    # setp( ax2.get_xticklabels(), visible=False)
    
    
    plt.legend(loc="upper left", bbox_to_anchor=[0, 1], ncol=1, shadow=True, title="Step Pattern", fancybox=True)
    ax2.get_legend().get_title().set_color("purple")
    ax2.set_xlim([0, x.size])
    ax2.set_ylim([0, y.size])
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax2.get_yticklabels(), visible=False)
        
    plt.suptitle("DTWPy: Time Series Alignment"+title)
    plt.show()

def plotalignment_with_window(window,path,title):
    plt.figure(num=None, figsize=(5, 5), dpi=80, facecolor='w', edgecolor='k')
    xlim,ylim = max(path,key=lambda item:item[1])
    wx = [x[0] for x in window]
    wy = [y[1] for y in window]
    x = [x[0] for x in path]
    y = [y[1] for y in path]
    plt.plot(wx,wy,color="gray")
    plt.plot(x,y,color="black")  
    plt.ylabel('Template index')
    plt.xlabel('Query index')
    plt.axis([0, xlim, 0, ylim])
    plt.suptitle("Itakura-window"+title)
    plt.show()
        
def plotalignment(path,title):
    plt.figure(num=None, figsize=(5, 5), dpi=80, facecolor='w', edgecolor='k')
    xlim,ylim = max(path,key=lambda item:item[1])
    x = [x[0] for x in path]
    y = [y[1] for y in path]
    plt.plot(x,y,color="black")  
    plt.ylabel('Template index')
    plt.xlabel('Query index')
    plt.axis([0, xlim, 0, ylim])
    plt.suptitle("Itakura-window"+title)
    plt.show()
        
def plotdetailalignment(query,template,path, title=''):
    coef=5
    xaxis = [x1[0] for x1 in path]
    yaxis = [y1[1] for y1 in path]
    for i in range(len(template)):
        template[i] += template[i]+coef

    plt.figure(1)
    plt.plot(template, lw=3)
    plt.plot(query, lw=3)
    for i in range(0,max(len(xaxis),len(yaxis))):
        plt.plot()
    #for val in path:
        #plt.plot([val[0], val[1]],[template[val[0]], query[val[1]]], 'k', lw=0.5)

    plt.axis('off')
    plt.show()
    

    
    