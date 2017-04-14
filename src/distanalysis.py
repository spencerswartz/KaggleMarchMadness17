import numpy as np
import csv
import pandas as pd
import random
from math import sqrt, pi, exp
import math

import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib inline')


# In[2]:

def pdf(mu=0, sigma=1):
    s2= sigma*sigma
    return lambda x: 1/sqrt(2*pi*(s2)) * exp(-(x-mu)**2/(2*s2))


# In[3]:

def plot1(frame,pdfW,pdfM,m_mu,m_sigma,w_mu,w_sigma, plot=False):
########################
    if plot:
        plt.figure(1)
        plt.plot(frame['theta'], frame['w'], '-')
        plt.plot(frame['theta'], frame['m'], '-')

        axis= [frame['theta'].min()
               ,frame['theta'].max()
               ,0
               ,min(max(frame['m'].max(),frame['w'].max()),1.5)
              ]

        w_atMu= pdfW(w_mu)
        m_atMu= pdfM(m_mu)

        w_atSigma= pdfW(w_sigma+w_mu)
        m_atSigma= pdfM(m_sigma+m_mu)

        plt.axvline(x=w_mu, ymax=w_atMu/axis[3], color='C1', linestyle='dashed')
        plt.axvline(x=m_mu, ymax=m_atMu/axis[3], color='C0', linestyle='dashed')

        plt.axvline(x=(w_sigma+w_mu), ymax=w_atSigma/axis[3], color='C0')
        plt.axvline(x=(m_sigma+m_mu), ymax=m_atSigma/axis[3], color='C1')

        plt.axis(axis)
        plt.title("Distribution")
        plt.legend()
        plt.show()
        ########################
        pass


# In[ ]:




# In[4]:

def plot2(frame,step,w_wing,m_wing, plot=False):
########################
    
    #frame['m_sum']= frame['m'].cumsum()*step+m_wing
    #frame['w_sum']= frame['w'].cumsum()*step+w_wing

    if plot:
        f, (ax1,ax2)= plt.subplots(2,1,sharex=True)
        ax1.plot(frame['theta'],frame['w_sum'])
        ax1.plot(frame['theta'],frame['m_sum'])
        plt.xticks([])
        

        ax1.legend()
        ax1.set_title('CDF')

    #frame['m_sum_inv']= 1-frame['m_sum']
    #frame['w_sum_inv']= 1-frame['w_sum']

    if plot:
        #plt.subplot(122)
        ax2.plot(frame['theta'],frame['w_sum_inv'])
        ax2.plot(frame['theta'],frame['m_sum_inv'])
        ax2.legend()
        ax2.set_title("Inverted-CDF")
        ########################

        plt.show()


# In[5]:

def plot3(frame, plot=False):
########################
    if plot:
        plt.figure(3)
        ax1= plt.subplot(211)
        plt.plot(frame['theta'],frame['w'])
        plt.plot(frame['theta'],frame['m'])

        axis= [frame['theta'].min()
               ,frame['theta'].max()
               ,0
               ,min(max(frame['m_sum'].max(),frame['w_sum'].max()),1.5)
              ]
        
        #ax1.axis(axis)
        plt.title("Density")
        plt.xticks([])
        plt.legend()
        #plt.show()

        plt.subplot(212)

        plt.plot(frame['theta'],frame['w_sum'])
        plt.plot(frame['theta'],frame['m_sum'])
        plt.title("Distribution")
        #plt.axis([-2,2,0,1])
        plt.legend()
        plt.show()
    ########################


# In[6]:

def plot4(frame, plot=False):
########################
    #frame['all']= frame['w_sum']+frame['m_sum']
    #frame['w_part']= frame['w_sum']/frame['all']
    #frame['m_part']= frame['m_sum']/frame['all']

    if plot:
        plt.figure(4)
        #f, (ax1, ax2)= plt.subplots(2, sharey=True)
        plt.subplot(211)
        plt.plot(frame['theta'],frame['w_part'])
        plt.plot(frame['theta'],frame['m_part'])
        plt.xticks([])
        plt.title("Integrate from Left")
        plt.legend()
        #plt.show()
    ########################

    ########################
    #frame['all_inv']= frame['w_sum_inv']+frame['m_sum_inv']
    #frame['w_part_inv']= frame['w_sum_inv']/frame['all_inv']
    #frame['m_part_inv']= frame['m_sum_inv']/frame['all_inv']

    if plot:
        #plt.figure(5)
        plt.subplot(212)
        plt.plot(frame['theta'],frame['w_part_inv'])
        plt.plot(frame['theta'],frame['m_part_inv'])
        plt.title("Integrate to Right")
        plt.legend()
        plt.show()
        ########################
        pass


# In[7]:

def plot5(frame, plot=False):
########################
    #frame['all']= frame['w_sum']+frame['m_sum']
    #frame['w_part']= frame['w_sum']/frame['all']
    #frame['m_part']= frame['m_sum']/frame['all']

    if plot:
        plt.figure(4)
        #f, (ax1, ax2)= plt.subplots(2, sharey=True)
        plt.subplot(211)
        plt.plot(frame['theta'],frame['w_part_kern'])
        plt.plot(frame['theta'],frame['m_part_kern'])
        plt.xticks([])
        plt.title("Peer Group Kernel Proportion (+/- 0.5)")
        plt.legend()
        #plt.show()
    ########################

    ########################
    #frame['all_inv']= frame['w_sum_inv']+frame['m_sum_inv']
    #frame['w_part_inv']= frame['w_sum_inv']/frame['all_inv']
    #frame['m_part_inv']= frame['m_sum_inv']/frame['all_inv']

    if plot:
        #plt.figure(5)
        plt.subplot(212)
        plt.plot(frame['theta'],frame['w_kernel'])
        plt.plot(frame['theta'],frame['m_kernel'])
        plt.title("Peer Group Kernel Density")
        plt.legend()
        plt.show()
        ########################
        pass


# In[8]:

def thetaIndex(theta=-7,step=.01):
    start= -abs(theta)
    stop= abs(theta)
    step = abs(step)
    nSteps= (stop-start)/step

    x= np.linspace(start, stop, int((nSteps-1)/2)*2+1)
    frame= pd.DataFrame({'theta':x})
    frame.reset_index(inplace=True)
    return frame


# In[9]:

def plotSamples(results):
    plt.figure(1)
    f, (ax1, ax2)= plt.subplots(2, sharex= True)
    
    if 'dMu' in results.columns:
        ax1.plot(results['dMu'],results['m_mu'])
        ax1.plot(results['dMu'],results['w_mu'])
        ax1.legend()
    else:
        ax1.plot(results['m_mu'],results['w_mu'],'.')
    #ax1.title('Mu')

    #plt.subplots(212)
    if 'dSigma' in results.columns:
        results.sort_values(by='dSigma', inplace=True)
        ax2.plot(results['dSigma'],results['m_sigma'])
        ax2.plot(results['dSigma'],results['w_sigma'])
        ax2.legend()
    else:
        ax2.plot(results['m_sigma'],results['w_sigma'],'.')
    #ax2.title('Sigma')
    plt.show()

    if 'dMu' in results.columns and 'dSigma' in results.columns:
        plt.figure(1)
        plt.plot(results['dMu'],results['dSigma'],'.')
        plt.xlabel('dMu')
        plt.ylabel('dSigma')
        plt.show()

    pass


# In[10]:

def simResults(baseMu=0,baseSigma=1,numRuns=420,seed=42,theta=-7,step=0.001):
    w_mu=baseMu
    w_sigma=baseSigma
    
    results= []#None
    random.seed(seed)
    for i in range(numRuns):
        #print(i)

        dMu= random.random()
        dSigma= random.random()
        dSigma_pctChange= math.copysign(abs(dSigma-0.5)+0.5,dSigma-0.5)

        #m_mu= w_mu+(dMu-0.5)
        #m_mu= w_mu+((dMu-0.5)/2)
        m_mu= w_mu+(dMu/2)*0
        m_sigma= w_sigma*(1+dSigma_pctChange)

        row= [w_mu,w_sigma,dMu,dSigma,dSigma_pctChange,m_mu,m_sigma,theta,step]
        #print(row)
        if type(results)!=type(None):
            results.append(row)
        else:
            results= pd.DataFrame(row)

    results= pd.DataFrame(results, columns=['w_mu','w_sigma','dMu','dSigma','dSigma_pctChange','m_mu','m_sigma','theta','step'])
    #results.head()
    return results


# In[11]:

#def plotData(m_mu,m_sigma,w_mu,w_sigma,frame,start,stop,step,plot=False):
def plotData(m_mu,m_sigma,w_mu,w_sigma,frame,step,plot=False):
    d= frame.copy()
    m_total,m_wing,pdfM,w_total,w_wing,pdfW,d= computeRunData(m_mu,m_sigma,w_mu,w_sigma,frame,step)
    #pdfM= pdf(m_mu,m_sigma)
    #pdfW= pdf(w_mu,w_sigma)
    #d['m']= d.applymap(pdfM)[['theta']]
    #d['w']= d.applymap(pdfW)[['theta']]
    #m_total= d['m'].sum()*step
    #w_total= d['w'].sum()*step
    #m_wing= (1-m_total)/2
    #w_wing= (1-w_total)/2
    
    d=computeStats(d,step,m_wing,w_wing)
    #if plot:
    plot1(d,pdfW,pdfM,m_mu,m_sigma,w_mu,w_sigma,plot=plot)
    plot2(d,step,w_wing,m_wing,plot=plot)
    plot3(d,plot=plot)
    plot4(d,plot=plot)
    
    return d


# In[12]:

def run(result, plot=False):
    #print("Running...", result.name)

    frame= thetaIndex(result.theta,result.step)
    frame= plotData(result.m_mu                    ,result.m_sigma
                    ,result.w_mu                    ,result.w_sigma
                    ,frame
                    ,result.step                    ,plot)
    frame['name']= result.name

    return frame


# In[13]:

def runAll(results, plot=False):
    cols= list(results.columns)
    if 'dMu' in cols:
        plotSamples(results)
    q= pd.DataFrame([])
    q= q.append([run(results.iloc[x]) for x in range(0,results.count()[0]) ])
    return q


# In[14]:

def plotProportion(z,col='m_part_inv',_printEvery=0):
    h,l = [],[]
    for i,r in enumerate(z['name'].unique()):
        if i % max(_printEvery,1) == (_printEvery==0): plt.show()

        #plt.axis([-5,5,0,1])
        plotDataFrame= z[z['name']==r]
        h0,= plt.plot(plotDataFrame['theta'],plotDataFrame[col])
        h.append(h0)
        l.append(r)
    return h,l


# In[15]:

def computeStats(data,step,m_wing,w_wing):
    frame = data.copy()
    frame.sort_values(by='theta', inplace=True)
    frame['m_sum']= frame['m'].cumsum()*step+m_wing
    frame['w_sum']= frame['w'].cumsum()*step+w_wing
    frame['all']= frame['w_sum']+frame['m_sum']
    frame['w_part']= frame['w_sum']/frame['all']
    frame['m_part']= frame['m_sum']/frame['all']

    frame['m_sum_inv']= 1-frame['m_sum']
    frame['w_sum_inv']= 1-frame['w_sum']
    frame['all_inv']= frame['w_sum_inv']+frame['m_sum_inv']
    frame['w_part_inv']= frame['w_sum_inv']/frame['all_inv']
    frame['m_part_inv']= frame['m_sum_inv']/frame['all_inv']

    #q[['theta']].applymap(lambda x: q[abs(q['theta']-q.ix[x]['theta']) <= 0.05]['m'].sum())
    #frame['m_kernel']= frame[['theta']].applymap(lambda x: frame[abs(frame['theta']-x) <= 0.5]['m'].sum()*step)
    frame['m_kernel']= frame[['theta']].applymap(lambda x: frame[abs(frame['theta']-x) <= 0.5]['m'].sum()*step)
    frame['w_kernel']= frame[['theta']].applymap(lambda x: frame[abs(frame['theta']-x) <= 0.5]['w'].sum()*step)
    frame['all_kernel']= frame['m_kernel']+frame['w_kernel']
    frame['m_part_kern']= frame['m_kernel']/frame['all_kernel']
    frame['w_part_kern']= frame['w_kernel']/frame['all_kernel']
    
    return frame


# In[16]:

def computeRunData(m_mu,m_sigma,w_mu,w_sigma,frame,step):
    d= frame.copy()
    pdfM= pdf(m_mu,m_sigma)
    pdfW= pdf(w_mu,w_sigma)
    d['m']= d.applymap(pdfM)[['theta']]
    d['w']= d.applymap(pdfW)[['theta']]
    m_total= d['m'].sum()*step
    w_total= d['w'].sum()*step
    m_wing= (1-m_total)/2
    w_wing= (1-w_total)/2
    
    return m_total,m_wing,pdfM,w_total,w_wing,pdfW,d
