#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 22:55:38 2021

@author: fadeli
"""
import numpy as np
#----------------------------------------------------------
def getamp(rfdata, tarray, t):
    amp = rfdata[(np.abs(tarray - t).argmin())]
    return amp
#----------------------------------------------------------

def hkZhu(rfstream,savePath, ext = "rf", H= np.linspace(20,60,201), 
          K = np.linspace(1.65,1.95,121), Vp = 6.5, rmneg = None, w1 = 0.7, 
          w2 = 0.2, w3 = 0.1, g = [75.,10., 15., 2.5], suff = "", 
          figleg = ["(A)","(B)","(C)"], legloc = [58,1.945,28,380,28,108] ):
    from os import chdir, getcwd
    from obspy import read
    from obspy.io.sac import SACTrace
    from glob import glob
    #from scipy.interpolate import griddata
    import matplotlib.pyplot as plt
    import matplotlib.transforms as mtransforms
    from matplotlib.legend_handler import HandlerLine2D as HL
    import matplotlib.gridspec as gridspec
    from obspy.geodetics.base import kilometer2degrees
    import seaborn
    seaborn.set_style("darkgrid", {"axes.facecolor": ".85"})

#    rfProc = "/home/islam/NARSPhD/RF/PROCESS/Pup/Gao/AF.MAUN/08_STACK/1.0_R"
#    savePath = "/home/islam/NARSPhD/RF/RESULTS"

    # if g is None: g = [75.,10., 15., 2.5]
    # if Vp is None: Vp = 6.5
    # if w1 is None: w1 = 0.6
    # if w2 is None: w2 = 0.3
    # if w3 is None: w3 = 0.1
    # if suff is None: suff = ""
    # if ext is None: ext = "rf"    
    # if rmneg is None: rmneg = False
    # if figleg is None: figleg=["(A)","(B)","(C)"]
    # if legloc is None: legloc = [58,1.945,28,380,28,108]    
    
    staName = rfstream[0].stats.station
    gau = rfstream[0].stats['gaussian']
    comp =  rfstream[0].stats.channel
    
    sta = "%s %s %s"%(staName, gau, comp) 
    print(sta)
    
    rfdata = rfstream[0].data
    delta = rfstream[0].stats.delta
    l = len(rfdata)
    t = np.arange(0, l)
    t = (delta *  t) - 10 
       
    # if  H is None: H = np.linspace(20,60,201)
    print(H)
    # if K is None: K = np.linspace(1.65,1.95,121)
    print(K)
    
    stk = np.zeros((len(K)*len(H),3))
    dbgvalues = np.zeros((len(K)*len(H)*len(rfstream),9))
    z = 0 ; q =0
    for i in range(len(K)):
        Ktemp = K[i]
        for j in range(len(H)):
            Htemp = H[j]
            s = 0.0
            for tr in rfstream:
                trdata = tr.data
                # b = tr.stats.sac.b
                b = tr.stats.starttime - tr.stats.onset
                delta= tr.stats.delta
                # rp = tr.stats.sac.user4
                rp = kilometer2degrees(1) * tr.stats['slowness']
                # gauw = tr.stats.sac.user0
                gauw = tr.stats['gaussian']
                
                #--------------------------------------
                term1= ((Ktemp/Vp)**2 - rp**2)**0.5
                term2= ((1/Vp)**2 - rp**2)**0.5
                tps = Htemp * (term1 - term2)
                tppps = Htemp * (term1 + term2)
                tpsps = Htemp * 2 * (term1) 
                #--------------------------------------
                ampps = getamp(trdata, t, tps)
                ampppps = getamp(trdata,t, tppps)
                amppsps = getamp(trdata, t, tpsps)
                #--------------------------------------
                stemp = (w1 * ampps) + (w2 * ampppps) - (w3 * amppsps)
                s = stemp + s
                dbgvalues[q, :] =[rp, Htemp, Ktemp, tps, tppps, tpsps,\
                ampps, ampppps, amppsps]
                q =q + 1
            stk[z,:] = [Htemp, Ktemp, s]
            z = z + 1
    bmod = stk[np.argmax(stk[:,2]),:]
    #----------------------------------------------------------
    # Add Ps, PpPs, PsPs + PbSs timing to the sac header
    #----------------------------------------------------------    
    Hbm = bmod[0]
    Kbm = bmod[1]
    Sbm = bmod[2]
    print("Best depth: ", Hbm, "Best Vp/Vs:", Kbm, "Max stack: ", Sbm)
    
    for tr in rfstream:
        # trhead = SACTrace.read(i, headonly = True)
        # rp = trhead.user4
        rp = kilometer2degrees(1) * tr.stats['slowness']
        term1= ((Kbm/Vp)**2 - rp**2)**0.5
        term2= ((1/Vp)**2 - rp**2)**0.5
        
        tps = Hbm * (term1 - term2)
        tppps = Hbm * (term1 + term2)
        tpsps = Hbm * 2 * (term1)
#        print tps, tppps, tpsps
        
        tr.stats['t6'] = tps
        tr.stats['t7'] = tppps
        tr.stats['t8'] = tpsps
#         trhead.write(i, headonly = True)
    
    #----------------------------------------------------------
    # plot H-K results
    #----------------------------------------------------------
    plt.figure(figsize=(18, 6)) 
    gs = gridspec.GridSpec(1, 3, width_ratios=[6,4,4])
    
    #----------------------------------------------------------
    # plot misfit
    #----------------------------------------------------------
    #-----Remove the stack values that is lower than zero------
    if rmneg == True:
        for i in range(len(stk)):
            if stk[i,2] <= 0:
                stk[i,2] = 0
    #----------------------------------------------------------
    ax1 = plt.subplot(gs[0, 0])
    plt.tricontourf(stk[:,0], stk[:,1], stk[:,2],50, cmap ="jet")
    #-----Other plotting options------
    #plt.scatter(stk[:,0], stk[:,1], c = stk[:,2], s = 10, lw = 0.0)
    #X,Y= np.meshgrid(H,K)
    #Z = griddata((stk[:,0],stk[:,1]), stk[:,2], (X, Y), method='cubic')
    #plt.contour(X,Y,Z,10)
    #----------------------------------------------------------
    cb=plt.colorbar(format='%.3f', orientation="vertical")
    cb.ax.set_xlabel('S', rotation=0)
    p1,=ax1.plot(bmod[0],bmod[1], 'k+', mew=2, ms=10, \
    label='%s Best Model %.2f km %.2f Vp/Vs'%(staName,bmod[0],bmod[1]+0.001))

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax1.text(legloc[0], legloc[1], figleg[0], fontsize=12, \
    bbox=props, fontweight ='bold', ha = 'center', va = 'top')

    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,\
        ncol=2, mode="expand", borderaxespad=0, \
        handler_map={p1:HL(numpoints=1)})
    ax1.set_ylabel('Vp/Vs')
    ax1.set_xlabel('Depth km')
    ax1.set_xlim(min(H), max(H))
    ax1.set_ylim(min(K), max(K))
    #----------------------------------------------------------
    # Plot the receiver function witht the estimated times
    #----------------------------------------------------------
    ax2 = plt.subplot(gs[0, 1]) 
    gcarcarr = np.zeros(len(rfstream))
    for rf in rfstream:
        # rf = SACTrace.read(rffiles[i])
        
        delta = rf.stats.delta
        baz = rf.stats.back_azimuth
        b = rf.stats.starttime - rf.stats.onset
        # a = rf.user0
    
        t6 = rf.stats.t6
        t7 = rf.stats.t7
        t8 = rf.stats.t8
    
        rfdata = ((rf.data) * g[0]) + baz
        
        t = (np.arange(0, len(rfdata), 1) * delta) + b
        
        major_ticks_x = np.arange(-10, 41, 5)                                              
        minor_ticks_x = np.arange(-10, 41, 1) 
        
        major_ticks_y = np.arange(0, 361, 45)                                              
        minor_ticks_y = np.arange(-30, 390, 5)
        
        ax2.set_xticks(major_ticks_x)                                                       
        ax2.set_xticks(minor_ticks_x, minor=True)                                           
        ax2.set_yticks(major_ticks_y)                                                       
        ax2.set_yticks(minor_ticks_y, minor=True)                                          
                                                           
    #    ax2.grid()                                                           
    #                                
    #    ax2.grid(which='minor', alpha=0.5)                                                
    #    ax2.grid(which='major', alpha=1.0) 
        
        ax2.set_xlim(-5, 30)
        ax2.set_ylim(-29, 390)
        
        ax2.plot(t, rfdata, "k-", lw=0.5, label = "RF Gaussian %s"%gau)
        ax2.fill_between(t, baz, rfdata, where=rfdata > baz, 
                        facecolor='red', alpha = 0.25)
        ax2.fill_between(t, baz, rfdata, where=rfdata <= baz, 
                        facecolor='blue', alpha = 0.25)
        mtransforms.blended_transform_factory(ax2.transData, ax2.transAxes)
    #    #----------------------------------------------------------
    #    # Plot H-K sequential formula times
    #    #----------------------------------------------------------
    #    for i in range(0,30,5):
    #        ax2.plot([i, i], [-30, 390], 'k-', lw=0.25)
        ampgain = g[1]
        ax2.plot([t6, t6], [baz+ampgain, baz-ampgain], 'g-', lw=1.5, \
                                                        label="Moho")
        ax2.plot([t7, t7], [baz+ampgain, baz-ampgain], 'g-', lw=1.5)
        ax2.plot([t8, t8], [baz+ampgain, baz-ampgain], 'g-', lw=1.5)
    
    #----------------------------------------------------------
    # Highlight specifiec time range on RF plot
    #----------------------------------------------------------
    ax2.axvspan(t6-1.5, t6+1.5, alpha=0.5, color='lightgreen')
    ax2.axvspan(t7-1.5, t7+1.5, alpha=0.5, color='lightgreen')
    ax2.axvspan(t8-1.5, t8+1.5, alpha=0.5, color='lightgreen')
    
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax2.text(legloc[2], legloc[3], figleg[1], fontsize=12, \
    bbox=props, fontweight ='bold', ha = 'center', va = 'top')
    
    ax2.set_xlabel("Time (sec)")
    ax2.set_ylabel("Back-Azimuth (deg)")
    #ax2.set_title("The H-K of station %s" %sta)
    
    
    #----------------------------------------------------------
    # Plot the receiver function witht the estimated times
    #----------------------------------------------------------
    ax3 = plt.subplot(gs[0, 2])
    gcarcarr = np.zeros(len(rfstream))   
    for rf, j in zip(rfstream, np.arange(len(rfstream))):
        # rf = SACTrace.read(rffiles[i])
        
        delta = rf.stats.delta
        gcarc = rf.stats.distance
        gcarcarr[j] = gcarc
        b = rf.stats.starttime - rf.stats.onset
        # a = rf.user0
    
        t6 = rf.stats.t6
        t7 = rf.stats.t7
        t8 = rf.stats.t8
        
        rfdata = ((rf.data) * g[2]) + gcarc
        
        t = (np.arange(0, len(rfdata), 1) * delta) + b    
        
        major_ticks_x = np.arange(-10, 41, 5)                                              
        minor_ticks_x = np.arange(-10, 41, 1) 
        
        major_ticks_y = np.arange(30, 110, 10)                                              
        minor_ticks_y = np.arange(30, 111, 5)
        
        ax3.set_xticks(major_ticks_x)                                                       
        ax3.set_xticks(minor_ticks_x, minor=True)                                           
        ax3.set_yticks(major_ticks_y)                                                       
        ax3.set_yticks(minor_ticks_y, minor=True)                                          
    
        #----------------------------------------------------------
        # Plot grid on the figure
        #----------------------------------------------------------                                                       
    #    ax3.grid()                                                           
    #                                
    #    ax3.grid(which='minor', alpha=0.5)                                                
    #    ax3.grid(which='major', alpha=1.0) 
        
        ax3.set_xlim(-5, 30)
        ax3.set_ylim(25,110)
        
        ax3.plot(t, rfdata, "k-", lw=0.5, label = "RF Gaussian %s"%gau)
        ax3.fill_between(t, gcarc, rfdata, where=rfdata > gcarc, 
                        facecolor='red', alpha = 0.25)
        ax3.fill_between(t, gcarc, rfdata, where=rfdata <= gcarc, 
                        facecolor='blue', alpha = 0.25)
        mtransforms.blended_transform_factory(ax3.transData, ax3.transAxes)
        #----------------------------------------------------------
        # Plot H-K sequential formula times
        #----------------------------------------------------------
    #    for i in range(0,30,5):
    #        ax3.plot([i, i], [20, 110], 'k-', lw=0.1)
        ampgain = g[3]
        ax3.plot([t6, t6], [gcarc+ampgain, gcarc-ampgain], 'g-', \
                                            lw=1.5, label="Moho")
        ax3.plot([t7, t7], [gcarc+ampgain, gcarc-ampgain], 'g-', \
                                                          lw=1.5)
        ax3.plot([t8, t8], [gcarc+ampgain, gcarc-ampgain], 'g-', \
                                                          lw=1.5)
    #----------------------------------------------------------
    # Highlight specifiec time range on RF plot
    #----------------------------------------------------------
    ax3.axvspan(t6-1.5, t6+1.5, alpha=0.5, color='lightgreen')
    ax3.axvspan(t7-1.5, t7+1.5, alpha=0.5, color='lightgreen')
    ax3.axvspan(t8-1.5, t8+1.5, alpha=0.5, color='lightgreen')    

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax3.text(legloc[4], legloc[5], figleg[2], fontsize=12, \
    bbox=props, fontweight ='bold', ha = 'center', va = 'top')
    
    ax3.set_xlabel("Time (sec)")
    ax3.set_ylabel("Distance (deg)")
    #ax3.set_title("The H-K results of station %s" %sta)
    
    filename = "%s/%s_%s_%s_%s.pdf"%(savePath, staName, gau, comp, suff)
    
    plt.savefig(filename , format='pdf', transparent=False,\
        dpi=250, bbox_inches = 'tight', pad_inches=0.1)
    
    # chdir(cwd)
    

#%% Test the function
# from rf import read_rf
# import os
# import matplotlib.pyplot as plt
# home = os.path.join('/Users/fadeli/Work/Codes/rf/test/', '99_test/')
# rfgood_visqc = home + '02_rfgood_vis.h5'
# rfst = read_rf(rfgood_visqc, 'H5')

# kw = {'trim': (-5, 40), 'fillcolors': ('black', 'gray'), 'trace_height': 0.1}
# rfst.select(component='R', station='NE208').sort(['back_azimuth']).plot_rf(**kw)
# plt.show()

# rfn08 = rfst.select(component='R', station='NE208')
# hkZhu(rfn08, home)
    
    
    
    