from bokeh.plotting import show, figure
from bokeh.io import output_file as outputfile
import pickle as pk
import numpy as np

def jwst_1d_spec(result_dict, model=True, title='Model + Data + Error Bars', output_file = 'data.html',legend = False, 
        R=False,  num_tran = False, plot_width=800, plot_height=400,x_range=[1,10]):
    """Plots 1d simulated spectrum and rebin or rescale for more transits
    
    Plots 1d data points with model in the background (if wanted). Designed to read in exact 
    output of run_pandexo. 
    
    Parameters 
    ----------
    result_dict : dict or list of dict
        Dictionary from pandexo output. If parameter space was run in run_pandexo 
        make sure to restructure the input as a list of dictionaries without they key words 
        that run_pandexo assigns. 
    model : bool 
        (Optional) True is default. True plots model, False does not plot model     
    title : str
        (Optional) Title of plot. Default is "Model + Data + Error Bars".  
    output_file : str 
        (Optional) name of html file for you bokeh plot. After bokeh plot is rendered you will 
        have the option to save as png. 
    legend : bool 
        (Optional) Default is False. True, plots legend. 
    R : float 
        (Optional) Rebin data from native instrument resolution to specified resolution. Dafult is False, 
        no binning. 
    num_tran : float
        (Optional) Scales data by number of transits to improve error by sqrt(`num_trans`)
    plot_width : int 
        (Optional) Sets the width of the plot. Default = 800
    plot_height : int 
        (Optional) Sets the height of the plot. Default = 400 
    x_range : list of int
        (Optional) Sets x range of plot. Default = [1,10]

    Returns
    -------
    x,y,e : list of arrays 
        Returns wave axis, spectrum and associated error in list format. x[0] will be correspond 
        to the first dictionary input, x[1] to the second, etc. 
        
    Examples
    --------
    
    >>> jwst_1d_data(result_dict, num_tran = 3, R = 35) #for a single plot 
    
    If you wanted to save each of the axis that were being plotted: 
    
    >>> x,y,e = jwst_1d_data([result_dict1, result_dict2], model=False, num_tran = 5, R = 100) #for multiple 
    """
    outx=[]
    outy=[]
    oute=[]
    TOOLS = "pan,wheel_zoom,box_zoom,resize,reset,save"
    outputfile(output_file)
    colors = ['black','blue','red','orange','yellow','purple','pink','cyan','grey','brown']
    #make sure its iterable
    if type(result_dict) != list: 
        result_dict = [result_dict]
        
    if type(legend)!=bool:
        legend_keys = legend
        legend = True
        if type(legend_keys) != list:
            legend_keys = [legend_keys]
      
    i = 0     
    for dict in result_dict: 
        ntran_old = dict['timing']['Number of Transits']
        #remove any nans 
        y = dict['FinalSpectrum']['spectrum_w_rand']
        x = dict['FinalSpectrum']['wave'][~np.isnan(y)]
        err = dict['FinalSpectrum']['error_w_floor'][~np.isnan(y)]
        y = y[~np.isnan(y)]

        
        if (R == False) & (num_tran == False): 
            x=x 
            y=y 
        elif (R != False) & (num_tran != False):     
            new_wave = Rspec(x, R)
            out = uniform_tophat_sum(new_wave,x, dict['RawData']['flux_out']*num_tran/ntran_old)
            inn = uniform_tophat_sum(new_wave,x, dict['RawData']['flux_in']*num_tran/ntran_old)
            vout = uniform_tophat_sum(new_wave,x, dict['RawData']['var_out']*num_tran/ntran_old)
            vin = uniform_tophat_sum(new_wave,x, dict['RawData']['var_in']*num_tran/ntran_old)
            if dict['input']['Primary/Secondary']=='fp/f*':
                fac = -1.0
            else:
                fac = 1.0
            rand_noise = np.sqrt((vin+vout))*(np.random.randn(len(new_wave)))
            sim_spec = fac*(out-inn + rand_noise)/out 
            x = new_wave
            y = sim_spec
            err = np.sqrt(vout+vin)/out
        elif (R == False) & (num_tran != False):     
            out = dict['RawData']['flux_out']*num_tran/ntran_old
            inn = dict['RawData']['flux_in']*num_tran/ntran_old
            vout = dict['RawData']['var_out']*num_tran/ntran_old
            vin = dict['RawData']['var_in']*num_tran/ntran_old
            if dict['input']['Primary/Secondary']=='fp/f*':
                fac = -1.0
            else:
                fac = 1.0
            rand_noise = np.sqrt((vin+vout))*(np.random.randn(len(x)))
            sim_spec = fac*(out-inn + rand_noise)/out 
            x = x
            y = sim_spec
            err = np.sqrt(vout+vin)/out
        elif (R != False) & (num_tran == False):     
            new_wave = Rspec(x, R)
            out = uniform_tophat_sum(new_wave,x, dict['RawData']['flux_out'])
            inn = uniform_tophat_sum(new_wave,x, dict['RawData']['flux_in'])
            vout = uniform_tophat_sum(new_wave,x, dict['RawData']['var_out'])
            vin = uniform_tophat_sum(new_wave,x, dict['RawData']['var_in'])
            if dict['input']['Primary/Secondary']=='fp/f*':
                fac = -1.0
            else:
                fac = 1.0
            rand_noise = np.sqrt((vin+vout))*(np.random.randn(len(new_wave)))
            sim_spec = fac*(out-inn + rand_noise)/out 
            x = new_wave
            y = sim_spec
            err = np.sqrt(vout+vin)/out
        else: 
            print "Something went wrong. Cannot enter both resolution and ask to bin to new wave"
            return
            
        #create error bars for Bokeh's multi_line
        y_err = []
        x_err = []
        for px, py, yerr in zip(x, y, err):
            np.array(x_err.append((px, px)))
            np.array(y_err.append((py - yerr, py + yerr)))
        #initialize figure
        if i == 0: 
            #Define units for x and y axis
            y_axis_label = dict['input']['Primary/Secondary']

            if y_axis_label == 'fp/f*': p = -1.0
            else: y_axis_label = '('+y_axis_label+')^2'

            if dict['input']['Calculation Type'] =='phase_spec':
                x_axis_label='Time (secs)'
            else:
                x_axis_label='Wavelength [microns]'
            
            ylims = [min(dict['OriginalInput']['model_spec'])- 0.1*min(dict['OriginalInput']['model_spec']),
                 0.1*max(dict['OriginalInput']['model_spec'])+max(dict['OriginalInput']['model_spec'])]
            xlims = [min(x), max(x)]
         
            fig1d = figure(x_range=x_range, y_range = ylims, 
               plot_width = plot_width, plot_height =plot_height,title=title,x_axis_label=x_axis_label,
              y_axis_label = y_axis_label, tools=TOOLS, background_fill_color = 'white')
        
              
        #plot model, data, and errors 
        if model:
            mxx = dict['OriginalInput']['model_wave']
            myy = dict['OriginalInput']['model_spec']
            
            my = uniform_tophat_mean(x, mxx,myy)
            fig1d.line(x,my, color='white',alpha=0.2, line_width = 4)
        if legend: 
            fig1d.circle(x, y, color=colors[i], legend = legend_keys[i])
        else: 
            fig1d.circle(x, y, color=colors[i])
        outx += [x]
        outy += [y]
        oute += [err]
        fig1d.multi_line(x_err, y_err,color=colors[i])
        i += 1 
    show(fig1d)
    return outx,outy,oute

    
def bin_wave_to_R(w, R):
    """Creates new wavelength axis at specified resolution
    
    Parameters
    ----------
    w : list of float or numpy array of float
        Wavelength axis to be rebinned 
    R : float or int 
        Resolution to bin axis to 
    
    Returns
    -------
    list of float
        New wavelength axis at specified resolution
    
    Examples
    --------
    
    >>> newwave = bin_wave_to_R(np.linspace(1,2,1000), 10)
    >>> print len(newwave)
    11
    """
    wave = []
    tracker = min(w)
    i = 1 
    ind= 0
    while(tracker<max(w)):
        if i <len(w)-1:
        
            dlambda = w[i]-w[ind]
            newR = w[i]/dlambda
            if newR < R:
                tracker = w[ind]+dlambda/2.0
                wave +=[tracker]
                ind = (np.abs(w-tracker)).argmin()
                i = ind
            else:            
                i+=1    
        else:
            tracker = max(w)
    return wave

def uniform_tophat_sum(xnew,x, y):
    """Adapted from Mike R. Line to rebin spectra
    
    Takes sum of group of points in bin of wave points 
    Parameters
    ----------
    xnew : list of float or numpy array of float
        New wavelength grid to rebin to 
    x : list of float or numpy array of float 
        Old wavelength grid to get rid of 
    y : list of float or numpy array of float 
        New rebinned y axis 
    
    Returns
    -------
    array of floats 
        new y axis

    Examples 
    --------
        
    >>> oldgrid = np.linspace(1,3,100)
    >>> y = np.zeros(100)+10.0
    >>> newy = uniform_tophat_sum(np.linspace(2,3,3), oldgrid, y)
    >>> newy
    array([ 240.,  250.,  130.])
    """
    xnew = np.array(newx)
    szmod=newx.shape[0]
    delta=np.zeros(szmod)
    ynew=np.zeros(szmod)
    delta[0:-1]=newx[1:]-newx[:-1]  
    delta[szmod-1]=delta[szmod-2] 
    #pdb.set_trace()
    for i in range(szmod-1):
        i=i+1
        loc=np.where((x >= newx[i]-0.5*delta[i-1]) & (x < newx[i]+0.5*delta[i]))
        ynew[i]=np.sum(y[loc])
    loc=np.where((x > newx[0]-0.5*delta[0]) & (x < newx[0]+0.5*delta[0]))
    ynew[0]=np.sum(y[loc])
    return ynew

def uniform_tophat_mean(xnew,x, y):
    """Adapted from Mike R. Line to rebin spectra
    
    Takes average of group of points in bin
    
    Parameters
    ----------
    xnew : list of float or numpy array of float
        New wavelength grid to rebin to 
    x : list of float or numpy array of float 
        Old wavelength grid to get rid of 
    y : list of float or numpy array of float 
        New rebinned y axis 

    Returns
    -------
    array of floats 
        new y axis 

    Examples 
    --------
        
    >>> oldgrid = np.linspace(1,3,100)
    >>> y = np.zeros(100)+10.0
    >>> newy = uniform_tophat_sum(np.linspace(2,3,3), oldgrid, y)
    >>> newy
    array([ 240.,  250.,  130.])
    """
    xnew = np.array(newx)
    szmod=newx.shape[0]
    delta=np.zeros(szmod)
    ynew=np.zeros(szmod)
    delta[0:-1]=newx[1:]-newx[:-1]  
    delta[szmod-1]=delta[szmod-2] 
    #pdb.set_trace()
    for i in range(szmod-1):
        i=i+1
        loc=np.where((x >= newx[i]-0.5*delta[i-1]) & (x < newx[i]+0.5*delta[i]))
        ynew[i]=np.mean(y[loc])
    loc=np.where((x > newx[0]-0.5*delta[0]) & (x < newx[0]+0.5*delta[0]))
    ynew[0]=np.mean(y[loc])
    return ynew
