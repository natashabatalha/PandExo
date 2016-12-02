from bokeh.plotting import show, figure
from bokeh.io import output_file 
import pickle as pk
import numpy as np

def data(result_dict, model=True, title='Model + Data + Error Bars', outputfile = 'data.html',legend = False, 
        R=False,  num_tran = False, plot_width=800, plot_height=400,x_range=[1,10], new_wave=False):
    """
    Plots 1d data points with model in the background (if wanted) 
    
    Inputs: 
        result_dict: dictionary from pandexo output or list of dictionaries 
        
        Optional: 
        output_file: html filename of your output (must be something.html). can save as png later. 
        title: title of plot (default = Model + Data + Error Bars) 
        legend: if you have multiple input dicts you can specify legends for them as a list ['thing1','thing2']
        fixed_bin: add in fixed bin width in units of micron 
        num_tran: add in number of transits 
    Outputs: 
        no outputs, renders html plot 
        
    Example: 
        import pandexo.engine.justplotit as jpi 
        jpi.data(result_dict) #for a single plot 
        jpi.data([result_dict1, result_dict2], color_data= ['red','blue']) #for multiple 
    """
    TOOLS = "pan,wheel_zoom,box_zoom,resize,reset,save"
    output_file(outputfile)
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
            new_wave = Rspec(x, R)
            out = dict['RawData']['flux_out']*num_tran/ntran_old
            inn = dict['RawData']['flux_in']*num_tran/ntran_old
            vout = new_wave,dict['RawData']['var_out']*num_tran/ntran_old
            vin = new_wave,dict['RawData']['var_in']*num_tran/ntran_old
            if dict['input']['Primary/Secondary']=='fp/f*':
                fac = -1.0
            else:
                fac = 1.0
            rand_noise = np.sqrt((vin+vout))*(np.random.randn(len(new_wave)))
            sim_spec = fac*(out-inn + rand_noise)/out 
            x = new_wave
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
            
            ylims = [min(dict['OriginalInput']['og_spec'])- 0.1*min(dict['OriginalInput']['og_spec']),
                 0.1*max(dict['OriginalInput']['og_spec'])+max(dict['OriginalInput']['og_spec'])]
            xlims = [min(x), max(x)]
         
            fig1d = figure(x_range=x_range, y_range = ylims, 
               plot_width = plot_width, plot_height =plot_height,title=title,x_axis_label=x_axis_label,
              y_axis_label = y_axis_label, tools=TOOLS, background_fill_color = 'white')
        
              
        #plot model, data, and errors 
        if model:
            mxx = dict['OriginalInput']['og_wave']
            myy = dict['OriginalInput']['og_spec'][mxx>1]
            mxx = mxx[mxx>1]
            mx, my = bin_data_smoothe(mxx,myy , 50)
            fig1d.line(mx,my, color='white',alpha=0.2, line_width = 4)
        if legend: 
            fig1d.circle(x, y, color=colors[i], legend = legend_keys[i])
        else: 
            fig1d.circle(x, y, color=colors[i])
        fig1d.multi_line(x_err, y_err,color=colors[i])
        i += 1 
    show(fig1d)
    
"""
def bin_data(x,y,R):
 
    Takes 2 arrays x and y and bins them into groups of blength.
    
    Parameters
	----------
        Inputs:     
            -x, y:                   1D lists or numpy arrays
        Outputs:    
            - xout, yout, yerrout,noise:    1D numpy arrays

    wlength = min(x)/R
    ii = 0
    start = 0
    ind = []
    for i in range(0, len(x)-1):
        if x[i+1] - x[start] >= wlength:
            ind.append(i+1)
            start = i 
    
    if ind[len(ind)-1] != (len(x)):
        ind.append(len(x))
    
    # convert to arrays if necessary
    x = np.array(x)
    y = np.array(y)


    xout,yout= [],[]
    first = 0
    for i in ind:
        xout.append(sum(x[first:i])/len(x[first:i]))
        yout.append(sum(y[first:i]))
        first = i 
    return np.array(xout),np.array(yout)
"""    
def bin_data_smoothe(x,y,R):
    """ 
    Takes 2 arrays x and y and bins them into groups of blength.
    
    Parameters
	----------
        Inputs:     
            -x, y:                   1D lists or numpy arrays
        Outputs:    
            - xout, yout, yerrout,noise:    1D numpy arrays
    """
    wlength = min(x)/R
    ii = 0
    start = 0
    ind = []
    for i in range(0, len(x)-1):
        if x[i+1] - x[start] >= wlength:
            ind.append(i+1)
            start = i 
    
    if ind[len(ind)-1] != (len(x)):
        ind.append(len(x))
    
    # convert to arrays if necessary
    x = np.array(x)
    y = np.array(y)


    xout,yout= [],[]
    first = 0
    for i in ind:
        xout.append(sum(x[first:i])/len(x[first:i]))
        yout.append(sum(y[first:i])/len(x[first:i]))
        first = i 
    return np.array(xout),np.array(yout)
    
def Rspec(a, R):
    """
    given a wavelength axis and a resolution, rebins just the axis 
    to the correct resolution
    """
    wave = []
    tracker = min(a)
    i = 1 
    ind= 0
    while(tracker<max(a)):
        if i <len(a)-1:
        
            dlambda = a[i]-a[ind]
            newR = a[i]/dlambda
            if newR < R:
                tracker = a[ind]+dlambda/2.0
                wave +=[tracker]
                ind = (np.abs(a-tracker)).argmin()
                i = ind
            else:            
                i+=1    
        else:
            tracker = max(a)
    return wave
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]   

def uniform_tophat_sum(wlgrid,wno, Fp):
    wlgrid = np.array(wlgrid)
    szmod=wlgrid.shape[0]

    delta=np.zeros(szmod)
    Fint=np.zeros(szmod)
    delta[0:-1]=wlgrid[1:]-wlgrid[:-1]  
    delta[szmod-1]=delta[szmod-2] 
    #pdb.set_trace()
    for i in range(szmod-1):
        i=i+1
        loc=np.where((wno >= wlgrid[i]-0.5*delta[i-1]) & (wno < wlgrid[i]+0.5*delta[i]))
        Fint[i]=np.sum(Fp[loc])
	    
    loc=np.where((wno > wlgrid[0]-0.5*delta[0]) & (wno < wlgrid[0]+0.5*delta[0]))
    Fint[0]=np.sum(Fp[loc])
    
    return Fint
     