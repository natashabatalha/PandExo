from bokeh.plotting import show, figure
from bokeh.io import output_file 
import pickle as pk
import numpy as np

def data(result_dict, model=True, title='Model + Data + Error Bars', outputfile = 'data.html',legend = False, 
        R=False,  num_tran = 1.0, plot_width=800, plot_height=400):
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
    colors = ['black','blue','red','orange','yellow','green','purple','pink','cyan','grey','brown']
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
    
        #remove any nans 
        y = dict['FinalSpectrum']['spectrum_w_rand']
        x = dict['FinalSpectrum']['wave'][~np.isnan(y)]
        err = dict['FinalSpectrum']['error_w_floor'][~np.isnan(y)]
        y = y[~np.isnan(y)]


        if R == False: 
            x=x 
            y=y 
        else:
            x, y ,err= binning(x,y,err,R)
            
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
         
            fig1d = figure(x_range=xlims, y_range = ylims, 
               plot_width = plot_width, plot_height =plot_height,title=title,x_axis_label=x_axis_label,
              y_axis_label = y_axis_label, tools=TOOLS)
        
              
        #plot model, data, and errors 
        if model:
            fig1d.line(dict['OriginalInput']['og_wave'], dict['OriginalInput']['og_spec'], color=colors[i],alpha=0.5)
        if legend: 
            fig1d.circle(x, y, color=colors[i], legend = legend_keys[i])
        else: 
            fig1d.circle(x, y, color=colors[i])
        fig1d.multi_line(x_err, y_err,color=colors[i])
        i += 1 
    show(fig1d)
    

def binning(x, y,e, R):
    """ 
    Takes 2 arrays x and y and bins them into groups of blength.
    only use for binning error bars. 
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
    e = np.array(e)

    xout,yout,eout= [],[],[]
    first = 0
    for i in ind:
        xout.append(sum(x[first:i])/len(x[first:i]))
        yout.append(sum(y[first:i])/len(y[first:i]))
        eout.append(np.sqrt(sum([jj**2 for jj in e[first:i]]))/len(e[first:i]))
        first = i 
    return np.array(xout),np.array(yout),np.array(eout)