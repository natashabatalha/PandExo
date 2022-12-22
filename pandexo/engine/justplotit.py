from bokeh.plotting import show
from bokeh.plotting import figure as Figure
from bokeh.io import output_file as outputfile
from bokeh.io import output_notebook  as outnotebook
import pickle as pk
import numpy as np
from bokeh.layouts import row
import pandas as pd
def jwst_1d_spec(result_dict, model=True, title='Model + Data + Error Bars', output_file = 'data.html',legend = False,
        R=False,  num_tran = False, plot_width=800, plot_height=400,x_range=[1,10],y_range=None, plot=True,
        output_notebook=True):
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
        no binning. Here I adopt R as w[1]/(w[2] - w[0]) to maintain consistency with `pandeia.engine`
    num_tran : float
        (Optional) Scales data by number of transits to improve error by sqrt(`num_trans`)
    plot_width : int
        (Optional) Sets the width of the plot. Default = 800
    plot_height : int
        (Optional) Sets the height of the plot. Default = 400
    y_range : list of int
        (Optional) sets y range of plot. Defaut is +- 10% of max and min
    x_range : list of int
        (Optional) Sets x range of plot. Default = [1,10]
    plot : bool
        (Optional) Supresses the plot if not wanted (Default = True)
    out_notebook : bool 
        (Optional) Output notebook. Default is false, if true, outputs in the notebook

    Returns
    -------
    x,y,e : list of arrays
        Returns wave axis, spectrum and associated error in list format. x[0] will be correspond
        to the first dictionary input, x[1] to the second, etc.

    Examples
    --------

    >>> jwst_1d_spec(result_dict, num_tran = 3, R = 35) #for a single plot

    If you wanted to save each of the axis that were being plotted:

    >>> x,y,e = jwst_1d_data([result_dict1, result_dict2], model=False, num_tran = 5, R = 100) #for multiple

    See Also
    --------
    jwst_noise, jwst_1d_bkg, jwst_1d_flux, jwst_1d_snr, jwst_2d_det, jwst_2d_sat

    """
    outx=[]
    outy=[]
    oute=[]
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save"
    if output_notebook & plot:
        outnotebook()
    elif plot:
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
    for dictt in result_dict:
        ntran_old = dictt['timing']['Number of Transits']
        to = dictt['timing']["Num Integrations Out of Transit"]
        ti = dictt['timing']["Num Integrations In Transit"]
        #remove any nans
        y = dictt['FinalSpectrum']['spectrum_w_rand']
        x = dictt['FinalSpectrum']['wave'][~np.isnan(y)]
        err = dictt['FinalSpectrum']['error_w_floor'][~np.isnan(y)]
        y = y[~np.isnan(y)]


        if (R == False) & (num_tran == False):
            x=x
            y=y
        elif (R != False) & (num_tran != False):
            new_wave = bin_wave_to_R(x, R)
            out = uniform_tophat_sum(new_wave,x, dictt['RawData']['electrons_out']*num_tran/ntran_old)
            inn = uniform_tophat_sum(new_wave,x, dictt['RawData']['electrons_in']*num_tran/ntran_old)
            vout = uniform_tophat_sum(new_wave,x, dictt['RawData']['var_out']*num_tran/ntran_old)
            vin = uniform_tophat_sum(new_wave,x, dictt['RawData']['var_in']*num_tran/ntran_old)
            var_tot = (to/ti/out)**2.0 * vin + (inn*to/ti/out**2.0)**2.0 * vout
            if dictt['input']['Primary/Secondary']=='fp/f*':
                fac = -1.0
            else:
                fac = 1.0
            rand_noise = np.sqrt((var_tot))*(np.random.randn(len(new_wave)))
            raw_spec = (out/to-inn/ti)/(out/to)
            sim_spec = fac*(raw_spec + rand_noise )
            x = new_wave
            y = sim_spec
            err = np.sqrt(var_tot)
        elif (R == False) & (num_tran != False):
            out = dictt['RawData']['electrons_out']*num_tran/ntran_old
            inn = dictt['RawData']['electrons_in']*num_tran/ntran_old
            vout = dictt['RawData']['var_out']*num_tran/ntran_old
            vin = dictt['RawData']['var_in']*num_tran/ntran_old
            var_tot = (to/ti/out)**2.0 * vin + (inn*to/ti/out**2.0)**2.0 * vout
            if dictt['input']['Primary/Secondary']=='fp/f*':
                fac = -1.0
            else:
                fac = 1.0
            rand_noise = np.sqrt((var_tot))*(np.random.randn(len(x)))
            raw_spec = (out/to-inn/ti)/(out/to)
            sim_spec = fac*(raw_spec + rand_noise )
            x = x
            y = sim_spec
            err = np.sqrt(var_tot)
        elif (R != False) & (num_tran == False):
            new_wave = bin_wave_to_R(x, R)
            out = uniform_tophat_sum(new_wave,x, dictt['RawData']['electrons_out'])
            inn = uniform_tophat_sum(new_wave,x, dictt['RawData']['electrons_in'])
            vout = uniform_tophat_sum(new_wave,x, dictt['RawData']['var_out'])
            vin = uniform_tophat_sum(new_wave,x, dictt['RawData']['var_in'])
            var_tot = (to/ti/out)**2.0 * vin + (inn*to/ti/out**2.0)**2.0 * vout
            if dictt['input']['Primary/Secondary']=='fp/f*':
                fac = -1.0
            else:
                fac = 1.0
            rand_noise = np.sqrt((var_tot))*(np.random.randn(len(new_wave)))
            raw_spec = (out/to-inn/ti)/(out/to)
            sim_spec = fac*(raw_spec + rand_noise )
            x = new_wave
            y = sim_spec
            err = np.sqrt(var_tot)
        else:
            print("Something went wrong. Cannot enter both resolution and ask to bin to new wave")
            return
        
        x = np.array(x,dtype=float)
        y = np.array(y,dtype=float)
        err= np.array(err,dtype=float)
        
        #create error bars for Bokeh's multi_line and drop nans
        data = pd.DataFrame({'x':x, 'y':y,'err':err}).dropna()

        y_err = []
        x_err = []
        for px, py, yerr in zip(data['x'], data['y'], data['err']):
            np.array(x_err.append((px, px)))
            np.array(y_err.append((py - yerr, py + yerr)))
        #initialize Figure
        if i == 0:
            #Define units for x and y axis
            y_axis_label = dictt['input']['Primary/Secondary']

            if y_axis_label == 'fp/f*': p = -1.0
            else: y_axis_label = y_axis_label

            if dictt['input']['Calculation Type'] =='phase_spec':
                x_axis_label='Time (secs)'
                x_range = [min(x), max(x)]
            else:
                x_axis_label='Wavelength [microns]'

            if y_range!=None:
                ylims = y_range
            else:
                ylims = [min(dictt['OriginalInput']['model_spec'])- 0.1*min(dictt['OriginalInput']['model_spec']),
                 0.1*max(dictt['OriginalInput']['model_spec'])+max(dictt['OriginalInput']['model_spec'])]

            fig1d = Figure(x_range=x_range, y_range = ylims,
               width = plot_width, height =plot_height,title=title,x_axis_label=x_axis_label,
              y_axis_label = y_axis_label, tools=TOOLS, background_fill_color = 'white')


        #plot model, data, and errors
        if model:
            mxx = dictt['OriginalInput']['model_wave']
            myy = dictt['OriginalInput']['model_spec']
            my = uniform_tophat_mean(x, mxx,myy)
            model_line = pd.DataFrame({'x':x, 'my':my}).dropna()
            fig1d.line(model_line['x'],model_line['my'], color=colors[i],alpha=0.2, line_width = 4)


        if legend:
            fig1d.circle(data['x'], data['y'], color=colors[i], legend = legend_keys[i])
        else:
            fig1d.circle(data['x'], data['y'], color=colors[i])
        outx += [data['x'].values]
        outy += [data['y'].values]
        oute += [data['err'].values]
        fig1d.multi_line(x_err, y_err,color=colors[i])
        i += 1
    if plot:
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
    >>> print(len(newwave))
    11
    """
    wave = []
    tracker = min(w)
    i = 1
    ind= 0
    firsttime = True
    while(tracker<max(w)):
        if i <len(w)-1:
            dlambda = w[i]-w[ind]
            newR = w[i]/dlambda
            if (newR < R) & (firsttime):
                tracker = w[ind]
                wave += [tracker]
                ind += 1
                i += 1
                firsttime = True
            elif newR < R:
                tracker = w[ind]+dlambda/2.0
                wave +=[tracker]
                ind = (np.abs(w-tracker)).argmin()
                i = ind+1
                firsttime = True
            else:
                firsttime = False
                i+=1
        else:
            tracker = max(w)
            wave += [tracker]
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
    xnew = np.array(xnew)
    szmod=xnew.shape[0]
    delta=np.zeros(szmod)
    ynew=np.zeros(szmod)
    delta[0:-1]=xnew[1:]-xnew[:-1]
    delta[szmod-1]=delta[szmod-2]
    #pdb.set_trace()
    for i in range(szmod-1):
        i=i+1
        loc=np.where((x >= xnew[i]-0.5*delta[i-1]) & (x < xnew[i]+0.5*delta[i]))
        ynew[i]=np.sum(y[loc])
    loc=np.where((x > xnew[0]-0.5*delta[0]) & (x < xnew[0]+0.5*delta[0]))
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
    xnew = np.array(xnew)
    szmod=xnew.shape[0]
    delta=np.zeros(szmod)
    ynew=np.zeros(szmod)
    delta[0:-1]=xnew[1:]-xnew[:-1]
    delta[szmod-1]=delta[szmod-2]
    #pdb.set_trace()
    for i in range(szmod-1):
        i=i+1
        loc=np.where((x >= xnew[i]-0.5*delta[i-1]) & (x < xnew[i]+0.5*delta[i]))
        ynew[i]=np.mean(y[loc])
    loc=np.where((x > xnew[0]-0.5*delta[0]) & (x < xnew[0]+0.5*delta[0]))
    ynew[0]=np.mean(y[loc])
    return ynew

def jwst_1d_flux(result_dict, plot=True, output_file= 'flux.html'):
    """Plot flux rate in e/s

    Parameters
    ----------
    result_dict : dict
        Dictionary from pandexo output. If parameter space was run in run_pandexo
        make sure to restructure the input as a list of dictionaries without they key words
        that run_pandexo assigns.
    plot : bool
        (Optional) True renders plot, Flase does not. Default=True
    output_file : str
        (Optional) Default = 'flux.html'
    Return
    ------
    x : numpy array
        micron
    y : numpy array
        1D flux rate in electrons/s

    See Also
    --------
    jwst_1d_spec, jwst_1d_bkg, jwst_noise, jwst_1d_snr, jwst_2d_det, jwst_2d_sat
    """
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save"
    out = result_dict['PandeiaOutTrans']

    # Flux 1d
    x, y = out['1d']['extracted_flux']
    x = x[~np.isnan(y)]
    y = y[~np.isnan(y)]

    plot_flux_1d1 = Figure(tools=TOOLS,
                         x_axis_label='Wavelength [microns]',
                         y_axis_label='Flux (e/s)', title="Out of Transit Flux Rate",
                         width=800, height=300)
    plot_flux_1d1.line(x, y, line_width = 4, alpha = .7)

    if plot:
        outputfile(output_file)
        show(plot_flux_1d1)
    return x,y

def jwst_1d_snr(result_dict, plot=True, output_file='snr.html'):
    """Plot SNR

    Parameters
    ----------
    result_dict : dict
        Dictionary from pandexo output. If parameter space was run in run_pandexo
        make sure to restructure the input as a list of dictionaries without they key words
        that run_pandexo assigns.
    plot : bool
        (Optional) True renders plot, Flase does not. Default=True
    output_file : str
        (Optional) Default = 'snr.html'

    Return
    ------
    x : numpy array
        micron
    y : numpy array
        1D SNR

    See Also
    --------
    jwst_1d_bkg, jwst_noise, jwst_1d_flux, jwst_1d_spec, jwst_2d_det, jwst_2d_sat
    """
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save"
    # Flux 1d
    x= result_dict['RawData']['wave']
    electrons_out = result_dict['RawData']['electrons_out']
    y = electrons_out/np.sqrt(result_dict['RawData']['var_out'])
    x = x[~np.isnan(y)]
    y = y[~np.isnan(y)]
    plot_snr_1d1 = Figure(tools=TOOLS,
                         x_axis_label='Wavelength (micron)',
                         y_axis_label='SNR', title="SNR Out of Trans",
                         width=800, height=300)
    plot_snr_1d1.line(x, y, line_width = 4, alpha = .7)
    if plot:
        outputfile(output_file)
        show(plot_snr_1d1)
    return x,y

def jwst_1d_bkg(result_dict, plot=True, output_file='bkg.html'):
    """Plot background

    Parameters
    ----------
    result_dict : dict
        Dictionary from pandexo output. If parameter space was run in run_pandexo
        make sure to restructure the input as a list of dictionaries without they key words
        that run_pandexo assigns.
    plot : bool
        (Optional) True renders plot, Flase does not. Default=True
    output_file : str
        (Optional) Default = bkt.html

    Return
    ------
    x : numpy array
        micron
    y : numpy array
        1D bakground e/s

    See Also
    --------
    jwst_1d_spec, jwst_noise, jwst_1d_flux, jwst_1d_snr, jwst_2d_det, jwst_2d_sat
    """
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save"
    # BG 1d
    out = result_dict['PandeiaOutTrans']
    x, y = out['1d']['extracted_bg_only']
    y = y[~np.isnan(y)]
    x = x[~np.isnan(y)]
    plot_bg_1d1 = Figure(tools=TOOLS,
                         x_axis_label='Wavelength [microns]',
                         y_axis_label='Flux (e/s)', title="Background",
                         width=800, height=300)
    plot_bg_1d1.line(x, y, line_width = 4, alpha = .7)
    if plot:
        outputfile(output_file)
        show(plot_bg_1d1)
    return x,y

def jwst_noise(result_dict, plot=True, output_file= 'noise.html'):
    """Plot background

    Parameters
    ----------
    result_dict : dict
        Dictionary from pandexo output. If parameter space was run in run_pandexo
        make sure to restructure the input as a list of dictionaries without they key words
        that run_pandexo assigns.
    plot : bool
        (Optional) True renders plot, Flase does not. Default=True
    output_file : str
        (Optional) Default = 'noise.html'
    Return
    ------
    x : numpy array
        micron
    y : numpy array
        1D noise (ppm)

    See Also
    --------
    jwst_1d_spec, jwst_1d_bkg, jwst_1d_flux, jwst_1d_snr, jwst_2d_det, jwst_2d_sat
    """
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save"    #saturation

    x = result_dict['FinalSpectrum']['wave']
    y = result_dict['FinalSpectrum']['error_w_floor']*1e6
    x = x[~np.isnan(y)]
    y = y[~np.isnan(y)]
    ymed = np.median(y)


    plot_noise_1d1 = Figure(tools=TOOLS,#responsive=True,
                         x_axis_label='Wavelength (micron)',
                         y_axis_label='Error on Spectrum (PPM)', title="Error Curve",
                         width=800, height=300, y_range = [0,2.0*ymed])
    ymed = np.median(y)
    plot_noise_1d1.circle(x, y, line_width = 4, alpha = .7)
    if plot:
        outputfile(output_file)
        show(plot_noise_1d1)
    return x,y

def jwst_2d_det(result_dict, plot=True, output_file='det2d.html'):
    """Plot 2d detector image

    Parameters
    ----------
    result_dict : dict
        Dictionary from pandexo output. If parameter space was run in run_pandexo
        make sure to restructure the input as a list of dictionaries without they key words
        that run_pandexo assigns.

    plot : bool
        (Optional) True renders plot, Flase does not. Default=True
    output_file : str
        (Optional) Default = 'det2d.html'

    Return
    ------
    numpy array
        2D array of out of transit detector simulation
    See Also
    --------
    jwst_1d_spec, jwst_1d_bkg, jwst_1d_flux, jwst_1d_snr, jwst_noise, jwst_2d_sat

    """
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save"
    out = result_dict['PandeiaOutTrans']
    data = out['2d']['detector']


    xr, yr = data.shape

    plot_detector_2d = Figure(tools="pan,wheel_zoom,box_zoom,reset,hover,save",
                         x_range=[0, yr], y_range=[0, xr],
                         x_axis_label='Pixel', y_axis_label='Spatial',
                         title="2D Detector Image",
                        width=800, height=300)

    plot_detector_2d.image(image=[data], x=[0], y=[0], dh=[xr], dw=[yr],
                      palette="Spectral11")
    if plot:
        outputfile(output_file)
        show(plot_detector_2d)
    return data


def jwst_2d_sat(result_dict, plot=True, output_file='sat2d.html'):
    """Plot 2d saturation profile

    Parameters
    ----------
    result_dict : dict
        Dictionary from pandexo output. If parameter space was run in run_pandexo
        make sure to restructure the input as a list of dictionaries without they key words
        that run_pandexo assigns.

    plot : bool
        (Optional) True renders plot, Flase does not. Default=True
    output_file : str
        (Optional) Default = 'sat2d.html'

    Return
    ------
    numpy array
        2D array of out of transit detector simulation

    See Also
    --------
    jwst_1d_spec, jwst_1d_bkg, jwst_1d_flux, jwst_1d_snr, jwst_2d_det, jwst_noise
    """
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save"    #saturation
    out = result_dict['PandeiaOutTrans']
    data = out['2d']['saturation']
    xr, yr = data.shape
    plot_sat_2d = Figure(tools=TOOLS,
                         x_range=[0, yr], y_range=[0, xr],
                         x_axis_label='Pixel', y_axis_label='Spatial',
                         title="Saturation",
                        width=800, height=300)

    plot_sat_2d.image(image=[data], x=[0], y=[0], dh=[xr], dw=[yr],
                      palette="Spectral11")
    if plot:
        outputfile(output_file)
        show(plot_sat_2d)
    return data

def hst_spec(result_dict, plot=True, output_file ='hstspec.html', model = True, output_notebook=True):
    """Plot 1d spec with error bars for hst

    Parameters
    ----------
    result_dict : dict
        Dictionary from pandexo output.

    plot : bool
        (Optional) True renders plot, False does not. Default=True
    model : bool
        (Optional) Plot model under data. Default=True
    output_file : str
        (Optional) Default = 'hstspec.html'
    output_notebook : bool 
        (Optional) Default true, plots in notebook

    Return
    ------
    x : numpy array
        micron
    y : numpy array
        1D spec fp/f* or rp^2/r*^2
    e : numpy array
        1D rms noise
    modelx : numpy array
        micron
    modely : numpy array
        1D spec fp/f* or rp^2/r*^2
    See Also
    --------
    hst_time
    """
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save"
    #plot planet spectrum
    mwave = result_dict['planet_spec']['model_wave']
    mspec = result_dict['planet_spec']['model_spec']

    binwave = result_dict['planet_spec']['binwave']
    binspec = result_dict['planet_spec']['binspec']

    error = result_dict['planet_spec']['error']
    error = np.zeros(len(binspec))+ error
    xlims = [result_dict['planet_spec']['wmin'], result_dict['planet_spec']['wmax']]
    ylims = [np.min(binspec)-2.0*error[0], np.max(binspec)+2.0*error[0]]

    plot_spectrum = Figure(width=800, height=300, x_range=xlims,
                               y_range=ylims, tools=TOOLS,#responsive=True,
                                 x_axis_label='Wavelength [microns]',
                                 y_axis_label='Ratio',
                               title="Original Model with Observation")

    y_err = []
    x_err = []
    for px, py, yerr in zip(binwave, binspec, error):
        np.array(x_err.append((px, px)))
        np.array(y_err.append((py - yerr, py + yerr)))
    if model:
        plot_spectrum.line(mwave,mspec, color= "black", alpha = 0.5, line_width = 4)
    plot_spectrum.circle(binwave,binspec, line_width=3, line_alpha=0.6)
    plot_spectrum.multi_line(x_err, y_err)

    if output_notebook & plot:
        outnotebook()
        show(plot_spectrum)
    elif plot:
        outputfile(output_file)
        
        show(plot_spectrum)

    return binwave, binspec, error, mwave, mspec

def hst_time(result_dict, plot=True, output_file ='hsttime.html', model = True, output_notebook=True):
    """Plot earliest and latest start times for hst observation

    Parameters
    ----------
    result_dict : dict
        Dictionary from pandexo output.

    plot : bool
        (Optional) True renders plot, False does not. Default=True
    model : bool
        (Optional) Plot model under data. Default=True
    output_file : str
        (Optional) Default = 'hsttime.html'

    Return
    ------
    obsphase1 : numpy array
        earliest start time
    obstr1 : numpy array
        white light curve
    obsphase2 : numpy array
        latest start time
    obstr2 : numpy array
        white light curve
    rms : numpy array
        1D rms noise

    See Also
    --------
    hst_spec
    """
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save"
    #earliest and latest start times
    obsphase1 = result_dict['calc_start_window']['obsphase1']
    obstr1 = result_dict['calc_start_window']['obstr1']
    rms = result_dict['calc_start_window']['light_curve_rms']
    obsphase2 = result_dict['calc_start_window']['obsphase2']
    obstr2 = result_dict['calc_start_window']['obstr2']
    phase1 = result_dict['calc_start_window']['phase1']
    phase2 = result_dict['calc_start_window']['phase2']
    trmodel1 = result_dict['calc_start_window']['trmodel1']
    trmodel2 = result_dict['calc_start_window']['trmodel2']

    if isinstance(rms, float):
        rms = np.zeros(len(obsphase1))+rms
    y_err1 = []
    x_err1 = []
    for px, py, yerr in zip(obsphase1, obstr1, rms):
        np.array(x_err1.append((px, px)))
        np.array(y_err1.append((py - yerr, py + yerr)))

    y_err2 = []
    x_err2 = []
    for px, py, yerr in zip(obsphase2, obstr2, rms):
        np.array(x_err2.append((px, px)))
        np.array(y_err2.append((py - yerr, py + yerr)))

    early = Figure(width=400, height=300,
                               tools=TOOLS,#responsive=True,
                                 x_axis_label='Orbital Phase',
                                 y_axis_label='Flux',
                               title="Earliest Start Time")

    if model: early.line(phase1, trmodel1, color='black',alpha=0.5, line_width = 4)
    early.circle(obsphase1, obstr1, line_width=3, line_alpha=0.6)
    early.multi_line(x_err1, y_err1)

    late = Figure(width=400, height=300,
                                tools=TOOLS,#responsive=True,
                                 x_axis_label='Orbital Phase',
                                 y_axis_label='Flux',
                               title="Latest Start Time")
    if model: late.line(phase2, trmodel2, color='black',alpha=0.5, line_width = 3)
    late.circle(obsphase2, obstr2, line_width=3, line_alpha=0.6)
    late.multi_line(x_err2, y_err2)

    start_time = row(early, late)


    if output_notebook & plot:
        outnotebook()
        show(start_time)
    elif plot:
        outputfile(output_file)
        show(start_time)


    return obsphase1, obstr1, obsphase2, obstr2,rms


def hst_simulated_lightcurve(result_dict, plot=True, output_file ='hsttime.html', model = True, output_notebook=True):
    """Plot simulated HST light curves (in fluece) for earliest and latest start times

    Parameters
    ----------
    result_dict : dict
        Dictionary from pandexo output.

    plot : bool
        (Optional) True renders plot, False does not. Default=True
    model : bool
        (Optional) Plot model under data. Default=True
    output_file : str
        (Optional) Default = 'hsttime.html'

    Return
    ------
    obsphase1 : numpy array
        earliest start time
    counts1 : numpy array
        white light curve in fluence (e/pixel)
    obsphase2 : numpy array
        latest start time
    counts2 : numpy array
        white light curve in fluence (e/pixel)
    rms : numpy array
        1D rms noise

    See Also
    --------
    hst_spec
    """
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save"
    # earliest and latest start times
    obsphase1 = result_dict['light_curve']['obsphase1']
    rms = result_dict['light_curve']['light_curve_rms']
    obsphase2 = result_dict['light_curve']['obsphase2']
    phase1 = result_dict['light_curve']['phase1']
    phase2 = result_dict['light_curve']['phase2']
    counts1 = result_dict['light_curve']['counts1']
    counts2 = result_dict['light_curve']['counts2']
    count_noise = result_dict['light_curve']['count_noise']
    ramp_included = result_dict['light_curve']['ramp_included']
    model_counts1 = result_dict['light_curve']['model_counts1']
    model_counts2 = result_dict['light_curve']['model_counts2']

    if isinstance(count_noise, float):
        rms = np.zeros(len(counts1)) + count_noise
    y_err1 = []
    x_err1 = []
    for px, py, yerr in zip(obsphase1, counts1, rms):
        np.array(x_err1.append((px, px)))
        np.array(y_err1.append((py - yerr, py + yerr)))

    y_err2 = []
    x_err2 = []
    for px, py, yerr in zip(obsphase2, counts2, rms):
        np.array(x_err2.append((px, px)))
        np.array(y_err2.append((py - yerr, py + yerr)))

    if ramp_included:
        title_description = " (Ramp Included)"
    else:
        title_description =" (Ramp Removed)"

    early = Figure(width=400, height=300,
                               tools=TOOLS,#responsive=True,
                                 x_axis_label='Orbital Phase',
                                 y_axis_label='Flux [electrons/pixel]',
                               title="Earliest Start Time" + title_description)

    if model:
        early.line(phase1, model_counts1, color='black', alpha=0.5, line_width=4)
    early.circle(obsphase1, counts1, line_width=3, line_alpha=0.6)
    early.multi_line(x_err1, y_err1)

    late = Figure(width=400, height=300,
                  tools=TOOLS,  # responsive=True,
                  x_axis_label='Orbital Phase',
                  y_axis_label='Flux [electrons/pixel]',
                  title="Latest Start Time" + title_description)
    if model:
        late.line(phase2, model_counts2, color='black', alpha=0.5, line_width=3)
    late.circle(obsphase2, counts2, line_width=3, line_alpha=0.6)
    late.multi_line(x_err2, y_err2)

    start_time = row(early, late)

    if plot:
        if output_notebook: 
            outnotebook()
        else:
            outputfile(output_file)
        show(start_time)

    return obsphase1, counts1, obsphase2, counts2, rms
