from bokeh.layouts import column,row
import numpy as np
from bokeh.plotting import output_file, show
from bokeh.plotting import figure as Figure

from bokeh.embed import components
from bokeh.models import Tabs
from bokeh.models import TabPanel as Panel
from bokeh.models import  ColumnDataSource, Slider,Select
from bokeh.models.callbacks import CustomJS
from bokeh.io import curdoc
from bokeh.layouts import row

def create_component_jwst(result_dict):
    """Generate front end plots JWST
    
    Function that is responsible for generating the front-end interactive plots for JWST.

    Parameters 
    ----------
    result_dict : dict 
        the dictionary returned from a PandExo run
    
    Returns
    -------
    tuple 
        A tuple containing `(script, div)`, where the `script` is the
        front-end javascript required, and `div` is a dictionary of plot
        objects.
    """  
    timing = result_dict['timing']
    noccultations = result_dict['timing']['Number of Transits']
    out = result_dict['PandeiaOutTrans']
    
    # select the tools we want
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save"

    #Define units for x and y axis
    punit = result_dict['input']['Primary/Secondary']
    p=1.0
    if punit == 'fp/f*': p = -1.0
    else: punit = punit
    
    if result_dict['input']['Calculation Type'] =='phase_spec':
        x_axis_label='Time (secs)'
        frac = 1.0
    else:
        x_axis_label='Wavelength [microns]'
        frac = result_dict['timing']['Num Integrations Out of Transit']/result_dict['timing']['Num Integrations In Transit']

    electrons_out = result_dict['RawData']['electrons_out']
    electrons_in = result_dict['RawData']['electrons_in']
    
    var_in = result_dict['RawData']['var_in']
    var_out = result_dict['RawData']['var_out']
    
    
    x = result_dict['FinalSpectrum']['wave']
    y = result_dict['FinalSpectrum']['spectrum_w_rand']
    err = result_dict['FinalSpectrum']['error_w_floor']

    y_err = []
    x_err = []
    for px, py, yerr in zip(x, y, err):
        np.array(x_err.append((px, px)))
        np.array(y_err.append((py - yerr, py + yerr)))

    source = ColumnDataSource(data=dict(x=x, y=y, y_err=y_err, x_err=x_err, err=err, 
                                electrons_out=electrons_out, electrons_in=electrons_in, var_in=var_in, var_out=var_out, 
                                p=var_in*0+p,nocc=var_in*0+noccultations, frac = var_in*0+frac))
    original = ColumnDataSource(data=dict(x=x, y=y, y_err=y_err, x_err=x_err, err=err, electrons_out=electrons_out, electrons_in=electrons_in, var_in=var_in, var_out=var_out))

    ylims = [min(result_dict['OriginalInput']['model_spec'])- 0.1*min(result_dict['OriginalInput']['model_spec']),
                 0.1*max(result_dict['OriginalInput']['model_spec'])+max(result_dict['OriginalInput']['model_spec'])]
    xlims = [min(result_dict['FinalSpectrum']['wave']), max(result_dict['FinalSpectrum']['wave'])]

    plot_spectrum = Figure(width=800, height=300, x_range=xlims,
                               y_range=ylims, tools=TOOLS,#responsive=True,
                                 x_axis_label=x_axis_label,
                                 y_axis_label=punit, 
                               title="Original Model with Observation")
    
    plot_spectrum.line(result_dict['OriginalInput']['model_wave'],result_dict['OriginalInput']['model_spec'], color= "black", alpha = 0.5, line_width = 4)
        
    plot_spectrum.circle('x', 'y', source=source, line_width=3, line_alpha=0.6)
    plot_spectrum.multi_line('x_err', 'y_err', source=source)

    callback = CustomJS(args=dict(source=source, original=original), code="""
            // Grab some references to the data
            var sdata = source.data;
            var odata = original.data;

            // Create copies of the original data, store them as the source data
            sdata['x'] = odata['x'].slice(0);
            sdata['y'] = odata['y'].slice(0);

            sdata['y_err'] = odata['y_err'].slice(0);
            sdata['x_err'] = odata['x_err'].slice(0);
            sdata['err'] = odata['err'].slice(0);

            sdata['electrons_out'] = odata['electrons_out'].slice(0);
            sdata['electrons_in'] = odata['electrons_in'].slice(0);
            sdata['var_in'] = odata['var_in'].slice(0);
            sdata['var_out'] = odata['var_out'].slice(0);

            // Create some variables referencing the source data
            var x = sdata['x'];
            var y = sdata['y'];
            var y_err = sdata['y_err'];
            var x_err = sdata['x_err'];
            var err = sdata['err'];
            var p = sdata['p'];
            var frac = sdata['frac'];
            var og_ntran = sdata['nocc'];

            var electrons_out = sdata['electrons_out'];
            var electrons_in = sdata['electrons_in'];
            var var_in = sdata['var_in'];
            var var_out = sdata['var_out'];

            var wbin = wbin.value;
            var ntran = ntran.value;

            var ind = [];
            ind.push(0);
            var start = 1;

            var i=1;
            for (i = 0; i < x.length-1; i++) {
                if (start == wbin) {
                    ind.push(i+1);
                    start =0;
                }
                start = start+1;
            }

            if (ind[ind.length-1] != x.length) {
                ind.push(x.length);
            }

            var xout = [];


            var foutout = [];
            var finout = [];
            var varinout = [];
            var varoutout = [];

            var xslice = []; 

            var foutslice = [];
            var finslice = [];
            var varoutslice = [];
            var varinslice = [];

            function add(a, b) {
                return a+b;
            }

            for (i = 0; i < ind.length-1; i++) {
                xslice = x.slice(ind[i],ind[i+1]);

                foutslice = electrons_out.slice(ind[i],ind[i+1]);
                finslice = electrons_in.slice(ind[i],ind[i+1]);
                
                varinslice = var_in.slice(ind[i],ind[i+1]);
                varoutslice = var_out.slice(ind[i],ind[i+1]);

                xout.push(xslice.reduce(add, 0)/xslice.length);
                foutout.push(foutslice.reduce(add, 0));
                finout.push(finslice.reduce(add, 0));
                
                varinout.push(varinslice.reduce(add, 0));
                varoutout.push(varoutslice.reduce(add, 0));

                xslice = [];
                foutslice = [];
                finslice = [];
                varinslice = [];
                varoutslice = [];
            }
            
            var new_err = 1.0;
            var rand = 1.0;

            for (i = 0; i < x.length; i++) {
                new_err = Math.pow((frac[i]/foutout[i]),2)*varinout[i] + Math.pow((finout[i]*frac[i]/Math.pow(foutout[i],2)),2)*varoutout[i];
                new_err = Math.sqrt(new_err)*Math.sqrt(og_ntran[i]/ntran);
                rand = new_err*(Math.random()-Math.random());
                y[i] = p[i]*((1.0 - frac[i]*finout[i]/foutout[i]) + rand); 
                x[i] = xout[i];
                x_err[i][0] = xout[i];
                x_err[i][1] = xout[i];
                y_err[i][0] = y[i] + new_err;
                y_err[i][1] = y[i] - new_err;            
            }

            source.change.emit();
        """)

    #var_tot = (frac/electrons_out)**2.0 * var_in + (electrons_in*frac/electrons_out**2.0)**2.0 * var_out

    sliderWbin =  Slider(title="# pixels to bin", value=1, start=1, end=50, step= 1)
    sliderWbin.js_on_change('value', callback)
    callback.args["wbin"] = sliderWbin
    
    sliderTrans =  Slider(title="Num Trans", value=noccultations, start=1, end=50, step= 1)
    sliderTrans.js_on_change('value', callback)
    callback.args["ntran"] = sliderTrans

    layout = column(row(sliderWbin,sliderTrans), plot_spectrum)


    #out of transit 2d output 
    raw = result_dict['RawData']
    
    # Flux 1d
    x, y = raw['wave'], raw['e_rate_out']*result_dict['timing']['Seconds per Frame']*(timing["APT: Num Groups per Integration"]-1)
    x = x[~np.isnan(y)]
    y = y[~np.isnan(y)]

    plot_flux_1d1 = Figure(tools=TOOLS,
                         x_axis_label='Wavelength [microns]',
                         y_axis_label='e-/integration', title="Flux Per Integration",
                         width=800, height=300)
    plot_flux_1d1.line(x, y, line_width = 4, alpha = .7)
    tab1 = Panel(child=plot_flux_1d1, title="Flux per Int")

    # BG 1d
    #x, y = out['1d']['extracted_bg_only']
    #y = y[~np.isnan(y)]
    #x = x[~np.isnan(y)]
    #plot_bg_1d1 = Figure(tools=TOOLS,
    #                     x_axis_label='Wavelength [microns]',
    #                     y_axis_label='Flux (e/s)', title="Background",
    #                     width=800, height=300)
    #plot_bg_1d1.line(x, y, line_width = 4, alpha = .7)
    #tab2 = Panel(child=plot_bg_1d1, title="Background Flux")

    # SNR 
    y = np.sqrt(y) #this is computing the SNR (sqrt of photons in a single integration)


    plot_snr_1d1 = Figure(tools=TOOLS,
                         x_axis_label=x_axis_label,
                         y_axis_label='sqrt(e-)/integration', title="SNR per integration",
                         width=800, height=300)
    plot_snr_1d1.line(x, y, line_width = 4, alpha = .7)
    tab3 = Panel(child=plot_snr_1d1, title="SNR per Int")


    # Error bars (ppm) 
    x = result_dict['FinalSpectrum']['wave']
    y = result_dict['FinalSpectrum']['error_w_floor']*1e6
    x = x[~np.isnan(y)]
    y = y[~np.isnan(y)]    
    ymed = np.median(y)

    source2 = ColumnDataSource(data=dict(x=x, y=y,nocc=x*0+noccultations))
    original2 = ColumnDataSource(data=dict(x=x, y=y,nocc=x*0+noccultations))

    plot_noise_1d1 = Figure(tools=TOOLS,#responsive=True,
                         x_axis_label=x_axis_label,
                         y_axis_label='Spectral Precision (ppm)', title="Spectral Precision",
                         width=800, height=300, y_range = [0,2.0*ymed])
    ymed = np.median(y)
    plot_noise_1d1.circle('x', 'y', line_width = 4, alpha = .7, source=source2)

    callback2 = CustomJS(args=dict(source=source2, original=original2), code="""
            // Grab some references to the data
            var sdata = source.data;
            var odata = original.data;

            // Create copies of the original data, store them as the source data
            sdata['x'] = odata['x'].slice(0);
            sdata['y'] = odata['y'].slice(0);

            // Create some variables referencing the source data
            var x = sdata['x'];
            var y = sdata['y'];

            var og_ntran = sdata['nocc'];

            var wbin = wbin.value;
            var ntran = ntran.value;

            var ind = [];
            ind.push(0);
            var start = 1;

            var i=1;
            for (i = 0; i < x.length-1; i++) {
                if (start == wbin) {
                    ind.push(i+1);
                    start =0;
                }
                start = start+1;
            }

            if (ind[ind.length-1] != x.length) {
                ind.push(x.length);
            }

            var xout = [];
            var yout = [];

            var xslice = []; 
            var yslice = [];

            function addpow(tot, num) {
                return tot+ Math.pow(num,2);
            }
            function add(a, b) {
                return a+b;
            }
            for (i = 0; i < ind.length-1; i++) {

                xslice = x.slice(ind[i],ind[i+1]);
                yslice = y.slice(ind[i],ind[i+1]);

                xout.push(xslice.reduce(add, 0)/xslice.length);
                yout.push(Math.sqrt(yslice.reduce(addpow, 0))/yslice.length);

                xslice = [];
                yslice = [];
            }
            
            var new_err = 1.0;
            for (i = 0; i < x.length; i++) {
                new_err = yout[i]*Math.sqrt(og_ntran[i]/ntran);
                x[i] = xout[i];
                y[i] = new_err       
            }

            source.change.emit();
        """)    

    sliderWbin2 =  Slider(title="# pixels to bin", value=1, start=1, end=50, step= 1)
    sliderWbin2.js_on_change('value', callback2)
    callback2.args["wbin"] = sliderWbin2  

    sliderTrans2 =  Slider(title="Num Trans", value=noccultations, start=1, end=50, step= 1)
    sliderTrans2.js_on_change('value', callback2)
    callback2.args["ntran"] = sliderTrans2
    
    noise_lay = column(row(sliderWbin2, sliderTrans2) , plot_noise_1d1)
    tab4 = Panel(child=noise_lay, title="Precision")

    #Not happy? Need help picking a different mode? 
    plot_spectrum2 = Figure(width=800, height=300, x_range=xlims,y_range=ylims, tools=TOOLS,
                             x_axis_label=x_axis_label,
                             y_axis_label=punit, title="Original Model",y_axis_type="log")

    plot_spectrum2.line(result_dict['OriginalInput']['model_wave'],result_dict['OriginalInput']['model_spec'],
                        line_width = 4,alpha = .7)
    tab5 = Panel(child=plot_spectrum2, title="Original Model")


    #create set of five tabs 
    tabs1d = Tabs(tabs=[ tab1,tab3, tab4, tab5])



    # Detector 2d
    data = out['2d']['detector']

    
    xr, yr = data.shape
    
    plot_detector_2d = Figure(tools="pan,wheel_zoom,box_zoom,reset,hover,save",
                         x_range=[0, yr], y_range=[0, xr],
                         x_axis_label='Pixel', y_axis_label='Spatial',
                         title="2D Detector Image",
                        width=800, height=300)
    
    plot_detector_2d.image(image=[data], x=[0], y=[0], dh=[xr], dw=[yr],
                      palette="Spectral11")


    #2d tabs 

    #2d snr 
    data = out['2d']['snr']
    data[np.isinf(data)] = 0.0
    xr, yr = data.shape
    plot_snr_2d = Figure(tools=TOOLS,
                         x_range=[0, yr], y_range=[0, xr],
                         x_axis_label='Pixel', y_axis_label='Spatial',
                         title="Signal-to-Noise Ratio",
                        width=800, height=300)
    
    plot_snr_2d.image(image=[data], x=[0], y=[0], dh=[xr], dw=[yr],
                      palette="Spectral11")
    
    tab1b = Panel(child=plot_snr_2d, title="SNR")

    #saturation
    
    data = out['2d']['saturation']
    xr, yr = data.shape
    plot_sat_2d = Figure(tools=TOOLS,
                         x_range=[0, yr], y_range=[0, xr],
                         x_axis_label='Pixel', y_axis_label='Spatial',
                         title="Saturation",
                        width=800, height=300)
    
    plot_sat_2d.image(image=[data], x=[0], y=[0], dh=[xr], dw=[yr],
                      palette="Spectral11")
    
    tab2b = Panel(child=plot_sat_2d, title="Saturation")

    tabs2d = Tabs(tabs=[ tab1b, tab2b])
    
 
    result_comp = components({'plot_spectrum':layout, 
                              'tabs1d': tabs1d, 'det_2d': plot_detector_2d,
                              'tabs2d': tabs2d})

    return result_comp
    


def create_component_hst(result_dict):
    """Generate front end plots HST
    
    Function that is responsible for generating the front-end spectra plots for HST.
    
    Parameters
    ----------
    result_dict : dict 
        The dictionary returned from a PandExo (HST) run

    Returns
    -------
    tuple
        A tuple containing `(script, div)`, where the `script` is the
        front-end javascript required, and `div` is a dictionary of plot
        objects.
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
    
    eventType = result_dict['calc_start_window']['eventType'] 

    if eventType=='transit':
        y_axis = '(Rp/R*)^2'
    elif eventType =='eclipse':
        y_axis='Fp/F*'

    plot_spectrum = Figure(width=800, height=300, x_range=xlims,
                               y_range=ylims, tools=TOOLS,#responsive=True,
                                 x_axis_label='Wavelength [microns]',
                                 y_axis_label=y_axis, 
                               title="Original Model with Observation")
    
    y_err = []
    x_err = []
    for px, py, yerr in zip(binwave, binspec, error):
        np.array(x_err.append((px, px)))
        np.array(y_err.append((py - yerr, py + yerr)))

    plot_spectrum.line(mwave,mspec, color= "black", alpha = 0.5, line_width = 4)
    plot_spectrum.circle(binwave,binspec, line_width=3, line_alpha=0.6)
    plot_spectrum.multi_line(x_err, y_err)
    
    
    #earliest and latest start times without ramp
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
    
    early.line(phase1, trmodel1, color='black',alpha=0.5, line_width = 4)
    early.circle(obsphase1, obstr1, line_width=3, line_alpha=0.6)
    early.multi_line(x_err1, y_err1)
     
    late = Figure(width=400, height=300, 
                                tools=TOOLS,#responsive=True,
                                 x_axis_label='Orbital Phase',
                                 y_axis_label='Flux', 
                               title="Latest Start Time")
    late.line(phase2, trmodel2, color='black',alpha=0.5, line_width = 3)
    late.circle(obsphase2, obstr2, line_width=3, line_alpha=0.6)
    late.multi_line(x_err2, y_err2)
        
    start_time = row(early, late)
    tab1b = Panel(child=start_time, title="Timing")
    
    #now create plot that has optional ramp
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


    early.line(phase1, model_counts1, color='black', alpha=0.5, line_width=4)
    early.circle(obsphase1, counts1, line_width=3, line_alpha=0.6)
    early.multi_line(x_err1, y_err1)

    late = Figure(width=400, height=300,
                  tools=TOOLS,  # responsive=True,
                  x_axis_label='Orbital Phase',
                  y_axis_label='Flux [electrons/pixel]',
                  title="Latest Start Time" + title_description)

    late.line(phase2, model_counts2, color='black', alpha=0.5, line_width=3)
    late.circle(obsphase2, counts2, line_width=3, line_alpha=0.6)
    late.multi_line(x_err2, y_err2)

    start_time = row(early, late)


    tab2b = Panel(child=start_time, title="Electrons")

    tabs2d = Tabs(tabs=[ tab1b, tab2b])
    
    
    result_comp = components({'plot_spectrum':plot_spectrum, 
                              'start_time':tabs2d})

    return result_comp


 
