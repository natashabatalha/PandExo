from bokeh.layouts import column,row
import numpy as np
from bokeh.plotting import Figure, output_file, show
from bokeh.embed import components
from bokeh.models.widgets import Panel, Tabs
from bokeh.models import CustomJS, ColumnDataSource, Slider,Select
from bokeh.io import curdoc


def create_component(result_dict):
    """
    Responsible for generating the front-end interactive plots.

    Parameters
    ----------
    result_dict : dict
        The dictionary returned from a Pandeia run.

    Returns
    -------
    result_comp : tuple
        A tuple containing `(script, div)`, where the `script` is the
        front-end javascript required, and `div` is a dictionary of plot
        objects.
    """
    
    noccultations = result_dict['timing']['Number of Transits']
    
    # select the tools we want
    TOOLS = "pan,wheel_zoom,box_zoom,resize,reset,save"

    #MAIN spectrum 
    punit = result_dict['input']['Primary/Secondary']
    p=1.0
    if punit == 'fp/f*': p = -1.0
    else: punit = '('+punit+')^2'

    flux_out = result_dict['RawData']['flux_out']
    flux_in = result_dict['RawData']['flux_in']
    var_tot = result_dict['RawData']['var_out'] + result_dict['RawData']['var_in']

    x = result_dict['FinalSpectrum']['wave']
    y = result_dict['FinalSpectrum']['spectrum_w_rand']
    err = result_dict['FinalSpectrum']['error_w_floor']

    y_err = []
    x_err = []
    for px, py, yerr in zip(x, y, err):
        np.array(x_err.append((px, px)))
        np.array(y_err.append((py - yerr, py + yerr)))

    source = ColumnDataSource(data=dict(x=x, y=y, y_err=y_err, x_err=x_err, err=err, flux_out=flux_out, flux_in=flux_in, var_tot=var_tot, p=var_tot*0+p))
    original = ColumnDataSource(data=dict(x=x, y=y, y_err=y_err, x_err=x_err, err=err, flux_out=flux_out, flux_in=flux_in, var_tot=var_tot))

    ylims = [min(result_dict['OriginalInput']['og_spec'])- 0.1*min(result_dict['OriginalInput']['og_spec']),
                 0.1*max(result_dict['OriginalInput']['og_spec'])+max(result_dict['OriginalInput']['og_spec'])]
    xlims = [min(result_dict['FinalSpectrum']['wave']), max(result_dict['FinalSpectrum']['wave'])]

    plot_spectrum = Figure(plot_width=800, plot_height=300, x_range=xlims,
                               y_range=ylims, tools=TOOLS,#responsive=True,
                                 x_axis_label='Wavelength [microns]',
                                 y_axis_label=punit, 
                               title="Original Model with Observation")
    
    plot_spectrum.line(result_dict['OriginalInput']['og_wave'],result_dict['OriginalInput']['og_spec'], color= "black", alpha = 0.5, line_width = 4)
        
    plot_spectrum.circle('x', 'y', source=source, line_width=3, line_alpha=0.6)
    plot_spectrum.multi_line('x_err', 'y_err', source=source)

    callback = CustomJS(args=dict(source=source, original=original), code="""
            // Grab some references to the data
            var sdata = source.get('data');
            var odata = original.get('data');

            // Create copies of the original data, store them as the source data
            sdata['x'] = odata['x'].slice(0);
            sdata['y'] = odata['y'].slice(0);

            sdata['y_err'] = odata['y_err'].slice(0);
            sdata['x_err'] = odata['x_err'].slice(0);
            sdata['err'] = odata['err'].slice(0);

            sdata['flux_out'] = odata['flux_out'].slice(0);
            sdata['flux_in'] = odata['flux_in'].slice(0);
            sdata['var_tot'] = odata['var_tot'].slice(0);

            // Create some variables referencing the source data
            var x = sdata['x'];
            var y = sdata['y'];
            var y_err = sdata['y_err'];
            var x_err = sdata['x_err'];
            var err = sdata['err'];
            var p = sdata['p'];

            var flux_out = sdata['flux_out'];
            var flux_in = sdata['flux_in'];
            var var_tot = sdata['var_tot'];

            var f = wbin.get('value');
            var ntran = ntran.get('value');

            var wlength = Math.pow(10.0,f);

            var ind = [];
            ind.push(0);
            var start = 0;


            for (i = 0; i < x.length-1; i++) {
                if (x[i+1] - x[start] >= wlength) {
                    ind.push(i+1);
                    start = i;
                }
            }

            if (ind[ind.length-1] != x.length) {
                ind.push(x.length);
            }

            var xout = [];


            var foutout = [];
            var finout = [];
            var varout = [];

            var xslice = []; 


            var foutslice = [];
            var finslice = [];
            var varslice = [];

            function add(a, b) {
                return a+b;
            }

            for (i = 0; i < ind.length-1; i++) {
                xslice = x.slice(ind[i],ind[i+1]);

                foutslice = flux_out.slice(ind[i],ind[i+1]);
                finslice = flux_in.slice(ind[i],ind[i+1]);
                varslice = var_tot.slice(ind[i],ind[i+1]);

                xout.push(xslice.reduce(add, 0)/xslice.length);
                foutout.push(foutslice.reduce(add, 0));
                finout.push(finslice.reduce(add, 0));
                varout.push(varslice.reduce(add, 0));

                xslice = [];

                foutslice = [];
                finslice = [];
                varslice = [];
            }

            for (i = 0; i < x.length; i++) {
                y[i] = p[i]*(foutout[i]-finout[i]+ (Math.sqrt(varout[i]/ntran)*(Math.random()-Math.random())))/foutout[i]; 
                x[i] = xout[i];
                x_err[i][0] = xout[i];
                x_err[i][1] = xout[i];
                y_err[i][0] = y[i] + (Math.sqrt(varout[i])/foutout[i]/Math.sqrt(ntran));
                y_err[i][1] = y[i] -(Math.sqrt(varout[i])/foutout[i]/Math.sqrt(ntran));            
            }

            source.trigger('change');
        """)

    sliderWbin =  Slider(title="binning", value=np.log10(x[1]-x[0]), start=np.log10(x[1]-x[0]), end=np.log10(max(x)/2.0), step= .05, callback=callback)
    callback.args["wbin"] = sliderWbin
    sliderTrans =  Slider(title="Num Trans", value=noccultations, start=1, end=50, step= 1, callback=callback)
    callback.args["ntran"] = sliderTrans
    layout = column(row(sliderWbin,sliderTrans), plot_spectrum)


    #out of transit 2d output 
    out = result_dict['PandeiaOutTrans']
    
    # Flux 1d
    x, y = out['1d']['extracted_flux']
    x = x[~np.isnan(y)]
    y = y[~np.isnan(y)]

    plot_flux_1d1 = Figure(tools=TOOLS,
                         x_axis_label='Wavelength [microns]',
                         y_axis_label='Flux (e/s)', title="Out of Transit Flux",
                         plot_width=800, plot_height=300)
    plot_flux_1d1.line(x, y, line_width = 4, alpha = .7)
    tab1 = Panel(child=plot_flux_1d1, title="Total Flux")

    # BG 1d
    x, y = out['1d']['bg']
    y = y[~np.isnan(y)]
    x = x[~np.isnan(y)]
    plot_bg_1d1 = Figure(tools=TOOLS,
                         x_axis_label='Wavelength [microns]',
                         y_axis_label='Flux (e/s)', title="Background",
                         plot_width=800, plot_height=300)
    plot_bg_1d1.line(x, y, line_width = 4, alpha = .7)
    tab2 = Panel(child=plot_bg_1d1, title="Background Flux")

    # SNR 1d accounting for number of occultations
    x, y = out['1d']['sn']
    x = x[~np.isnan(y)]
    y = y[~np.isnan(y)]
    y = y*np.sqrt(noccultations)
    plot_snr_1d1 = Figure(tools=TOOLS,
                         x_axis_label='Wavelength [microns]',
                         y_axis_label='SNR', title="Signal-to-Noise Ratio",
                         plot_width=800, plot_height=300)
    plot_snr_1d1.line(x, y, line_width = 4, alpha = .7)
    tab3 = Panel(child=plot_snr_1d1, title="SNR")


    # Error bars (ppm) 

    x = result_dict['FinalSpectrum']['wave']
    y = result_dict['FinalSpectrum']['error_w_floor']*1e6
    x = x[~np.isnan(y)]
    y = y[~np.isnan(y)]    
    ymed = np.median(y)


    plot_noise_1d1 = Figure(tools=TOOLS,#responsive=True,
                         x_axis_label='Wavelength [microns]',
                         y_axis_label='Error on Spectrum (PPM)', title="Error Curve",
                         plot_width=800, plot_height=300, y_range = [0,2.0*ymed])
    ymed = np.median(y)
    plot_noise_1d1.circle(x, y, line_width = 4, alpha = .7)
    tab4 = Panel(child=plot_noise_1d1, title="Error")

    #Not happy? Need help picking a different mode? 
    plot_spectrum2 = Figure(plot_width=800, plot_height=300, x_range=xlims,y_range=ylims, tools=TOOLS,
                             x_axis_label='Wavelength [microns]',
                             y_axis_label=punit, title="Original Model",y_axis_type="log")

    plot_spectrum2.line(result_dict['OriginalInput']['og_wave'],result_dict['OriginalInput']['og_spec'],
                        line_width = 4,alpha = .7)
    tab5 = Panel(child=plot_spectrum2, title="Original Model")


    #create set of five tabs 
    tabs1d = Tabs(tabs=[ tab1, tab2,tab3, tab4, tab5])



    # Detector 2d
    data = out['2d']['detector']

    
    xr, yr = data.shape
    
    plot_detector_2d = Figure(tools="pan,wheel_zoom,box_zoom,resize,reset,hover,save",
                         x_range=[0, yr], y_range=[0, xr],
                         x_axis_label='Pixel', y_axis_label='Spatial',
                         title="2D Detector Image",
                        plot_width=800, plot_height=300)
    
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
                        plot_width=800, plot_height=300)
    
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
                        plot_width=800, plot_height=300)
    
    plot_sat_2d.image(image=[data], x=[0], y=[0], dh=[xr], dw=[yr],
                      palette="Spectral11")
    
    tab2b = Panel(child=plot_sat_2d, title="Saturation")

    tabs2d = Tabs(tabs=[ tab1b, tab2b])
    
 
    result_comp = components({'plot_spectrum':layout, 
                              'tabs1d': tabs1d, 'det_2d': plot_detector_2d,
                              'tabs2d': tabs2d})

    return result_comp
    
def create_component2(result_dict):
    num = -1
    color = ["red", "blue", "green", "purple", "black", "yellow", "orange", "pink","cyan","brown"]

    TOOLS = "pan,wheel_zoom,box_zoom,resize,reset,save"
        
    plot1 = Figure(plot_width=800, plot_height=350,  tools=TOOLS, responsive =False,
                                 x_axis_label='Wavelength [um]', x_axis_type="log",
                                 y_axis_label='Alpha Lambda', y_range = [min(result_dict['alpha']), max(result_dict['alpha'])]) 
    plot1.line(result_dict['w'], result_dict['alpha'], alpha = 0.5, line_width = 3)    
    
    plot2 = Figure(plot_width=800, plot_height=350,  tools=TOOLS, responsive =False, y_axis_type="log",
                                 x_axis_label='Wavelength [um]', x_axis_type="log",
                                 y_axis_label='Weighted Cross Section', y_range = [1e-29, 1e-17])
    alpha = result_dict['alpha']
    squig =  result_dict['mols']    
    xsec = result_dict['xsec'] 
    waves = result_dict['w']                       
    for i in squig.keys():  
        if i=="H2":
            num +=1
            plot2.line(waves,squig[i]*xsec[i][alpha> 0.0]*squig["H2"], color = color[num], legend = i)
        elif i=="He":
            num +=1
            plot2.line(waves,squig[i]*xsec[i][alpha> 0.0]*squig["H2"], color = color[num], legend = i)
        else:
            num +=1
            plot2.line(waves,squig[i]*xsec[i][alpha> 0.0], color = color[num], legend = i)            
                             
                          
    result_comp =   components({'plot1':plot1, 
                              'plot2': plot2})
    return result_comp 

                                 




