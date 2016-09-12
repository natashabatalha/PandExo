import numpy as np


def loopingL(cen, signal_col, noise_col, bkg_col):
#create function to find location where SNR is the highest
#loop from the highest value of the signal downward 
    sn_old= 0
    for ii in range(0,cen+1): #0,1,2,3,4,5 range(0,6).. center=5
        sig_sum = sum(signal_col[cen-ii:cen+1])
        noi_sum = np.sqrt(sum(noise_col[cen-ii:cen+1]))
        bkg_sum = sum(bkg_col[cen-ii:cen+1])
        sn_new = sig_sum/noi_sum
        if sn_old >=sn_new:               
            return cen-ii+1                     
        else:
            sn_old = sn_new
    return 0


def loopingU(cen, signal_col, noise_col, bkg_col):
#create function to find location where SNR is the highest
#loop from the highest value of the signal upward
    sn_old= 0
    for ii in range(1,(len(signal_col)-cen+1)): #1,2,3,4,5,6 range(1,6).. center=5, edge =10
        sig_sum = sum(signal_col[cen:cen+ii])
        noi_sum = np.sqrt(sum(noise_col[cen:cen+ii]))
        bkg_sum = sum(bkg_col[cen:cen+ii])
        sn_new = sig_sum/noi_sum
           
        if sn_old >=sn_new:               
            return cen+ii-1                     
        else:
            sn_old = sn_new
    return len(signal_col)+1 


def sum_spatial(extract_info, nocc, nint_in, nint_out):    #last 
    #extract info 
    LBout = extract_info['bounds']['LBout']
    LBin = extract_info['bounds']['LBin']
    UBin = extract_info['bounds']['UBin']
    UBout = extract_info['bounds']['UBout']
    photon_sig_in = extract_info['photons']['photon_sig_in']    
    photon_sig_out = extract_info['photons']['photon_sig_out']
    var_pix_in = extract_info['photons']['var_pix_in']
    var_pix_out = extract_info['photons']['var_pix_out']
    photon_sky_in = extract_info['noise']['photon_sky_in']
    photon_sky_out = extract_info['noise']['photon_sky_out']
    rn_var_in = extract_info['noise']['rn_var_in']
    rn_var_out = extract_info['noise']['rn_var_out']
    

    rr, lenw = photon_sig_in.shape

    # sum spectrum in the spatial direction to create 1d spectrum 
    photon_out_1d = range(0,lenw)
    photon_in_1d = range(0,lenw)
    var_out_1d = range(0,lenw)
    var_in_1d = range(0,lenw)
    sky_in_1d = range(0,lenw)
    sky_out_1d = range(0,lenw)
    rn_in_1d = range(0,lenw)
    rn_out_1d = range(0,lenw)
    

    #sum 2d spectrum in extraction region in spatial direciton to create 1d spec
    for i in range(0,lenw):
        photon_out_1d[i] = sum(photon_sig_out[LBout[i]:UBout[i],i])*nint_out*nocc
        photon_in_1d[i] = sum(photon_sig_in[LBout[i]:UBout[i],i])*nint_in*nocc
        var_out_1d[i] = sum(var_pix_out[LBout[i]:UBout[i],i])*nint_out*nocc
        var_in_1d[i] = sum(var_pix_in[LBout[i]:UBout[i],i])*nint_in*nocc
        sky_out_1d[i] = sum(photon_sky_out[LBout[i]:UBout[i],i])*nint_out*nocc
        sky_in_1d[i] = sum(photon_sky_in[LBout[i]:UBout[i],i])*nint_in*nocc
        rn_out_1d[i] = sum(rn_var_out[LBout[i]:UBout[i],i])*nint_out*nocc
        rn_in_1d[i] = sum(rn_var_in[LBout[i]:UBout[i],i])*nint_in*nocc

        
    photon_out_1d = np.array(photon_out_1d)
    photon_in_1d = np.array(photon_in_1d)
    var_out_1d = np.array(var_out_1d)
    var_in_1d = np.array(var_in_1d)
    sky_in_1d = np.array(sky_in_1d)
    sky_out_1d = np.array(sky_out_1d)
    rn_out_1d = np.array(rn_out_1d)
    rn_in_1d = np.array(rn_in_1d)

    return {'photon_out_1d':photon_out_1d, 'photon_in_1d':photon_in_1d, 'var_in_1d':var_in_1d, 'var_out_1d':var_out_1d, 'rn_in_1d':rn_in_1d,'rn_out_1d':rn_out_1d,'sky_in_1d': sky_in_1d, 'sky_out_1d': sky_out_1d, 'extract_info':extract_info}

def extract_region(inn, out, exptime_per_int, ngroups_per_int): #second to last 

    #Full variance including all noise sources
    detector_var=out.noise.var_pix
    #Standard deviation (i.e. sqrt(detector_var))
    detector_stdev=out.noise.stdev_pix

    #signal (no noise no background)
    s_out = out.signals[0].rate 
    s_in = inn.signals[0].rate 

    #size of wavelength direction 
    rr,lenw = s_out.shape

    #background 
    bkgd_out = out.signals[0].rate_plus_bg - out.signals[0].rate
    bkgd_in = inn.signals[0].rate_plus_bg - inn.signals[0].rate


    #define psf center based on max value on detector
    cenRo,cenCo = np.unravel_index(s_out.argmax(), s_out.shape)
    cenRi,cenCi = np.unravel_index(s_out.argmax(), s_in.shape)

    #define out of transit parameters to calculate extraction region
    
    #multiaccum sample data factor 
    factor_flux = 6.0/5.0*(ngroups_per_int**2.0+1.0)/(ngroups_per_int**2.0+ngroups_per_int)
    factor_rn = 12.0*(ngroups_per_int-1.0)/(ngroups_per_int**2.0 + ngroups_per_int)

    photon_sig_out = s_out*exptime_per_int*factor_flux #total photons per pixel in signal
    photon_sky_out = bkgd_out*exptime_per_int*factor_flux #total photons per pixel in background 
    rn_var_out= out.noise.var_rn_pix*factor_rn #variance stricly due to detector readnoise
    var_pix_out = photon_sig_out + photon_sky_out + rn_var_out # variance of noise per pixel

    #define parameters for IN transit 
    photon_sig_in = s_in*exptime_per_int*factor_flux #total photons per pixel in signal
    photon_sky_in = bkgd_in*exptime_per_int*factor_flux #total photons per pixel in background 
    rn_var_in= inn.noise.var_rn_pix*factor_rn #variance stricly due to detector readnoise
    var_pix_in = photon_sig_in + photon_sky_in + rn_var_in # variance of noise per pixel 

    UBout = range(0,lenw)
    LBout = range(0,lenw)
    UBin = range(0,lenw)
    LBin = range(0,lenw)


    #start new loop over column
    for j in range(0,lenw):

        #define column of interest
        noise_col_out = var_pix_out[:, j]
        signal_col_out = photon_sig_out[:,j]
        bkg_col_out = photon_sky_out[:,j]
        rn_col_out = rn_var_out[:,j]
        
        noise_col_in =  var_pix_in[:, j]
        signal_col_in = photon_sig_in[:,j]
        bkg_col_in = photon_sky_in[:,j]
        
 
        #store lower and upper bound extraction regions for in and out of transit 
        LBout[j] = loopingL(cenRo, signal_col_out, noise_col_out, bkg_col_out)
        UBout[j] = loopingU(cenRo, signal_col_out, noise_col_out, bkg_col_out)
        LBin[j] = loopingL(cenRi, signal_col_in, noise_col_in, bkg_col_in    )
        UBin[j] = loopingU(cenRi, signal_col_in, noise_col_in, bkg_col_in)
        
    #this could be made more elegant later... not very efficient 
    noise = {'rn_var_out':rn_var_out, 'rn_var_in':rn_var_in, 'photon_sky_in':photon_sky_in, 'photon_sky_out': photon_sky_out }
    bounds = {'LBout':LBout, 'UBout':UBout, 'LBin':LBin, 'UBin':UBin}
    photons = {'photon_sig_out': photon_sig_out, 'photon_sig_in':photon_sig_in, 'var_pix_in':var_pix_in,'var_pix_out': var_pix_out}
    extract_info ={'bounds':bounds, 'photons':photons, 'noise':noise}
    return extract_info

