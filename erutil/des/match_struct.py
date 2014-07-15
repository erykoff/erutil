"""

need FITS_LDAC ...

"""

import fitsio
import esutil
import numpy as np
import sys

from erutil import header

def make_match_struct(infile, outfile, bands=['g','r','i','z']):
    """


    """
    
    inst=fitsio.read(infile,ext=1)

    new_hdr = False
    if ('HDRPATH' in inst.dtype.names) :
        new_hdr = True

    nmag=len(bands)
    nobs=inst.size
    
    
    # this is object by object
    elt = [('RA','>f8',nobs),\
               ('DEC','>f8',nobs),\
               ('EXPNUM','i4',nobs),\
               ('CCDNUM','i2',nobs),\
               ('BAND','a1',nobs),\
               ('MAG_AUTO','f4',nobs),\
               ('MAGERR_AUTO','f4',nobs),\
               ('MAG_PSF','f4',nobs),\
               ('MAGERR_PSF','f4',nobs),\
               ('NGOOD','i4',nmag),\
               ('MAG_AUTO_MEAN','f4',nmag),\
               ('MAG_AUTO_RMS','f4',nmag),\
               ('MAG_PSF_MEAN','f4',nmag),\
               ('MAG_PSF_RMS','f4',nmag),\
               ('RA_MEAN','f8'),\
               ('RA_RMS','f8'),\
               ('DEC_MEAN','f8'),\
               ('DEC_RMS','f8')]

    # this is image-by-image
    selt=[('EXPNUM','i4'),\
              ('CCDNUM','i4'),\
              ('BAND','a2'),\
              ('AIRMASS','f4'),\
              ('NGOOD','i4'),\
              ('MAG_AUTO_OFFSET','f4'),\
              ('MAG_AUTO_OFFSET_SIGMA','f4'),\
              ('MAG_PSF_OFFSET','f4'),\
              ('MAG_PSF_OFFSET_SIGMA','f4'),\
              ('RA_OFFSET','f8'),\
              ('RA_OFFSET_RMS','f8'),\
              ('DEC_OFFSET','f8'),\
              ('DEC_OFFSET_RMS','f8')]

    stats=np.recarray((nobs,),dtype=selt)    
    
    for i in range(inst.size):
        sys.stdout.write('\b'*20)
        sys.stdout.write('%d of %d' % (i, inst.size))
        sys.stdout.flush()

        if (new_hdr):
            hdr = header.read_astromatic_header(inst[i]['PATH'], inst[i]['HDRPATH'])
        else :
            hdr = header.read_astromatic_header(inst[i]['PATH'], None)

        filt=hdr['FILTER']
        band=filt[0]

        stats[i]['BAND'] = band
        stats[i]['EXPNUM'] = hdr['EXPNUM']
        stats[i]['CCDNUM'] = hdr['CCDNUM']
        stats[i]['AIRMASS'] = hdr['AIRMASS']

        data=fitsio.read(inst[i]['PATH'],ext=2)

        mag_psf_corr = data['MAG_PSF'] - 25.0 + inst['MAG_ZERO'][i]

        stars,=np.where((data['SPREAD_MODEL'] < 0.002) & (mag_psf_corr < 22.0) & (data['FLAGS'] == 0))
        data=data[stars]

        if (new_hdr):
            # need to translate wcs
            wcs = esutil.wcsutil.WCS(hdr)
            ra, dec = wcs.image2sky(data['XWIN_IMAGE'],data['YWIN_IMAGE'])
        else :
            ra = data['ALPHAWIN_J2000']
            dec = data['DELTAWIN_J2000']

        if (i == 0) :
            # this is the first one...
            mcat = np.zeros((data.size,),dtype=elt)

            mcat['RA'][:,0] = ra
            mcat['DEC'][:,0] = dec
            mcat['EXPNUM'][:,0] = hdr['EXPNUM']
            mcat['BAND'][:,0] = band
            mcat['CCDNUM'][:,0] = hdr['CCDNUM']
            mcat['MAG_AUTO'][:,0] = data['MAG_AUTO'] - 25.0 + inst['MAG_ZERO'][i]
            mcat['MAGERR_AUTO'][:,0] = data['MAGERR_AUTO']
            mcat['MAG_PSF'][:,0] = data['MAG_PSF'] - 25.0 + inst['MAG_ZERO'][i]
            mcat['MAGERR_PSF'][:,0] = data['MAGERR_PSF']            

            mcat['RA_MEAN'][:] = ra
            mcat['DEC_MEAN'][:] = dec

        else:            
            # need to match these...
            h=esutil.htm.HTM(13)

            m1,m2,d12=h.match(ra,dec,mcat['RA_MEAN'],mcat['DEC_MEAN'],1./3600.0,maxmatch=1)
            if (m1.size > 0):
                # we have matches
                mcat['RA'][m2,i] = ra[m1]
                mcat['DEC'][m2,i] = dec[m1]
                mcat['EXPNUM'][m2,i] = hdr['EXPNUM']
                mcat['BAND'][m2,i] = band
                mcat['CCDNUM'][m2,i] = hdr['CCDNUM']
                mcat['MAG_AUTO'][m2,i] = data['MAG_AUTO'][m1] - 25.0 + inst['MAG_ZERO'][i]
                mcat['MAGERR_AUTO'][m2,i] = data['MAGERR_AUTO'][m1]
                mcat['MAG_PSF'][m2,i] = data['MAG_PSF'][m1] - 25.0 + inst['MAG_ZERO'][i]
                mcat['MAGERR_PSF'][m2,i] = data['MAGERR_PSF'][m1]
                
                
                for kk in m2:                
                    # need to recompute mean positions...
                    gd,=np.where(mcat['EXPNUM'][kk,:] > 0)
                    mcat['RA_MEAN'][kk] = np.mean(mcat['RA'][kk,gd])
                    mcat['RA_RMS'][kk] = np.std(mcat['RA'][kk,gd])
                    mcat['DEC_MEAN'][kk] = np.mean(mcat['DEC'][kk,gd])
                    mcat['DEC_RMS'][kk] = np.std(mcat['DEC'][kk,gd])
                    
            if (m1.size < data.size) :
                # we have misses
                indices=np.arange(data.size)
                miss = np.delete(indices,m1)

                tempcat=np.zeros(miss.size,dtype=elt)
                tempcat['RA'][:,i] = ra[miss]
                tempcat['DEC'][:,i] = dec[miss]
                tempcat['EXPNUM'][:,i] = hdr['EXPNUM']
                tempcat['BAND'][:,i] = band
                tempcat['CCDNUM'][:,i] = hdr['CCDNUM']
                tempcat['MAG_AUTO'][:,i] = data['MAG_AUTO'][miss] - 25.0 + inst['MAG_ZERO'][i]
                tempcat['MAGERR_AUTO'][:,i] = data['MAGERR_AUTO'][miss]
                tempcat['MAG_PSF'][:,i] = data['MAG_PSF'][miss] - 25.0 + inst['MAG_ZERO'][i]
                tempcat['MAGERR_PSF'][:,i] = data['MAGERR_PSF'][miss]
                
                tempcat['RA_MEAN'][:] = ra[miss]
                tempcat['DEC_MEAN'][:] = dec[miss]

                # and append
                mcat=np.append(mcat,tempcat)

    # compute mean positions and magnitudes...this is object by object
    for i in range(nmag):
        for j in range(mcat.size):
            gd,=np.where((mcat['EXPNUM'][j,:] > 0) & (mcat['MAG_AUTO'][j,:] < 25.0) & (mcat['MAG_PSF'][j,:] < 25.0) & (mcat['BAND'][j,:] == bands[i]))

            mcat['NGOOD'][j,i] = gd.size

            if (gd.size >= 3) :
                wt=1./mcat['MAGERR_AUTO'][j,gd]**2
                mcat['MAG_AUTO_MEAN'][j,i] = np.sum(mcat['MAG_AUTO'][j,gd]*wt)/np.sum(wt)
                V1 = np.sum(wt)
                V2 = np.sum(wt**2)
                mcat['MAG_AUTO_RMS'][j,i] = np.sqrt((V1/(V1**2.-V2))*np.sum(wt*(mcat['MAG_AUTO'][j,gd] - mcat['MAG_AUTO_MEAN'][j,i])**2.))

                wt=1./mcat['MAGERR_PSF'][j,gd]**2
                mcat['MAG_PSF_MEAN'][j,i] = np.sum(mcat['MAG_PSF'][j,gd]*wt)/np.sum(wt)
                V1 = np.sum(wt)
                V2 = np.sum(wt**2)
                mcat['MAG_PSF_RMS'][j,i] = np.sqrt((V1/(V1**2.-V2))*np.sum(wt*(mcat['MAG_PSF'][j,gd] - mcat['MAG_PSF_MEAN'][j,i])**2.))

    # and finally look at the offsets relative to the means, image by image
    for i in range(inst.size):
        # get all the magnitudes from the ith observation...
        mag_auto = mcat['MAG_AUTO'][:,i]
        magerr_auto = mcat['MAGERR_AUTO'][:,i]
        mag_psf = mcat['MAG_PSF'][:,i]
        magerr_psf = mcat['MAGERR_PSF'][:,i]

        ok,=np.where(mcat['BAND'][:,i] != '')
        if (ok.size == 0):
            continue

        test,=np.where(mcat['BAND'][ok[0],i] == np.array(bands))
        test=test[0]
        
        mag_auto_mean = mcat['MAG_AUTO_MEAN'][:,test]
        mag_psf_mean = mcat['MAG_PSF_MEAN'][:,test]

        ra = mcat['RA'][:,i]
        dec = mcat['DEC'][:,i]

        ra_mean = mcat['RA_MEAN']
        dec_mean = mcat['DEC_MEAN']

        expnum = mcat['EXPNUM'][:,i]

        imag_psf = mcat['MAG_PSF_MEAN'][:,2]

        gd,=np.where((expnum > 0) & (mag_auto < 25.0) & (mag_psf < 25.0) & (np.abs(mag_auto - mag_auto_mean) < 0.2) & (np.abs(mag_psf - mag_psf_mean) < 0.2) & (imag_psf > 0.0) & (imag_psf < 21.0))

        stats['NGOOD'][i] = gd.size

        if (gd.size < 10):
            continue

        wt=1./(magerr_auto[gd]**2. + 0.003**2.)
        dm=mag_auto[gd] - mag_auto_mean[gd]
        stats['MAG_AUTO_OFFSET'][i] = np.sum(dm*wt)/np.sum(wt)
        stats['MAG_AUTO_OFFSET_SIGMA'][i] = np.sqrt(1./np.sum(wt))

        wt=1./(magerr_psf[gd]**2. + 0.003**2.)
        dm=mag_psf[gd] - mag_psf_mean[gd]
        stats['MAG_PSF_OFFSET'][i] = np.sum(dm*wt)/np.sum(wt)
        stats['MAG_PSF_OFFSET_SIGMA'][i] = np.sqrt(1./np.sum(wt))

        dr = ra[gd] - ra_mean[gd]
        stats['RA_OFFSET'][i] = np.mean(dr)
        stats['RA_OFFSET_RMS'][i] = np.std(dr)

        dd = dec[gd] - dec_mean[gd]
        stats['DEC_OFFSET'][i] = np.mean(dd)
        stats['DEC_OFFSET_RMS'][i] = np.std(dd)

    fitsio.write(outfile, mcat, clobber=True)
    fitsio.write(outfile, stats)
