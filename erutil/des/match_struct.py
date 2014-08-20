"""

need FITS_LDAC ...

"""

import fitsio
import esutil
import numpy as np
import sys

from erutil import header


def make_match_struct(infile, outfile, bands=['g','r','i','z'],matchrad=1.0,statsonly=False,maximag=21.0,quiet=False,meanonly=False):

    inst=fitsio.read(infile,ext=1)

    if (inst.size == 0):
        print "Error! Nothing in info file " + infile
        return


    new_hdr = False
    if ('HDRPATH' in inst.dtype.names) :
        new_hdr = True

    nmag=len(bands)
    nimage=inst.size

    # first, a master list...

    elt=[('RA','f8'),\
             ('DEC','f8'),\
             ('EXPNUM','i4'),\
             ('CCDNUM','i2'),\
             ('BAND','a1'),\
             ('MAG_AUTO','f4'),\
             ('MAGERR_AUTO','f4'),\
             ('MAG_PSF','f4'),\
             ('MAGERR_PSF','f4')]

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

    stats=np.recarray((nimage,),dtype=selt)
    mark=0

    for i in range(inst.size):
        if (not quiet):
            sys.stdout.write('\b'*20)
            sys.stdout.write('%d of %d' % (i, inst.size))
            sys.stdout.flush()

        if (new_hdr):
            hdr = header.read_astromatic_header2(inst[i]['PATH'].strip(), inst[i]['HDRPATH'].strip())
        else :
            hdr = header.read_astromatic_header2(inst[i]['PATH'].strip(), None)

        filt=hdr['FILTER']
        band=filt[0]

        stats[i]['BAND'] = band
        stats[i]['EXPNUM'] = hdr['EXPNUM']
        stats[i]['CCDNUM'] = hdr['CCDNUM']
        stats[i]['AIRMASS'] = hdr['AIRMASS']

        data=fitsio.read(inst[i]['PATH'].strip(),ext=2)

        mag_psf_corr = data['MAG_PSF'] - 25.0 + inst['MAG_ZERO'][i]

        stars,=np.where((data['SPREAD_MODEL'] < 0.002) & (mag_psf_corr < 22.0) & (data['FLAGS'] == 0) & (data['MAGERR_PSF'] < 0.1))
        data=data[stars]

        if (new_hdr):
            # need to translate wcs
            wcs = esutil.wcsutil.WCS(hdr)
            ra, dec = wcs.image2sky(data['XWIN_IMAGE'],data['YWIN_IMAGE'])
        else :
            ra = data['ALPHAWIN_J2000']
            dec = data['DELTAWIN_J2000']

        if (i == 0):
            nstart = data.size*4
            cat=np.zeros((nstart,),dtype=elt)
        elif (cat.size < (mark+data.size)):
            # we need to expand
            tempcat=np.zeros(cat.size,dtype=elt)
            cat=np.append(cat,tempcat)
            tempcat=0
            
        cat['RA'][mark:mark+data.size] = ra
        cat['DEC'][mark:mark+data.size] = dec
        cat['EXPNUM'][mark:mark+data.size] = hdr['EXPNUM']
        cat['BAND'][mark:mark+data.size] = band
        cat['CCDNUM'][mark:mark+data.size] = hdr['CCDNUM']
        cat['MAG_AUTO'][mark:mark+data.size] = data['MAG_AUTO'] - 25.0 + inst['MAG_ZERO'][i]
        cat['MAGERR_AUTO'][mark:mark+data.size] = data['MAGERR_AUTO']
        cat['MAG_PSF'][mark:mark+data.size] = data['MAG_PSF'] - 25.0 + inst['MAG_ZERO'][i]
        cat['MAGERR_PSF'][mark:mark+data.size] = data['MAGERR_PSF']

        mark=mark+data.size

    # crop it down...
    cat=cat[0:mark]

    if not quiet:
        print ""
        print "Matching objects..."
        
    # now we need to match...
    htm=esutil.htm.HTM(13)
    ind1,ind2,d12=htm.match(cat['RA'],cat['DEC'],cat['RA'],cat['DEC'],matchrad/3600.,maxmatch=0)

    # and break out the histogram...

    fake_id=np.arange(cat.size)
    hist,rev=esutil.stat.histogram(fake_id[ind1],rev=True)

    if (hist.max() == 1) :
        print "Error: only a single unique list of objects!"
        fitsio.write(outfile,stats, clobber=True)
        return

    maxobs=hist.max()

    # build the match structure...
    elt = [('RA','>f8',maxobs),\
               ('DEC','>f8',maxobs),\
               ('EXPNUM','i4',maxobs),\
               ('CCDNUM','i2',maxobs),\
               ('BAND','a1',maxobs),\
               ('MAG_AUTO','f4',maxobs),\
               ('MAGERR_AUTO','f4',maxobs),\
               ('MAG_PSF','f4',maxobs),\
               ('MAGERR_PSF','f4',maxobs),\
               ('NGOOD','i4',nmag),\
               ('MAG_AUTO_MEAN','f4',nmag),\
               ('MAG_AUTO_RMS','f4',nmag),\
               ('MAG_PSF_MEAN','f4',nmag),\
               ('MAG_PSF_RMS','f4',nmag),\
               ('RA_MEAN','f8'),\
               ('RA_RMS','f8'),\
               ('DEC_MEAN','f8'),\
               ('DEC_RMS','f8')]

    # this needs to be the number of unique objects!
    
    hist_temp = hist.copy()
    count=0
    for j in range(hist_temp.size):
        jj=fake_id[j]
        if (hist_temp[jj] > 0):
            i1a=rev[rev[jj]:rev[jj+1]]
            hist_temp[ind2[i1a]] = 0
            count=count+1

    
    mcat=np.zeros(count,dtype=elt)

    if not quiet:
        print "Measuring object stats..."
        
    index=0
    for j in range(hist.size):
        jj=fake_id[j]
        if (hist[jj] > 0):
            i1a=rev[rev[jj]:rev[jj+1]]
            
            gind=ind2[i1a]

            # and zero out the other entries in the histogram
            hist[gind] = 0
            
            mcat['RA'][index,0:gind.size] = cat['RA'][gind]
            mcat['DEC'][index,0:gind.size] = cat['DEC'][gind]
            mcat['EXPNUM'][index,0:gind.size] = cat['EXPNUM'][gind]
            mcat['CCDNUM'][index,0:gind.size] = cat['CCDNUM'][gind]
            mcat['BAND'][index,0:gind.size] = cat['BAND'][gind]
            mcat['MAG_AUTO'][index,0:gind.size] = cat['MAG_AUTO'][gind]
            mcat['MAGERR_AUTO'][index,0:gind.size] = cat['MAGERR_AUTO'][gind]
            mcat['MAG_PSF'][index,0:gind.size] = cat['MAG_PSF'][gind]
            mcat['MAGERR_PSF'][index,0:gind.size] = cat['MAGERR_PSF'][gind]

            # all bands together...
            mcat['RA_MEAN'][index] = np.mean(mcat['RA'][index,0:gind.size])
            mcat['RA_RMS'][index] = np.std(mcat['RA'][index,0:gind.size])
            mcat['DEC_MEAN'][index] = np.mean(mcat['DEC'][index,0:gind.size])
            mcat['DEC_RMS'][index] = np.std(mcat['DEC'][index,0:gind.size])
            
            for i in range(nmag):
                gd,=np.where((mcat['MAG_AUTO'][index,0:gind.size] < 25) &
                             (mcat['MAG_PSF'][index,0:gind.size] < 25) &
                             (mcat['BAND'][index,0:gind.size] == bands[i]))

                mcat['NGOOD'][index,i] = gd.size

                if (gd.size >= 3):
                    wt=1./mcat['MAGERR_AUTO'][index,gd]**2
                    mcat['MAG_AUTO_MEAN'][index,i] = np.sum(mcat['MAG_AUTO'][index,gd]*wt)/np.sum(wt)
                    V1 = np.sum(wt)
                    V2 = np.sum(wt**2)
                    mcat['MAG_AUTO_RMS'][index,i] = np.sqrt((V1/(V1**2.-V2))*np.sum(wt*(mcat['MAG_AUTO'][index,gd] - mcat['MAG_AUTO_MEAN'][index,i])**2.))
                                                    
                    wt=1./mcat['MAGERR_PSF'][index,gd]**2
                    mcat['MAG_PSF_MEAN'][index,i] = np.sum(mcat['MAG_PSF'][index,gd]*wt)/np.sum(wt)
                    V1 = np.sum(wt)
                    V2 = np.sum(wt**2)
                    mcat['MAG_PSF_RMS'][index,i] = np.sqrt((V1/(V1**2.-V2))*np.sum(wt*(mcat['MAG_PSF'][index,gd] - mcat['MAG_PSF_MEAN'][index,i])**2.))

            index=index+1

    # need to go image by image

    if not quiet:
        print "Measuring image stats..."

    # make a hash and histogram it...
    exphash = (mcat['EXPNUM']*100 + mcat['CCDNUM']).ravel()
    stathash = stats['EXPNUM']*100 + stats['CCDNUM']

    # hopefully doesn't use too much memory...
    smin=stathash.min()
    smax=stathash.max()
    
    hist,rev=esutil.stat.histogram(exphash,rev=True,min=smin,max=smax)

    for i in range(stats.size):
        if (hist[stathash[i]-smin] > 0):
            i1a=rev[rev[stathash[i]-smin]:rev[stathash[i]-smin+1]]

            xx=i1a % maxobs
            yy=i1a / maxobs

            # get the magnitudes of all the stars in the image
            mag_auto=mcat['MAG_AUTO'][yy,xx]
            magerr_auto=mcat['MAGERR_AUTO'][yy,xx]
            mag_psf=mcat['MAG_PSF'][yy,xx]
            magerr_psf=mcat['MAGERR_PSF'][yy,xx]

            # we need to know the band...
            test,=np.where(stats['BAND'][i] == np.array(bands))
            test=test[0]

            # and the mean magnitudes for each of these objects
            mag_auto_mean = mcat['MAG_AUTO_MEAN'][yy,test]
            mag_psf_mean = mcat['MAG_PSF_MEAN'][yy,test]

            # individual ra/dec
            ra=mcat['RA'][yy,xx]
            dec=mcat['DEC'][yy,xx]

            # mean ra/dec
            ra_mean = mcat['RA_MEAN'][yy]
            dec_mean = mcat['DEC_MEAN'][yy]

            imag_psf = mcat['MAG_PSF_MEAN'][yy,2]

            gd,=np.where((mag_auto < 25.0) & (mag_psf < 25.0) & (np.abs(mag_auto-mag_auto_mean) < 0.2) & (np.abs(mag_psf - mag_psf_mean) < 0.2) & (imag_psf > 0.0) & (imag_psf < maximag))

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

    if not quiet:
        print "Saving..."
        
    if (not statsonly):
        fitsio.write(outfile, mcat, clobber=True)
        fitsio.write(outfile, stats)
    else :
        fitsio.write(outfile,stats, clobber=True)
                             


def make_match_struct_gals(infile, outfile, bands=['g','r','i','z'],matchrad=1.0,statsonly=False,maximag=22.0,quiet=False):

    inst=fitsio.read(infile,ext=1)

    if (inst.size == 0):
        print "Error! Nothing in info file " + infile
        return


    new_hdr = False
    if ('HDRPATH' in inst.dtype.names) :
        new_hdr = True

    nmag=len(bands)
    nimage=inst.size

    # first, a master list...

    elt=[('RA','f8'),\
             ('DEC','f8'),\
             ('EXPNUM','i4'),\
             ('CCDNUM','i2'),\
             ('BAND','a1'),\
             ('MAG_AUTO','f4'),\
             ('MAGERR_AUTO','f4')]

    # this is image-by-image
    selt=[('EXPNUM','i4'),\
              ('CCDNUM','i4'),\
              ('BAND','a2'),\
              ('NGOOD','i4'),\
              ('MAG_AUTO_OFFSET','f4'),\
              ('MAG_AUTO_OFFSET_SIGMA','f4'),\
              ('RA_OFFSET','f8'),\
              ('RA_OFFSET_RMS','f8'),\
              ('DEC_OFFSET','f8'),\
              ('DEC_OFFSET_RMS','f8')]

    stats=np.recarray((nimage,),dtype=selt)
    mark=0

    for i in range(inst.size):
        if (not quiet):
            sys.stdout.write('\b'*20)
            sys.stdout.write('%d of %d' % (i, inst.size))
            sys.stdout.flush()

        if (new_hdr):
            hdr = header.read_astromatic_header2(inst[i]['PATH'].strip(), inst[i]['HDRPATH'].strip())
        else :
            hdr = header.read_astromatic_header2(inst[i]['PATH'].strip(), None)

        filt=hdr['FILTER']
        band=filt[0]

        stats[i]['BAND'] = band
        stats[i]['EXPNUM'] = hdr['EXPNUM']
        stats[i]['CCDNUM'] = hdr['CCDNUM']

        data=fitsio.read(inst[i]['PATH'].strip(),ext=2)

        mag_auto_corr = data['MAG_AUTO'] - 25.0 + inst['MAG_ZERO'][i]

        gals,=np.where((data['SPREAD_MODEL'] > 0.002) & (mag_auto_corr < 24.0) & (data['FLAGS'] <= 3))
        data=data[gals]

        if (new_hdr):
            # need to translate wcs
            wcs = esutil.wcsutil.WCS(hdr)
            ra, dec = wcs.image2sky(data['XWIN_IMAGE'],data['YWIN_IMAGE'])
        else :
            ra = data['ALPHAWIN_J2000']
            dec = data['DELTAWIN_J2000']

        if (i == 0):
            nstart = data.size*4
            cat=np.zeros((nstart,),dtype=elt)
        elif (cat.size < (mark+data.size)):
            # we need to expand
            tempcat=np.zeros(cat.size,dtype=elt)
            cat=np.append(cat,tempcat)
            tempcat=0
            
        cat['RA'][mark:mark+data.size] = ra
        cat['DEC'][mark:mark+data.size] = dec
        cat['EXPNUM'][mark:mark+data.size] = hdr['EXPNUM']
        cat['BAND'][mark:mark+data.size] = band
        cat['CCDNUM'][mark:mark+data.size] = hdr['CCDNUM']
        cat['MAG_AUTO'][mark:mark+data.size] = data['MAG_AUTO'] - 25.0 + inst['MAG_ZERO'][i]
        cat['MAGERR_AUTO'][mark:mark+data.size] = data['MAGERR_AUTO']

        mark=mark+data.size

    # crop it down...
    cat=cat[0:mark]

    if not quiet:
        print ""
        print "Matching objects..."
        
    # now we need to match...
    htm=esutil.htm.HTM(13)
    ind1,ind2,d12=htm.match(cat['RA'],cat['DEC'],cat['RA'],cat['DEC'],matchrad/3600.,maxmatch=0)

    # and break out the histogram...

    fake_id=np.arange(cat.size)
    hist,rev=esutil.stat.histogram(fake_id[ind1],rev=True)

    if (hist.max() == 1) :
        print "Error: only a single unique list of objects!"
        fitsio.write(outfile,stats, clobber=True)
        return

    maxobs=hist.max()

    # build the match structure...
    elt = [('RA','>f8',maxobs),\
               ('DEC','>f8',maxobs),\
               ('EXPNUM','i4',maxobs),\
               ('CCDNUM','i2',maxobs),\
               ('BAND','a1',maxobs),\
               ('MAG_AUTO','f4',maxobs),\
               ('MAGERR_AUTO','f4',maxobs),\
               ('NGOOD','i4',nmag),\
               ('MAG_AUTO_MEAN','f4',nmag),\
               ('MAG_AUTO_RMS','f4',nmag),\
               ('RA_MEAN','f8'),\
               ('RA_RMS','f8'),\
               ('DEC_MEAN','f8'),\
               ('DEC_RMS','f8')]

    # this needs to be the number of unique objects!
    
    hist_temp = hist.copy()
    count=0
    for j in range(hist_temp.size):
        jj=fake_id[j]
        if (hist_temp[jj] > 0):
            i1a=rev[rev[jj]:rev[jj+1]]
            hist_temp[ind2[i1a]] = 0
            count=count+1

    
    mcat=np.zeros(count,dtype=elt)

    if not quiet:
        print "Measuring object stats..."
        
    index=0
    for j in range(hist.size):
        jj=fake_id[j]
        if (hist[jj] > 0):
            i1a=rev[rev[jj]:rev[jj+1]]
            
            gind=ind2[i1a]

            # and zero out the other entries in the histogram
            hist[gind] = 0
            
            mcat['RA'][index,0:gind.size] = cat['RA'][gind]
            mcat['DEC'][index,0:gind.size] = cat['DEC'][gind]
            mcat['EXPNUM'][index,0:gind.size] = cat['EXPNUM'][gind]
            mcat['CCDNUM'][index,0:gind.size] = cat['CCDNUM'][gind]
            mcat['BAND'][index,0:gind.size] = cat['BAND'][gind]
            mcat['MAG_AUTO'][index,0:gind.size] = cat['MAG_AUTO'][gind]
            mcat['MAGERR_AUTO'][index,0:gind.size] = cat['MAGERR_AUTO'][gind]

            # all bands together...
            mcat['RA_MEAN'][index] = np.mean(mcat['RA'][index,0:gind.size])
            mcat['RA_RMS'][index] = np.std(mcat['RA'][index,0:gind.size])
            mcat['DEC_MEAN'][index] = np.mean(mcat['DEC'][index,0:gind.size])
            mcat['DEC_RMS'][index] = np.std(mcat['DEC'][index,0:gind.size])
            
            for i in range(nmag):
                gd,=np.where((mcat['MAG_AUTO'][index,0:gind.size] < 25) &
                             (mcat['BAND'][index,0:gind.size] == bands[i]))

                mcat['NGOOD'][index,i] = gd.size

                if (gd.size >= 3):
                    wt=1./mcat['MAGERR_AUTO'][index,gd]**2
                    mcat['MAG_AUTO_MEAN'][index,i] = np.sum(mcat['MAG_AUTO'][index,gd]*wt)/np.sum(wt)
                    V1 = np.sum(wt)
                    V2 = np.sum(wt**2)
                    mcat['MAG_AUTO_RMS'][index,i] = np.sqrt((V1/(V1**2.-V2))*np.sum(wt*(mcat['MAG_AUTO'][index,gd] - mcat['MAG_AUTO_MEAN'][index,i])**2.))
                                                    

            index=index+1

    # need to go image by image

    if not quiet:
        print "Measuring image stats..."

    # make a hash and histogram it...
    exphash = (mcat['EXPNUM']*100 + mcat['CCDNUM']).ravel()
    stathash = stats['EXPNUM']*100 + stats['CCDNUM']

    # hopefully doesn't use too much memory...
    smin=stathash.min()
    smax=stathash.max()
    
    hist,rev=esutil.stat.histogram(exphash,rev=True,min=smin,max=smax)

    for i in range(stats.size):
        if (hist[stathash[i]-smin] > 0):
            i1a=rev[rev[stathash[i]-smin]:rev[stathash[i]-smin+1]]

            xx=i1a % maxobs
            yy=i1a / maxobs

            # get the magnitudes of all the gals in the image
            mag_auto=mcat['MAG_AUTO'][yy,xx]
            magerr_auto=mcat['MAGERR_AUTO'][yy,xx]

            # we need to know the band...
            test,=np.where(stats['BAND'][i] == np.array(bands))
            test=test[0]

            # and the mean magnitudes for each of these objects
            mag_auto_mean = mcat['MAG_AUTO_MEAN'][yy,test]

            # individual ra/dec
            ra=mcat['RA'][yy,xx]
            dec=mcat['DEC'][yy,xx]

            # mean ra/dec
            ra_mean = mcat['RA_MEAN'][yy]
            dec_mean = mcat['DEC_MEAN'][yy]

            imag_auto = mcat['MAG_AUTO_MEAN'][yy,2]

            gd,=np.where((mag_auto < 25.0) & (np.abs(mag_auto-mag_auto_mean) < 0.2) & (imag_auto > 0.0) & (imag_auto < maximag))

            stats['NGOOD'][i] = gd.size
            
            if (gd.size < 10):
                continue

            wt=1./(magerr_auto[gd]**2. + 0.003**2.)
            dm=mag_auto[gd] - mag_auto_mean[gd]
            stats['MAG_AUTO_OFFSET'][i] = np.sum(dm*wt)/np.sum(wt)
            stats['MAG_AUTO_OFFSET_SIGMA'][i] = np.sqrt(1./np.sum(wt))
            
            
            dr = ra[gd] - ra_mean[gd]
            stats['RA_OFFSET'][i] = np.mean(dr)
            stats['RA_OFFSET_RMS'][i] = np.std(dr)
            
            dd = dec[gd] - dec_mean[gd]
            stats['DEC_OFFSET'][i] = np.mean(dd)
            stats['DEC_OFFSET_RMS'][i] = np.std(dd)

    if not quiet:
        print "Saving..."
        
    if (not statsonly):
        fitsio.write(outfile, mcat, clobber=True)
        fitsio.write(outfile, stats)
    else :
        fitsio.write(outfile,stats, clobber=True)
                             



