import healpy as hp
import fitsio
import numpy as np

def read_map2(filename,field=0,dtype=np.float64,nest=False,hdu=1,h=False,verbose=True,memmap=False):
    hdr=fitsio.read_header(filename,ext=hdu)
    
    fullsky = False
    try:
        if (hdr['OBJECT'].strip() == 'PARTIAL') :
            # partial sky format
            fullsky=False
        else:
            fullsky=True
    except:
        # if no OBJECT in header, assume full sky
        fullsky=True

    if fullsky:
        m=hp.read_map(filename,field=field,dtype=dtype,nest=nest,hdu=hdu,h=h,verbose=verbose,memmap=memmap)
    else:
        # partial sky
        st=fitsio.read(filename,ext=1)
        nside=hdr['NSIDE']

        m=np.zeros(12*nside*nside,dtype=dtype) + hp.UNSEEN

        if ((hdr['ORDERING'].strip() == 'NESTED') and (not nest)) :
            # change from nest to ring...
            m[hp.nest2ring(nside,st['PIXEL'])] = st['SIGNAL']
        elif ((hdr['ORDERING'].strip() == 'RING') and (nest)):
            # change from ring to nest...
            m[hp.ring2nest(nside,st['PIXEL'])] = st['SIGNAL']
        else :
            # straight up
            m[st['PIXEL']] = st['SIGNAL']

    return m
