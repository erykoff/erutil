#!/usr/bin/env python

import healpy as hp
import numpy as np
import fitsio
import esutil
from .read_map2 import *


def deres_map(nside_out,minfrac=0.8,minsub=1,nest=False,range=[hp.UNSEEN,np.abs(hp.UNSEEN)],fracfile=None,frac_in=None,return_errmap=False,errfile=None,infile=None,outfile=None,map_in=None):

    if (map_in is not None and infile is not None):
        raise RuntimeError("Cannot set both map_in and infile")
    if (map_in is None and infile is None):
        raise RuntimeError("Must set one of map_in or infile")
    if (frac_in is not None and fracfile is not None):
        raise RuntimeError("Cannot set both frac_in and fracfile")

    if (infile is not None):
        _map_in = read_map2(infile,hdu=1,dtype=np.float32)
    else:
        _map_in = map_in.copy()
    nside_in = hp.npix2nside(_map_in.size)

    if (fracfile is not None):
        frac_in = read_map2(fracfile,hdu=1,dtype=np.float32)

    if (frac_in is not None):
        bad,=np.where((frac_in < minfrac) & (frac_in > 0.0))
        if (bad.size > 0):
            _map_in[bad] = hp.UNSEEN
            
    # generate pixels at old nside
    theta_in,phi_in = hp.pix2ang(nside_in, np.arange(12*nside_in*nside_in),nest=nest)

    use_in, = np.where((_map_in > range[0]) & (_map_in < range[1]))

    # and put in out nside; note that we may want nest output
    hpix_out = hp.ang2pix(nside_out,theta_in[use_in],phi_in[use_in],nest=nest)

    # what are the out pixels?
    #st=np.argsort(hpix_out)
    hpix_out_u=np.unique(hpix_out)
    npix = hpix_out_u.size

    # match together
    suba,subb=esutil.numpy_util.match_multi(hpix_out_u,hpix_out)

    # make a compressed index
    fake_id=np.arange(npix)+1
    # and split with histogram
    h,rev=esutil.stat.histogram(fake_id[suba],min=0,rev=True)

    # finally, prep the out map
    map_out = np.zeros(12*nside_out*nside_out,dtype=np.float32) + hp.UNSEEN

    do_err=False
    if (errfile is not None or return_errmap):
        map_err = np.zeros(12*nside_out*nside_out,dtype=np.float32)
        do_err=True        
        
        
    for i in xrange(npix):
        id=fake_id[i]
        if (rev[id] < rev[id+1]):
            i1a=rev[rev[id]:rev[id+1]]

            pind=subb[i1a]

            if (pind.size >= minsub):
                map_out[hpix_out_u[i]] = np.mean(_map_in[use_in[pind]])

                if do_err:
                    map_err[hpix_out_u[i]] = np.std(_map_in[use_in[pind]])


    if (outfile is not None):
        hp.write_map(outfile,map_out,nest=nest,coord='C')

    if (errfile is not None):
        hp.write_map(errfile,map_err,nest=nest,coord='C')
        

    if (return_errmap):
        return map_out,map_err
    else :
        return map_out
            
    
                  
                  
    
    
