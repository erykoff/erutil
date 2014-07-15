"""
Matching functions

close_match

close_match_radec

"""

license="""
  Copyright (C) 2014 Eli Rykoff

    This program is free software; you can redistribute it and/or modify it
    under the terms of version 2 of the GNU General Public License as
    published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
"""

import numpy as np
import esutil.coords as coords



def close_match_radec(ra1,dec1,ra2,dec2,ep,allow,silent=False,box=False) :
    """
    Find the nearest neighbors between two arrays of ra/dec

    parameters
    ----------
    ra1,dec1: scalar or array
         coordinates of a set of points (degrees).  Must be same length.
    ra2,dec2: scalar or array
         coordinates of a second set of points (degrees).  Must be same length.
    ep: scalar
         maximum match distance between pairs (degrees)
    allow: scalar
         maximum number of matches in second array to each element in first array.
    silent: boolean
         make quiet
    box: boolean
         use square of size ep (default circle):

    Original by Dave Johnston, University of Michigan, 1997
    
    Translated from IDL by Eli Rykoff, SLAC
    """
    ra1=np.atleast_1d(ra1)
    dec1=np.atleast_1d(dec1)
    ra2=np.atleast_1d(ra2)
    dec2=np.atleast_1d(dec2)
    
    epdec = ep

    n1=ra1.size
    n2=ra2.size

    matcharr=np.zeros([n1,allow],dtype=np.int32)
    matcharr.fill(-1)
    ind=np.arange(n2,dtype=np.int32)
    sor=ra2.argsort()
    ra2sort=ra2[sor]
    dec2sort=dec2[sor]
    ind=ind[sor]
    runi=0
    endra2=ra2sort[n2-1]

    for i in range(n1) :
        epra=ep/np.cos(dec1[i]*0.01745329)
        ra1minus = ra1[i]-epra
        ra1plus = ra1[i]+epra
        in1=_binary_search(ra2sort,ra1minus)

        if (in1 == -1) :
            if (ra1minus < endra2): in1=0
        if (in1 != -1):
            in2=in1
            jj=in2+1
            while (jj < n2):
                if (ra2sort[in2+1] < ra1plus):
                    in2+=1
                    jj+=1
                else:
                    jj=n2
            if (n2 == 1):
                in2 = 0

            if (in1 <= in2) :
                dec2check=dec2sort[in1:in2+1]
                ra2check=ra2sort[in1:in2+1]

                decoffby=np.abs(dec2check-dec1[i])
                raoffby=np.abs(ra2check-ra1[i])
                good=np.where((decoffby < epdec) & \
                              (raoffby < epra))[0]+in1
                ngood=good.size
                if (ngood != 0):
                    if (not box):
                        offby = coords.sphdist(ra1[i],dec1[i],\
                                                      ra2sort[good],dec2sort[good],\
                                                      units=['deg','deg'])
                        good_offby=np.where(offby <= ep)[0]
                        ngood = good_offby.size
                    else :
                        good_offby = np.arange(ngood)
                        offby = raoffby

                    if (ngood != 0) :
                        good = good[good_offby]
                        offby=offby[good_offby]
                        if (ngood > allow) :
                            good=good[offby.argsort()]
                            ngood = allow
                            good=good[0:allow]

                        matcharr[i,0:ngood]=good
                        runi=runi+ngood
    
        
    if (not silent): print "total put in bytarr:",runi
    matches=np.where(matcharr != -1)
    if (matches[0].size == 0):
        if (not silent):print "no matches found"
        m1=np.array([-1])
        m2=np.array([-1])
        return (m1,m2)
    m1=matches[0] % n1
    m2=matcharr[matches]
    m2=ind[m2].flatten()
    if (not silent):print m1.size,' matches'
    return (m1,m2)

def close_match(t1,s1,t2,s2,ep,allow,silent=False,circle=False) :
    """
    Find the nearest neighbors between two arrays of x/y

    parameters
    ----------
    ra1,dec1: scalar or array
         coordinates of a set of points.  Must be same length.
    ra2,dec2: scalar or array
         coordinates of a second set of points.  Must be same length.
    ep: scalar
         maximum match distance between pairs (pixels)
    allow: scalar
         maximum number of matches in second array to each element in first array.
    silent: boolean
         make quiet
    circle: boolean
         use circle of size ep (default square):

    Original by Dave Johnston, University of Michigan, 1997
         
    Translated from IDL by Eli Rykoff, SLAC

    """
    t1=np.atleast_1d(t1)
    s1=np.atleast_1d(s1)
    t2=np.atleast_1d(t2)
    s2=np.atleast_1d(s2)

    n1=t1.size
    n2=t2.size

    matcharr=np.zeros([n1,allow],dtype=np.int32)
    matcharr.fill(-1)
    ind=np.arange(n2,dtype=np.int32)
    sor=t2.argsort()
    t2s=t2[sor]
    s2s=s2[sor]
    ind=ind[sor]
    runi=0
    endt=t2s[n2-1]

 
    for i in range(n1):
        t=t1[i]
        tm=t-ep
        tp=t+ep
        in1=_binary_search(t2s,tm)  # I can improve this?
        
        if in1 == -1:
            if (tm < endt) : in1=0
        if in1 != -1:
            in1=in1+1
            in2=in1-1
            jj=in2+1
            while (jj < n2):
                if (t2s[in2+1] < tp):
                    in2+=1
                    jj+=1
                else :
                    jj=n2
            if (n2 == 1) :
                in2=0  # hmmm

            if (in1 <= in2):
                if (n2 != 1) :
                    check = s2s[in1:in2+1]
                    tcheck = t2s[in1:in2+1]
                else :
                    check = s2s[0]
                    tcheck=t2s[0]
                s=s1[i]
                t=t1[i]
                offby=abs(check-s)
                toffby=abs(tcheck-t)
                good=np.where(np.logical_and(offby < ep,toffby < ep))[0]+in1
                ngood=good.size
                if (ngood != 0) :
                    if (ngood > allow) :
                        offby=offby[good-in1]
                        toffby=toffby[good-in1]
                        dist=np.sqrt(offby**2+toffby**2)
                        good=good[dist.argsort()]
                        ngood=allow
                    good=good[0:ngood]
                    matcharr[i,0:ngood]=good
                    run=runi+ngood
        

    if (not silent): print "total put in bytarr:",runi
    #matches=np.where(matcharr != -1)[0]
    matches=np.where(matcharr != -1)
    #if (matches.size == 0):
    if (matches[0].size == 0):
        if (not silent):print "no matches found"
        m1=np.array([-1])
        m2=np.array([-1])
        return (m1,m2)
    m1 = matches[0] % n1
    #m2 = matcharr[matches].flatten()
    m2 = matcharr[matches]
    m2 = ind[m2].flatten()
    if (not silent): print m1.size,' matches'
    return (m1,m2)



def _binary_search(arr,x,edgedefault=False,round=False):
    n=arr.size
    if (x < arr[0]) or (x > arr[n-1]):
        if (edgedefault):
            if (x < arr[0]): index = 0
            elif (x > arr[n-1]): index = n-1
        else: index = -1
        return index

    down=-1
    up=n
    while (up-down) > 1:
        mid=down+(up-down)/2
        if x >= arr[mid]:
            down=mid
        else:
            up=mid

    index=down

    if (round) and (index != n-1):
        if (abs(x-arr[index]) >= abs(x-arr[index+1])): index=index+1

    return index

