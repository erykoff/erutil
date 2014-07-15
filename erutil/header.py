"""
Read astromatic header files

read_astromatic_header

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


import pyfits
import numpy as np



# determine pyfits version
vers = pyfits.__version__.split('.')
if (vers[0] >= 3) and (int(vers[1]) >= 1)  :
    new_pyfits = True
else:
    new_pyfits = False

def read_astromatic_header(ofile,headfile,ext=0):
    """
    Read in an astromatic header (.head) text file and combine with original header

    parameters
    ----------
    ofile: original file.  May be FITS_LDAC or image.
    headfile: astromatic header text file
    ext: extension of ofile; only used if not FITS_LDAC

    returns
    -------
    pyfits header

    Written by Eli Rykoff, SLAC, 2014
    
    """

    # determine if ofile is a FITS_LDAC file
    is_ldac = False
    
    hdulist=pyfits.open(ofile)
    if (len(hdulist) > 0):
        # definitely not 
        hdunames=[]
        for hdu in hdulist[1:len(hdulist)]:
            hdunames.append(hdu.name)
        if ('LDAC_IMHEAD' in hdunames):
            is_ldac = True

    if (is_ldac):
        odata=hdulist['LDAC_IMHEAD'].data[0][0]

        hdr=pyfits.Header()
        for j in range(odata.size):
            card=pyfits.Card.fromstring(odata[j])
            try:
                test = card.value
            except:
                card.verify('fix')

            if (card.key == 'END'):
                continue

            if (new_pyfits):
                hdr.append(card)
            else:
                if (len(card.key) > 8) :
                    hdr.update('hierarch '+card.key, card.value, comment=card.comment)
                else:
                    hdr.update(card.key, card.value, comment=card.comment)

    else:
        hdr=hdulist[ext].header.copy()
        
    hdulist.close()

    fs=open(headfile)
    for text in fs:
        card=pyfits.Card.fromstring(text.strip())
        if (card.key == 'END' or card.key == 'COMMENT'):
            continue
        if (new_pyfits):
            hdr[card.key] = (card.value, card.comment)
        else :
            if (len(card.key) > 8):
                hdr.update('hierarch '+card.key, card.value, comment=card.value)
            else:
                hdr.update(card.key, card.value, comment=card.comment)
    fs.close()

    return hdr
    
        
