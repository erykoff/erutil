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
import fitsio
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
    headfile: astromatic header text file.  May be None
        If this is None, useful to read in FITS_LDAC into pyfits header.
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

    if (headfile is not None):    
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
    
        
def read_astromatic_header2(ofile,headfile,ext=0):
    """
    Read in an astromatic header (.head) text file and combine with original header

    parameters
    ----------
    ofile: original file.  May be FITS_LDAC or image.
    headfile: astromatic header text file.  May be None
        If this is None, useful to read in FITS_LDAC into pyfits header.
    ext: extension of ofile; only used if not FITS_LDAC

    returns
    -------
    fitsio header

    Written by Eli Rykoff, SLAC, 2014
    
    """

    # determine if ofile is a FITS_LDAC file
    is_ldac = False

    hdulist=fitsio.FITS(ofile)
    
    if 'LDAC_IMHEAD' in hdulist:
        is_ldac = True

    if (is_ldac):
        odata=hdulist['LDAC_IMHEAD']['Field Header Card'][0]

        hdr=fitsio.FITSHDR()
        for j in range(odata.size):
            line = odata[j].strip()
            if line != 'END' and line[0:8] != 'CONTINUE':
                hdr.add_record(line)

    else:
        hdr=hdulist[ext].read_header()

    hdulist.close()

    if (headfile is not None):
        fs=open(headfile)
        for text in fs:
            r=fitsio.FITSRecord(text.strip())
            if (r['name'] == 'COMMENT') or (r['name'] == 'HISTORY'):
                hdr.add_record(r)
            else :
                if (r['name'] not in hdr) :
                    hdr.add_record(r)
                else :
                    hdr[r['name']] = r['value']
            
        fs.close()

    return hdr
