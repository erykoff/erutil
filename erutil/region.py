import numpy as np
import os

class Regions:
    def __init__(self,x,y,index=None,color=None,radius=None,width=None,fk5=True,labelvals=None):
        self.x = x
        self.y = y

        if (len(self.x) != len(self.y)) :
            raise ValueError("x and y must be same length")

        self.labelvals=labelvals

        self.uselabels = False
        if (self.labelvals is not None):
            self.uselabels = True

        if (index is None) :
            index = range(len(x))
        self.index = index
        if (color is None) :
            color = "green"
        self.color = color
        if (radius is None) :
            if (not fk5) :
                radius = 10
            else :
                radius = 3
        self.radius=np.array(radius)

        if (self.radius.size == 1) :
            self.singlerad = True
        elif (self.radius.size == x.size):
            self.singlerad = False
        else:
            raise ValueError("Must either have one radius or same as x/y")
        
        
        if (width is None) :
            width = 1
        self.width = width

        self.fk5 = fk5

    def output(self,outfile,clobber=False,append=False):
        if (os.path.isfile(outfile)) :
            if (clobber) :
                print "Note: clobbering existing "+outfile
            else :
                print outfile+" already exists.  Exiting."
                return

        if (not append):
            f=open(outfile,"w")
            #f.write("global color="+self.color+" width=%02d" % self.width+"\n")
        else :
            f=open(outfile,"a")
            
        if (not self.fk5) :
            if (not append):
                f.write("physical\n")
            for i in self.index :
                if (self.uselabels) :
                    if (self.singlerad) :
                        f.write("circle(%.2f,%.2f,%.2f) # color = %s width = %02d text={%d}\n" \
                                    % (self.x[i],self.y[i],self.radius,self.color,self.width,self.labelvals[i]))
                    else :
                        f.write("circle(%.2f,%.2f,%.2f) # color = %s width = %02d text={%d}\n" \
                                    % (self.x[i],self.y[i],self.radius[i],self.color,self.width,self.labelvals[i]))
                else :
                    if (self.singlerad):                        
                        f.write("circle(%.2f,%.2f,%.2f) # color = %s width = %02d\n" % (self.x[i],self.y[i],self.radius,self.color,self.width))
                    else :
                        f.write("circle(%.2f,%.2f,%.2f) # color = %s width = %02d\n" % (self.x[i],self.y[i],self.radius[i],self.color,self.width))

        else :
            if (not append):
                f.write("fk5\n")
            for i in self.index :
                if (self.uselabels) :
                    if (self.singlerad):
                        f.write("circle(%.7f,%.7f,%.2f\") # color = %s width = %02d text={%d}\n" \
                                    % (self.x[i],self.y[i],self.radius,self.color,self.width,self.labelvals[i]))
                    else:
                        f.write("circle(%.7f,%.7f,%.2f\") # color = %s width = %02d text={%d}\n" \
                                    % (self.x[i],self.y[i],self.radius[i],self.color,self.width,self.labelvals[i]))
                else :
                    if (self.singlerad):
                        f.write("circle(%.7f,%.7f,%.2f\") # color = %s width = %02d\n" % (self.x[i],self.y[i],self.radius,self.color,self.width))
                    else :
                        f.write("circle(%.7f,%.7f,%.2f\") # color = %s width = %02d\n" % (self.x[i],self.y[i],self.radius[i],self.color,self.width))

        f.close()
        
