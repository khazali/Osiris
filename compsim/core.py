import math
import scipy.integrate as integrate

def function1(x, ystart, yend, xmean, R):
    if x<xmean:
        return math.sqrt(R*R-x*x)-ystart
    else:
        return yend-ystart


def function2(x, ystart, R):
    return math.sqrt(R*R-x*x)-ystart


def function3(x, ystart, yend, xmean, R):
    if x<xmean:        
        return math.sqrt(R*R-x*x)-ystart
    else:
        return 0


def function4(x, ystart, yend, xmean, R):
    if x<xmean:        
        return yend-ystart
    else:
        return math.sqrt(R*R-x*x)-ystart
        

def function5(x, ystart, yend, xmean1, xmean2, R):
    if xmean1<xmean2:
        if x<xmean1:
            return yend-ystart
        elif x<xmean2:
            return math.sqrt(R*R-x*x)-ystart
        else:
            return 0
    else:
        if x<xmean2:
            return 0
        elif x<xmean1:
            return math.sqrt(R*R-x*x)-ystart
        else:
            return yend-ystart


def function6(x, ystart, yend, xmean, R):
    if x<xmean:        
        return 0
    else:
        return math.sqrt(R*R-x*x)-ystart


def function7(x, ystart, yend, xmean1, xmean2, xmean3, xmean4, R):
    if xmean1>xmean4:
        temp=xmean1
        xmean1=xmean4
        xmean4=temp

    if xmean2>xmean3:
        temp=xmean2
        xmean2=xmean3
        xmean3=temp

    if x<xmean1:        
        return 0
    elif x<xmean2:
        return math.sqrt(R*R-x*x)-ystart
    elif x<xmean3:
        return yend-ystart
    elif x<xmean4:
        return math.sqrt(R*R-x*x)-ystart
    else:
        return 0


def function8(x, ystart, yend, xmean1, xmean2, xmean3, R):
    if xmean2>xmean3:
        temp=xmean2
        xmean2=xmean3
        xmean3=temp

    if x<xmean1:        
        return 0
    elif x<xmean2:
        return math.sqrt(R*R-x*x)-ystart
    elif x<xmean3:
        return yend-ystart
    else:
        return math.sqrt(R*R-x*x)-ystart


def function9(x, ystart, yend, xmean1, xmean2, xmean3, R):
    if xmean1>xmean2:
        temp=xmean2
        xmean2=xmean1
        xmean1=temp

    if x<xmean1:
        return math.sqrt(R*R-x*x)-ystart        
    elif x<xmean2:
        return yend-ystart        
    elif x<xmean3:
        return math.sqrt(R*R-x*x)-ystart
    else:
        return 0

def function10(x, ystart, yend, xmean1, xmean2, R):
    if xmean1>xmean2:
        temp=xmean2
        xmean2=xmean1
        xmean1=temp

    if x<xmean1:
        return math.sqrt(R*R-x*x)-ystart        
    elif x<xmean2:
        return yend-ystart        
    else:
        return math.sqrt(R*R-x*x)-ystart


def function11(x, ystart, yend, xmean1, xmean2, R):
    if xmean1>xmean2:
        temp=xmean2
        xmean2=xmean1
        xmean1=temp

    if x<xmean1:
        return 0
    elif x<xmean2:
        return math.sqrt(R*R-x*x)-ystart
    else:
        return 0





class Grid:
    x=[0, 0]
    y=[0, 0]
    z=[0, 0]

    
    
    Permeability=0
    Porosity=0
    R=0

    mynone=-5e30
       

    def DetermineInsideStat(self):
        if self.x[0]>self.x[1]:
            temp=self.x[0]
            self.x[0]=self.x[1]
            self.x[1]=temp

        if self.y[0]>self.y[1]:
            temp=self.y[0]
            self.y[0]=self.y[1]
            self.y[1]=temp

        if self.z[0]>self.z[1]:
            temp=self.z[0]
            self.z[0]=self.z[1]
            self.z[1]=temp
        

        ExpectedY=[0, 0]
        ExpectedX=[0, 0, 0, 0]
        ConvertedX=[0, 0]
        ConvertedY=[0, 0]
        side=[False, False, False, False, False, False]
        integ=[0, 0]
        

        for i in range(0, 2):
            ConvertedX[i]=self.z[i]-self.R
            ConvertedY[i]=self.y[i]
            
            if abs(ConvertedX[i])<=self.R:
                ExpectedY[i]=math.sqrt(self.R*self.R-ConvertedX[i]*ConvertedX[i])
            else:
                ExpectedY[i]=self.mynone

            if abs(ConvertedY[i])<=self.R:
                ExpectedX[i]=math.sqrt(self.R*self.R-ConvertedY[i]*ConvertedY[i])
                ExpectedX[i+2]=-math.sqrt(self.R*self.R-ConvertedY[i]*ConvertedY[i])
            else:
                ExpectedX[i]=self.mynone
                ExpectedX[i+2]=self.mynone

        area=(ConvertedX[1]-ConvertedX[0])*(ConvertedY[1]-ConvertedY[0])

        if (ExpectedY[0]>(self.mynone)) and (ExpectedY[0]>=ConvertedY[0]) and (ExpectedY[0]<=ConvertedY[1]):
            side[0]=True
        if (ExpectedY[1]>(self.mynone)) and (ExpectedY[1]>=ConvertedY[0]) and (ExpectedY[1]<=ConvertedY[1]):
            side[2]=True
        if (ExpectedY[0]>(self.mynone)) and (ExpectedY[1]>(self.mynone)) and (ExpectedY[0]>=ConvertedY[1]) and (ExpectedY[1]>=ConvertedY[1]):
            integ[0]=area
            

       
        if ((ExpectedX[1]>(self.mynone)) and (ExpectedX[3]>(self.mynone))) and (((ExpectedX[1]>=ConvertedX[0]) and (ExpectedX[1]<ConvertedX[1])) and ((ExpectedX[3]>=ConvertedX[0]) and (ExpectedX[3]<ConvertedX[1]))):
            side[4]=True
        elif (ExpectedX[3]>(self.mynone)) and (ExpectedX[3]>=ConvertedX[0]) and (ExpectedX[3]<ConvertedX[1]):
            side[1]=True
            xx1=ExpectedX[3]
        elif (ExpectedX[1]>(self.mynone)) and (ExpectedX[1]>=ConvertedX[0]) and (ExpectedX[1]<ConvertedX[1]):
            side[1]=True
            xx1=ExpectedX[1]

        if ((ExpectedX[0]>(self.mynone)) and (ExpectedX[2]>(self.mynone))) and (((ExpectedX[0]>=ConvertedX[0]) and (ExpectedX[0]<ConvertedX[1])) and ((ExpectedX[2]>=ConvertedX[0]) and (ExpectedX[2]<ConvertedX[1]))):
            side[5]=True
        elif (ExpectedX[0]>(self.mynone)) and (ExpectedX[0]>=ConvertedX[0]) and (ExpectedX[0]<ConvertedX[1]):
            side[3]=True
            xx0=ExpectedX[0]
        elif (ExpectedX[2]>(self.mynone)) and (ExpectedX[2]>=ConvertedX[0]) and (ExpectedX[2]<ConvertedX[1]):
            side[3]=True
            xx0=ExpectedX[2]
        
        
        if side[0]:
            if side[1]:                
                integ=integrate.quad(function1, ConvertedX[0], ConvertedX[1], args=(ConvertedY[0], ConvertedY[1], xx1, self.R))
                
            elif side[2]:
                integ=integrate.quad(function2, ConvertedX[0], ConvertedX[1], args=(ConvertedY[0], self.R))               
            
            elif side[3]:
                integ=integrate.quad(function3, ConvertedX[0], ConvertedX[1], args=(ConvertedY[0], ConvertedY[1], xx0, self.R))

        elif side[1]:
            if side[2]:
                integ=integrate.quad(function4, ConvertedX[0], ConvertedX[1], args=(ConvertedY[0], ConvertedY[1], xx1, self.R))
            elif side[3]:
                integ=integrate.quad(function5, ConvertedX[0], ConvertedX[1], args=(ConvertedY[0], ConvertedY[1], xx1, xx0, self.R))
               
        elif side[2]:
            if side[3]:
                integ=integrate.quad(function6, ConvertedX[0], ConvertedX[1], args=(ConvertedY[0], ConvertedY[1], xx0, self.R))
        
        
        elif side[4]:
            if side[5]:
                integ=integrate.quad(function7, ConvertedX[0], ConvertedX[1], args=(ConvertedY[0], ConvertedY[1], ExpectedX[0], ExpectedX[1], ExpectedX[3], ExpectedX[2], self.R))
            
            elif side[2] and side[3]:
                integ=integrate.quad(function8, ConvertedX[0], ConvertedX[1], args=(ConvertedY[0], ConvertedY[1], xx0, ExpectedX[1], ExpectedX[3], self.R))

            elif side[0] and side[3]:
                integ=integrate.quad(function9, ConvertedX[0], ConvertedX[1], args=(ConvertedY[0], ConvertedY[1], xx0, ExpectedX[1], ExpectedX[3], self.R))

            elif side[0] and side[3]:
                integ=integrate.quad(function9, ConvertedX[0], ConvertedX[1], args=(ConvertedY[0], ConvertedY[1], ExpectedX[1], ExpectedX[3], xx0, self.R))

            elif side[0] and side[2]:
                integ=integrate.quad(function10, ConvertedX[0], ConvertedX[1], args=(ConvertedY[0], ConvertedY[1], ExpectedX[1], ExpectedX[3], self.R))

        elif side[5] and (not side[0]) and (not side[1]) and (not side[2]) and (not side[3]) and (not side[4]):
            integ=integrate.quad(function11, ConvertedX[0], ConvertedX[1], args=(ConvertedY[0], ConvertedY[1], ExpectedX[0], ExpectedX[2], self.R))




        
        if (integ[0]/area)<0.5:
            self.Porosity=0
            self.Permeability=0
        
            






########################################################################################################
R=38e-3/2
Lx=85.5e-3
Ly=R
Lz=2*R
#Ly=2*R
#Lz=4*R

Nx=81
Ny=20
Nz=23

Hx=Lx/Nx
Hy=Ly/Ny
Hz=Lz/Nz

Fperm=3.5
Mperm=3

#Fpor=0.03242
Fpor=1
Mpor=46.5e-2
#Mpor=1
#ZeroPore=0.01

fp=open('e:\\filefile.txt', 'w')

fp.write('GRID '+repr(Nx)+' '+repr(Ny)+' '+repr(Nz)+'\n\n')
fp.write('DI CON '+repr(Hx)+'\n')
fp.write('DJ CON '+repr(Hy)+'\n')
fp.write('DK CON '+repr(Hz)+'\n')

Mesh = [[[Grid() for k in range(Nz)] for j in range(Ny)] for i in range(Nx)]

fp.write('\n\nPOR VAR\n')
for k in range(Nz):
    for j in range(Ny):
        for i in range (Nx):
            (Mesh[i][j][k]).x[0]=i*Hx
            (Mesh[i][j][k]).x[1]=(i+1)*Hx
            (Mesh[i][j][k]).y[0]=j*Hy
            (Mesh[i][j][k]).y[1]=(j+1)*Hy
            (Mesh[i][j][k]).z[0]=k*Hz
            (Mesh[i][j][k]).z[1]=(k+1)*Hz
            (Mesh[i][j][k]).R=R
            (Mesh[i][j][k]).Permeability=Mperm
            (Mesh[i][j][k]).Porosity=Mpor
            (Mesh[i][j][k]).DetermineInsideStat()

            str=repr((Mesh[i][j][k]).Porosity)
            fp.write(str+ ' ')
        fp.write('\n')
    fp.write('\n')

fp.write('\n\nPERMI VAR\n')
for k in range(Nz):
    for j in range(Ny):
        for i in range (Nx):
            str=repr((Mesh[i][j][k]).Permeability)
            fp.write(str+ ' ')
        fp.write('\n')
    fp.write('\n')

fp.write('\n\nPERMJ VAR\n')
for k in range(Nz):
    for j in range(Ny):
        for i in range (Nx):
            str=repr((Mesh[i][j][k]).Permeability)
            fp.write(str+ ' ')
        fp.write('\n')
    fp.write('\n')

fp.write('\n\nPERMK VAR\n')
for k in range(Nz):
    for j in range(Ny):
        for i in range (Nx):
            str=repr((Mesh[i][j][k]).Permeability)
            fp.write(str+ ' ')
        fp.write('\n')
    fp.write('\n')

fp.write('\n\nFPOR VAR\n')
for k in range(Nz):
    for j in range(Ny):
        for i in range (Nx):
#            (Mesh[i][j][k]).x[0]=i*Hx
#            (Mesh[i][j][k]).x[1]=(i+1)*Hx
#            (Mesh[i][j][k]).y[0]=j*Hy
#            (Mesh[i][j][k]).y[1]=(j+1)*Hy
#            (Mesh[i][j][k]).z[0]=k*Hz
#            (Mesh[i][j][k]).z[1]=(k+1)*Hz
#            (Mesh[i][j][k]).R=R
#            (Mesh[i][j][k]).Permeability=Fperm
#            (Mesh[i][j][k]).Porosity=Fpor
#            (Mesh[i][j][k]).DetermineInsideStat()

            if (i==0) or (i==(Nx-1)):
#                str=repr((Mesh[i][j][k]).Porosity/2)
                str=repr(Fpor)
            elif (j==0):
                str=repr(Fpor)
            else:
                str=repr(0)
            fp.write(str+ ' ')
        fp.write('\n')
    fp.write('\n')

fp.write('\n\nFPERMI VAR\n')
for k in range(Nz):
    for j in range(Ny):
        for i in range (Nx):
            if (i==0) or (i==(Nx-1)) or (j==0):
                str=repr(Fperm)
            else:
                str=repr(0)
            fp.write(str+ ' ')
        fp.write('\n')
    fp.write('\n')

fp.write('\n\nFPERMJ VAR\n')
for k in range(Nz):
    for j in range(Ny):
        for i in range (Nx):
            if (i==0) or (i==(Nx-1)) or (j==0):
                str=repr(Fperm)
            else:
                str=repr(0)
            fp.write(str+ ' ')
        fp.write('\n')
    fp.write('\n')

fp.write('\n\nFPERMK VAR\n')
for k in range(Nz):
    for j in range(Ny):
        for i in range (Nx):
            if (i==0) or (i==(Nx-1)) or (j==0):
                str=repr(Fperm)
            else:
                str=repr(0)
            fp.write(str+ ' ')
        fp.write('\n')
    fp.write('\n')

fp.close()
