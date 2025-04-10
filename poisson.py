import numpy
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D
%matplotlib inline
# Parameters
nx = 50
ny = 50
nt  = 200000
xmin = 0
xmax = numpy.pi
ymin = 0
ymax = numpy.pi

dx = (xmax - xmin) / (nx - 1)
dy = (ymax - ymin) / (ny - 1)

# Initialization
p  = numpy.zeros((ny, nx))
pd = numpy.zeros((ny, nx))
b  = numpy.zeros((ny, nx))
x  = numpy.linspace(xmin, xmax, nx)
y  = numpy.linspace(xmin, xmax, ny)
X, Y = numpy.meshgrid(x, y)

# Source
cx=8
cy=8
b[int(ny / 4), int(nx / 4)]  = 100
b[int(3 * ny / 4), int(3 * nx / 4)] = -100
b=2.*((Y*Y-Y)+(X*X-X))
pexact=(X-X*X)*(Y-Y*Y)
pexact=numpy.sin(cx*X)*numpy.sin(cy*Y)
b=-1.*(cx*cx*pexact+cy*cy*pexact)
for it in range(nt):

    pd = p.copy()

    p[1:-1,1:-1] = (((pd[1:-1, 2:] + pd[1:-1, :-2]) * dy**2 +
                    (pd[2:, 1:-1] + pd[:-2, 1:-1]) * dx**2 -
                    b[1:-1, 1:-1] * dx**2 * dy**2) / 
                    (2 * (dx**2 + dy**2)))

    p[0, :] = 0
    p[ny-1, :] = 0
    p[:, 0] = 0
    p[:, nx-1] = 0
    
err=p-pexact
totalerr=sum(err)
meanerr=numpy.mean(err)
maxerr=numpy.max(err)
maxp=numpy.max(p)
maxpexact=numpy.max(pexact)
print(maxerr,meanerr,maxp,maxpexact)

