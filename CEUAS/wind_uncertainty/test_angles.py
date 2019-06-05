""" Script to test if angles extracted from u,v components are correct """


import numpy as np


# x components
x = [0, 0, 1, 1, 1, -1, -1, -1 , 1000, 1.0001]
y = [1,-1, 0, 1,-1,  0,  1, -1 , 1000, 0]





x = [-1,  1 , 1 ]
y = [-1, -1 , 1 ]
z = np.arctan2(y,x)


zz = 270 -(180/np.pi)*np.arctan2(y,x)

print ('zz is', zz)
for a,b,c,d in zip (x,y,z,zz):
    if c < 0:
       c  = 2*np.pi - abs(c)

    ang_deg = 360*c/(2*np.pi) # just converts radians to degree

    if d > 360:
        d = d -360 
    print('components u and v : ', b, a,  'arctan2: ', c, ' angle degree: ',  ang_deg, ' old: ', d)





 
