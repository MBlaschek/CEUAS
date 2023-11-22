import numpy

#@njit(cache=True)
# pythran export esplit(float32[][][][],float32[][][][],int,float32,float32[][][][][],int32[],int)
def esplit(tgrid,fv,offset,missval,corr,tidx,im):

    vs=fv.shape
    if im==-1:
        for i in range(offset,vs[0]*2,2):
            for j in range(vs[1]):
                for k in range(vs[2]):
                    for l in range(vs[3]):
                        if fv[i/2,j,k,l]!=missval:
                            tgrid[i,vs[1]-j-1,k,l]=fv[i/2,j,k,l]
                        else:
                            tgrid[i,vs[1]-j-1,k,l]=numpy.nan
    else:
        for i in range(offset,vs[0]*2,2):
            for j in range(tidx.shape[0]):
                if i/2>=tidx[j]:
                    ik=j
            for j in range(vs[1]):
                for k in range(vs[2]):
                    for l in range(vs[3]):
                        if fv[i/2,j,k,l]!=missval:
                            tgrid[i,vs[1]-j-1,k,l]=fv[i/2,j,k,l]-corr[im+ik,offset,vs[1]-j-1,k,l]
                        else:
                            tgrid[i,vs[1]-j-1,k,l]=numpy.nan

    return

