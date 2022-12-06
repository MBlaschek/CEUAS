import os,sys
import gzip
import time
from numba import * #njit
import numpy
import datetime

#@njit(cache=True)
def dswitch(itime,idx):
    
#    print(itime)
#    print(idx)
    j=0
    itold=-1
    for i in range(itime.shape[0]):
        if itold!=itime[i]:
            idx[j]=i
            itold=itime[i]
            j+=1
    return j

@njit(cache=True)
def findtabs(x,idxt,idxn,tab,cr):
    c=0
    j=0
    k=0
    d=0
    for i in range(x.shape[0]):
        if x[i]==tab[0]:
            if c==8 or c==9:
                
                idxt[j]=i
                if c==9 and j>2:
                    for m in range(1,i-idxt[j-1]):
                        if x[idxt[j-3]+m]!=x[idxt[j-1]+m]:
                            idxt[d]=k-1
                            d+=1
                            break
                j+=1
            c=c+1
        if x[i]==cr[0]:
            idxn[k]=i
            k+=1
            c=0
    idxt[d]=k-1
    d+=1
    return k,d
        
    
def filter_odbgz(ifile,ofile):

    #countlist=glob.glob(opath+'/*.count')
    alldata=''
    alldict=dict()


    t=time.time()
    try:
        #with gzip.open(os.path.expandvars('$RSCRATCH/era5/odbs/1/era5.conv._10393.gz'), 'rb') as f:
            #rodata = f.readline()   
        with open(ifile, 'rb') as f:
            rodata = f.readline()   
            rdata = f.read()#.split()#
            print('read gzip',time.time()-t)
            rl=len(rdata)
            idxt=numpy.empty(rl//100,dtype=numpy.int)
            idxn=numpy.empty(rl//100,dtype=numpy.int)
            tab=numpy.frombuffer(b'\t',dtype='S1')
            cr=numpy.frombuffer(b'\n',dtype='S1')
            k,d=findtabs(numpy.frombuffer(rdata,dtype='S1'),idxt,idxn,tab,cr)
            idxt=idxt[:d]
            idxn=idxn[:k]
            
    except MemoryError:
        print('odb failed!:'+' '+ifile)
        return
    
    #y=numpy.frombuffer(rdata,dtype='S1')
    #idt=numpy.zeros(len(y)//20,dtype=numpy.int32)
    #idt1=numpy.zeros(len(y)//20,dtype=numpy.int32)
    #idn=numpy.zeros(len(y)//20,dtype=numpy.int32)
    #k=findtabs(rdata,idt,idt1,idn)
    print('read gzip',time.time()-t)
    
    header=rodata.split()
    itime=numpy.empty(d,dtype=numpy.int32)
    idate=numpy.empty(d,dtype=numpy.int32)
    rtype=numpy.empty(d,dtype=numpy.int32)
    #plev=numpy.empty(d)
    tstamp=numpy.empty(d,dtype='datetime64[s]')
    for i in range(d):
        list0=rdata[idxn[idxt[i]-1]:idxn[idxt[i]]].split()
        itime[i]=int(list0[19])
        idate[i]=int(list0[18])
        #plev[i]=float(list0[37])
        #tstamp.append(datetime.datetime.strptime((list0[18]+list0[19]).decode(),'%Y%m%d%H%M%S'))
        tstamp[i]=numpy.datetime64((list0[18][:4]+b'-'+list0[18][4:6]+b'-'+list0[18][6:]).decode()+\
                                   'T{:0>2}:{:0>2}:{:0>2}'.format(itime[i]//10000,(itime[i]%10000)//100,itime[i]%100))
        rtype[i]=list0[6]
    #print(time.time()-t)
    #iseqno=numpy.int32(rlist[9::len(header)])
    
    #itime=numpy.int32(rlist[19::len(header)])
#    idate=numpy.int32(rdata[18::len(header)])
#    adate=numpy.int32(rdata[4::len(header)])
#    atime=numpy.int32(rdata[5::len(header)])
    #rtype=numpy.int32(rlist[6::len(header)])
    
    
    print('preanalyze',time.time()-t)
    wdata=[]
    idxnold=0
    diffs=numpy.empty_like(idxt)
    diffs[1:]=idxt[1:]-idxt[:-1]
    diffs[0]=idxt[0]
    
    try:
        with open(ofile, 'wb') as f:
            f.write(rodata)
            j=0
            for i in range(d):
#                idup=numpy.where(numpy.logical_and(idate[i]==idate,itime[i]==itime))[0]
                idup=numpy.where(abs((tstamp[i]-tstamp)/numpy.timedelta64(1, 's'))<7200)[0]
                if diffs[i]==numpy.max(diffs[idup]):
                    if i>0:
                        wdata=rdata[idxn[idxt[i-1]]+1:idxn[idxt[i]]+1]
                    else:
                        wdata=rdata[:idxn[idxt[i]]+1]
                    f.write(wdata)
                    j+=1
                    #f.write(b'\n'.join(wdata))
        #            wdata.append(rdata[idxnold:idxn[idxt[i]]])
                    
    except MemoryError:
        print('odb failed!:'+' '+ifile)
    #if itime.shape[0]<10000:
        #idx=numpy.zeros(itime.shape[0],dtype=numpy.int32)
    #else:
        #idx=numpy.zeros(itime.shape[0]//5,dtype=numpy.int32)
    #j=dswitch(itime, idx)
    #idx=idx[:j]
    #diff=numpy.empty_like(idx)
    #diff[:-1]=idx[1:]-idx[:-1]
    #diff[-1]=iseqno.shape[0]-idx[-1]
    #j=dswitch(iseqno, idx)
    #idx=numpy.concatenate([idx,numpy.array([len(rdata)//63],dtype=numpy.int32)])

    #allgood=[]
    #tgood=[]
    #goodl=[]
    #good=[]
    #atold=-1
    #adold=-1
    #for i in range(idx.shape[0]-1):
        #if atold!=atime[idx[i]] or adold!=adate[idx[i]]:
            #atold=atime[idx[i]]
            #adold=adate[idx[i]]
            #for it in range(len(tgood)-1,-1,-1):
                #diffs=itime[tgood[it][0]]-itime[good]
                #diffs[diffs<-120000]+=240000
                #if any(abs(diffs)<20000):
                    #del tgood[it]
            #allgood=allgood+tgood+goodl
            #good=[]
            #goodl=[]
            #tgood=[]
        #if rtype[idx[i]] in [16045,16068] and itime[idx[i]] not in itime[good]:
            #goodl.append([idx[i],idx[i+1]])
            #good.append(idx[i])
        #if rtype[idx[i]] in [16022,16013]:
            #tgood.append([idx[i],idx[i+1]])
            
    #for it in range(len(tgood)-1,-1,-1):
        #diffs=itime[tgood[it][0]]-itime[good]
        #diffs[diffs<-120000]+=240000
        #if any(abs(diffs)<20000):
            #del tgood[it]
    #allgood=allgood+tgood+goodl
    
    #print('analyze',len(wdata),time.time()-t)
   
    
    #try:
        #with open(ofile, 'wb') as f:
            #f.writelines([rodata])
            #f.write(b'\n'.join(wdata))
            ##for i in range(len(allgood)):
                ##chunk=rdata[63*allgood[i][0]:63*allgood[i][1]]
                ##for j in range(0,len(chunk),63):
                    ##f.writelines([b'\t'.join(chunk[j:j+63])+b'\n'])
    #except MemoryError:
        #print('odb failed!:'+' '+ifile)
        
    print('wrote',ofile,',',j,'out of',d,'records',time.time()-t)
    return
           
    
                      
    
if __name__ == '__main__':

    print(sys.argv[1])
    filter_odbgz(sys.argv[1]+'.txt',sys.argv[1]+'.txt')
    #filter_odbgz(os.path.expandvars('$RSCRATCH/era5/odbs/1/era5.conv.201901.10393.txt'), 
                 #os.path.expandvars('$RSCRATCH/era5/odbs/1/era5.conv.201901.10393._txt'))
    #filter_odbgz(os.path.expandvars('$RSCRATCH/era5/odbs/1/era5.conv.201901.10393._txt'), 
                 #os.path.expandvars('$RSCRATCH/era5/odbs/1/era5.conv.201901.10393.__txt'))
    
    print('ready')