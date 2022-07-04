import os,sys

if len(sys.argv)>1:
    print sys.argv[1]


ifname='anomaly.py'
#f=open(sys.argv[1],'r')
f=open(ifname,'r')
lines=f.readlines()
imports=[]
sfile=[]
for l in lines:
    if 'import' in l:
        imports.append(l)
    if 'pythran' in l :
        sub=l.split('(')
        sub=sub[0].split()
        
        sfile.append(sub[3])

dname=ifname.split('.')[0]
try:
    os.mkdir(dname)
except:
    pass
j=0
sf='xxx'
for l in lines:
    
    if 'pythran' in l and sfile[j] in l:
        fo=open(dname+'/'+sfile[j]+'.py','w')
        fo.writelines(imports)
#        fo.writelines(['\n','# pythran export '+sfile[j]+'()\n','\n'])
        fo.writelines([l])
        lsave=l
        sf=sfile[j]
        j=j+1
    elif '@njit' in l:
        pass
    elif 'import' in l:
        pass
    elif 'def' in l and sf in l:
        t1=l.split(',')
        t2=lsave.split(',')
        if len(t1)!=len(t2):
            print l,lsave 
            print 'argument mismatch'
            exit()
        else:
            print sf+': argument match'
        fo.writelines([l])
        fo.writelines(['\n','    miss_val=-1.e30\n','\n'])
    elif j>0:
        fo.writelines([l])


pass
