""" Simple script to plot the distributions of the size of the datasets files """

import os,sys
import numpy as np
import matplotlib.pyplot as plt



""" To create the dictionaries with the data: 
#l = [ f for f in os.listdir('.') if 'era5.conv._' in f ]   
#ll = [ i for i in l if '.gz' not in i and '.nc' not in i ]
s = [ os.path.getsize(f) for f in  ll ]
a = { 'l' : s }
#c = np.save('/raid60/scratch/federico/era5_1_histo_dic', a  )




ll = [ f for f in os.listdir('.')  if '.conv.' in f and '.nc' not in f ]
s = [ os.path.getsize(f) for f in  ll ]
a = { 'l' : s }
c = np.save('/raid60/scratch/federico/era5_1759_histo_dic', a  )





                         

ll = [ f for f in os.listdir('.')  if '.conv.' in f and '.nc' not in f ]
s = [ os.path.getsize(f) for f in  ll ]
a = { 'l' : s }
c = np.save('/raid60/scratch/federico/era5_1761_histo_dic', a  )

                                                                                                                                                                                                 

ll = [ f for f in os.listdir('.')  if '.C:' in f and '.nc' not in f ]
s = [ os.path.getsize(f) for f in  ll ]
a = { 'l' : s }
c = np.save('/raid60/scratch/federico/era5_3188_histo_dic', a  )



ll = [ f for f in os.listdir('.')  if '.bfr' in f and '.nc' not in f ]
s = [ os.path.getsize(f) for f in  ll ]
a = { 'l' : s }
c = np.save('/raid60/scratch/federico/bufr_histo_dic', a  )



ll = [ f for f in os.listdir('.')  if '.txt' in f and '.nc' not in f ]
s = [ os.path.getsize(f) for f in  ll ]
a = { 'l' : s }
c = np.save('/raid60/scratch/federico/igra2_histo_dic', a  )




ll = [ f for f in os.listdir('.')  if '.txt' in f and '.nc' not in f ]
s = [ os.path.getsize(f) for f in  ll ]
a = { 'l' : s }
c = np.save('/raid60/scratch/federico/ncar_histo_dic', a  )
"""



ds = ['igra2' , 'era5_1' , 'ncar' , 'bufr' , 'era5_1759' , 'era5_1761' , 'era5_3188']


for d in ds:
    c = np.load( d + '_histo_dic.npy', allow_pickle = True ).item()
    h = c['l']
    h = [ float(j)/10**9 for j in h ]
    plt.hist(h , histtype = 'step' , label = d + ' ' + str(len(h)) )
    plt.xlabel('Dimension [GB]')
    plt.legend()

plt.title('Size of the original source files', y = 1.03 )
plt.savefig('size.png' , dpi = 250)


