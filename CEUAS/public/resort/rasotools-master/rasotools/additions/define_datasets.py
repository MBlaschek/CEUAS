import numpy
import copy
import os

def make_andeplist():
    andeplist=['jra55_andep','jra55_fgdep','e20c_andep','n20c_andep','erapresat_andep','bgdep',
               'eijra_fgdep','eice20c_andep','jrace20c_andep','fg_dep','an_dep',
               'comp-bg','comp2-bg','comp4-bg', 'era5_andep','era5_fgdep','e5_fg_dep','e5_an_dep','e5_bias','e5_3116_bias','oper_bias',
               'oper_fg_dep','jracepre_fgdep','jrapresat_andep',
     'ce20c0_andep','eice20c0_andep','ce20c1_andep','ce20c2_andep','ce20c3_andep','ce20c4_andep','ce20c5_andep',
     'ce20c6_andep','ce20c7_andep','ce20c8_andep','ce20c9_andep','ce20cens_spread','ce20cens_mean','erai_fggpsdep','erai_fggpswetdep']
    andepstdlist=[]
    andeprmslist=[]
    for x in andeplist:
        andepstdlist.append(x+'_std')
    for x in andeplist:
        andeprmslist.append(x+'_rms')
    andepstdlist.append('ce20cens_spread')    
    return andeplist,andepstdlist,andeprmslist   

def make_e5difflist():
    e5difflist=['e5_ei','e5_2502_ei','e5_2504_ei','e5_2506_ei','e5_oper_ei','e5_presat_ei','e5_1759_ei','e5_1770_ei']
    e5difflist+=['e5obs_ei','e5obs_2502_ei','e5obs_2504_ei','e5obs_2506_ei',
                'e5obs_oper_ei','e5obs_presat_ei','e5obs_1759_ei','e5obs_1770_ei']
    
def define_ds(dslist,rfpar,plotproperties,stnames,path):
    
    pathtrunk='/'.join(path.split('/')[:path.split('/').index('rise')])+'/'
    initial_size=plotproperties['initial_size']
    pindex=plotproperties['pindex']
    jpindex=plotproperties['jpindex']
    msupindex=plotproperties['msupindex']
    
    andeplist,andepstdlist,andeprmslist=make_andeplist()
    
    dailyfile='feedbackmerged'
    dailyobsvar=['temperatures']
    dailysuff=['']
    for k in range(len(plotproperties['parameters'])):
        plotproperties['parameters'][k]=str(plotproperties['parameters'][k])
    if plotproperties['parameters'][0]=='temperatures' and plotproperties['exps'][0] in ('3188','1759','1761','ai_bfr','igra','erai'):
        dailyfile='ERA5_'+plotproperties['exps'][0]+'_'
        #if plotproperties['exps'][0]=='odbei':
            #dailyfile=''
        dailysuff=['_t']
    if plotproperties['parameters'][0]=='uwind':
        dailysuff=['_u']
        dailyobsvar=['uwind']
    if plotproperties['parameters'][0]=='vwind':
        dailysuff=['_v']
        dailyobsvar=['vwind']
    if plotproperties['parameters'][0]=='wspeed':
        dailysuff=['_u','_v']
        dailyobsvar=['uwind','vwind']
    if plotproperties['parameters'][0]=='wdir':
        dailysuff=['_u','_v']
        dailyobsvar=['uwind','vwind']
    
        
    def make_ddict(
        shortname='tm',
        name='unadjusted',
        file="feedbackglobbincorrsave",
        #msufile="rttov_",
        #msufile="feedbackglobbincorrmon",
        msufile="feedbackglobbinmon",
        var="rasocorr",
        suff="",    
        dvar="temperatures",
        msuvar="montemp",
        ens=[0],
        dfile=['feedbackmerged'],
        dsuff=[''],
        fpath=[''],
        startdate=rfpar["startdate"],
        color="blue",
        index=numpy.zeros([initial_size,1,rfpar['brmax']],numpy.int32),
        dindex=[0],
        maxindex=numpy.zeros(initial_size,numpy.int32),
        minindex=numpy.zeros(initial_size,numpy.int32),
        data=[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['brmax']],
        msudata=[initial_size,1,rfpar['parmax'],4,rfpar['mmax']],
        ddata=[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
        zslopes=numpy.empty([1,18,len(jpindex)]),
        beltslopes=numpy.empty([1,6,len(jpindex)]),
        beltanomalies=numpy.empty([1,6,len(jpindex),rfpar['mmax']]),
        type='REA'
        ):
        
        x=locals()
        return x

    dsdict=[]
    suff=""
    l=0
    for ds in dslist:
        if ds=='tm':
            dsdict.append(make_ddict(type='ROB'))
        elif ds in ['uah','rss','star']:
            name={'uah':['UAH MSU brightness temperatures','orange'],
                      'rss':['RSS MSU brightness temperatures','red'],
                      'star':['STAR MSU brightness temperatures','orange_yellow']}
            
            dsdict.append(make_ddict(**{"shortname":ds,
                "name":name[ds][0],
                "file":"",
                "var":"",
                "suff":"",    
                "index":[0],
                'data':[],
                'color':name[ds][1],
                #'msudata':[initial_size,1,rfpar['parmax'],len(msupindex),rfpar['mmax']],
                #'zslopes':numpy.empty([1,18,len(msupindex)]),
                #'beltslopes':numpy.empty([1,6,len(msupindex)]),
                #'beltanomalies':numpy.empty([1,6,len(msupindex),rfpar['mmax']]),
                'msudata':[initial_size,1,rfpar['parmax'],4,rfpar['mmax']],
                'zslopes':numpy.empty([1,18,4]),
                'beltslopes':numpy.empty([1,6,4]),
                'beltanomalies':numpy.empty([1,6,4,rfpar['mmax']]),
                'type':'SAT'
                }))

        elif ds in ['wegc']:
            name={'wegc':['Wegener Center GPS dry temperatures','charcoal']}
            
            dsdict.append(make_ddict(**{"shortname":ds,
                "name":name[ds][0],
                "file":"",
                "var":"",
                "suff":"",    
                "index":[0],
                'data':[],
                'msudata':[initial_size,1,rfpar['parmax'],4,rfpar['mmax']],
                'color':name[ds][1],
                'zslopes':numpy.empty([1,18,len(jpindex)]),
                'beltslopes':numpy.empty([1,6,len(jpindex)]),
                'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']]),
                'type':'SAT'
                }))

        elif ds in ['merra2']:
            name={'merra2':['MERRA reanalysis','turquoise']}
            
            dsdict.append(make_ddict(**{"shortname":ds,
                "name":name[ds][0],
                "file":"",
                "var":"",
                "suff":"",    
                "index":[0],
                'data':[],
                'msudata':[initial_size,1,rfpar['parmax'],4,rfpar['mmax']],
                'color':name[ds][1],
                'zslopes':numpy.empty([1,18,len(jpindex)]),
                'beltslopes':numpy.empty([1,6,len(jpindex)]),
                'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']]),
                'type':'REA'
                }))

        elif ds in ['20CRv3']:
            name={'20CRv3':['20CRv3','turquoise']}
            
            dsdict.append(make_ddict(**{"shortname":ds,
                "name":name[ds][0],
                "file":"",
                "var":"",
                "suff":"",    
                "index":[0],
                'data':[],
                'msudata':[initial_size,1,rfpar['parmax'],4,rfpar['mmax']],
                'color':name[ds][1],
                'zslopes':numpy.empty([1,18,len(jpindex)]),
                'beltslopes':numpy.empty([1,6,len(jpindex)]),
                'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']]),
                'type':'REA'
                }))

        elif ds in ['hadatx']:
            name={'hadat':['HadAT surface temperatures','black']}
            
            dsdict.append(make_ddict(**{"shortname":ds,
                "name":name[ds][0],
                "file":"feedbackglobbincorrmon",
                "var":"montemp",
                "suff":"",    
                "index":[0],
                'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                'color':name[ds][1],
                'msudata':[initial_size,1,rfpar['parmax'],len(msupindex),rfpar['mmax']],
                'zslopes':numpy.empty([1,18,len(jpindex)]),
                'beltslopes':numpy.empty([1,6,len(jpindex)]),
                'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']]),
                'type':'SAT'
                }))

        elif ds in ['hadcrut5','hadcrut4','hadcrutem4']:
            name={ds:[ds.upper()+' surface temperatures','black']}
            if ds=='hadcrutem4':
                name[ds][1]='tan'
            
            dsdict.append(make_ddict(**{"shortname":ds,
                "name":name[ds][0],
                "file":"feedbackglobbincorrmon",
                "var":"montemp",
                "suff":"",    
                "index":[0],
                'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                'color':name[ds][1],
                'msudata':[initial_size,1,rfpar['parmax'],4,rfpar['mmax']],
                'zslopes':numpy.empty([1,18,len(jpindex)]),
                'beltslopes':numpy.empty([1,6,len(jpindex)]),
                'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']]),
                'type':'SAT'
                }))

        elif ds=='tmcorr':
                dsdict.append(make_ddict(**{"shortname":ds,
                    "name":'RAOBCORE-adjusted',
                    "file":"feedbackglobbincorrmon",
                    "msufile":"feedbackglobbincorrmon",
                    #"msufile":"RAOBCORE_rttov_",
                    "dfile":["feedbackglobbincorrsave"],
                    "var":"montemp",
                    "dvar":"rasocorr",
                    "msuvar":"montemp",
                    "suff":"",    
                    "dsuff":[''],
                    "color":"cyan",
                    "index":[0],
                    "dindex":numpy.zeros([len(stnames),1,rfpar['brmax']],numpy.int32),
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['brmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),
                    'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']]),
                    'type':'ROB'
                   }))

        elif ds in ['rharm_h','radrharm_h']:
                dsdict.append(make_ddict(**{"shortname":ds,
                    "name":'RHARM-adjusted',
                    "file":"rharm_h_mon_",
                    "msufile":"rharm_h",
                    "dfile":["rharm_h_"],
                    "var":"temperatures_h",
                    "dvar":"ta_h",
                    "msuvar":"montemp",
                    "suff":"",    
                    "dsuff":[''],
                    "color":"burgundy",
                    "index":[0],
                    "dindex":[0],
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),
                    'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']]),
                    'type':'ROB'
                   }))

        elif ds in ['rharmbc','radrharmbc']:
                dsdict.append(make_ddict(**{"shortname":ds,
                    "name":'RHARM-adjustment',
                    "file":"rharm_h_mon_",
                    "msufile":"rharm_h_mon_",
                    "dfile":["rharm_h_"],
                    "var":"rharmbc",
                    "dvar":"rharmbc",
                    "msuvar":"rharmbc",
                    "suff":"",    
                    "dsuff":[''],
                    "color":"cyan",
                    "index":[0],
                    "dindex":[0],
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),
                    'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']]),
                    'type':'ROB'
                   }))

        elif ds in ['igra','rharm','radrharm']:
                dsdict.append(make_ddict(**{"shortname":ds,
                    "name":'RHARM-IGRA',
                    "file":"rharm_h_mon_",
                    "msufile":"IH_rttov_",
                    "dfile":["rharm_h_"],
                    "var":"temperatures",
                    "dvar":"ta",
                    "msuvar":"temperatures",
                    "suff":"",    
                    "dsuff":[''],
                    "color":"cyan",
                    "index":[0],
                    "dindex":[0],
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),
                    'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']]),
                    'type':'ROB'
                   }))

        elif ds in ('uwind','vwind'):
                dsdict.append(make_ddict(**{"shortname":ds,
                    "name":'RAOBCORE-adjusted',
                    "file":ds+"mon_u",
                    "msufile":"feedbackglobbincorrmon",
                    "dfile":["feedbackglobbincorrsave"],
                    "var":ds,
                    "dvar":"rasocorr",
                    "msuvar":"montemp",
                    "suff":"",    
                    "dsuff":[''],
                    "color":"cyan",
                    "index":[0],
                    "dindex":numpy.zeros([len(stnames),1,rfpar['brmax']],numpy.int32),
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['brmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),
                    'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']]),
                    'type':'ROB'
                   }))

        elif ds=='ctmcorr':
                dsdict.append(make_ddict(**{"shortname":ds,
                    "name":'RAOBCORE-adjusted',
                    "file":"feedbackglobbincorrmon",
                    "fpath":pathtrunk+"rise/1.0/Control/",
                    "msufile":"feedbackglobbincorrmon",
                    "dfile":["feedbackglobbincorrsave"],
                    "var":"montemp",
                    "dvar":"rasocorr",
                    "msuvar":"montemp",
                    "suff":"",    
                    "dsuff":[''],
                    "color":"turquoise",
                    "index":[0],
                    "dindex":numpy.zeros([len(stnames),1,rfpar['brmax']],numpy.int32),
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['brmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),
                    'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']]),
                    'type':'ROB'
                   }))

        elif ds=='milancorr':
                dsdict.append(make_ddict(**{"shortname":ds,
                    "name":'VPROF-adjusted',
                    "file":"feedbackglobbincorrmon",
                    "dfile":["feedbackglobbincorrsave"],
                    "var":"rasocorrmon",
                    "dvar":"rasocorr",
                    "suff":"",    
                    "dsuff":[''],
                    "fpath":'/fio/srvx7/leo/fastscratch/sabiner/datidl/leoformatclstatmonthangnum',
                #    "fpath":'/fio/srvx7/leo/fastscratch/sabiner/datidl/leoformatstat',
                #    "fpath":'/fio/srvx7/leo/fastscratch/sabiner/datidl/leoformatst',
                   #"fpath":'/fio/srvx7/leo/fastscratch/sabiner/datidl/leoformatstatmonthcl1Kanglenumlin',
                   #"fpath":'/fio/srvx7/leo/fastscratch/sabiner/datidl/leoformatstatmonthsinglelin',
                    "color":"orange",
                    "index":[0],
                    "dindex":numpy.zeros([len(stnames),1,rfpar['brmax']],numpy.int32),
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['brmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),
                   }))

        elif ds  in ['suny','sunyhom','sunycorr','rharm']:
                nam={'suny':'SUNY raw','sunyhom':'SUNY hom','rharm':'RHARM','sunycorr':'SUNY corr'}
                col={'suny':'burgundy','sunyhom':'avocado','rharm':'rose','sunycorr':'rose'}
                dspath={'suny':'SUNY','sunyhom':'SUNY','rharm':'RHARM','sunycorr':'SUNY'}
                dsdict.append({"shortname":ds,
                    "name":nam[ds],
                    "file":pathtrunk+"SUNY/",
                    "dfile":[pathtrunk+"SUNY/"],
                    "relpath":["./"],
                    "var":suff+"bias",
                    "dvar":suff+"bias",
                    "suff":"",    
                    "dsuff":[""],    
                    "color":col[ds],
                    "startdate":rfpar["startdate"],
                    "ens":[0],
                    "index":[0],
                    "dindex":[0],
                    "maxindex":numpy.zeros(initial_size,numpy.int32),
                    "minindex":numpy.zeros(initial_size,numpy.int32),
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),
                    'msudata':[initial_size,1,rfpar['parmax'],4,rfpar['mmax']],
                    'msufile':ds,
                    "msuvar":"montemp",'type':'ROB'
                    })
        elif ds in ['rcorr','bgdeprcorr','obsrcorr']:
                adpref=''
                if ds!='rcorr':
                    adpref=ds.split('rcorr')[0]+'-'
                dsdict.append({"shortname":ds,
                    "name":adpref+'RAOBCORE-adjustment',
                    "file":"feedbackglobbincorrsave",
                    "dfile":["feedbackglobbincorrsave"],
                    "msufile":"",
                    "st5":True,
                    "var":"rasocorr",
                    "dvar":"rasocorr",
                    "ens":[0],
                    "suff":"",    
                    "dsuff":[""],    
                    "startdate":rfpar["startdate"],
                    "color":"cyan",
                    "index":numpy.zeros([initial_size,1,rfpar['brmax']],numpy.int32),
                    "dindex":numpy.zeros([len(stnames),1,rfpar['brmax']],numpy.int32),
                    "maxindex":numpy.zeros(initial_size,numpy.int32),
                    "minindex":numpy.zeros(initial_size,numpy.int32),
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['brmax']],
                    'msudata':[initial_size,1,rfpar['parmax'],len(msupindex),rfpar['mmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),
                    'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']]),'type':'ROB',
                   })
        elif ds in ['rcorrann']:
                adpref=''
                if ds!='rcorr':
                    adpref=ds.split('rcorr')[0]+'-'
                dsdict.append({"shortname":ds,
                    "name":adpref+'RASE-adjustment',
                    "file":"feedbackglobbincorrsave",
                    "dfile":["CEUAS_merged_v1"],
                    "msufile":"",
                    "relpath":[os.path.expandvars('$RSCRATCH/converted_v7')],    
                    "st5":True,
                    "var":"rasocorr",
                    "dvar":"advanced_homogenisation/RASE_bias_estimate",
                    "ens":[0],
                    "suff":"",    
                    "dsuff":[""],    
                    "startdate":rfpar["startdate"],
                    "color":"cyan",
                    "index":numpy.zeros([initial_size,1,rfpar['brmax']],numpy.int32),
                    "dindex":numpy.zeros([len(stnames),1,rfpar['brmax']],numpy.int32),
                    "maxindex":numpy.zeros(initial_size,numpy.int32),
                    "minindex":numpy.zeros(initial_size,numpy.int32),
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['brmax']],
                    'msudata':[initial_size,1,rfpar['parmax'],len(msupindex),rfpar['mmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),
                    'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']]),'type':'ROB',
                   })
        elif ds in ('rio','rit','riocorr','ritcorr'):
                adpref=''
                if 'rio' in ds:
                    col='grey'
                    adpref='RICH-obs'
                else:
                    col='lavender'
                    adpref='RICH-tau'
                dsdict.append(make_ddict(**{"shortname":ds,
                    "name":adpref+' adjustments',
                    "file":"feedbackglobbincorrsave_"+ds,
                    "msufile":"feedbackglobbincorrsave_"+ds,
                    "dfile":["feedbackglobbincorrsave_"+ds],
                    "var":"rasocorr",
                    "dvar":"rasocorr",
                    "msuvar":"montemp",
                    "ens":numpy.arange(32,dtype=numpy.int32),
                    "suff":"",    
                    "dsuff":[""],    
                    "startdate":rfpar["startdate"],
                    "color":col,
                    "index":numpy.zeros([initial_size,32,rfpar['brmax']],numpy.int32),
                    "dindex":numpy.zeros([len(stnames),1,rfpar['brmax']],numpy.int32),
                    'data':[initial_size,32,rfpar['parmax'],rfpar['pmax'],rfpar['brmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['brmax']],
                    'msudata':[initial_size,32,rfpar['parmax'],4,rfpar['mmax']],
                    'zslopes':numpy.empty([32,18,len(jpindex)]),
                    'beltslopes':numpy.empty([32,6,len(jpindex)]),
                    'beltanomalies':numpy.empty([32,6,len(jpindex),rfpar['mmax']]),'type':'ROB'
                   }))
        elif ('rio' in ds or 'rit' in ds) and 'mean' not in ds:
            adpref=''
            if 'bgdep' in ds:
                adpref='bgdep - '
            if 'rio' in ds:
                col='grey'
                rv='RICH-obs'
            else:
                col='lavender'
                rv='RICH-tau'
            
            for i in range(32):
                if ds in ['rio{:0>2}'.format(i),'rit{:0>2}'.format(i),'riocorr{:0>2}'.format(i),'ritcorr{:0>2}'.format(i),'bgdeprio{:0>2}'.format(i),'bgdeprit{:0>2}'.format(i)]:
                    dsfile=ds[-5:]
                    if 'corr' in ds:
                        dsfile=''.join(ds.split('corr'))
                    dsdict.append(make_ddict(**{"shortname":ds,
                        "name":adpref+' '+rv+'{:0>2}-adjustment'.format(i),
                        "file":"feedbackglobbincorrsave_"+dsfile+"_",
                        "dfile":["feedbackglobbincorrsave_"+dsfile+"_"],
                        "var":"rasocorr",
                        "dvar":"rasocorr",
                        "msufile":"feedbackglobbincorrmon_"+dsfile,
                        "msuvar":"montemp",
                        'msudata':[initial_size,1,rfpar['parmax'],4,rfpar['mmax']],
                        "ens":[0],
                        "suff":"",    
                        "dsuff":[""],    
                        "startdate":rfpar["startdate"],
                        "color":col,
                        "index":numpy.zeros([initial_size,1,rfpar['brmax']],numpy.int32),
                        "dindex":numpy.zeros([len(stnames),1,rfpar['brmax']],numpy.int32),
                        'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                        'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['brmax']],
                        'zslopes':numpy.empty([1,18,len(jpindex)]),
                        'beltslopes':numpy.empty([1,6,len(jpindex)]),
                        'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']]),'type':'ROB'
                       }))

        elif ds in ('riomean','ritmean'):
                dsdict.append(make_ddict(**{"shortname":ds,
                    "name":'RICH ensemble mean adjustment',
                    "file":"feedbackglobbincorrsave_"+ds,
                    "dfile":["feedbackglobbincorrsave_rio"],
                    "var":"rasocorr",
                    "dvar":"rasocorr",
                    "msufile":"feedbackglobbincorrsave_"+ds+'_',
                    "msuvar":"rasocorrmon",
                    'msudata':[initial_size,1,rfpar['parmax'],4,rfpar['mmax']],
                    "ens":numpy.arange(32,dtype='int'),
                    "suff":"",    
                    "dsuff":[""],    
                    "startdate":rfpar["startdate"],
                    "color":"brown",
                    "index":numpy.zeros([initial_size,1,rfpar['brmax']],numpy.int32),
                    "dindex":numpy.zeros([len(stnames),32,rfpar['brmax']],numpy.int32),
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),32,rfpar['parmax'],len(pindex),rfpar['brmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),
                    'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']]),'type':'ROB'
                   }))


        elif ds in ['bgdep']:
                dsdict.append({"shortname":ds,
                    "name":'Background departures',
                    "file":"feedbackglobbgmon",
                    "dfile":[dailyfile],
                    "var":"montemp",
                    "dvar":plotproperties['fgdepname'], #"fg_dep",
                    "suff":"",   
                    "dsuff":dailysuff,
                    "color":"green",
                    "startdate":rfpar["startdate"],
                    "ens":[0],
                    "index":[0],
                    "dindex":[0],
                    "maxindex":numpy.zeros(initial_size,numpy.int32),
                    "minindex":numpy.zeros(initial_size,numpy.int32),
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),
                    'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']])
                    })
        elif ds in ['radbg']:
                dsdict.append({"shortname":ds,
                    "name":'Background',
                    "file":"feedbackglobbgmon",
                    "dfile":[dailyfile],
                    "var":"montemp",
                    "dvar":plotproperties['fgdepname'], #"fg_dep",
                    "suff":"",   
                    "dsuff":dailysuff,
                    "color":"green",
                    "startdate":rfpar["startdate"],
                    "ens":[0],
                    "index":[0],
                    "dindex":[0],
                    "maxindex":numpy.zeros(initial_size,numpy.int32),
                    "minindex":numpy.zeros(initial_size,numpy.int32),
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),
                    'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']])
                    })
        elif ds=='bgpresat':
                dsdict.append({"shortname":ds,
                    "name":'ERA-PreSAT background departures',
                    "file":"feedbackglobbgmon",
                    "dfile":[""],
                    "var":"montemp",
                    "dvar":"fg_dep",
                    "suff":"_t",    
                    "dsuff":["_t"],    
                    "relpath":["/fio/srvx7/leo/fastscratch/presat/"],    
                    "color":"green",
                    "startdate":rfpar["startdate"],
                    "ens":[0],
                    "index":[0],
                    "dindex":[0],
                    "maxindex":numpy.zeros(initial_size,numpy.int32),
                    "minindex":numpy.zeros(initial_size,numpy.int32),
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),
                    'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']])
                    })
        elif ds=='e5_fg_dep':
                dsdict.append({"shortname":ds,
                    "name":'Background departures',
                    "file":"feedbackglobbgmon",
                    "msufile":"",
                    'msudata':[initial_size,1,rfpar['parmax'],4,rfpar['mmax']],
                    "dfile":["ERA5_1_"],
                    "relpath":["/fio/srvx7/leo/fastscratch/ei6/"],
                    "var":"montemp",
                    "dvar":"fg_dep",
                    "suff":"_t",    
                    "dsuff":["_t"],    
                    "color":"green",
                    "startdate":rfpar["startdate"],
                    "ens":[0],
                    "index":[0],
                    "dindex":[0],
                    "maxindex":numpy.zeros(initial_size,numpy.int32),
                    "minindex":numpy.zeros(initial_size,numpy.int32),
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),
                    })
        elif ds=='e5_bias':
                dsdict.append({"shortname":ds,
                    "name":'Bias Estimate',
                    "file":"feedbackglobbgmon",
                    "dfile":["feedbackmerged"],
                    "relpath":["./"],
                    "var":"montemp",
                    "dvar":"bias_estimate",
                    "suff":"_t",    
                    "dsuff":[""],    
                    "color":"green",
                    "startdate":rfpar["startdate"],
                    "ens":[0],
                    "index":[0],
                    "dindex":[0],
                    "maxindex":numpy.zeros(initial_size,numpy.int32),
                    "minindex":numpy.zeros(initial_size,numpy.int32),
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),'type':'REA',
                    })
        elif ds=='e5_3116_bias':
                dsdict.append({"shortname":ds,
                    "name":'Bias Estimate',
                    "file":"feedbackglobbgmon",
                    "dfile":["ERA5_3116_"],
                    "relpath":["../../../ei6/"],
                    "var":"montemp",
                    "dvar":"bias",
                    "suff":"_t",    
                    "dsuff":["_t"],    
                    "color":"green",
                    "startdate":rfpar["startdate"],
                    "ens":[0],
                    "index":[0],
                    "dindex":[0],
                    "maxindex":numpy.zeros(initial_size,numpy.int32),
                    "minindex":numpy.zeros(initial_size,numpy.int32),
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),
                    })
        elif ds=='e5_an_dep':
                dsdict.append({"shortname":ds,
                    "name":'Background departures',
                    "file":"feedbackglobbgmon",
                    "dfile":["ERA5_e5_"],
                    "relpath":["../../../ei6/"],
                    "var":"montemp",
                    "dvar":"an_dep",
                    "suff":"_t",    
                    "dsuff":["_t"],    
                    "color":"gold",
                    "startdate":rfpar["startdate"],
                    "ens":[0],
                    "index":[0],
                    "dindex":[0],
                    "maxindex":numpy.zeros(initial_size,numpy.int32),
                    "minindex":numpy.zeros(initial_size,numpy.int32),
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),'type':'REA',
                    })
        elif ds=='oper_fg_dep':
                dsdict.append({"shortname":ds,
                    "name":'Background departures',
                    "file":"feedbackglobbgmon",
                    "dfile":["ERA5_oper_"],
                    "relpath":["../../../ei6/"],
                    "var":"montemp",
                    "dvar":"fg_dep",
                    "suff":"_t",    
                    "dsuff":["_t"],    
                    "color":"green",
                    "startdate":rfpar["startdate"],
                    "ens":[0],
                    "index":[0],
                    "dindex":[0],
                    "maxindex":numpy.zeros(initial_size,numpy.int32),
                    "minindex":numpy.zeros(initial_size,numpy.int32),
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),'type':'REA',
                    })
        elif ds=='oper_bias':
                dsdict.append({"shortname":ds,
                    "name":'Bias Estimate',
                    "file":"feedbackglobbgmon",
                    "dfile":["ERA5_oper_"],
                    "relpath":["../../../ei6/"],
                    "var":"montemp",
                    "dvar":"bias",
                    "suff":"_t",    
                    "dsuff":["_t"],    
                    "color":"green",
                    "startdate":rfpar["startdate"],
                    "ens":[0],
                    "index":[0],
                    "dindex":[0],
                    "maxindex":numpy.zeros(initial_size,numpy.int32),
                    "minindex":numpy.zeros(initial_size,numpy.int32),
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),
                    })
        elif ds=='oper_an_dep':
                dsdict.append({"shortname":ds,
                    "name":'Background departures',
                    "file":"feedbackglobbgmon",
                    "dfile":["ERA5_oper_"],
                    "relpath":["../../../ei6/"],
                    "var":"montemp",
                    "dvar":"an_dep",
                    "suff":"_t",    
                    "dsuff":["_t"],    
                    "color":"gold",
                    "startdate":rfpar["startdate"],
                    "ens":[0],
                    "index":[0],
                    "dindex":[0],
                    "maxindex":numpy.zeros(initial_size,numpy.int32),
                    "minindex":numpy.zeros(initial_size,numpy.int32),
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),'type':'REA',
                    })
        elif ds=='andep':
                dsdict.append({"shortname":ds,
                    "name":'Background departures',
                    "file":"feedbackglobanmon",
                    "dfile":["feedbackmerged"],
                    "var":"montemp",
                    "dvar":"an_dep",
                    "suff":"",    
                    "color":"gold",
                    "startdate":rfpar["startdate"],
                    "ens":[0],
                    "index":[0],
                    "dindex":[0],
                    "maxindex":numpy.zeros(initial_size,numpy.int32),
                    "minindex":numpy.zeros(initial_size,numpy.int32),
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),'type':'REA'
                    })
        elif ds  in  ['era5v2','era5v2rich','era5v2corr','era5v2richcorr']:
                color="yellow_green"
                if 'rich' in ds:
                    suff='rich'
                    color='yellow'
                dsdict.append({"shortname":ds,
                    "name":'ERA5 '+suff.upper()+' adjustments',
                    "file": pathtrunk+"rise/1.0/ERA5_v2//ERA5bc_RAOBCORE_v1.5_",
                    "dfile":[pathtrunk+"rise/1.0/ERA5_v2/ERA5bc_RAOBCORE_v1.5_"],
                    "var":suff+"bias",
                    "dvar":suff+"bias",
                    "suff":"",    
                    "dsuff":[""],    
                    "color":color,
                    "startdate":rfpar["startdate"],
                    "ens":[0],
                    "index":[0],
                    "dindex":[0],
                    "maxindex":numpy.zeros(initial_size,numpy.int32),
                    "minindex":numpy.zeros(initial_size,numpy.int32),
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),
                    'msudata':[initial_size,1,rfpar['parmax'],4,rfpar['mmax']],
                    'msufile':'',
                    "msuvar":"rasocorrmon",'type':'ROB'
                    })
        elif ds  in  ['era5v2429']:
                color="orange"
                dsdict.append({"shortname":ds,
                    "name":'ERA5 '+suff.upper()+' adjustments',
                    "file": pathtrunk+"rise/1.0/ERA5_v2429/ERA5_2429_",
                    "dfile":[pathtrunk+"rise/1.0/ERA5_v2429/ERA5_2429_"],
                    "var":suff+"bias",
                    "dvar":suff+"bias",
                    "suff":"_t",    
                    "dsuff":["_t"],    
                    "color":color,
                    "startdate":rfpar["startdate"],
                    "ens":[0],
                    "index":[0],
                    "dindex":[0],
                    "maxindex":numpy.zeros(initial_size,numpy.int32),
                    "minindex":numpy.zeros(initial_size,numpy.int32),
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),
                    'msudata':[initial_size,1,rfpar['parmax'],4,rfpar['mmax']],
                    'msufile':'',
                    "msuvar":"rasocorrmon",'type':'ROB'
                    })
        elif ds=='era5v3':
                dsdict.append({"shortname":ds,
                    "name":'ERA5 adjustments',
                    #"file":pathtrunk+"rise/1.0/ERA5_v3/ERA5bc_RAOBCORE_v1.6_",
                    #"dfile":[pathtrunk+"rise/1.0/ERA5_v3/ERA5bc_RAOBCORE_v1.6_"],
                    "file":pathtrunk+"rise/1.0/exp01/ERA5bc_RAOBCORE_v1.6_",
                    "dfile":["ERA5bc_RAOBCORE_v1.6_"],
                    "var":"richbias",
                    "dvar":"bias",
                    "suff":"",    
                    "dsuff":[""],    
                    "color":"yellow_green",
                    "startdate":rfpar["startdate"],
                    "ens":[0],
                    "index":[0],
                    "dindex":[0],
                    "maxindex":numpy.zeros(initial_size,numpy.int32),
                    "minindex":numpy.zeros(initial_size,numpy.int32),
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),
                    'msudata':[initial_size,1,rfpar['parmax'],4,rfpar['mmax']],
                    'msufile':'',
                    "msuvar":"rasocorrmon",'type':'ROB'
                    })
        elif ds  in ['era5v4','era5v4rich']:
                if 'rich' in ds:
                    suff='rich'
                dsdict.append({"shortname":ds,
                    "name":'ERA5 adjustments',
                    "file":pathtrunk+"rise/1.0/ERA5_v4/ERA5bc_RAOBCORE_v1.6_",
                    "dfile":[pathtrunk+"rise/1.0/ERA5_v4/ERA5bc_RAOBCORE_v1.6_"],
                    "relpath":["./"],
                    "var":suff+"bias",
                    "dvar":suff+"bias",
                    "suff":"",    
                    "dsuff":[""],    
                    "color":"yellow_green",
                    "startdate":rfpar["startdate"],
                    "ens":[0],
                    "index":[0],
                    "dindex":[0],
                    "maxindex":numpy.zeros(initial_size,numpy.int32),
                    "minindex":numpy.zeros(initial_size,numpy.int32),
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),
                    'msudata':[initial_size,1,rfpar['parmax'],4,rfpar['mmax']],
                    'msufile':'',
                    "msuvar":"rasocorrmon",'type':'ROB'
                    })
        elif ds in ['era5v5','era5v5rich']:
                if 'rich' in ds:
                    suff='rich'
                dsdict.append({"shortname":ds,
                    "name":'ERA5 RICH adjustments',
                    "file":pathtrunk+"rise/1.0/ERA5_v5/ERA5bc_RAOBCORE_v1.6_",
                    "dfile":[pathtrunk+"rise/1.0/ERA5_v5/ERA5bc_RAOBCORE_v1.6_"],
                    "relpath":["./"],
                    "var":suff+"bias",
                    "dvar":suff+"bias",
                    "suff":"",    
                    "dsuff":[""],    
                    "color":"yellow_green",
                    "startdate":rfpar["startdate"],
                    "ens":[0],
                    "index":[0],
                    "dindex":[0],
                    "maxindex":numpy.zeros(initial_size,numpy.int32),
                    "minindex":numpy.zeros(initial_size,numpy.int32),
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),
                    'msudata':[initial_size,1,rfpar['parmax'],4,rfpar['mmax']],
                    'msufile':'',
                    "msuvar":"rasocorrmon",'type':'ROB'
                    })
        elif ds in ['era5v7','era5v7rich']:
                color='yellow_green'
                if 'rich' in ds:
                    suff='rich'
                    color='yellow'
                dsdict.append({"shortname":ds,
                    "name":'ERA5 RICH adjustments',
                    "file":pathtrunk+"rise/1.0/ERA5_v7/ERA5bc_RAOBCORE_v1.7_",
                    "dfile":[pathtrunk+"rise/1.0/ERA5_v7/ERA5bc_RAOBCORE_v1.7_"],
                    "relpath":["./"],
                    "var":suff+"bias",
                    "dvar":suff+"bias",
                    "suff":"",    
                    "dsuff":[""],    
                    "color":color,
                    "startdate":rfpar["startdate"],
                    "ens":[0],
                    "index":[0],
                    "dindex":[0],
                    "maxindex":numpy.zeros(initial_size,numpy.int32),
                    "minindex":numpy.zeros(initial_size,numpy.int32),
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),
                    'msudata':[initial_size,1,rfpar['parmax'],4,rfpar['mmax']],
                    'msufile':'',
                    "msuvar":"rasocorrmon",'type':'ROB'
                    })
        elif ds=='eraibc':
                dsdict.append({"shortname":ds,
                    "name":'ERAI adjustments',
                    "relpath":[pathtrunk+"ei6/"],
                    "file":"" , #"../exp00/biasmon",
                    "dfile":[""],
                    "var":"bias",
                    "dvar":"bias",
                    "dsuff":["_t"],    
                    "suff":"_t",    
                    "color":"yellow",
                    "startdate":rfpar["startdate"],
                    "ens":[0],
                    "index":[0],
                    "dindex":[0],
                    "maxindex":numpy.zeros(initial_size,numpy.int32),
                    "minindex":numpy.zeros(initial_size,numpy.int32),
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),
                    'msudata':[initial_size,1,rfpar['parmax'],4,rfpar['mmax']],
                    'msufile':'',
                    "msuvar":"rasocorrmon",'type':'ROB'
                    })
        elif ds=='oper_bc':
                dsdict.append({"shortname":ds,
                    "name":'operational adjustments',
                    "file":"feedbackglobbgmon",
                    "dfile":["ERA5_oper_"],
                    "relpath":["../../../ei6/"],
                    "dvar":"bias",
                    "dsuff":["_t"],    
                    "var":"bias",
                    "suff":"_t",    
                    "color":"yellow_green",
                    "startdate":rfpar["startdate"],
                    "ens":[0],
                    "index":[0],
                    "dindex":[0],
                    "maxindex":numpy.zeros(initial_size,numpy.int32),
                    "minindex":numpy.zeros(initial_size,numpy.int32),
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),
                    'msudata':[initial_size,1,rfpar['parmax'],4,rfpar['mmax']],
                    'msufile':'',
                    "msuvar":"rasocorrmon",'type':'ROB'
                    })
 #               andeplist=['jra55_andep','jra55_fgdep','e20c_andep','n20c_andep','erapresat_andep','bgdep',
 #                          'eijra_fgdep','eice20c_andep','jrace20c_andep','fg_dep','an_dep',
 #                          'comp-bg','comp2-bg','comp4-bg', 'era5_andep','era5_fgdep','e5_fg_dep','e5_an_dep','e5_bias','e5_3116_bias','oper_bias','oper_fg_dep',
 #                'ce20c0_andep','eice20c0_andep','eice20c1_andep','eice20c2_andep','eice20c3_andep','eice20c4_andep','eice20c5_andep',
 #                'eice20c6_andep','eice20c7_andep','eice20c8_andep','eice20c9_andep','ce20cens_spread','ce20cens_mean','erai_fggpsdep','erai_fggpswetdep']
        elif ds in andeplist:
            col=['yellow_green','olive','blue_purple','olive','avocado','gold','brown','yellow_green','peach']+['green']*7+['peach']*10+['orange_yellow','turquoise','khaki']
            while len(col)<len(andeplist):
                col=col+['black']
            col=col[andeplist.index(ds)]
            dsdict.append({"shortname":ds,
                "name":ds[:-3]+' departures',
                "file":ds+"mon",
                "prefix":"./",
                "dfile":[dailyfile],
                "dsuff":dailysuff,
                "var":ds,
                "dvar":ds,
                "suff":"",    
                "color":col,
                'msudata':[initial_size,1,rfpar['parmax'],4,rfpar['mmax']],
                'msufile':ds+"mon",
                "msuvar":"montemp",
                "startdate":rfpar["startdate"],
                "ens":[0],
                "index":[0],
                "dindex":[0],
                "maxindex":numpy.zeros(initial_size,numpy.int32),
                "minindex":numpy.zeros(initial_size,numpy.int32),
                'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                'zslopes':numpy.empty([1,18,len(jpindex)]),
                'beltslopes':numpy.empty([1,6,len(jpindex)]),
                'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']]),'type':'REA'
                })
            if ds=='comp-bg':
                dsdict[-1]["dfile"]=['feedbackglobbincomp']
                dsdict[-1]["file"]='feedbackglobbincompmon'
                dsdict[-1]["var"]='rasocorrmon'
                dsdict[-1]["name"]='RAOBCORE composite-bg'
            if ds=='comp2-bg':
                dsdict[-1]["dfile"]=['feedbackglobbincomp_2_']
                dsdict[-1]["file"]='feedbackglobbincompmon_2_'
                dsdict[-1]["var"]='rasocorrmon'
                dsdict[-1]["dvar"]='comp-bg'
                dsdict[-1]["name"]='RICH-obs composite-bg'
            if ds=='comp4-bg':
                dsdict[-1]["dfile"]=['feedbackglobbincomp_4_']
                dsdict[-1]["file"]='feedbackglobbincompmon_4_'
                dsdict[-1]["var"]='rasocorrmon'
                dsdict[-1]["dvar"]='comp-bg'
                dsdict[-1]["name"]='RICH-obs composite-bg'
            if 'gps' in ds:
                dsdict[-1]["type"]='SAT'
               
                
        elif ds in andepstdlist:
            col=['yellow_green','olive','blue_purple','olive','avocado','gold']+['peach']*10+['orange_yellow','turquoise','kahki']+['green']*30
            col=col[andepstdlist.index(ds)]
            dsdict.append({"shortname":ds,
                "name":ds[:-6]+' Analysis departure standard deviations',
                "file":ds[:-4]+"mon",
                "prefix":"./",
                "dfile":['feedbackmerged'],
                "dsuff":[''],
                "var":ds,
                "dvar":ds,
                "suff":"",    
                "color":col,
                'msudata':[initial_size,1,rfpar['parmax'],4,rfpar['mmax']],
                'msufile':'',
                "msuvar":"rasocorrmon",
                "startdate":rfpar["startdate"],
                "ens":[0],
                "index":[0],
                "dindex":[0],
                "maxindex":numpy.zeros(initial_size,numpy.int32),
                "minindex":numpy.zeros(initial_size,numpy.int32),
                'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                'zslopes':numpy.empty([1,18,len(jpindex)]),
                'beltslopes':numpy.empty([1,6,len(jpindex)]),
                'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']]),'type':'REA'
                })
        elif ds in andeprmslist:
            col=['yellow_green','olive','blue_purple','olive','avocado','gold']+['peach']*10+['orange_yellow','turquoise','kahki']+['green']*30
            col=col[andeprmslist.index(ds)]
            dsdict.append({"shortname":ds,
                "name":ds[:-6]+' Analysis departure rms',
                "file":ds[:-4]+"mon",
                "prefix":"./",
                "dfile":['feedbackmerged'],
                "dsuff":[''],
                "var":'_std'.join(ds.split('_rms')),
                "dvar":ds,
                "suff":"",    
                "color":col,
                'msudata':[initial_size,1,rfpar['parmax'],4,rfpar['mmax']],
                'msufile':'',
                "msuvar":"rasocorrmon",
                "startdate":rfpar["startdate"],
                "ens":[0],
                "index":[0],
                "dindex":[0],
                "maxindex":numpy.zeros(initial_size,numpy.int32),
                "minindex":numpy.zeros(initial_size,numpy.int32),
                'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                'zslopes':numpy.empty([1,18,len(jpindex)]),
                'beltslopes':numpy.empty([1,6,len(jpindex)]),
                'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']]),'type':'REA'
                })
        elif ds in ['jra55_fgdep','eijra_fgdep']:
            dsdict.append({"shortname":ds,
                "name":ds[:-6]+' Background departures',
                "file":ds+"mon",
                "prefix":"./",
                "dfile":['feedbackmerged'],
                "dsuff":[''],
                "var":ds,
                "dvar":ds,
                "suff":"",    
                "color":"green",
                'msudata':[initial_size,1,rfpar['parmax'],4,rfpar['mmax']],
                'msufile':'',
                "msuvar":"rasocorrmon",
                "startdate":rfpar["startdate"],
                "ens":[0],
                "index":[0],
                "dindex":[0],
                "maxindex":numpy.zeros(initial_size,numpy.int32),
                "minindex":numpy.zeros(initial_size,numpy.int32),
                'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                'zslopes':numpy.empty([1,18,len(jpindex)]),
                'beltslopes':numpy.empty([1,6,len(jpindex)]),
                'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']]),'type':'REA'
                })

        elif ds in ['era5v5-rcorr']:
            dsdict.append({"shortname":ds,
                "name":ds[:-6],
                "file":"",
                "prefix":"./",
                "dfile":[''],
                "dsuff":[''],
                "var":ds,
                "dvar":ds,
                "suff":"",    
                "color":"green",
                'msudata':[initial_size,1,rfpar['parmax'],4,rfpar['mmax']],
                'msufile':'',
                "msuvar":"rasocorrmon",
                "startdate":rfpar["startdate"],
                "ens":[0],
                "index":[0],
                "dindex":[0],
                "maxindex":numpy.zeros(initial_size,numpy.int32),
                "minindex":numpy.zeros(initial_size,numpy.int32),
                'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                'zslopes':numpy.empty([1,18,len(jpindex)]),
                'beltslopes':numpy.empty([1,6,len(jpindex)]),
                'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']]),'type':'REA'
                })

        elif ds in ['bgdep-era5v5','bgdep-era5v4','bgdep-era5v3','bgdep-erai','bgdep-era5v5rich','bgdep-era5v4rich','bgdep-era5v3rich']:
            dsdict.append({"shortname":ds,
                "name":ds,
                "file":"",
                "prefix":"./",
                "dfile":[''],
                "dsuff":[''],
                "var":ds,
                "dvar":ds,
                "suff":"",    
                "color":"green",
                'msudata':[initial_size,1,rfpar['parmax'],4,rfpar['mmax']],
                'msufile':'',
                "msuvar":"rasocorrmon",
                "startdate":rfpar["startdate"],
                "ens":[0],
                "index":[0],
                "dindex":[0],
                "maxindex":numpy.zeros(initial_size,numpy.int32),
                "minindex":numpy.zeros(initial_size,numpy.int32),
                'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                'zslopes':numpy.empty([1,18,len(jpindex)]),
                'beltslopes':numpy.empty([1,6,len(jpindex)]),
                'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']]),'type':'REA'
                })

        elif ds in ['ce20cens_andep']:
                dsdict.append({"shortname":ds,
                "prefix":"./",
                'msufile':'',
                "msuvar":"rasocorrmon",
                "startdate":rfpar["startdate"],
                "maxindex":numpy.zeros(initial_size,numpy.int32),
                "minindex":numpy.zeros(initial_size,numpy.int32),
                'beltanomalies':numpy.empty([10,6,len(jpindex),rfpar['mmax']]),'type':'REA'})
                dsdict[-1]["ens"]=numpy.arange(10,dtype=numpy.int32)
                dsdict[-1]["file"]="ce20c"
                dsdict[-1]["dvar"]='ce20c'
                dsdict[-1]["var"]='ce20c'
                dsdict[-1]["suff"]='_andepmon'
                dsdict[-1]["dsuff"]='_andep'
                dsdict[-1]["shortname"]=ds
                dsdict[-1]["name"]='CE20C Analysis departures'
                dsdict[-1]["color"]='ochre'
                dsdict[-1]["index"]=[0]
                dsdict[-1]["dindex"]=[0]
                dsdict[-1]['data']=[initial_size,10,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']]
                dsdict[-1]['ddata']=[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']]
                dsdict[-1]['msudata']=[initial_size,10,rfpar['parmax'],4,rfpar['mmax']]
                dsdict[-1]['zslopes']=numpy.empty([10,18,len(jpindex)])
                dsdict[-1]['beltslopes']=numpy.empty([10,6,len(jpindex)])
                
        elif ds in ['ce20c_andep']:
                dsdict.append({"shortname":ds,
                "prefix":"./",
                'msufile':'',
                "msuvar":"rasocorrmon",
                "startdate":rfpar["startdate"],
                "maxindex":numpy.zeros(initial_size,numpy.int32),
                "minindex":numpy.zeros(initial_size,numpy.int32),
                'beltanomalies':numpy.empty([10,6,len(jpindex),rfpar['mmax']]),'type':'REA'})
                dsdict[-1]["ens"]=numpy.arange(10,dtype=numpy.int32)
                dsdict[-1]["file"]="ce20c"
                dsdict[-1]["dvar"]='ce20c'
                dsdict[-1]["var"]='ce20c'
                dsdict[-1]["suff"]='_andepmon'
                dsdict[-1]["dsuff"]='_andep'
                dsdict[-1]["shortname"]=ds
                dsdict[-1]["name"]='CE20C Analysis departures'
                dsdict[-1]["color"]='ochre'
                dsdict[-1]["index"]=[0]
                dsdict[-1]["dindex"]=[0]
                dsdict[-1]['data']=[initial_size,10,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']]
                dsdict[-1]['ddata']=[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']]
                dsdict[-1]['msudata']=[initial_size,10,rfpar['parmax'],4,rfpar['mmax']]
                dsdict[-1]['zslopes']=numpy.empty([10,18,len(jpindex)])
                dsdict[-1]['beltslopes']=numpy.empty([10,6,len(jpindex)])
                
        elif ds in ['e20c_pandep']:
                dsdict.append({"shortname":ds,
                'msufile':'',
                "msuvar":"rasocorrmon",
                "startdate":rfpar["startdate"],
                "maxindex":numpy.zeros(initial_size,numpy.int32),
                "minindex":numpy.zeros(initial_size,numpy.int32),
                'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']])})
                dsdict[-1]["ens"]=[0]
                dsdict[-1]["file"]=""
                dsdict[-1]["dfile"]=[""]
                dsdict[-1]["dvar"]='e20c_0'
                dsdict[-1]["var"]='e20c_0'
                dsdict[-1]["suff"]='_andepmon'
                dsdict[-1]["dsuff"]=['_20c']
                dsdict[-1]["relpath"]="/fio/srvx7/leo/fastscratch/presat/",    
                dsdict[-1]["shortname"]=ds
                dsdict[-1]["name"]='E20C Analysis departures'
                dsdict[-1]["color"]='ochre'
                dsdict[-1]["index"]=[0]
                dsdict[-1]["dindex"]=[0]
                dsdict[-1]['data']=[initial_size,10,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']]
                dsdict[-1]['ddata']=[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']]
                dsdict[-1]['msudata']=[initial_size,10,rfpar['parmax'],4,rfpar['mmax']]
                dsdict[-1]['zslopes']=numpy.empty([10,18,len(jpindex)])
                dsdict[-1]['beltslopes']=numpy.empty([10,6,len(jpindex)])
                
        elif ds in ['bg']:
            dsdict.append(make_ddict(**{"shortname":'bg',
                "name":'Background',
                "file":"feedbackglobbgmon",
                "msufile":"feedbackglobbgmon",
                "dfile":["",'feedbackmerged'],
                "dsuff":["_t",''],
                "var":"montemp",
                "dvar":"fg_dep",
                "msuvar":"montemp",
                "suff":"",    
                "color":"reddish_purple",
                "startdate":rfpar["startdate"],
                "ens":[0],
                "index":[0],
                "dindex":[0],
                "maxindex":numpy.zeros(initial_size,numpy.int32),
                "minindex":numpy.zeros(initial_size,numpy.int32),
                'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                'msudata':[initial_size,1,rfpar['parmax'],4,rfpar['mmax']],
                'zslopes':numpy.empty([1,18,len(jpindex)]),
                'beltslopes':numpy.empty([1,6,len(jpindex)]),
                'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']]),'type':'REA'
                }))
        elif ds in ['obs','rad','radc']:
            if plotproperties['exps'][0] in ['3188','1759','1761','ai_bfr','igra','erai']:
                dsdict.append({"shortname":ds,
                    "name":'Observations',
                    "file":"feedbackglobbgmon",
                    "dfile":[dailyfile],
                    "var":"temperatures",
                    "dvar":dailyobsvar[0],
                    "suff":"",
                    "dsuff":dailysuff,
                    "color":"blue",
                    "startdate":rfpar["startdate"],
                    "ens":[0],
                    "index":[0],
                    "dindex":[0],
                    "maxindex":numpy.zeros(initial_size,numpy.int32),
                    "minindex":numpy.zeros(initial_size,numpy.int32),
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),'type':'ROB',
                    })
            else:
                dsdict.append({"shortname":ds,
                    "name":'Observations',
                    "file":"feedbackglobbgmon",
                    "dfile":["feedbackmerged"],
                    "var":"temperatures",
                    "dvar":dailyobsvar[0],
                    "suff":"",
                    "dsuff":dailysuff,
                    "color":"blue",
                    "startdate":rfpar["startdate"],
                    "ens":[0],
                    "index":[0],
                    "dindex":[0],
                    "maxindex":numpy.zeros(initial_size,numpy.int32),
                    "minindex":numpy.zeros(initial_size,numpy.int32),
                    'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
                    'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
                    'zslopes':numpy.empty([1,18,len(jpindex)]),
                    'beltslopes':numpy.empty([1,6,len(jpindex)]),'type':'ROB',
                    })
            l+=1
        else:
            print(('Key '+ds+' not found'))
                
    return dsdict



#e20c_andep["color"]='orange'

#n20c_andep["color"]='brown'

#erapresat_andep=copy.deepcopy(ce20c0_andep)
#erapresat_andep["file"]="erapresat_andepmon"
#erapresat_andep["dvar"]='erapresat_andep'
#erapresat_andep["var"]='erapresat_andep'
#erapresat_andep["shortname"]='erapresat_andep'
#erapresat_andep["name"]='ERA-preSAT Analysis departures'

#e20cminusce20c=copy.deepcopy(ce20c0_andep)
#e20cminusce20c["dvar"]='e20c_andep'
#e20cminusce20c["shortname"]='e20c-ce20c0'
#e20cminusce20c["name"]='E20C minus CE20C0'



#eijra_fgdeprcorr=rcorr.copy()
#eijra_fgdeprcorr["shortname"]='eijradeprcorr'
#eijra_fgdeprcorr["name"]='Adjusted ei+jra background departures'

#eijra_fgdepriocorr=rio.copy()
#eijra_fgdepriocorr["shortname"]='eijradepriocorr'
#eijra_fgdepriocorr["name"]='RICH-adjusted ei+jra background departures'


#fg_dep={"shortname":'fg_dep',
    #"name":'Background',
    #"file":"fg_depmon",
    #"dfile":["feedbackmerged"],
    #"var":"fg_dep",
    #"dvar":"fg_dep",
    #"suff":"",    
    #"color":"green",
    #"startdate":rfpar["startdate"],
    #"ens":[0],
    #"index":[0],
    #"dindex":[0],
    #"maxindex":numpy.zeros(initial_size,numpy.int32),
    #"minindex":numpy.zeros(initial_size,numpy.int32),
    #"msufile":"",
    #"msuvar":"",
    #'msudata':[initial_size,1,rfpar['parmax'],4,rfpar['mmax']],
    #'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
    #'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
    #'msudata':[initial_size,1,rfpar['parmax'],4,rfpar['mmax']],
    #'zslopes':numpy.empty([1,18,len(jpindex)]),
    #'beltslopes':numpy.empty([1,6,len(jpindex)]),
    #'beltanomalies':numpy.empty([1,6,len(jpindex),rfpar['mmax']])
    #}
#jra55={"shortname":'jra55',
    #"name":'JRA55 analysis',
    #"file":"feedbackglobbgmon",
    #"dfile":["feedbackmerged"],
    #"var":"jra55_andep",
    #"dvar":"jra55_andep",
    #"suff":"",    
    #"color":"green",
    #"startdate":rfpar["startdate"],
    #"ens":[0],
    #"index":[0],
    #"dindex":[0],
    #"maxindex":numpy.zeros(initial_size,numpy.int32),
    #"minindex":numpy.zeros(initial_size,numpy.int32),
    #'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
    #'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
    #'zslopes':numpy.empty([1,18,len(pindex)]),
    #'beltslopes':numpy.empty([1,6,len(jpindex)]),
    #}
#obs={"shortname":'obs',
    #"name":'Observations',
    #"file":"feedbackglobbgmon",
    #"dfile":["feedbackmerged"],
    #"var":"jra55_andep",
    #"dvar":"temperatures",
    #"suff":"",    
    #"color":"green",
    #"startdate":rfpar["startdate"],
    #"ens":[0],
    #"index":[0],
    #"dindex":[0],
    #"maxindex":numpy.zeros(initial_size,numpy.int32),
    #"minindex":numpy.zeros(initial_size,numpy.int32),
    #'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
    #'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
    #'zslopes':numpy.empty([1,18,len(pindex)]),
    #'beltslopes':numpy.empty([1,6,len(jpindex)]),
    #}

#bgdiff=copy.deepcopy(bgpresat)
#bgdiff['shortname']='bgdiff'
#bgdiff['name']='Background departures'

#anpresat={"shortname":'an',
    #"name":'Analysis departures',
    #"file":"feedbackglobanmon",
    #"dfile":["","feedbackmerged"],
    #"dsuff":["_t",""],
    #"var":"montemp",
    #"dvar":"an_dep",
    #"suff":"",    
    #"color":"green",
    #"startdate":rfpar["startdate"],
    #"ens":[0],
    #"index":[0],
    #"dindex":[0],
    #"maxindex":numpy.zeros(initial_size,numpy.int32),
    #"minindex":numpy.zeros(initial_size,numpy.int32),
    #'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
    #'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
    #'zslopes':numpy.empty([1,18,len(jpindex)]),
    #'beltslopes':numpy.empty([1,6,len(jpindex)]),
    #}
#presatminuse20c=copy.deepcopy(bgpresat)
#presatminuse20c["shortname"]='presat-e20c'
#presatminuse20c["name"]='bgpresat-e20c'
#anpresat20c={"shortname":'an20c',
    #"name":'E20C Analysis departures',
    #"file":"feedbackglobanmon",
    #"prefix":"../../presat",
    #"dfile":["","feedbackmergede20c"],
    #"dsuff":["_20c",""],
    #"var":"montemp",
    #"dvar":"e20c_0",
    #"suff":"",    
    #"color":"green",
    #"startdate":rfpar["startdate"],
    #"ens":[0],
    #"index":[0],
    #"dindex":[0],
    #"maxindex":numpy.zeros(initial_size,numpy.int32),
    #"minindex":numpy.zeros(initial_size,numpy.int32),
    #'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
    #'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
    #'zslopes':numpy.empty([1,18,len(jpindex)]),
    #'beltslopes':numpy.empty([1,6,len(jpindex)]),
    #}
#anpresatn20cr={"shortname":'an20cr',
    #"name":'NOAA 20CR Analysis departures',
    #"file":"feedbackglobanmon",
    #"prefix":"../../presat",
    #"dfile":["","feedbackmergedn20cr"],
    #"dsuff":["_n20cr",""],
    #"var":"montemp",
    #"dvar":"e20c_0",
    #"suff":"",    
    #"color":"green",
    #"startdate":rfpar["startdate"],
    #"ens":[0],
    #"index":[0],
    #"dindex":[0],
    #"maxindex":numpy.zeros(initial_size,numpy.int32),
    #"minindex":numpy.zeros(initial_size,numpy.int32),
    #'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
    #'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
    #'zslopes':numpy.empty([1,18,len(jpindex)]),
    #'beltslopes':numpy.empty([1,6,len(jpindex)]),
    #}
#aninc={"shortname":'aninc',
    #"name":'Analysis increments',
    #"file":"feedbackglobanmon",
    #"dfile":["","feedbackmerged"],
    #"dsuff":["_t",""],
    #"var":"montemp",
    #"dvar":"an_dep",
    #"suff":"",    
    #"color":"green",
    #"startdate":rfpar["startdate"],
    #"ens":[0],
    #"index":[0],
    #"dindex":[0],
    #"maxindex":numpy.zeros(initial_size,numpy.int32),
    #"minindex":numpy.zeros(initial_size,numpy.int32),
    #'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
    #'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
    #'zslopes':numpy.empty([1,18,len(jpindex)]),
    #'beltslopes':numpy.empty([1,6,len(jpindex)]),
    #}
#n20crminuse20c=copy.deepcopy(anpresatn20cr)
#n20crminuse20c["shortname"]='20cr-e20c'
#n20crminuse20c["name"]='20cr-e20c'
#anince20c={"shortname":'anince20c',
    #"name":'Analysis increments',
    #"file":"feedbackglobanmon",
    #"dfile":["","feedbackmergede20c"],
    #"dsuff":["_20c",""],
    #"var":"montemp",
    #"dvar":"e20c_0",
    #"suff":"",    
    #"color":"green",
    #"startdate":rfpar["startdate"],
    #"ens":[0],
    #"index":[0],
    #"dindex":[0],
    #"maxindex":numpy.zeros(initial_size,numpy.int32),
    #"minindex":numpy.zeros(initial_size,numpy.int32),
    #'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
    #'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
    #'zslopes':numpy.empty([1,18,len(jpindex)]),
    #'beltslopes':numpy.empty([1,6,len(jpindex)]),
    #}

#bgrcorr=rcorr.copy()
#bgrcorr["shortname"]='bgdeprcorr'
#bgrcorr["name"]='Adjusted background departures'

#bgcorrection={"shortname":'bgcorrection',
    #"name":'Background adjustment',
    #"file":"feedbackglobbincorrmon",
    #"dfile":["feedbackbgcorr"],
    #"st5":True,
    #"var":"montemp",
    #"dvar":"bg correctio",
    #"ens":[0],
    #"suff":[""],    
    #"dsuff":[""],    
    #"startdate":rfpar["startdate"],
    #"color":"cyan",
    #"index":[0],
    #"dindex":[0],
    #'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
    #'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
    #'zslopes':numpy.empty([1,18,len(jpindex)]),
    #'beltslopes':numpy.empty([1,6,len(jpindex)]),
   #}


#bgriocorr=rio.copy()
#bgriocorr["shortname"]='bgdepriocorr'
#bgriocorr["name"]='RICH-adjusted background departures'


#bgdeppresat=make_ddict(**{"shortname":'bgdep',
    #"name":'Background departures',
    #"file":"feedbackglobbgmon",
    #"dfile":["",'feedbackmerged'],
    #"dsuff":["_t",''],
    #"var":"montemp",
    #"dvar":"fg_dep",
    #"suff":"",    
    #"color":"green",
    #"startdate":rfpar["startdate"],
    #"ens":[0],
    #"index":[0],
    #"dindex":[0],
    #"maxindex":numpy.zeros(initial_size,numpy.int32),
    #"minindex":numpy.zeros(initial_size,numpy.int32),
    #'data':[initial_size,1,rfpar['parmax'],rfpar['pmax'],rfpar['mmax']],
    #'ddata':[len(stnames),1,rfpar['parmax'],len(pindex),rfpar['nmax']],
    #'zslopes':numpy.empty([1,18,len(jpindex)]),
    #'beltslopes':numpy.empty([1,6,len(jpindex)]),
    #})
