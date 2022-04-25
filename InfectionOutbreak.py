from numpy import *
import numpy as np
import matplotlib.pyplot as plt
from tracerecorder import TraceRecorder
import matplotlib.gridspec as gridspec
%matplotlib inline
def step1(itrial,t):
    global nusers,nmachines
    global mv,mt,uv,ut,ui
    global tr
    
    # create views on next day's array slices
    MV = mv[itrial,t+1,:,:]
    MT = mt[itrial,t+1,:,:]
    UV = uv[itrial,t+1,:,:]
    UT = ut[itrial,t+1,:,:]
    UI = ui[itrial,t+1,:,:]

    # copy current state to day t+1
    MV[:,:] = mv[itrial,t,:,:]
    MT[:,:] = mt[itrial,t,:,:]
    UV[:,:] = uv[itrial,t,:,:]
    UT[:,:] = ut[itrial,t,:,:]
    UI[:,:] = ui[itrial,t,:,:]
    
    # then modify state on day t+1
    MT[MV] += 1 # age the germs on the infected machines by 1
    UT[UV] += 1 # advance the infections in the users by 1
     
    MV[ MT >= viral_lifetime ]     = False # let the old viruses on the machines die
    MT[ MT >= viral_lifetime ]     = 0     # not necessary but cleaner
    # let the long-infected users recover
    UV[ UT >= infection_duration ] = False
    UT[ UT >= infection_duration ] = 0  # not necessary but cleaner

    #random_infect(UV,UT,.1) # TODO: user infections from external sources

    # assign users to machines for this day
    # Will want to play with this
    um = random.choice(range(nmachines),nusers,replace=True)

    # transfer user germs to machines and set contamination times to 0
    # transfer machine germs to users and set infection times to 0
    # DONE Need to make infection time only set to zero if new infection
    # TODO Want to make viral transfer probabilistic: not certain
    # TODO Introduce random infection of users from external sources
    # DONE Make a immunity array, so that users gain immunity to virus after infection

    #print(UV.shape)
    for iu in range(nusers):
        user_viruses    = UV[iu].copy()       # save state to allow overwriting
        machine_viruses = MV[ um[iu] ].copy() # save state to allow overwriting
        MV[ um[iu] ][user_viruses]    = True  # contaminate machine with viruses on user
        MT[ um[iu] ][user_viruses]    = 0     # set age of contamination to 0
        # infect non-immune user with viruses on machine
        #print(UV.shape,user_viruses.shape,machine_viruses.shape,UI.shape)
        UV[ iu     ][machine_viruses  & logical_not(UI[iu])] = True 
        UT[ iu     ][machine_viruses  & logical_not(user_viruses)] = 0 # set infection age to 0
        UI[ iu, UV[iu] ] = True  # set user iu immune to future infection by current viruses

TRIAL_AXIS,DAY_AXIS,USER_AXIS,VIRUS_AXIS = 0,1,2,3
MACHINE_AXIS = 2
def bunch_of_trials1(nm=50):
    global uv,mv,ui,ut,mt
    global tr
    global nusers,nmachines
    nmachines = nm
    days = range(ndays+1)
    imax = 0

    uv = zeros((ntrials,ndays+1,nusers   ,nspecies),dtype=bool)  # uv[l,i,j,k] = True if in trial l on day i user j infected with species k
    mv = zeros((ntrials,ndays+1,nmachines,nspecies),dtype=bool)  # mv[l,i,j,k] = True if in trial l on day i machine j contaminated with species k
    ui = zeros(uv.shape,dtype=bool) # on day i user j immune to species k in trial l 
    ut = zeros(uv.shape,dtype=int)  # days since infection
    mt = zeros(mv.shape,dtype=int)  # days since contamination


    tr = TraceRecorder(30,ndays)

    for it in range(ntrials):

        uv[it,0,0,0] = True  # initial infection

        for t in range(ndays):
            step1(it,t)

        for iv in range(nspecies):
            tr.record(uv[it,:,:,iv].sum(axis=1) )
        imax = max(imax,uv.sum(axis=USER_AXIS).max())

    plt.figure(figsize=(20,8))
    plt.subplot(1,2,1)
    tr.imageplot(.5)
    plt.xlabel('day')
    plt.ylabel('number of users with virus')
    plt.title('darkess of dot indicates frequency of occurrence in trials');
    avg_infecteduserdays = uv.sum(axis=(TRIAL_AXIS,DAY_AXIS,USER_AXIS,VIRUS_AXIS))/ntrials
    print('Avg person-days of infection',avg_infecteduserdays)

    plt.subplot(1,2,2)
    ninfected = (uv.sum(axis=(DAY_AXIS,VIRUS_AXIS))>0).sum(axis=1)
    # histogram
    inffreqs = zeros(nusers+1)
    for n in ninfected:
        inffreqs[n] += 1
    plt.bar(range(nusers+1),inffreqs)
    plt.xlabel('number of users catching virus over entire time')
    plt.ylabel('frequency');
  nspecies = 1
nusers = 35
ndays = 21
viral_lifetime = 3
infection_duration = 2
ntrials = 500
vcolors = 'rgbmcy'

bunch_of_trials1(nm=60)
plt.show()
def step2(itrial,t,externalp):
    global nusers,nmachines,mv,mt,uv,ut,ui,tr

    # create views on next day's array slices
    MV = mv[itrial,t+1,:,:]
    MT = mt[itrial,t+1,:,:]
    UV = uv[itrial,t+1,:,:]
    UT = ut[itrial,t+1,:,:]
    UI = ui[itrial,t+1,:,:]
    
    # copy current state to day t+1
    MV[:,:] = mv[itrial,t,:,:]
    MT[:,:] = mt[itrial,t,:,:]
    UV[:,:] = uv[itrial,t,:,:]
    UT[:,:] = ut[itrial,t,:,:]
    UI[:,:] = ui[itrial,t,:,:]
    
    # then modify state on day t+1
    MT[MV] += 1 # age the germs on the infected machines by 1
    UT[UV] += 1 # advance the infections in the users by 1
     
    MV[ MT >= viral_lifetime ]     = False # let the old viruses on the machines die
    MT[ MT >= viral_lifetime ]     = 0     # not necessary but cleaner
    # let the long-infected users recover
    UV[ UT >= infection_duration ] = False
    UT[ UT >= infection_duration ] = 0  # not necessary but cleaner

    um = random.choice(range(nmachines),nusers,replace=True)
    external = random.rand(nusers,nspecies)< externalp

    for iu in range(nusers):
        user_viruses    = UV[iu].copy()       # save state to allow overwriting
        machine_viruses = MV[ um[iu] ].copy() # save state to allow overwriting
        if random.rand(1) > 0.3:
            MV[ um[iu] ][user_viruses]    = True  # contaminate machine with viruses on user, randomly, only 70% of the time
        MT[ um[iu] ][user_viruses]    = 0     # set age of contamination to 0
        # infect non-immune user with viruses on machine
        #print(UV.shape,user_viruses.shape,machine_viruses.shape,UI.shape)
        UV[external & logical_not(UI[iu])] = True 
        UV[ iu     ][machine_viruses  & logical_not(UI[iu])] = True
        UT[ iu     ][machine_viruses  & logical_not(user_viruses)] = 0 # set infection age to 0
        UI[ iu, UV[iu] ] = True  # set user iu immune to future infection by current viruses
TRIAL_AXIS,DAY_AXIS,USER_AXIS,VIRUS_AXIS = 0,1,2,3
MACHINE_AXIS = 2
def bunch_of_trials2(nm,externalp,doplot = False):
    global uv,mv,ui,ut,mt,tr,nusers,nmachines,q,qx
    global ax1
    global ax2
    global ax3
    nmachines = nm
    days = range(ndays+1)
    imax = 0

    uv = zeros((ntrials,ndays+1,nusers   ,nspecies),dtype=bool)
    mv = zeros((ntrials,ndays+1,nmachines,nspecies),dtype=bool)
    ui = zeros(uv.shape,dtype=bool) # on day i user j immune to species k in trial l 
    ut = zeros(uv.shape,dtype=int)  # days since infection
    mt = zeros(mv.shape,dtype=int)  # days since contamination
    q = []

    tr = TraceRecorder(30,ndays)

    for it in range(ntrials):

        uv[it,0,0,0] = True  # initial infection

        for t in range(ndays):
            step2(it,t,externalp)

        for iv in range(nspecies):
            tr.record(uv[it,:,:,iv].sum(axis=1) )
        imax = max(imax,uv.sum(axis=USER_AXIS).max())
    
    avg_infecteduserdays = uv.sum(axis=(TRIAL_AXIS,DAY_AXIS,USER_AXIS,VIRUS_AXIS))/ntrials
    ninfected = (uv.sum(axis=(DAY_AXIS,VIRUS_AXIS))>0).sum(axis=1)
    q.append(ninfected)
    
    if doplot:
        if nms == 30:
            
            plt.figure(figsize=(20,8))
            plt.subplot(1,2,1)
            tr.imageplot(.5)
            plt.xlabel('day')
            plt.ylabel('number of users with virus')
            plt.title.set_text('darkess of dot indicates frequency of occurrence in trials');
            print('Avg person-days of infection',avg_infecteduserdays)   
            
            plt.subplot(1,2,2)
            inffreqs = zeros(nusers+1)
            for n in ninfected:
                inffreqs[n] += 1
            plt.bar(range(nusers+1),inffreqs)
            plt.xlabel('number of users catching virus over entire time')
            plt.ylabel('frequency');

            plt.show()
nspecies = 1
nusers = 35
ndays = 21
viral_lifetime = 3
infection_duration = 2
ntrials = 5000
vcolors = 'gbmrck'
q = []                    # empty list of, after X trials, the list of how many N Infected for X trials
e = []                    # empty list of number of machines
qx = []                   # empty list of average of number of people who got sick, mean of q, so mean of 500 trials
ebu = []
qxu = []

plt.clf()
fig = plt.figure(figsize = (20,8))

gs = gridspec.GridSpec(1, 2)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[0,1])
ax1.title.set_text('Plot 1, Number of Machines VS Total Number of Users Infected')
ax2.title.set_text('Plot 2, Number of Machines Per Total Users VS Total Number of Users Infected Per Total Users')
ax1.legend(['External Sickness Prob. = 1%','External Sickness Prob. = 0.1%','External Sickness Prob. = 0.01%'])
ax2.legend(['External Sickness Prob. = 1%','External Sickness Prob. = 0.1%','External Sickness Prob. = 0.01%'])

extp = [0.01,0.001,0.0001]

for i,color in zip(extp,vcolors):
    for nms in [10,12,14,16,20,22,24,26,28,30,35,40,50,60,70,90,100]:
        bunch_of_trials2(nms,i,doplot = False)
        e.append(nms)
        ebu.append(nms/nusers)
        qx.append(np.mean(q))
        qxu.append(np.mean(q)/nusers)
        ax1.plot(array(e),np.transpose(qx),color)
        ax2.plot(array(ebu),np.transpose(qxu))
    e = []
    q = []
    qx = []
    ebu = []
    qxu = []

plt.show();
