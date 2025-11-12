#!/usr/bin/python

#Author: Juan de Vicente
#17/01/2015
#photoz.py
#v 1.0
#Usage: python photoz.py

#Program to test the functions of nf (Neigborhood Fit) module

import math
import numpy as np
import os
import sys
import fnmatch
#import pylab
import sys
import scipy.stats as st
#from sklearn import neighbors
import dnf
#import knnz
#import nipknnz
##import anf
#import ncf
#import matplotlib.pyplot as plt

#from sklearn import svm,datasets
#from sklearn.cross_validation import train_test_split
from sklearn.metrics import confusion_matrix
#import matplotlib.cm as cm
from scipy.stats import ks_2samp,entropy

import astropy

from sklearn import neighbors
import numpy.lib.recfunctions as rec
#import weighted_kde as kde
#import pandas as pd

from astropy.table import Table


#parameters
nfilters=4


#####reading training file ######################

GALAXY=Table.read(sys.argv[1]) #'Y1SPEC_MATCH2as_Y1Goldv103.fits')
#GALAXY=GALAXY[:10000]
Ngalaxies=len(GALAXY)
#print 'Ngalaxies=',Ngalaxies
Ngalaxies=len(GALAXY)
print('Ngalaxies=',Ngalaxies)

G=np.zeros((Ngalaxies,nfilters),dtype='double')
Gerr=np.zeros((Ngalaxies,nfilters),dtype='double')

G[:,0]=GALAXY['BDF_MAG_G_CORRECTED']
G[:,1]=GALAXY['BDF_MAG_R_CORRECTED']
G[:,2]=GALAXY['BDF_MAG_I_CORRECTED']
G[:,3]=GALAXY['BDF_MAG_Z_CORRECTED']
#G[:,4]=GALAXY['BDF_MAG_Y_CORRECTED']


Gerr[:,0]=GALAXY['BDF_MAG_ERR_G']
Gerr[:,1]=GALAXY['BDF_MAG_ERR_R']
Gerr[:,2]=GALAXY['BDF_MAG_ERR_I']
Gerr[:,3]=GALAXY['BDF_MAG_ERR_Z']
#Gerr[:,4]=GALAXY['BDF_MAG_ERR_Y']
#Gerr[:,5]=GALAXY['BDF_T_ERR']



Ntrain=Ngalaxies
TRAIN=GALAXY
T=G
Terr=Gerr

#read valid
GALAXY=Table.read(sys.argv[2]) 
#GALAXY=GALAXY[:10000]
#print 'GALAXY=',GALAXY
Ngalaxies=len(GALAXY)
print('Ngalaxies=',Ngalaxies)

G=np.zeros((Ngalaxies,nfilters),dtype='double')
Gerr=np.zeros((Ngalaxies,nfilters),dtype='double')

G[:,0]=GALAXY['BDF_MAG_G_CORRECTED']
G[:,1]=GALAXY['BDF_MAG_R_CORRECTED']
G[:,2]=GALAXY['BDF_MAG_I_CORRECTED']
G[:,3]=GALAXY['BDF_MAG_Z_CORRECTED']
#G[:,4]=GALAXY['BDF_MAG_Y_CORRECTED']
#G[:,5]=GALAXY['BDF_T']


#GM[:,0]=-2.5*np.log10(np.abs(GALAXY['BDF_FLUX_CORRECTED_G']))
#GM[:,1]=-2.5*np.log10(np.abs(GALAXY['BDF_FLUX_CORRECTED_R']))
#GM[:,2]=-2.5*np.log10(np.abs(GALAXY['BDF_FLUX_CORRECTED_I']))
#GM[:,3]=-2.5*np.log10(np.abs(GALAXY['BDF_FLUX_CORRECTED_Z']))
#GM[:,4]=-2.5*np.log10(np.abs(GALAXY['BDF_FLUX_CORRECTED_Y']))


Gerr[:,0]=GALAXY['BDF_MAG_ERR_G']
Gerr[:,1]=GALAXY['BDF_MAG_ERR_R']
Gerr[:,2]=GALAXY['BDF_MAG_ERR_I']
Gerr[:,3]=GALAXY['BDF_MAG_ERR_Z']
#Gerr[:,4]=GALAXY['BDF_MAG_ERR_Y']
#Gerr[:,5]=GALAXY['BDF_T_ERR']


Nvalid=Ngalaxies
VALID=GALAXY
V=G
Verr=Gerr
#################################### 
      


#bins
#start=0.0
#stop=0.8
#step=0.1

#start
#step=0.066

start=0.0
stop=1.6
step=0.01

#start=0.0
#stop=2.0
#step=0.1

#start=0.0
#stop=0.9
#step=0.1

#start=0.1
#stop=0.7
#step=0.0375

#start=np.double(raw_input('start:'))
#stop=np.double(raw_input('stop:'))
#step=np.double(raw_input('step:'))

zbins=np.arange(start,stop,step)
nbins=len(zbins)-1 
bincenter=(np.double(zbins[1:])+np.double(zbins[:-1]))/2.0

#reponer
#zbins = np.linspace(0, 2.0, 50) #201)
#bincenter = (zbins[0:-1] + zbins[1:])/2.0
#nbins=len(zbins)-1

binning = zbins  #np.linspace(0, 2.0, 201)
bin_centers = (binning[0:-1] + binning[1:])/2.0


algorithm=sys.argv[3]
#algorithm=raw_input('\nenf\ndnf\nanf\nEnter an option:')

#names=('z_photo','z1','zerr','zerrabs','zerr_e','mode_z','mean_z','sample_z','std_z','SNmean','SNmax'

#PHOTOZ CALL
z_photo,zerr_e,photozerr_param,photozerr_fit,Vpdf,z1,nneighbors,de1,d1,id1,C=dnf.dnf(T,TRAIN['Z'],V,Verr,zbins,pdf=False,Nneighbors=80,bound=False,radius=2,magflux='mag',metric=algorithm,coeff=True) 

print("mean Nneighbors=",np.mean(nneighbors))


#SAVE RESULTS
#*********point prediction file**********
#f.writeto('test.fits')
#create the test fits files 
from astropy.table import Table
#d = {} #dictionary
d=VALID

d['DNF_Z']=z_photo #+algorithm.lower()]=z_photo
d['DNF_ZN']=z1 #+algorithm.lower()]=z1
d['DNF_ZSIGMA']=zerr_e
d['DNF_ZERR_PARAM']=photozerr_param
d['DNF_ZERR_FIT']=photozerr_fit
d['DNF_D1']=d1
d['DNF_NNEIGHBORS']=nneighbors
#d['Z']=VALID['Z']
d['DNF_DE1']=de1
d['DNF_ID1']=id1
#d['M']=V
#d['M1']=T[id1]
d['C']=C
# d['Vpdf']=Vpdf

#d['Vpdf_start']=str(start)
#d['Vpdf_stop']=str(stop)
#d['Vpdf_step']=str(step)

#d['numberOfNeighgors']=nneighbors
#d['closestDistance']=closestDistance
#d['zerr_'+algorithm.lower()]=zerr_e
#d['nneighbors_'+algorithm.lower()]=nneighbors
#d['closestDistance_'+algorithm.lower()]=closestDistance
#print 'd=',d
#fit = Table(d)
#fit.write('test1wlenf/jvicente'+'_'+algorithm+'_'+testfile, format='fits',overwrite=True)
d.write(sys.argv[4], format='fits',overwrite=True)


#raw_input()


sys.exit()




#####Signal to noise analysis
Msnm=np.zeros(Nvalid,dtype='double')
Msnstd=np.zeros(Nvalid,dtype='double')
Msnsum=np.zeros(Nvalid,dtype='double')
Msnmax=np.zeros(Nvalid,dtype='double')
Msnmin=np.zeros(Nvalid,dtype='double')
Msnmod=np.zeros(Nvalid,dtype='double')

####z1 results
#photoz error
zerr1=z1-VALID['REDSHIFT']
zerrabs1=np.abs(z1-VALID['REDSHIFT'])  
print('z1err results')
print('mean=',zerr1.mean())
print('median=',np.median(zerr1))
print('std=',zerr1.std())
print('mad=',zerrabs1.mean())
print('biasNorm=',(zerr1/zerr_e).mean())
print('sigmaNorm=',(err1/zerr_e).std())

zerr1Sort=np.sort(zerrabs1)
sigma68=zerr1Sort[Nvalid*68/100]
print('sigma68=',sigma68)
aux=np.where(np.sqrt(z_true_hist)==0.0,0.0,(z_true_hist-z_photo_hist)/np.sqrt(z_true_hist))
print('Npoisson_Nz=',np.linalg.norm(aux)/np.sqrt(nbins))
print('kS=',ks_2samp(z_true_hist,z_photo_hist)[0])

aux=np.where(np.sqrt(z_true_hist)==0.0,0.0,(z_true_hist-z_1_hist)/np.sqrt(z_true_hist))
print('Npoisson_Nz_1=',np.linalg.norm(aux)/np.sqrt(nbins))
print('kS_1=',ks_2samp(z_true_hist,z_1_hist)[0])

#print 'Npoisson_Nz=',(np.linalg.norm((z_true_hist-z_photo_hist)/np.sqrt(z_true_hist)))/np.sqrt(nbins)
print('Chi2_Nz=',np.linalg.norm((z_true_hist-z_photo_hist)/np.sqrt(z_photo_hist))) #/np.sqrt(nbins)

#END
sys.exit()
