# EEG-Artifact-Reduction
In this repository, some blind and semi-blind source separation methods, such as PCA, GEVD, ICA, and DSS, have been implemented to eliminate the EEG artifacts.
## Question 1  
- Ex1.mat data was scattered.
- the PCA algorithm was implemented on the data, and the principal components were extracted. Moreover, the whitened data was plotted. 
- PCA algorithm was implemented with the Matlab function, and the result was similar to the results obtained in the previous part.
- SVD decomposition was used, and checked the connections of PCA and SVD.
## Question 2  
- Ex2.m contains non-epileptic simulated data, which we contaminated with background EEG signal and muscle noise.
- Sources were extracted using PCA and Com2 (on of ICA algorithms) and undesirable sources were eliminated.
- Denoised data was transfered to sensor space. Furthermore, RRMSE was calculated for the two methods with different SNRs. 
## Question 3
- Ex3.m includes seizure/non-seizure periods of recorded EEG signals from epileptic patients.
- Main artifacts were detected and Com2 was applied to detect sources.
- Fequency and spatial charactristics (topoplots) were obtained.
- Undesired sources were removed based on the temporal, spatial and frequency characteristics, and denoised signal was constructed.
## Question 4

## Question 5
- pure.mat contains a 19-channel EEG signal of a subject with closed eyes. 
- contaminted.m includes artificially contaminated version of pure.mat. 
- The main artifact is electrooculogram signal that is mixed up with EEG signal. 
    -   




