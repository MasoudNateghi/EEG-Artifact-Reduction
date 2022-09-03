X_org  % 32-channel data
load('Electrodes') ;

% Plot Data
% Use function disp_eeg(X,offset,feq,ElecName)
offset = max(abs(X_org(:))) ;
feq = 250 ;
ElecName = Electrodes.labels ;
disp_eeg(X_org,offset,feq,ElecName);