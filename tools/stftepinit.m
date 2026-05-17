function p=stftepinit()
% this function initializes all the stfte pameters to their default values
%
p.pitchlim=[40 50 400];     % gs_frames:    [min target max] pitch in Hz
p.gcifrac=0.3;              % gs_frames:    fractional position of GCI in analysis frame
p.GCImethod='DYPSA';        % gs_frames:    GCI detection algorithm {'DYPSA','YAGA','SEDREAMS'}
p.smoothalg='none';         % smoothframes: smoothing algorithm {'none','lin','quadlin','quadlog'}
p.voithresh=0.6;            % smoothframes: cross-correlation threshold for voiced frames
p.voilenrat=1.1;            % smoothframes: maximum length ratio of successive voiced frames
p.offset='none';            % stfte:        offset removal: {'none','mean','ends'}
p.scale='none';             % stfte:        scaling method: {'none','peakabs','rms','len','sqlen'}
p.pad='none';               % stfte:        zero-padding method: {'none','zero','ends'}
p.groupdelay='none';        % stfte:        linear phase component: {'none','dct','ewgd','cplx','phgr','gcif','gpdf','fmnb' + optional 'int' suffix}
p.gcifrac=0;                % stfte:        fraction of frame length for p.groupdelay='gpdf' option
p.gpdfrac=0.3;              % stfte:        fraction of frame length for p.groupdelay='gpdf' option
p.fmbound=[0.3 0.5];        % stfte:        bounds for group delay option p.groupdelay='fmbd' as fraction of frame length
p.window='r';               % stfte:        window: 'r'=rectangular, 'n'=hanning, 'm'=hamming
p.windowmode='Es';          % stfte:        window mode (see v_windows.m)
p.interpdom='magcph';       % stftegridz:   Interpolatione domain: {'cplx','magcph','crmcph'}
p.interpext='rep';          % stftegridz:   Handling of extrapolated frames: {'omit','zero','rep','refl'}
p.interpstft='indep';       % stftegridz:   Interpolation method for call to griddata: {'none','indep','nearest','linear','natural','cubic','v4'}
p.interpfsps=10e-6;         % stftegridz:   Dimensionless interpolation scale factor: fs^-2 multiplied by the distance in Hz that is equivalent to a distance of one second
p.interpgd='none';          % stftegridz:   Interpolation of frame group delay: {'none','lin','linrep'}
p.interpof='none';          % stftegridz:   Interpolation of frame offset: {'none','lin'}
p.interpsc='none';          % stftegridz:   Interpolation of frame scale factor: {'none','lin','log'}
p.interptz='none';          % stftegridz:   Compensation of phase for the frame starting sample: {'none','origin'}
par.interpseq=  1;          % stfteresft:   0=interpolate T and F jointly, 1=interpolate T and F sequentially with an intermediate complex spectrum
par.interpp=    5;          % stfteresft:   length of intermediate filter