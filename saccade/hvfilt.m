function hvfilt(datfiles,outfiles, Fcutoff,Fsample)
% Digital low-pass filtering on HVFILES
%
% HVFILT(HVFILES,OUTFILES,FCUTOFF)
%
% Digital low-pass filtering on HVFILES (with extension HV).
% File names have to be entered as ['AAA';'BBB';'CCC'], with or without the
% '.hv'-extension. If HVFILES are not entered, HVFILT will by default
% filter all hvfiles in the current directory.
% Optional input includes:
% OUTFILES - by default the HVFILES are overwritten!!! If you do
% not want this to occur you have to supply the new name in OUTFILES.
% FCUTOFF - Cut-off frequency (Hz), by default this will amount to 80 Hz.
%
% See also CALIBRATE, SACDET, SAC2MAT

%% Initialization
if nargin<1
    d                   = dir('*.hv');
    datfiles            = char(d.name);
end
if nargin <2
    outfiles            = datfiles;
end
nchan                   = 2;        % Azimuth channel and elevation channel = 2 channels
h                       = 1;        % Channel number 1 = azimuth
v                       = 2;        % Channel number 2 = elevation

%% Design filter
if nargin<3
    Fcutoff             = 80;       % Cutoff Frequency  (Hz)
end
Order                   = 50;
if nargin<4
	Fsample                 = 6104;      % Sample Frequency
% 	Fsample				= 1017; % Hoop
end
ext                     = '.hv';

lpFilt = designfilt('lowpassfir', 'FilterOrder', Order, 'CutoffFrequency', ...
                    Fcutoff, 'SampleRate', Fsample, 'Window', 'hamming');
%% Looping all files
for i                   = 1 : size(datfiles,1)
    fname               = datfiles(i,:);
    fname               = pa_fcheckext(fname,ext);
    csvfile             = pa_fcheckext(fname,'csv');
    [~,chaninfo]		= pa_readcsv(csvfile);
    nsample             = chaninfo(1,6);
    if nsample/3>Order
        disp(['   Reading ' fname]);
        fname               = pa_fcheckexist(fname);
        fid                 = fopen(fname,'r','l');
        if fid==-1
            disp(['   Error reading ' fname ]);
            return;
        end;
        [mtx,n]             = fread(fid,[nchan,inf],'float');
        fclose(fid);
        % filter the signals trial by trial
        disp(['   Low-Pass Fc = ' int2str(Fcutoff) ' Hz   Order = ' int2str(Order)] );
        nblock      = n/(nsample*nchan);
        for j       = 1:nblock
            indx        = (j-1)*nsample+1:j*nsample;
            fs          = filtfilt(lpFilt,mtx(h,indx));
            mtx(h,indx) = fs;
            fs          = filtfilt(lpFilt,mtx(v,indx));
            mtx(v,indx) = fs;
			
		end
        %% save filtered data
        fname       = outfiles(i,:);
        fname       = pa_fcheckext(fname,ext);
        disp(['   Writing ' fname]);
        fid         = fopen(fname,'w','l');
        if fid==-1
            disp(['   Error writing ' fname ]);
            return;
		end
        MTX         = [mtx(h,:); mtx(v,:)];
        MTX         = MTX(:);
        fwrite(fid,MTX,'float');
        fclose(fid);
    elseif nsample/3<Order
        disp(['skip ' fname])
        warning('HVfilt:OrderotEnoughSamples','Data is not filtered due to small number of samples');
    end
end



