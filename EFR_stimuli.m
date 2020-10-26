function [cal_stimulus, calibs] = EFR_stimuli(stim)

% this function creates the stimulus
%
% Input:
%
% stim -(string): stimulus to generate:

%                           'tone_4_sin',
%                           'square_4_50_same_rms','square_4_50_same_ptp',
%                           'square_4_25_same_rms', 'square_4_25_same_ptp',
%                           'multi_harmonics_4','transposed'
%
% Output:
%
%   cal_stimulus: one instance of requested stimulus
%   calibs: structure with relevant calibration information for the
%           stim_test.m file
%
%
%
%

%% Set the general stimulus parameters

pars = struct();

pars.stim = stim;
pars.fs = 48000; % stimulus sampling rate
pars.mf = 120;   % modulation frequency
pars.avs = 1000; % repetitions 
pars.loops = 10; % number of loops to devide the stimuli in --> rem(pars.avs,pars.loops) needs to be 0!
pars.md =  0.95; % modulation depth in percent; 0.95 = 95% depth; in dB: -0.4455 --> %: 10^(-0.4455/20) ~ dB: 20*log10(0.95); exception: multi_harmonics_4: fixed at 1
%you may further increase the EFR strength by setting md to 1.
pars.jit = 0.1; % percentage of the silence gap to jitter
pars.silence = 100e-3; %silence gap in sec
pars.dur = 0.4; % duration of epoch in sec
pars.phi = 3*pi/2; % starting phase of modulator (sin); --> carrier (sin) phases are 0; exception: multi_harmonics_4 phase can not be changed here
pars.order = 1024; % Order of filter in stim_filter.m function
pars.level = 70; % Target output level
pars.f = 4000; % carrier frequency
% THE ONLY PARAMETER THAT NEEDS TO BE CHANGED IN THE LAB:
% The tone_4_sin condition was the reference condition to which all other conditions were calibrated.
% It needs to be calibrated to be at the level of pars.level_total
%%%%%%%%%%%%%%%%%%%%%%%% Calibrate reference condition
pars.level_total = 70; % Target level of overall output of reference condition
pars.lcal_total = 70;%; % Change this, so the reference is played at target dB SPL level in the lab
% --> if pars.level_total is equal to pars.lcal_total no amplification is applied
pars.ref_rms = 0.0087485551; % Calibrated on 180316 - Markus; Note the rms value of the reference and plug it in here
pars.ref_ptp = 0.0201700786; % Calibrated on 180316 - Markus;% Note the max(abs() value of the reference and plug it in here (for symetric signal)
%--> 60 second stimulus no silence at beginning or end
pars.corfactor = 0; % dB - correction factor for sound meter (if needed)
%%%%%%%%%%%%%%%%%%%%%%%%

switch stim    % select the desired modulator shape
    case 'tone_4_sin'
        pars.lcal_stim_mod = 109.7;% Calibrated on 180316 - Markus; Value to calibrate level of stimulus
        pars.modulator = 'sin'; % type of modulator to use

    case 'transposed'
        pars.modulator = 'transposed';
        
    case 'square_4_50_same_rms'
        pars.modulator = 'square'; % type of modulator to use
        
    case 'square_4_50_same_ptp'
        pars.modulator = 'square'; % type of modulator to use
        
    case 'square_4_25_same_rms'
        pars.modulator = 'square_25'; % type of modulator to use
        
    case 'square_4_25_same_ptp'
        pars.modulator = 'square_25'; % type of modulator to use
        
    case 'multi_harmonics_4'
        pars.modulator = 'multi_harmonics'; % type of modulator to use

    otherwise
        error('Non existing stimulus condition')       
end    


%% Generate the modulator
n  = round(pars.dur*pars.fs); % number of samples needed

switch pars.modulator
    
    case 'sin'
        %================ Modulator ==============================
        modsine = sin([0:n-1]'/pars.fs * 2*pi* pars.mf + pars.phi);
        %=========================================================
    case 'square_25' % 25 % duty cycle: This is the percent of the period in which the signal is positive.
        %================ Square wave Modulator ==================
        modsine = square([0:n-1]'/pars.fs * 2*pi* pars.mf + pars.phi,25);
        %=========================================================
    case 'square' % 50 % duty cycle: This is the percent of the period in which the signal is positive.
        %================ Square wave Modulator ==================
        modsine = square([0:n-1]'/pars.fs * 2*pi* pars.mf + pars.phi,50);
        %=========================================================
        %================ Alternative Square wave Modulator ==================
        %modsine = sign(sin([0:n-1]'/pars.fs * 2*pi* pars.mf + pars.phi));
        %=========================================================
    case 'multi_harmonics' % add up the first 10 harmonics of the modulation frequency with a start phase (cos) = 0 (starts at 1)
        %================ Modulator ==============================
        modsine =   cos([0:n-1]'/pars.fs * 2*pi*1*pars.mf)+...
                    cos([0:n-1]'/pars.fs * 2*pi*2*pars.mf)+...
                    cos([0:n-1]'/pars.fs * 2*pi*3*pars.mf)+...
                    cos([0:n-1]'/pars.fs * 2*pi*4*pars.mf)+...
                    cos([0:n-1]'/pars.fs * 2*pi*5*pars.mf)+...
                    cos([0:n-1]'/pars.fs * 2*pi*6*pars.mf)+...
                    cos([0:n-1]'/pars.fs * 2*pi*7*pars.mf)+...
                    cos([0:n-1]'/pars.fs * 2*pi*8*pars.mf)+...
                    cos([0:n-1]'/pars.fs * 2*pi*9*pars.mf)+...
                    cos([0:n-1]'/pars.fs * 2*pi*10*pars.mf);    
        %=========================================================
    case 'transposed' %implements the Van de Par & Kohlrausch (1997) version of the transposed tone
        %==========================================================
        % "modulator"   
        hrsine = sin(2*pi*pars.mf * [0:n-1]'/pars.fs);         
        % half-wave rectify the modulator
        hrsine = max(0,hrsine);             
        % Compute coefficients for low-pass filtering
        fcuth = .2*pars.f/(pars.fs/2);
        [b,a] = butter(4,fcuth,'low');
        % Low-pass filter the modulator
        modsine = filter(b,a,hrsine);  
        
    case 'none'
        % no modulation
        modsine = 'none';
        
    otherwise
        
        error('Non existing stimulus condition')
end

calibs.modsine = modsine; % store modulator in output variable 'calibs'

%% generate the stimuli
switch stim   
    
    case {'tone_4_sin','square_4_50_same_rms','square_4_50_same_ptp'}
               
        % Create pure tone stimulus
        stim_raw = sin([0:n-1]'/pars.fs * 2*pi* pars.f); % sine with phase 0 (starts at 0)
        
        % Modulate the pure tone
        stim_mod = stim_raw .* (1 + pars.md * modsine);
        
    case {'transposed'}
        % Create pure tone stimulus
        stim_raw = sin([0:n-1]'/pars.fs * 2*pi* pars.f); % sine with phase 0 (starts at 0)
        
        % Add the carrier*modulator to the transposed tone
        stim_mod = stim_raw .* (1 + pars.md * (2 * modsine-1));
        %the modsine is half wave rectified (0-1), in stead of SAM (-1 to +1) so we need to do the
        %modulation differently than for a pure tone.
        
        
    case 'multi_harmonics_4'
        
        % Create pure tone stimulus
        stim_raw = sin([0:n-1]'/pars.fs * 2*pi* pars.f); % sine with phase 0 (starts at 0)
        
        % Modulate the pure tone
        stim_mod = stim_raw .* (1 + 1 * modsine); % changing the modulation depth makes no sense because we 
        % overmodulate. So it is set to 1
        
    case {'square_4_25_same_rms','square_4_25_same_ptp'}
        
        % Create pure tone stimulus
        stim_raw = sin([0:n-1]'/pars.fs * 2*pi* pars.f); % sine with phase 0 (starts at 0)
        
        % Modulate the pure tone
        stim_mod = stim_raw .* (2 + (pars.md*2) * modsine);
        
    otherwise
        
        error('Non existing stimulus condition')
        
end


%% Calibrate the stimuli

window = tukeywin(length(stim_mod), 0.025); % Create a Tukey (tapered cosine) window
   
wind_stim_mod = stim_mod.*window;

switch stim
 
% Reference stimulus
    case {'tone_4_sin'} % reference stimulus
    
        % Calibrate modulated stimulus
        calfact_stim_mod = 10^((pars.level-pars.lcal_stim_mod)/20); % Computer needed RMS to achieve a target outputlevel of 'level' using 'lcal' to calibrate
        cal_stimulus = wind_stim_mod * calfact_stim_mod; % multiply stimulus with target RMS
        
        
% same RMS
    case {'square_4_50_same_rms','square_4_25_same_rms', 'multi_harmonics_4','transposed'}
        stim_mod_norm = wind_stim_mod/rms(wind_stim_mod); % normalize RMS to 1
        cal_stimulus = stim_mod_norm * pars.ref_rms; % multiply with reference RMS
        
% same PTP
    case {'square_4_50_same_ptp','square_4_25_same_ptp'}
        
        stim_mod_norm = wind_stim_mod/max(abs(wind_stim_mod)); % normalize the peak to be 1
        cal_stimulus = stim_mod_norm* pars.ref_ptp; % multiply by reference peak value
        
end
 
%% Save data in output variable 'calibs' to be used in the 'stim_test.m' file for calibration and plotting

calibs.stim_raw = stim_raw; % raw stimulus

calibs.stim_mod = stim_mod; % Modulate stimulus

% --> the 'cal_stimulus' output stimulus is its own output variable
end

%Code written by Markus Garrett, Viacheslav Vasilkov, and
%beloging to the recordings presented in the publication: "Enhancing the sensitivity of the
%envelope-following response for cochlear synaptopathy screening in humans:
%the role of stimulus envelope" by Viacheslav Vasilkov, Markus Garrett,
%Manfred Mauermann and Sarah Verhulst.
%see license file in the code repository


