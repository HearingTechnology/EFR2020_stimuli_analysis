function [ EFR ] = EFR_extractor( epochs, mf, varargin)
% This function takes an array of size: epochs x time or channels x epochs
% x time and extracts different temporal envelope measures using a
% bootstrapping approach. If the epochs alternate in polarity this will be
% accounted for in the epochs drawing procedure.The random number generator
% is reset for every bootstrap draw, which leads to reproducable results:
%
% - Spectral magnitude
% - Spectral noise floor estimate
% - Peak-to-noise floor estimate
% - Phase locking value
% - Noise floor correctd frequency-time domain reconstruction of the neural envelope
% - Uncorrected frequency-time domain reconstruction of the neural envelope
% - Average cycle of modulation frequency
%
% Metrics representing the envelope coding at the modulation frequency
% are based on the energy of the modulation frequeny and its
% harmonics.
%
% Required input arguments:
% epochs:   EEG data set of shape:
%           epochs x time or
%           channels x epochs x time
%
% mf:       Modulation frequency of the recorded signal in Hz (e.g. 120)
%
% Optional input arguments:
% channels: integer of vector specifying the channels to use e.g. 32 or
%           [20:32]. Waring: Very computationally expensive for large
%           number of channel; default: use all channels
%
% fs:       Integer - Sampling rate in Hz; default: 16384

% harms:    Integer - Number of harmonics including mf to use for the metrics;
%           default: 5
%
% bs_draws: Integer - Number of epochs to use per bootstrap draw (should be even
%           number); default is the number of epochs in the input signal
%
% bs_sig:   Integer - Number of bootstrap iterations for the signal estimate;
%           defautl: 200
%
% bs_noise: Integer - Number of bootstrap iterations for the noise floor estimate;
%           defautl: 1000
%
% zeropad: Use number of samples for fft (0) or zeropad signal (1) so nfft
%          is equal to 2^nextpow2(samples in signal); default: 0
%
% tukey_width: Width of the window (0 to 1) applied to the epochs befor fft (see
%             'tukeywin'); default: 0.02 (2% of signal length)
%
% visual:  Integer - Shows a visual represnetation of the results for one of
%          the channel in 'channels'. e.g. if channels 14 and 15 were specified
%          one can choose to plot either input channel 14 or 15;
%          default: no plotting
%
% Output:
% EFR - structure including all extracted information
%
%
% For further information contact:
%
% sarah.verhulst@ugent.be
%
% Authors: Markus Garrett, markus.garrett@uol.de; Viacheslav Vasilkov, Viacheslav.Vasilkov@UGent.be
% A UGent-noncommercial license applies: details in the License.txt file.  


%% Check inputs
tStart = tic; % start timer

p = inputParser;

checkepoch = @(epochs) (ndims(epochs) == 2 || ndims(epochs) == 3); % data set
addRequired(p, 'epochs', checkepoch);

check_mf = @(mf) (mf == floor(mf)); % modulation frequency
addRequired(p,'mf', check_mf);

default_channels = []; % What channels to use
check_channels = @(channels) (isvector(channels));
addParameter(p,'channels',default_channels, check_channels);

default_fs = 16384; % Sampling rate
check_fs = @(fs) (fs == floor(fs));
addParameter(p,'fs',default_fs, check_fs);

default_harms = 5; % Number of harmonics to include in computation
check_harms = @(harms) (harms == floor(harms));
addParameter(p,'harms',default_harms, check_harms);

default_bs_draws = []; % Number of epochs per bootstrap draw
check_bs_draws = @(bs_draws) (bs_draws == floor(bs_draws));
addParameter(p,'bs_draws',default_bs_draws, check_bs_draws);

default_bs_sig = 200; % number of signal bootstrap runs
check_bs_sig = @(bs_sig) (bs_sig == floor(bs_sig));
addParameter(p,'bs_sig',default_bs_sig, check_bs_sig);

default_bs_noise = 1000; % number of noise floor bootstrap runs
check_bs_noise = @(bs_noise) (bs_noise == floor(bs_noise));
addParameter(p,'bs_noise',default_bs_noise, check_bs_noise);

default_zeropad = 0; % zeropad or not
check_zeropad = @(zeropad) (zeropad == 0 || zeropad == 1);
addParameter(p,'zeropad',default_zeropad, check_zeropad);

default_tukey_width = 0.02; % width of tukey window
check_tukey_width = @(tukey_width) (tukey_width >= 0 && tukey_width <= 1);
addParameter(p,'tukey_width',default_tukey_width, check_tukey_width);

default_visual = []; % plot channel or not
check_visual = @(visual) (visual == floor(visual));
addParameter(p,'visual',default_visual, check_visual);

parse(p, epochs, mf, varargin{:});
epochs = p.Results.epochs;
mf = p.Results.mf;
channels = p.Results.channels;
fs = p.Results.fs;
harms = p.Results.harms;
bs_draws = p.Results.bs_draws;
bs_sig = p.Results.bs_sig;
bs_noise = p.Results.bs_noise;
zeropad = p.Results.zeropad;
tukey_width = p.Results.tukey_width;
visual = p.Results.visual;

%% Adjust input if necessary


% add addtional singelton dimension if epochs is of shape epochs x time
if ndims(epochs)  == 2
    epochs = reshape(epochs,1,size(epochs,1),size(epochs,2));
end

if isempty(bs_draws)
    bs_draws = size(epochs,2); % set number of epochs if not specified
end

if ~isempty(channels)
    epochs = epochs(channels,:,:); % reduce data to channels specified
else
    channels = 1:size(epochs,1); % Note that all channels were used
end

% Make sure the number of draws per epoch is an even number
if rem(bs_draws,2) == 1
    bs_draws = bs_draws - 1;
    warning('Number of epochs per bootstrap draw was not an even number it has been changed from %d to %d.',bs_draws-1,bs_draws)
end


%% Extract basic parameters
chans = size(epochs,1); % number of channels
eps = size(epochs,2); % number of epochs
N = size(epochs,3); % time points in epoch

if zeropad == 0
    NFFT = N; % number of samples in signal
elseif zeropad == 1
    NFFT = 2^nextpow2(N); % zeropad the time signal to the next power of 2
end

freq = 0:fs/NFFT:fs-1/NFFT;   % Frequency vector
time =(0:N-1)/fs;               % Time vector
freq_res = fs/NFFT; % frequency resolution
samples_per_cycle = fs/mf; % exact number of samples per cycle (most likely not an integer)
number_of_cycles = floor(mf*(N/fs)); % number of whole cycles per epoch
cut_off = 0.2; % Amount of the time domain signal to cut off both sides to find the max and min peaks

% design butterworth bandpass filter of 4th order to smoothen the average cycle
d = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',round(mf-mf/3),'HalfPowerFrequency2',round(mf*harms+mf), ...
    'SampleRate',fs);

closesed_match = zeros(1,harms); % frequency index that most closely matches harmonics
for i = 1:harms
    [ ~ , closesed_match(i) ] = min( abs( freq-i*mf ) );
end

%% EEG preprocessing:

baseline = mean(epochs,3); % extract the average value of every epoch - baseline
epochs = bsxfun(@minus,epochs,baseline); % substract baseline from every epoch

%% Extract indices of trials with positive and negative polarities
% given that they alternated in the experimental design

pos_index = 2:2:eps; % positive polarity epochs
neg_index = 1:2:eps; % negative polarity epochs

%% Estimate noise floor
tmp_nf = zeros(size(epochs,1),bs_noise,ceil(NFFT/2)+1); % preallocate memory to speed things up

for j = 1:bs_noise % for all noise floor draws
    fprintf('\n:: Noise floor - Draw %d out of %d',j,bs_noise)
    rng(bs_noise-j) % change random number generator for each bootstrap iteration --> keeps data reproducable
    pos_draws = randi(size(pos_index,2),1, bs_draws/2); % choose epochs with replacement
    neg_draws = randi(size(neg_index,2),1, bs_draws/2); % choose epochs with replacement
    tmp_epochs = epochs(:,[pos_index(pos_draws), neg_index(neg_draws)],:);  % draw epochs (same amount of positive and negative epochs)
    tmp_epochs(:,1:floor(size(pos_draws,2)/2),:) = tmp_epochs(:,1:floor(size(pos_draws,2)/2),:)*-1; % multiply half of the positive polarity epochs by -1 --> flip phase
    tmp_epochs(:,(end-floor(size(neg_draws,2)/2)+1):end,:) = tmp_epochs(:,(end-floor(size(neg_draws,2)/2)+1):end,:)*-1; % multiply half of negative polarity epochs by -1 --> flip phase
    tmp_epochs = reshape(mean(tmp_epochs,2),chans,N); % Average the epochs in the time domain and remove epoch dimnesion
    tmp_fft = fft(tmp_epochs,NFFT,2); % compute FFT
    tmp_nf(:,j,:) = tmp_fft(:,1:ceil(NFFT/2)+1); %discard second half of spectrum
end

nf_spec = tmp_nf./N; % devide the spectrum by the number of DATA POINTS PER EPOCH



%% Estimate signal

% preallocate memory to speed things up
ptn_efr = zeros(size(epochs,1),bs_sig); % collect the time domain EFRs (amplitude)
ptn_efr_uncor = zeros(size(epochs,1),bs_sig); % collect the uncorrected time domain EFRs (amplitude)
ptn_td = zeros(size(epochs,1),bs_sig,N); % collect the time domain reconstruction (whole signal)
ptn_td_uncor = zeros(size(epochs,1),bs_sig,N); % collect the uncorrected time domain reconstruction (whole signal)
plv = zeros(size(epochs,1),bs_sig,ceil(NFFT/2)+1); % collect the phase locking values per frequency (whole signal)
spec = zeros(size(epochs,1),bs_sig,ceil(NFFT/2)+1); % collect the frequency domain spectrum (whole signal)
td = zeros(size(epochs,1),bs_sig,N); % collect the time domain signal (whole signal)
cycles = zeros(size(epochs,1),bs_sig, ceil(samples_per_cycle)); % Collect average cycles of the epoch
cycle_amp = zeros(size(epochs,1),bs_sig); % collect each (max-mix)/2 cycle amplitude per channel
nf_bins =  zeros(size(epochs,1),bs_sig,harms); % collect instentanious noise floor per harmonic
plv_nf_bins = zeros(size(epochs,1),bs_sig,harms); % collect instentanious plv noise floor per harmonic
plv_ptn = zeros(size(epochs,1),bs_sig,harms); % % collect instentanious plv ptn values of the harmonics
ftd_sum_harm_uncor = zeros(size(epochs,1),bs_sig); % collect uncorrected summed harmonic magnitudes of freq/time domain approach
ftd_sum_harm = zeros(size(epochs,1),bs_sig); % collect noise-floor corrected summed harmonic magnitudes of freq/time domain approach

fprintf('\n') % display one empty rows

for j = 1:bs_sig % for all bootstrap draws ...
    fprintf('\n:: Signal - Draw %d out of %d',j,bs_sig)
    rng(j) % change random number generator for each bootstrap iteration --> keeps data reproducable
    pos_draws = randi(size(pos_index,2),1, bs_draws/2); % choose epochs with replacement
    neg_draws = randi(size(neg_index,2),1, bs_draws/2); % choose epochs with replacement
    tmp_epochs = epochs(:,[pos_index(pos_draws), neg_index(neg_draws)],:);  % draw epochs (same amount of positive and negative epochs)
    
    %% Spectral domain and phase locking value extraction
    
    tmp_wave = reshape(mean(tmp_epochs,2),chans,N); % average over epochs in time domain and remove epoch dimension
    
    % multiply time signals per channel and epoch with tukey window
    % --> first and last sample of time signal will be zero
    tmp_win_epochs = reshape(bsxfun(@times,reshape(tmp_epochs,chans*bs_draws,N),tukeywin(N,tukey_width)'),[chans,bs_draws,N]);
    tmp_fft = fft(tmp_win_epochs,NFFT,3); % do the fft for every epoch
    
    tmp_phase = angle(tmp_fft); % extract phase for all channels, epochs and frequencies
    tmp_phase = tmp_phase(:,:,1:ceil(NFFT/2)+1); % discard second half of spectrum
    tmp_phase = abs(sum(exp(1i*tmp_phase(:,1:bs_draws/2,:)),2) + sum(exp(1i*tmp_phase(:,bs_draws/2+1:end,:)),2))./bs_draws; % compute plv values
    
    tmp_plv_ptn = zeros(size(epochs,1),harms); % collect temporary PLV ptn values
    
    for h = 1:harms % for all harmonics ...
        plv_nf_index = [closesed_match(h)-5:closesed_match(h)-1,closesed_match(h)+1:closesed_match(h)+5]; % extract 5 bins around signal each side
        plv_nf_bins(:,j,h) = mean(tmp_phase(:,plv_nf_index),2); % average the bins per channel for a noise floor estimate
        tmp_plv_ptn(:,h) = tmp_phase(:,closesed_match(h)) - plv_nf_bins(:,j,h); % substract instentanious noise floor from PLV at harmonic frequencies
    end
    
    %%%%%%%%%%%%%% I droped this line, so negative values are not set to 0 anymore
     %plv_ptn(:,j,:) = tmp_plv_ptn .* (tmp_plv_ptn > 0); % set negative PLV PtN peaks to 0
    
    tmp_fft = reshape(mean(tmp_fft,2),chans,NFFT); % average the complex values
    
    tmp_spec = tmp_fft./N; % devide the spectrum by the NUMBER OF DATA POINTS per epoch
    tmp_spec = tmp_spec(:,1:ceil(NFFT/2)+1); % discard second half of spectrum
    
    tmp_signal = abs(tmp_fft); % extract real part of signal
    tmp_mean_phase = angle(tmp_fft); % extract phase of averaged spectrum
    
    for i = [1,2] % for noise floor uncorrected and corrected estimate ...
        % Substract the instentanious noise floor (5 bins above and below mf)
        if i == 2 % if corrected run ...
            for h = 1:harms % for all harmonics
                nf_index = [closesed_match(h)-5:closesed_match(h)-1,closesed_match(h)+1:closesed_match(h)+5]; % extract 5 bins around signal each side
                nf_bin = mean(tmp_signal(:,nf_index),2); % average the bins per channel
                % substract noise floor only at harmonic
                tmp_signal(:,closesed_match(h)) = tmp_signal(:,closesed_match(h)) - nf_bin;
                nf_bins(:,j,h) = 2*(nf_bin./N); % store the instentanious noise floor and convert to magnitude to make it comparable with spectral estimates
            end
        end
        
        % Re-compute complex array after noise floor correction
        tmp_z = tmp_signal.*exp(1i*tmp_mean_phase); % recomputed complex numbers
        
        tmp_mask = zeros(size(tmp_z)); % create mask with zeros matching tmp_z
        tmp_mask(:,closesed_match) = 1; % set the harmonic frequencies in all channels of the mask to one
        tmp_z = tmp_z.*tmp_mask; % set all values except harmonics to zero
        
        if i == 1 % if uncorrected run ...
            tmp_ftd_sum_harm_uncor = sum(2*abs(tmp_z/N),2); % sum up magnitudes of harmonics in spectrum
            tmp_ifft_uncor = ifft(tmp_z, NFFT, 2 ,'symmetric'); % Doing the ifft for each channel (row wise) treating tmp_z as conjugate symmetric
            tmp_ifft_uncor = tmp_ifft_uncor(:,1:N); % if the fft was zero padded reduce time signal to original time signal length
            
        elseif i == 2 % if corrected run ...
            %%%%% I dropped this line as it had no function
            %tmp_signal(:,closesed_match) = tmp_signal(:,closesed_match) .* (tmp_signal(:,closesed_match) > 0); % set negative PtN peaks to 0
            tmp_ftd_sum_harm = sum(2*abs(tmp_z/N),2); % sum up noise floor corrected magnitudes of harmonics in spectrum
            tmp_ifft = ifft(tmp_z, NFFT, 2 ,'symmetric'); % Doing the ifft for each channel (row wise) treating tmp_z as conjugate symmetric
            tmp_ifft = tmp_ifft(:,1:N); % if the fft was zero padded reduce time signal to original time signal length
        end
    end
    
    %% Time domain extraction approach based on averaging cycles
    
    tmp_wave_filt = (filtfilt(d,tmp_wave'))'; % filter the mean epoched data per channel to smoothen the signal
    
    tmp_cycle = zeros(size(epochs,1),number_of_cycles, ceil(samples_per_cycle)); % collect all individual cycles per epoch for each channel
    for cnt = 1:number_of_cycles % for the number of whole cycles per epoch ...
        % cut the epoch into cycle chunks
        tmp_cycle(:,cnt,:)= tmp_wave_filt(:,(1+ceil((cnt-1)*samples_per_cycle)):(ceil((cnt-1)*samples_per_cycle)+ceil(samples_per_cycle)));
    end
    
    %% Save the results of the bootstrap draws to the preallocated variables
    
    % Find the max and the min in the middle of the reconstructed time
    % domain signal for each channel (avoiding any possible irregularities
    % at the edges of the signal).
    % Amount of epoch to cut off each side is controlled by cut_off
    ptn_efr(:,j) = (max(tmp_ifft(:,floor(N*cut_off):ceil(N*(1-cut_off))),[],2) - min(tmp_ifft(:,floor(N*cut_off):ceil(N*(1-cut_off))),[],2))/2;
    % do the same for uncorrected estimate
    ptn_efr_uncor(:,j) = (max(tmp_ifft_uncor(:,floor(N*cut_off):ceil(N*(1-cut_off))),[],2) - min(tmp_ifft_uncor(:,floor(N*cut_off):ceil(N*(1-cut_off))),[],2))/2;
    ptn_td(:,j,:) = tmp_ifft; % reconstructed time domain signal (whole)
    ptn_td_uncor(:,j,:) = tmp_ifft_uncor; % uncorrected reconstructed time domain signal (whole)
    plv(:,j,:)= tmp_phase; % PLV
    plv_ptn(:,j,:) = tmp_plv_ptn; % noise floor corrected PLV
    spec(:,j,:) = tmp_spec; % frequency domain spectrum (whole signal) --> devided by number of data points in signal
    td(:,j,:) = tmp_wave; % collect the time domain signal (whole signal)
    cycles(:,j,:) = mean(tmp_cycle,2); % average filtered cycle in time domain
    cycle_amp(:,j) = (max(reshape(mean(tmp_cycle,2),chans,ceil(samples_per_cycle)),[],2)- min(reshape(mean(tmp_cycle,2),chans,ceil(samples_per_cycle)),[],2))/2; % (max-min)/2 amplitude in the time domain using time domain approach
    ftd_sum_harm_uncor(:,j) = tmp_ftd_sum_harm_uncor; % summed up uncorrected magnitudes of harmonics in spectrum
    ftd_sum_harm(:,j) = tmp_ftd_sum_harm; % summed up magnitudes of harmonics in spectrum
    
end % finish boostrap runs

%% Extract the average/std measures from boostrap runs

spec = abs(spec); % discard complex part of spectrum
nf_spec = abs(nf_spec); % discard complex part of noise floor spectrum
spec(:,:,2:end-1) = 2*spec(:,:,2:end-1); % extract magnitude (do not multiply DC and Nyquist frequency)
nf_spec(:,:,2:end-1) = 2*nf_spec(:,:,2:end-1); % extract magnitude (do not multiply DC and Nyquist frequency)

% Spectral approach:
EFR.peak_freq = freq(closesed_match); % save the frequencies of the extracted harmonic peaks
EFR.peak_index = closesed_match; % frequency indicies of the extracted harmonic peaks
EFR.wave_mean =  reshape(mean(td,2),chans,N); % Average bootstrapped time domain waveform and..
EFR.wave_std =  reshape(std(td,0,2),chans,N); % standard deviation
EFR.spec_mean = reshape(mean(spec,2),chans,ceil(NFFT/2)+1); % Average bootstrapped spectrum and ...
EFR.spec_std = reshape(std(spec,0,2),chans,ceil(NFFT/2)+1); %  standard deviation
EFR.nf_spec_mean = reshape(mean(nf_spec,2),chans,ceil(NFFT/2)+1); % Average bootstrapped noise floor spectrum and ...
EFR.nf_spec_std = reshape(std(nf_spec,0,2),chans,ceil(NFFT/2)+1); % standard deviation
EFR.ptn_spec = EFR.spec_mean - EFR.nf_spec_mean; % Average peak-to-noise floor spectrum
EFR.ptn_spec_harms = EFR.ptn_spec(:,closesed_match); % extract harmonics per channel
EFR.ptn_spec_sum = sum(EFR.ptn_spec_harms,2); % compute trigonometric sume of harmonics per channel

% Phase locking value:
EFR.plv_mean =  reshape(mean(plv,2),chans,ceil(NFFT/2)+1); % Average phase locking spectrum and ...
EFR.plv_std =  reshape(std(plv,0,2),chans,ceil(NFFT/2)+1); % ... standard deviation
EFR.plv_harms = EFR.plv_mean(:,closesed_match); % extract harmonics per channel
EFR.plv_nf_bins_mean = reshape(mean(plv_nf_bins,2),chans,harms); % Average instentanoius noise floor and ...
EFR.plv_nf_bins_std = reshape(std(plv_nf_bins,0,2),chans,harms); % ... standard deviation
EFR.ptn_plv_mean = reshape(mean(plv_ptn,2),chans,harms); % mean noise floor corrected PLV estimates of harmonics and ...
EFR.ptn_plv_std = reshape(std(plv_ptn,0,2),chans,harms); % ... standard deviation

% Frequency/time domain extraction approach
EFR.ptn_ftd_mean_ind = mean(ptn_efr,2); % Average amplitude peaks of reconstructed time domain signal (this seems to be strongly effected by outliers, use ptn_ftd_mean instead) and ...
EFR.ptn_ftd_std_ind = std(ptn_efr,0,2); % ... standard deviation
EFR.ptn_ftd_wave_mean = reshape(mean(ptn_td,2),chans,N); % Average reconstructed time domain waveform and
EFR.ptn_ftd_wave_std = reshape(std(ptn_td,0,2),chans,N); % ... standard deviation
EFR.ptn_ftd_mean = (max(EFR.ptn_ftd_wave_mean(:,floor(N*cut_off):ceil(N*(1-cut_off))),[],2) - min(EFR.ptn_ftd_wave_mean(:,floor(N*cut_off):ceil(N*(1-cut_off))),[],2))/2; % EFR amplitude extracted from the averaged bootstrapped waveforms


EFR.nf_bins_mean = reshape(mean(nf_bins,2),chans,harms); % Average instentanious noise floor and ...
EFR.nf_bins_std = reshape(std(nf_bins,0,2),chans,harms); % ... standard deviation
EFR.ptn_ftd_sum_harm_mean = reshape(mean(ftd_sum_harm,2),chans,1); % mean summed corrected harmonics before ifft of ftd approach and ...
EFR.ptn_ftd_sum_harm_std = reshape(std(ftd_sum_harm,0,2),chans,1); % ... standard deviation

% Uncorrected Frequency/time domain extraction approach
EFR.ftd_uncor_mean = mean(ptn_efr_uncor,2); % Average amplitude peaks of reconstructed time domain signal and ...
EFR.ftd_uncor_std = std(ptn_efr_uncor,0,2); % ... standard deviation
EFR.ftd_uncor_wave_mean =  reshape(mean(ptn_td_uncor,2),chans,N); % Average reconstructed time domain wave form ... and
EFR.ftd_uncor_wave_std =  reshape(std(ptn_td_uncor,0,2),chans,N); % ... standard deviation
EFR.ftd_sum_harm_uncor_mean = mean(ftd_sum_harm_uncor,2); % mean summed uncorrected harmonics before ifft of ftd approach and ...
EFR.ftd_sum_harm_uncor_std = std(ftd_sum_harm_uncor,0,2); % ... standard deviation

% Time domain extraction approach
EFR.td_amp_mean = mean(cycle_amp,2); % Average amplitude peaks of average cycle in time domain signal and ...
EFR.td_amp_std = std(cycle_amp,0,2); % ... standard deviation
EFR.td_cycle_mean = reshape(mean(cycles,2),chans,ceil(samples_per_cycle)); % Average cycle waveform and ...
EFR.td_cycle_std = reshape(std(cycles,0,2),chans,ceil(samples_per_cycle)); % ... standard deviation

% Genearl information:
EFR.time = time; % time domain vector
EFR.freq = freq(1:ceil(NFFT/2)+1); % frequency vector
EFR.nchannels = chans; % number of channels
EFR.channels = channels; % used channels
EFR.nepochs = eps; % number of epochs in input signal
EFR.nsamples = N; % time points in epoch
EFR.nfft = NFFT; % number of points used for fft
EFR.freq_res = freq_res; % frequency resolution
EFR.samples_per_cycle = samples_per_cycle; % exact number of samples per cycle (most likely not an integer)
EFR.ncycles = number_of_cycles; % number of whole cycles per epoch
EFR.fs = fs; % sampling rate
EFR.mf = mf; % modulation frequency
EFR.nharms = harms; % number of harmonics
EFR.bs_sig = bs_sig;  % Number of bootstrap iterations for the signal estimate
EFR.bs_noise = bs_noise; %  Number of bootstrap iterations for the nose-floor estimate
EFR.bs_draws = bs_draws; % Number of epochs use per bootstrap draw
EFR.cut_off = cut_off; % Amount of time signal to cut of at each side for peak extraction (0 to 1 = 0% to 100%)
EFR.tukey_width = tukey_width; % Width of the tukey window (0 to 1 = 0% to 100%)

EFR =  orderfields(EFR); % order structure alphabetically

%% Add explanation structure

EFR.info.bs_draws = 'Number of epochs use per bootstrap draw';
EFR.info.bs_noise = 'Number of bootstrap iterations for the noise-floor estimate';
EFR.info.bs_sig = 'Number of bootstrap iterations for the signal estimate';
EFR.info.channels = 'Used channels';
EFR.info.cut_off = 'Amount of time signal to cut off at each side for peak extraction (0 to 1 = 0% to 100%)';
EFR.info.freq = 'Frequency vector';
EFR.info.freq_res = 'Frequency resolution';
EFR.info.fs = 'Sampling rate';
EFR.info.ftd_sum_harm_uncor_mean = 'Mean summed uncorrected magnitude of harmonics before ifft of freq/time domain approach';
EFR.info.ftd_sum_harm_uncor_std = 'Standard deviation of ...mean variable';
EFR.info.ftd_uncor_mean = 'Average amplitude (max-min)/2 of uncorrected reconstructed time domain signal';
EFR.info.ftd_uncor_std = 'Standard deviation of ...mean variable';
EFR.info.ftd_uncor_wave_mean = 'Average uncorrected reconstructed time domain signal';
EFR.info.ftd_uncor_wave_std = 'Standard deviation of ...mean variable';
EFR.info.mf = 'Modulation frequency';
EFR.info.nchannels = 'Number of channels in data set';
EFR.info.ncycles = 'Number of whole mf cycles per epoch';
EFR.info.nepochs = 'Number of epochs in data set';
EFR.info.nf_bins_mean = 'Average instentanious noise floor for all harmonics';
EFR.info.nf_bins_std = 'Standard deviation of ...mean variable';
EFR.info.nf_spec_mean = 'Average bootstrapped magnitude noise floor spectrum';
EFR.info.nf_spec_std = 'Standard deviation of ...mean variable';
EFR.info.nfft = 'Number of samples used for fft';
EFR.info.nharms = 'Number of extracted harmonics';
EFR.info.nsamples = 'Number of samples per epoch in data set';
EFR.info.peak_freq = 'Frequencies of extracted harmonics';
EFR.info.peak_index = 'Frequeny indicies of extracted harmonics';
EFR.info.plv_harms = 'Average PLV estimate for all harmonics';
EFR.info.plv_mean = 'Average phase locking value spectrum';
EFR.info.plv_nf_bins_mean = 'Average instentanious PLV noise floor for all harmonics';
EFR.info.plv_nf_bins_std = 'Standard deviation of ...mean variable';
EFR.info.plv_std = 'Standard deviation of ...mean variable';
EFR.info.ptn_ftd_mean_ind = 'Average amplitude (max-min)/2 of noise-floor corrected reconstructed time domain signal';
EFR.info.ptn_ftd_std_ind = 'Standard deviation of ...mean variable';
EFR.info.ptn_ftd_mean = 'EFR amplitude extracted from the averaged bootstrapped waveforms';
EFR.info.ptn_ftd_sum_harm_mean = 'Mean summed noise-floor corrected magnitude of harmonics before ifft of freq/time domian approach';
EFR.info.ptn_ftd_sum_harm_std = 'Standard deviation of ...mean variable';
EFR.info.ptn_ftd_wave_mean = 'Average reconstructed time domain waveform';
EFR.info.ptn_ftd_wave_std = 'Standard deviation of ...mean variable';
EFR.info.ptn_plv_mean = 'Mean noise floor corrected PLV estimates of harmonics';
EFR.info.ptn_plv_std = 'Standard deviation of ...mean variable';
EFR.info.ptn_spec = 'Average peak-to-noise floor spectrum';
EFR.info.ptn_spec_harms = 'Average peak-to-noise estimates for all harmonics';
EFR.info.ptn_spec_sum = 'Summed up magnitudes of spectral PtN harmonics';
EFR.info.samples_per_cycle = 'Exact number of samples per cycle';
EFR.info.spec_mean = 'Average bootstrapped magnitude spectrum';
EFR.info.spec_std = 'Standard deviation of ...mean variable';
EFR.info.td_amp_mean = 'Average peak amplitude (max-min)/2 of average cycle in time domain signal';
EFR.info.td_amp_std = 'Standard deviation of ...mean variable';
EFR.info.td_cycle_mean = 'Average cycle waveform';
EFR.info.td_cycle_std = 'Standard deviation of ...mean variable';
EFR.info.time = 'Time vector';
EFR.info.tukey_width = 'Width of the tukey window (0 to 1 = 0% to 100%)';
EFR.info.wave_mean = 'Average bootstrapped time domain waveform';
EFR.info.wave_std = 'Standard deviation of ...mean variable';

%% Plot the results
if ~isempty(visual)
    
    ch = find(EFR.channels == visual); % find the right index of the channel
    figure;
    hold all
    plot(EFR.time, EFR.wave_mean(ch,:))
    plot(EFR.time, EFR.ptn_ftd_wave_mean(ch,:))
    plot(EFR.time, EFR.ftd_uncor_wave_mean(ch,:))
    plot(EFR.time(1:length(EFR.td_cycle_mean(ch,:))), EFR.td_cycle_mean(ch,:))
    hold off
    title(['Amplitudes - PtN: ',num2str(EFR.ptn_ftd_mean(ch)),', Uncorrected: ',num2str(EFR.ftd_uncor_mean(ch)), ', Cycle: ', num2str(EFR.td_amp_mean(ch))])
    legend('Time signal', 'PtN Frequency/time domain reconstruction', 'Frequency/time domain reconstruction','Mean cycle')
    xlabel('Time re: stimulus onset (sec)')
    ylabel('Amplitude in uV')
    axis tight
    
    %%
    
    figure;
    hold all
    plot(EFR.freq, EFR.spec_mean(ch,:))
    plot(EFR.freq, EFR.nf_spec_mean(ch,:))
    plot(EFR.freq, EFR.ptn_spec(ch,:))
    plot(EFR.peak_freq,EFR.ptn_spec_harms(ch,:),'*','linewidth',2)
    plot(EFR.peak_freq,EFR.nf_bins_mean(ch,:),'x','linewidth',2)
    legend('Spectral magnitudes', 'Noise floor', 'PtN magnitudes','PtN Peaks','Instantaneous noise floor')
    title(['Sum of harmonics - Spec-PtN: ',num2str(EFR.ptn_spec_sum(ch)),', FTD-PtN: ',num2str(EFR.ptn_ftd_sum_harm_mean(ch)),', Spec: ',num2str(sum(EFR.spec_mean(ch,EFR.peak_index))),', FTD: ',num2str(EFR.ftd_sum_harm_uncor_mean(ch))])
    xlabel('Frequency (Hz)')
    ylabel('Magnitude in uV')
    xlim([0,EFR.nharms*EFR.mf+100])
    
    %%
    
    figure;
    hold all
    plot(EFR.freq, EFR.plv_mean(ch,:))
    plot(EFR.peak_freq,EFR.plv_harms(ch,:),'o','linewidth',2)
    plot(EFR.peak_freq,EFR.ptn_plv_mean(ch,:),'*','linewidth',2)
    plot(EFR.peak_freq,EFR.plv_nf_bins_mean(ch,:),'x','linewidth',2)
    legend('Phase locking values','PLV Peaks', 'PtN PLVs','PLV noise floor')
    xlabel('Frequency (Hz)')
    ylabel('PLV (-)')
    title('Phase locking values')
    xlim([0,EFR.nharms*EFR.mf+100])
    
end

tEnd = toc(tStart);
fprintf('\n\n:: Function finished after %d minutes and %.0f seconds\n', floor(tEnd/60), rem(tEnd,60));

end

