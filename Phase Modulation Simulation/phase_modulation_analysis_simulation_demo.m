%% Create synthetic spiking data for a number of neurons
% We create a distribution of average spiking rates (0 to 20Hz) and we generate
% random number of spikes for each neuron.
% If Statistics Toolbox is available, it will calculate a poisson
% distribution for each neuron (random lambda) otherwise uniform.

N_neurons = 500;
n_spikes = 1e3;
t_max = 600; % 10 min of recording
SR = 1000; % sampling rate 1kHz

spike_rates = 20*rand(N_neurons,1);
spike_times = {};
for n=1:N_neurons;    
    if license('test','statistics_toolbox')
        tmp = cumsum(poissrnd(SR/spike_rates(n),1e5,1)); % poisson
        tmp = tmp(tmp<t_max*SR);
    else
        tmp = randi(t_max*SR,round(spike_rates(n)*t_max),1); % uniform
    end
    spike_times{n} = tmp;
end

%% Create random reference LFP signal 
% Cverlay multiple sinusoids with phase differences, on top of a carrier (respiratory) frequency.

min_f = 5;
max_f = 100;
n_f = 10;
freqs = max_f * rand(n_f,1) + min_f;
amplitudes = rand(n_f,1);
phase = pi*rand(n_f,1);

t = linspace(0,t_max,t_max*SR);

lfp = 2*sin(2*pi*4*t);

for f=1:length(freqs);
   lfp = lfp + amplitudes(f)*sin(2*pi*freqs(f)*t+phase(f));
end

%% Filter LFP in the 2-6Hz band
frequency_band = [2 6];
order = 6;
ripple = 20;
nyquist = SR/2;

[z, p, k] = butter(order,frequency_band/nyquist,'bandpass');
[sos_var,g] = zp2sos(z, p, k);
Hd = dfilt.df2sos(sos_var, g);
reset(Hd) % reset initial states
resp =  filtfilt(sos_var,g,lfp);

%% Extract phase in the 2-6Hz band
resp_ph = angle(hilbert(resp));
resp_ph_ctr = (randperm(length(resp_ph)));

%% Calculate phase modulation of each unit
phase_mod = [];
phase_mod_uniform_corrected = [];
phase_mod_ctr = [];

for n=1:length(spike_times);
    spike_phase = resp_ph(spike_times{n});
    if ~isempty(spike_phase);
        phase_mod(n) = phase_modulation(spike_phase,0);
        phase_mod_uniform_corrected(n) = phase_modulation(spike_phase,1);
        spike_phase = resp_ph_ctr(spike_times{n});
        phase_mod_ctr(n) = phase_modulation(spike_phase,0);
    end;
end;

%% Plot the spike firing for each neuron and LFP
figure('units','normalized','outerposition',[0 0 1 1]);

subplot(4,1,1);
plot(t,lfp,'k');
hold on;
plot(t,resp,'r');
xlim([10 11]);
set(gca,'XTick',[10:0.25:11])
xlabel('Time (s)');
ylabel('Amplitude (a.u.)');
title('LFP trace');
legend({'Raw LFP','Filtered LFP'},'box','off');

% LFP spectrogram
subplot(4,1,2)

window = 2*SR;
overlap = 0.8 * window;
[S,f_s,t_s,P] = spectrogram(lfp,window,overlap,512,SR,'power');
imagesc(t_s,f_s,20*log10(P)); axis xy; %colorbar
ylim([0 20]);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Spectral decomposition of LFP');

% Binned spike count color coded
subplot(4,1,3);
bins = 0:0.1:t_max;
binned_spikes = [];
for n=1:length(spike_times);
    binned_spikes(n,:) = histc(spike_times{n}/SR,bins);    
end

%imagesc(bins,1:length(spike_times),zscore(binned_spikes,[],2)); %plot normalized histogram
imagesc(bins,1:length(spike_times),binned_spikes);
title('Binned spike count');
xlabel('Time (s)');
ylabel('Neurons (#)');


% Plot the modulation statistics for the neuronal population
subplot(4,1,4);
[f,x]=ecdf(phase_mod,'function','cdf');
h1 = plot(x,100*(1-f),'k');
[f,x]=ecdf(phase_mod_uniform_corrected,'function','cdf');
hold on;
h2 = plot(x,100*(1-f),'r');

%[f,x]=ecdf(phase_mod_ctr,'function','cdf');
%h3 = plot(x,100*(1-f),'r--');
%legend({'','',''});
sig_levels = log(-log([0.05 0.01 0.001]));
verticalline(sig_levels(1),'k--')
%xlim([0 max(get(gca,'xlim'))])
title('Phase modulation simulation');
ylabel('Neurons (%)');
xlabel('Phase mod. (logZ)');
legend({'Raw','Uniform transformed phases','Significant modulation (p<0.05)'},'box','off','Location','southwest');