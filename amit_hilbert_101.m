%%
clear
clc
close all

%% Load data
filename = "pini_12_sec_a.wav";
[record,fs] = audioread(filename);
record= 0.99*record/max(abs(record));

% Crop 5 seconds
% T = 5;
% record(round(T*fs):end)=[];

%% Create demo signal
dt=1/fs;
freq1 = 81;
freq2 = 89;
freq3 = 100;
T = 5;
t = 0:dt:T;

a=1;
b=1;
c=1;
demo_data = a*sin(2*pi*freq1*t)+b*sin(2*pi*freq2*t)+c*sin(2*pi*freq3*t);

% Add noise
% demo_data = demo_data + 0.5*randn(size(t));

demo_data = 0.99*demo_data /max(abs(demo_data));

%% Choose data to work on, real or demo
x = demo_data;  % demo_data, record
figure
subplot(2,1,1)
plot(demo_data)
subplot(2,1,2)
plot(record)

%% delete part of the signal
% Set
n  = 1024;
IX = round(length(x)/2) + [1:n];
x(IX)=[];

dt = 1/fs;
t  = [1:length(x)].*dt;
plot(t,x)

%% Simple example with a single band-pass
fpass = [85 105];
y = bandpass(x,fpass,fs);
ph = angle(hilbert(y));
plot(diff(ph),'.')
ylim([-1 1].*0.3)

%% More Complex example with a filter bank
bw  = 5;
hop = bw+1;
filter_bank = [78:hop:100]'+[0 bw];
% filter_bank = [75 85; 85 95; 95 105];
% filter_bank = [79 81; 89 91; 99 101];
% filter_bank = [99 101];

nBands = size(filter_bank,1);

Y = zeros(length(x),nBands);
t2 = 0:dt:1;
estimated_ph_diff = zeros(nBands,1);
for ii_band = 1:nBands
    fpass = filter_bank(ii_band,:);
    f1 = fpass(1);
    f2 = fpass(2);
    filterOrder = 1000;
    cutoffFreq = [f1/(fs/2), f2/(fs/2)]; % Normalize cutoff frequencies
    b = fir1(filterOrder, cutoffFreq, 'bandpass');
%     Y(:,ii_band) = bandpass(x,fpass,fs);
    Y(:,ii_band) = filtfilt(b,1,x);
    sdf = diff(angle(hilbert(sin(2*pi*t2*mean(fpass)))));
    estimated_ph_diff(ii_band) = median(sdf);
end

%% Analyze the results
% get the diff of the phase
ph = angle(hilbert(Y));
ph_diff = diff(ph);

% center the results, and take abs
z = abs(ph_diff - median(ph_diff)); 

% Remove artifacts: high values, and start+end of the signal
z(z>pi)=nan;
z(1:3*filterOrder)=nan;
z(end-3*filterOrder:end)=nan;

% Remove noisy pass' that cross the threshold too much
ph_diff_factor = 0.8;  % low number allow more crossing. will be filtered out later
w = sum(z > ph_diff_factor*estimated_ph_diff') % for each band, the number of crossing the threshold
% w(w==0)=1;

% #crossing the threshold allowed. Low number filter out more data
cross_factor=10;
z(:,w>cross_factor)=nan;

% z = z./w;
plot(z,'.')
xline(IX(1))
% legend(filter_bank)

%%

%%

%%
fs = 1000;
dt = 1/fs;
f = 100;
T = 10;
t = 0:dt:T;
x = sin(2*pi*f*t);
% x = x+randn(size(t));
IX = 1000+[100:117];
% IX = [];
t(IX)=[];
x(IX)=[];
plot(x)

%%
ph = angle(hilbert(x));
plot(diff(ph),'.')
