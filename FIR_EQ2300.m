clear
close all
clc

%% Task 2: low pass filter design
%%%% LP filter
fc = 1/16;  %cutoff freq
fs = 1/8;   %supressed freq
N = 51;     % filter length
N_pad = N + 1e4;    %pad zeros afterwards

df = 1/N;
df_pad = 1/N_pad;

f = (0:df:1-df)';
f_pad = (0:df_pad:1-df_pad)';


%window selection
% w = window(@rectwin, N);
w = window(@hamming, N);
% w = window(@hann, N);
% w = window(@triang, N);
% w = window(@chebwin, N);
% w = window(@kaiser, N);


% get ideal filter frequency response

H = ones(N,1);
H(abs(f)>fc & abs(1-f)>fc)=0;

h = real(ifft(H));
h = fftshift(h);   

% get filter impulse response
hd = h .* w;
hd_pad = cat(1, hd, zeros(N_pad-N,1));  %padding zeros

% hd = hd / max(hd);
% n = (0:N-1)';
% 
% OUTPUT = [n,hd]';;
% fID = fopen('lowpass.txt', 'w');
% fprintf(fID, '%.6f %.12f\n', OUTPUT);
% fclose(fID);

figure()
fig = stem(0:N-1, hd, 'k', 'linewidth', 1.5);
pbaspect([3 1 1])

% find frequency response
Hd = abs(fft(hd));
Hd_pad = abs(fft(hd_pad));

% normalise frequency response
Hd = Hd / max(Hd);
Hd_pad = Hd_pad / max(Hd_pad);
figure()
hold on
plot(f_pad(f_pad<0.5), Hd_pad(f_pad<0.5), 'linewidth', 1.5)
ylim([0 1.2])
plot([fc fc],[0 1.2], '-.k', 'linewidth', 1.5 )
grid on; grid minor; box on


% frequency response in dB scale
Hd_dB = mag2db(Hd);
Hd_pad_dB = mag2db(Hd_pad);


% clip only half frequency
f_half = f(f<0.5);
f_pad_half = f_pad(f_pad<0.5);
Hd_dB_half = Hd_dB(f<0.5);
Hd_pad_dB_half = Hd_pad_dB(f_pad<0.5);
figure();
hold on
p4=plot(f_pad_half, Hd_pad_dB_half, 'linewidth',2.5);
p5=plot(f_half, Hd_dB_half, '.k', 'markersize', 30);
ylim([-120,20])
xax = get(gca,'xlim');
yax = get(gca,'ylim');
p1=plot([fc fc], yax, '-.k', 'linewidth', 2);
p2=plot([fs fs], [-120 -40], '-.r', 'linewidth', 2);
p3=plot([fs 0.5],[-40 -40], '-.g', 'linewidth', 2);
set(gca, 'fontsize', 20)
xlabel('Normalised frequency \nu') % cycles/sample
ylabel('Magnitude, [dB]')
% title('Frequency respo')
legend([p4,p5,p1, p2,p3],{'frequncy response with padding zeros', ...
    'frequency response using 51 data points', 'cutoff frequency', ...
    'stop band edge frequency', 'supression level'})
grid on; grid minor; box on;



%% Task 3: high pass filter design

m = (N-1)/2;
delta = zeros(N,1);
delta(m+1) = 1;

hH = delta - hd;
hH_pad = cat(1, hH, zeros(N_pad-N,1));

figure()
stem(0:N-1, hH, 'k', 'linewidth', 1.5);
pbaspect([2 1 1])

% n = (0:N-1)';
% hH = hH/max(hH);
% 
% OUTPUT = [n,hH]';
% fID = fopen('highpass.txt', 'w');
% fprintf(fID, '%.6f %.12f\n', OUTPUT);
% fclose(fID);

% frequency response
Hh = abs(fft(hH));
Hh_pad = abs(fft(hH_pad));

%normalise the frequency response
Hh_pad = Hh_pad / max(Hh_pad);
Hh = Hh/max(Hh);

% magnitude in dB scale
Hh_pad_dB = mag2db(Hh_pad(f_pad<0.5));
Hh_dB = mag2db(Hh(f<0.5));
figure();
hold on
plot(f_pad(f_pad<0.5), Hh_pad_dB,'r','linewidth',2.5)
plot(f(f<0.5), Hh_dB, '.k','markersize',30);
set(gca, 'fontsize', 20)
xlabel('Normalised frequency \nu') % cycles/sample
ylabel('Magnitude, [dB]')
ylim([-120 10])
plot([fs fs], [-120 10], '-.g', 'linewidth', 2);
grid on; grid minor; box on;
legend('Frequency response with padding zeros', ...
    'Frequency response w/o padding zeros'...
    ,'Stop band edge frequency')
% p1=plot([fc fc], yax, '-.k', 'linewidth', 2);
% p2=plot([fs fs], [-120 -40], '-.r', 'linewidth', 2);


%% task 4, quantization

F = 10;

% quantized impulsee response
hd_Q = round(hd*2^F);
hd_Q = hd_Q*2^-F;
hd_pad_Q = cat(1, hd_Q, zeros(N_pad-N,1));

%frequency response
Hd_Q = abs(fft(hd_Q));
Hd_pad_Q = abs(fft(hd_pad_Q));

%normalise frequency response
Hd_Q = Hd_Q / max(Hd_Q);
Hd_pad_Q = Hd_pad_Q / max(Hd_pad_Q);

% frequency resonsee in dB scale
Hd_Q_dB = mag2db(Hd_Q);
Hd_pad_Q_dB = mag2db(Hd_pad_Q);

%clip the half part of the frequency response
Hd_Q_dB = Hd_Q_dB(f<0.5);
Hd_pad_Q_dB = Hd_pad_Q_dB(f_pad<0.5);

figure()
hold on
p4=plot(f_pad(f_pad<0.5), Hd_pad_Q_dB,'b','linewidth',2.5);
p5=plot(f(f<0.5),Hd_Q_dB,'.k','markersize',30);
ylim([-120,20])
xax = get(gca,'xlim');
yax = get(gca,'ylim');
p1=plot([fc fc], yax, '-.k', 'linewidth', 2);
p2=plot([fs fs], [-120 -40], '-.r', 'linewidth', 2);
p3=plot([fs 0.5],[-40 -40], '-.g', 'linewidth', 2);
set(gca, 'fontsize', 20)
xlabel('Normalised frequency \nu') % cycles/sample
ylabel('Magnitude, [dB]')
% title('Frequency respo')
legend([p4,p5,p1, p2,p3],{'frequncy response with padding zeros', ...
    'frequency response using 51 data points', 'cutoff frequency', ...
    'stop band edge frequency', 'supression level'})
grid on; grid minor; box on;

%% task 5 SQNR

E_hw = sum(hd.^2);
E_he = sum((hd - hd_Q).^2);

SQNR = 10*log10(E_hw / E_he)

% end



