clc;
clear;
close all;

% Init
nb=100;  % Number of bits
nsb=100;  % Number of samples per bit
ns=nsb*nb;  % Number of samples
Tb=100;  % Bit transmission time (s)
Eb=(1/2)*Tb;  % Engery per pit
Rb=1/Tb;  % Bit rate (bps)
Ac1=sqrt(2*Eb/Tb);  % Carrier 1 amplitude
fc1=2/Tb;  % Carrier 1 frequency (Hz)
pc1=0;  % Carrier 1 phase (rad)
Ac2=-Ac1;
fc2=fc1;
pc2=0;
t=linspace(0,Tb*nb,ns);  % Discrete time sequence [0; Tb*nb] (ns samples)
snr_db=1;  % Signal-to-noise ratio (dB)
snr=10^(snr_db/10)
N0=2*Eb/snr;  % Variance = N0/2

% Generate carrier signal
c1=Ac1*sin(2*pi*fc1*t+pc1);
c2=Ac2*sin(2*pi*fc2*t+pc2);

% Generate binary data bits
b=randi([0,1],1,nb);

% Generate data signal
d=repelem(b,nsb);  % Repeat each bit nsb times
d_inv=1-d;  % Inverse d

% Perform ASK modulation
ask1=c1.*d;
ask2=c2.*d_inv;
ask=ask1+ask2;

% Add white noise to signal
% ask=awgn(ask,snr_db,'measured');
n = sqrt(N0/2)*randn([1,length(ask)]);
ask=ask+n;

% Perform ASK demodulation
for i=1:nb
    ask_i=ask((i-1)*nsb+1:i*nsb);
    c1_i=c1((i-1)*nsb+1:i*nsb);
    c2_i=c2((i-1)*nsb+1:i*nsb);

    % Correlator
    I1=sum(ask_i.*c1_i);
    I2=sum(ask_i.*c2_i);
    E1=sum(c1_i.^2);
    E2=sum(c2_i.^2);

    % Decision device
    if I1-E1/2>I2-E2/2
        demod(i)=1;
    else
        demod(i)=0;
    end
end

% Init figs
fig1 = figure;
fig2 = figure;
figure(fig1);

% Plot binary data bits
subplot(6,1,1);stem(b);
title('Binary data bits');ylabel('b(n)');grid on;

% Plot binary data signal
subplot(6,1,2);plot(t,d);
title('Data signal');ylabel('d(t)');grid on;

% Plot carrier 1 signal
subplot(6,1,3);plot(t,c1);
title('Carrier 1 signal');ylabel('c(t)');grid on;

% Plot carrier 2 signal
subplot(6,1,4);plot(t,c2);
title('Carrier 2 signal');ylabel('c(t)');grid on;

% Plot ASK signal
subplot(6,1,5);plot(t, ask);
title('ASK signal');ylabel('ask(t)');grid on;

% Plot demodulated
subplot(6,1,6);stem(demod);
title('ASK demodulated');ylabel('b(n)');grid on;

% Calculate bit error rate
[numBitErrs, BER] = biterr(b, demod);
fprintf("Number of Bit Errors: %d\nBit Error Rate: %f%%\n", numBitErrs, BER*100);

% Calculate theoretical bit error rate
fprintf("Bit Error Rate (Theory): %f%%\n", 0.5*erfc(sqrt(Eb/N0))*100);

% Plot BER
figure(fig2);
x_db = -5:0.5:15;
y = 0.5*erfc(sqrt(10.^(x_db/10)));
semilogy(x_db,y);
title('BER');xlabel('SNR (dB)');ylabel('BER');grid on; hold on;
plot(snr_db,BER,'rx');
