clc;
clear;
close all;

% Init
nb=100;  % Number of bits
nsb=100;  % Number of samples per bit
ns=nsb*nb;  % Number of samples
Tb=100;  % Bit transmission time
fc1=1;  % Carrier 1 frequency
fc2=2;  % Carrier 2 frequency
t=linspace(0,Tb*nb,ns);  % Discrete time sequence [0; Tb*nb] (ns samples)
snr=1; % Signal-to-noise ratio in dB

% Generate carrier signal
c1=sqrt(2/Tb)*sin(2*pi*fc1*t);
c2=sqrt(2/Tb)*sin(2*pi*fc2*t);

% Generate binary data bits
b=randi([0,1],1,nb);

% Generate data signal
d=repelem(b,nsb);  % Repeat each bit nsb times
d_inv=1-d;  % Inverse d

% Perform FSK modulation
fsk1=c1.*d;
fsk2=c2.*d_inv;
fsk=fsk1+fsk2;

% Add white noise to signal
fsk=awgn(fsk,snr,'measured');

% Perform FSK demodulation
for i=1:nb
    fsk_i=fsk((i-1)*nsb+1:i*nsb);
    c1_i=c1((i-1)*nsb+1:i*nsb);
    c2_i=c2((i-1)*nsb+1:i*nsb);

    % Correlator
    I1=sum(fsk_i.*c1_i);
    I2=sum(fsk_i.*c2_i);
    E1=sum(c1_i.^2);
    E2=sum(c2_i.^2);

    % Decision device
    if I1-E1/2>I2-E2/2
        demod(i)=1;
    else
        demod(i)=0;
    end
end

% Plot binary data bits
subplot(6,1,1);stem(b);
xticks(1:nsb);title('Binary data bits');ylabel('b(n)');grid on;

% Plot binary data signal
subplot(6,1,2);plot(t,d);
xticks(1:nsb);title('Data signal');ylabel('d(t)');grid on;

% Plot carrier 1 signal
subplot(6,1,3);plot(t,c1);
xticks(1:nsb);title('Carrier 1 signal');ylabel('c1(t)');grid on;

% Plot carrier 2 signal
subplot(6,1,4);plot(t,c2);
xticks(1:nsb);title('Carrier 2 signal');ylabel('c2(t)');grid on;

% Plot FSK signal
subplot(6,1,5);plot(t, fsk);
xticks(1:nsb);title('FSK signal');ylabel('fsk(t)');grid on;

% Plot demodulated
subplot(6,1,6);stem(demod);
xticks(1:nsb);title('FSK demodulated');ylabel('b(n)');grid on;

% Calculate bit error rate
[numBitErrs, BER] = biterr(b, demod);
fprintf("Number of Bit Errors: %d\nBit Error Rate: %f\n", numBitErrs, BER * 100);
