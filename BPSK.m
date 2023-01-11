%this code is written by Lakshmi Sudha Rani Nimmakayala with
%RollNo:22GS61R08
%wireless test bed on BLUE ESTIMATION and BLUE EQUALISATION FOR BPSK
% Date: 12-11-2022 
%Version: 1

%Transmitter
clc;
clear all;
n_p=100; %number of packets
bits_packet=500; %number of bits per packet
len=bits_packet*n_p; %length of bits to be generated
LTS=[ 1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 0, 1, -1, -1, 1,1, -1, 1, -1, 1, -1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, 1, 1, 1]; 
h=readmatrix("channel.csv");
h_Channel=h(:,1);%rayleigh channel coefficients(1st column)
ts=[]; 
%taking the Long training sequence 5-times
for j = 1:5 
 ts = [ts LTS]; 
end
len_ts=length(ts); 
ip = randn(1,len)<0;%binary sequence 

%BPSK modulation
x = 2*ip -1; %binary sequence mapping to +1 and -1
Tx=[]; 
%creating transmit matrix for 100 packets
for i=1:n_p 
 Tx(i,:)=[ts x((((len/n_p)*i)-(len/n_p)+1):(len/n_p)*i)]; 
end

SNRdB=0:15; 
for j=1:length(SNRdB) 
    SNR=10^(SNRdB(j)/10); 
    w=(1/sqrt(2*SNR))*(randn(n_p,len_ts+bits_packet)+1j*randn(n_p,len_ts+bits_packet));%complex noise matrix
    RX_AWGN=Tx+w; %received Signal in AWGN channel
    RX_RAY=h_Channel.*Tx+w; %received Signal in Rayleigh channel 
    Rx_out_blue=[]; 
    Rx_out_AWGN=[]; 
    Rx_out_RAY=[]; 
    Rx_out_noe=[]; 
    temp=zeros(1,bits_packet); 
 
    %RECEIVER side

    for k=1:n_p
        %AWGN channel
        for l=1:bits_packet 
            if real(RX_AWGN(k,l+len_ts))>0 
                temp(l)=1; 
            else 
                temp(l)=0; 
            end
        end
        Rx_out_AWGN=[Rx_out_AWGN temp]; 
        % detection of received bits without Channel estimation 
        for l=1:bits_packet 
            if real(RX_RAY(k,l+len_ts))>0 
                temp(l)=1; 
            else 
                temp(l)=0; 
            end
        end
        Rx_out_noe=[Rx_out_noe temp]; 
        % dection of received bits with perfect CSI
        for l=1:bits_packet
            if real(RX_RAY(k,l+len_ts)/h_Channel(k))>0 
                temp(l)=1; 
            else 
                temp(l)=0; 
            end
        end
        Rx_out_RAY=[Rx_out_RAY temp]; 
        %Blue estimation
        h_blue(k)=sum(Tx(k,1:len_ts).*RX_RAY(k,1:len_ts))/sum(Tx(k,1:len_ts).^2);
        %blue equilisation using blue estimated channel coefficients
        x_est(1:bits_packet)=(conj(h_blue(k)).*(RX_RAY(k,1+len_ts:end)))/(conj(h_blue(k)).*(h_blue(k)));
        for l=1:bits_packet
            if real(x_est(l))>0 
                temp(l)=1; 
            else 
                temp(l)=0; 
            end
        end
        Rx_out_blue=[Rx_out_blue temp];
    end
    BER1(j)=sum(xor(ip,Rx_out_AWGN))/len; %BER for awgn
    BER2(j)=sum(xor(ip,Rx_out_noe))/len; %BER for no channel estimation
    BER3(j)=sum(xor(ip,Rx_out_RAY))/len; %BER for perfect CSI
    BER4(j)=sum(xor(ip,Rx_out_blue))/len; %BER for blue estimation,blue equalisation
    BER_th(j)=0.5*erfc(sqrt(SNR)); %Theoretical BER for AWGN
    BER_rth(j)=0.5*(1-sqrt(SNR/(1+SNR))); %Theoretical BER for rayleigh 
    mse(j)=mean(abs(h_blue-h_Channel.').^2); %mean square error for h and h_est
end

% graphs for BER vs SNR
figure; 
semilogy(SNRdB,BER1,'r-*'); 
hold on; 
semilogy(SNRdB,BER_th,'c'); 
grid on; 
semilogy(SNRdB,BER_rth,'g'); 
semilogy(SNRdB,BER2,'k->'); 
semilogy(SNRdB,BER3,'b'); 
semilogy(SNRdB,BER4,'m+'); 
legend('AWGN simulation results','AWGN theoretical','rayleigh theoretical','No channel estimation','equalisation with perfect CSI','BLUE estimation&equalisation'); 
title ('BER vs SNR(dB) for BPSK'); 
xlabel('SNR(dB)'); 
ylabel('BER'); 
axis([0 15 10^-5 1]); 

figure; 
semilogy(SNRdB,mse,'k'); 
legend('mse'); 
grid on; 
title ('mse vs SNR(dB)'); 
xlabel('SNR(dB)'); 
ylabel('mse'); 