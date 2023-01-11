%this code is written by Lakshmi Sudha Rani Nimmakayala with
%RollNo:22GS61R08
%wireless test bed on BLUE ESTIMATION and BLUE EQUALISATION FOR QPSK
% Date: 12-11-2022 
%Version: 2

%Transmitter
clc;
clear all;
n_p=100; %number of packets
symbols_packet=500; %number of symbols per packet
len=symbols_packet*2*n_p; %length of binary sequence to be generated
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

%QPSK modulation
for n =1:length(ip)/2
     two_bits = ip(2*n-1:2*n);
     if (two_bits(1) == 0 && two_bits(2) == 0)
           symbol_temp = -0.7071-1i*0.7071;
     elseif(two_bits(1) == 0 && two_bits(2) == 1)
           symbol_temp = -0.7071+1i*0.7071;
     elseif(two_bits(1) == 1 && two_bits(2) == 1)
           symbol_temp = 0.7071+1i*0.7071;        
     else
           symbol_temp = 0.7071-1i*0.7071;  
     end
           x(n) = symbol_temp;
end
len_sym=length(x);
Tx=[]; 
%creating transmit matrix for 100 packets
for i=1:n_p 
 Tx(i,:)=[ts x((((len_sym/n_p)*i)-(len_sym/n_p)+1):(len_sym/n_p)*i)]; 
end

SNRdB=0:15; 
for j=1:length(SNRdB) 
    SNR=10^(SNRdB(j)/10); 
    w=(1/sqrt(2*SNR))*(randn(n_p,len_ts+symbols_packet)+1j*randn(n_p,len_ts+symbols_packet));%complex noise matrix
    RX_AWGN=Tx+w; %received Signal in AWGN channel
    RX_RAY=h_Channel.*Tx+w; %received Signal in Rayleigh channel 
    Rx_out_blue=[]; 
    Rx_out_AWGN=[]; 
    Rx_out_RAY=[]; 
    Rx_out_noe=[];  
 
    %RECEIVER side

    for k=1:n_p
        %AWGN channel
        for l=1:symbols_packet 
            if (real(RX_AWGN(k,l+len_ts))>0 && imag(RX_AWGN(k,l+len_ts))>0)
                rx_bits=[1 1]; 
            elseif (real(RX_AWGN(k,l+len_ts))>0 && imag(RX_AWGN(k,l+len_ts))<0)
                rx_bits=[1 0]; 
            elseif (real(RX_AWGN(k,l+len_ts))<0 && imag(RX_AWGN(k,l+len_ts))>0)
                rx_bits=[0 1];
            else
                rx_bits=[0 0];
            end
            Rx_out_AWGN=[Rx_out_AWGN rx_bits];
        end
        %No channel estimation
        for l=1:symbols_packet 
            if (real(RX_RAY(k,l+len_ts))>0 && imag(RX_RAY(k,l+len_ts))>0)
                rx_bits=[1 1]; 
            elseif (real(RX_RAY(k,l+len_ts))>0 && imag(RX_RAY(k,l+len_ts))<0)
                rx_bits=[1 0]; 
            elseif (real(RX_RAY(k,l+len_ts))<0 && imag(RX_RAY(k,l+len_ts))>0)
                rx_bits=[0 1];
            else
                rx_bits=[0 0];
            end
            Rx_out_noe=[Rx_out_noe rx_bits];
        end
        %equalisation with perfect CSI
        true_est=RX_RAY(k,len_ts+1:end)/h_Channel(k);
        for l=1:symbols_packet 
            if (real(true_est(l))>0 && imag(true_est(l))>0)
                rx_bits=[1 1]; 
            elseif (real(true_est(l))>0 && imag(true_est(l))<0)
                rx_bits=[1 0]; 
            elseif (real(true_est(l))<0 && imag(true_est(l))>0)
                rx_bits=[0 1];
            else
                rx_bits=[0 0];
            end
            Rx_out_RAY=[Rx_out_RAY rx_bits];
        end
        %---------------------------------------------------------
         %Blue estimation
        h_blue(k)=sum(Tx(k,1:len_ts).*RX_RAY(k,1:len_ts))/sum(Tx(k,1:len_ts).^2);
        %blue equilisation using blue estimated channel coefficients
        x_est(1:symbols_packet)=(conj(h_blue(k)).*(RX_RAY(k,1+len_ts:end)))/(conj(h_blue(k)).*(h_blue(k)));
        for l=1:symbols_packet 
            if (real(x_est(l))>0 && imag(x_est(l))>0)
                rx_bits=[1 1]; 
            elseif (real(x_est(l))>0 && imag(x_est(l))<0)
                rx_bits=[1 0]; 
            elseif (real(x_est(l))<0 && imag(x_est(l))>0)
                rx_bits=[0 1];
            else
                rx_bits=[0 0];
            end
            Rx_out_blue=[Rx_out_blue rx_bits];
        end
        %---------------------------------------------------------
    end
    BER1(j)=sum(xor(ip,Rx_out_AWGN))/len; %BER for awgn
    BER_th(j)=0.5*erfc(sqrt(SNR/2)); %Theoretical BER for AWGN
    BER2(j)=sum(xor(ip,Rx_out_noe))/len; %BER for no channel estimation
    BER3(j)=sum(xor(ip,Rx_out_RAY))/len; %BER for perfect CSI
    BER4(j)=sum(xor(ip,Rx_out_blue))/len; %BER for blue estimation,blue equalisation
end
% graphs for BER vs SNR
figure; 
semilogy(SNRdB,BER1,'r-*'); 
hold on; 
semilogy(SNRdB,BER_th,'c'); 
grid on;
semilogy(SNRdB,BER2,'k->'); 
semilogy(SNRdB,BER3,'b'); 
semilogy(SNRdB,BER4,'m+'); 
legend('AWGN simulation results','AWGN theoretical','No channel estimation','equalisation using perfect CSI','BLUE estimation&equalisation'); 
title ('BER vs SNR(dB) for QPSK'); 
xlabel('SNR(dB)'); 
ylabel('BER'); 
axis([0 15 10^-5 1]); 

  