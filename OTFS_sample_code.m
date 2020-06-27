%
% Copyright (c) 2018, Raviteja Patchava, Yi Hong, and Emanuele Viterbo, Monash University
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution.
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
%ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
%WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
%ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
%(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
%ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%    - Latest version of this code may be downloaded from: https://ecse.monash.edu/staff/eviterbo/
%    - Freely distributed for educational and research purposes
%%

clc
clear all
close all
tic
%% OTFS parameters%%%%%%%%%%
% number of symbol
N = 8;
% number of subcarriers
M = 8;
% size of constellation
M_mod = 4;
M_bits = log2(M_mod);
% average energy per data symbol
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));  %为什么是这样子计算
% number of symbols per frame
N_syms_perfram = N*M;
% number of bits per frame
N_bits_perfram = N*M*M_bits;

SNR_dB = 20:2:20;
SNR = 10.^(SNR_dB/10);
noise_var_sqrt = sqrt(1./SNR);
sigma_2 = abs(eng_sqrt*noise_var_sqrt).^2;
%%
rng(1)
N_fram = 10^4;
err_ber = zeros(length(SNR_dB),1);
for iesn0 = 1:length(SNR_dB)
    for ifram = 1:N_fram
        %% random input bits generation%%%%%
        data_info_bit = randi([0,1],N_bits_perfram,1);  %128×1的0,1序列
        data_temp = bi2de(reshape(data_info_bit,N_syms_perfram,M_bits));%reshape将其变成64×2的矩阵，bi2de将其变成10进制的64×1的向量
        x1 = qammod(data_temp,M_mod,0,'gray');
        x = reshape(x1,N,M);
%         test1=ifft(x1);
%         test2=ifft(x);
%         X1= fft(ifft(x1).').'/sqrt(M/N);
%         X = fft(ifft(x).').'/sqrt(M/N); %%%ISFFT
%         s_mat = ifft(X.')*sqrt(M); % Heisenberg transform
%         s = s_mat(:);
        %% OTFS modulation%%%%
        s = OTFS_modulation(N,M,x);
        
        %% OTFS channel generation%%%%
        [taps,delay_taps,Doppler_taps,chan_coef] = OTFS_channel_gen(N,M);
        
        %% OTFS channel output%%%%%
        r = OTFS_channel_output(N,M,taps,delay_taps,Doppler_taps,chan_coef,sigma_2(iesn0),s);
        
        %% OTFS demodulation%%%%
        y = OTFS_demodulation(N,M,r);
        
        %% message passing detector%%%%
        x_est = OTFS_mp_detector(N,M,M_mod,taps,delay_taps,Doppler_taps,chan_coef,sigma_2(iesn0),y);
        
        %% output bits and errors count%%%%%
        data_demapping = qamdemod(x_est,M_mod,0,'gray');
        data_info_est = reshape(de2bi(data_demapping,M_bits),N_bits_perfram,1);
        errors = sum(xor(data_info_est,data_info_bit));
        err_ber(iesn0) = errors + err_ber(iesn0);
        %ifram
    end
end
err_ber_fram = err_ber/N_bits_perfram./N_fram;
semilogy(SNR_dB, err_ber_fram,'-*','LineWidth',2);
title(sprintf('OTFS'))
ylabel('BER'); xlabel('SNR in dB');grid on

toc