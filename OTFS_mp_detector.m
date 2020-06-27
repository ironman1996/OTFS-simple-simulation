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
function x_est = OTFS_mp_detector(N,M,M_mod,taps,delay_taps,Doppler_taps,chan_coef,sigma_2,y)

yv = reshape(y,N*M,1);
n_ite = 200;
delta_fra = 0.6;
alphabet = qammod(0:M_mod-1,M_mod,0,'gray');

mean_int = zeros(N*M,taps);
var_int = zeros(N*M,taps);
p_map = ones(N*M,taps,M_mod)*(1/M_mod);

conv_rate_prev = -0.1;
for ite=1:n_ite
    %% Update mean and var
    for ele1=1:1:M
        for ele2=1:1:N
            mean_int_hat = zeros(taps,1);
            var_int_hat = zeros(taps,1);
            for tap_no=1:taps
                m = ele1-1-delay_taps(tap_no)+1;
                add_term = exp(1i*2*(pi/M)*(m-1)*(Doppler_taps(tap_no)/N));
                add_term1 = 1;
                if ele1-1<delay_taps(tap_no)
                    n = mod(ele2-1-Doppler_taps(tap_no),N) + 1;
                    add_term1 = exp(-1i*2*pi*((n-1)/N));
                end
                new_chan = add_term * (add_term1) * chan_coef(tap_no);
                
                for i2=1:1:M_mod
                    mean_int_hat(tap_no) = mean_int_hat(tap_no) + p_map(N*(ele1-1)+ele2,tap_no,i2) * alphabet(i2);
                    var_int_hat(tap_no) = var_int_hat(tap_no) + p_map(N*(ele1-1)+ele2,tap_no,i2) * abs(alphabet(i2))^2;
                end
                mean_int_hat(tap_no) = mean_int_hat(tap_no) * new_chan;
                var_int_hat(tap_no) = var_int_hat(tap_no) * abs(new_chan)^2;
                var_int_hat(tap_no) = var_int_hat(tap_no) - abs(mean_int_hat(tap_no))^2;
            end
            
            mean_int_sum = sum(mean_int_hat);
            var_int_sum = sum(var_int_hat)+(sigma_2);
            
            for tap_no=1:taps
                mean_int(N*(ele1-1)+ele2,tap_no) = mean_int_sum - mean_int_hat(tap_no);
                var_int(N*(ele1-1)+ele2,tap_no) = var_int_sum - var_int_hat(tap_no);
            end
            
        end
    end
    %% Update probabilities
    sum_prob_comp = zeros(N*M,M_mod);
    dum_eff_ele1 = zeros(taps,1);
    dum_eff_ele2 = zeros(taps,1);
    for ele1=1:1:M
        for ele2=1:1:N
            dum_sum_prob = zeros(M_mod,1);
            log_te_var = zeros(taps,M_mod);
            for tap_no=1:taps
                
                if ele1+delay_taps(tap_no)<=M
                    eff_ele1 = ele1 + delay_taps(tap_no);
                    add_term = exp(1i*2*(pi/M)*(ele1-1)*(Doppler_taps(tap_no)/N));
                    int_flag = 0;
                else
                    eff_ele1 = ele1 + delay_taps(tap_no)- M;
                    add_term = exp(1i*2*(pi/M)*(ele1-1-M)*(Doppler_taps(tap_no)/N));
                    int_flag = 1;
                end
                add_term1 = 1;
                if int_flag==1
                    add_term1 = exp(-1i*2*pi*((ele2-1)/N));
                end
                eff_ele2 = mod(ele2-1+Doppler_taps(tap_no),N) + 1;
                new_chan = add_term * add_term1 * chan_coef(tap_no);
                
                dum_eff_ele1(tap_no) = eff_ele1;
                dum_eff_ele2(tap_no) = eff_ele2;
                for i2=1:1:M_mod
                    dum_sum_prob(i2) = abs(yv(N*(eff_ele1-1)+eff_ele2)- mean_int(N*(eff_ele1-1)+eff_ele2,tap_no) - new_chan * alphabet(i2))^2;
                    dum_sum_prob(i2)= -(dum_sum_prob(i2)/var_int(N*(eff_ele1-1)+eff_ele2,tap_no));
                end
                dum_sum = dum_sum_prob - max(dum_sum_prob);
                dum1 = sum(exp(dum_sum));
                log_te_var(tap_no,:) = dum_sum - log(dum1);
            end
            for i2=1:1:M_mod
                ln_qi(i2) = sum(log_te_var(:,i2));
            end
            dum_sum = exp(ln_qi - max(ln_qi));
            dum1 = sum(dum_sum);
            sum_prob_comp(N*(ele1-1)+ele2,:) = dum_sum/dum1;
            for tap_no=1:1:taps
                eff_ele1 = dum_eff_ele1(tap_no);
                eff_ele2 = dum_eff_ele2(tap_no);
                
                dum_sum = log_te_var(tap_no,:);
                ln_qi_loc = ln_qi - dum_sum;
                dum_sum = exp(ln_qi_loc - max(ln_qi_loc));
                dum1 = sum(dum_sum);
                p_map(N*(eff_ele1-1)+eff_ele2,tap_no,:) = (dum_sum/dum1)*delta_fra + (1-delta_fra)*reshape(p_map(N*(eff_ele1-1)+eff_ele2,tap_no,:),1,M_mod);
            end
            
        end
    end
    conv_rate =  sum(max(sum_prob_comp,[],2)>0.99)/(N*M);
    if conv_rate==1
        sum_prob_fin = sum_prob_comp;
        break;
    elseif conv_rate > conv_rate_prev
        conv_rate_prev = conv_rate;
        sum_prob_fin = sum_prob_comp;
    elseif (conv_rate < conv_rate_prev - 0.2) && conv_rate_prev > 0.95
        break;
    end
end
x_est = zeros(N,M);
for ele1=1:1:M
    for ele2=1:1:N
        [~,pos] = max(sum_prob_fin(N*(ele1-1)+ele2,:));
        x_est(ele2,ele1) = alphabet(pos);
    end
end
end