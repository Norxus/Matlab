function [X] = ourTest(Y,H,M)
% M为调制方式（进制），NT为发射端个数，NR为接收端个数，IterationTimes为迭代次数，
% Eb为设定的SNR范围（目前只支持5个不同的信噪比进行比较）
NR = size(H);
NT = size(H,2);
        
if M == 2
    H_hat = [real(H);imag(H)];
else
    H_hat = [real(H),-imag(H);imag(H),real(H)];
end


if M == 2
    GR = zeros(2*NR,NT);
else
    GR = zeros(2*NR,2*NT);
end
YR = [sign(real(Y)); sign(imag(Y))];
for kk = 1:2*NR
    GR(kk,:)= H_hat(kk,:)*YR(kk);
end
        
        if M == 2
               Constraints = [y<=GR*x_yalmip,y<=0];
        elseif M == 4
               Constraints = [y<=GR*x_yalmip/sqrt(2),y<=0];
        elseif M == 16
               Constraints = [y<=GR*((P*x_yalmip)/sqrt(10)),y<=0];
        end
        
        Objective = -sum(y);
        sol = optimize(Constraints, Objective, options);
        
        if M == 16
               solution = P*round(value(x_yalmip));
        else 
               solution = round(value(x_yalmip));
        end
        
        if M == 2
               X = solution;
        else 
               X = [solution(1:NT) + 1i*solution(NT+1:2*NT)];
        end      
end