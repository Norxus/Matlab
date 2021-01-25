clear all
close all
NT=4;      %发送天线数
NR=32;      %接受天线数
Eb = 10;

A = [1+i,1-i,-1+i,-1-i]/sqrt(2);
% A = [1;-1;-3;3]/sqrt(10);
[x4,x3,x2,x1] = ndgrid(A,A,A,A);
X_candidate = [x1(:) x2(:) x3(:) x4(:)].'; %所有排列组合情况
X_candidate = [real(X_candidate);imag(X_candidate)];
% [x4,x3,x2,x1] = ndgrid(A,A,A,A);
% X_candidate = [x1(:) x2(:) x3(:) x4(:)]';

ip = randn(2*NT,1)>0; %等概率产生0和1
% ip = randn(4*NT,1)>0;
X_bpsk = 2*ip-1; % 0 -> -1; 1 -> 1   

% P = zeros(2*NT,4*NT);
% P_1 = zeros(1,4*NT);
% P_1(1) = 2;
% P_1(2) = 1;
% P(1,:) = P_1;
% for cc = 2:2*NT
%     P_temp = circshift(P_1,2*(cc-1));
%     P(cc,:) = P_temp;
% end
    
% X_hat = P*X_bpsk;
X_hat = X_bpsk;
% X = [X_hat(1:NT)+1i*X_hat(NT+1:2*NT)]/sqrt(10);
X = [X_hat(1:NT)+1i*X_hat(NT+1:2*NT)]/sqrt(2);
% X = X_hat;
N = [randn(NR,1)+1i*randn(NR,1)]/sqrt(2); % 0均值高斯白噪声
H = [randn(NR,NT)+1i*randn(NR,NT)]/sqrt(2); % 瑞利衰落信道

Y = H*X+(10^(-Eb/20))*N;%利用上面的信号x，信道H，noise n计算出输出信号Y
Y_hat = [real(Y);imag(Y)];
H_hat = [real(H),-imag(H);imag(H),real(H)];
% H_hat = [real(H);imag(H)];
N_hat = [real(10^(-Eb/20)*N);imag(10^(-Eb/20)*N)];

GR = zeros(2*NR,2*NT);
% GR = zeros(2*NR,NT);
YR = [sign(real(Y)); sign(imag(Y))];
for kk = 1:2*NR
    GR(kk,:)= H_hat(kk,:)*YR(kk);
end
    
for i = 1:size(X_candidate,2)
   Z(:,i) = GR*X_candidate(:,i);
end

dist = sum(Z-abs(Z),1);

[minvalue indest] = max(dist);
% x = round(X_candidate(:,indest)*sqrt(10))
x = round(X_candidate(:,indest)*sqrt(2))
F = x;
   