% function [SER] = our(M,NT,NR,IterationTimes,Eb)
Eb = -10 : 5 : 10;
M = 16;
NT = 4;
NR = 32;
IterationTimes = 1000;
E_f10 = [];
E_f5 = [];
E_0 = [];
E_5 = []; 
E_10 = [];
errorNumber =0;
errorNumber_5 =0;
errorNumber_10 =0;
errorNumber_f5 =0;
errorNumber_f10 =0;
l=0;
l1=0;
l2=0;
l3=0;
l4=0;

y = sdpvar(2*NR,1)
options= sdpsettings;
options.solver = 'intlinprog';

for jj=1:IterationTimes
    if M == 2
        ip = randn(NT,1)>0; %等概率产生0和1
        X_bpsk = 2*ip-1; % 0 -> -1; 1 -> 1
        x_yalmip = 2*binvar(NT,1)-1;
        X = X_bpsk;
        X_transmit = X_bpsk;
    elseif M == 4
             ip = randn(2*NT,1)>0; %等概率产生0和1
             X_bpsk = 2*ip-1; % 0 -> -1; 1 -> 1
             X_hat = X_bpsk;
             x_yalmip = 2*binvar(2*NT,1)-1;
             X = [X_hat(1:NT)+1i*X_hat(NT+1:2*NT)]/sqrt(2);
             X_transmit = [X_hat(1:NT)+1i*X_hat(NT+1:2*NT)];
    elseif M == 16
            ip = randn(4*NT,1)>0; %等概率产生0和1
            X_bpsk = 2*ip-1; % 0 -> -1; 1 -> 1
            P = zeros(2*NT,4*NT);
            P_1 = zeros(1,4*NT);
            P_1(1) = 2;
            P_1(2) = 1;
            P(1,:) = P_1;
            for cc = 2:2*NT
                P_temp = circshift(P_1,2*(cc-1));
                P(cc,:) = P_temp;
            end
            X_hat = P*X_bpsk;
           x_yalmip = 2*binvar(4*NT,1)-1;
           X = [X_hat(1:NT)+1i*X_hat(NT+1:2*NT)]/sqrt(10);
           X_transmit = [X_hat(1:NT)+1i*X_hat(NT+1:2*NT)];
      else
             error('不支持这种调制方式');
      end
        
    N = [randn(NR,1)+1i*randn(NR,1)]/sqrt(2); % 0均值高斯白噪声
    H = [randn(NR,NT)+1i*randn(NR,NT)]/sqrt(2); % 瑞利衰落信道
    
    for ii = 1:length(Eb) %计算多个信噪比情况下的输出Y    
        Y = H*X+(10^(-Eb(ii)/20))*N;%利用上面的信号x，信道H，noise n计算出输出信号Y
        
%         Y_hat = [real(Y);imag(Y)];
%         if M == 2
%             H_hat = [real(H);imag(H)];
%         else
%             H_hat = [real(H),-imag(H);imag(H),real(H)];
%         end
%         
%         N_hat = [real(10^(-Eb(ii)/20)*N);imag(10^(-Eb(ii)/20)*N)];
%         
%         if M == 2
%             GR = zeros(2*NR,NT);
%         else
%             GR = zeros(2*NR,2*NT);
%         end
%         YR = [sign(real(Y)); sign(imag(Y))];
%         for kk = 1:2*NR
%             GR(kk,:)= H_hat(kk,:)*YR(kk);
%         end
%         
%         if M == 2
%                Constraints = [y<=GR*x_yalmip,y<=0];
%         elseif M == 4
%                Constraints = [y<=GR*x_yalmip/sqrt(2),y<=0];
%         elseif M == 16
%                Constraints = [y<=GR*((P*x_yalmip)/sqrt(10)),y<=0];
%         end
%         
%         Objective = -sum(y);
%         sol = optimize(Constraints, Objective, options);
%         
%         if M == 16
%                solution = P*round(value(x_yalmip));
%         else 
%                solution = round(value(x_yalmip));
%         end
%         
%         if M == 2
%                Y_mid_OneBit = solution;
%         else 
%                Y_mid_OneBit = [solution(1:NT) + 1i*solution(NT+1:2*NT)];
%         end

       Y_mid_OneBit = ourTest(Y,H,M);
        
   
      
       if(ii == 1)
        E_f10 = [E_f10,Y_mid_OneBit(1:NT)];  
      end
      
      if(ii == 2)
        E_f5 = [E_f5,Y_mid_OneBit(1:NT)];  
      end
      
      if(ii == 3)
        E_0 = [E_0,Y_mid_OneBit(1:NT)];  
      end
      
      if(ii == 4)
        E_5 = [E_5,Y_mid_OneBit(1:NT)];  
      end
      
      if(ii == 5)
        E_10 = [E_10,Y_mid_OneBit(1:NT)];  
       end
  end
    
        
        l =NT-length(find(E_0(:,jj)==X_transmit));
        errorNumber = errorNumber + l;


       
        l1 =NT-length(find(E_5(:,jj)==X_transmit));
        errorNumber_5 = errorNumber_5 + l1;
            
        
                                                                                                                                                                                                                                                                                 
        l2 = NT-length(find(E_10(:,jj)==X_transmit));
        errorNumber_10 = errorNumber_10 + l2;

        
        
        l3 = NT-length(find(E_f5(:,jj)==X_transmit));
        errorNumber_f5 = errorNumber_f5 + l3;

        
        
        l4 = NT-length(find(E_f10(:,jj)==X_transmit));
        errorNumber_f10 = errorNumber_f10 + l4;
    
    end


SER_0 = errorNumber/(NT*IterationTimes);
SER_5 = errorNumber_5/(NT*IterationTimes);
SER_10 = errorNumber_10/(NT*IterationTimes);
SER_f5 = errorNumber_f5/(NT*IterationTimes);
SER_f10 = errorNumber_f10/(NT*IterationTimes);
% 
SER = [SER_f10,SER_f5,SER_0,SER_5,SER_10];
% end