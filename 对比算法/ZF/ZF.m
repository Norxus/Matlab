clc
clear

IterationTimes = 1000; % ���͵ķ�����Ŀ
NT=4;      %����������
NR=32;      %����������
Eb=-10:5:10;     %�����
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
 M=4;
 X_S = [1+1i,1-1i,-1+1i,-1-1i];

for jj=1:IterationTimes
    jj
    
    ip = randn(2*NT,1)>0; %�ȸ��ʲ���0��1
    X_hat = 2*ip-1; % 0 -> -1; 1 -> 1
    
    X_transmit = [X_hat(1:NT)+1i*X_hat(NT+1:2*NT)];
    X = [X_hat(1:NT)+1i*X_hat(NT+1:2*NT)]/sqrt(2);
    N = [randn(NR,1)+1i*randn(NR,1)]/sqrt(2); % 0��ֵ��˹������
    H = [randn(NR,NT)+1i*randn(NR,NT)]/sqrt(2); % ����˥���ŵ�
    
    for ii = 1:length(Eb) %���������������µ����Y
        p = 10^(Eb(ii)/10);
        Y = sqrt(p)*H*X+N;%����������ź�x���ŵ�H��noise n���������ź�Y
%       Y = sqrt(2p)*H*X_transmit+(10^(-Eb(ii)/20))*N;
        Y_hat = [sign(real(Y));sign(imag(Y))];
        Y = [Y_hat(1:NR)+1i*Y_hat(NR+1:2*NR)];
        H_hat = [real(H),-imag(H);imag(H),real(H)];
        N_hat = [real(10^(-Eb(ii)/20)*N);imag(10^(-Eb(ii)/20)*N)];
        


      Y_mid_OneBit = pinv(H)*Y;  
      Y_mid_OneBit = sqrt(NT)*Y_mid_OneBit/norm(Y_mid_OneBit);
      for i_mid = 1:NT
%           Y_mid_OneBit(i_mid) = Y_mid_OneBit(i_mid)+1i*Y_mid_OneBit(i_mid+NT);
          mid = [];
          for i_inmid = 1:M
              midNumber = abs(Y_mid_OneBit(i_mid)-X_S(i_inmid));
              mid=[mid,midNumber];
          end
          Y_mid_OneBit(i_mid) = X_S(find(mid == min(mid)));
      end
      
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
x = 5:5:25;
semilogy(x,SER,'-o');
grid on
% hold on
% ZF_SER = [0.129015486725664,0.0461946902654867,0.0128871681415929,0.00554203539823009,0.00246681415929204];
% semilogy(x,ZF_SER,'-o');
% hold on
% ZF_MCD_SER = [0.0653761061946903,0.0274115044247788,0.0108075221238938,0.00217920353982301,0.00119469026548673];
% semilogy(x,ZF_MCD_SER,'-o');
