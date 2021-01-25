clc
close all
clear all

IterationTimes = 10000; % ���͵ķ�����Ŀ
NT=2;      %����������
NR=6;      %����������
Eb=0:5:30;     %�����


%�ܴ������ָ���
errorNumber_4 =0;
errorNumber_6 =0;
errorNumber_8 =0;
errorNumber_10 =0;
errorNumber_12 =0;
errorNumber_14 =0;
errorNumber_16 =0;

%ÿ�ε����Ĵ�����������
 l4=0;
 l6=0;
 l8=0;
 l10=0;
 l12=0;
 l14=0;
 l16=0;
 M=4; %���Ʒ�ʽM-QAM
 X_S = [1+1i,1-1i,-1+1i,-1-1i]; %���Է��͵�������Ԫ
 sigma = 30;
 T = 8;  %ѵ�����ݸ���
 Td = 512;%������ݸ���
 Tt = T*(M)^NT;%ѵ�����ݸ���
 Tu = 10*Tt;


 
for jj=1:IterationTimes
      jj
    ip = randn(2*NT,Td)>0; %�ȸ��ʲ���0��1
    X_hat = 2*ip-1; % 0 -> -1; 1 -> 1
    X = [X_hat(1:NT,:)+1i*X_hat(NT+1:2*NT,:)]/sqrt(2);%��������
    N = [randn(NR,1)+1i*randn(NR,1)]/sqrt(2); % 0��ֵ��˹������
    H = [randn(NR,NT)+1i*randn(NR,NT)]/sqrt(2); % ����˥���ŵ�
    H_hat = [real(H),-imag(H);imag(H),real(H)];


    for ii = 1:length(Eb) %���������������µ����Y
%%
        Y = H*X + (10^(-Eb(ii)/20))*N;
        YR = [sign(real(Y));sign(imag(sign(Y)))];%���������
        
        %����ѵ���׶�
        [x2,x1] = ndgrid(X_S,X_S);
        x_train_all = [x2(:),x1(:)]'/sqrt(2);%���п��ܵĴ�������
        x_D_all = zeros(2*NR,(M)^NT*T);%���п��ܴ�����������õ������ս�������
        x_D_prosibility = zeros(2*NR,(M)^NT);%������ʲ���
        x_D_exception = zeros(2*NR,(M)^NT);%���ڶԱȵ�c����
        x_DL_temp = zeros(2*NR,T);
        x_D_prosibility_calculate = zeros(2*NR,(M)^NT*T);
       
       %����źż�D
            for j = 1 : (M)^(NT)%�����������п��ܵ����
                for l = 1 : T %ÿ�ֿ��ܵ���������ظ�����T��
                    Z_train = (10^(-Eb(ii)/20))*[randn(NR,1)+1i*randn(NR,1)]/sqrt(2);%ÿһ�������Ӧ��ͬ����������.
                    y_temp = H*x_train_all(:,j) + Z_train;
                    y_1_kjl = [sign(real(y_temp));sign(imag(y_temp))];
                    x_D_all(:,l+(j-1)*T) = y_1_kjl;
                end %���forѭ�������󣬾��൱��ѵ����һ�������
                x_DL_temp = x_D_all(:,(j-1)*T+1:j*T);
                x_D_exception(:,j) = sign(sum( x_DL_temp,2));%�������ڶԱȵ�c����
                temp_exception = x_D_exception(:,j);
                temp_exception(find(temp_exception == 0)) = 1;
                x_D_exception(:,j) = temp_exception;%������������Ҫ�ǽ�sum�õ���0����1����Ϊsign�����޷����д����滻
                
                for n = 1 : 2*NR %Ȼ�������forѭ������������
                    errnumber_train = 0;
                    for l = 1 : T
                        if(temp_exception(n) ~= x_DL_temp(n,l))
                            errnumber_train = errnumber_train + 1;
                        end
                    end
                x_D_prosibility(n,j) = (errnumber_train+1)/(T+2);
                end
            end
            

     
      X_finalDetection = zeros(2*NT,Td);
      theta_exception_update = x_D_exception;
%       
      for k = 1 : Td %�õ���������֮�󣬽��м��
          detection_different = theta_exception_update - YR(:,k);%���������0�����ݾ;���Ϊ���Ӧλ����ȣ���ÿһ����������������п��ܵ�������жԱ�
          detection_posibility = zeros(2*NR,(M)^NT);
          theta_posibility_update = x_D_prosibility;
          for i = 1 : (M)^NT
              for j = 1 : 2*NR
                  if(detection_different(j,i) ~= 0)
                      detection_posibility(j,i) = theta_posibility_update(j,i);
                  else
                      detection_posibility(j,i) = 1 - theta_posibility_update(j,i);
                  end
              end
          end

         [value,detectionIndex] = max(prod(detection_posibility));%���˲������ҵ���Ӧ������
          X_finalDetection_mid = [real(x_train_all(:,detectionIndex));imag(x_train_all(:,detectionIndex))];
          %����˳�����ѵ������ô�Ա����ֵ�˳��Ҳ����ͬ�ģ�
          X_finalDetection(:,k) = X_finalDetection_mid;
%         X_finalDetection(:,k) = x_train_all(:,detectionIndex);
      end
      
       if(ii == 1)
        E_4 = X_finalDetection;  
      end
      
      if(ii == 2)
        E_6 = X_finalDetection;  
      end
      
      if(ii == 3)
        E_8 = X_finalDetection;  
      end
      
      if(ii == 4)
        E_10 = X_finalDetection;  
      end
      
      if(ii == 5)
        E_12 = X_finalDetection;  
      end
       
      if(ii == 6)
        E_14 = X_finalDetection;  
      end
      
      if(ii == 7)
        E_16 = X_finalDetection;  
      end
 end
    
        
        l4 =length(find(E_4~=X_hat/sqrt(2)));
        errorNumber_4 = errorNumber_4 + l4;


       
        l6 =length(find(E_6~=X_hat/sqrt(2)));
        errorNumber_6 = errorNumber_6 + l6;
            
        
                                                                                                                                                                                                                                                                                 
        l8 =length(find(E_8~=X_hat/sqrt(2)));
        errorNumber_8 = errorNumber_8 + l8;

        
        
        l10 =length(find(E_10~=X_hat/sqrt(2)));
        errorNumber_10 = errorNumber_10 + l10;

        
        
        l12 =length(find(E_12~=X_hat/sqrt(2)));
        errorNumber_12 = errorNumber_12 + l12;
        
        l14 =length(find(E_14~=X_hat/sqrt(2)));
        errorNumber_14 = errorNumber_14 + l14;
        
        l16 =length(find(E_16~=X_hat/sqrt(2)));
        errorNumber_16 = errorNumber_16 + l16;
    
    end

SER_4 = errorNumber_4/(2*NT*IterationTimes*Td);
SER_6 = errorNumber_6/(2*NT*IterationTimes*Td);
SER_8 = errorNumber_8/(2*NT*IterationTimes*Td);
SER_10 = errorNumber_10/(2*NT*IterationTimes*Td);
SER_12 = errorNumber_12/(2*NT*IterationTimes*Td);
SER_14 = errorNumber_14/(2*NT*IterationTimes*Td);
SER_16 = errorNumber_16/(2*NT*IterationTimes*Td);
% 
SER = [SER_4,SER_6,SER_8,SER_10,SER_12,SER_14,SER_16];
x = 0:5:30;
semilogy(x,SER,'-o');
grid on
