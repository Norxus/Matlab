clc
close all
clear all

IterationTimes = 10000; % ���͵ķ�����Ŀ
NT= 1;      %����������
NR= 4;      %����������
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
% X_S = [-1,1];
 sigma = 1;
 L = 4;  %ѵ�����ݸ���
 Td = 512;
 Tt = L*(M)^NT;
 Tu = 10*Tt;


 
for jj=1:IterationTimes
      jj
    ip = randn(2*NT,Td)>0; %�ȸ��ʲ���0��1
    X_hat = 2*ip-1; % 0 -> -1; 1 -> 1
    X = [X_hat(1:NT,:)+1i*X_hat(NT+1:2*NT,:)]/sqrt(2);
    N = [randn(NR,1)+1i*randn(NR,1)]/sqrt(2); % 0��ֵ��˹������
    H = [randn(NR,NT)+1i*randn(NR,NT)]/sqrt(2); % ����˥���ŵ�
    H_hat = [real(H),-imag(H);imag(H),real(H)];

    for ii = 1:length(Eb) %���������������µ����Y
%%
        Y = H*X + (10^(-Eb(ii)/20))*N;
        YR = [sign(real(Y));sign(imag(sign(Y)))];  
        
        %����ѵ���׶�
        [x1] = ndgrid(X_S);
        x_train_all = [x1(:)]'/sqrt(2);
        x_D_all = zeros(2*NR,(M)^NT*L);
        x_D_prosibility = zeros(2*NR,(M)^NT);
        x_D_exception = zeros(2*NR,(M)^NT);
        x_DL_temp = zeros(2*NR,L);
        x_D_prosibility_calculate = zeros(2*NR,(M)^NT*L);
       
       %����źż�D
            for j = 1 : (M)^(NT)
                for l = 1 : L
                    Z_train = (10^(-Eb(ii)/20))*[randn(NR,1)+1i*randn(NR,1)]/sqrt(2);
                    y_temp = H*x_train_all(:,j) + Z_train;
                    y_1_kjl = [sign(real(y_temp)); sign(imag(y_temp))];
                    x_D_all(:,l+(j-1)*L) = y_1_kjl;
                end
                x_DL_temp = x_D_all(:,(j-1)*L+1:j*L);
                x_D_exception(:,j) = sign(sum( x_DL_temp,2));
                temp = x_D_exception(:,j);
                temp(find(temp == 0)) = 1;
                x_D_exception(:,j) = temp;
                
                for n = 1 : 2*NR
                    errnumber_train = 0;
                    for l = 1 : L
                        if(temp(n) ~= x_DL_temp(n,l))
                            errnumber_train = errnumber_train + 1;
                        end
                    end
                x_D_prosibility(n,j) = (errnumber_train+1)/(L+2);
                end
                temp = x_DL_temp - x_D_exception(:,j);
               for l = 1 : L
                   for n = 1 : 2*NR
                       if(temp(n,l)~= 0)
                          x_D_prosibility_calculate(n,l+(j-1)*L) =  x_D_prosibility(n,j);
                       else
                          x_D_prosibility_calculate(n,l+(j-1)*L) = 1 - x_D_prosibility(n,j);
                       end
                   end
               end
            end
            
          %δ������ݼ�L
          x_L_prosibility_calculate = zeros(2*NR,Tu*(M)^(NT));
          x_L_p_sum = zeros(1,Tu);
          x_L_p_temp = zeros(1,Tu);
          x_L = YR(:,1:Tu);
          for i = 1 : Tu
              for j = 1 : (M)^(NT)
                   temp = x_L(:,i) - x_D_exception(:,j);
                  for n = 1 : 2*NR
                      if(temp(n) ~= 0)
                          x_L_prosibility_calculate(n,(i-1)*(M)^(NT)+j) = x_D_prosibility(n,j);
                      else
                          x_L_prosibility_calculate(n,(i-1)*(M)^(NT)+j) = 1 - x_D_prosibility(n,j);
                      end
                  end
              end
              
              x_L_p_sum(i) = sum(prod(x_L_prosibility_calculate(:,(i-1)*(M)^(NT)+1:i*(M)^NT)));
              x_L_p_temp(i) = x_L_p_sum(i)/(M^(NT)); 
          end
          
          %�����0�εĺ������
          x_D_p = prod(x_D_prosibility_calculate)/(M^(NT));
          x_L_p = prod(x_L_prosibility_calculate);
         
          logp_theta = sum(log(x_D_p)) +sum(log(x_L_p_temp));
          
%           ���ϴ������ȷ
 %         
          %����Tt��Tu�ĸ���gamma
          x_L_p_gamma = zeros(1,Tu*(M)^NT);
          for i = 1 : Tu
                 x_L_p_gamma((i-1)*(M)^NT+1:i*(M)^NT) = x_L_p((i-1)*(M)^NT+1:i*(M)^NT)/x_L_p_sum(i); 
          end
          
          
          
          %����theta
          theta_posibility = zeros(1,L+Tu);
          e_numerator = zeros(2*NR,(M)^NT);
          theta_posibility_update = zeros(2*NR,(M)^NT);
          theta_exception_update = zeros(2*NR,(M)^NT);
          theta_temp = zeros(1,Tt+Tu);
          theta_posibility(1:L) = ones(1,L);
          
              for j = 1 : (M)^NT
                  temp = [x_D_all(:,(j-1)*L+1:j*L),x_L];                  
                  for k = 1 : Tu
                      theta_posibility(L+k) = x_L_p_gamma((k-1)*(M)^NT+j);
                  end
                    theta_temp = sign(sum(theta_posibility.*temp,2));                    
                    theta_temp(find(theta_temp == 0)) = 1;
                    theta_exception_update(:,j) = theta_temp;
                
                   e_numerator_number = temp - theta_exception_update(:,j);
                   for n = 1 : 2*NR
                         e_numerator_Index = find(e_numerator_number(n,:) ~= 0);
                          e_numerator(n,j) = sum(theta_posibility(e_numerator_Index));
                   end
                   theta_posibility_update(:,j) = e_numerator(:,j)./(sum(theta_posibility));
              end
          
          x_D_prosibility_calculate = zeros(2*NR,(M)^NT*L);
          for i = 1 : (M)^NT
           temp = x_D_all(:,(i-1)*L+1:i*L) - theta_exception_update(:,i);
              for l = 1 : L
                  for j = 1 : 2*NR
                      if(temp(j,l) ~= 0)
                         temp(j,l) = theta_posibility_update(j,i);
                      else
                         temp(j,l) = 1 - theta_posibility_update(j,i);
                      end
                  end
              end
              x_D_prosibility_calculate(:,(i-1)*L+1:i*L) = temp;
          end
          
          x_L_prosibility_calculate = zeros(2*NR,Tu*(M)^(NT));
          x_L_p_temp = zeros(1,Tu);
          x_L_p_sum = zeros(1,Tu);
          for i = 1 : Tu
              for j = 1 : (M)^(NT)
                   temp = x_L(:,i) - theta_exception_update(:,j);
                  for n = 1 : 2*NR
                      if(temp(n) ~= 0)
                          x_L_prosibility_calculate(n,(i-1)*(M)^(NT)+j) = theta_posibility_update(n,j);
                      else
                          x_L_prosibility_calculate(n,(i-1)*(M)^(NT)+j) = 1 - theta_posibility_update(n,j);
                      end
                  end
              end
              x_L_p_sum(i) = sum(prod(x_L_prosibility_calculate(:,(i-1)*(M)^(NT)+1:i*(M)^NT)));
              x_L_p_temp(i) = x_L_p_sum(i)/(M^(NT)); 
          end
        
          %�����1�εĺ������
          x_D_p = prod(x_D_prosibility_calculate)/(M^(NT));%������������ô���и���������
          x_L_p = prod(x_L_prosibility_calculate);
        
          logp_theta_update = sum(log(x_D_p)) +sum(log(x_L_p_temp));
           logp_theta_update - logp_theta;
%            
     while (logp_theta_update - logp_theta) < sigma
         logp_theta_update - logp_theta;
         logp_theta  = logp_theta_update;
                   %����Tt��Tu�ĸ���gamma
          x_L_p_gamma = zeros(1,Tu*(M)^NT);
          for i = 1 : Tu
                 x_L_p_gamma((i-1)*(M)^NT+1:i*(M)^NT) = x_L_p((i-1)*(M)^NT+1:i*(M)^NT)/x_L_p_sum(i); 
          end
          
%           gamma = zeros((M)^NT,Tt+Tu);
%           e_numerator = zeros(2*NR,(M)^NT);
%           theta_posibility_update = zeros(2*NR,(M)^NT);
%           theta_exception_update = zeros(2*NR,(M)^NT);
%           theta_temp = zeros(1,Tt+Tu);
%           X_D = [x_D_all,x_L];
%           
%           for j = 1 : (M)^NT
%               gamma(j,(j-1)*L+1:j*L) = ones(1,L);
%               for i = 1 : Tu
%                   gamma(j,Tt+i) = x_L_p_gamma((i-1)*(M)^NT+j);
%               end
%               theta_exception_update(:,j) = sign(sum(gamma(j,:).*X_D,2));
%               temp = X_D - theta_exception_update(:,j);
%               for n = 1 : 2*NR
%                   temp_in = temp(n,:);
%                   gamma_temp = gamma(j,:);
%                   e_numerator(n,j) = sum(gamma_temp(find(temp_in ~= 0)));
%               end
%               theta_posibility_update(:,j) = e_numerator(:,j)/sum(gamma(j,:));
%           end
          
          
          %����theta
          theta_posibility = zeros(1,L+Tu);
          e_numerator = zeros(2*NR,(M)^NT);
          theta_posibility_update = zeros(2*NR,(M)^NT);
          theta_exception_update = zeros(2*NR,(M)^NT);
          theta_temp = zeros(1,Tt+Tu);
          theta_posibility(1:L) = ones(1,L);
          
              for j = 1 : (M)^NT
                  temp = [x_D_all(:,(j-1)*L+1:j*L),x_L];                  
                  for k = 1 : Tu
                      theta_posibility(L+k) = x_L_p_gamma((k-1)*(M)^NT+j);
                  end
                    theta_temp = sign(sum(theta_posibility.*temp,2));                    
                    theta_temp(find(theta_temp == 0)) = 1;
                    theta_exception_update(:,j) = theta_temp;
                
                   e_numerator_number = temp - theta_exception_update(:,j);
                   for n = 1 : 2*NR
                         e_numerator_Index = find(e_numerator_number(n,:) ~= 0);
                          e_numerator(n,j) = sum(theta_posibility(e_numerator_Index));
                   end
                   theta_posibility_update(:,j) = e_numerator(:,j)./(sum(theta_posibility));
              end
          
          x_D_prosibility_calculate = zeros(2*NR,(M)^NT*L);
          for i = 1 : (M)^NT
           temp = x_D_all(:,(i-1)*L+1:i*L) - theta_exception_update(:,i);
              for l = 1 : L
                  for j = 1 : 2*NR
                      if(temp(j,l) ~= 0)
                         temp(j,l) = theta_posibility_update(j,i);
                      else
                         temp(j,l) = 1 - theta_posibility_update(j,i);
                      end
                  end
              end
              x_D_prosibility_calculate(:,(i-1)*L+1:i*L) = temp;
          end
          
          x_L_prosibility_calculate = zeros(2*NR,Tu*(M)^(NT));
          x_L_p_temp = zeros(1,Tu);
          x_L_p_sum = zeros(1,Tu);
          for i = 1 : Tu
              for j = 1 : (M)^(NT)
                   temp = x_L(:,i) - theta_exception_update(:,j);
                  for n = 1 : 2*NR
                      if(temp(n) ~= 0)
                          x_L_prosibility_calculate(n,(i-1)*(M)^(NT)+j) = theta_posibility_update(n,j);
                      else
                          x_L_prosibility_calculate(n,(i-1)*(M)^(NT)+j) = 1 - theta_posibility_update(n,j);
                      end
                  end
              end
              x_L_p_sum(i) = sum(prod(x_L_prosibility_calculate(:,(i-1)*(M)^(NT)+1:i*(M)^NT)));
              x_L_p_temp(i) = x_L_p_sum(i)/(M^(NT)); 
          end
        
          %�����1�εĺ������
          x_D_p = prod(x_D_prosibility_calculate)/(M^(NT));%������������ô���и���������
          x_L_p = prod(x_L_prosibility_calculate);
        
          logp_theta_update = sum(log(x_D_p)) +sum(log(x_L_p_temp));
           logp_theta_update - logp_theta;     
     end
     
      X_finalDetection = zeros(2*NT,Td);
%       theta_exception_update = x_D_exception;
      for m = 1 : Tu
         [value,detectionIndex] = max(x_L_p_gamma((m-1)*(M)^NT+1 : m*(M)^NT));
         X_finalDetection_mid = [real(x_train_all(:,detectionIndex));imag(x_train_all(:,detectionIndex))];
         X_finalDetection(:,m) = X_finalDetection_mid;
      end
%       
      for k = Tu+1 : Td
          detection_different = theta_exception_update - YR(:,k);
          detection_posibility = zeros(2*NR,(M)^NT);
%           theta_posibility_update = x_D_prosibility;
          for i = 1 : (M)^NT
              for j = 1 : 2*NR
                  if(detection_different(j,i) ~= 0)
                      detection_posibility(j,i) = theta_posibility_update(j,i);
                  else
                      detection_posibility(j,i) = 1 - theta_posibility_update(j,i);
                  end
              end
          end

         [value,detectionIndex] = max(prod(detection_posibility));
          X_finalDetection_mid = [real(x_train_all(:,detectionIndex));imag(x_train_all(:,detectionIndex))];
          X_finalDetection(:,k) = X_finalDetection_mid;
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
