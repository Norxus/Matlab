clc
clear
IterationTimes = 10000; % ���͵ķ�����Ŀ
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
 c_TwoStage_Final = 1.3;

for jj=1:IterationTimes
    
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
        Y_hat = [real(Y);imag(Y)];
        H_hat = [real(H),-imag(H);imag(H),real(H)];
        N_hat = [real(10^(-Eb(ii)/20)*N);imag(10^(-Eb(ii)/20)*N)];
        
        GR_T = [];
        GR = [real(H.'),imag(H.');-imag(H.'),real(H.')]';
        GR_NEW = mat2cell(GR,ones(1,2*NR),2*NT);
        GR_OLD = GR_NEW;
        YR = [sign(real(Y)); sign(imag(Y))];
        for kk = 1:2*NR
            GR_NEW{kk,1}= GR_NEW{kk,1}*YR(kk);
            GR_T = cat(2,GR_T,GR_NEW{kk,1}');
        end
        G_R = GR_NEW;

%         Final = 0;
%         X_RML = [];
%          for k = 1:2^(2*NT)
%             Final = 0;
%             X_transmit = 2*bitget(k,2*NT:-1:1)'-1;
%             for n = 1:2*NR
%               Final = Final + log(normcdf(sqrt(2*p)*G_R{1,n}'*X_transmit,0,1));
%             end
%           X_RML = [X_RML,Final];
%         end
%         [X_max,X_index] = max(X_RML);
%         X_reserve = 2*bitget(X_index,2*NT:-1:1)-1;
%         if(ii == 1)   
%             E_f10 = [E_f10,X_reserve'];
%         end
%         if(ii ==2)
%             E_f5 = [E_f5, X_reserve'];
%         end
%         if(ii == 3)
%             E_0 = [E_0, X_reserve'];
%         end
%         if(ii == 4)
%             E_5 = [E_5, X_reserve'];
%         end
%         if(ii == 5)
%             E_10 = [E_10, X_reserve'];
%         end
        X_0 =  sqrt(NT)*(GR_T*ones(2*NR,1)/norm(GR_T*ones(2*NR,1)));
        f = [];
        k_final = 0.01;
        c_final = 10^(-3);
        for TT = 1:2*NR
             temp = (1/sqrt(2*pi))*(exp((-p)*abs(G_R{TT,1}*X_0)^2))/normcdf(sqrt(2*p)*G_R{TT,1}*X_0,0,1);
             f = [f,temp];
        end
        X_1 = X_0 + k_final*GR_T*f';
        if(norm(X_1)^2 > NT)
               X_1 = sqrt(NT)*(X_1/norm(X_1));
        end
        R(:,1) = X_0;
        R(:,2) = X_1;
        CN = 2;
        while(norm(R(:,CN) - R(:,CN-1))>= c_final*norm(R(:,CN-1)))
          f_while = [];
          for T = 1:2*NR
             temp = (1/sqrt(2*pi))*(exp((-p)*abs(G_R {T,1}*R(:,CN-1))^2))/normcdf(sqrt(2*p)*G_R{T,1}*R(:,CN-1),0,1);
             f_while = [f_while,temp];
          end
            CN=CN+1;
            R(:,CN) = R(:,CN-1) +k_final*GR_T*f_while';
            if(norm(R(:,CN))^2 > NT)
                R(:,CN) = sqrt(NT)*(R(:,CN)/norm(R(:,CN)));
            end
        end
%     Y_OneBit = sign(Y_hat);
%     Y_OneBit = sign(R);
      Y_mid_OneBit = sqrt(NT)*(1/norm(R(:,CN)))*(R(:,CN));
       X_twoSatgenML_variableIndex = [];

      for i_mid = 1:NT
          Y_mid_OneBit(i_mid) = Y_mid_OneBit(i_mid)+1i*Y_mid_OneBit(i_mid+NT);       
          X_nML = Y_mid_OneBit(1:NT);
         
          mid = [];
          for i_inmid = 1:M
              midNumber = abs(Y_mid_OneBit(i_mid)-X_S(i_inmid));
              mid=[mid,midNumber];
          end
          Y_mid_OneBit(i_mid) = X_S(find(mid == min(mid)));
          Y_Denominator = abs(X_nML(i_mid) - Y_mid_OneBit(i_mid));
          
          for i_twoStage_inmid = 1:M
             Judge =  abs(mid(i_twoStage_inmid))/ Y_Denominator;
             if(Judge < c_TwoStage_Final)
                 Judge_twoStage(i_twoStage_inmid) = 1;
             else
                 Judge_twoStage(i_twoStage_inmid) = 0;
             end
          end
          
          X_twoSatgenML_variableIndex =  [ X_twoSatgenML_variableIndex,Judge_twoStage'];
              
      end
      
     X1_Index = find(X_twoSatgenML_variableIndex(:,1)==1);
     X2_Index = find(X_twoSatgenML_variableIndex(:,2)==1);
     X3_Index = find(X_twoSatgenML_variableIndex(:,3)==1);
     X4_Index = find(X_twoSatgenML_variableIndex(:,4)==1);
     
     X_CandidateMatrix = [];
     
     for i1 = 1:length(X1_Index)
         X1_temp_TwoStage = X_S(X1_Index(i1));
         for i2 = 1:length(X2_Index)
             X2_temp_TwoStage = X_S(X2_Index(i2));
             for i3 = 1:length(X3_Index)
                 X3_temp_TwoStage = X_S(X3_Index(i3));
                 for i4 = 1:length(X4_Index)
                     X4_temp_TwoStage = X_S(X4_Index(i4));
                     X_Candidate = [X1_temp_TwoStage,X2_temp_TwoStage,X3_temp_TwoStage,X4_temp_TwoStage];
                     X_CandidateMatrix = [X_CandidateMatrix,X_Candidate.'];
                 end
             end
         end
     end
    
    X_CandidateMatrix_R = [real(X_CandidateMatrix);imag(X_CandidateMatrix)];
    X_estimate = [];
    for k = 1:size(X_candidateMatrix_R,2)
        X_Two_nML = 1;
     for i = 1:2*NR
            X_Two_nML = X_Two_nML* normcdf(sqrt(2*p)*G_R{i,1}*X_CandidateMatrix_R(:,k),0,1);
     end
      X_estimate = [X_estimate,X_Two_nML];
    end        
             
     X_TwoStage_FinalOutput = X_CandidateMatrix(:,find(X_estimate == max(X_estimate)));
%       if(ii == 1)
%         E_f10 = [E_f10,FinalOutput];  
%       end
%       
%       if(ii == 2)
%         E_f5 = [E_f5,FinalOutput];  
%       end
%       
%       if(ii == 3)
%         E_0 = [E_0,FinalOutput];  
%       end
%       
%       if(ii == 4)
%         E_5 = [E_5,FinalOutput];  
%       end
%       
%       if(ii == 5)
%         E_10 = [E_10,FinalOutput];  
%       end
       if(ii == 1)
        E_f10 = [E_f10,X_TwoStage_FinalOutput];  
      end
      
      if(ii == 2)
        E_f5 = [E_f5,X_TwoStage_FinalOutput];  
      end
      
      if(ii == 3)
        E_0 = [E_0,X_TwoStage_FinalOutput];  
      end
      
      if(ii == 4)
        E_5 = [E_5,X_TwoStage_FinalOutput];  
      end
      
      if(ii == 5)
        E_10 = [E_10,X_TwoStage_FinalOutput];  
       end
  end
    
% 
%         l =2*NT-length(find(E_0(:,jj)==X_hat));
%         errorNumber = errorNumber + l;
% 
% 
%        
%         l1 =2*NT-length(find(E_5(:,jj)==X_hat));
%         errorNumber_5 = errorNumber_5 + l1;
%             
%         
%                                                                                                                                                                                                                                                                                  
%         l2 = 2*NT-length(find(E_10(:,jj)==X_hat));
%         errorNumber_10 = errorNumber_10 + l2;
% 
%         
%         
%         l3 =2*NT-length(find(E_f5(:,jj)==X_hat));
%         errorNumber_f5 = errorNumber_f5 + l3;
% 
%         
%         
%         l4 =2*NT-length(find(E_f10(:,jj)==X_hat));
%         errorNumber_f10 = errorNumber_f10 + l4;
        
        
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


% 
% SER_0 = (1/NT)*errorNumber/(2*NT*IterationTimes);
% SER_5 = (1/NT)*errorNumber_5/(2*NT*IterationTimes);
% SER_10 = (1/NT)*errorNumber_10/(2*NT*IterationTimes);
% SER_f5 = (1/NT)*errorNumber_f5/(2*NT*IterationTimes);
% SER_f10 = (1/NT)*errorNumber_f10/(2*NT*IterationTimes);

SER_0 = errorNumber/(NT*IterationTimes);
SER_5 = errorNumber_5/(NT*IterationTimes);
SER_10 = errorNumber_10/(NT*IterationTimes);
SER_f5 = errorNumber_f5/(NT*IterationTimes);
SER_f10 = errorNumber_f10/(NT*IterationTimes);
% 
SER = [SER_f10,SER_f5,SER_0,SER_5,SER_10];
x = [-10,-5,0,5,10];
semilogy(x,SER,'-o');
grid on
