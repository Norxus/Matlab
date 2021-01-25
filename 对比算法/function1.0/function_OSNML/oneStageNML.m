function [X] = oneStageNML(Y,H,M,Eb)
NR = size(H,1);
NT = size(H,2);
p = 10^(Eb/10);
if M == 2
      X_S = [1,-1];
elseif M == 4
     X_S = [1+1i,1-1i,-1+1i,-1-1i]/sqrt(2);
elseif M ==16
    X_S = [3+1i,3-1i,3-3i,3+3i,-3-3i,-3+3i,-3-1i,-3+1i,-1+1i,-1-1i,-1+3i,-1-3i,1+3i,1-3i,1+1i,1-1i]/sqrt(10);
else
    error('不支持该调制方式');
end

        if M == 2
            H_hat = [real(H);imag(H)];
        else
            H_hat = [real(H),-imag(H);imag(H),real(H)];
        end
        
        
        GR_T = [];
%         if M==2
%             GR = [real(H);imag(H)];
%         else
%             GR = [real(H.'),imag(H.');-imag(H.'),real(H.')]';
%         end
        GR = H_hat;
        if M == 2
            GR_NEW = mat2cell(GR,ones(1,2*NR),NT);
        else
            GR_NEW = mat2cell(GR,ones(1,2*NR),2*NT);
        end
        
        GR_OLD = GR_NEW;
        YR = [sign(real(Y)); sign(imag(Y))];
        
        for kk = 1:2*NR
            GR_NEW{kk,1}= GR_NEW{kk,1}*YR(kk);
            GR_T = cat(2,GR_T,GR_NEW{kk,1}');
        end
        G_R = GR_NEW;

        X_0 =  sqrt(NT)*(GR_T*ones(2*NR,1)/norm(GR_T*ones(2*NR,1)));
        f = [];
        if M == 2
            R = zeros(NT,600);
        else
            R = zeros(2*NT,600);
        end
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
%             norm(R(:,CN) - R(:,CN-1))
          f_while = [];
          for T = 1:2*NR
             temp = (1/sqrt(2*pi))*(exp((-p)*abs(G_R {T,1}*R(:,CN))^2))/normcdf(sqrt(2*p)*G_R{T,1}*R(:,CN),0,1);
             f_while = [f_while,temp];
          end
            CN=CN+1;
            R(:,CN) = R(:,CN-1) +k_final*GR_T*f_while';
%             if(CN>100||prod(isnan(R(:,CN))))
%                 break;
%             end
            if(norm(R(:,CN))^2 > NT)
                R(:,CN) = sqrt(NT)*(R(:,CN)/norm(R(:,CN)));
            end
        end

%             if(CN>100||prod(isnan(R(:,CN))))
%                 break;
%             end 
      Y_mid_OneBit = sqrt(NT)*(1/norm(R(:,CN)))*(R(:,CN));

      for i_mid = 1:NT
          if M ~= 2
              Y_mid_OneBit(i_mid) = Y_mid_OneBit(i_mid)+1i*Y_mid_OneBit(i_mid+NT);
          end
          mid = [];
          for i_inmid = 1:M
              midNumber = abs(Y_mid_OneBit(i_mid)-X_S(i_inmid));
              mid=[mid,midNumber];
          end
          Y_mid_OneBit(i_mid) = X_S(find(mid == min(mid)));
      end
      X = Y_mid_OneBit;
end