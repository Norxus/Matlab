function [X] = ZF(Y,H,M)
NR = size(H);
NT = size(H,2);

if M == 2
      X_S = [1,-1];
elseif M == 4
     X_S = [1+1i,1-1i,-1+1i,-1-1i]/sqrt(2);
elseif M ==16
    X_S = [3+1i,3-1i,3-3i,3+3i,-3-3i,-3+3i,-3-i,-3+i,-1+i,-1-i,-1+3i,-1-3i,1+3i,1-3i,1+i,1-i]/sqrt(10);
else
    error('不支持该调制方式');
end

        Y_hat = [sign(real(Y));sign(imag(Y))];
        
        if M == 2
            H_hat = [real(H);imag(H)];
        else
            H_hat = [real(H),-imag(H);imag(H),real(H)];
        end
        
        


      Y_mid_OneBit = pinv(H)*Y;  
      Y_mid_OneBit = sqrt(NT)*Y_mid_OneBit/norm(Y_mid_OneBit);
      for i_mid = 1:NT
          mid = [];
          for i_inmid = 1:M
              midNumber = abs(Y_mid_OneBit(i_mid)-X_S(i_inmid));
              mid=[mid,midNumber];
          end
          Y_mid_OneBit(i_mid) = X_S(find(mid == min(mid)));
      end
      
      X = Y_mid_OneBit;

end

