function [QQnext] = return_mapping(Qnow,K,q_dot1,Qy,dt)
      Qtrail=Qnow+(-K*(Qnow'*q_dot1)/Qy^2*Qnow+K*q_dot1)*dt;
      ntrail=Qtrail./norm(Qtrail);
      QQnext=Qy*ntrail;                    
end