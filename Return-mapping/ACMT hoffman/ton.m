function [t_on,Qt,qt,d_ton] = ton(Q_dot1,Q_off,q_dot1,r,tton,Y,b, q_off)
syms x
  eqn=0.5*(Q_off+x.*Q_dot1)'*Y*(Q_off+x.*Q_dot1)+b'*(Q_off+x.*Q_dot1).*r-0.5.*r.^2==0;
        d_ton=double(max(solve(eqn,x)));
        t_on=double(max(solve(eqn,x)))+tton;
         Qt=Q_off+Q_dot1.*(t_on-tton);
         qt=q_off+(q_dot1)*d_ton;
         
end