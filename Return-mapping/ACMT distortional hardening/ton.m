function  [t_on,Qton,qton,d_ton] = ton(Q_dot1,Q_now,Qb_now,c, q_now,q_dot1,tau,tton,P)
syms x
  eqn=1.5.*(1-c*(Qb_now'*P*(Q_now-Qb_now+x.*Q_dot1))/(sqrt(Q_now-Qb_now+x.*Q_dot1)'*P*(Q_now-Qb_now+x.*Q_dot1)))*((Q_now-Qb_now+x.*Q_dot1)'*P*(Q_now-Qb_now+x.*Q_dot1))-tau.^2==0;
        d_ton=double(max(solve(eqn,x)));
        t_on=double(max(solve(eqn,x)))+tton;
         Qton=Q_now+Q_dot1.*(t_on-tton);
         qton=q_now+(q_dot1)*d_ton;
         
end