function [nextQ,nextq,nextdq,nextdQ,nextQb,nextQa,nextyf,nextqp,nextdqp]=on_phase(q,ke,kp,dt,Qy,Q,Qb,qp,Qa,qy,beta,ee,ii)
nextQ=Q(ii,:);
nextdQ=(Q(ii,:)-Q(ii-1,:))/dt;
 a=cosh((norm(nextdQ)*dt)/Qy);
 b=sinh((norm(nextdQ)*dt)/Qy);
 G=round([eye(6)+((a-1)*(nextdQ')*(nextdQ))/(norm(nextdQ))^2 (b*nextdQ')/norm(nextdQ) ; b*(nextdQ)/norm(nextdQ) a],1);
 g=[eye(6) zeros(6,1) ; zeros(1,6) ,-1];
 a1=G'*g*G;
 X=[ Qa'/Qy ; 1 ];
 X0=G*X; 
 nextQa=((X0(1:6)/X0(7))*Qy)';
%%
  nextQa=Qa;
  dlambda=(nextQa*nextdQ')/(kp*Qy);
  nextdq=nextdQ/ke+dlambda*nextQa/Qy;
  nextq=(nextdq*dt+q);
  nextdqp=dlambda*nextQa'/Qy;
  nextqp=(nextdqp*dt+qp);
  nextQb=kp*nextqp;
  nextyf=norm(nextQa)-Qy;
%   nextyf(abs(nextyf)<ee)=0;

 
 





