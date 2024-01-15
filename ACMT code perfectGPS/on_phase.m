function [nextQ,nextq,nextdq,nextdQ,nextyf,nextqp,nextdqp]=on_phase(q,ke,dt,Qy,Q,qp,dq,qy,ee,ii)
nextq=q(ii,:); 
nextdq=(q(ii,:)-q(ii-1,:))/dt;
 a=cosh((norm(nextdq)*dt)/qy);
 b=sinh((norm(nextdq)*dt)/qy);
 G=[eye(6)+((a-1)*(nextdq')*(nextdq))/(norm(nextdq))^2 (b*nextdq')/norm(nextdq) ; b*(nextdq)/norm(nextdq) a];
 g=[eye(6) zeros(6,1) ; zeros(1,6) ,-1];
 a1=G'*g*G
 X=[ Q'/Qy ; 1 ];
 X0=G*X; 
 nextQ=((X0(1:6)/X0(7))*Qy)';
%%
% nextdq=(q(ii,:)-q(ii-1,:))/dt;
  dlambda=(nextQ*nextdq')/Qy;
  nextdQ=nextdq*ke-ke*dlambda*nextQ/Qy;
  nextdqp=((nextQ*nextdq')/Qy)*nextQ/Qy^2;
  nextqp=(nextdqp*dt+qp);
 nextyf=norm(nextQ)-Qy;
%  nextyf=norm((ke*q(ii,:)+ke*nextdq*dt))-Qy;
%   nextyf(abs(nextyf)<ee)=0;

 
 





