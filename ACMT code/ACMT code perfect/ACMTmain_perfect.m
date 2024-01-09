clc;
clear all;
close all;
%% parameters
dt=0.005;
ke=10;
Qy=4000;
ee=1e-10;
%Qy=1;
%% input data
t=0:dt:20;
q11=(t.^2)';
q12=(2*t.^2)';
q21=(3*t.^2)';
q22=(4*t.^2)';
q31=(5*t.^2)';
q32=(6*t.^2)';
q=[ q11  q12  q21  q22  q31  q32 ];
q(1,:)=[0 ; 0 ; 0 ; 0 ; 0 ; 0 ];
%% A matrix
qy=Qy/ke;
%% initial values
Q=[0 ; 0 ; 0 ; 0 ; 0 ; 0 ]';
dq=[0 ; 0 ; 0 ; 0 ; 0 ; 0 ];
dQ=[0 ; 0 ; 0 ; 0 ; 0 ; 0 ];
qp=[0 ; 0 ; 0 ; 0 ; 0 ; 0 ];
qe=[0 ; 0 ; 0 ; 0 ; 0 ; 0 ];
dqp=[0 ; 0 ; 0 ; 0 ; 0 ; 0 ];
% straining condition
s=0;
% switching index
iflag=[0];
state1=[0 dq'];
state2=[0 Q dQ' norm(Q)-Qy qp' q11(1) q12(1) q21(1) q22(1) q31(1) q32(1) iflag];

%%
for ii=2:length(q11)
    %ggg=ii;
       if iflag==[0];
        % off_phase
        [nextQ,nextq,nextqe,nextdq,nextdQ,nextyf,nextqp]=off_phase(Q,ke,dt,Qy,q,qe,ee,ii);
       if nextyf>0
            % off to on
            [ts,nextQ,nextq,nextqe,nextdq,nextdQ,nextyf,nextqp] = pull_back (iflag,dQ,q,qe,ke,dt,dq,Qy,Q,ee,ii);
           
             q1=q(ii-1,:)+dq*ts; Q=nextQ; dq=nextdq; dQ=nextdQ; qp=nextqp; qe=nextqe; nextyf=0;
            state1=[state1;(ii-2)*dt+ts dq ];
            state2=[state2; [(ii-2)*dt+ts Q dQ nextyf qp q1 iflag]];
   
            iflag=[1];
            [nextQ,nextq,nextdq,nextdQ,nextyf,nextqp,nextdqp]=on_phase(q,ke,dt-ts,Qy,Q,qp,dq,qy,ee,ii);
       else
           iflag=[0];
       end
       elseif iflag==[1];
             s=(Q*dq')/(Qy);
             if s>0 
             iflag=[1];
            [nextQ,nextq,nextdq,nextdQ,nextyf,nextqp,nextdqp]=on_phase(q,ke,dt,Qy,Q,qp,dq,qy,ee,ii);
               else 
                   iflag=[0];

                     [nextQ,nextq,nextqe,nextdq,nextdQ,nextyf,nextqp]=off_phase(Q,ke,dt,Qy,q,qe,ee,ii);
            
                   if nextyf>0
            % off to on
                    [ts,nextQ,nextq,nextqe,nextdq,nextdQ,nextyf,nextqp] = pull_back (iflag,dQ,q,qe,ke,dt,dq,Qy,Q,ee,ii);
                       q1=q(ii-1,:)+dq*ts; Q=nextQ; dq=nextdq; dQ=nextdQ; qp=nextqp; qe=nextqe; nextyf=0;
                    state1=[state1;(ii-2)*dt+ts dq ];
                    state2=[state2; [(ii-2)*dt+ts Q dQ nextyf qp q1 iflag]];
                   
                    iflag=[1];
                    [nextQ,nextq,nextdq,nextdQ,nextyf,nextqp,nextdqp]=on_phase(q,ke,dt-ts,Qy,Q,qp,dq,qy,ee,ii);
      
           
                   else
                    iflag=[0];
                   end

                   
               end
       end
    q1=nextq;
    Q=nextQ;
    qe=nextqe;
    dq=nextdq;
    dQ=nextdQ;
    qp=nextqp;
    state1=[state1;(ii-1)*dt dq];
    state2=[state2;[(ii-1)*dt Q dQ nextyf qp q1 iflag]];

end

figure(1)
axes('FontSize',14,'FontName','TimesNewRoman');
plot(state2(:,21),state2(:,2),'linewidth',1.0)
title('q-Q');
xlabel('q');
ylabel('Q');


% figure(2)
% axes('FontSize',14,'FontName','TimesNewRoman');
% plot(state2(:,1),state2(:,5),'linewidth',1.0)
% title('X-direction Displacement Response');
% xlabel('Time (s)');
% ylabel('Displacement  q_x (cm)');
% 
% figure(3)
% axes('FontSize',14,'FontName','TimesNewRoman');
% plot(state2(:,1),state2(:,8),'linewidth',1.0)
% title('X-direction Restoring-Force Response');
% xlabel('Time (s)');
% ylabel(' Restoring-force  Q_x (kgf)');
% 
% figure(4)
% axes('FontSize',14,'FontName','TimesNewRoman');
% plot(state2(:,1),state2(:,11),'linewidth',1.0)
% title('X-direction Back-Force Response');
% xlabel('Time (s)');
% ylabel('Back_Force  Qb_x (kgf)');
% 
% figure(5)
% axes('FontSize',14,'FontName','TimesNewRoman');
% plot(state2(:,1),state2(:,14),'linewidth',1.0)
% title('X-direction Active-Force Response');
% xlabel('Time (s)');
% ylabel('Active_Force  Qa_x (kgf)');
% 
% 
% 
% 
% figure(6)
% axes('FontSize',14,'FontName','TimesNewRoman');
% plot(state2(:,5),state2(:,8)-c*state2(:,2)/m,'linewidth',1.0)
% title('X-direction Q_x-q_x Response');
% xlabel('Displacement  q_x (cm)');
% ylabel('Total force  Q_x (kgf)');
% 
% figure(7)
% axes('FontSize',14,'FontName','TimesNewRoman');
% plot(state2(:,5),state2(:,8),'linewidth',1.0)
% title('X-direction Q_x-q_x Response');
% xlabel('Displacement  q_x (cm)');
% ylabel('Total force  Q_x (kgf)');
% 
%  figure(8)
%  axes('FontSize',14,'FontName','TimesNewRoman');
%  plot(state2(:,1),state2(:,21),'linewidth',1.0)
%  title('Iflag');
%  xlabel('Time(s)');
%  ylabel('Iflag');
% t1=[0:dt:dt*(length(ag)-1)]';
% 
%  figure(9)
%  axes('FontSize',14,'FontName','TimesNewRoman');
%  plot(t1,ag(:,1),'linewidth',1.0)
%  title('Input force');
%  xlabel('Time(s)');
%  ylabel('Force');
% 
% figure(10)
% axes('FontSize',14,'FontName','TimesNewRoman');
% plot(state2(:,1),state2(:,2)/m,'linewidth',1.0)
% title('X-direction velocity Response');
% xlabel('Time (s)');
% ylabel('velocity  v_x (kg*cm/s)');
% 
% figure (11)
% axes('FontSize',14,'FontName','TimesNewRoman');
% plot(state1(:,1),state1(:,2)/m,'linewidth',1.0)
% title('X-direction acceleration Response');
% xlabel('Time (s)');
% ylabel('acceleration  a_x (kg*cm/s^2)');
% 
% 
% % figure (12)
% % axes('FontSize',14,'FontName','TimesNewRoman');
% % plot(state1(:,1),state1(:,5),'linewidth',1.0)
% % title('X-direction acceleration Response');
% % xlabel('Time (s)');
% % ylabel('acceleration  a_x (kg*cm/s^2)');
% peakvalue1=max(abs(state2(:,5)))
% peakvalue2=max(abs(state2(:,6)))
% peakvalue3=max(abs(state2(:,8)))
% peakvalue4=max(abs(state2(:,9)))













































