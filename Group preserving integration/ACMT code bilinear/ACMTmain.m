clc;
clear all;
close all;
%% parameters
dt=0.005;
ke=100;
kp=10;
Qy=40;
ee=1e-10;
%Qy=1;
%% input data
t=0:dt:20;
Q11=(t)';
Q12=(t)';
Q21=(t)';
Q22=(t)';
Q31=(t)';
Q32=(t)';
Q=[ Q11  Q12  Q21  Q22  Q31  Q32 ];
Q(1,:)=[0 ; 0 ; 0 ; 0 ; 0 ; 0 ];
%% A matrix
beta=(ke)/(ke+kp);
qy=Qy/ke;
%% initial values
q=[0 ; 0 ; 0 ; 0 ; 0 ; 0 ]';
dq=[0 ; 0 ; 0 ; 0 ; 0 ; 0 ];
dQ=[0 ; 0 ; 0 ; 0 ; 0 ; 0 ];
qp=[0 ; 0 ; 0 ; 0 ; 0 ; 0 ];
dqp=[0 ; 0 ; 0 ; 0 ; 0 ; 0 ];
Qa=[0 ; 0 ; 0 ; 0 ; 0 ; 0 ];
Qb=[0 ; 0 ; 0 ; 0 ; 0 ; 0 ];
% straining condition
s=0;
% switching index
iflag=[0];
state1=[0 dq'];
state2=[0 q dQ' Qb' Qa' norm(Qa)-Qy qp' Q11(1) Q12(1) Q21(1) Q22(1) Q31(1) Q32(1) iflag];

%%
for ii=2:length(Q11)
    %ggg=ii;
       if iflag==[0];
        % off_phase
        [nextQ,nextq,nextdq,nextdQ,nextQb,nextQa,nextyf,nextqp]=off_phase(q,ke,kp,dt,Qy,Q,Qb,ee,ii);
       if nextyf>0
            % off to on
            [ts,nextq,nextdq,nextdQ,nextQb,nextQa,nextyf,nextqp,nextQ] = pull_back (iflag,dQ,q,ke,kp,dt,Qy,Q,Qb,ee,ii);
           
           Q1=Q(ii-1,:)+dQ*ts; q=nextq; dq=nextdq; dQ=nextdQ; Qb=nextQb; Qa=Q1; qp=nextqp; nextyf=0;
            state1=[state1;(ii-2)*dt+ts dq ];
            state2=[state2; [(ii-2)*dt+ts q dQ Qb' Qa nextyf qp' Q1 iflag]];
          
            iflag=[1];
            [nextQ,nextq,nextdq,nextdQ,nextQb,nextQa,nextyf,nextqp,nextdqp]=on_phase(q,ke,kp,dt-ts,Qy,Q,Qb,qp,Qa,qy,beta,ee,ii);
       else
           iflag=[0];
       end
       elseif iflag==[1];
             s=(Qa*dQ')/(kp*Qy);
             if s>0 
             iflag=[1];
            [nextQ,nextq,nextdq,nextdQ,nextQb,nextQa,nextyf,nextqp,nextdqp]=on_phase(q,ke,kp,dt,Qy,Q,Qb,qp,Qa,qy,beta,ee,ii); 
               else 
                   iflag=[0];

                     [nextQ,nextq,nextdq,nextdQ,nextQb,nextQa,nextyf,nextqp]=off_phase(q,ke,kp,dt,Qy,Q,Qb,ee,ii);
            
                   if nextyf>0
            % off to on
                    [ts,nextq,nextdq,nextdQ,nextQb,nextQa,nextyf,nextqp,nextQ] = pull_back (iflag,dQ,q,ke,kp,dt,Qy,Q,Qb,ee,ii);
                     Q1=Q(ii-1,:)+dQ*ts;q=nextq; dq=nextdq; dQ=nextdQ; Qb=nextQb; Qa=Q1; qp=nextqp;
                    state1=[state1;(ii-2)*dt+ts dq ];
                    state2=[state2; [(ii-2)*dt+ts q dQ Qb' Qa nextyf qp' Q1 iflag]];
                   
                    iflag=[1];
                    [nextq,nextdq,nextdQ,nextQb,nextQa,nextyf,nextqp,nextdqp]=on_phase(q,ke,kp,dt-ts,Qy,Q,Qb,qp,Qa,qy,beta,ee,ii);
      
           
                   else
                    iflag=[0];
                   end

                   
               end
       end
    Q1=nextQ;
    q=nextq;
    dq=nextdq;
    dQ=nextdQ;
    Qb=nextQb;
    Qa=nextQa;
    qp=nextqp;
    state1=[state1;(ii-1)*dt dq];
    state2=[state2;[(ii-1)*dt q dQ Qb' Qa nextyf qp' Q1 iflag]];

end

figure(1)
axes('FontSize',14,'FontName','TimesNewRoman');
plot(state2(:,2),state2(:,33),'linewidth',1.0)
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













































