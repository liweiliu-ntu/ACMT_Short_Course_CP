clc
clear all
close all
syms c11 c22 c33 c12 c23 c13 t11 t22 t33 t12 t13 t23 Qy d11 d22 d33  d12 d13 d23 x xx yy zz
load("Propor strain loading.csv");
% E1111=30910;
% E2222=21240;
% E3333=21240;
% E1122=14350;
% E1133=14350;
% E2233=11900;
% E2323=5630;
% E1313=5630;
% E1212=4670;
% theta1=7.20721e-05;
% theta2=2.75482e-04;
% theta3=1.31273e-05;
% theta4=3.34595e-05;
% theta5=2.77469e-05;
% theta6=3.96396e-03;
% theta7=1.06061e-02;


% material constant
% Es=200e3;
% v=0.3;
% 
% % K= [ E1111 E1122 E1133 0      0      0;
% %      E1122 E2222 E2233 0      0      0;
% %      E1133 E2233 E3333 0  0      0;
% %      0     0     0     E2323  0      0;
% %      0     0     0     0      E1313  0;
% %      0     0     0     0      0      E1212];
% 
% Y= [ 2/3 -1/3 -1/3 0      0      0;
%      -1/3 2/3 -1/3 0      0      0;
%      -1/3 -1/3 2/3  0     0      0;
%      0     0     0     1   0      0;
%      0     0     0     0      1  0;
%      0     0     0     0      0      1];
% P1=Y/((1+v)*(1-2*v));
% K=Es*P1;
% % Y= [ theta1 theta4 theta4 0      0      0      ;
% %      theta4 theta2 theta5 0      0      0      ;
% %      theta4 theta5 theta2 0      0      0      ;
% %      0     0     0     theta3  0      0      ;
% %      0     0     0     0      theta3  0      ;
% %      0     0     0     0      0      0.5*theta1-0.5*theta4  ];
% 
% % b=[ theta6 ;theta7 ; theta7   ;    0    ;  0     ;  0 ];
% C=1/Es*[ 1    -v     -v    0      0      0;
%         -v     1     -v    0      0      0;
%         -v     -v     1    0      0      0;
%          0     0     0      2*(1+v)      0      0;
%          0     0     0     0       2*(1+v)     0;
%          0     0     0     0      0       2*(1+v)];
% b=[ 0;0 ; 0   ;    0    ;  0     ;  0 ];
% material constant
Es=200e3;
v=0.3;
tau_y=sqrt(0.5)*150;
Et=0e3;

% operator matrix
P=zeros(6,6);
Pd=P;
Pd(1:3,1:3)=-1/3*ones(3,3);
Pd=Pd+eye(6);
P(1:3,1:3)=v*ones(3,3)-(2*v-1)*eye(3);
P(4:6,4:6)=(1-2*v)/2*eye(3);
P=P/((1+v)*(1-2*v));

% stiffness modulus matreix
K=Es*P;
%Kp=Et*P;
YY=   [ 1    -v     -v    0      0      0;
        -v     1     -v    0      0      0;
        -v     -v     1    0      0      0;
         0     0     0      2*(1+v)      0      0;
         0     0     0     0       2*(1+v)     0;
         0     0     0     0      0       2*(1+v)];
C=1/Es*YY
% yield parameter
Y=1/(tau_y)^2*YY;
b=zeros(6,1);

q=transpose(Propor_strain_loading(:,2:7));
t=transpose(Propor_strain_loading(:,1));
Q = zeros(size(q));
ff=[];
dt=0.01;
%intial stress
Q(:,1)=[0;0;0;0;0;0];
%yield function
r=1;
dd=[d11;d22;d33;0;0;0];

yf=0.5.*dd'*Y*dd+b'*dd.*r-0.5.*r^2;
err=10e-8;
Iflag(1)=[0];
for i = 1:length(q)-1
    q_dot(:,i)=(q(:,i+1)-q(:,i))/dt;
    Q_dot(:,i)=K*(q(:,i+1)-q(:,i))/dt;
    u(:,i)=[(q(:,i+1)-q(:,i))./dt; 0];
end
tStart = cputime;
for i = 1:length(q)-1
     ff(i)=0.5.*(Q(:,i)'*Y*Q(:,i))+(b'*Q(:,i)).*r-0.5.*r.^2;
     Qnext = Q(:,i) + K*(q(:,i+1)-q(:,i));
     
     if Iflag==[0]  

      [fnext] = Yield_function(Qnext,Y,b,r);
          if fnext>err 
              Iflag=[1];
              iflag(i)=1;
           %elastic to plastic

           Q_off=Q(:,i); q_off=q(:,i);tton=t(i);Q_dot1=Q_dot(:,i);q_dot1=q_dot(:,i);
           %time_on
            [t_on,Qt,qt,d_ton] = ton(Q_dot1,Q_off,q_dot1,r,tton,Y,b, q_off);
            uu=u(:,i);
            Qnow=Qt;
            qnow=qt;
            qnext=q(:,i+1);

            [QQnext] = return_mapping(Qnow,Y,K,qnext,b,r,qnow,C);

             Q(:,i+1)=QQnext;
             Qnext=QQnext;
             [fnext] = Yield_function(Qnext,Y,b,r);
             f(i)=fnext;
          else  
iflag(i)=2;
          %elastic to elastic
          Iflag=[0];

          Q(:,i+1)=Qnext;
          [fnext] = Yield_function(Qnext,Y,b,r);
          f(i)=fnext;
       
          end

     elseif Iflag==[1] 
          q_dot1=q_dot(:,i);
          Qnow=Q(:,i);
          [Sg] = Straining_condition(Y,K,r,b,q_dot1,Qnow);
           if Sg<=0
         %plastic to elastic
iflag(i)=3;
          Q(:,i+1)=Qnext;
          [fnext] = Yield_function(Qnext,Y,b,r);
          f(i)=fnext;
          Iflag=[0];
         else
iflag(i)=4;  
          %plastic to plastic

             Qnow=Q(:,i);
             qnow=q(:,i);
             qnext=q(:,i+1);
             [QQnext] = return_mapping(Qnow,Y,K,qnext,b,r,qnow,C);
              Q(:,i+1)=QQnext;
              Qnext=QQnext;
             [fnext] = Yield_function(Qnext,Y,b,r);
             f(i)=fnext;
              Iflag=[1];
     end
   
     

      
         
     end

% f(i)=0.5.*Q(:,i)'*Y*Q(:,i)+b'*Q(:,i).*r-0.5.*r^2 ;
end

tEnd = cputime - tStart


figure(1)
 plot(q(1,:),Q(1,:))
 xlabel('Strain q1 ');ylabel('Stress Q1 (Mpa)')
 
figure(2)
  plot(t(1,13:length(f)),f(1,13:length(f)))
% grid on
%  xlabel('Stress Q1 (Mpa)');ylabel('Stress Q2 (Mpa)');zlabel('Stress Q3 (Mpa)')
%  yfz=solve(yf,d33);
%     yfz1=yfz(1,1);
%     yfz2=yfz(2,1);
% yfz11=subs(yfz1,[d12,d13,d23,d11,d22],[Q(4,1),Q(5,1),Q(6,1),xx,yy]);
% yfz12=subs(yfz2,[d12,d13,d23,d11,d22],[Q(4,1),Q(5,1),Q(6,1),xx,yy]);

% 
% figure(5)
%   plot(t(1,8:length(f)),f(1,8:length(f)))
% 
 figure(6)
plot3(Q(1,:),Q(2,:),Q(3,:))
grid on
 xlabel('Stress Q1 (Mpa)');ylabel('Stress Q2 (Mpa)');zlabel('Stress Q3 (Mpa)')

% 
% figure(7)
% pp=fsurf(yfz11,[-5000000 5000000 -5000000 5000000])
% pp.FaceAlpha = 0.2;
% hold on
% ppp=fsurf(yfz12,[-5000000 5000000 -5000000 5000000] )
% ppp.FaceAlpha = 0.2;
% hold on
% pppp=plot3(Q(1,:),Q(2,:),Q(3,:))
% pppp.LineWidth = 2
% grid on
%  xlabel('Stress Q1 (Mpa)');ylabel('Stress Q2 (Mpa)');zlabel('Stress Q3(Mpa)')