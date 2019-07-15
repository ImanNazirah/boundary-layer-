
function MHD
format long g
global M Pr Sc St s
M=0.1;
Pr=0.7;
Sc=0.6;
St=0.2;
s=1;
 %beta
%==========INITIAL CONDITION=======
solinit = bvpinit(linspace(0,8,16),@MHD_init);
%==========SOLVE USING BVP4C=======
options = bvpset('stats','on','RelTol',1e-10);
sol = bvp4c(@MHD_ode,@MHD_bc,solinit,options);
x=linspace(0,8,16);
y=deval(sol,x);
%Result will be display at excel file, where x,f,f',f'',theta,theta'
xlswrite('ProjectFluid2019.xlsx',[(sol.x).',(sol.y(1,:)).',(sol.y(2,:)).',(sol.y(3,:)).',(sol.y(4,:)).',(sol.y(5,:)).',(sol.y(6,:)).',(sol.y(7,:)).']);
%==========OUTPUT==================
figure(1)
lines={'blue','red','black'};
plot(sol.x,sol.y(2,:),lines{1})
xlabel('\eta')
ylabel('f ''(\eta)')
title('Velocity profiles for different value of M')
hold on
fprintf('\nFirst solution:\n');
fprintf('f"(0)=%0.8f\n',y(3));
hold on
for i=2:3
  M = M*2;
  sol = bvp4c(@MHD_ode,@MHD_bc,solinit,options);
  %fprintf('For M = %5i, A = %4.2f.\n',M,sol.parameters);
  plot(sol.x,sol.y(2,:),lines{i});
  drawnow
end
legend('M =    0.1','M =   0.2','M = 0.4');
hold off
 
figure(2)
lines={'blue','red','black'};
plot(sol.x,sol.y(4,:),lines{1})
xlabel('\eta')
ylabel('\theta(\eta)')
title('Temperature profiles for different value of M')
hold on
fprintf('\nFirst solution:\n');
fprintf('\theta''(0)=%0.8f\n',y(5));
hold on
M=0.1;
for i=2:3
  M = M*2;
  sol = bvp4c(@MHD_ode,@MHD_bc,solinit,options);
  %fprintf('For M = %5i, A = %4.2f.\n',M,sol.parameters);
  plot(sol.x,sol.y(4,:),lines{i});
  drawnow
end
legend('M =    0.1','M =   0.2','M = 0.4');
hold off
 
figure(3)
lines={'blue','red','black'};
plot(sol.x,sol.y(1,:),lines{1})
xlabel('\eta')
ylabel('f(\eta)')
legend('M=1')
title('f for different value of M')
hold on
fprintf('\nFirst solution:\n');
fprintf('\theta''(0)=%0.8f\n',y(2));
hold on
M=0.1;
for i=2:3
  M = M*2;
  sol = bvp4c(@MHD_ode,@MHD_bc,solinit,options);
  %fprintf('For M = %5i, A = %4.2f.\n',M,sol.parameters);
  plot(sol.x,sol.y(1,:),lines{i});
  drawnow
end
legend('M =    0.1','M =   0.2','M = 0.4');
hold off

figure(4)
lines={'blue','red','black'};
plot(sol.x,sol.y(6,:),lines{1})
xlabel('\eta')
ylabel('\phi(\eta)')
legend('M=1')
title('Concentration profile for different value of M')
hold on
fprintf('\nFirst solution:\n');
fprintf('\tphi(0)=%0.8f\n',y(6));
hold on
M=0.1;
for i=2:3
  M = M*2;
  sol = bvp4c(@MHD_ode,@MHD_bc,solinit,options);
  %fprintf('For M = %5i, A = %4.2f.\n',M,sol.parameters);
  plot(sol.x,sol.y(6,:),lines{i});
  drawnow
end
legend('M =    0.1','M =   0.2','M = 0.4');
hold off
%============INPUT ODEs===================
function dydx = MHD_ode(x,y,M,Pr,Sc,St)
global  M Pr Sc St
dydx = [ y(2)
         y(3)
         -y(1)*y(3)+2*y(2)^2+M*y(2)
         y(5)
         -Pr*(y(1)*y(5)-y(2)*y(4)-St*y(2))
         y(7)
         Sc*(-y(1)*y(7)+y(2)*y(6))];
end
%===========INPUT BCs=====================
function res = MHD_bc(ya,yb,St,Sc,s)
global St Sc s
res = [ya(1)-s
       ya(2)-1
       yb(3) 
       ya(4)-1+St
       yb(4)
       ya(6)-1
       yb(6)
       ];end
%===========INPUT INITIAL VALUEs============
function v = MHD_init(x,St,Sc,s)
global St Sc s
v = [0
0
0
0
0
0
0];end
end
