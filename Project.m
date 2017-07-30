%% Projeto Final de ES927 - 29/05/2015
% Grupo:
%   Alejo Rayo Alarcon              RA: 161267
%   Iago Carvalho                   RA: 091543
%   Rafael Fernandes Ayres Fonseca	RA: 103844
%   Romeu Andrade                   RA: 094425
%   Thiago Ribeiro                  RA: 104202


% Neste arquivo serao efetuados calculos relativos as:
%   - Definicao do robo Kuka KR 6 R900, usando toolbox rcvtools
%   - Definicao das equacoes de movimento do sistema.
%   - Preparacao para controle PD do robo

%% Definicao do robo Kuka KR 6 R900, visto pela toolbox rcvtools

% Definicao de thetas simbolicos
theta1 = sym('theta1');
theta2 = sym('theta2');
theta3 = sym('theta3');
theta4 = sym('theta4');
theta5 = sym('theta5');
theta6 = sym('theta6');

% Definicao dos links do sistema, escolha dos eixos seguindo o padrao
%  proposto por Denavit Hartenberg

% Definicao com 6 links - FUNCIONA
L(1) = Link([0 4 0.25 pi/2]);
L(2) = Link([0 0 4.55 0]);
L(3) = Link([0 0 0.35 -pi/2]);
L(4) = Link([0 4.20 0 pi/2]);
L(5) = Link([0 0 0 -pi/2]);
L(6) = Link([0 1 0 0]);

%% Matrizes de rotacao
T01 = trotz(theta1)*transl(0, 0, 4)*transl(0.25,0,0)*round(trotx(pi/2));
T12 = trotz(theta2)*transl(0, 0, 0)*transl(4.55,0,0)*round(trotx(0));
T02 = T01 * T12;
T23 = trotz(theta3)*transl(0, 0, 0)*transl(0.35,0,0)*round(trotx(-pi/2));
T03 = T02 * T23;
T34 = trotz(theta4)*transl(0, 0, 4.2)*transl(0,0,0)*round(trotx(pi/2));
T04 = T03 * T34;
T45 = trotz(theta5)*transl(0, 0, 0)*transl(0,0,0)*round(trotx(-pi/2));
T05 = T04 * T45;
T56 = trotz(theta6)*transl(0, 0, 1)*transl(0,0,0)*round(trotx(0));
T06 = T05 * T56;

% Referenciando os links em Objeto usado do toolbox
kukaKR6 = SerialLink(L); 
% Ilustracao do robo em sua configuracao padrao
kukaKR6.name = 'Kuka KR 6 R900'; 
% figure(1);
% 6 thetas
% ThreeLink.plot([-pi pi/2 0 -pi/2 pi/3 0]);


%% Cinematica direta do robo, tendo como base configuracao padrao do robo
% 6 Thetas
TL = kukaKR6.fkine([-pi pi/2 0 -pi/2 pi/3 0]);

%% Cinematica inversa, tendo como ponto final o resultado da cinematica
% direta calculado acima
QF = kukaKR6.ikunc(TL);

% Ilustracao da configuracao obtida pela cinematica inversa 
% figure(2)
% ThreeLink.plot(QF);

%% Definicao da equacao de movimento do sistema
% Equacao (forma matricial): D(q)q''+ C(q, q')q' + g(q) = t
%   - q, q', q'': Coordenadas generalizadas, aqui sendo os 6 thetas de 
%       nosso sistema com suas reespectivas derivas;
%   - D: Matriz de inercia do sistema;
%   - C: Matriz baseada nos simbolos de Christoffel;
%   - g: Derivada da Energia potencial em funcao de cada theta. No calculo
%       futuro do controle PD tal matriz sera igual a zero;

%% Calculo dos jacobianos 
%thetas = [theta1 theta2 theta3 theta4 theta5 theta6 theta7];
thetas = [-pi pi/2 0 -pi/2 pi/3 0];
% J6 = ThreeLink.jacob0(thetas);
% J6 = tr2jac(T06);
% J1 = J6;
% J2 = J6;
% J3 = J6;
% J4 = J6;
% J5 = J6;
% %J6 = ThreeLink.jacob0(thetas);
% J1(:,2:6) = zeros(6,1:5);
% J2(:,3:6) = zeros(6,1:4);
% J3(:,4:6) = zeros(6,1:3);
% J4(:,5:6) = zeros(6,1:2);
% J5(:,6) = zeros(6,1);

%% Massa de cada elo, em kg, total mass = 52Kg
% Aproximacoes foram feitas com o objetivo de calcular 
M1 = 15;
M2 = 17.5;
M3 = 2.5;
M4 = 12.5;
M5 = 2.5;
M6 = 2;

kukaKR6.links(1).m = M1;
kukaKR6.links(1).m = M2;
kukaKR6.links(1).m = M3;
kukaKR6.links(1).m = M4;
kukaKR6.links(1).m = M5;
kukaKR6.links(1).m = M6;

% Raio de cada elo
r1 = 7.5;
r2 = 7.5; 
r3 = 5;
r4 = 5;
r5 = 2.5;
r6 = 2;

% Altura de cada elo
h1 = 4.00;
h2 = 4.55;
h3 = 0.35;
h4 = 4.20;
h5 = 0.01;
h6 = 1.00;

% Centro de massa de cada elo [x y z]
kukaKR6.links(1).r = [0 0 h1/2];
kukaKR6.links(2).r = [0 0 h2/2];
kukaKR6.links(3).r = [0 0 h3/2];
kukaKR6.links(4).r = [0 0 h4/2];
kukaKR6.links(5).r = [0 0 h5/2];
kukaKR6.links(6).r = [0 0 h6/2];

%% Momento de inercia, aproximando os elos a cilindros.
%  Todos os os cilindros possuem seu eixo Z passando por seu centro
%  Iz = M * r^2 / 2  
%  Ix = Iy = M*(3*r^2 + h^2)/12
Iz = 0.5 * [M1 M2 M3 M4 M5 M6]' * [r1^2 r2^2 r3^2 r4^2 r5^2 r6^2];
Ixy = (1/12) * [M1 M2 M3 M4 M5 M6]' * [(3*r1^2 + h1^2) (3*r2^2 + h2^2) (3*r3^2 + h3^2) (3*r4^2 + h4^2) (3*r5^2 + h5^2) (3*r6^2 + h6^2)];

I1 = [ Ixy(1) 0 0; 0 Ixy(1) 0; 0 0 Iz(1)];
I2 = [ Ixy(2) 0 0; 0 Ixy(2) 0; 0 0 Iz(2)];
I3 = [ Ixy(3) 0 0; 0 Ixy(3) 0; 0 0 Iz(3)];
I4 = [ Ixy(4) 0 0; 0 Ixy(4) 0; 0 0 Iz(4)];
I5 = [ Ixy(5) 0 0; 0 Ixy(5) 0; 0 0 Iz(5)];
I6 = [ Ixy(6) 0 0; 0 Ixy(6) 0; 0 0 Iz(6)];

kukaKR6.links(1).I = I1;
kukaKR6.links(2).I = I2;
kukaKR6.links(3).I = I3;
kukaKR6.links(4).I = I4;
kukaKR6.links(5).I = I5;
kukaKR6.links(6).I = I6;

% ...................
%% Matriz de inercia 
% D = M1 * J1(1:3,:)' * J1(1:3,:) + J1(4:6,:)' * T01(1:3,1:3) * I1 * T01(1:3,1:3)' * J1(4:6,:);
% D = D + M2 * J2(1:3,:)' * J2(1:3,:) + J2(4:6,:)' * T02(1:3,1:3) * I2 * T02(1:3,1:3)' * J2(4:6,:); 
% D = D + M3 * J3(1:3,:)' * J3(1:3,:) + J3(4:6,:)' * T03(1:3,1:3) * I3 * T03(1:3,1:3)' * J3(4:6,:); 
% D = D + M4 * J4(1:3,:)' * J4(1:3,:) + J4(4:6,:)' * T04(1:3,1:3) * I4 * T04(1:3,1:3)' * J4(4:6,:); 
% D = D + M5 * J5(1:3,:)' * J5(1:3,:) + J5(4:6,:)' * T05(1:3,1:3) * I5 * T05(1:3,1:3)' * J5(4:6,:); 
% D = D + M6 * J6(1:3,:)' * J6(1:3,:) + J6(4:6,:)' * T06(1:3,1:3) * I6 * T06(1:3,1:3)' * J6(4:6,:); 

%% Posicao e velocidade de referencia
qf = [pi/2 4*pi/3 -pi/2 0 -pi/4 0]
qd = [pi/6 pi/2 0 0 -pi pi/2]

% Calculos do Jacobiano 
kukaKR6.jacob0(qf)
% Calculo da matriz de inercia
kukaKR6.inertia(qf)
% Calculo da matriz de Coriolis
kukaKR6.coriolis(qf,qd)
% Calculo do torque
kukaKR6.rne(qf,qd,0*qd)












