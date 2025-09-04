%% Linealizacion
% Defino el orden del sistema (tamaño del vector de estados)

orden = 3;
x = sym('x',[orden 1],'real');
u = sym('u','real');

% Defino las ecuaciones de estado (la cantidad de estas ecuaciones son el
% orden del sistema)
m = 0.1;
R = 4;
g = 10;
L = 0.01;

% Variables de estado
% x1 = y (posicion vertical de la bolita)
% x2 = y´ (velocidad vertical de la bolita)
% x3 = i (intensidad de corriente)
% u = v (tension de la fuente)

F1 = x(2);
F2 = g - x(3) / (m * x(1));
F3 = u/L - R * x(3)/L;
H1 = x(1);

% Linealizo mediante el metodo del Jacobiano alrededor del punto de
% equilibrio hallado

f = [F1; F2; F3];
h = [H1];

A = jacobian(f,x);
B = jacobian(f,u);
C = jacobian(h,x);
D = jacobian(h,u);
%%
x0 = [1; 0; m*g];
u0 = m*g*R;

% Evaluo las matrices en el equilibrio
A=subs(A, {'x1','x2', 'x3', 'u'}, {x0(1), x0(2), x0(3), u0});
B=subs(B, {'x1','x2', 'x3', 'u'}, {x0(1), x0(2), x0(3), u0});
C=subs(C, {'x1','x2', 'x3', 'u'}, {x0(1), x0(2), x0(3), u0});
D=subs(D, {'x1','x2', 'x3', 'u'}, {x0(1), x0(2), x0(3), u0});

A=double(A);
B=double(B);
C=double(C);
D=double(D);

P0=zpk(ss(A,B,C,D));

%% Configuracion del Bode 
optionss=bodeoptions;
optionss.MagVisible='on';
optionss.PhaseMatching='on';
optionss.PhaseMatchingValue=-180;
optionss.PhaseMatchingFreq=1;
optionss.Grid='on';
%optionss.Xlim={[1, 1e4]};

bode(P0, optionss)
%% Compensacion

% Tengo la planta P0, la divido en parte pasa todo (Pap) y de fase
% minima (Pmp)

s = tf('s');
p1 = 3.162;
p2 = 400;
Pap = zpk(-p1, p1, 1);
Pmp = zpk(minreal(P0 / Pap));

%%
figure();bode(Pap, optionss); legend('Pap'); set(findall(gcf,'type','line'),'linewidth',2);
figure();bode(Pmp, optionss); legend('Pmp'); set(findall(gcf,'type','line'),'linewidth',2);

%%

% Busco Ts de tal forma que la digitalización me reste 5º de fase. Para
% ello busco la frecuencia de corte wgc para la cual el resto de la parte
% pasa todo resta 25º. 

%fase_ap = 2*m*atan(wgc/z) -> wgc <= z*tan(fase_ap/(2*m))
%fase_ap = 2*m*atan(p/wgc) -> wgc >= p/tan(fase_ap/(2*m))
% m es la multiplicidad de la singularidad

fase_ap = deg2rad(25);
m = 1;
wgc = p1/tan(fase_ap/(2*m));
%%
%fase_ap = 2*m*atan(p/wgc) + 2*m*atan(wgc/z)
%wgc = sqrt(p*z)


% Despejo Ts de tal forma que reste 5º.
fase_dig = deg2rad(5);
%fase_dig = 2*atan(wgc*Ts/4) -> Ts = 4*tan(fase_dig/2)/wgc
Ts = 4*tan(fase_dig/2)/wgc;

%%
% La digitalizacion del controlador se puede modelar como un retardo de
% media muestra, que se puede aproximar con el aproximante de Pade de
% primer orden
%OJO NUNCA PONER EL PADE CON ZPK PORQUE TE RETRASA 180°
Pade = (1-s*Ts/4)/(1+s*Ts/4);
P = P0 * Pade;
close all;
figure();bode(P, optionss); legend('P'); set(findall(gcf,'type','line'),'linewidth',2);
%%

% Propongo el compensador y ajusto kc de tal forma de obtener el wgc
% hallado anteriormente

kc = db2mag(104);
Control = -zpk([-p1 -p1], [0 -1e4], kc);
L = minreal(Control*P);
optionss.PhaseMatching='on';
optionss.PhaseMatchingValue=-180;
optionss.PhaseMatchingFreq=1;
close all; figure();bode(L, optionss); legend('C * P (con kc)'); set(findall(gcf,'type','line'),'linewidth',2);
figure; nyqlog(L)


%%
% Digitalizo el controlador
Cdig = c2d(Control, Ts, 'Tustin');
Pdig = c2d(P, Ts, 'Tustin');
L = minreal(Cdig*Pdig);


%% Grupo de los 4
S=1/(1+L);              % funcion de sensibilidad
T=1-S;                  % transferencia a lazo cerrado
% Bodes
freqrange={1,1e5};
figure();
h1=subplot(2,1,1); bode(L,optionss,freqrange);title('Ld')
optionss.PhaseVisible='off';
h2=subplot(2,1,2); bode(S,T,optionss,freqrange);title('S & T')



%% Respuesta al escalon
figure();
time=10;
step(T,time);title('S & T')



%% Margen de estabilidad (Sm)
% Para calcular el margen de estabilidad, me fijo en el bode de la funcion
% de sensibilidad S. El pico maximo So de dicha funcion (habiendolo pasado de
% dB a veces) es la inversa del Sm. Es decir:
% 20 * log(1/Sm) = So -> 1/Sm = 10^(So/20) -> Sm = 10^(-So/20)
bode(S, optionss)

So = 2.68;
Sm = 10^(-So/20);           % margen de estabilidad

%% Realimentacion en Espacio de Estados
% Defino el sistema de espacio de estados
sys = ss(A, B, C, D);
% Polos de lazo abierto
E = eig(A);

% Polos que quiero
P = [-1000 -3.1623 -400];
% Solucion de matriz K para los polos que quiero
K = place(A, B, P);
% Chequeo
A_cl = A - B*K;
E_cl = eig(A_cl);

% Defino el nuevo sistema
sys_cl = ss(A_cl, B, C, D);
% Respuesta al escalon
step(sys_cl);

% Ajusto kr para obtener la unidad a la salida
k_dc = dcgain(sys_cl);
kr = 1/k_dc;
sys_cl_kr = ss(A_cl, kr * B, C, D);
step(sys_cl_kr);