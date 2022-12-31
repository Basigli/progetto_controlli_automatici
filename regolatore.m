clear; close all;clc;

% SPECIFICHE
% 
% -) Errore a regime in risposta a un gradino w(t) = 20*1(t) e d(t) = 20*1(t) pari a 0.01
%
% -) Attenuazione di almeno 50dB per d(t)
%       con [omega_d_min, omega_d_MAX] = [0,0.1]
%
% -) Attenuazione di almeno 35dB per n(t)
%       con [omega_n_min, omega_n_MAX] = [10^3,10^6]
%
% -) S% <= 5%
% -) Ta,1 <= 0.075 s
% -) Mf >= 45°

%%Specifiche

% ampiezze gradini
WW = 20;
DD = 0.1;
NN = 0.6;

% errore a regime

% attenuazione disturbo sull'uscita
omega_d_MAX = 0.1;
A_d = 50;

% attenuazione disturbo di misura
omega_n_min = 10^3;
omega_n_MAX = 10^6;
A_n = 35;

% Sovraelongazione massima e tempo d'assestamento al 5%
S_p = 5;
T_a5_spec = 0.075;

%parametri
bb = 50;
g = 9.8;
mi = 5;
ei = 0.1;
ie = 50;

%% Creazione sistema
A = [0 1; -g*mi*ei*sqrt(3)/(2*(mi*ei*ei+ie)) -bb/(mi*ei*ei+ie)];
B = [0; 1/(mi*ei*ei+ie)];
C = [0 1];
D = 0;

sys = ss(A,B,C,D);

GG = tf(sys);
s = tf('s');

%% Diagramma di Bode

figure(1);
bode(GG);
zoom on;
%% Regolatore statico 

RR_s =  0.3 * 10^4/s;
GG_e = RR_s * GG;
GG_e = minreal(GG_e);
%% Diagrammi di Bode di Ge con specifiche
% [Gm, Mf_curr, Wcg, Wcp]=margin(GG_e)
figure(2);
hold on;

%scrivere calcolo specifiche S% => Margine di fase
xi = log(S_p)/(sqrt(log(S_p)^2+pi*pi)); %smorzamento
Mf_spec = xi * 100; 

% Specifiche su d
omega_d_min = 0.0001; % lower bound per il plot
Bnd_d_x = [omega_d_min; omega_d_MAX; omega_d_MAX; omega_d_min];
Bnd_d_y = [A_d; A_d; -150; -150];
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche su n
Bnd_n_x = [omega_n_min; omega_n_MAX; omega_n_MAX; omega_n_min];
Bnd_n_y = [-A_n; -A_n; 100; 100];
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche tempo d'assestamento (minima pulsazione di taglio)
omega_Ta_min = 1e-4; % lower bound per il plot
omega_Ta_MAX = 300/(Mf_spec * T_a5_spec); % omega_c >= 300/(Mf*T^*) ~ 4.6
Bnd_Ta_x = [omega_Ta_min; omega_Ta_MAX; omega_Ta_MAX; omega_Ta_min];
Bnd_Ta_y = [0; 0; -150; -150];
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);

% Legenda colori
Legend_mag = ["A_d", "A_n", "\omega_{c,min}", "G(j\omega)"];
legend(Legend_mag);

% Plot Bode con margini di stabilità
margin(GG_e);
grid on; zoom on

% Specifiche sovraelongazione (margine di fase)
omega_c_min = omega_Ta_MAX;
omega_c_MAX = 10^3;

phi_spec = Mf_spec - 180;
phi_low = -270; % lower bound per il plot

Bnd_Mf_x = [omega_c_min; omega_c_MAX; omega_c_MAX; omega_c_min];
Bnd_Mf_y = [phi_spec; phi_spec; phi_low; phi_low];
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

% Legenda colori
Legend_arg = ["G(j\omega)"; "M_f"];
legend(Legend_arg);

%% Design del regolatore dinamico

w_c_star = 60;
M_f_star = 85 * pi/180;
M_star = 1/abs(evalfr(GG_e, i *w_c_star));
phi_star = M_f_star-pi-angle(evalfr(GG_e,i * w_c_star));

tau = (M_star-cos(phi_star))/(w_c_star*sin(phi_star));
alpha = (cos(phi_star)-1/M_star)/(w_c_star*sin(phi_star)*tau);
R_d = ((1+tau*s))/((1+alpha*tau*s) * (1+s*1/(5.5*10^2)));
% R_d = (1+tau*s)/((1+alpha*tau*s)*(1+s*1/(7*10^2))^2*(1+s*1/(7.5*10^2)));

%%funzione ad anello
LL = R_d * GG_e;
figure(3);
hold on;

% Specifiche di ampiezza
[Gm, Mf_curr, Wcg, Wcp] = margin(LL);
omega_Ta_MAX = 300/(Mf_curr * T_a5_spec);
Bnd_Ta_x = [omega_Ta_min; omega_Ta_MAX; omega_Ta_MAX; omega_Ta_min];

patch(Bnd_d_x, Bnd_d_y, 'r', 'FaceAlpha', 0.2, 'EdgeAlpha', 0);
patch(Bnd_n_x, Bnd_n_y, 'g', 'FaceAlpha', 0.2, 'EdgeAlpha', 0);
patch(Bnd_Ta_x, Bnd_Ta_y, 'b', 'FaceAlpha', 0.2, 'EdgeAlpha', 0);
legend(Legend_mag);

% Plot Bode con margini di stabilità
margin(LL);
grid on; zoom on;

%Specifiche su fase
omega_c_min = omega_Ta_MAX;
Bnd_Mf_x = [omega_c_min; omega_c_MAX; omega_c_MAX; omega_c_min];

patch(Bnd_Mf_x, Bnd_Mf_y, 'g', 'FaceAlpha', 0.2, 'EdgeAlpha', 0);
hold on;
legend(Legend_arg);

%% Check prestazioni in anello chiuso

% Funzione di sensitività complementare
FF = LL/(1+LL);
FF = minreal(FF);
FF
figure(4);
hold on;
bode(FF);
margin(LL);
poles = pole(FF);
fprintf('Poli dei F(s):\n');
disp(poles); 

% Risposta ad un gradino

figure(5);
T_simulation = 0.5;
[y_step, t_step] = step(WW*FF, T_simulation);
plot(t_step, y_step, 'b');
grid on; zoom on; hold on;

% vincolo di sovraelongazione
S_100_spec = 5;
patch([0, T_simulation, T_simulation, 0], [WW*(1+S_100_spec), WW*(1+S_100_spec), WW+1, WW+1], 'r', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.5);
ylim([0,WW+1]);

% vincolo tempo di assestamento al 5%
LV = abs(evalfr(WW*FF, 0)); %valore limite gradino: W*F(0)
patch([T_a5_spec, T_simulation, T_simulation, T_a5_spec], [LV*(1-0.05), LV*(1-0.05),0,0], 'g', 'FaceAlpha', 0.1, 'EdgeAlpha',0.5);
patch([T_a5_spec, T_simulation, T_simulation, T_a5_spec], [LV*(1+0.05), LV*(1+0.05), LV+1, LV+1], 'g', 'FaceAlpha', 0.1, 'EdgeAlpha',0.1);

Legend_step=["Risposta al gradino"; "Vincolo sovraelongazione"; "Vincolo tempo di assestamento"];
legend(Legend_step);

%% Check disturbo in uscita

% Funzione di sensitività
SS=1/(1+LL);
figure(7);

% Simulazione del disturbo a pulsazione 0.025 (sono 4 a multipli di 0.025)
omega_d = 0.0025;
tt = (0:1e-2: 5*1e3);
dd = DD*(sin(omega_d*tt) + sin(omega_d*2*tt) + sin(omega_d*3*tt) + sin(omega_d*4*tt));
y_d = lsim(SS, dd, tt);
d_attenuato = dd / 10^(50/20); 
hold on, grid on, zoom on;
plot(tt, dd,'m');
plot(tt,y_d, 'b');
plot(tt,d_attenuato, 'g');
grid on
legend('disturbo d(t)', 'uscita del sistema', 'attenuazione minima');

%Simulazione dell'errore di misura n 
figure(8);
omega_n = 10^3;
tt = (0:1e-4:1e-1);
nn = NN*(sin(omega_n*tt) + sin(omega_n*2*tt) + sin(omega_n*3*tt) + sin(omega_n*4*tt));
n_attenuato = nn / 10^(35/20); 
y_n = lsim(FF, nn, tt);
hold on, grid on, zoom on;
plot(tt, nn, 'm');
plot(tt, y_n, 'b');
plot(tt, n_attenuato, 'g'); 
grid on
legend('disturbo n(t)', 'uscita del sistema', 'attenuazione minima');