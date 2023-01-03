clear ; close all; clc



%parametri
bb=50;
g=9.8;
mi=5;
ei=0.1;
ie=50;

x0=[7*pi/6;0];
interv=0:0.2:100;

width=0.1;
heigth=1;
radius= .1;

%% risoluzione equazione differenziale
uu= @(t) 0;
% dinamica: f(t, x)
dyn = @(t, x) [x(2); (-bb*x(2)-g*mi*ei*sin(x(1))+uu(t))/(mi*ei*ei+ie)];

% risolviamo l'equazione differenziale
[time, traj] = ode45(dyn, interv, x0);


%% animazione
xx=traj';

%creo 8 pale che aggiungo al rotore
for i=0:7
    pala=antenna.Polygon('Vertices',[-width/2 0 0; width/2 0 0; ...
        width/2 0.8 0; 0 heigth 0; -width/2 heigth 0]);
    pala=rotate(pala, 45*i, [0 0 0], [0 0 1]);
    if i==0
        rotore=pala;
    else
        rotore=add(rotore, pala);
    end
    
    
end

%imposto la condizione iniziale
rotore=rotate(rotore, x0(1)*180/pi, [0 0 0], [0 0 -1]); 

for tt=1:length(interv)
    if tt==1
        x_curr=x0(1)*180/pi;
    end
   
    %ruoto il rotore basandomi sulla soluzione dell'eq. diff. trovata
    rotore=rotate(rotore, x_curr-(xx(1,tt)*180/pi), [0 0 0], [0 0 0]);
    x_curr=xx(1,tt)*180/pi;
    pl=rotore.plot;
    
    patch(pl.XData, pl.YData, 'black', 'FaceColor','#92B4A7', 'FaceAlpha', 0.4);

    %cerchio utile per visualizzare l'andamento della posizione angolare
    %nel tempo
    rectangle('Position',[-0.8*sin(xx(1,tt))-radius/2 -0.8*cos(xx(1,tt))-radius/2 0.1 0.1], ...
        'Curvature',[1 1], 'FaceColor','#8C8A93');
    rectangle('Position', [-width -width 2*width 2*width], 'Curvature', [0.8 0.8],'FaceColor','#81667A');
    axis ([-1.1 1.1 -1.1 1.1]);

    %pause(0.02);
    anim(i)=getframe;
    if ishghandle(pl)~=1
        break;
    end

end 
