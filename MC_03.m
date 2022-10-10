function ris = MC_03(alpha)

h = 0.015;
a1 = 145; a2 = 225; a3 = 235; a4 = 295; a5 = 360;
in = [0 a1 a2 a3 a4 a5];
xv_1.v = 0.5;
xv_1.w = 0.9;

if alpha>= in(1) && alpha < in(2)   %rise
    da = in(2) - in(1);
    alpha_da=(alpha-in(1))/da;

    ris.pos = 0;
    ris.vel = 0;
    ris.acc = 0;
elseif alpha>=in(2) && alpha < in(3)
    da = (in(3)-in(2));
    alpha_ad = (alpha-in(2))/da;
    out = MCMtraposoidal(alpha_ad);
    %out = MCM_sshape(alpha_ad,xv_1);
    ris.pos = h*out.pos;
    ris.vel = h/deg2rad(da)*out.vel;
    ris.acc = h/(deg2rad(da))^2*out.acc;
elseif alpha>=in(3) && alpha < in(4)
    da = (in(4)-in(3));
    alpha_ad=(alpha-in(3))/da;
    ris.pos = h;
    ris.vel = 0;
    ris.acc = 0;
elseif alpha>=in(4) && alpha < in(5)
    da = (in(5)-in(4));
    alpha_ad=(alpha-in(4))/da;
    out = MCMtraposoidal(alpha_ad);
    %out = MCM_sshape(alpha_ad,xv_1);
    ris.pos = h- h*out.pos;
    ris.vel = -h/deg2rad(da)*out.vel;
    ris.acc = -h/(deg2rad(da))^2*out.acc;

elseif alpha>=in(5) && alpha < in(6)
    da = (in(6)-in(5));
    alpha_ad=(alpha-in(5))/da;    
    ris.pos=0;
    ris.vel=0;
    ris.acc=0;
else
    da = 0;
    ris.pos=0;
    ris.vel=0;
    ris.acc=0;    
end
end