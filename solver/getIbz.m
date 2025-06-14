function [ptx, pty] = getIbz (A1, A2)
        D1 = [A1 A2];
        B1 = inv(D1);
        piB = pi*B1;
        pib1 = piB(1,:); %pt 1
        pib2 = piB(2,:);
        
        %%%pt 2 calculation
        pt1 = pib1; %pt1
        r1 = norm(pt1);
        r2 = r1/ tan(pi/3);
        r3 = r1/sin(pi/3);
        
        p1 = (r3^2 - r2^2 + (pt1(1))^2 + (pt1(2))^2)/(2*pt1(1));
        q1 = -(pt1(2))/pt1(1);
        x = (1+q1*q1);
        y = 2*p1*q1;
        z = p1*p1-r3*r3;
        b = roots([x y z]);
        a1 = p1 + q1*b;
        pt2 = [a1 b];%pt 2 %each row one coordinate

        ptx = [pt1(1) pt2(1,1)];
        pty = [pt1(2) pt2(1,2)];
    end