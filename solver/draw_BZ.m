function draw_BZ(a1,a2)

D =  [a1' a2'];
B = D^-1;

b1 = pi*B(1,:)';
b2 = pi*B(2,:)';
b3 = b1 + b2;

[xb1,yb1] = perpendicular_coord(b1);
[xb2,yb2] = perpendicular_coord(b2);
[xb3,yb3] = perpendicular_coord(b3);
[minxb1,minyb1] = perpendicular_coord(-b1);
[minxb2,minyb2] = perpendicular_coord(-b2);
[minxb3,minyb3] = perpendicular_coord(-b3);

a = [xb3(1) yb3(1)] - [xb3(end) yb3(end)];

dot(a,b3)


plot([0 b1(1)],[0 b1(2)],"k")
hold on
plot(xb1,yb1,"r")
hold on
plot([0 b2(1)],[0 b2(2)],"k")
hold on
plot(xb2,yb2,"r")
hold on
plot([0 b3(1)],[0 b3(2)],"k")
hold on
plot(xb3,yb3,"r")
hold on
plot([0 -b1(1)],[0 -b1(2)],"k")
hold on
plot(minxb1,minyb1,"r")
hold on
plot([0 -b2(1)],[0 -b2(2)],"k")
hold on
plot(minxb2,minyb2,"r")
hold on
plot([0 -b3(1)],[0 -b3(2)],"k")
hold on
plot(minxb3,minyb3,"r")
hold on
xlabel("$k_x$","interpreter","latex","fontsize",20)
ylabel("$k_y$","interpreter","latex","fontsize",20)

% 
% print("-dpng", "-r300","../result/IBZ")
end


function [x,y]  = perpendicular_coord(point)

if point(1) ~= 0
x = (point(1)-2*pi): 0.0001:(point(1)+2*pi);
m = point(2)/point(1);

for p = 1:length(x)
    y(p) = (x(p)+(-point(1)-m*point(2)))/(-m);
end
else
   x =  -2*pi:0.0001:2*pi;
   y = point(2)*ones(length(x),1);
end

end
