function FBZline = getFBZ(lattice)

a1 =  lattice.ax;
a2 =  lattice.ay;

D =  [a1' a2'];
B =  inv(D);
b1 =  B(1,:);
b2 =  B(2,:);

firstLine_x = linspace(0,b1(1), 100);
firstLine_y = zeros(1,100);

firstLine = [firstLine_x' firstLine_y'];

scondLine_x = ones(1,100)* b1(1);
secondLine_y =  linspace(0,b2(2), 100);

secondLine = [scondLine_x' secondLine_y'];

thirdLine_x =  linspace(b1(1),0, 100);
thirdLine_y = ones(1,100) * b2(2);

thirdLine =  [thirdLine_x' thirdLine_y'];

fourthLine_x = zeros(1,100);
fourthLine_y = linspace(b2(2),0, 100);

fourthLine = [fourthLine_x' fourthLine_y'];

FBZline =  pi * [firstLine;secondLine;thirdLine;fourthLine];


end