%%% routine to determine strip nodal coordinates in fem model

starting_point

he_length = 21.09;
re_length = 18.70;

%%% hexagon part

a1 = he_length*[2*sqrt(3)*cos(theta) 0];
a2 = he_length*[0 -2*cos(theta)];

x1  =  starting_point;


