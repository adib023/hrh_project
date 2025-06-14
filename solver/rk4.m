function rk4(kiri,startTime,dt,finalTime,w_ext,w,A)

result_folder = strcat("../result/rk/bic");

if not(isfolder(result_folder))
    mkdir(result_folder)
end

f_file = strcat(result_folder,"/frequency.txt");

f_fid = fopen(f_file,'w');
write(f_file,[w_ext w]);
fclose(f_fid);


coord_file =  strcat(result_folder,"/def_coord.txt");
if w_ext ==  0
    coord_file =  strcat(result_folder,"/def_coord_zeroExt.txt");
end
fid = fopen(coord_file,'w');

count = 0;
for t = startTime : dt :finalTime
     
    if mod(count,1/dt/10) == 0
t
    write1by1(coord_file,kiri.xCoord);
    end
    
    k1 = get_u_dot(kiri,t,zeros(4*kiri.nNode,1),w,w_ext,A);
    k2 = get_u_dot(kiri,t + 0.5*dt, 0.5*k1*dt,w,w_ext,A);
    k3 = get_u_dot(kiri,t + 0.5*dt, 0.5*k2*dt,w,w_ext,A);
    k4 = get_u_dot(kiri,t+ dt, k3*dt,w,w_ext,A);
    
    du =  dt/6 *( k1 + 2* k2 + 2* k3 + k4 );
    kiri =  updateKiri(kiri,du);
    count =  count + 1;
end

fclose(fid);

end

function u_dot = get_u_dot(kiri,t,add,w,w_ext,A)
temp_kiri = vary_xCoord_vCoord ( kiri, add );


[K,f_int] = getStiffness_matrix(temp_kiri.xCoord,temp_kiri.indx,temp_kiri.rotIndx,0,1);
[~,f_int] = apply_PBC(temp_kiri.xCoord,temp_kiri.axialPBC,temp_kiri.rotPBC,K,f_int,0,1);

v  =  temp_kiri.vCoord;
v =  v';
v =  v(:);

%%%% apply damping


for p = 1: size(kiri.dmp,1)
    
    n1 = kiri.dmp(p,1);
    n2 = kiri.dmp(p,2);
    c =  kiri.dmp(p,3);
    
    f_int(2*n1-1) = f_int(2*n1-1) + c*(v(2*n1-1)-v(2*n2-1));
    f_int(2*n1) = f_int(2*n1) + c*(v(2*n1)-v(2*n2));
    f_int(2*n2-1) = f_int(2*n2-1) + c*(v(2*n2-1)-v(2*n1-1));
    f_int(2*n2) = f_int(2*n2) + c*(v(2*n2)-v(2*n1));
end


u_dot = [v ; - f_int];

num =  length(kiri.presDispID);
% if t <= 20*pi/w
%     alpha =  w/20/pi*t;
% else 
%     alpha = 1; 
% end


if t>=100 && t<= 200
    pres_v = 1/100;
else
    pres_v = 0;
end

u_dot(2*kiri.presDispID-1) = 0; 
u_dot(2*kiri.nNode + 2*kiri.presDispID-1) = 0;
u_dot(2*kiri.presDispID) = kiri.presDispY*pres_v;
%u_dot(2*kiri.presDispID(1:num/2)) = 0;
u_dot(2*kiri.nNode + 2*kiri.presDispID) = 0;


% % % % %  important section for future, as now I am trying to have
% deformed configuration without excitation, but when excitation will be
% applied, I need to uncomment the section
%%%%%% define excitation
%if w_ext > 0
u_dot(2*kiri.frcID) =  0.001*w_ext*sin(w_ext*t);
u_dot(2*kiri.nNode +2*kiri.frcID) =  0;
%end

end




function temp_kiri =  vary_xCoord_vCoord(kiri, add)
temp_kiri = kiri;
count =  0;
for p =  1:temp_kiri.nNode
    temp_kiri.xCoord(p,1) = temp_kiri.xCoord(p,1) + add(count+1);
    temp_kiri.xCoord(p,2) = temp_kiri.xCoord(p,2) + add(count+2);
    count = count + 2;
end

for p =  1:temp_kiri.nNode
    temp_kiri.vCoord(p,1) = temp_kiri.vCoord(p,1) + add(count+1);
    temp_kiri.vCoord(p,2) = temp_kiri.vCoord(p,2) + add(count+2);
    count = count + 2;
end

end





function kiri = updateKiri( kiri, du )

count =  0;

for p  =  1: kiri.nNode
    kiri.xCoord(p,1) = kiri.xCoord(p,1) + du(count+1);
    kiri.xCoord(p,2) = kiri.xCoord(p,2) + du(count+2);
    count = count + 2 ;
end

for p  =  1: kiri.nNode
    kiri.vCoord(p,1) = kiri.vCoord(p,1) + du(count+1);
    kiri.vCoord(p,2) = kiri.vCoord(p,2) + du(count+2);
    count = count + 2 ;
end

end





