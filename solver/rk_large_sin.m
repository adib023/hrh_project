function rk_large_sin(kiri,startTime,finalTime,dt,pres_d,w,w2,w_count)

count = 0;

if w2 == 0
    result_folder =  '../result/rk_result/transLargeSin/';
else
    result_folder =  '../result/rk_result/transLargeSin/reflect_spec/';
end

if not(isfolder(result_folder))
    mkdir(result_folder)
end

deformFile = strcat(result_folder,'deformCoord',num2str(w_count),'.txt');


fid1 = fopen(deformFile,'w');

for t = startTime : dt : 2*pi/w*4
    
    
    if mod(count, (1/dt)/2) == 0
    %%%%%%%%%%%%%%% here I store the result for time steps with interval of 0.5s    
    %%%%%%%%%%%%%%% for all the time steps deformed coordinates are stored in a single file    
    t    
    %fid = fopen(filename,'w');
    %filename = strcat(result_folder,num2str(count),'.txt');
    write1by1(deformFile,kiri.xCoord)
    end    
    
    k1 = get_u_dot(kiri,zeros(4*kiri.nNode,1),pres_d,t,w,w2);
    k2 = get_u_dot(kiri, 0.5*k1*dt,pres_d,t,w,w2);
    k3 = get_u_dot(kiri, 0.5*k2*dt,pres_d,t,w,w2);
    k4 = get_u_dot(kiri, k3*dt,pres_d,t,w,w2);
    
    du =  dt/6 *( k1 + 2* k2 + 2* k3 + k4 );
    kiri =  updateKiri(kiri,du);
    count = count + 1;
end

if w2 == 0
otherInfo = strcat(result_folder,'otherInfo.txt');
fid2 = fopen(otherInfo,'w');
write(otherInfo, [pres_d w kiri.dmp(1,3)]);
fclose(fid2);
end
%%%%% write damping ratio, freq, and pres_d

fclose(fid1);

end




function u_dot = get_u_dot (kiri,add,pres_d,t,w,w2)

temp_kiri = vary_xCoord_vCoord ( kiri, add );
[K,f_int] = getStiffness_matrix(temp_kiri.xCoord,temp_kiri.indx,temp_kiri.rotIndx,0,1);
[~,f_int] = apply_PBC(temp_kiri.xCoord,temp_kiri.axialPBC,temp_kiri.rotPBC,K,f_int,0,1);

v  =  temp_kiri.vCoord;
v =  v';
v =  v(:);

for p = 1: size(kiri.dmp,1)
    %     n = kiri.dmp(p,1);
    %     c = kiri.dmp(p,2);
    %
    %     f_int(2*n-1) = f_int(2*n-1) + c*v(2*n-1);
    %     f_int(2*n) = f_int(2*n) + c*v(2*n);
    
    n1 = kiri.dmp(p,1);
    n2 = kiri.dmp(p,2);
    c =  kiri.dmp(p,3);
    
    f_int(2*n1-1) = f_int(2*n1-1) + c*(v(2*n1-1)-v(2*n2-1));
    f_int(2*n1) = f_int(2*n1) + c*(v(2*n1)-v(2*n2));
    f_int(2*n2-1) = f_int(2*n2-1) + c*(v(2*n2-1)-v(2*n1-1));
    f_int(2*n2) = f_int(2*n2) + c*(v(2*n2)-v(2*n1));
end



u_dot =  [v ;  - f_int];

% this section is only compatable for fixedNodes
% as in our problem there is no other condision rather than
% fixed conditions


num =  length(kiri.presDispID);

% if t <= 20*pi/w
%     alpha =  w/20/pi*t;
% else 
    alpha = 1; 
% end

u_dot(2*kiri.presDispID-1) = 0; 
u_dot(2*kiri.nNode + 2*kiri.presDispID-1) = 0;
u_dot(2*kiri.presDispID(num/2+1:num)) = pres_d*w*sin(w*t)+ 0.0001*w2*sin(w2*t);
u_dot(2*kiri.presDispID(1:num/2)) = 0;
% u_dot(2*kiri.presDispID(num/2+1:num)) = 0;
% u_dot(2*kiri.presDispID(1:num/2)) = -pres_v;
u_dot(2*kiri.nNode + 2*kiri.presDispID) = 0;
   


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




function plot_lattice(xCoord,indx)

figure

for i = 1 : size ( indx, 1 )
      obj1 =  indx ( i, 1);
      obj2 =  indx ( i, 2);
      x1 =  xCoord(obj1, 1);
      y1 =  xCoord(obj1, 2);
      x2= xCoord(obj2, 1);
      y2= xCoord(obj2, 2);
      k = indx(i,3);
      k_ref = indx(1,3);
      
if k == k_ref
    line = "-k";
else
    line = "-r";
end
plot([x1 x2], [y1 y2],line)
      
      hold on;
      
    
end



axis equal;
axis off
end

