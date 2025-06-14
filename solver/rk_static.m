function rk_static(kiri,startTime,finalTime,dt,pres_v)

count = 0;

result_folder =  strcat('../result/rk_result/transStat_v',num2str(pres_v*1000),'.txt');
fid = fopen(result_folder,'w');
% 
% if not(isfolder(result_folder))
%     mkdir(result_folder)
% end


for t = startTime : dt : finalTime
    
    
    if mod(count, (1/dt)/2) == 0
    %%%%%%%%%%%%%%% here I store the result for time steps with interval of 0.5s    
    %%%%%%%%%%%%%%% for all the time steps deformed coordinates are stored in a single file    
    t    
    %fid = fopen(filename,'w');
    %filename = strcat(result_folder,num2str(count),'.txt');
    write1by1(result_folder,kiri.xCoord)
    end    
    
    k1 = get_u_dot(kiri,zeros(4*kiri.nNode,1),pres_v);
    k2 = get_u_dot(kiri, 0.5*k1*dt,pres_v);
    k3 = get_u_dot(kiri, 0.5*k2*dt,pres_v);
    k4 = get_u_dot(kiri, k3*dt,pres_v);
    
    du =  dt/6 *( k1 + 2* k2 + 2* k3 + k4 );
    kiri =  updateKiri(kiri,du);
    count = count + 1;
end

fclose(fid);
end




function u_dot = get_u_dot (kiri,add,pres_v)

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



u_dot(2*kiri.presDispID-1) = 0; 
u_dot(2*kiri.nNode + 2*kiri.presDispID-1) = 0;
u_dot(2*kiri.presDispID(num/2+1:num)) = pres_v;
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

