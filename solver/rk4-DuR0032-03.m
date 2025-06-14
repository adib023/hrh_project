function rk4(kiri,appliedForce,startTime,dt,finalTime,w,count,pres_v)

%%%%%%%%%%%%%%%%%%%%%% must note here\
%%%%%%%%%%% to run the parfor in linux
%%%%%%%%%%% we set w as input argument instead of seting it as a property
%%%%%%% of 'applied force'

result_folder = strcat("../result/rk_result/transTxtFiles_with_",num2str(pres_v*10000),"/",num2str(count));

if not(isfolder(result_folder))
    mkdir(result_folder)
end

f_file = strcat(result_folder,"/frequency.txt");

f_fid = fopen(f_file,'w');
write(f_file,w);
fclose(f_fid);


%count = 0;
count2 = 0;
% t_record = [];
% xCoord_t = [];
% %vCoord_t = [];

for t = startTime : dt :finalTime
    
    %if mod(count, 100) == 0
    
    if t == 0
        xc_file = strcat("../result/rk_result/transTxtFiles1/",num2str(count),"/10001.txt");
        xc =  read(xc_file);
        kiri.xCoord = xc;
    end
    
    %if t>=900
        count2 = count2+1;
        file =  strcat(result_folder,"/",num2str(count2),'.txt'); 
        fid =  fopen(file,"w");
        write(file,kiri.xCoord);
        fclose(fid);
    %end
    %end
    
    k1 = get_u_dot(kiri,appliedForce,t,zeros(4*kiri.nNode,1),w,pres_v);
    k2 = get_u_dot(kiri,appliedForce,t + 0.5*dt, 0.5*k1*dt,w,pres_v);
    k3 = get_u_dot(kiri,appliedForce,t + 0.5*dt, 0.5*k2*dt,w,pres_v);
    k4 = get_u_dot(kiri,appliedForce,t+ dt, k3*dt,w,pres_v);
    
    du =  dt/6 *( k1 + 2* k2 + 2* k3 + k4 );
    kiri =  updateKiri(kiri,du);
    
    %count = count + 1 ;
end


end

function u_dot = get_u_dot(kiri,appliedForce,t,add,w,pres_v)
temp_kiri = vary_xCoord_vCoord ( kiri, add );
f_app  =  zeros(2*kiri.nNode,1);
%f_app  =  force(temp_kiri,appliedForce,t);
[K,f_int] = getStiffness_matrix(temp_kiri.xCoord,temp_kiri.indx,temp_kiri.rotIndx);
[~,f_int] = apply_PBC(temp_kiri.xCoord,temp_kiri.axialPBC,temp_kiri.rotPBC,K,f_int);

v  =  temp_kiri.vCoord;
v =  v';
v =  v(:);

%%%% apply damping

%if t<50
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
%end


u_dot = [v ; f_app - f_int];

num =  length(kiri.presDispID);
%w = appliedForce.w;
A =  appliedForce.A;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% please change this alpha when doing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calculation from t = 0, and when there
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% is no constant velocity applied on top
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% boundary
% if t <= 20*pi/w
%     alpha =  w/20/pi*t;
% else 
    alpha = 1; 
% end

u_dot(2*kiri.presDispID-1) = 0; 
u_dot(2*kiri.nNode + 2*kiri.presDispID-1) = 0;
u_dot(2*kiri.presDispID(num/2+1:num)) = alpha*A*w*sin(w*t)+pres_v;
u_dot(2*kiri.presDispID(1:num/2)) = 0;
u_dot(2*kiri.nNode + 2*kiri.presDispID) = 0;

%u_dot(2*kiri.nNode + 2*kiri.presDispID(num/2+1:num)) = 0;


% if t>=bndryStart & t<=bndryFinish
%     u_dot(2*kiri.presDispID) = u_dot(2*kiri.presDispID) + kiri.presDispY/...
%         (bndryFinish - bndryStart);
% end


% this section is only compatable for fixedNodes
% as in our problem there is no other condision rather than
% fixed conditions
% if t == bndryStart
%     kiri.vCoord(kiri.presDispID,1:2) = kiri.vCoord(kiri.presDispID,1:2)+[kiri.presDispX kiri.presDispY]/(-bndryStart+bndryFinish);
% end
% 
% w = appliedForce.w;
% A =  appliedForce.A;
% 
% if t == bndryFinish
%     num =  length(kiri.presDispID);
%     kiri.vCoord(kiri.presDispID(num/2+1:num),1) = 0;
%     kiri.vCoord(kiri.presDispID(num/2+1:num),2) = A*w*cos(w*t);
% end
% 
% 
% % 
% % if t>bndryStart && t < bndryFinish
% for p =  1:length(kiri.presDispID)/2
%     %         u_dot(2*kiri.presDispID(p)-1) = 0;
%     u_dot(2*kiri.nNode + 2*kiri.presDispID(p)-1) = 0;
%     
%     %         u_dot(2*kiri.presDispID(p)) = -0.1;
%     u_dot(2*kiri.nNode + 2*kiri.presDispID(p)) = 0;
%     
% end
% for p =  length(kiri.presDispID)/2+1:length(kiri.presDispID)
%     
%     %         u_dot(2*kiri.presDispID(p)-1) = 0;
%     u_dot(2*kiri.nNode + 2*kiri.presDispID(p)-1) = 0;
%     
%     %         u_dot(2*kiri.presDispID(p)) = 0.1;
%     u_dot(2*kiri.nNode + 2*kiri.presDispID(p)) = A*w^2*sin(w*t);
%     
% end
% end




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

function f =  force(kiri,appliedForce,t)
nodeId  = appliedForce.id;
omega = appliedForce.w;
N = appliedForce.N;
% A =  appliedForce.inAmp;

f  = zeros(2*kiri.nNode,1);
%if t <= N*pi/omega 
 f(2*nodeId(1)-1) = 0.001*omega^2*sin(omega*t);
 f(2*nodeId(1)) = 0*sin(omega*t);
 f(2*nodeId(2)-1) = -0.001*omega^2*sin(omega*t);
 f(2*nodeId(2)) = 0*sin(omega*t);
%end
% if appliedForce.Direc == "x"
%     f(2*nodeId-1) = A*sin(omega*t);
% elseif appliedForce.Direc == "y"
%     f(2*nodeId) = A*sin(omega*t);
% elseif appliedForce.Direc == "-x"
%     f(2*nodeId-1) = -A*sin(omega*t);
% elseif appliedForce.Direc == "-y"
%     f(2*nodeId) = -A*sin(omega*t);
% end

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

% 
% scatter(xCoord(1:numNode_re1,1),xCoord(1:numNode_re1,2),20,"r","filled");
% hold on
% 
% scatter(xCoord(numNode_re1+1:numNode_re1+numNode_h,1),xCoord(numNode_re1+1:numNode_re1+numNode_h,2),20,"b","filled");
% hold on
% 
% scatter(xCoord(numNode_re1+numNode_h+1:numNode,1),xCoord(numNode_re1+numNode_h+1:numNode,2),20,"r","filled");
% hold on


% 
% for n = 1:size(xCoord,1)
%     idxNumber = num2str(n);
%     text(xCoord(n,1),xCoord(n,2),idxNumber);
%     hold on
% end



axis equal;
axis off
% fileDir =  "./lattice_fig/";
% filename =  strcat(fileDir,fig_name);
% print("-dpng","-r300",filename);
% print("-depsc","-r300",filename);
end


function [] =  write1by1(fname, data_to_be_stored)
  
  [nr,nc] = size (data_to_be_stored);

  file = fopen(fname, "a");
  
  for i = 1: nr
    for j = 1 : nc
      fprintf(file,"%4.10f\t", data_to_be_stored(i,j));
    end
    fprintf(file,"\n");
  end
  %fprintf(file,"\n");
  
  fclose ( file );

end
