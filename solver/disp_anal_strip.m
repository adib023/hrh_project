function disp_anal_strip(lattice, strpID,strpElem,fix)

mainStrp =  strpID(:,3:4);


% K =  K([2*strpID(:,2:3)-1,2*strpID(:,2:3)],[2*strpID-1,2*strpID]);
n1 = strpID(1,1);
n2 = strpID(1,3);
a =  lattice.xCoord(n1,:) - lattice.xCoord(n2,:);

count = 0;
for k = -pi/abs(a(1)):0.01:pi/abs(a(1))  %pi/abs(a(1))/2
    [K,~] = getStiffness_matrix (lattice.xCoord,lattice.indx,lattice.rotIndx);
    count = count + 1;
    for p = 1: size(strpID,1)
        n1 = strpID(p,1);
        n2 = strpID(p,3);
        a =  lattice.xCoord(n1,:) - lattice.xCoord(n2,:);
        K(:,[2*n1-1,2*n1]) = K(:,[2*n1-1,2*n1]) * exp(i*k*a(1));
        K(:,[2*n2-1,2*n2]) = K(:,[2*n2-1,2*n2]) + K(:,[2*n1-1,2*n1]);
        
        n1 = strpID(p,2);
        n2 = strpID(p,4);
        a =  lattice.xCoord(n1,:) - lattice.xCoord(n2,:);
        K(:,[2*n1-1,2*n1]) = K(:,[2*n1-1,2*n1]) * exp(i*k*a(1));
        K(:,[2*n2-1,2*n2]) = K(:,[2*n2-1,2*n2]) + K(:,[2*n1-1,2*n1]);
        
        n1 = strpID(p,5);
        n2 = strpID(p,3);
        a =  lattice.xCoord(n1,:) - lattice.xCoord(n2,:);
        K(:,[2*n1-1,2*n1]) = K(:,[2*n1-1,2*n1]) * exp(i*k*a(1));
        K(:,[2*n2-1,2*n2]) = K(:,[2*n2-1,2*n2]) + K(:,[2*n1-1,2*n1]);
        
        n1 = strpID(p,6);
        n2 = strpID(p,4);
        a =  lattice.xCoord(n1,:) - lattice.xCoord(n2,:);
        K(:,[2*n1-1,2*n1]) = K(:,[2*n1-1,2*n1]) * exp(i*k*a(1));
        K(:,[2*n2-1,2*n2]) = K(:,[2*n2-1,2*n2]) + K(:,[2*n1-1,2*n1]);
        
    end
    mainStrp = sort(mainStrp(:));
%     mainStrp(mainStrp == fix(1)) = [];
%     mainStrp(mainStrp == fix(2)) = [];
    K_red = K(sort([2*mainStrp-1;2*mainStrp]),sort([2*mainStrp-1;2*mainStrp]));
    w2 =  eig(K_red);
    [A,w2] = eig(K_red,"vector");
    [w2,ind] = sort(w2);
    A = A(:,ind);
    w = sqrt(w2);
    w = real(w);
    w = sort(w);
	result(count,:) = [k;w];
end
% delete('../result/modeshape_strp/*')
% plot_mode(w,A,lattice.xCoord,strpElem,mainStrp,'../result/modeshape_strp/')
% % 
f = fopen('../result/strip_disp.txt','w');
write('../result/strip_disp.txt',result);
fclose(f);
% % % % % 

end

function plot_mode(w,A,xCoord,indx,strp,save_dir)

for p = 1:length(w)
    def_xCoord = xCoord;
    def_xCoord(strp(:),1) = def_xCoord(strp(:),1) + 10*A(1:2:size(A,1),p);
    def_xCoord(strp(:),2) = def_xCoord(strp(:),2) + 10*A(2:2:size(A,1),p);
    plot_lattice(def_xCoord,xCoord,indx,p,w(p),save_dir);
end

end

function plot_lattice(defcoord,coord,indx,count,w,save_dir)

% f= figure('visible','off');
figure
for i = 1 : size ( indx, 1 )
      obj1 =  indx ( i, 1);
      obj2 =  indx ( i, 2);
      x1 =  defcoord(obj1, 1);
      y1 =  defcoord(obj1, 2);
      x2= defcoord(obj2, 1);
      y2= defcoord(obj2, 2);
      plot([x1 x2], [y1 y2],'-k');      
      hold on;
      
    
end


for i = 1 : size ( indx, 1 )
      obj1 =  indx ( i, 1);
      obj2 =  indx ( i, 2);
      x1 =  coord(obj1, 1);
      y1 =  coord(obj1, 2);
      x2= coord(obj2, 1);
      y2= coord(obj2, 2);
      plot([x1 x2], [y1 y2],'--b');
%       xlim([1,2]);
%       ylim([0,1])
      
      hold on;
      
    
end


%axis equal;
axis off
title(strcat("$\omega =$",num2str(w)),"interpreter","latex","fontsize",20)
print("-dpng","-r300",strcat(save_dir,num2str(count)))


end