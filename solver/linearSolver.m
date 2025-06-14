function [u, def_xCoord] = linearSolver(lattice)
nDof = 2*lattice.nNode;

F = zeros(nDof,1);
K = getStiffness_matrix(lattice.xCoord,lattice.indx,lattice.rotIndx,lattice.def_xCoord);

K= applyPBC(lattice.xCoord,lattice.bndryID,lattice.pbcID,K,[0 0]);


F =  applyPresDisp(F,K,lattice.disp_bndry);
F = applyPresFrc(F,lattice.frc_bndry);




illiminateIDs = [2*lattice.disp_bndry(:,1)-1;2*lattice.disp_bndry(:,1);
    2*lattice.pbcID(:)-1;2*lattice.pbcID(:)];



K(illiminateIDs,:) = [];
K(:,illiminateIDs) = [];
F(illiminateIDs,:) = []; %%%% for applying PBC

u =  K^-1*F;

T = eye(nDof);
T(:,illiminateIDs) = [];
u = T*u;

%lattice.disp_bndry(:,1)
u(2*lattice.disp_bndry(:,1)) = u(2*lattice.disp_bndry(:,1)) + lattice.disp_bndry(:,3); 

u = setDispForPBC(lattice.bndryID,lattice.pbcID,u);

lattice.def_xCoord(:,1) = lattice.def_xCoord(:,1) + u(1:2:nDof);
lattice.def_xCoord(:,2) = lattice.def_xCoord(:,2) + u(2:2:nDof);
def_xCoord = lattice.def_xCoord;

file = "../file/linearResult.txt";

fid =  fopen(file,"w");

u_main = u;
u_main([2*lattice.pbcID(:)-1;2*lattice.pbcID(:)]) = [];
write(file,u_main);
fclose(fid);

plot_def_lattice(lattice.def_xCoord,lattice.xCoord,lattice.indx);
end


function F =  applyPresDisp(F,K,disp_bndry)

id = disp_bndry(:,1);
disp = disp_bndry(:,2:3);


for p = 1:length(id)
    F = F - K(:,2*id(p)-1)*disp(p,1);
    F = F - K(:,2*id(p))*disp(p,2);
end

end

function F =  applyPresFrc(F,frc_bndry)
id = frc_bndry(:,1);
frc = frc_bndry(:,2:3);

for p = 1:length(id)
    F(2*id(p)-1) = F(2*id(p)-1) + frc(p,1);
    F(2*id(p)) = F(2*id(p)) + frc(p,2); 
end
end

function u = setDispForPBC(bndryID,pbcID,u)

for p = 1: size(bndryID,1)
    
    nb = bndryID(p,1);
    n =  pbcID(p,3);
    u(2*n-1) = u(2*n-1)+u(2*nb-1);
    u(2*n) = u(2*n)+u(2*nb);
    
    nb = bndryID(p,2);
    n =  pbcID(p,4);
    u(2*n-1) = u(2*n-1)+u(2*nb-1);
    u(2*n) = u(2*n)+u(2*nb);
    
    nb = bndryID(p,3);
    n =  pbcID(p,1);
    u(2*n-1) = u(2*n-1)+u(2*nb-1);
    u(2*n) = u(2*n)+u(2*nb);
    
    nb = bndryID(p,4);
    n =  pbcID(p,2);
    u(2*n-1) = u(2*n-1)+u(2*nb-1);
    u(2*n) = u(2*n)+u(2*nb);
    
end

end


function plot_def_lattice(defcoord,coord,indx)

%f= figure('visible','off');
figure
for i = 1 : size ( indx, 1 )
      obj1 =  indx ( i, 1);
      obj2 =  indx ( i, 2);
      x1 =  defcoord(obj1, 1);
      y1 =  defcoord(obj1, 2);
      x2= defcoord(obj2, 1);
      y2= defcoord(obj2, 2);
%       plot([x1 x2], [y1 y2],'-ok',...
%                             'LineWidth', 1,...
%                             'MarkerFaceColor','r','MarkerEdgeColor','r',...
%                             'MarkerSize',5);
%       xlim([1,2]);
%       ylim([0,1])
plot([x1 x2], [y1 y2],'-k');      
      hold on;
      
    
end

% 
% for i = 1 : size ( indx, 1 )
%       obj1 =  indx ( i, 1);
%       obj2 =  indx ( i, 2);
%       x1 =  coord(obj1, 1);
%       y1 =  coord(obj1, 2);
%       x2= coord(obj2, 1);
%       y2= coord(obj2, 2);
%       plot([x1 x2], [y1 y2],'--b');
% %       xlim([1,2]);
% %       ylim([0,1])
%       
%       hold on;
%       
%     
% end

% 
%axis equal;
axis off
% title(strcat("$\omega =$",num2str(w)),"interpreter","latex","fontsize",20)
% print("-dpng","-r300",strcat("../result/modeshape/",num2str(count)))


end