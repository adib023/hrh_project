function deformed_config(lattice,fVal)
nDof = 2*lattice.nNode;

F = zeros(nDof,1);

K = getStiffness_matrix(lattice.xCoord,lattice.indx,lattice.rotIndx);

for p = 1:size(lattice.frcID,1)
    F(2*lattice.frcID(p,1)) = -fVal;
    F(2*lattice.frcID(p,2)) = fVal;
end
M = eye(nDof);
[K,~,F] = applyPBC(lattice.xCoord,lattice.bndryID,lattice.pbcID,K,M,[0 0],F);
F
u = K^-1*F;

file_wo_pbc_info = "../file/re_hexa_re.txt";
[nNode,xCoord,nElem,indx]= read(file_wo_pbc_info);

def_xCoord =  zeros(size(xCoord));

for p = 1:nNode
    def_xCoord(p,1) =  xCoord(p,1)+ u(2*p-1);
    def_xCoord(p,2) =  xCoord(p,2)+ u(2*p-1);
end
% 
% xCoord
% def_xCoord
% indx

plot_def_lattice(def_xCoord,xCoord,indx,fVal)

end


function plot_def_lattice(defcoord,coord,indx,fVal)

%f= figure('visible','off');
figure
% for i = 1 : size ( indx, 1 )
%     
%       obj1 =  indx ( i, 1);
%       obj2 =  indx ( i, 2);
%       x1 =  defcoord(obj1, 1);
%       y1 =  defcoord(obj1, 2);
%       x2= defcoord(obj2, 1);
%       y2= defcoord(obj2, 2);
% %       plot([x1 x2], [y1 y2],'-ok',...
% %                             'LineWidth', 1,...
% %                             'MarkerFaceColor','r','MarkerEdgeColor','r',...
% %                             'MarkerSize',5);
% %       xlim([1,2]);
% %       ylim([0,1])
% plot([x1 x2], [y1 y2],'-k');      
%       hold on;
%       
%     
% end


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
title(strcat("$f =$",num2str(fVal)),"interpreter","latex","fontsize",20)
%print("-dpng","-r300",strcat("../result/modeshape/",num2str(count)))


end

