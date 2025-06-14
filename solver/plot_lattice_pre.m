function plot_lattice_pre(xCoord,indx,rotIndx,fig_name)

numNode = size(xCoord,1);

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
plot([x1 x2], [y1 y2],'k-')
      
      hold on;
      
    
end




% scatter(xCoord(1:2:numNode,1),xCoord(1:2:numNode,2),30,"r","filled");
% hold on
% scatter(xCoord(2:2:numNode,1),xCoord(2:2:numNode,2),30,"r","filled");
% 
% scatter(xCoord(numNode_hexa+1:numNode_hexa+numNode_re,1),...
%     xCoord(numNode_hexa+1:numNode_hexa+numNode_re,2),10,"b","filled");
% hold on
% 
% scatter(xCoord(numNode_hexa+numNode_re+1:numNode,1),...
%     xCoord(numNode_hexa+numNode_re+1:numNode,2),10,"r","filled");
% hold on

% 


for n = 1:size(xCoord,1)
    idxNumber = num2str(n);
    text(xCoord(n,1),xCoord(n,2),idxNumber);
    hold on
end
% set(gca, 'fontsize', 10)

% 
% for p =  1:size(rotIndx,1)
%     
%     n = rotIndx(p,:);
%     x = xCoord(n(1:3),:);
%     x1 = (x(1,:)+x(2,:))/2;
%     x2 = (x(1,:)+x(3,:))/2;
%     
%     plot([x1(1) x2(1)], [x1(2) x2(2)],'b--','LineWidth', 0.5);
%                         
%     hold on
% end

axis equal;
axis off
% fileDir =  "./lattice_fig/";
% filename =  strcat(fileDir,fig_name);
% print("-dpng","-r300",filename);
% print("-depsc","-r300",filename);
end

