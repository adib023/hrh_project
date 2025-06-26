function plotDisp(probType,xNode,elemNode,elemMatIndx,elemType, ... 
		uDisp,scatNode)

if probType == 6

	plotScatterDisp(xNode,elemNode,elemMatIndx,elemType,uDisp,scatNode) ; 
else 
	error('not yet implemented');
end

end

function plotScatterDisp(xNode,elemNode,elemMatIndx,elemType,uDisp,scatNode) 

nNode = size(xNode,1);
count = 0 ;

for p = 1:nNode
	uv(p,1:3) = real(uDisp(3*p-2:3*p)); 
end

count = 0 ; 
for p = 1:size(elemNode,1)

	if elemMatIndx(p) == 0 
		continue;
	end

	count = count + 1 ; 

	en4(count,1:4) = [elemNode(p,1) elemNode(p,3) ... 
			elemNode(p,9) elemNode(p,7)] ; 
end


figure
patch('Vertices',xNode(:,1:2),'Faces',en4,'FaceVertexCData',uv(:,1), ...
		'FaceColor','interp','EdgeColor','none') ; 
colorbar
axis equal
print('-dpng','-r300','m1') ; 

figure
patch('Vertices',xNode(:,1:2),'Faces',en4,'FaceVertexCData',uv(:,2), ...
		'FaceColor','interp','EdgeColor','none') ; 
colorbar
axis equal
print('-dpng','-r300','m2') ; 

figure
patch('Vertices',xNode(:,1:2),'Faces',en4,'FaceVertexCData',abs(uv(:,3)), ...
		'FaceColor','interp','EdgeColor','none') ; 
colorbar
axis equal
print('-dpng','-r300','m3') ; 

end
