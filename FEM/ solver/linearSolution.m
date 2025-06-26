function linearSolution()



[KFF,PF] = assembleMat(nElem,nForce,nPresDisp, ...
		xcoords,elemNode,elemStiff,elemLoad, ...
		forceNode,forceValue,dispNode,dispValue,eqNum,nDof);

end
