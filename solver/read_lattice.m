function lattice = read_lattice(filename)

[a1, a2, N_list,numNode,xCoord,numElem,indx,numRot,rotIndx,axialPBC,rotPBC,disp_bndry,frcID, dmp,id1, id2, id3, id4]= read(filename);

lattice =  lattice_file;
lattice.ax = a1;
lattice.ay = a2;
lattice.N1 = N_list(1);
lattice.N2_list = N_list(2:end);
lattice.nNode =  numNode;
lattice.xCoord = xCoord;
lattice.def_xCoord = xCoord;
lattice.vCoord = zeros(size(xCoord,1),2);
lattice.nElem =  numElem;
lattice.indx = indx;
lattice.nRot = numRot;
lattice.rotIndx = rotIndx;
lattice.axialPBC = axialPBC;
lattice.rotPBC = rotPBC;
lattice.presDispID = disp_bndry(:,1);
lattice.presDispX = disp_bndry(:,2);
lattice.presDispY = disp_bndry(:,3);
lattice.frcID = frcID;
lattice.dmp = dmp;
lattice.id1 = id1;
lattice.id2 = id2;
lattice.id3 = id3;
lattice.id4 = id4;

end 