%%%%function  finiteDiff ()
clc
clear all
close all
filename = "../file/re_hexa_re_pbc.txt";
hexa = read_lattice(filename);
nDof =  hexa.nNode *2;

u = 0.1*rand(nDof, 1);

u(2*hexa.presDispID -1) = hexa.presDispX;
u(2*hexa.presDispID) = hexa.presDispY;

def_xCoord = updateLattice(u,hexa.xCoord,hexa.pbc,hexa.bndry);  

elimID = [hexa.pbc(:)' hexa.presDispID'];

f = cal_dp_dx(def_xCoord, hexa.indx, hexa.rotIndx,hexa.bndry,hexa.pbc)
K = getStiffness_matrix(hexa.xCoord, hexa.indx, hexa.rotIndx, def_xCoord);



plot_lattice(hexa.xCoord,hexa.indx,hexa.bndry,hexa.pbc,'k-',f)




function def_xCoord =  updateLattice(x,xCoord,pbc,bndry)

%%%%% here xCoord means the initial xCoord before updating latest
%%%%% deformation

def_xCoord =  xCoord;

for n =  1:size(def_xCoord,1)
    p = n;
    if any(pbc(:) == n)
        [row,col] = find(pbc ==  n);
        p = bndry(row, col);
    end
    
    def_xCoord(n,1) = def_xCoord(n,1) + x(2*p-1);
    def_xCoord(n,2) = def_xCoord(n,2) + x(2*p);
end


end

function [f,K] = applyBC(f,K,lattice,def_xCoord,elimID)

%%%%%% considering wave_vector =  [0 0]

K = applyPBC(def_xCoord,lattice.bndry,lattice.pbc,K,[0 0]);

%elimID = [lattice.pbc(:)' lattice.presDispID'];

K([2*elimID-1, 2*elimID],:) = [];
K(:,[2*elimID-1, 2*elimID]) = [];
f([2*elimID-1, 2*elimID]) = [];
x([2*elimID-1, 2*elimID]) = [];

end

function plot_lattice(xCoord,indx,bndryID,pbcID,line,f)

%figure

for i = 1 : size ( indx, 1 )
      obj1 =  indx ( i, 1);
      obj2 =  indx ( i, 2);
      x1 =  xCoord(obj1, 1);
      y1 =  xCoord(obj1, 2);
      x2= xCoord(obj2, 1);
      y2= xCoord(obj2, 2);
      k = indx(i,3);
      k_ref = indx(1,3);
      
      
      plot([x1 x2], [y1 y2],line)
      
      hold on;
      
    
end

a = bndryID(:,1);
b = bndryID(:,4);
 
plot(xCoord(a,1),xCoord(a,2),"b--","linewidth",1);
hold on
plot(xCoord(b,1),xCoord(b,2),"r--","linewidth",1);
hold on

for n = 1:size(xCoord,1)
    idxNumber = num2str(n);
    text(xCoord(n,1),xCoord(n,2),idxNumber);
    hold on
end


axis equal;
axis off
end


