clear all
clc
close all


filename = "../file/re_hexa_re_pbc.txt";
hexa = read_lattice(filename);

nNode = hexa.nNode;
xCoord = hexa.xCoord;
elem =  hexa.indx(1,1:2);
k = hexa.indx(1,3);
l_0 = hexa.indx(1,4);
u =  0.01*rand(2*nNode,1);
def_xCoord = updateLattice(xCoord,u);

xc = xCoord(elem,:);
def_xc =  def_xCoord(elem,:);
deb_xc = def_xc;


f  =  getAxialForce(deb_xc(1,:),deb_xc(2,:),k,l_0);
debK =  zeros(4);

count = 0;
TOL =  1e-6;

for p = 1:2
    for q = 1:2
        count = count + 1;
        deb_xc =  def_xc;
        deb_xc(p,q) = deb_xc(p,q)+ TOL;
        f1 = getAxialForce(deb_xc(1,:),deb_xc(2,:),k,l_0);
        debK(count,:) = (f1-f)/TOL;
    end
end


K =  getAxialSpringK(deb_xc(1,:),deb_xc(2,:),xc(1,:),xc(2,:),k)
debK
K-debK


function def_xCoord = updateLattice(xCoord,u)

def_xCoord = xCoord;

for n = 1:size(xCoord,1)
    def_xCoord(n,1) = def_xCoord(n,1) + u(2*n-1);
    def_xCoord(n,2) = def_xCoord(n,1) + u(2*n);
end
end

