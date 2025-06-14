clear all
clc
close all
addpath("../pre")

filename = "../file/re_hexa_re_pbc.txt";
hexa = read_lattice(filename);

elemID =  5;
nNode = hexa.nNode;
xCoord = hexa.xCoord;
elem =  hexa.rotIndx(elemID,1:3);
k_rot = hexa.rotIndx(elemID,4);
th_0 = hexa.rotIndx(elemID,5);
u =  0.1*rand(2*nNode,1)
def_xCoord = updateLattice(xCoord,u);

xc =  xCoord(elem,:);
def_xc = def_xCoord(elem,:);

xc_debug =  def_xc;

f_t =  getTorForce(xc_debug,k_rot,th_0);

TOL =  1e-6;

K_debug =   zeros(6);
count =  0;
for p = 1:3

    for q = 1:2
        count =  count + 1;
        xc_debug =  def_xc;
        xc_debug(p,q) = xc_debug(p,q)+TOL;
        f_t1 =  getTorForce(xc_debug,k_rot,th_0);
        K_debug(count,:) =  (f_t1-f_t)/TOL;
    end
end



K =  getRotK(xc,def_xc,k_rot)
K_debug
K-K_debug

function def_xCoord = updateLattice(xCoord,u)

def_xCoord = xCoord;

for n = 1:size(xCoord,1)
    def_xCoord(n,1) = def_xCoord(n,1) + u(2*n-1);
    def_xCoord(n,2) = def_xCoord(n,1) + u(2*n);
end
end

