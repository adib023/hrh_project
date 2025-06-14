function strip_modeshape(lattice,mu)

N1 =  lattice.N1;
N2 = sum(lattice.N2_list);

%%%% node 1 list

lattice.id1
list1(:,1) =  lattice.id1(3,1:N2);
list1(:,2) =  lattice.id1(2,1:N2);
list1(:,3) =  lattice.id1(4,1:N2);

list2(:,1) =  lattice.id2(3,1:N2);
list2(:,2) =  lattice.id2(2,1:N2);
list2(:,3) =  lattice.id2(4,1:N2);

list3(:,1) =  lattice.id3(3,1:N2);
list3(:,2) =  lattice.id3(2,1:N2);
list3(:,3) =  lattice.id3(4,1:N2);

list4(:,1) =  lattice.id4(3,1:N2);
list4(:,2) =  lattice.id4(2,1:N2);
list4(:,3) =  lattice.id4(4,1:N2);

list = [list1;list2;list3;list4];


[K,~] = getStiffness_matrix(lattice.xCoord, lattice.indx, lattice.rotIndx,1,0);

mainStripNode = sort(list(:,1));

main_strip_id = sort([2*list(:,1)-1;2*list(:,1)]);
left_strip_id = sort([2*list(:,2)-1;2*list(:,2)]);
right_strip_id =sort([2*list(:,3)-1;2*list(:,3)]);

fixedNode =  [lattice.id1(3,1) lattice.id4(3,N2)];
fixedDOF_ID = [2*fixedNode-1 2*fixedNode];

does_ilim =  ismember(main_strip_id,fixedDOF_ID);
ilimID = find(does_ilim == 1);

red_main_strip = main_strip_id;
red_main_strip(ilimID)= []; 


K_disp = K;
K_disp(:,main_strip_id) = K_disp(:,main_strip_id) + exp(-i*mu)* K_disp(:,left_strip_id);
K_disp(:,main_strip_id) = K_disp(:,main_strip_id) + exp(i*mu)* K_disp(:,right_strip_id);

K_red = K_disp(red_main_strip,red_main_strip);
%K_red = K_disp(main_strip_id,main_strip_id);
[A,w2] = eig(K_red,'vector');
[w2,ind] = sort(w2);
w =  sqrt(w2);
w = real(w);

A = A(:,ind);
T = eye(length(main_strip_id));
T(:,ilimID) = [];
A = T*A;



%%%%%%%%% plot mode shape %%%%%%%%%%%
figure
plot(w,'r.')
print('-dpng','-r300','../result/strip_modeshape/frequency')

for p = 1:length(w)
    plot_strip_modeshape(lattice,mainStripNode, A(:,p), mu, w(p),p)
end
end

function plot_strip_modeshape(lattice,mainStripNode, A, mu, w,p)

def_coord = lattice.xCoord;

for n = 1:length(mainStripNode)
    node  =  mainStripNode(n);
    def_coord(node,1:2) = def_coord(node,1:2) + A(2*n-1:2*n)'; 
end

a = figure('Visible','off');
for elem = 1:size(lattice.indx,1)
    node1 = lattice.indx(elem,1);
    node2 = lattice.indx(elem,2);

    if ismember(node1, mainStripNode) & ismember(node2, mainStripNode)
        x1 = def_coord(node1,1);
        x2 = def_coord(node2,1);
        y1 = def_coord(node1,2);
        y2 = def_coord(node2,2);

        plot([x1 x2],[y1 y2],'k-');
        hold on

        x1 = lattice.xCoord(node1,1);
        x2 = lattice.xCoord(node2,1);
        y1 = lattice.xCoord(node1,2);
        y2 = lattice.xCoord(node2,2);

        plot([x1 x2],[y1 y2],'r--');
        hold on
    end
end

scatter(def_coord(mainStripNode,1),def_coord(mainStripNode,2),10,"r","filled"); %%%%% 10 
hold on

axis equal
axis off
values  = strcat('(',num2str(mu),',',num2str(w),')');
title(strcat('$(\mu,\omega) =',values,'$'),'Interpreter','latex')

print('-dpng','-r300',strcat('../result/strip_modeshape/',num2str(p)))


end