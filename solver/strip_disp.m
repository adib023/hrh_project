function strip_disp(lattice)

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

list = [list1;list2;list3;list4]


[K,~] = getStiffness_matrix(lattice.xCoord, lattice.indx, lattice.rotIndx,1,0);

main_strip_id = sort([2*list(:,1)-1;2*list(:,1)]);
left_strip_id = sort([2*list(:,2)-1;2*list(:,2)]);
right_strip_id =sort([2*list(:,3)-1;2*list(:,3)]);

fixedNode =  [lattice.id1(3,1) lattice.id4(3,N2)];
fixedDOF_ID = [2*fixedNode-1 2*fixedNode];


does_ilim =  ismember(main_strip_id,fixedDOF_ID);
ilimID = find(does_ilim == 1);
% ilimID
red_main_strip = main_strip_id;
%red_main_strip(ilimID)= [];

 mu_list = -2*pi:0.01:2*pi;


for p = 1: length(mu_list)
    mu = mu_list(p);
    K_disp = K;
    K_disp(:,main_strip_id) = K_disp(:,main_strip_id) + exp(-i*mu)* K_disp(:,left_strip_id);
    K_disp(:,main_strip_id) = K_disp(:,main_strip_id) + exp(i*mu)* K_disp(:,right_strip_id);
    
    K_red = K_disp(red_main_strip,red_main_strip);
    %K_red = K_disp(main_strip_id,main_strip_id);
    w = sqrt(eig(K_red));
    w = sort(real(w));
    w_list(:,p) = w;
end

figure 
plot(w,'r.')

figure
for p = 1:size(w_list,1)
    plot(mu_list./pi,w_list(p,:),'k-')
    hold on
end
xlabel('$\mu/\pi$','Interpreter','latex');
ylabel('$\omega$','Interpreter','latex');
set(gca,'fontsize',15)
set(gca,'fontname','times new roman')

end