function defrmdConfig(lattice)

tolx = 1e-4;
tolf = 1e-1;
nDof  =  lattice.nNode *2;

x = zeros(nDof,1);
xCoord = lattice.xCoord;
nStep = abs(lattice.presDispY(end))/0.001;

T = eye(lattice.nNode*2);
T(:,[2*lattice.presDispID-1,2*lattice.presDispID]) = [];

for step = 1:nStep
    step
    
    x(2*lattice.presDispID -1) = lattice.presDispX/nStep;
    x(2*lattice.presDispID) = lattice.presDispY/nStep;
    
    
    [xCoord, rotIndx] = updateLattice(x,xCoord,lattice.rotIndx);
    %elimID = [lattice.pbc(:)' lattice.presDispID'];
    
    for itr =  1:200
        
        [K,f] = getStiffness_matrix(xCoord, lattice.indx, rotIndx,1,1);
        [K,f] = apply_PBC(xCoord,lattice.axialPBC,lattice.rotPBC,K,f,1,1);
        
        K([2*lattice.presDispID-1,2*lattice.presDispID],:)= [];
        K(:,[2*lattice.presDispID-1,2*lattice.presDispID])= [];
        f([2*lattice.presDispID-1,2*lattice.presDispID],:)= [];
        
        x = - K^-1*f;
        
        if itr == 1
            norm_f0 = norm(f);
        end
        %err(itr)  =  norm(x-x_i);
        
        
        %%%%%%% I also need to check the err in force, but I am not doing
        %%%%%%% it right now. When I will do it, will check the code sent
        %%%%%%% by professor
        
        x = T*x;
        xc_i =  xCoord(:);
        [xCoord, rotIndx] = updateLattice(x,xCoord,rotIndx);
        
        err_f(itr) = norm(f)/norm_f0;
        err_x(itr) = norm(xCoord(:)-xc_i)/norm(xc_i);

        
        if err_f(itr)<=tolf & err_x(itr)<=tolx
            %             figure
            %             plot_lattice(xCoord,lattice.indx,'-k')
            %             plot_lattice(lattice.xCoord,lattice.indx,'--b')
            %             title(num2str(step))
            %               print('-r300','-dpng',strcat('../result/NRresult/',num2str(step)))
            itr
%             figure
% plot_lattice(xCoord,lattice.indx,'-k')
% plot_lattice(lattice.xCoord,lattice.indx,'--b')
% print('-r300','-dpng',strcat('../result/NRresult/',num2str(step)))
            break
        end

        if itr == 200
            figure
            plot_lattice(lattice.xCoord,lattice.indx,'r-')
            plot_lattice(xCoord,lattice.indx,'-k')
            title(strcat('Step ',num2str(step)),'(cant converge)')
            print('-r300','-dpng',strcat('../result/NRresult/',num2str(step)))
            figure
            plot(err_f)
            title(strcat('Force Error step ',num2str(step)))
            print('-r300','-dpng','../result/NRresult/forceError')
            figure
            plot(err_x)
            title(strcat('Displacement Error step ',num2str(step)))
            print('-r300','-dpng','../result/NRresult/dispError')
            error(strcat('code does not converges at step', num2str(step)))
        end
        
    end
    
end
figure

%plot_lattice(lattice.xCoord,lattice.indx,'r-')
plot_lattice(xCoord,lattice.indx,'-k')
%title(Deform)
print('-r300','-dpng',strcat('../result/NRresult/','dfrmLattice'))
print('-r300','-depsc',strcat('../result/NRresult/','dfrmLattice'))
%%%%% calculate modeshapes for deformed lattice
[K,f] = getStiffness_matrix(xCoord, lattice.indx, rotIndx,1,0);
[K,~] = apply_PBC(xCoord,lattice.axialPBC,lattice.rotPBC,K,f,1,0);

K([2*lattice.presDispID-1,2*lattice.presDispID],:)= [];
K(:,[2*lattice.presDispID-1,2*lattice.presDispID])= [];

[A,w2] = eig(K,'vector');
[w2,ind] = sort(w2);
A = A(:,ind);
A = T*A;
w = sqrt(w2);


filename =  '../result/NRresult/eigenResult.txt';
f = fopen(filename,'w');
write(filename,w);
write(filename,real(A));
write(filename,imag(A));

fclose(f);

fid =  fopen('../result/NRresult/def_xCoord.txt','w');
write('../result/NRresult/def_xCoord.txt',xCoord);
fclose(fid);
end


function plot_lattice(xCoord,indx,line)



for i = 1 : size ( indx, 1 )
      obj1 =  indx ( i, 1);
      obj2 =  indx ( i, 2);
      x1 =  xCoord(obj1, 1);
      y1 =  xCoord(obj1, 2);
      x2= xCoord(obj2, 1);
      y2= xCoord(obj2, 2);
      k = indx(i,3);
      k_ref = indx(1,3);
%       
% if k == k_ref
%     line = "-k";
% else
%     line = "-r";
% end
plot([x1 x2], [y1 y2],line)
      
      hold on;
      
    
end



% a = bndryID(:,3);
% b = pbcID(:,3);

% 
% plot(xCoord(a,1),xCoord(a,2),"b--","linewidth",1);
% hold on
% plot(xCoord(b,1),xCoord(b,2),"r--","linewidth",1);
% hold on
%plot(xCoord(c,1),xCoord(c,2),"k--","linewidth",1);
% hold on
% plot(xCoord(e,1),xCoord(e,2),"g--","linewidth",1);
% hold on



% for n = 1:size(xCoord,1)
%     idxNumber = num2str(n);
%     text(xCoord(n,1),xCoord(n,2),idxNumber);
%     hold on
% end



axis equal;
axis off
end

