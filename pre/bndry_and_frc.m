function [disp_bndry, frcID] = bndry_and_frc(N1,N2,ucID, TB_disp)

%%%%% TB_disp = [top_disp;bottom_disp]
node =  [ucID(1:N1,N2,3);ucID(1:N1,N2,4)];

pres_disp = TB_disp(1,:);
topDispBndry(:,1) =  node;
topDispBndry(:,2) =  pres_disp(1);
topDispBndry(:,3) =  pres_disp(2);

pres_disp = TB_disp(2,:);

node =  [ucID(1:N1,1,1);ucID(1:N1,1,2)];

bottomDispBndry(:,1) =  node;
bottomDispBndry(:,2) =  pres_disp(1);
bottomDispBndry(:,3) =  pres_disp(2);

disp_bndry = [bottomDispBndry;topDispBndry];


% frcID = [];
% for p = fix(N1/2)+1
%     frcID = [frcID ucID(p,fix(N2/2)+1,1) ucID(p,fix(N2/2)+1,2)];
% end

frcID = ucID(1,fix(N2/2)+1,1);




end