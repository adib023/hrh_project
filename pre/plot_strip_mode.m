function plot_strip_mode(xCoord, indx, modeshape, frequnecy, stripNode, problemType,freqID)

xCoord =  xCoord(stripNode,1:2);
originalIndx = indx(:,1:2);


check =  ismember(originalIndx, stripNode);
doRowSum = sum(check');
targetRows = find(doRowSum == 2);

stripIndx =  originalIndx(targetRows,:)

for p = 1:length(stripNode)
    for q = 1:size(stripIndx,1)
        loc = find(stripIndx(q,:) == stripNode(p));
        if length(loc) == 1                    
            stripIndx(q,loc) = p;
        end
    end
end

% plot_lattice(xCoord,stripIndx,'r-')

alpha =  2;

for p = freqID %1: length(frequnecy)
    mode =  modeshape(:,p);
    disp = []

    mode_xCoord = xCoord;

    mode_xCoord(:,1) = mode_xCoord(:,1) + alpha * mode(1:2:end);
    mode_xCoord(:,2) = mode_xCoord(:,2) + alpha * mode(2:2:end);
    
    disp =  sqrt(mode(1:2:end).^2 + mode(2:2:end).^2);
    % 
    % h1 = figure('Visible','off');
    % plot_lattice(xCoord, stripIndx, 'r-');
    % plot_lattice(mode_xCoord, stripIndx, 'k-');
    % 
    % if problemType == 2
    %     print('-r300','-dpng',strcat('../result/undef_strip_mode/MODE',num2str(p)))
    % else
    %     print('-r300','-dpng',strcat('../result/def_strip_mode/MODE',num2str(p)))
    % end
    % 
    % h2 = figure('Visible','on');
    plot(1:length(stripNode),disp,'md-')
    hold on
    %legend('x','y','disp')
    %title(num2str(frequnecy(p)))
    % if problemType == 2
    %     print('-r300','-dpng',strcat('../result/undef_strip_mode/', num2str(p)))
    % else
    %     print('-r300','-dpng',strcat('../result/def_strip_mode/', num2str(p)))
    % end
end



end

