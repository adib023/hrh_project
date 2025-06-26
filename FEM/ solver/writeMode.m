function writeMode(alpha,UUR,freq,modeshape,fileName) 

fid = fopen(strcat('../Results/',fileName),'w');


mod_UUR =  zeros(size(UUR));

for p = 1:length(freq)
    mode = modeshape(:,p);
    % mode =  zeros(size(mode));

    count = 0;
    for n = 1:size(UUR,1)
        for d = 1:3
            count = count + 1;

            mod_UUR(n,d) = mode(count);
        end
    end

    fid = write(fid,mod_UUR);

end

fclose(fid);

end
