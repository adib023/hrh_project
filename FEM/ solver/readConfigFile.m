function [meshFile,pbc,a1,a2,a3,nFreq] = readConfigFile(configFile)

a1 = [0 0 0] ; a2 = [0 0 0] ; a3 = [0 0 0] ; 

fid = fopen(configFile); 

dumy = fgetl(fid) ;  
while (length(dumy) == 0) ||  (strcmp(dumy(1),'%') == 1) 
	dumy = fgetl(fid) ;
end

dumy = fgetl(fid) ; 
while (length(dumy) == 0) ||  (strcmp(dumy(1),'%') == 1) 
	dumy = fgetl(fid) ;
end

meshFile = dumy 

dumy = fgetl(fid) ; 
while (length(dumy) == 0) ||  (strcmp(dumy(1),'%') == 1) 
	dumy = fgetl(fid) ;
end
str1 = split(dumy) ; 

for p = 1:3
	pbc(p) = str2num(str1{p}) ; 
end
pbc

dumy = fgetl(fid) ; 
while (length(dumy) == 0) ||  (strcmp(dumy(1),'%') == 1) 
	dumy = fgetl(fid) ;
end
str1 = split(dumy) ; 

if (pbc(1) == 1) && (pbc(2) == 0) && (pbc(3) == 0)
	for p = 1:3
		a1(p) = str2double(str1{p}) ; 
	end
elseif (pbc(1) == 1) && (pbc(2) == 1) && (pbc(3) == 0)
	for p = 1:3
		a1(p) = str2double(str1{p}) ; 
	end

	dumy = fgetl(fid) ; 
	while (length(dumy) == 0) ||  (strcmp(dumy(1),'%') == 1) 
		dumy = fgetl(fid) ;
	end
	str1 = split(dumy) ; 
	for p = 1:3
		a2(p) = str2double(str1{p}) ; 
	end
elseif  (pbc(1) == 1) && (pbc(2) == 1) && (pbc(3) == 1)
	error('not yet implemented');
end

dumy = fgetl(fid) ; 
while (length(dumy) == 0) ||  (strcmp(dumy(1),'%') == 1) 
	dumy = fgetl(fid) ;
end
nFreq = str2num(dumy); 

fclose(fid);

end
