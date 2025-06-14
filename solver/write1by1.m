
function [] =  write1by1(fname, data_to_be_stored)
  
  [nr,nc] = size (data_to_be_stored);

  file = fopen(fname, "a");
  
  for i = 1: nr
    for j = 1 : nc
      fprintf(file,"%4.10f\t", data_to_be_stored(i,j));
    end
    fprintf(file,"\n");
  end
  %fprintf(file,"\n");
  
  fclose ( file );

end