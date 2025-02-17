function write(a,name)
        fid   = fopen(name,'w');
        [m,n] = size(a);
     for i=1:m
        for j=1:n
            if j==n
               fprintf(fid,'%g\n',a(i,j));
            else
               fprintf(fid,'%g\t',a(i,j));
            end
        end    
     end
    fclose(fid);
    judge    = exist('Results');
    if judge ~= 7
        system('mkdir Results');
    end
    file_path   = strcat(cd,'\Results');
    movefile(name,file_path);
end