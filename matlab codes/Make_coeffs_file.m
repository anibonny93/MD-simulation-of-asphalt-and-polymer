path='PET.lmp';




files=dir(path);
files=files(3:end);

for i=1:length(files)
    f=fopen([path '/' files(i).name],'r')
    fout=fopen(['coeffs_' files(i).name],'w')
    testflag=1
    while testflag==1
        x=fgets(f);
        if isnumeric(x)
            break
        end



        if startsWith(x,'Bond Coeffs')
            xtest=fgets(f);
            testflag2=1;
            while testflag2==1
                x=fgets(f);
                if ~strcmp(x,xtest)
                    fwrite(fout,['bond_coeff ' x]);
                elseif isnumeric(x)
                    testflag2=-1;
                else
                    break
                end
            end
        end

        if startsWith(x,'Pair Coeffs')
            xtest=fgets(f);
            testflag2=1;
            while testflag2==1
                x=fgets(f);
                if ~strcmp(x,xtest)
                    fwrite(fout,['pair_coeff ' x]);
                elseif isnumeric(x)
                    testflag2=-1;
                else
                    break
                end
            end
        end

        if startsWith(x,'Angle Coeffs')
            xtest=fgets(f);
            testflag2=1;
            while testflag2==1
                x=fgets(f);
                if ~strcmp(x,xtest)
                    fwrite(fout,['angle_coeff ' x]);
                elseif isnumeric(x)
                    testflag2=-1;
                else
                    break
                end
            end
        end

        if startsWith(x,'Dihedral Coeffs')
            xtest=fgets(f);
            testflag2=1;
            while testflag2==1
                x=fgets(f);
                if ~strcmp(x,xtest)
                    fwrite(fout,['dihedral_coeff ' x]);
                elseif isnumeric(x)
                    testflag2=-1;
                else
                    break
                end
            end
        end


        if startsWith(x,'Improper Coeffs')
            xtest=fgets(f);
            testflag2=1;
            while testflag2==1
                x=fgets(f);
                if ~strcmp(x,xtest)
                    fwrite(fout,['improper_coeff ' x]);
                elseif isnumeric(x)
                    testflag2=-1;
                else
                    break
                end
            end
        end


    end
end


