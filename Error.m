function [ error_1 ] = Error(U_Mean,N_m )
error_1=zeros(1,1);

%%Calculate the relative error value
    for i=1:4
        for j=1:N_m-1
            error_1(i,j)=U_Mean(i,j+1)-U_Mean(i,j);
        end
    end
end