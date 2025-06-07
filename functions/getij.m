
function [iindex,jindex] = getij(Nn,Nge,type)

if nargin == 2
iindex = reshape((1:Nn*Nge)',Nn,[]);
iindex = repmat(iindex,Nn,1);
iindex = iindex(:);
jindex = repmat(1:Nn*Nge,Nn,1);
jindex = jindex(:);
else
    switch lower(rmspaces(type))
        case 'l'
            iindex = 1:(Nn*(Nn+1)/2);
        case 'u'  
            % length of indices vectors
            Niv = (Nn*(Nn+1)/2);
            iindex = zeros(Niv,1);
            jindex = zeros(Niv,1);
            index = 0;
            for i = 1:Nn
                for j = i:Nn
                    index = index + 1;
                    iindex(index) = i;
                    jindex(index) = j;
                end
            end
            
            iindex = repmat(iindex,Nge,1);
            jindex = repmat(jindex,Nge,1);
            
            temp = 1:Nn:Nn*Nge;
            temp = repmat(temp-1,Niv,1);
            temp = temp(:);
            
            iindex = iindex + temp;
            jindex = jindex + temp;
            
            
    end
end

end