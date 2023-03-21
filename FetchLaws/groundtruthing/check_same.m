function [result,Index] = check_same(old,new)
% checks if old and new 4D matrices are the same and finds indices which are different
% not using this function at the moment but could be used as a function call during checking
problem = 0;
for o = 1:25
    for p = 1:287
        A = old(:,:,o,p);
        B = new(:,:,o,p);
        Index{o,p} = find(A~=B);
        if ~isempty(Index{o,p})
            problem = problem + 1;
        end
%         surf(xplot',yplot',abs(nwn(:,:,o,p)-own(:,:,o,p)))
%         title("o "+o+"p "+p)
%         drawnow;
%         pause(0.1)
    end
end
result = problem;

end

