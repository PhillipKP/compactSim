function [cosOfAng] = computeCosine(a, b )

% This function takes vector column vector a, normalizes a so it's L2 Norm = 1
% Then it takes matrix b, and normalizes each of the columns so each of the
% respective columns = 1


    
    aNorm = a / norm(a);

    % Normalize the columns of the matrix b
    bNorm = bsxfun( @rdivide, b, sqrt(diag(b'*b))' );

    cosOfAng = bNorm.' * aNorm;

    cosOfAng = cosOfAng.';



    

    
    
end