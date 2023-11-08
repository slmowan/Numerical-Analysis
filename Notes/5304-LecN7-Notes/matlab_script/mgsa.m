% modified Gram-Schmidt script
function [Q,R] = mgsa (A)
    % [Q,R] = cgsa (A)
    % modified Gram Schmidt QR factorization of A
     [m,n] = size(A);
     for j=1:n
         q = A(:,j); 
         for i=1:j-1 
             r = q'*Q(:,i);         
             q = q - r*Q(:,i);
             %fprintf(1,'----- j = %d i = %d\n',j,i)
             %disp('------')
             %disp(r); disp(q);          pause
             R(i,j) = r;
         end
    %%---------- error exit for case rjj == 0
         r = norm(q) ;
         if (r==0.0), error('** zero column'), end
         Q(:,j) = q / r; 
         R(j,j) = r; 
    %%     Q
    %%     R
    %%     pause 
end