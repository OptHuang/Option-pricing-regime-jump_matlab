function x_star = Newton_perMerton(max_iter,eps,sigma,r,d,lam,kappaM,mu,delt)
    
   x0 = -1;
   [y,deriv_y] = Merton_func(sigma,r,d,lam,kappaM,mu,delt,x0);
   iter = 0;
   while((iter<max_iter)&&(abs(y)>eps))
       x0 = x0-(1/deriv_y)*y;
       [y,deriv_y] = Merton_func(sigma,r,d,lam,kappaM,mu,delt,x0);
       iter = iter+1;
   end

   x_star = x0;
    
end