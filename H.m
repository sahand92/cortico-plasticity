function Htau=H(tau)
    if tau > 0
        Htau=exp(-tau);
    elseif tau <= 0
        Htau=-exp(tau);
    end
    
end
