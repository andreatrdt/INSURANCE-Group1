function liabilities = Liabilities(F0, P_death, lt, RD, COMM, discounts, expenses,dt,F, benefit_commission,T)
    % computes the liabilities
    % INPUT
    %

    V = zeros(size(dt));
    Expenses = zeros(size(dt));
    Contract_prob =cumprod([1; (1-P_death(1:end-1)).*(1-lt(1:end-1))]);
    for i=1:length(dt)

        % cash flows at each year
        death_cf = (max(F0, F(i+1))-benefit_commission)*P_death(i)*(1-lt(i));
    
            if i==T   % condition because in the last year it's mandatory to lapse
                lapse_cf = (F(i+1)-benefit_commission)*(1-P_death(i));
            else
                lapse_cf = (F(i+1)-benefit_commission)*lt(i)*(1-P_death(i));
            end
    
        Expenses(i)=Contract_prob(i)*expenses(i);
    
        V(i)=Contract_prob(i)*(lapse_cf+death_cf+expenses(i)+F(i+1)/(1-RD)*COMM);
    end
    
    liabilities = V' * discounts;

end % function Liabilities