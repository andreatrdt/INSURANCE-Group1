function liabilities = Liabilities(C0, qx, lt, RD, COMM, discounts, expenses)
    % computes the liabilities
    % INPUT
    %

    V = zeros(size(time_grid));
    Expenses = zeros(size(time_grid));
    Contract_prob =cumprod([1; (1-qx(1:end-1)).*(1-lt(1:end-1))]);
    for i=1:T

        % cash flows at each year
        death_cf = (max(1.1*C0, F(i+1))-benefit_commission)*qx(i)*(1-lt(i));
    
            if i==T   % condition because in the last year it's mandatory to lapse
                lapse_cf = (F(i+1)-benefit_commission)*(1-qx(i));
            else
                lapse_cf = (F(i+1)-benefit_commission)*lt(i)*(1-qx(i));
            end
    
        Expenses(i)=Contract_prob(i)*expenses(i);
    
        V(i)=Contract_prob(i)*(lapse_cf+death_cf+expenses(i)+F(i+1)/(1-RD)*COMM);
    end
    
    liabilities = V' * discounts;

end % function Liabilities