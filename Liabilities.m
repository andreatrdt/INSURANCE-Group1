function [liabilities, Lapse_BEL, Death_BEL, Expenses_BEL,Commissions_BEL] = Liabilities(F0, P_death, lt, RD, COMM, discounts, expenses,dt,F, benefit_commission,T)
    % computes the liabilities
    % INPUT
    % F0: initial fund
    % P_death: probability of death
    % lt: probability of lapse
    % RD: regular deduction
    % COMM: commission
    % discounts: discounts
    % expenses: expenses
    % dt: time steps
    % F: fund
    % benefit_commission: benefit commission
    % T: time to maturity
    %
    % OUTPUT
    % liabilities: liabilities of the insurance company


    Expenses = zeros(size(dt));
    Contract_prob = cumprod([1; (1-P_death(1:end-1)).*(1-lt(1:end-1))]);

    for i=1:length(dt)

        % cash flows at each year
        death_cf = (max(F0, F(:,i+1)))*P_death(i)*(1-lt(i));
    
            if i==T   % condition because in the last year it's mandatory to lapse
                lapse_cf = (F(:,i+1)-benefit_commission)*(1-P_death(i));
            else
                lapse_cf = (F(:,i+1)-benefit_commission)*lt(i)*(1-P_death(i));
            end
    
        Expenses(i)=Contract_prob(i)*expenses(i);
        Lapse_benefits(:,i)=Contract_prob(i)*lapse_cf;
        Death_benefits(:,i)=Contract_prob(i)*death_cf;
        Commissions(:,i)=Contract_prob(i)*F(:,i)/(1-RD)*COMM;

        Val(:,i) = Lapse_benefits(:,i) +  Death_benefits(:,i) + Expenses(i) + Commissions(:,i);

    end

    size(Lapse_benefits)
    size(Death_benefits)
    size(Expenses)
    
    % computation of the death benefits 
    Death_BEL=(mean(Death_benefits)*discounts);

 % computation of the lapse benefits
    Lapse_BEL=(mean(Lapse_benefits)*discounts);

    % computation of expenses
    Expenses_BEL=(mean(Expenses)*discounts);

    % compuation of commissions
    Commissions_BEL=(mean(Commissions)*discounts);

    liabilities = sum(mean(Val)* discounts);

end % function Liabilities