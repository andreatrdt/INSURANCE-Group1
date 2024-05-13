function [liabilities, Lapse_BEL, Death_BEL, Expenses_BEL,Commissions_BEL] = Liabilities(F0, P_death, lt, regular_deduction, COMM, discounts, expenses,dt,F, lapse_pay,T)
    % Computes the liabilities and split the BEL value into its main PV components: 
    % death benefits, lapse benefits, expenses and commissions
    %       
    % INPUTS
    % F0: initial fund
    % P_death: probability of death
    % lt: probability of lapse
    % RD: regular deduction
    % COMM: commission
    % discounts: discounts
    % expenses: expenses per year
    % dt: time steps
    % F: fund value
    % benefit_commission: benefit commission
    % T: time to maturity
    %
    % OUTPUTS
    % liabilities: liabilities of the insurance company
    % main PV components of the BEL value:
    % Lapse_BEL 
    % Death_BEL
    % Expenses_BEL
    % Commissions_BEL

    % Initialize variables
    Expenses = zeros(size(dt));
    Contract_prob = cumprod([1; (1-P_death(1:end-1)).*(1-lt(1:end-1))]);
    
    % Initialize waitbar
    h = waitbar(0, 'Calculating...');

    for i = 1:length(dt)
        % cash flows at each year
        death_cf = (max(F0, F(:,i+1))) * P_death(i);     
        
        % at the end of the 50 years all the people leave the contract with a massive surrender
        if i == T   % condition because in the last year it's mandatory to lapse
            lapse_cf = (F(:,i+1) - lapse_pay) * (1 - P_death(i));
        else
            lapse_cf = (F(:,i+1) - lapse_pay) * lt(i) * (1 - P_death(i));
        end

        Expenses(i) = Contract_prob(i) * expenses(i);
        Lapse_benefits(:,i) = Contract_prob(i) * lapse_cf;
        Death_benefits(:,i) = Contract_prob(i) * death_cf;
        Commissions(:,i) =  COMM * Contract_prob(i) * F(:,i+1) / (1 - regular_deduction);

        Val(:,i) = Lapse_benefits(:,i) +  Death_benefits(:,i) + Expenses(i) + Commissions(:,i);
        
        % Update waitbar
        waitbar(i / length(dt), h, sprintf('Calculating... %d%%', round(i / length(dt) * 100)));
    end

    % Close waitbar after loop completes
    close(h);

    % computation of the death_BEL 
    Death_BEL=(mean(Death_benefits)*discounts);
     
    % computation of the lapse_BEL
    Lapse_BEL=(mean(Lapse_benefits)*discounts);

    % computation of expenses_BEL
    Expenses_BEL=(Expenses'*discounts);

    % computation of commissions_BEL
    Commissions_BEL=(mean(Commissions)*discounts);
    
    % computation of the liabilities
    liabilities = sum(mean(Val)* discounts);

    % figure
    % bar ( mean(Val)' .* discounts )
    % xlabel (' years ')
    % title ('Discounted cash flows')

    end % function Liabilities