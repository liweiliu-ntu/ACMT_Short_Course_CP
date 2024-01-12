function [QQnext, alphanext, taunext,fnext] = flowRule(Qnow, Qbnow, taunow, Q_dot1, c, alpha2,kappa1, kappa2, a1, a2, err,dt,K,P)
    % Compute df_dalpha
    s=P*Qnow;
    nr=(s-alpha2)/sqrt((s-alpha2)'*(s-alpha2));
    df_dalpha = sqrt(3/8) * ((-2 + c.*(nr' * alpha2)) * nr + 2 * c * alpha2- c * s) / sqrt(1 - c.*(nr'* alpha2));
    df_dsigma = sqrt(3/2) * ((2 - c.*(nr' * alpha2)) * nr - c * alpha2) / sqrt(1 - c.*(nr'* alpha2));
    df_dsigma_norm=norm(df_dsigma);
    n=df_dsigma/df_dsigma_norm;
    % Compute Kp
    Kp = 0.5 * kappa1 * (1 - kappa2 * taunow) - df_dalpha' * (df_dsigma - a2 * sqrt(df_dsigma'*df_dsigma) * alpha2) * a1;

    dQtrail=Q_dot1*dt;

    Qtrail=Qnow+dQtrail;
    
    % Compute lambda
    lambda = df_dsigma' * dQtrail * (Kp + df_dsigma' * K* df_dsigma)^-1;

    % Update variables

    Qnow1 = Qnow + dQtrail - lambda * K * df_dsigma;

    snow=P*Qnow1;

%     epsilon_p = epsilon_p + sqrt(2/3) * lambda * df_dsigma_norm;

    alphanext = alpha2+ lambda * a1 * df_dsigma_norm * (n - a2 * alpha2);


    taunext = taunow + 0.5 * lambda * kappa1 * (1 - kappa2 * taunow);


     tau1=taunext;

     Qnext=Qnow1;

     alpha1=alphanext;

    [fnext] = Yield_function(Qnext,P,tau1,c, alpha1)  ;
    f1=fnext;

      if fnext < -err
        Qnext= Qtrail;

        [fnext] = Yield_function(Qnext,P,tau1,c, alpha1);
        f3=fnext;

     end

    if fnext > err
        % Update sigma using Newton's method
        QQnext = Qnow1 + (sqrt(3/2) * ( taunext  / sqrt(1 - c.*(nr'* alphanext)))- sqrt((snow-alphanext)'*P*(snow-alphanext)))*nr;
        
        Qnext=QQnext;
        % Yield function
        [fnext] = Yield_function(Qnext,P,tau1,c, alpha1) ;
        f2=fnext;

        if n'*df_dsigma < 0
            n=dQtrail;
        end

    end

        while fnext > err
          QQnext= Qnext-fnext.*n*(n'*df_dsigma)^(-1);
          Qnext=QQnext;
          [fnext] = Yield_function(Qnext,P,tau1,c, alpha1) ;
          f4=fnext;
     
        end
       

end





