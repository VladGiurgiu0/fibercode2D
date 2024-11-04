function [yfit,gof]=fit_data(x,y,Model,kernel)
%%% fits a curve (x,y) defined by the model and give back fit, derivative
%%% and goodness of fit
    switch Model
        case "Poly1"
            opts=fitoptions('Method','LinearLeastSquares'); warning('off');
            [fitresult, gof] = fit(x,y,fittype( 'poly1' ),opts);
            p=coeffvalues(fitresult);
            yfit=p(1)*x.^1 + p(2);
            p_prime=polyder(p);
            %yfit_prime=p_prime(1);
        case "Poly2"
            opts=fitoptions('Method','LinearLeastSquares'); warning('off');
            [fitresult, gof] = fit(x,y,fittype( 'poly2' ),opts);
            p=coeffvalues(fitresult);
            yfit=p(1)*x.^2 + p(2)*x.^1 + p(3);
            p_prime=polyder(p);
            %yfit_prime=p_prime(1)*x.^1 + p_prime(2);
        case "Poly3"
            opts=fitoptions('Method','LinearLeastSquares'); warning('off');
            [fitresult, gof] = fit(x,y,fittype( 'poly3' ),opts);
            p=coeffvalues(fitresult);
            yfit=p(1)*x.^3 + p(2)*x.^2 + p(3)*x.^1 + p(4);
            p_prime=polyder(p);
%             if p_prime==0
%                 yfit_prime=0;
%             else
%                 yfit_prime=p_prime(1)*x.^2 + p_prime(2)*x.^1 + p_prime(3);
%             end
        case "Poly4"
            opts=fitoptions('Method','LinearLeastSquares'); warning('off');
            [fitresult, gof] = fit(x,y,fittype( 'poly4' ),opts);
            p=coeffvalues(fitresult);
            yfit=p(1)*x.^4 + p(2)*x.^3 + p(3)*x.^2 + p(4)*x.^1 + p(5);
            p_prime=polyder(p);
            %yfit_prime=p_prime(1)*x.^3 + p_prime(2)*x.^2 + p_prime(3)*x.^1 + p_prime(4);
        case "Poly5"
            opts=fitoptions('Method','LinearLeastSquares'); warning('off');
            [fitresult, gof] = fit(x,y,fittype( 'poly5' ),opts); 
            p=coeffvalues(fitresult);
            yfit=p(1)*x.^5 + p(2)*x.^4 + p(3)*x.^3 + p(4)*x.^2 + p(5)*x.^1 + p(6);
            p_prime=polyder(p);
            %yfit_prime=p_prime(1)*x.^4 + p_prime(2)*x.^3 + p_prime(3)*x.^2 + p_prime(4)*x.^1 + p_prime(5);
       case "Poly9"
            opts=fitoptions('Method','LinearLeastSquares'); warning('off');
            [fitresult, gof] = fit(x,y,fittype( 'poly9' ),opts);
            p=coeffvalues(fitresult);
            yfit=p(1)*x.^9 + p(2)*x.^8 + p(3)*x.^7 + p(4)*x.^6 + p(5)*x.^5 + p(6)*x.^4 + p(7)*x.^3 + p(8)*x.^2 + p(9)*x.^1 + p(10);
            p_prime=polyder(p);
            %yfit_prime=p_prime(1)*x.^7 + p_prime(2)*x.^6 + p_prime(3)*x.^5 + p_prime(4)*x.^4 + p_prime(5)*x.^3 + p_prime(6)*x.^2 + p_prime(7)*x.^3;
       case "Marco Sgolay2 - 1" %%% Savitzky-Golay filtering - order 2 length kern_L
%         yfit = sgolayfilt(y,2,min(floor(numel(y)/2)*2-1,kernel));
%         yfit = sgolayfilt(yfit,2,min(floor(numel(yfit)/2)*2-1,kernel));
        gof=[];
        yfit = sgolayfilt(y,2,round_odd(floor(numel(y)/4)));
        yfit = sgolayfilt(yfit,2,round_odd(floor(numel(y)/4)));
        case "Marco Smooth 2x" %%% Savitzky-Golay filtering - order 2 length kern_L
        yfit = smooth(x',y,min(floor(numel(y)/2)*2-1,kernel),'sgolay',2);
        yfit = smooth(x',yfit,min(floor(numel(y)/2)*2-1,kernel),'sgolay',2);

        gof=[];

        case "Automatic"
           Model_Auto=["Poly1"; "Poly2"; "Poly3"; "Poly4";"Poly5"];
           for i=1:size(Model_Auto)
                [yfits(i,:),gof(i)]=fit_data(x,y,Model_Auto(i),kernel);
                sse(i)=gof(i).sse;
           end
           [~,idx_min]=min(sse);
           yfit=yfits(idx_min,:)';

        case 'movmean'
            yfit=smoothdata(y,1,Model,kernel);
            gof=[];
        case 'movmedian'
            yfit=smoothdata(y,1,Model,kernel);
            gof=[];
        case 'gaussian'
            yfit=smoothdata(y,1,Model,kernel);
            gof=[];
        case 'lowess'
            yfit=smoothdata(y,1,Model,kernel);
            gof=[];
        case 'loess'
            yfit=smoothdata(y,1,Model,kernel);
            gof=[];
        case 'rlowess'
            yfit=smoothdata(y,1,Model,kernel);
            gof=[];
        case 'rloess'
            yfit=smoothdata(y,1,Model,kernel);
            gof=[];
        case 'sgolay'
            yfit=smoothdata(y,1,Model,kernel);
            gof=[];
        case 'spline'
            [spline, gof] = fit(x,y,'smoothingspline');
            yfit=spline(x);
    end
end