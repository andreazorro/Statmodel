function [sNet] = Statmodel(genes, regulators, expressiondata)

    expressiondata = expressiondata';

    distribution = {'Beta','Binomial','BirnbaumSaunders','Burr','Exponential',...
    'ExtremeValue','Gamma','GeneralizedExtremeValue','GeneralizedPareto',...
    'HalfNormal','InverseGaussian','Kernel','Logistic','Loglogistic',...
    'Lognormal','Nakagami','NegativeBinomial','Normal','Poisson','Rayleigh',...
    'Rician','Stable','tLocationScale','Weibull'};

    nr = size(regulators,1);
    tfs = zeros(1,nr);

    parfor i = 1:nr
        try
            tfs(1,i) = find(strcmp(genes,regulators{i}));
        catch
            continue
        end
    end

    tfs = tfs(tfs ~= 0);

    ngenes = size(expressiondata,1);
    ntf = size(tfs,2);
    ndist = size(distribution,2);

    pv = zeros(ngenes,ntf);

    for j = 1:ngenes
        
        Y = expressiondata(j,:)';
        meanY = mean(Y);
        maxY = max(Y);
        minY = min(Y);
        stY = (Y - meanY)/(maxY - minY); 
        stYm = stY+abs(min(stY))+1;
        
        parfor k = 1:ndist    
            warning ('off','all')
            try 
                pd = fitdist(stYm,distribution{k}); 
                [h(k),p(k)] = chi2gof(stYm,'CDF',pd); 
            catch
                h(k)=-Inf;
                p(k)=-Inf; 
            end
        end
        
        if find(h==0) == 1
            bestfit = 1;
        else
            [~,bestfit] = max(p);
        end
        
        try 
            bestpdY = fitdist(stYm,distribution{bestfit}); 
        catch
            fprintf('Error fitting gene #%d \n',j); 
        end
        
        parfor n = 1:ntf
            
            if tfs(1,n) == j 
                continue
            else
                X = expressiondata(tfs(n),:)';
                meanX = mean(X);           
                maxX = max(X); 
                minX = min(X); 
                stX = (X - meanX)/(maxX - minX);         
                
                Xmp = find(stX >= mean(stX));
                Xmn = find(stX <= mean(stX));
            
                stYp = stYm(Xmp);
                stYn = stYm(Xmn); 

                if size(stYp)<=size(stYn) 
                    stYt=stYp; 
                else
                    stYt=stYn; 
                end
            
                if abs(mean(stYp) - mean(stYn)) >= 0.05 
                    [h0,p0] = chi2gof(stYt,'CDF',bestpdY); 
                    if  h0 == 1 
                        pv(j,n) = p0; 
                    else
                        continue 
                    end
                end
            end
        end
    end

    r = 1;
    Net = cell(nnz(pv),3);

    for i = 1:ntf 
        for j = 1:ngenes 
            if pv(j,i) == 0 
                continue
            else
                Net{r,1} = genes{tfs(i)};
                Net{r,2} = genes{j}; 
                Net{r,3} = -log(pv(j,i)); 
                r = r+1; 
            end
        end
    end

    sNet = sortrows(Net,3,'descend');
    
 end


