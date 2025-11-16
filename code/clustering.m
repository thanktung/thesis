function cluster_output = cluster_analysis(QRSm, threshold)
    m = sort(cell2mat(QRSm));

    if length(m)>=2
        Y = pdist(m'); 
        Z = linkage(Y);
 
        T = cluster(Z,'cutoff', threshold,'criterion','distance')'; 
                       
        [uT,oT] = unique(T,'first'); 
        [~,pos] = sort(oT); 
        cluster_output = cell(1,length(uT));
        
        for i = uT
            cluster_output{i} = sort(m(T==pos(i)));
        end
    else
        cluster_output = {m};
    end
end
