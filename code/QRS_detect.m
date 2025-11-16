function [QRS] = QRS_detect(data,fs)    
    extent = fs*0.12;

    beginning(1:extent) = data(1);          
    ending(1:extent) = data(end);   

    data = [beginning data ending];   
    [filt_dat,diff_dat,mov_dat] = usr_fil(data,fs,round(fs*0.12));

    mov_dat = mov_dat(extent+1 : end-extent);
	
    [FP_pks,FP_lcs] = findpeaks(mov_dat);
    [FP_max FP_lc] = max(FP_pks);
    FP_max_in_mov = find(mov_dat == FP_max);
	
    kc_to_FP_max_lc = 0;
    cc_to_FP_max_vl = 0;

    for i = 1:length(FP_lcs)
        kc_to_FP_max_lc = [kc_to_FP_max_lc abs(FP_lcs(i) - FP_max_in_mov)];
        cc_to_FP_max_vl = [cc_to_FP_max_vl FP_max - FP_pks(i)];        
    end
    
    where = find(kc_to_FP_max_lc==0);    
    kc_to_FP_max_lc(where) = eps;   
    kc_to_FP_max_lc(1) = [];
    
    where = find(cc_to_FP_max_vl==0);
    cc_to_FP_max_vl(where) = eps; 
    cc_to_FP_max_vl(1) = []; 
           
    agl = cc_to_FP_max_vl./kc_to_FP_max_lc;  
	
    [s_y,s_x] = findpeaks(-agl);
	
    s_x = [1 s_x length(FP_lcs)];   
    s_x = sort([FP_lcs(s_x) FP_lcs(FP_lc)]);               
        
    [tf,loc] = ismember(s_x, FP_lcs);    
    k3 = 0;
   
    for i = 1:length(loc)-1
        k = mov_dat(s_x(i)+1 : s_x(i+1)-1);
        thres = 0.55*std(k,1);     
        k1 = FP_lcs(loc(i)+1:loc(i+1)-1);
        mov_dat = mov_dat - mean(mov_dat);
        k2 = find(mov_dat(k1)>=thres);
        k3 = [k3 k1(k2)];
    end    
    
    QRS = unique(sort([s_x k3(2:end)]));
    
    elim = diff(QRS);             
    where = find(elim < fs*0.12);       
    QRS(where+1) = [];          
    
    ii = 0;
    for i = 2:length(QRS)-1
        if mov_dat(QRS(i)) < 0.97*mean(mov_dat(QRS(i-1):QRS(i+1)))         
            ii = [ii,i];
        end
    end

	QRS(ii(2:end)) = [];                
    
    QRS = real_r_peak_detection(mov_dat, fs, QRS, 0.2*fs);
end