function [QRSonset] = QRSonsetdetect(QRSpos, x, fs)        
    [x_filter] = filter_user(x,fs);

    x_diff = diff(x_filter);
    x_ht = hilbert(x_diff);
    x_ht = real(x_ht);
 
    Q_loc = 0;
    for i = 1:length(QRSpos)
        if QRSpos(i)-0.12*fs <=0
            x_ht_dat = x_ht(1:QRSpos(i));
        else
            x_ht_dat = x_ht(QRSpos(i)-0.12*fs:QRSpos(i));
        end
                        
        [a,b] = max(x_ht_dat);
        Q_loc = [Q_loc, QRSpos(i)-0.12*fs + b-1];                       
    end    
    
    Q_loc(1) = [];   

    QRSon_uv = 0;
    for i = 1:length(Q_loc)
        if Q_loc(i)-0.12*fs <=0
            x_ht_dat = x_ht(1:Q_loc(i));
        else
            x_ht_dat = x_ht(Q_loc(i)-0.12*fs:Q_loc(i));
        end
        
        k = find(diff(x_ht_dat)<0,1,'last');
        QRSon_uv = [QRSon_uv, Q_loc(i)-0.12*fs + k]; 
    end
    
    where = find(QRSon_uv <= 0);
    QRSon_uv(where) = [];
    
    where = find(QRSon_uv >= length(x));
    QRSon_uv(where) = [];
    
    QRSonset = QRSon_uv;   
end


function [QRSoffset] = QRSoffsetdetect(QRSpos, x, fs)             
    [x_filter] = filter_user(x,fs);
    x_diff = diff(x_filter);    
    x_ht = hilbert(x_diff);
    
    x_ht = real(x_ht);
    x_ht_in = x_ht;
    
    S_loc = 0;
    for i = 1:length(QRSpos)
        if QRSpos(i)+0.12*fs >=length(x_ht)
            x_ht_dat = x_ht(QRSpos(i):length(x_ht));
        else
            x_ht_dat = x_ht(QRSpos(i):QRSpos(i)+0.12*fs);
        end
        
        [a,b] = min(x_ht_dat);
        S_loc = [S_loc, QRSpos(i) + b-1];                       
    end    
    
    where = find(S_loc <= 0);
    S_loc(where) = [];
    
    where = find(S_loc >= length(x));
    S_loc(where) = [];
    
    QRSoff_uv = 0;
    for i = 1:length(S_loc)
        if S_loc(i)+0.12*fs >=length(x_ht)
            x_ht_dat = x_ht(S_loc(i):length(x_ht));
        else
            x_ht_dat = x_ht(S_loc(i):S_loc(i)+0.12*fs);
        end
        
        [a,b] = max(x_ht_dat);
        QRSoff_uv = [QRSoff_uv, S_loc(i) + b-1];     
    end
           
    where = find(QRSoff_uv <= 0);
    QRSoff_uv(where) = [];
    
    where = find(QRSoff_uv >= length(x));
    QRSoff_uv(where) = [];

    QRSoffset = 0;
    for i = 1:length(QRSoff_uv)
        if QRSoff_uv(i)+0.12*fs >= length(x_ht)
            x_ht_dat = x_ht(QRSoff_uv(i):length(x_ht));
        else
            x_ht_dat = x_ht(QRSoff_uv(i):QRSoff_uv(i)+0.12*fs);
        end
        
        k = find(diff(x_ht_dat)>0,1,'first');
        QRSoffset = [QRSoffset, QRSoff_uv(i) + k];
    end
    
    QRSoffset(1) = [];  
    
    QRSoffset1 = 0;
    for i = 1:length(QRSoffset)
        x_filt_dat = x_filter(QRSoffset(i)-10:QRSoffset(i));
        [a b] = max(x_filt_dat);
        
        QRSoffset1 = [QRSoffset1 QRSoffset(i) - b];
    end
    
    where = find(QRSoffset1 <= 0);
    QRSoffset1(where) = [];
    
    where = find(QRSoffset1 >= length(x_ht_in));
    QRSoffset1(where) = [];
    
    QRSoffset = QRSoffset1;
end


function [Ponset] = Ponsetdetect(Ppos,x,fs)
    filt1 = medfilt1(x, 0.2*fs);
    filt2 = medfilt1(filt1, 0.6*fs);
    xin = x - filt2;

    [x_filter] = filter_user(x,fs);
    x_diff = diff(x_filter);    
    x_ht = hilbert(x_diff);
    
    x_ht = real(x_ht);
    x_ht_in = x_ht;    
    Ponset = 0;
	
    for i = 1:length(Ppos)
        if Ppos(i)-0.12*fs <= 0
            x_dat_Ponloc = x_ht_in(1:Ppos(i));        
        else
            x_dat_Ponloc = x_ht_in(Ppos(i)-0.12*fs:Ppos(i));        
        end
        
        [zero_loc, pks_val, pks_loc] = find_zero_extremes(x_dat_Ponloc);
        zero_loc = Ppos(i)-0.12*fs + zero_loc;
        elim = Ppos(i) - zero_loc(end);
        if elim <= 5 && length(zero_loc)>=2
            Ponset = [Ponset zero_loc(end-1)];
        else
            Ponset = [Ponset zero_loc(end)];
        end      
    end
    
    where = find(Ponset >= length(x));
    Ponset(where) = [];
    
    where = find(Ponset <= 0);
    Ponset(where) = [];    
end


function [Poffset] = Poffsetdetect(Ppos,x,fs)
    filt1 = medfilt1(x, 0.2*fs);
    filt2 = medfilt1(filt1, 0.6*fs);
    xin = x - filt2;

    [x_filter] = filter_user(x,fs);
    x_diff = diff(x_filter);    
    x_ht = hilbert(x_diff);
    
    x_ht = real(x_ht);
    x_ht_in = x_ht;

    Poffset = 0;

    for i = 1:length(Ppos)        
        x_dat_Ponloc = x_ht_in(Ppos(i):Ppos(i)+0.1*fs);        
        
        [zero_loc, pks_val, pks_loc] = find_zero_extremes(x_dat_Ponloc);
        zero_loc = Ppos(i) + zero_loc;
        elim = zero_loc(1) - Ppos(i);
        
        if elim <= 5 && length(zero_loc)>=2
            Poffset = [Poffset zero_loc(2)];
        else
            Poffset = [Poffset zero_loc(1)];
        end      
    end
    
    where = find(Poffset >= length(x));
    Poffset(where) = [];
    
    where = find(Poffset <= 0);
    Poffset(where) = [];    
end


function [Toffset] = Tenddetect(Tpos, QRS, x, fs)
    [x_filter] = filter_user(x,fs);
    x_diff = diff(x_filter);    
    x_ht = hilbert(x_diff);
    
    x_ht = real(x_ht);
    x_ht_in = x_ht;
    
    Toffset = 0;

    for i = 1:length(Tpos)-1
        x_T_end_dat = x_ht_in(Tpos(i):QRS(i+1));
        [zero_loc, pks_val, pks_loc] = find_zero_extremes(x_T_end_dat);
        zero_loc = zero_loc + Tpos(i) - 1;
        
        if length(zero_loc) > 1
            Toffset = [Toffset zero_loc(2)];
        else
            Toffset = [Toffset zero_loc(1)];
        end
    end
    
    x_T_end_dat = x_ht_in(Tpos(end):length(x_ht_in));
    [zero_loc, pks_val, pks_loc] = find_zero_extremes(x_T_end_dat);
    zero_loc = zero_loc + Tpos(end) - 1;
        
    if length(zero_loc) > 1
        Toffset = [Toffset zero_loc(2)];
    else
        Toffset = [Toffset zero_loc(1)];
    end

    where = find(Toffset <= 0);
    Toffset(where) = [];
    
    where = find(Toffset >= length(x_ht_in));
    Toffset(where) = [];
end