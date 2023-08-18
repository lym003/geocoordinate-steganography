% The significant digit of geo-coordinate is more than 8
% construct cover element sequences
% s_p is the geo-coordinate sequence; point_num is the number of digits; msg is the secret messages after encrypting
% dig_num is used to store the digit number of each part in the geo-coordinate
    stc_len = 1;
    lsb_len = fra_len - 7;
    dig_num(,1) = int_len; % int_len is the digit number of integer part
    dig_num(,2) = fra_len; % fra_len is the digit number of fraction part
    stc_len = 1; % stc_len is the digit number of SiSDPs
    dig_num(,3) = stc_len; 
    stc_sum = stc_sum + stc_len;
    lsb_len = fra_len - 7; % lsb_len is the digit number of LSDPs
    dig_num(,4) = lsb_len;
    lsb_sum = lsb_sum + lsb_len;

% the embedding process
% the MLSBR algorithm
    dec_str = int2str(bin2dec(msg(m_start:m_start+52)));
    dec_str_len = length(dec_str);
    if dec_str_len<16   %53位二进制表示十进制位数不足16位的，前面补足0
        add_len = 16 - dec_str_len;
        for i=1:add_len
            dec_str = ['0',dec_str];
        end
    end
    lsb_embstr = [lsb_embstr,dec_str];    %将lsb每一组载密（16位）合成总的字符串
        
    lsb_remain = mod(lsb_sum,16);  %16个十进制一组嵌入后，lsb部分不足一组的十进制位数
    lsb_remain_bin = floor(log2(10^lsb_remain-1));
    lsb_remain_dec = int2str(bin2dec(msg(m_start:m_start+lsb_remain_bin-1)));
    if length(lsb_remain_dec)<lsb_remain
        add_len = lsb_remain - length(lsb_remain_dec);
        for i=1:add_len
            lsb_remain_dec = ['0',lsb_remain_dec];
        end
    end
    lsb_embstr = [lsb_embstr,lsb_remain_dec];

% the STC adaptive algorithm
% bp is the embedding ratio;
    emb_len = stc_emb_len;
    [cover,costs] = extra_c(p,dig_num,cover_p,row_num,col_num,stc_sum); % extract cover sequence of SiSDPs
    [d stego n_msg_bits l] = stc_pm1_pls_embed(cover, costs, stc_m, h); % embed message with STC algorithm
    for i=1:stc_sum
        stc_embstr = [stc_embstr,int2str(stego(i))];
    end

% update geo-coordinate after embedding process
stego_str = [lsb_embstr,stc_embstr];
stego_len = length(stego_str);
while point<=stego_len
    point = 1;  %载密元素（lsb+stc）序列中元素序号
    for i=1:row_num
        for j=1:col_num
            coord_str = s_p{i,j};
            coord_str(end_plane-lay+1) = stego_str(point);  % end_plane is the last digit plane of coord_str; lay is the number of RDPs
            s_p{i,j} = coord_str;
            point = point + 1;
        end
    end
    lay = lay + 1;
end
        
% the extracting process
% stego_point is the processed point sequences from geographic data file
% extracting messages from LSDPs
k = 1; 
n = 1;  
while k <= lay_mlsb  % lay_mlsb is the number of LSDPs
    for i=1:row_num
        for j=1:col_num
            ext_point = stego_point{i,j};
            dec_msg(n) = ext_point(end_plane-k+1);
            n = n + 1;
        end
    end
    k = k + 1;
end

% convert each 16-digit to 53-bit
dec_point = 1;  
while length(dec_msg(dec_point:))>15
    dec = str2num(dec_msg(dec_point:dec_point+15));
    bin_str = dec2bin(dec,53);
    str_lsb = [str_lsb,bin_str];
    dec_point = dec_point + 16;
end
bin_str = dec2bin(str2num(dec_msg(dec_point:)));
str_lsb = [str_lsb,bin_str];

% extracting messages from SiSDPs
n = 1;
for i=1:row_num
    for j=1:col_num
        stego_stc(n) = stego_point{i,j}(stc_plane) % stc_plane is the number of the SiSDPs
        n = n +1;
    end
end
bin_stcmsg = stc_ml_extract(stego_stc, n_msg_bits, h); 
str_lsb = [str_lsb,bin_stcmsg];

% extracting elements from SiSDPs and computing distortion cost
function [cover,costs] = extra_c(cell_p,dig_num,p,row_num,col_num,stc_sum)
    s = 1;
    for i=1:row_num
        for j=1:col_num
            value = str2double(cell_p{i,j}(stc_plane)); % stc_plane is the number of the SiSDPs
            cover(1,s) = int32(value); 
            if value==0 || value==9    
                costs(:,s) = [1e+5 0 1e+5];
            else
                costs(1,s) = emb_distortion(p,7,i,j,-1);   %+1嵌入失真代价
                costs(3,s) = emb_distortion(p,7,i,j,1);   %-1嵌入失真代价
            end
            s = s + 1;
        end
    end
end
       
% emb_distortion is the distortion function
function cost = emb_distortion(p,stc_k,r_num,c_num,sg)
    bp = 0.8;    % bp is an adjustment parameter
    if(p(r_num,c_num)<0)
        abs_p = abs(p(r_num,c_num));
        sign = -1;
    else
        abs_p = p(r_num,c_num);
        sign = 1;
    end
    new_coor = sign*(abs_p + sg * 10^(-stc_k)); % the new geo-coordinate after "±1" embedding
    b = p(r_num-1,:) - p(r_num,:);
    c = p(r_num-1,:) - point_emb;  % point_emb is the new point by updating new_coor
    sigma_front = acos(dot(b,c)/sqrt(sum(b.^2)*sum(c.^2))*bp); 
    d = p(r_num+1,:) - p(r_num,:);
    e = p(r_num+1,:) - point_emb;
    sigma_back = acos(dot(d,e)/sqrt(sum(d.^2)*sum(e.^2))*bp); 
    cost = rad2deg(sigma_front + sigma_back);
end
