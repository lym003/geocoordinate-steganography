%stc嵌入小数点后第7位，第8位至最后用mlsb算法

% clc;
% clear;

function ed = stc_lsb_ed(s_digit,p,msg)

switch s_digit %坐标的有效数字决定最大嵌入比率；
    case '9'
        bp = 0.51;
    case '10'
        bp = 1.00;
    case '11'
        bp = 1.00;
    case '12'
        bp = 1.00;
    case '13'
        bp = 1.00;
    case '14'
        bp = 0.99;
    case '15'
        bp = 0.99;
    case '16'
        bp = 1.00;
end
% if s_digit
%     bp = 0.50;
% end

[row_num,col_num] = size(p);
s_p = p;    %cell类型stego坐标点
cover_p = zeros(row_num,col_num);
[~,m_len] = size(msg);
h = 10;
dig_num = zeros(row_num*col_num,4);
s = 1;  %坐标点数量
lay_lsb = 0;   %载密矩阵中lsb的列数
lsb_sum = 0;    %lsb位平面总和
for i=1:row_num
    for j=1:col_num
        point_str = s_p{i,j};
        point_len = length(point_str);
        point_dec = str2double(s_p{i,j});
        cover_p(i,j) = str2double(s_p{i,j});
        if(point_dec < 0)
           point_num = point_len - 2;  %有效位数减去正负号和小数点
           point_abs = abs(point_dec);
        else
           point_num = point_len - 1;
           point_abs = point_dec;
        end
        int_len = length(mat2str(floor(point_abs)));
        fra_len = point_num - int_len;
        if fra_len - 6>0
            stc_len = 1;
            lsb_len = fra_len - 7;
        else
            stc_len = 0;
            lsb_len = 0;
        end
        dig_num(s,1) = int_len;
        dig_num(s,2) = fra_len;
        dig_num(s,3) = stc_len;
        dig_num(s,4) = lsb_len;
        if lay_lsb<lsb_len
            lay_lsb = lsb_len;
        end
        lsb_sum = lsb_sum + lsb_len;
        s = s + 1;
    end
end

m_start = 1;    %秘密消息嵌入移动指针
% lay_back = 1;  %逆序位平面
% s_num = 1;  %第s_num个载体
stc_sum = row_num*col_num;    %stc位平面总和
% flag = 1;
%每次用53位二进制数，最高可表示16位十进制数。将十进制按高位到低位替换至lsb位平面。
r = floor(lsb_sum/16);  %得到lsb总位平面含有多少组16个十进制；
lsb_remain = mod(lsb_sum,16);  %16个十进制一组嵌入后，lsb部分不足一组的十进制位数
lsb_embstr = '';
flag = 1;
while flag&&m_start<m_len
    if r
        if m_start+52<=m_len
            dec_str = int2str(bin2dec(msg(m_start:m_start+52)));
            m_start = m_start + 53;
        else
            dec_str = int2str(bin2dec(msg(m_start:m_len)));
            m_start = m_len + 1;
        end
        dec_str_len = length(dec_str);
        if dec_str_len<16   %53位二进制表示十进制位数不足16位的，前面补足0
            add_len = 16 - dec_str_len;
            for i=1:add_len
                dec_str = ['0',dec_str];
            end
        end
        lsb_embstr = [lsb_embstr,dec_str];    %将lsb每一组载密（16位）合成总的字符串
        r = r - 1;
    elseif r==0&&length(lsb_embstr)<lsb_sum
        lsb_remain_bin = floor(log2(10^lsb_remain-1));
        if m_start+lsb_remain_bin-1>m_len
            lsb_remain_bin = m_len-m_start+1;
            lsb_remain_dec = int2str(bin2dec(msg(m_start:m_len)));
        else
            lsb_remain_dec = int2str(bin2dec(msg(m_start:m_start+lsb_remain_bin-1)));
        end
        if length(lsb_remain_dec)<lsb_remain
            add_len = lsb_remain - length(lsb_remain_dec);
            for i=1:add_len
                lsb_remain_dec = ['0',lsb_remain_dec];
            end
        end
        lsb_embstr = [lsb_embstr,lsb_remain_dec];
        m_start = m_start + lsb_remain_bin;
    else  %在stc部分嵌入
        stc_emb_len = floor(stc_sum*bp);   %stc部分可以嵌入的秘密信息长度
        if m_start+stc_emb_len-1<=m_len
            emb_len = stc_emb_len;
        else
            emb_len = m_len-m_start+1;  %剩余位数小于可嵌秘密信息长度
        end
        [cover,costs] = extra_c(p,dig_num,cover_p,row_num,col_num,stc_sum);   %提取stc部分载体
        for j=1:emb_len
            stc_m(j) = uint8(str2double(msg(m_start+j-1)));
        end
        [d stego n_msg_bits l] = stc_pm1_pls_embed(cover, costs, stc_m, h); % embed message
        for j=1:emb_len
            if stego(j)<0 || stego(j)>9
                fprintf('stego超出0-9！\n');
            end
        end
        m_start = m_start + emb_len;
        extr_msg(m_start-emb_len:m_start-1) = stc_ml_extract(stego, n_msg_bits, h); % extract message
        flag = 0;
    end
end
%fprintf('实际嵌入%d位，',m_start-1);
fprintf('单位坐标嵌入%.2f位，',(m_start-1)/row_num/col_num);

%更新载密坐标点
stc_embstr = '';
if exist('stego','var')
    for i=1:stc_sum
        stc_embstr = [stc_embstr,int2str(stego(i))];
    end
end
stego_str = [lsb_embstr,stc_embstr];
stego_len = length(stego_str);
point = 1;  %载密元素（lsb+stc）序列中元素序号
lay = 1;
w = 1;
while point<=stego_len&&lay<=lay_lsb
    for i=1:row_num
        for j=1:col_num
            if point<=stego_len&&lay<=dig_num(w,4)
                coord_str = s_p{i,j};
                coord_str(end-lay+1) = stego_str(point);
                s_p{i,j} = coord_str;
                point = point + 1;
                w = w + 1;
            elseif lay>dig_num(w,4)
                w = w + 1;
            else 
                break;
            end
        end
    end
    w = 1;
    lay = lay + 1;
end
if point<=stego_len %更新stc部分
    w = 1;
    for i=1:row_num
        for j=1:col_num
            if dig_num(w,3) == 1
                coord_str = s_p{i,j};
                coord_str(end-dig_num(w,4)) = stego_str(point);
                s_p{i,j} = coord_str;
            end
            point = point + 1;
            w = w + 1;
        end
    end
end
stego_p = str2double(s_p);

%提取lsb部分嵌入的秘密消息
ext_len = 1;  %提取秘密消息长度
str_lsb = '';
dec_msg = '';
k = 1;  %位平面序号
n = 1;  %lsb提取载密序号
c = 1;  %载体序号
while k<=lay_lsb&&n<=stego_len
    for i=1:row_num
        for j=1:col_num
            if k<=dig_num(c,4)&&n<=stego_len
                ext_point = s_p{i,j};
                dec_msg(n) = ext_point(end-k+1);
                n = n + 1;
                c = c + 1;
            elseif n<=stego_len
                c = c + 1;
            else 
                 break;   
            end
        end
    end
    c = 1;
    k = k + 1;
end
dec_point = 1;  %提取的lsb部分十进制数序列，每16位转为53位二进制。
while dec_point<n
    if dec_point+15<n&&ext_len+52<m_start
        dec = str2num(dec_msg(dec_point:dec_point+15));
        bin_str = dec2bin(dec,53);
        str_lsb = [str_lsb,bin_str];
        dec_point = dec_point + 16;
        if ~strcmp(str_lsb(ext_len:ext_len+52),msg(ext_len:ext_len+52))
            fprintf('本段秘密消息提取错误！');
        end
        ext_len = ext_len + 53;
    elseif dec_point+15<n
%         t = m_start - ext_len;
        dec = str2num(dec_msg(dec_point:dec_point+15));
        bin_str = dec2bin(dec,m_start - ext_len);
        str_lsb = [str_lsb,bin_str];
        dec_point = dec_point + 16;
        if ~strcmp(str_lsb(ext_len:m_start-1),msg(ext_len:m_start-1))
            fprintf('本段秘密消息提取错误！');
        end
        ext_len = m_start;
    else
        dec_rem = n - dec_point;
        bin_rem = floor(log2(10^dec_rem-1));
        if ext_len+bin_rem>m_start
            bin_rem = m_start - ext_len;
        end
        bin_str = dec2bin(str2num(dec_msg(dec_point:dec_point+dec_rem-1)),bin_rem);
        str_lsb = [str_lsb,bin_str];
        dec_point = n;
        if ~strcmp(str_lsb(ext_len:ext_len + bin_rem-1),msg(ext_len:ext_len + bin_rem-1))
            fprintf('本段秘密消息提取错误！');
        end
        ext_len = ext_len + bin_rem;
    end
end

for j=1:ext_len-1
    extr_msg(j) = uint8(str2double(str_lsb(j)));
end
%提取lsb部分嵌入的秘密消息

error_bit = 0;
for j=1:(m_start-1)
    if extr_msg(j)~= uint8(str2double(msg(j)))
        error_bit = error_bit + 1;
    end
end
if error_bit==0
    fprintf('秘密消息提取正确！');
else
    error_rate = error_bit/(m_start-1)*100;
    fprintf('秘密消息提取错误比特率为%.2f%%,',error_rate);
end

psnr = coord_psnr(cover_p,stego_p,row_num,col_num,180,90);
fprintf('PSNR为%.2f,',psnr);
ed = Tra_ed(cover_p,stego_p,row_num);

end

function [cover,costs] = extra_c(cell_p,dig_num,p,row_num,col_num,stc_sum)
    costs = zeros(3,stc_sum,'single');
    s = 1;
    for i=1:row_num
        for j=1:col_num
            if dig_num(s,2)>6
                c = dig_num(s,1)+8;
                value = str2double(cell_p{i,j}(c));
            else
                value = 0;
            end
%             if p(i,j)<0
%                abs_p = abs(p(i,j));
%             else
%                abs_p = p(i,j);
%             end
%             value = mod(floor(abs_p*10^7),10);  %直接提取小数点后第7位
            cover(1,s) = int32(value); %数字放入载体库
            if value==0 || value==9    %数字0和9，不嵌入。
               costs(:,s) = [1e+5 0 1e+5];
            else
               costs(1,s) = emb_distortion(p,7,i,j,-1);   %+1嵌入失真代价
               costs(3,s) = emb_distortion(p,7,i,j,1);   %-1嵌入失真代价
            end
            s = s + 1;
        end
    end
end
       

function cost = emb_distortion(p,stc_k,r_num,c_num,sg)
    a = 0.8;    %调节参数
    if(p(r_num,c_num)<0)
        abs_p = abs(p(r_num,c_num));
        sign = -1;
    else
        abs_p = p(r_num,c_num);
        sign = 1;
    end
    new = sign*(abs_p + sg * 10^(-stc_k));
    if(mod(c_num,2) == 1)
        point_emb = [new p(r_num,2)];   %-1嵌入后的新坐标点
    else
        point_emb = [p(r_num,1) new];
    end
    
    [p_len,~] = size(p);
    if(r_num == 1)
        d = p(r_num+1,:) - p(r_num,:);
        e = p(r_num+1,:) - point_emb;
        x = dot(d,e);
        y = sqrt(sum(d.^2)*sum(e.^2));
        sigma = acos(x/y*a);  %嵌入后新旧向量夹角（与后坐标）弧度制
        cost = rad2deg(sigma);
    elseif(r_num == p_len)
        b = p(r_num-1,:) - p(r_num,:);
        c = p(r_num-1,:) - point_emb;
        sigma = acos(dot(b,c)/sqrt(sum(b.^2)*sum(c.^2))*a); %嵌入后新旧向量夹角（与前坐标）弧度制
        cost = rad2deg(sigma);
    else
        b = p(r_num-1,:) - p(r_num,:);
        c = p(r_num-1,:) - point_emb;
        sigma_front = acos(dot(b,c)/sqrt(sum(b.^2)*sum(c.^2))*a); %嵌入后新旧向量夹角（与前坐标）弧度制
        d = p(r_num+1,:) - p(r_num,:);
        e = p(r_num+1,:) - point_emb;
        sigma_back = acos(dot(d,e)/sqrt(sum(d.^2)*sum(e.^2))*a);  %嵌入后新旧向量夹角（与后坐标）弧度制
        cost = rad2deg(sigma_front + sigma_back);
    end
end
