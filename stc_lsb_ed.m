%stcǶ��С������7λ����8λ�������mlsb�㷨

% clc;
% clear;

function ed = stc_lsb_ed(s_digit,p,msg)

switch s_digit %�������Ч���־������Ƕ����ʣ�
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
s_p = p;    %cell����stego�����
cover_p = zeros(row_num,col_num);
[~,m_len] = size(msg);
h = 10;
dig_num = zeros(row_num*col_num,4);
s = 1;  %���������
lay_lsb = 0;   %���ܾ�����lsb������
lsb_sum = 0;    %lsbλƽ���ܺ�
for i=1:row_num
    for j=1:col_num
        point_str = s_p{i,j};
        point_len = length(point_str);
        point_dec = str2double(s_p{i,j});
        cover_p(i,j) = str2double(s_p{i,j});
        if(point_dec < 0)
           point_num = point_len - 2;  %��Чλ����ȥ�����ź�С����
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

m_start = 1;    %������ϢǶ���ƶ�ָ��
% lay_back = 1;  %����λƽ��
% s_num = 1;  %��s_num������
stc_sum = row_num*col_num;    %stcλƽ���ܺ�
% flag = 1;
%ÿ����53λ������������߿ɱ�ʾ16λʮ����������ʮ���ư���λ����λ�滻��lsbλƽ�档
r = floor(lsb_sum/16);  %�õ�lsb��λƽ�溬�ж�����16��ʮ���ƣ�
lsb_remain = mod(lsb_sum,16);  %16��ʮ����һ��Ƕ���lsb���ֲ���һ���ʮ����λ��
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
        if dec_str_len<16   %53λ�����Ʊ�ʾʮ����λ������16λ�ģ�ǰ�油��0
            add_len = 16 - dec_str_len;
            for i=1:add_len
                dec_str = ['0',dec_str];
            end
        end
        lsb_embstr = [lsb_embstr,dec_str];    %��lsbÿһ�����ܣ�16λ���ϳ��ܵ��ַ���
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
    else  %��stc����Ƕ��
        stc_emb_len = floor(stc_sum*bp);   %stc���ֿ���Ƕ���������Ϣ����
        if m_start+stc_emb_len-1<=m_len
            emb_len = stc_emb_len;
        else
            emb_len = m_len-m_start+1;  %ʣ��λ��С�ڿ�Ƕ������Ϣ����
        end
        [cover,costs] = extra_c(p,dig_num,cover_p,row_num,col_num,stc_sum);   %��ȡstc��������
        for j=1:emb_len
            stc_m(j) = uint8(str2double(msg(m_start+j-1)));
        end
        [d stego n_msg_bits l] = stc_pm1_pls_embed(cover, costs, stc_m, h); % embed message
        for j=1:emb_len
            if stego(j)<0 || stego(j)>9
                fprintf('stego����0-9��\n');
            end
        end
        m_start = m_start + emb_len;
        extr_msg(m_start-emb_len:m_start-1) = stc_ml_extract(stego, n_msg_bits, h); % extract message
        flag = 0;
    end
end
%fprintf('ʵ��Ƕ��%dλ��',m_start-1);
fprintf('��λ����Ƕ��%.2fλ��',(m_start-1)/row_num/col_num);

%�������������
stc_embstr = '';
if exist('stego','var')
    for i=1:stc_sum
        stc_embstr = [stc_embstr,int2str(stego(i))];
    end
end
stego_str = [lsb_embstr,stc_embstr];
stego_len = length(stego_str);
point = 1;  %����Ԫ�أ�lsb+stc��������Ԫ�����
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
if point<=stego_len %����stc����
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

%��ȡlsb����Ƕ���������Ϣ
ext_len = 1;  %��ȡ������Ϣ����
str_lsb = '';
dec_msg = '';
k = 1;  %λƽ�����
n = 1;  %lsb��ȡ�������
c = 1;  %�������
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
dec_point = 1;  %��ȡ��lsb����ʮ���������У�ÿ16λתΪ53λ�����ơ�
while dec_point<n
    if dec_point+15<n&&ext_len+52<m_start
        dec = str2num(dec_msg(dec_point:dec_point+15));
        bin_str = dec2bin(dec,53);
        str_lsb = [str_lsb,bin_str];
        dec_point = dec_point + 16;
        if ~strcmp(str_lsb(ext_len:ext_len+52),msg(ext_len:ext_len+52))
            fprintf('����������Ϣ��ȡ����');
        end
        ext_len = ext_len + 53;
    elseif dec_point+15<n
%         t = m_start - ext_len;
        dec = str2num(dec_msg(dec_point:dec_point+15));
        bin_str = dec2bin(dec,m_start - ext_len);
        str_lsb = [str_lsb,bin_str];
        dec_point = dec_point + 16;
        if ~strcmp(str_lsb(ext_len:m_start-1),msg(ext_len:m_start-1))
            fprintf('����������Ϣ��ȡ����');
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
            fprintf('����������Ϣ��ȡ����');
        end
        ext_len = ext_len + bin_rem;
    end
end

for j=1:ext_len-1
    extr_msg(j) = uint8(str2double(str_lsb(j)));
end
%��ȡlsb����Ƕ���������Ϣ

error_bit = 0;
for j=1:(m_start-1)
    if extr_msg(j)~= uint8(str2double(msg(j)))
        error_bit = error_bit + 1;
    end
end
if error_bit==0
    fprintf('������Ϣ��ȡ��ȷ��');
else
    error_rate = error_bit/(m_start-1)*100;
    fprintf('������Ϣ��ȡ���������Ϊ%.2f%%,',error_rate);
end

psnr = coord_psnr(cover_p,stego_p,row_num,col_num,180,90);
fprintf('PSNRΪ%.2f,',psnr);
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
%             value = mod(floor(abs_p*10^7),10);  %ֱ����ȡС������7λ
            cover(1,s) = int32(value); %���ַ��������
            if value==0 || value==9    %����0��9����Ƕ�롣
               costs(:,s) = [1e+5 0 1e+5];
            else
               costs(1,s) = emb_distortion(p,7,i,j,-1);   %+1Ƕ��ʧ�����
               costs(3,s) = emb_distortion(p,7,i,j,1);   %-1Ƕ��ʧ�����
            end
            s = s + 1;
        end
    end
end
       

function cost = emb_distortion(p,stc_k,r_num,c_num,sg)
    a = 0.8;    %���ڲ���
    if(p(r_num,c_num)<0)
        abs_p = abs(p(r_num,c_num));
        sign = -1;
    else
        abs_p = p(r_num,c_num);
        sign = 1;
    end
    new = sign*(abs_p + sg * 10^(-stc_k));
    if(mod(c_num,2) == 1)
        point_emb = [new p(r_num,2)];   %-1Ƕ�����������
    else
        point_emb = [p(r_num,1) new];
    end
    
    [p_len,~] = size(p);
    if(r_num == 1)
        d = p(r_num+1,:) - p(r_num,:);
        e = p(r_num+1,:) - point_emb;
        x = dot(d,e);
        y = sqrt(sum(d.^2)*sum(e.^2));
        sigma = acos(x/y*a);  %Ƕ����¾������нǣ�������꣩������
        cost = rad2deg(sigma);
    elseif(r_num == p_len)
        b = p(r_num-1,:) - p(r_num,:);
        c = p(r_num-1,:) - point_emb;
        sigma = acos(dot(b,c)/sqrt(sum(b.^2)*sum(c.^2))*a); %Ƕ����¾������нǣ���ǰ���꣩������
        cost = rad2deg(sigma);
    else
        b = p(r_num-1,:) - p(r_num,:);
        c = p(r_num-1,:) - point_emb;
        sigma_front = acos(dot(b,c)/sqrt(sum(b.^2)*sum(c.^2))*a); %Ƕ����¾������нǣ���ǰ���꣩������
        d = p(r_num+1,:) - p(r_num,:);
        e = p(r_num+1,:) - point_emb;
        sigma_back = acos(dot(d,e)/sqrt(sum(d.^2)*sum(e.^2))*a);  %Ƕ����¾������нǣ�������꣩������
        cost = rad2deg(sigma_front + sigma_back);
    end
end
